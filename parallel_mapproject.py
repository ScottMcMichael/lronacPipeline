#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

import sys

import os, glob, re, shutil, subprocess, string, time, errno, optparse, math

import IrgFileFunctions, IrgIsisFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Calls mapproject in parallel with ISIS cameras.
'''
    sys.exit()

# List of currently running jobs
job_pool  = []


def add_job( cmd, suppressTileOutput=True, num_working_threads=4 ):
    
    FNULL = open(os.devnull, 'w')
    
    # If we already have too many running processes:
    if ( len(job_pool) >= num_working_threads):
        job_pool[0].wait(); # Wait until the first one finishes
        job_pool.pop(0); # Pull the finished job off the list
        
    # Add the new job to the end of the list
    if suppressTileOutput:
        job_pool.append( subprocess.Popen(cmd, stdout=FNULL, stderr=subprocess.STDOUT) )
    else:
        job_pool.append( subprocess.Popen(cmd) )

def wait_on_all_jobs():
    print "Waiting for jobs to finish";
    while len(job_pool) > 0:
        job_pool[0].wait();
        job_pool.pop(0);


def isOption(arg):
    """Returns True if the string is an argument, False otherwise"""
    
    # An option must start with '-' and not consist of all numbers
    if ( arg.startswith('-') and not re.match('^-[0-9.]+$', arg) ):
        return True
    else:
        return False

def generateTileList(fullWidth, fullHeight, tileSize):
    """Generate a full list of tiles for this image"""

    numTilesX = int(math.ceil(fullWidth  / float(tileSize)))
    numTilesY = int(math.ceil(fullHeight / float(tileSize)))

    tileList = []
    for r in range(0, numTilesY):
        for c in range(0, numTilesX):
            
            # Starting pixel positions for the tile
            tileStartY = r * tileSize
            tileStartX = c * tileSize
            
            # Determine the size of this tile
            thisWidth  = tileSize
            thisHeight = tileSize
            if (r == numTilesY-1): # If the last row
                thisHeight = fullHeight - tileStartY # Height is last remaining pixels
            if (c == numTilesX-1): # If the last col
                thisWidth  = fullWidth  - tileStartX # Width is last remaining pixels
            
            # Get the end pixels for this tile
            tileStopY  = tileStartY + thisHeight # Stop values are exclusive
            tileStopX  = tileStartX + thisWidth
            
            # Create a name for this tile
            # - Tile format is tile_col_row_width_height_.tif
            tileString = 'tile_' + str(c) + '_' + str(r) + '_' + str(thisWidth) + '_' + str(thisHeight) + '_.tif'
            
            tileList.append(tileStartX, tileStartY, thisWidth, thisHeight, tileString)
    
    return (numTilesX, numTilesY, tileList)

def handleArguments(args):
    """Split up arguments into required and optional lists which will be passed to subprocess"""

    requiredList = []
    optionsList  = []
       
    # Loop through all entries.
    iterable = iter(range(0, len(args)))
    for i in iterable:
        a = args[i]
        if (i < len(args)-1): # Don't load the next value when we are at the end!
            n = args[i+1]
        else:
            n = '-' # This will just cause us to get out of the loop
        
        if isOption(a):      # This is the start of an option.
            optionsList.append(a)  # Record this entry.
            
            if isOption(n):  # The next entry is the start of another option so this one has no values.
                continue
            
            optionsList.append(n)  # Otherwise record the next entry as a value.
            iterable.next()              # Skip the next entry in the loop.

            if (a == '--t_projwin') or (a == '--t_pixelwin'):  # These arguments have four values, not just one.
                optionsList.append(args[i+2])              # Add the additional three arguments and skip them in the loop.
                optionsList.append(args[i+3])
                optionsList.append(args[i+4])
                iterable.next()
                iterable.next()
                iterable.next()
        
        else: # This is one of the three positional arguments
            requiredList.append(a)
    
    # Return the two lists
    return (requiredList, optionsList)


#------------------------------------------------------------------------------

def main(argsIn):


    try:
        usage = "usage: parallel_mapproject.py [options] <dem> <camera-image> <output>"
        parser = IrgFileFunctions.PassThroughOptionParser(usage=usage) # Use parser that ignores unknown options

        parser.set_defaults(numThreads=8)
        parser.set_defaults(keep=False)
        parser.set_defaults(suppressOutput=False)

        parser.add_option("--num-threads",  dest="numThreads",  help="Number of threads to use")
        
        parser.add_option("--convert-tiles",  action="store_true",
                          dest="convertTiles",  help="Generate a uint8 version of each tile")
        parser.add_option("--suppress-output",  action="store_true",
                          dest="suppressOutput",  help="Suppress output of sub-calls.")

        parser.add_option("--manual", action="callback", callback=man,
                          help="Read the manual.")
        parser.add_option("--keep", action="store_true", dest="keep",
                          help="Do not delete the temporary files.")
        
        # This call handles all the parallel_mapproject specific options.
        (options, args) = parser.parse_args(argsIn)

        # This will parse all the mapproject options.
        requiredList, optionsList = handleArguments(args)

        # Check the required positional arguments.
        if len(requiredList) < 1: 
            parser.error("Need path to input image")
        if len(requiredList) < 2: 
            parser.error("Need path to DEM")
        if len(requiredList) < 3:
            parser.error("Need output path")
        
        options.imagePath  = requiredList[0]
        options.demPath    = requiredList[1]
        options.outputPath = requiredList[2]

        # Any additional arguments need to be forwarded to the mapproject function
        options.extraArgs = optionsList

    except optparse.OptionError, msg:
        raise Usage(msg)

    startTime = time.time()
    
    # If the input image is NOT an ISIS image then the normal map_project call
    #  can operate in parallel without this wrapper.
    if not IrgIsisFunctions.isIsisFile(options.imagePath):
        cmd = ['mapproject',  options.imagePath, options.demPath, options.outputPath]    
        cmd = cmd + extraArgs
        subprocess.call(cmd)
        return 0
    

    # Call mapproject on the input data using subprocess and record output
    cmd = ['mapproject',  '--query-projection', options.imagePath, options.demPath, options.outputPath]
    cmd = cmd + options.extraArgs # Append other options
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    projectionInfo, err = p.communicate()
    if not options.suppressOutput:
        print projectionInfo
        
        
    # Now find the image size in the output
    startPos    = projectionInfo.find('Output image bounding box:')
    widthStart  = projectionInfo.find('width:', startPos)
    heightStart = projectionInfo.find('height:', widthStart)
    heightEnd   = projectionInfo.find(')', heightStart)
    # Extract values
    fullWidth   = int(projectionInfo[widthStart+7  : heightStart-1])
    fullHeight  = int(projectionInfo[heightStart+8 : heightEnd])
    print 'Output image size is ' + str(fullWidth) + ' by ' + str(fullHeight) + ' pixels.'

    # Figure out how to break up the image into tiles.
    # - For now just do something simple.
    TILE_SIZE = 1000
    
    numTilesX, numTilesY, tileList = generateTileList(fullWidth, fullHeight, TILE_SIZE)
    
    
    print 'Splitting into ' + str(numTilesX) + ' by ' + str(numTilesY) + ' tiles.'

    # Make a temporary directory to store the tiles
    outputFolder = os.path.dirname(options.outputPath)
    if outputFolder == '':
        outputFolder = './' # Handle calls in same directory
    outputName   = os.path.basename(options.outputPath)
    IrgFileFunctions.createFolder(outputFolder)
    tempFolder   = os.path.join(outputFolder, outputName.replace('.', '_') + '_tiles/')
    IrgFileFunctions.createFolder(tempFolder)
    
    
    # Queue up one mapproject call for each file
    print 'Writing tiles...'
    tilePathList = []
    index = 0    
    for r in range(0, numTilesY):
        for c in range(0, numTilesX):
            
            # Starting pixel positions for the tile
            tileStartY = tiles[index][0]
            tileStartX = tiles[index][1]
                        
            # Get the end pixels for this tile
            tileStopY  = tileStartY + tiles[index][2] # Stop values are exclusive
            tileStopX  = tileStartX + tiles[index][3]
            
            # Get the output path for this tile
            tilePath = os.path.join(tempFolder, tiles[index][4])
    
            print 'Writing tile: ' + tiles[index][4]
            
            # Call mapproject on the input data using subprocess and record output
            cmd = ['mapproject',  '--t_pixelwin', str(tileStartX), str(tileStartY), str(tileStopX), str(tileStopY),
                                   options.imagePath, options.demPath, tilePath]
            cmd = cmd + options.extraArgs # Append other options
            if not os.path.exists(tilePath):
                add_job(cmd, options.suppressOutput, int(options.numThreads)) # Send to parallel job queue
    
            index = index + 1
    
    # Wait for all of the tiles to finish processing
    wait_on_all_jobs()
    
    if options.convertTiles: # Make uint8 version of all tiles for debugging
        print 'Writing out uint8 version of all tiles...'
        for t in tiles:
            tilePath   = os.path.join(tempFolder, t[4])
            tilePathU8 = os.path.splitext(tilePath)[0] + 'U8.tif'
            cmd = ['gdal_translate', '-ot', 'byte', '-scale', tilePath, tilePathU8]
            add_job(cmd, options.suppressOutput, options.numThreads) # Send to parallel job queue
    
        # Wait for all of the tiles to finish processing
        wait_on_all_jobs()
    
    # Build a gdal VRT file which is composed of all the processed tiles
    vrtPath = os.path.join(tempFolder, 'mosaic.vrt')
    cmd = "gdalbuildvrt  -resolution highest " + vrtPath + " " + tempFolder + "*_.tif";
    print cmd
    os.system(cmd)

    # Convert VRT file to final output file
    cmd = "gdal_translate -co compress=lzw " + vrtPath + " " + options.outputPath;
    print cmd
    os.system(cmd)

    # Clean up temporary files
    if not options.keep:
        IrgFileFunctions.removeFolderIfExists(tempFolder)


    endTime = time.time()

    print "Finished in " + str(endTime - startTime) + " seconds."

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
