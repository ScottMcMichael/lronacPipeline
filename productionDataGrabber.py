#!/usr/bin/env python
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

from BeautifulSoup import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing file for LRONAC batch processing '''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------

#TODO: Move this code!
# Retrieves LOLA data from the WUSTL REST web interface
def retrieveLolaFile(minLat, maxLat, minLon, maxLon, outputFolder, padAmount=0.25):

    # Build a query to the WUSTL REST interface for LOLA data
    lolaUrl        = 'http://oderest.rsl.wustl.edu/test/'
    baseQuery      = '?query=lolardr&results=v&'
    locationParams = ('maxlat='      + str(maxLat+padAmount) +
                      '&minlat='     + str(minLat-padAmount) +
                      '&westernlon=' + str(minLon-padAmount) +
                      '&easternlon=' + str(maxLon+padAmount))

    outputPath = os.path.join(outputFolder, 'lolaRdrPoints.csv')
    if os.path.exists(outputPath):
        return True

    queryUrl = lolaUrl + baseQuery + locationParams
    print queryUrl
    
    # Parse the response
    parsedPage = BeautifulSoup(urllib2.urlopen((queryUrl)).read())
    #print parsedPage.prettify()

    # Find the link containing '_pts_csv.csv' and download it
    found = False
        
    for url in parsedPage.findAll('url'):
        if (url.string.find('_pts_csv.csv') >= 0):
            cmd = "wget --output-document=" + outputPath + "  " + url.string
            print cmd
            os.system(cmd)
            found = True
            
    return found



"""
TODO: FUNCTIONS

get input csv file

decide which rows to fetch (all of them?)

for each row:
    
    come up with a name for the pair (name_name?)
     -- Must be reproducible so no pair is done more than once
    
    get path to left and right images (easy)
    
    download all four images
    
    *** How to determine the bounding box?  Start processing one frame and use campt?
    
    Request and download the LOLA file

-- OR --

Pick out only one uncompleted row from the file and grab that data.
 - The submission script will al


"""

#TODO: Move this code!
def getCubeBoundingBox(cubePath, workDir):
    """Returns (minLon, maxLon, minLat, maxLat)"""
    
    # Get the cube size, then request the positions of the four corners
    cubeSize = IsisTools.getCubeSize(cubePath)
    
    # Note that the underlying ISIS tool is one-based
    points  = []
    firstPt =     IsisTools.getPixelLocInCube(cubePath, 1,           1,           workDir)['gdc']
    points.append(IsisTools.getPixelLocInCube(cubePath, cubeSize[0], 1,           workDir)['gdc'])
    points.append(IsisTools.getPixelLocInCube(cubePath, 1,           cubeSize[1], workDir)['gdc'])
    points.append(IsisTools.getPixelLocInCube(cubePath, cubeSize[0], cubeSize[1], workDir)['gdc'])

    # Go through the four corners and get the bounding box
    minLon = firstPt[0]
    maxLon = firstPt[0]
    minLat = firstPt[1]
    maxLat = firstPt[1]
    
    for p in points:
        if p[0] < minLon:
            minLon = p[0]
        if p[0] > maxLon:
            maxLon = p[0]
        if p[1] < minLat:
            minLat = p[1]
        if p[1] > maxLat:
            maxLat = p[1]
            
    return (minLon, maxLon, minLat, maxLat)
    


def retrieveData(inputFile, outputFolder, startLine=0):
    
    # Make sure output folder exists
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    
    # Go through the lines in the file and look for the next unused line
    i = 0
    f = open(inputFile, 'r')
    for line in f:
        
        # Skip the first N lines of the file if requested
        # - Also always skip the first header line
        if i <= startLine:
            i = i + 1
            continue
        i = i + 1
        
        # Read the fields from this line
        strings = line.split(',')
        
        # Generate the name for this pair
        imageA      = strings[0]
        imageB      = strings[1]
        dataSetName = 'pair_' + imageA[:-2] + '_' + imageB[:-2]
        subFolder   = os.path.join(outputFolder, dataSetName+'/')
        
        # If the final log exists in this folder than it is finished and we can skip it
        logPath = os.path.join(subFolder, 'downloadLog.txt')
        if os.path.exists(logPath):
            continue

        IsisTools.createFolder(subFolder) # Create the output folder
        
        # Download the images
        IMAGE_BASE_URL = 'http://lroc.sese.asu.edu/data/'
        imagePathList  = strings[13:17]
        lastDiskPath   = ''
        for image in imagePathList:
            
            # Get the output path
            startPos     = image.rfind('/')
            nameOnDisk   = image[startPos+1:]
            diskPath     = os.path.join(subFolder, nameOnDisk)
            lastDiskPath = diskPath
            
            # Download the image if we don't have it
            if not os.path.exists(diskPath):
                fullUrl = IMAGE_BASE_URL + image
                cmd = "wget --output-document=" + diskPath + "  " + fullUrl
                print cmd
                os.system(cmd)
        
        # Finished downloading the images
        
        # Init the nav data from one cube
        cubePath = lastDiskPath + '.cub'
        cmd = 'lronac2isis from=' + lastDiskPath + ' to=' + cubePath
        os.system(cmd)
        cmd = 'spiceinit from=' + cubePath
        os.system(cmd)
        
        # Get the bounding box of the cube's footprint
        cubeBB = getCubeBoundingBox(cubePath, subFolder)
         
        # Remove the temporary file
        os.remove(cubePath)
         
        # Now fetch the LOLA data, expanding the bounds by a constant to provide wiggle room
        LOLA_PAD_AMOUNT = 0.5
        retrieveLolaFile(cubeBB[2], cubeBB[3], cubeBB[0], cubeBB[1], subFolder, LOLA_PAD_AMOUNT)
        
        # Make a small log file to indicate that all five input files were downloaded.
        logFile = open(logPath, 'w')
        logFile.write(dataSetName + '\n') # Why not store the data set name?
        for image in imagePathList:       # Record the partial paths for the image files
            logFile.write(image + '\n')
        logFile.write('minLon: ' + str(cubeBB[0]) + '\n')  # Record the LOLA bounding box
        logFile.write('maxLon: ' + str(cubeBB[1]) + '\n')
        logFile.write('minLat: ' + str(cubeBB[2]) + '\n')
        logFile.write('maxLat: ' + str(cubeBB[3]) + '\n')
        logFile.write('Finished downloading data: ' + time.strftime("%H:%M:%S"))     # Record the current date/time
        logFile.close()
        
        # Quit after processing a line
        break
    
    f.close()

    # Return a flag value if we did not process any lines from the file
    if not cubeBB:
        return -1
    else: # Return the index of the last line processed plus the folder we wrote to
        return str(i) + ' ' + subFolder 

#--------------------------------------------------------------------------------


def main():

    print "Started productionDataGrabber.py"

    try:
        try:
            usage = "usage: productionDataGrabber.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-i", "--input-file", dest="inputFile",
                              help="Specifies the data list file to read from.",
                              default='')
            parser.add_option("-o", "--output-folder", dest="outputFolder",
                              help="Specifies the folder to copy the data to.",
                              default='./')
            parser.add_option("-s", "--starting-line", dest="startLine",
                              help="Hint to the program to skip ahead past this line.",
                              default=0)
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

        except optparse.OptionError, msg:
            raise Usage(msg)

        #print "Beginning processing....."

        #startTime = time.time()

        outputString = retrieveData(options.inputFile, options.outputFolder, options.startLine)

        #endTime = time.time()

        #print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())