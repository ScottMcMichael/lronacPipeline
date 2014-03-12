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

import os, glob, re, shutil, subprocess, string, time, errno


#TODO: Move some of these to a non-ISIS python file!

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Contains utilities for working with ISIS data files.
'''
    sys.exit()

def createFolder(path):
    """Creates a folder if it does not already exist"""
    if path == '':
        return
    if not os.path.exists(path):
        os.mkdir(path)

def removeIfExists(path):
    """Removes a file if it exists"""
    try:
        os.remove(path)
    except OSError as e: 
        if e.errno != errno.ENOENT: # Continue if the error is "no such file or directory"
            raise # Re-raise the exception if a different error occured

def removeFolderIfExists(directory):
    """Removes a directory and everything in it"""
    try:
        shutil.rmtree(directory)
    except OSError as e: 
        if e.errno != errno.ENOENT: # Continue if the error is "no such file or directory"
            raise # Re-raise the exception if a different error occured

def removeIntermediateStereoFiles(stereoPrefix):
    """Deletes intermediate files from a stereo output directory"""

    # List of all non-final-output files
    fileList = ['-align-L.exr', \
                '-align-R.exr', \
                '-DEM.tif.aux.xml', \
                '-D_sub.tif', \
                '-D.tif', \
                '-F.tif', \
                '-GoodPixelMap.tif', \
                '-lMask_sub.tif', \
                '-lMask.tif', \
                '-L_sub.tif', \
                '-L.tif', \
                '-RD.tif' ,\
                '-rMask_sub.tif', \
                '-rMask.tif', \
                '-R_sub.tif' ,\
                '-R.tif']
    # Remove each of those files
    for f in fileList:
        path = stereoPrefix + f
        removeIfExists(path)

def checkIfToolExists(toolName):
    """Returns true if the system knows about the utility with this name (it is on the PATH)"""

    # Determine the percentage of good pixels   
    cmd = ['which', toolName]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    translateOut, err = p.communicate()

    # Parse the output
    failString = 'no ' + toolName + ' in ('
    if translateOut.find(failString) >= 0:
        raise Exception('Missing requested tool ' + toolName)
    else:
        return True


def parseHeadOutput(textPath, cubePath):
    """Parses the output from head [cube path] and returns a dictionary containing all kernels"""

    kernelDict = dict()

    isisDataFolder = os.environ['ISIS3DATA']
    cubeFolder     = os.path.dirname(cubePath)

    # Search each line in the folder for a required kernel file
    dataFile = open(textPath, 'r')
    lastLine = ''
    kernelsStarted = False
    currentKernelType = 'ERROR!'
    for line in dataFile:
        # Append leftovers from last line and clear left/right whitespace
        workingLine = lastLine + line.strip()
        lastLine = ''

#        print 'workingLine =' + workingLine

        # Skip lines until we find the start of the kernel section
        if (not kernelsStarted ) and (workingLine.find('Group = Kernels') < 0):
            continue
        kernelsStarted = True

        # Quit when we reach the end of the kernel section
        if (workingLine.find('End_Group') >= 0):
            return kernelDict

        # Check if the current line is cut off with an append character
        if (workingLine[-1] == '-'): # This means ISIS has done a weird truncation to the next line
            lastLine = workingLine[:-1] # Strip trailing - and append next line to it next pass
#            print '===   ' + lastLine
            continue

        # Maintain the current kernel type
        if (workingLine.find('LeapSecond') >= 0):
            currentKernelType = 'LeapSecond'
        elif (workingLine.find('TargetAttitudeShape') >= 0):
            currentKernelType = 'TargetAttitudeShape'
        elif (workingLine.find('TargetPosition') >= 0):
            currentKernelType = 'TargetPosition'
        elif (workingLine.find('InstrumentPointing') >= 0):
            currentKernelType = 'InstrumentPointing'
        elif (workingLine.find('InstrumentPosition') >= 0):
            currentKernelType = 'InstrumentPosition'
        elif (workingLine.find('InstrumentAddendum') >= 0):
            currentKernelType = 'InstrumentAddendum'
        elif (workingLine.find('Instrument') >= 0): # This must be check after the other instrument lines!
            currentKernelType = 'Instrument'
        elif (workingLine.find('SpacecraftClock') >= 0):
            currentKernelType = 'SpacecraftClock'
        elif (workingLine.find('ShapeModel') >= 0):
            currentKernelType = 'ShapeModel'

        # Now look for any kernel files on the line 
        # TODO: This will fail if one kernel ends on a line and the next gets a continuation!
        remainingSearchLine = workingLine
        while (len(remainingSearchLine) > 3):
        
            # Look through the line for the next kernel
            m = re.search('[$a-zA-Z0-9/._\-]*'+
                          '((\.tls)|(\.tpc)|(\.tf)|(\.bpc)|(\.bsp)|(\.bc)|(\.tf)|(\.ti)|(\.tsc)|(\.cub))', 
                          remainingSearchLine) 
            
            if not m: # Did not find a kernel
                #print 'Failed to find kernel in line: ' + workingLine
                break # If we did not find a match move on to the next line

            # Found a kernel, handle abbreviations
            if m.group(0)[0] == '$': # Located in ISIS data folder
                kernelPath = os.path.join(isisDataFolder, m.group(0)[1:])
            else: # Path relative to the file location, make it an absolute path
                kernelPath = os.path.join(cubeFolder, m.group(0))

            # Handle special case where two different kinds of files are in the same category
            if (currentKernelType == 'InstrumentPointing') and (kernelPath.find('.tf') > 0):
                currentKernelType = 'Frame'

            # Store the kernel in the dictionary
            if not (currentKernelType in kernelDict):
                kernelDict[currentKernelType] = [kernelPath]
            else:
                kernelDict[currentKernelType].append(kernelPath)

#            print 'In type: ' + currentKernelType + ' Found kernel ' + kernelPath

            # Set up whatever is left of the line for more searching
            remainingSearchLine = remainingSearchLine[m.end()+1:]
#            print remainingSearchLine
#            print len(remainingSearchLine)
#            print '\n'

    # Return the list of kernels
    return kernelDict


def getKernelsFromCube(cubePath, tempFolder):
    """Returns a list of all the SPICE kernels needed by a cube """

    # Call head -120 on file, write to a temp file for parsing
    tempTextPath = os.path.join(tempFolder, "headOutput.txt")
    cmd = "head -120 "+cubePath+" > "+tempTextPath
    #print cmd
    os.system(cmd)
    if not os.path.exists(tempTextPath):
        raise Exception('Failed to extract cube kernel data!')

    # Parse output looking for all the kernel files
    #print 'Looking for source frame file...'
    kernelList = parseHeadOutput(tempTextPath, cubePath)
    if not kernelList:
        raise Exception('Unable to find any kernel files in ' + cubePath)

    return kernelList # Success!


def readPositions(positionFilePath):
    """Reads in a list of GDC coordinates from a pc_align or LOLA RDR file"""

    if not os.path.exists(positionFilePath):
        print 'File ' + positionFilePath + ' is missing!'
        return []

    pointList = []

    #TODO: Read this from the file?
    MEAN_MOON_RADIUS = 1737400

    isLolaFile           = False
    isPcAlignErrorFileKm = False
    isPcAlignErrorFileM  = False
    f = open(positionFilePath, 'r')
    i = 0
    for line in f:
        # On first line check if this is a LOLA RDR file
        if (i == 0):
            if (line.find('Coordinated_Universal_Time') == 0):
                isLolaFile = True
                print 'Detected LOLA RDR file'
                continue # Skip this header line
            if line.find('radius (meters)') > 0:
                isPcAlignErrorFileM = True
                print 'Detected pc_align meters error file'
                continue # Skip this header line
            if  line.find('radius (km)') > 0:
                isPcAlignErrorFileKm = True
                print 'Detected pc_align kilometers error file'
                continue # Skip this header line

        if isLolaFile: # Pick out the correct fields

                strings = line.split(',')
                pointList.append(float(strings[1])) # lon
                pointList.append(float(strings[2])) # lat
                pointList.append(float(strings[3])*1000 - MEAN_MOON_RADIUS) # alt
        
        elif isPcAlignErrorFileM or isPcAlignErrorFileKm: # pc_align error file

                strings = line.split(',')
                pointList.append(float(strings[0])) # lon
                pointList.append(float(strings[1])) # lat
                if isPcAlignErrorFileM:
                    pointList.append(float(strings[2]) - MEAN_MOON_RADIUS) # alt
                else: # Units are km
                    pointList.append(float(strings[2])*1000.0 - MEAN_MOON_RADIUS) # alt

        else: # Default handling
            if line.find('#') < 0: # Skip lines containing the comment symbol
                strings = line.split(',')
                #print strings
                pointList.append(float(strings[1])) # lon
                pointList.append(float(strings[0])) # lat
                pointList.append(float(strings[2])) # alt
        i = i + 1
    f.close()

    #print pointList
    return pointList


def readJitregFile(filePath):
    """Reads the output file from a lronacjitreg call and returns [meanSampleOffset, meanLineOffset]"""

    # Fail if the input file is not present
    if not os.path.isfile(filePath):
        raise Exception('File ' + filePath + ' is missing!')

    averages = [0.0, 0.0]

    f = open(filePath,'r')
    for line in f:
        if ( line.rfind("Average Sample Offset:") >= 0 ):
            index       = line.rfind("Offset:");
            index_e     = line.rfind("StdDev:");
            crop        = line[index+7:index_e];
            if crop == " NULL ": # Check for null value
                raise Exception('Null sample offset in file ' + flat)
            averages[0] = float(crop);
        elif ( line.rfind("Average Line Offset:") >= 0 ):
            index       = line.rfind("Offset:");
            index_e     = line.rfind("StdDev:");
            crop        = line[index+7:index_e];
            if crop == "   NULL ": # Check for null value
                raise Exception('Null sample offset in file ' + flat)
            averages[1] = float(crop);
        elif ( line.rfind("Using IpFind result only:") >= 0 ):
            index       = line.rfind("only:");
            if (line[index + 7] == 1):
                print "Warning: This result based only on IpFind search."
    print str(averages)
    return averages



def writeLronacPvlFile(outputPath, isHalfRes):
    """
    Generates a .pvl file needed to use noproj with an LRONAC camera pair.
    - Can generate a version for either full or half sample resolution files.
    """

    if os.path.exists(outputPath):
        print outputPath + ' already exists, using existing file.'
        return True
    else: # Need to write the file
        print 'Generating LRONAC compatible .pvl file ' + outputPath

    f = open(outputPath, 'w')

    f.write('Object = IdealInstrumentsSpecifications\n');
    f.write('  UserName     = auto\n');
    f.write('  Created      = 2013-07-18T13:42:00\n');
    f.write('  LastModified = 2013-07-18T13:42:00\n\n');
    f.write('  Group = "LUNAR RECONNAISSANCE ORBITER/NACL"\n');

    if not isHalfRes: # Full resolution camera
        f.write('     TransY = 16.8833\n')
        f.write('     ItransS = -2411.9\n')
        f.write('     TransX = 0.6475\n')
        f.write('     ItransL = -92.5\n')
        f.write('     DetectorSamples = 10000\n')
    else: # Half resolution camera
        f.write('     TransY = 16.8833\n')
        f.write('     ItransS = -4823.8\n')     # Halved
        f.write('     TransX = 0.6475\n')
        f.write('     ItransL = -185\n')       # Halved
        f.write('     DetectorSamples = 5000\n') # Halved

    f.write('  End_Group\n\n')
    f.write('End_Object\n')
    f.write('End')

    f.close()

def getCubeCenterLatitude(cubePath, workDir='tmp'):
    """Calls caminfo on a mosaic and returns the CenterLatitude value"""

    # Make sure the requested file is present
    if not os.path.exists(cubePath):
        raise Exception('File ' + cubePath + ' does not exist!')

    # Call caminfo (from ISIS) on the input cube to find out the CenterLatitude value
    camInfoOuputPath = workDir + "/camInfoOutput.txt"
    cmd = 'caminfo from=' + cubePath + ' to=' + camInfoOuputPath
    os.system(cmd)

    if not os.path.exists(camInfoOuputPath):
        raise Exception('Call to caminfo failed on file ' + cubePath)

    # Read in the output file to extract the CenterLatitude value
    centerLatitude = -9999
    infoFile       = open(camInfoOuputPath, 'r')
    for line in infoFile:
        if (line.find('CenterLatitude') >= 0):
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            centerLatitude = float(numStr)
            break
    # Make sure we found the desired value
    if (centerLatitude == -9999):          
        raise Exception("Unable to find CenterLatitude from file " + cubePath)

    # Clean up temporary file
    os.remove(camInfoOuputPath)

    return centerLatitude # Return the latitude we found



def makeSpkSetupFile(leapSecondFilePath, outputPath):
    """Creates the required mkspk setup file if it does not already exist"""

    # If the file already exists, delete it and rewrite it.
    if os.path.exists(outputPath):
        os.remove(outputPath)

#        print 'Generating LRONAC compatible .pvl file ' + halfResFilePath
    f = open(outputPath, 'w')
    f.write("\\begindata\n")
    f.write("INPUT_DATA_TYPE   = 'STATES'\n")
    f.write("OUTPUT_SPK_TYPE   = 13\n")
    f.write("OBJECT_ID         = -85\n") # LRO
    f.write("CENTER_ID         = 301\n") # Moon
    f.write("REF_FRAME_NAME    = 'J2000'\n")
    f.write("PRODUCER_ID       = 'Lronac Pipeline'\n")
    f.write("DATA_ORDER        = 'epoch x y z vx vy vz'\n")
    f.write("DATA_DELIMITER    = ','\n")
    f.write("LEAPSECONDS_FILE  = '" + leapSecondFilePath + "'\n")
    f.write("LINES_PER_RECORD  = 1\n")
    f.write("TIME_WRAPPER      = '# ETSECONDS'\n")
    #f.write("EPOCH_STR_LENGTH  = 16\n")
    f.write("INPUT_DATA_UNITS  = ('ANGLES=DEGREES' 'DISTANCES=km')\n")
    f.write("POLYNOM_DEGREE    = 11\n")
    f.write("SEGMENT_ID        = 'SPK_STATES_13'\n")
#        f.write("INPUT_DATA_FILE   = 'spkDataFile.txt'")
#        f.write("OUTPUT_SPK_FILE   = '/home/smcmich1/testSpkFile.bsp'")
    f.write("\\begintext\n")
    f.close()


# TODO: Merge this with IsisTools functions
def getPixelLocInCube(cubePath, sample, line, workDir=''):
    """Returns the BodyFixedCoordinate of a pixel from a cube"""

    DEFAULT_MOON_RADIUS = 1737400 # In meters

    # Make sure the input file exists
    if not os.path.exists(cubePath):
        raise Exception('Cube file ' + cubePath + ' not found!')

    # Default working directory is the cubePath folder
    outputFolder = workDir
    if workDir == '':
        outputFolder = os.path.dirname(cubePath)
       
    if (len(outputFolder) > 1) and not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    # Call ISIS campt function to compute the pixel location
    tempTextPath = os.path.join(outputFolder, 'camptOutput.txt')
    if os.path.exists(tempTextPath):
        os.remove(tempTextPath) # Make sure any existing file is removed!
        
    # Use subprocess to suppress the command output
    cmd = ['campt', 'from=', cubePath, 'to=', tempTextPath, 'sample=', str(sample), 'line=', str(line)]
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
    #subprocess.call(cmd)

    # Check that we created the temporary file
    if not os.path.exists(tempTextPath):
        raise Exception('campt failed to create temporary file ' + tempTextPath)
    
    infoFile = open(tempTextPath, 'r')
        
    # Read in the output file to extract the pixel coordinates
    gccLine       = ''
    latLine       = ''
    lonLine       = ''
    radiusLine    = ''
    lineAfterBody = False
    for line in infoFile:
        
        # GCC stuff
        if lineAfterBody: # BodyFixedCoordinate takes up two lines
            gccLine       = gccLine + line
            lineAfterBody = False
            
        if (gccLine == ''): # Look for start of the info (this check must come second)
            if (line.find('BodyFixedCoordinate') >= 0):
                gccLine     = line
                lineAfterBody = True
        
        # GDC stuff
        if line.find('PlanetocentricLatitude') >= 0:
            latLine = line
            #print line
        if line.find('PositiveEast180Longitude') >= 0:
            lonLine = line
            #print line
        if line.find('LocalRadius') >= 0:
            radiusLine = line
            #print line

    os.remove(tempTextPath) # Remove the file to clean up

    # Make sure we found the desired lines
    if (gccLine == ''):
        raise Exception("Unable to find BodyFixedCoordinate in file " + tempTextPath)
    if (latLine == ''):
        raise Exception("Unable to find PlanetocentricLatitude in file " + tempTextPath)
    if (lonLine == ''):
        raise Exception("Unable to find PositiveEast180Longitude in file " + tempTextPath)
    if (radiusLine == ''):
        raise Exception("Unable to find LocalRadius in file " + tempTextPath)

    # Extract GCC coordinates
    startParen = gccLine.find('(')
    stopParen  = gccLine.find(')')
    numString  = gccLine[startParen+1:stopParen]
    x,y,z = numString.split(',')

    # Convert output from kilometers to meters
    pixelLocationGcc = [0, 0, 0]
    pixelLocationGcc[0] = float(x) * 1000.0
    pixelLocationGcc[1] = float(y) * 1000.0
    pixelLocationGcc[2] = float(z) * 1000.0
    
    # Extract GDC coordinates
    latStart     = latLine.find('=')+2
    lonStart     = lonLine.find('=')+2
    radiusStart  = radiusLine.find('=')+2
    radiusEnd    = radiusLine.find('<') - 1
    latNumStr    = latLine[latStart:]
    lonNumStr    = lonLine[lonStart:]
    radiusNumStr = radiusLine[radiusStart:radiusEnd]
    pixelLocationGdc = [float(lonNumStr), float(latNumStr), float(radiusNumStr)-1737400]
                        
    pixelInformation = dict()
    pixelInformation['gcc'] = pixelLocationGcc
    pixelInformation['gdc'] = pixelLocationGdc
    return pixelInformation



def getCubeElevationEstimate(cubePath, workDir=''):
    """Returns the surface elevation at the center of a cube"""

    DEFAULT_MOON_RADIUS = 1737400 # In meters

    # TODO: Get these values from the file!
    sample = 2500
    line   = 25000

    # Make sure the input file exists
    if not os.path.exists(cubePath):
        raise Exception('Cube file ' + cubePath + ' not found!')

    # Default working directory is the cubePath folder
    outputFolder = workDir
    if workDir == '':
        outputFolder = os.path.dirname(cubePath)
       
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    # Call ISIS campt function to compute the pixel location
    tempTextPath = os.path.join(outputFolder, 'camptOutput.txt')
    if os.path.exists(tempTextPath):
        os.remove(tempTextPath) # Make sure any existing file is removed!
        
    # Use subprocess to suppress the command output
    cmd = ['campt', 'from=', cubePath, 'to=', tempTextPath, 'sample=', str(sample), 'line=', str(line)]
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)

    # Check that we created the temporary file
    if not os.path.exists(tempTextPath):
        raise Exception('campt failed to create temporary file ' + tempTextPath)
        
    # Read in the output file to extract the pixel coordinates
    foundLine   = ''
    infoFile    = open(tempTextPath, 'r')
    for line in infoFile:
        if (line.find('LocalRadius') >= 0):
            foundLine = line
            break

    os.remove(tempTextPath) # Remove the file to clean up

    # Make sure we found the desired lines
    if (foundLine == ''):
        raise Exception("Unable to find LocalRadius in file " + tempTextPath)

    # ExtractfoundLine the desired coordinates
    eqPos     = foundLine.find('=')
    endPos    = foundLine.find('<')
    numString = foundLine[eqPos+2:endPos-2]

    # Convert the absolute radius into a height relative to the mean radius of the moon
    localRadius = float(numString) - DEFAULT_MOON_RADIUS
    print 'found local radius ' + str(localRadius)

    return localRadius


def getStereoGoodPixelPercentage(inputPrefix, workDir=''):
    """Returns the percentage of good pixels in a stereo output"""

    # Set up input folder
    inputFolder = os.path.dirname(inputPrefix)
    if not os.path.exists(inputFolder):
        raise Exception('Input folder ' + inputFolder + ' not found!')    
    if workDir == '':
        workDir = inputFolder

    
    #TODO: Look for goodPixelMap file!
    
    #TODO: Look for later stage estimates!
    
    # If the later stage files were not found, use the integer correlation file 
    
    # Extract the third band of the D_sub.tif image which contains a good pixel map
    inputPath = inputPrefix + '-D_sub.tif'
    if not os.path.exists(inputPath):
        raise Exception('Could not find file ' + inputPath)
    convertedImagePath = os.path.join(workDir,     'goodPixelMap-D_sub.tif')
    cmd = 'gdal_translate -of GTiff -ot BYTE -b 3 ' + inputPath + ' ' + convertedImagePath
    print cmd
    os.system(cmd)
    
    # Determine the percentage of good pixels   
    cmd = ['gdalinfo', '-hist', convertedImagePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    translateOut, err = p.communicate()

    # Parse the gdalinfo output
    bucket  = translateOut.find('buckets')
    colon   = translateOut.find(':', bucket)
    start   = translateOut.find('\n', colon)
    end     = translateOut.find('\n', start+1)
    buckets = translateOut[start+1:end] # Pick off the string containing the buckets
    numbers = buckets.strip().split(' ')
    
    numBad      = int(numbers[0]) # All pixels are in the first (bad) or last (good) buckets
    numGood     = int(numbers[-1])
    percentGood = float(numGood) / float(numGood + numBad)

    return percentGood
  


def getCubeSize(cubePath):
    """Returns the size [samples, lines] in a cube"""

    # Make sure the input file exists
    if not os.path.exists(cubePath):
        raise Exception('Cube file ' + cubePath + ' not found!')
       
    # Use subprocess to suppress the command output
    cmd = ['gdalinfo', cubePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    # Extract the size from the text
    sizePos    = textOutput.find('Size is')
    endPos     = textOutput.find('\n', sizePos+7)
    sizeStr    = textOutput[sizePos+7:endPos]
    sizeStrs   = sizeStr.strip().split(',')
    numSamples = int(sizeStrs[0])
    numLines   = int(sizeStrs[1])
    
    size = [numSamples, numLines]
    return size

def getCubeBoundingBox(cubePath, workDir):
    """Returns (minLon, maxLon, minLat, maxLat)"""
    
    # Get the cube size, then request the positions of the four corners
    cubeSize = getCubeSize(cubePath)
    
    # Note that the underlying ISIS tool is one-based
    points  = []
    firstPt =     getPixelLocInCube(cubePath, 1,           1,           workDir)['gdc']
    points.append(getPixelLocInCube(cubePath, cubeSize[0], 1,           workDir)['gdc'])
    points.append(getPixelLocInCube(cubePath, 1,           cubeSize[1], workDir)['gdc'])
    points.append(getPixelLocInCube(cubePath, cubeSize[0], cubeSize[1], workDir)['gdc'])

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


def modifyPixelPairs(inputPath, outputPath, leftOffsetX, leftOffsetY, rightOffsetX, rightOffsetY):
    """Adds offsets to the pixel pairs in a file on disk"""
    
    # Make sure the input file exists
    if not os.path.exists(inputPath):
        raise Exception('Pixel pair file ' + inputPath + ' not found!')
    
    
    # Open files
    fI = open(inputPath,  'r')
    fO = open(outputPath, 'w')
    
    
    for line in fI:
        # Get and modify input values
        values = line.split(',')
        leftX  = float(values[0]) + leftOffsetX
        leftY  = float(values[1]) + leftOffsetY
        rightX = float(values[2]) + rightOffsetX
        rightY = float(values[3]) + rightOffsetY
        
        # Write output values
        fO.write("%f,%f,%f,%f\n" % (leftX, leftY, rightX, rightY))
        
    # Close files
    fI.close()
    fO.close()
    
    return True


# Reads the results from a single diff_stats.txt file
def readLolaCompareFile(filePath):

    # Read in the output file to extract the CenterLatitude value
    meanValue  = -32768
    stdDev     = -32768
    percentile = []
    histogram  = []

    statsFile = open(filePath, 'r')
    for line in statsFile:
        if (line.find('Mean') >= 0):         # Parse out the mean
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            meanValue = float(numStr)
        elif (line.find('deviation') >= 0):  # Parse out the standard deviation
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            stdDev = float(numStr)
        elif (line.find('Percentile') >= 0): # Add a line to the percentile distribtion
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            percentile.append(float(numStr))
        elif (line.find('<-->') >= 0): # Add a line to the histogram
            eqPt   = line.rfind('=') # Get the second equal sign
            perPt  = line.find('%')
            numStr = line[eqPt+1:perPt-1]
            histogram.append(float(numStr))

    # Make sure we found the desired values
    if (meanValue == -32768) or (stdDev == -32768):
        raise Exception("Unable to find statistics in file " + filePath)
    
    return (meanValue, stdDev, percentile, histogram)

def makeDataSetName(fileNameA, fileNameB):
    """Given two of the input images, determines a data set name"""

    if fileNameA < fileNameB: # Always put the image with the lower number first
        dataSetName = 'NAC_DTM_' + fileNameA[:-2] + '_' + fileNameB[:-2]
    else:
        dataSetName = 'NAC_DTM_' + fileNameB[:-2] + '_' + fileNameA[:-2]

    return dataSetName

#TODO: Make this one general?
# Copies one file from the supercomputer - Caller needs to wait for the job to finish!
def grabSupercomputerFile(supercomputerPath, localPath):

    cmd = 'sup scp smcmich1@bridge3.nas.nasa.gov:' + supercomputerPath  + ' ' + localPath
    os.system(cmd)


def getLastGitTag(codePath):
    """Returns the last brief git tag of the repository containing the file""" 

    # Get path to git folder
    codeFolder = os.path.dirname(codePath)
    gitFolder  = os.path.join(codeFolder, '.git/')
    
    # Get the tag using a subprocess call
    cmd = ['git', '--git-dir', 'gitFolder', 'describe', '--abbrev=0']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    return textOutput


#==================================================
# This class implements a variant of OptionParser which ignores unknown options.

# TODO: Replace existing implementation in stereo_utils.py?
#       --> This version does not make the user treat the arguments differently than normal.

from optparse import (OptionParser,BadOptionError,AmbiguousOptionError)


class PassThroughOptionParser(OptionParser):

    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                self._process_args2(largs,rargs,values) 
            except (BadOptionError,AmbiguousOptionError), e:  # On failure, pass option to output list
                if sys.version_info < (2, 6, 0):
                    # Port to Python 2.4
                    p = re.match("^.*?no such option:\s*(.*?)$", e.msg)
                    if p:
                        largs.append(p.group(1))
                else:
                    largs.append(e.opt_str)
    

    # This version of the function successfully passes through negative numbers    
    def _process_args2(self, largs, rargs, values):
        """_process_args(largs : [string],
                         rargs : [string],
                         values : Values)

        Process command-line arguments and populate 'values', consuming
        options and arguments from 'rargs'.  If 'allow_interspersed_args' is
        false, stop at the first non-option argument.  If true, accumulate any
        interspersed non-option arguments in 'largs'.
        """
        while rargs:
            arg = rargs[0]

            p = re.match('^-[0-9.]+$', arg)
            if p:
                del rargs[0]
                raise BadOptionError(arg)
                #self.error(_("%s unrecognized number in arguments") % arg)
            
            # We handle bare "--" explicitly, and bare "-" is handled by the
            # standard arg handler since the short arg case ensures that the
            # len of the opt string is greater than 1.
            if arg == "--":
                del rargs[0]
                return
            elif arg[0:2] == "--":
                # process a single long option (possibly with value(s))
                OptionParser._process_long_opt(self, rargs, values)
            elif arg[:1] == "-" and len(arg) > 1:
                # process a cluster of short options (possibly with
                # value(s) for the last one only)
                OptionParser._process_short_opts(self, rargs, values)
            elif self.allow_interspersed_args:
                largs.append(arg)
                del rargs[0]
            else:
                return                  # stop now, leave this arg in rargs

        # Say this is the original argument list:
        # [arg0, arg1, ..., arg(i-1), arg(i), arg(i+1), ..., arg(N-1)]
        #                            ^
        # (we are about to process arg(i)).
        #
        # Then rargs is [arg(i), ..., arg(N-1)] and largs is a *subset* of
        # [arg0, ..., arg(i-1)] (any options and their arguments will have
        # been removed from largs).
        #
        # The while loop will usually consume 1 or more arguments per pass.
        # If it consumes 1 (eg. arg is an option that takes no arguments),
        # then after _process_arg() is done the situation is:
        #
        #   largs = subset of [arg0, ..., arg(i)]
        #   rargs = [arg(i+1), ..., arg(N-1)]
        #
        # If allow_interspersed_args is false, largs will always be
        # *empty* -- still a subset of [arg0, ..., arg(i-1)], but
        # not a very interesting subset!


