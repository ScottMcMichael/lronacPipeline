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

import IrgFileFunctions, IrgGeoFunctions, IrgIsisFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Contains utilities required for the LRONAC mass pipeline.
'''
    sys.exit()

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


def makeSpkSetupFile(leapSecondFilePath, outputPath):
    """Creates the required mkspk setup file if it does not already exist"""

    # If the file already exists, delete it and rewrite it.
    IrgFileFunctions.removeIfExists(outputPath)

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

    justFileA = os.path.splitext(os.path.basename(fileNameA))[0] # Strip paths and extensions
    justFileB = os.path.splitext(os.path.basename(fileNameB))[0]

    if justFileA < justFileB: # Always put the image with the lower number first
        dataSetName = 'NAC_DTM_' + justFileA[:-2] + '_' + justFileB[:-2]
    else:
        dataSetName = 'NAC_DTM_' + justFileB[:-2] + '_' + justFileA[:-2]

    return dataSetName

#TODO: Make this one general?
# Copies one file from the supercomputer - Caller needs to wait for the job to finish!
def grabSupercomputerFile(supercomputerPath, localPath):

    cmd = 'sup scp smcmich1@pfe.nas.nasa.gov:' + supercomputerPath  + ' ' + localPath
    os.system(cmd)
    
    
    
def cam2mapWithDem(inputCubePath, demPath, projectedCubePath):
    """Calls cam2map with a user DEM"""
    
    # TODO: Make sure this does not overwrite any custom NAV data!
    
    # Load the user DEM in to the input cube
    cmd = ['spiceinit', 'from=', inputCubePath, 'shape=', 'USER', 'model=', demPath]
    
    # Call the ASP wrapper of the ISIS cam2map function
    cmd = ['cam2map4stereo.py']
    
    

    
    
    


