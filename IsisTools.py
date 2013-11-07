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

import os, glob, optparse, re, shutil, subprocess, string, time


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Contains utilities for working with ISIS data files.
'''
    sys.exit()


def parseHeadOutput(textPath, cubePath):
    """Parses the output from head [cube path] and returns a list of all kernels"""

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

        #print 'workingLine =' + workingLine

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

        # Now look for any kernel files on the line (should never be more than one per line)
        m = re.search('[$a-zA-Z0-9/._\-]*((\.tls)|(\.tpc)|(\.tf)|(\.bpc)|(\.bsp)|(\.bc)|(\.tf)|(\.ti)|(\.tsc)|(\.cub))', workingLine) 
        if not m:
            #print 'Failed to find kernel in line: ' + workingLine
            continue # If we did not find a match move on to the next line

        if m.group(0)[0] == '$': # Located in ISIS data folder
            kernelPath = os.path.join(isisDataFolder, m.group(0)[1:])
        else: # Path relative to the file location, make it an absolute path
            kernelPath = os.path.join(cubeFolder, m.group(0))

        # Handle special case where two different kinds of files are in the same category
        if (currentKernelType == 'InstrumentPointing') and (kernelPath.find('.tf') > 0):
            currentKernelType = 'Frame'

        #print 'In type: ' + currentKernelType + ' Found kernel ' + kernelPath

        if not (currentKernelType in kernelDict):
            kernelDict[currentKernelType] = [kernelPath]
        else:
            kernelDict[currentKernelType].append(kernelPath)

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

# TODO: Move this to another module!
def readPositions(positionFilePath):
    """Reads in a list of GDC coordinates from a pc_align or LOLA RDR file"""

    if not os.path.exists(positionFilePath):
        print 'File ' + positionFilePath + ' is missing!'
        return []

    pointList = []

    MEAN_MOON_RADIUS = 1737400

    isLolaFile = False
    f = open(positionFilePath, 'r')
    i = 0
    for line in f:
        # On first line check if this is a LOLA RDR file
        if (i == 0) and (line.find('Coordinated_Universal_Time') == 0):
            isLolaFile = True
            print 'Detected LOLA RDR file'
            continue # Skip this header line

        if isLolaFile: # Pick out the correct fields

                strings = line.split(',')
                pointList.append(float(strings[1])) # lon
                pointList.append(float(strings[2])) # lat
                pointList.append(float(strings[3])*1000 - MEAN_MOON_RADIUS) # alt

        else: # Not a LOLA RDR file
            if line.find('#') < 0: # Skip lines containing the comment symbol
                strings = line.split(',')
                #print strings
                pointList.append(float(strings[1])) # lon
                pointList.append(float(strings[0])) # lat
    #            if (pcAlign):
    #                pointList.append( float(strings[2])*1000 - MEAN_MOON_RADIUS ) # alt
    #            else:
                pointList.append(float(strings[2])) # alt
                #pointList.append(float(line)) # Add all numbers to the list
        i = i + 1
    f.close()

    #print pointList
    return pointList


