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
Solves for the best rotation angle between the LROC cameras and updates the rotation in the FK file.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

# Makes a modified copy of an FK file using the new rotation
def modifyFrameFile(frameFilePath, outputPath, newRotData):

    # If the file already exists, delete it and rewrite it.
    if os.path.exists(outputPath):
        os.remove(outputPath)
    if not os.path.exists(frameFilePath):
        print 'Error, file ' + frameFilePath + ' not found!'
        return false

    # Find and modify this line: TKFRAME_-85610_ANGLES    = ( 0.0, 0.0, 0.0 )


    print 'reading file ' + frameFilePath
    i = open(frameFilePath, 'r')
    f = open(outputPath, 'w')
    for line in i:
        #print 'input --> ' + line
        if (line.find('TKFRAME_-85610_ANGLES') >= 0): # Line where rotation amounts are specified
            f.write('      TKFRAME_-85610_MATRIX    = ( ') # Replace with matrix
            for d in newRotData:
                f.write('\n' + '                                  ' + str(d))
            f.write(' )\n')
        elif (line.find('TKFRAME_-85610_SPEC') >= 0): # Line where data type is specified
            newString = "      TKFRAME_-85610_SPEC    =  'Matrix'\n" # Replace data type
            f.write(newString)
        elif (line.find('TKFRAME_-85610_AXES') >= 0): # Line where rotation axes are specified
            f.write('') # Skip the line
        elif (line.find('TKFRAME_-85610_UNITS') >= 0): # Line where angular units are specified
            f.write('') # Skip the line
        else: # Normal line, just copy it
            #print 'normal -->' + line
            f.write(line)
    i.close()
    f.close()

# Parses the output from head [cube path]
def readRotationFile(rotFilePath):

    if not os.path.exists(rotFilePath):
        print 'Error: file ' + rotFilePath + ' does not exist!'
        return False

    # Read through file one line at a time
    f = open(rotFilePath)
    rotData = []
    for line in f:
        rotData.append(float(line)) # Each line just contains the floating point value of the rotation matrix
    f.close()

    return rotData

# TODO: Replace this with call to IsisTools!
# Parses the output from head [cube path]
def parseHeadOutput(textPath):

    isisDataFolder = os.environ['ISIS3DATA']

    # Search each line in the folder for a required kernel file
    dataFile = open(textPath, 'r')
    lastLine = ''
    for line in dataFile:
        # Append leftovers from last line and clear left/right whitespace
        workingLine = lastLine + line.strip()
        #print 'workingLine =' + workingLine
        if (workingLine.find('/kernels/fk/lro_frames_') >= 0): # This should week out all other kernels
            m = re.search('\$[a-zA-Z0-9/._\-]*', workingLine)
            if m: # Path found
                if (m.group(0)[-1] == '-'): # This means ISIS has done a weird truncation to the next line
                    lastLine = m.group(0)[:-1] # Strip trailing - and append next line to it
                else: # Valid match
                    #print 'found kernel path ' + m.group(0)
                    kernelPath = os.path.join(isisDataFolder, m.group(0)[1:])

                    if not os.path.exists(kernelPath): # Make sure the kernel file exists
                        print 'Error! Specified kernel file ' + kernelPath + ' does not exist!'
                        return [] # Fail if we get a miss
                    
                    return kernelPath # We are only looking for this one kernel file
            else: 
               print 'Failed to find kernel in line: ' + line
               
    # Failed to find the frame kernel!
    return []

#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    print "Started lronacCameraRotationCorrector.py"

    try:
        try:
            usage = "usage: rotationCorrector.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--left",  dest="leftPath",  help="Path to LE .cub file")
            parser.add_option("--right", dest="rightPath", help="Path to RE .cub file")

            parser.add_option("-g", "--gdcLogPath", dest="gdcLogPath",
                              help="Optional path to save computed GDC points to")

            parser.add_option("-s", "--spk", dest="spkPath",
                              help="Path to optional specified SPK (position) file to use.")
            parser.add_option("-c", "--ck", dest="ckPath",
                              help="Path to optional specified CK (orientation) file to use.")


            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output (RE) file.")
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.outputPath: 
                parser.error("Need output path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder = os.path.dirname(options.outputPath)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_rotCorrectTemp/'
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

        # File must already have had spiceinit called


        # Copy the input file to the output location (only the RE image is modified
        cmd = "cp " + options.rightPath + " " + options.outputPath
        print cmd
        os.system(cmd)

        # Call head -120 on file
        tempTextPath = os.path.join(tempFolder, "headOutput.txt")
        cmd = "head -120 "+options.outputPath+" > "+tempTextPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempTextPath):
            print 'Error! Failed to extract cube kernel data!'
            return 1

        # Parse output looking for the IK frame file
        print 'Looking for source frame file...'
        inputFramePath = parseHeadOutput(tempTextPath)
        if not inputFramePath:
            print 'Error! Unable to find any IK kernel file in ' + tempTextPath
            return 1

        # Make sure the output path does not already exist
        rotationAnglePath = os.path.join(tempFolder, "solvedRotationAngles.txt")
        if os.path.exists(rotationAnglePath):
            os.remove(rotationAnglePath)
        
        # Call lronacSpkParser to generate modified text file
        gdcText = ''
        if options.gdcLogPath: # Handle GDC point logging option
            gdcText = ' --gdcPointsOutPath ' + options.gdcLogPath
        cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath '+ \
                  rotationAnglePath + gdcText + ' ' + options.leftPath + ' ' + options.rightPath
        print cmd
        os.system(cmd)
        if not os.path.exists(rotationAnglePath):
            print 'Error! Failed to solve for rotation angles!'
            return 1
        
        # Read the rotation angles
        newRotation = readRotationFile(rotationAnglePath)

        # Generate a modified frame file
        tempIkPath = os.path.join(tempFolder, "angleCorrectedIkFile.tf")
        modifyFrameFile(inputFramePath, tempIkPath, newRotation)
        
        # Re-run spiceinit (on the copied RE file) using the new frame file
        cmd = "spiceinit from=" + options.outputPath + " fk=" + tempIkPath;
        if (options.spkPath): # Add forced SPK path if needed
            cmd = cmd + " spk=" + options.spkPath
        if (options.ckPath): # Add forced CK path if needed
            cmd = cmd + " ck=" + options.ckPath
        print cmd
        os.system(cmd)

        # Clean up temporary files
        if not options.keep:
            os.remove(tempTextPath)
            os.remove(rotationAnglePath)
            os.remove(tempIkPath)
            os.remove(options.spkPath)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
