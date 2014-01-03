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

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tool for applying a rotation and position to a camera using mkspk and msopck
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------

# Creates the required msopck setup file if it does not already exist
def makeCkSetupFile(leapSecondFilePath, clockFilePath, frameFilePath, outputPath, tempFolder):

    # If the file already exists, delete it and rewrite it.
    if os.path.exists(outputPath):
        os.remove(outputPath)

    # The program can't handle long paths so we need to replace them with short symlinks
    leapFileName  = os.path.basename(leapSecondFilePath)
    clockFileName = os.path.basename(clockFilePath)
    frameFileName = os.path.basename(frameFilePath)
    leapSymPath   = os.path.join('/tmp', leapFileName)
    clockSymPath  = os.path.join('/tmp', clockFileName)
    frameSymPath  = os.path.join('/tmp', frameFileName)
    os.system('ln -s ' + leapSecondFilePath + ' ' + leapSymPath)
    os.system('ln -s ' + clockFilePath      + ' ' + clockSymPath)
    os.system('ln -s ' + frameFilePath      + ' ' + frameSymPath)

    f = open(outputPath, 'w')
    f.write("\\begindata\n")
    f.write("CK_TYPE              = 3\n") # 3 is usually the reccomended type
    f.write("INPUT_DATA_TYPE      = 'MATRICES'\n")
    f.write("INPUT_TIME_TYPE      = 'ET'\n")
    f.write("INSTRUMENT_ID        = -85000\n") # LRO
    f.write("REFERENCE_FRAME_NAME = 'J2000'\n")
    f.write("ANGULAR_RATE_PRESENT = 'NO'\n")
    f.write("PRODUCER_ID          = 'Lronac Pipeline'\n")
    f.write("LSK_FILE_NAME        = '" + leapSymPath + "'\n")
    f.write("SCLK_FILE_NAME       = '" + clockSymPath + "'\n")
    f.write("FRAMES_FILE_NAME     = '" + frameSymPath + "'\n")
    f.write("CK_SEGMENT_ID        = 'CK_MATRICES'\n")
    f.write("\\begintext\n")
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
        # Each line just contains the floating point value of the rotation matrix
        rotData.append(float(line)) 
    f.close()

    return rotData

#--------------------------------------------------------------------------------

def main():

# TODO: Rename this file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print "Started rotationCorrector.py"

    try:
        try:
            usage = "usage: rotationCorrector.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--input",  dest="inputPath",  help="Path to input .cub file to modify")

            parser.add_option("-s", "--spk", dest="spkPath", help="Path to write new SPK file to.")
            parser.add_option("-c", "--ck",  dest="ckPath",  help="Path to write new CK file to.")
            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output file.")
            parser.add_option("--transformPath",  dest="transformPath",  
                              help="Path to input file containing transform to apply")

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()

            if not options.inputPath: 
                parser.error("Need input path")

            if not options.outputPath: 
                parser.error("Need output path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder  = os.path.dirname(options.outputPath)
        inputBaseName = os.path.basename(options.inputPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_rotPosCorrectTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

        # File must already have had spiceinit called


        # Copy the input file to the output location (only the RE image is modified
        cmd = "cp " + options.inputPath + " " + options.outputPath
        print cmd
        os.system(cmd)


        # Retrieve a list of all the kernels needed by the input cube file
        kernelDict = IsisTools.getKernelsFromCube(options.outputPath, tempFolder)

        # Locate required kernels
        if not ('LeapSecond' in kernelDict):
            raise Exception('Error! Unable to find leap second file!')
        leapSecondFilePath = kernelDict['LeapSecond'][0] # Only deal with a single file

        if not ('SpacecraftClock' in kernelDict):
            raise Exception('Error! Unable to find clock kernel file!')
        clockFilePath = kernelDict['SpacecraftClock'][0] # Only deal with a single file

        if not ('Frame' in kernelDict):
            raise Exception('Error! Unable to find frame kernel file!')
        frameFilePath = kernelDict['Frame'][0] # Only deal with a single file

        # Convert the kernels into a space delimited string to pass as arguments
        kernelStringList = ""
        for k, v in kernelDict.iteritems(): # Iterate through dictionary entries
            # Add everything except the gigantic shape model (DEM)
            if k != 'ShapeModel':
                for i in v: # Iterate through type lists
                    kernelStringList = kernelStringList + ' ' + str(i)

        # Make sure the SPK and CK data paths do not already exist
        tempDataPrefix = os.path.join(tempFolder, "tempNavData")
        spkDataPath    = tempDataPrefix + "-spkData.txt"
        ckDataPath     = tempDataPrefix + "-ckData.txt"
        if os.path.exists(spkDataPath):
            os.remove(spkDataPath)
        if os.path.exists(ckDataPath):
            os.remove(ckDataPath)

        # Call lronac spice editor tool to generate modified text file
        cmd = ('spiceEditor --transformType 0 --transformFile ' + options.transformPath + 
                          ' --outputPrefix ' + tempDataPrefix + 
                          ' --kernels '      + kernelStringList + 
                          ' --sourceCube '   + options.inputPath)
        print cmd
        os.system(cmd)
        if not os.path.exists(spkDataPath):
            raise Exception('Error! Failed to create modified SPK data!')
        if not os.path.exists(ckDataPath):
            raise Exception('Error! Failed to create modified CK data!')

        # Write the config file needed for the mkspk function
        print 'Writing mkspk config file...'
        mkspkConfigPath = os.path.join(tempFolder, "spkConfig.txt")
        IsisTools.makeSpkSetupFile(leapSecondFilePath, mkspkConfigPath)

        # If the file already exists, delete it and rewrite it.
        if options.spkPath:
            tempSpkPath = options.spkPath
            print 'Storing modified SPK file ' + tempSpkPath
        else:
            tempSpkPath = os.path.join(tempFolder, "modifiedLrocSpk.bsp")
        if os.path.exists(tempSpkPath):
            os.remove(tempSpkPath)

        # Create new SPK file using modified data
        cmd = 'mkspk -setup ' + mkspkConfigPath + ' -input ' + spkDataPath + ' -output ' + tempSpkPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempSpkPath):
            raise Exception('Error! Failed to create modified SPK file!')

        # Write the config file needed for the msopck function
        print 'Writing msopck config file...'
        msopckConfigPath = os.path.join(tempFolder, "ckConfig.txt")
        makeCkSetupFile(leapSecondFilePath, clockFilePath, frameFilePath, msopckConfigPath, tempFolder)

        # If the file already exists, delete it and rewrite it.
        if options.ckPath:
            tempCkPath = options.ckPath
            print 'Storing modified CK file ' + tempCkPath
        else:
            tempCkPath = os.path.join(tempFolder, "modifiedLrocCk.bc")
        if os.path.exists(tempCkPath):
            os.remove(tempCkPath)

        # For unknown reasons there seems to be serious restrictions on directories this will actually write to.
        # We work around this by writing to a 'safe' path and copying the output to the desired location.
        reallyTempCkPath = os.path.join('/tmp/', inputBaseName + '_modifiedLrocCk.bc')  
        if os.path.exists(reallyTempCkPath):
            os.remove(reallyTempCkPath)

        # Create new CK file using modified data
        cmd = 'msopck ' + msopckConfigPath + ' ' + ckDataPath + '  ' + reallyTempCkPath
        print cmd
        os.system(cmd)
        if not os.path.exists(reallyTempCkPath):
            raise Exception('Error! Failed to create modified CK file!')

        # Copy the CK file to the actual desired location
        cmd = 'cp ' + reallyTempCkPath + ' ' + tempCkPath
        print cmd
        os.system(cmd)


        # Re-run spiceinit using the new SPK and CK file
        cmd = ("spiceinit attach=true from=" + options.outputPath + 
                           " spk=" + tempSpkPath + " ck=" + tempCkPath) # + " fk=../../lroFrameZero.tf"
        print cmd
        os.system(cmd)



        # Clean up temporary files
        if not options.keep:
            os.remove(tempTextPath)
            os.remove(rotationAnglePath)
            os.remove(tempIkPath)
            os.remove(options.spkPath)
            os.remove(options.ckPath)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
