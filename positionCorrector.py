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
Applies the LROC offset from spacecraft position to an LROC cube's spice data
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

# Creates the required mkspk setup file if it does not already exist
def makeSpkSetupFile(leapSecondFilePath, outputPath):

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

#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    print "Started positionCorrector.py"

    try:
        try:
            usage = "usage: positionCorrector.py [--input <cube file>][--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-i", "--input", dest="inputPath",
                              help="Path to the input file.")
            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output file.")

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true",
                              dest="keep",
                              help="Do not delete the temporary files.")
            parser.add_option("--spk", dest="spkPath", help='Output SPK path (always keep)')
            (options, args) = parser.parse_args()

            if not options.inputPath:  parser.error("Need input path")
            if not options.outputPath: parser.error("Need output path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder  = os.path.dirname(options.outputPath)
        inputBaseName = os.path.basename(options.inputPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_posCorrectTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

        # Copy the input file to the output location
        cmd = "cp " + options.inputPath + " " + options.outputPath
        print cmd
        os.system(cmd)

        #TODO: Remove or make optional?
        # Call spiceinit on input file
        cmd = "spiceinit attach=true from=" + options.outputPath
        print cmd
        os.system(cmd)

        #DEBUG
        #campPath = os.path.join(outputFolder, "campt_orig.txt") 
        #cmd = "campt from=" + options.outputPath + " > "+campPath
        #print cmd
        #os.system(cmd)

        # Retrieve a list of all the kernels needed by the input cube file
        kernelDict = IsisTools.getKernelsFromCube(options.outputPath, tempFolder)

        # Find the leap second file
        if not ('LeapSecond' in kernelDict):
            print 'Error! Unable to find leap second file!'
            return 1
        else:
            leapSecondFilePath = kernelDict['LeapSecond'][0] # Only deal with a single file
            

        # Convert the kernels into a space delimited string to pass as arguments
        kernelStringList = ""
        for k, v in kernelDict.iteritems(): # Iterate through dictionary entries
            for i in v: # Iterate through type lists
                kernelStringList = kernelStringList + ' ' + str(i)

        # Determine if the input file is LE or RE
        lePos = options.inputPath.rfind('LE')
        rePos = options.inputPath.rfind('RE')
        if (lePos > rePos):
            sideCode = '0' # LE
        else: 
            sideCode = '1' # RE


        # Make sure the SPK data path does not already exist
        tempDataPrefix = os.path.join(tempFolder, "tempNavData")
        spkDataPath    = tempDataPrefix + "-spkData.txt"
        if os.path.exists(spkDataPath):
            os.remove(spkDataPath)

        # Call lronac spice editor tool to generate modified text file
        cmd = 'spiceEditor --offsetCode ' + sideCode + ' --outputPrefix ' + tempDataPrefix + ' --kernels ' + kernelStringList
        print cmd
        os.system(cmd)
        if not os.path.exists(spkDataPath):
            print 'Error! Failed to create modified SPK data!'
            return 1

        # Write the config file needed for the mkspk function
        print 'Writing mkspk config file...'
        mkspkConfigPath = os.path.join(tempFolder, "spkConfig.txt")
        makeSpkSetupFile(leapSecondFilePath, mkspkConfigPath)

        # If the file already exists, delete it and rewrite it.
        if options.spkPath:
          tempSpkPath = options.spkPath
          print 'Storing modified SPK file ' + tempSpkPath
        else:
          tempSpkPath = os.path.join(tempFolder, "modifiedLrocSpk.bsp")
        if os.path.exists(tempSpkPath):
            os.remove(tempSpkPath)

        #TODO: Need to make sure this is stored as an absolute path!!!!!!!!!!!!!

        # Create new SPK file using modified data
        cmd = 'mkspk -setup ' + mkspkConfigPath + ' -input ' + spkDataPath + ' -output ' + tempSpkPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempSpkPath):
            print 'Error! Failed to create modified SPK file!'
            return 1

        # Re-run spiceinit using the new SPK file
        cmd = "spiceinit attach=true from=" + options.outputPath + " spk=" + tempSpkPath
        print cmd
        os.system(cmd)


        #DEBUG
        #campPath = os.path.join(outputFolder, "campt_modified.txt") 
        #cmd = "campt from=" + options.outputPath + " > "+campPath
        #print cmd
        #os.system(cmd)

        # Clean up temporary files
        if not options.keep:
            #os.remove(tempTextPath)
            os.remove(spkDataPath)
            os.remove(mkspkConfigPath)
            #if not options.spkPath: # Really need to keep the SPK file around
            #    os.remove(tempSpkPath)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
