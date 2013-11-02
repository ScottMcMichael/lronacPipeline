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

# Parses the output from head [cube path]
def parseHeadOutput(textPath):

    kernelList = []

    isisDataFolder = os.environ['ISIS3DATA']

    # Search each line in the folder for a required kernel file
    dataFile = open(textPath, 'r')
    lastLine = ''
    for line in dataFile:
        # Append leftovers from last line and clear left/right whitespace
        workingLine = lastLine + line.strip()
#        print 'workingLine =' + workingLine
        if (workingLine.find('/kernels/') >= 0):
            m = re.search('\$[a-zA-Z0-9/._\-]*', workingLine)
            if m: # Path found
                if (m.group(0)[-1] == '-'): # This means ISIS has done a weird truncation to the next line
                    lastLine = m.group(0)[:-1] # Strip trailing - and append next line to it
                else: # Valid match
#                    print 'found kernel path ' + m.group(0)
                    kernelPath = os.path.join(isisDataFolder, m.group(0)[1:])
                    
                    if not os.path.exists(kernelPath): # Make sure the kernel file exists
                        print 'Error! Specified kernel file ' + kernelPath + ' does not exist!'
                        return [] # Fail if we get a miss
                    kernelList.append(kernelPath)
                    lastLine = ''
            else: 
               print 'Failed to find kernel in line: ' + line
               
    # Return the list of kernels
    return kernelList

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
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

        # Copy the input file to the output location
        cmd = "cp " + options.inputPath + " " + options.outputPath
        print cmd
        os.system(cmd)

        # Call spiceinit on input file
        cmd = "spiceinit from=" + options.outputPath
        print cmd
        os.system(cmd)

        #DEBUG
        #campPath = os.path.join(outputFolder, "campt_orig.txt") 
        #cmd = "campt from=" + options.outputPath + " > "+campPath
        #print cmd
        #os.system(cmd)

        # Call head -120 on file
        tempTextPath = os.path.join(tempFolder, "headOutput.txt")
        cmd = "head -120 "+options.outputPath+" > "+tempTextPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempTextPath):
            print 'Error! Failed to extract cube kernel data!'
            return 1

        # Parse output
        print 'Parsing cube metadata for kernel locations...'
        kernelList = parseHeadOutput(tempTextPath)
        if not kernelList:
            print 'Error! Unable to find any kernel files in ' + tempTextPath
            return 1
#        print 'Found ' + str(len(kernelList)) + ' kernel files'

        # Find the leap second file
        for k in kernelList:
            if (k.find('/kernels/lsk/naif') >= 0): # This should week out all other kernels
                leapSecondFilePath = k
        if not leapSecondFilePath:
            print 'Error! Unable to find leap second file!'
            return 1


        # Convert the kernels into a space delimited string to pass as arguments
        kernelStringList = ""
        for i in kernelList:
          kernelStringList = kernelStringList + ' ' + str(i)
#        print kernelStringList

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
        cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/spiceEditor --offsetCode ' + sideCode + ' --outputPrefix ' + tempDataPrefix + ' --kernels ' + kernelStringList
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

        # Create new SPK file using modified data
        cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/mkspk -setup ' + mkspkConfigPath + ' -input ' + spkDataPath + ' -output ' + tempSpkPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempSpkPath):
            print 'Error! Failed to create modified SPK file!'
            return 1

        # Re-run spiceinit using the new SPK file
        cmd = "spiceinit from=" + options.outputPath + " spk=" + tempSpkPath
        print cmd
        os.system(cmd)


        #DEBUG
        #campPath = os.path.join(outputFolder, "campt_modified.txt") 
        #cmd = "campt from=" + options.outputPath + " > "+campPath
        #print cmd
        #os.system(cmd)

        # Clean up temporary files
        if not options.keep:
            os.remove(tempTextPath)
            os.remove(spkDataPath)
            os.remove(mkspkConfigPath)
            if not options.spkPath:
                os.remove(tempSpkPath)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
