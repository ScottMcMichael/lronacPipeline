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

        # Retrieve a list of all the kernels needed by the input cube file
        kernelDict = IsisTools.getKernelsFromCube(options.outputPath, tempFolder)

        # Find the leap second file
        if not ('LeapSecond' in kernelDict):
            raise Exception('Error! Unable to find leap second file!')
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

        # Set the offsets from LRO spacecraft to the LRONAC cameras in meters
        if (lePos > rePos):
            lronacOffset = [1.3462, 0.8890, -0.1778] # LE
        else: 
            lronacOffset = [1.0160, 0.8890, -0.1778] # RE

        # Write out a file containing the LRONAC camera offset
        lronacOffsetPath = os.path.join(tempFolder, "lronacCameraOffset.csv")
        f = open(lronacOffsetPath, 'w')
        f.write('1 0 0 ' + str(lronacOffset[0]) + '\n')
        f.write('0 1 0 ' + str(lronacOffset[1]) + '\n')
        f.write('0 0 1 ' + str(lronacOffset[2]) + '\n')
        f.write('0 0 0 1')
        f.close()


        # Make sure the SPK data path does not already exist
        tempDataPrefix = os.path.join(tempFolder, "tempNavData")
        spkDataPath    = tempDataPrefix + "-spkData.txt"
        if os.path.exists(spkDataPath):
            os.remove(spkDataPath)

        # Call lronac spice editor tool to generate modified text file
        cmd = 'spiceEditor --transformType 1 --transformFile ' + lronacOffsetPath + ' --outputPrefix ' + tempDataPrefix + ' --kernels ' + kernelStringList +  ' --sourceCube ' + options.inputPath
        print cmd
        os.system(cmd)
        if not os.path.exists(spkDataPath):
            raise Exception('Failed to create modified SPK data!')

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

        #TODO: Need to make sure this is stored as an absolute path!!!!!!!!!!!!!

        # Create new SPK file using modified data
        cmd = 'mkspk -setup ' + mkspkConfigPath + ' -input ' + spkDataPath + ' -output ' + tempSpkPath
        print cmd
        os.system(cmd)
        if not os.path.exists(tempSpkPath):
            raise Exception('Failed to create modified SPK file!')

        # Re-run spiceinit using the new SPK file
        cmd = "spiceinit attach=true from=" + options.outputPath + " spk=" + tempSpkPath
        print cmd
        os.system(cmd)

        # Clean up temporary files
        if not options.keep:
            #os.remove(tempTextPath)
            os.remove(spkDataPath)
            os.remove(mkspkConfigPath)

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
