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

import IsisTools, IrgIsisFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Applies the LROC offset from spacecraft position to an LROC cube's nav data
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main(argsIn):

    print "Started positionCorrectorV2.py"

    try:
        try:
            usage = "usage: positionCorrectorV2.py [--input <cube file>][--output <path>][--manual]\n  "
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
            (options, args) = parser.parse_args(argsIn)

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
        #print cmd
        os.system(cmd)

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


        # Call lronac spice editor tool to generate modified text file
        # - Call is silent unless there is an error
        cmd = ['spiceEditorV2', '--transformType', '1', '--transformFile', lronacOffsetPath, 
                                '--sourceCube', options.outputPath]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        outputText, err = p.communicate()       
        print outputText 

        # TODO: Tool should create the output file itself and check for error
#        if not os.path.exists(spkDataPath):
#            os.remove(options.outputPath)
#            print cmd
#            print outputText 
#            raise Exception('Failed to create modified SPK data!')

 
#        if (outputText.find('ERROR') >= 0):
#            print cmd
#            print outputText
#            raise Exception('Found error when calling spiceinit!') 

        # Clean up temporary files
        if not options.keep:
            os.remove(lronacOffsetPath)

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
