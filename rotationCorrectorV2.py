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

import os, glob, optparse, re, shutil, subprocess, string, time

import IsisTools, IrgIsisFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tool for applying a rotation and position to a cube file
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

def main(argsIn):

    print "Started rotationCorrectorV2.py"

    try:
        try:
            usage = "usage: rotationCorrectorV2.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.set_defaults(keep=False)
            parser.set_defaults(pcAlignTrans=False)

            parser.add_option("--input",  dest="inputPath",  help="Path to input .cub file to modify")

            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output file.")
            parser.add_option("--transformPath",  dest="transformPath",  
                              help="Path to input file containing transform to apply")

            # This flag means that the rotation also affects the position!
            parser.add_option("--pcAlignTrans", action="store_true", dest="pcAlignTrans",
                              help="The transform is from pc_align.")

#            # The default working directory path is kind of ugly...
#            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
#            parser.add_option("--keep", action="store_true", dest="keep",
#                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args(argsIn)

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
#        tempFolder    = outputFolder + '/' + inputBaseName + '_rotPosCorrectTemp/'
#        if (options.workDir):
#            tempFolder = options.workDir
#        if not os.path.exists(tempFolder):
#            os.mkdir(tempFolder) 

        # File must already have had spiceinit called


        # Copy the input file to the output location (only the RE image is modified
        cmd = "cp " + options.inputPath + " " + options.outputPath
        #print cmd
        os.system(cmd)

        transformType = 0 # GLOBAL, rotation does not affect position, from lronacAngleSolver
        if options.pcAlignTrans:
            transformType = 2 #PC_ALIGN, rotation affects position

        # Call lronac spice editor tool to generate modified text file
        # - Call is silent unless an error occurs
        cmd = ['spiceEditorV2', '--transformType', str(transformType),
                              '--transformFile', options.transformPath,
                              '--sourceCube',    options.outputPath]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        outputText, err = p.communicate()
        print outputText
#        if not os.path.exists(spkDataPath):
#            os.remove(options.outputPath)
#            print cmd
#            print outputText
#            raise Exception('Error! Failed to create modified SPK data!')

#        if (outputText.find('ERROR') >= 0):
#            print cmd
#            print outputText
#            raise Exception('Found error when calling spiceinit!')

#        # Clean up temporary files
#        if not options.keep:
#            os.remove(tempTextPath)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
