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

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging, threading

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two LRONAC pairs using SBA and LOLA for increased accuracy.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------


# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('lronac2refinedMosaics.py')
    IsisTools.checkIfToolExists('makeDemAndCompare.py')
    return True

#--------------------------------------------------------------------------------

def main():

    print '#################################################################################'
    print "Running lronac2dem.py"

    try:
        try:
            usage = "usage: lronac2dem.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)

            parser.set_defaults(keep=False)

            inputGroup = optparse.OptionGroup(parser, 'Input Paths')
            inputGroup.add_option("--left",  dest="leftPath",  help="Path to LE .IMG file")
            inputGroup.add_option("--right", dest="rightPath", help="Path to RE .IMG file")            
            inputGroup.add_option("--stereo-left",  dest="stereoLeft", 
                                  help="Path to LE .IMG file with overlapping view of --left file")
            inputGroup.add_option("--stereo-right", dest="stereoRight", 
                                  help="Path to RE .IMG file with overlapping view of --right file")

            inputGroup.add_option("--lola",    dest="lolaPath", help="Path to LOLA DEM")
            inputGroup.add_option("--asu",     dest="asuPath",  help="Path to ASU DEM")
            
            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--prefix",  dest="prefix",   help="Output prefix.")

            parser.add_option("--crop",  dest="cropAmount", 
                              help="Crops the output image to reduce processing time.")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.stereoLeft: 
                parser.error("Need stereo left input path")
            if not options.stereoRight: 
                parser.error("Need stereo right input path")
            if not options.lolaPath: 
                parser.error("Need LOLA data path")
            if not options.prefix: 
                parser.error("Need output prefix")

        except optparse.OptionError, msg:
            raise Usage(msg)

        startTime = time.time()

        # Set up the output folders
        outputFolder  = os.path.dirname(options.prefix)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder) 
        hadToCreateTempFolder = not os.path.exists(tempFolder)
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 


        # Set up logging
        logPath = options.prefix + '-Log.txt'
        logging.basicConfig(filename=logPath,level=logging.INFO)

        keepString = ''
        if options.keep:
            keepString = ' --keep'
  

        # Call lronac2refinedMosaics.py
        refineTemp = os.path.join(tempFolder, 'refinement')
        cmd = ('lronac2refinedMosaics.py --left '          + options.leftPath + 
                                       ' --right '         + options.rightPath + 
                                       ' --stereo-left '   + options.stereoLeft + 
                                       ' --stereo-right '  + options.stereoRight + 
                                       ' --lola '          + options.lolaPath + 
                                       ' --workDir '       + refineTemp + 
                                       ' --output-folder ' + tempFolder + 
                                       ' --log-path '      + logPath + 
                                       keepString)
        print cmd
        os.system(cmd)
        
        # Copy the pc_align log file to the output folder
        pcAlignLogPath = os.path.join(tempFolder, 'pcAlignLog.txt')
        shutil.copyfile(pcAlignLogPath, os.path.join(outputFolder, 'pcAlignLog.txt'))


        # These are the output mosaic paths
        filename         = os.path.splitext(options.leftPath)[0] + '.correctedMosaic.cub'
        mainMosaicPath   = os.path.join(tempFolder, os.path.basename(filename))
        filename         = os.path.splitext(options.stereoLeft)[0] + '.correctedMosaic.cub'
        stereoMosaicPath = os.path.join(tempFolder, os.path.basename(filename))

        if (not os.path.exists(mainMosaicPath) or not os.path.exists(stereoMosaicPath)):
            raise Exception('lronac2refinedMosaics failed to produce mosaics!')

        # Call makeDemAndCompare.py
        cmd = ('makeDemAndCompare.py --left '     + mainMosaicPath + 
                                   ' --right '    + stereoMosaicPath + 
                                   ' --lola '     + options.lolaPath + 
                                   ' --asu '      + options.asuPath + 
                                   ' --workDir '  + tempFolder + 
                                   ' --prefix '   + options.prefix + 
                                   ' --log-path ' + logPath + 
                                   keepString)
        if options.cropAmount:
            cmd = cmd + ' --crop ' + str(options.cropAmount)
        print cmd
        os.system(cmd)


        if not options.keep:
            print 'Deleting temporary files'
            IsisTools.removeIfExists(mainMosaicPath)
            IsisTools.removeIfExists(stereoMosaicPath)


        endTime = time.time()

        logging.info('lronac2dem finished in %f seconds', endTime - startTime)
        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
