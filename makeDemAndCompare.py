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

import sys, os, glob, optparse, re, subprocess, string, time, logging, threading

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two input mosaics and optionally compares to ASU and LOLA DEMs.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------


# Generate a comparison of a DEM to the LOLA point cloud
def compareDemToLola(lolaPath, demPath, outputPath, csvPath, force):

    if force or not os.path.exists(outputPath):
        cmd = ('lola_compare --absolute 1 --limit-hist=2 -l "' + lolaPath + 
                             '" -d ' + demPath + ' -o ' + outputPath + ' -c ' + csvPath)
        print cmd
        os.system(cmd)
        
    return True


# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('lola_compare')
    IsisTools.checkIfToolExists('parallel_stereo')
    IsisTools.checkIfToolExists('point2dem')
    IsisTools.checkIfToolExists('hillshade')
    IsisTools.checkIfToolExists('crop')

    return True

#--------------------------------------------------------------------------------

def main(argsIn):

    print '#################################################################################'
    print "Running makeDemAndCompare.py"

    try:
        try:
            usage = "usage: makeDemAndCompare.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
            inputGroup = optparse.OptionGroup(parser, 'Input Paths')
            inputGroup.add_option("--left",  dest="leftPath",  help="Path to left  cube file")
            inputGroup.add_option("--right", dest="rightPath", help="Path to right cube file")            


            inputGroup.add_option("--lola",    dest="lolaPath", help="Path to LOLA DEM")
            inputGroup.add_option("--asu",     dest="asuPath",  help="Path to ASU DEM")
            
            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--prefix",  dest="prefix",   help="Output prefix.")

            parser.add_option("--log-path",  dest="logPath",        
                              help="Where to write the output log file.")

            parser.add_option("--crop",  dest="cropAmount", 
                              help="Crops the output image to reduce processing time.")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args(argsIn)

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.prefix: 
                parser.error("Need output prefix")            

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        # Make sure we have all the functions we need
        functionStartupCheck()

        # Set this to true to force steps after it to activate
        carry = False

        # Verify input files are present
        if not os.path.exists(options.leftPath):
            raise Exception('Input file ' + options.leftPath + ' not found!')
        if not os.path.exists(options.rightPath):
            raise Exception('Input file ' + options.rightPath + ' not found!')

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
        if not options.logPath:
            options.logPath = options.prefix + '-Log.txt'
        logging.basicConfig(filename=options.logPath,level=logging.INFO)

        # If specified, crop the inputs that will be passed into the stereo function to reduce processing time
        mainMosaicCroppedPath   = os.path.join(tempFolder, 'mainMosaicCropped.cub')
        stereoMosaicCroppedPath = os.path.join(tempFolder, 'stereoMosaicCropped.cub')
        if options.cropAmount and (options.cropAmount > 0):
            if (not os.path.exists(mainMosaicCroppedPath)) or carry:
                cmd = ('crop from= ' + options.leftPath   + ' to= ' + mainMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))
                print cmd
                os.system(cmd)
            if (not os.path.exists(stereoMosaicCroppedPath) or carry):
                cmd = ('crop from= ' + options.rightPath + ' to= ' + stereoMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))
                print cmd
                os.system(cmd)
            options.leftPath  = mainMosaicCroppedPath
            options.rightPath = stereoMosaicCroppedPath

        print '\n-------------------------------------------------------------------------\n'

        # Call stereo to generate a point cloud from the two images
        # - This step takes a really long time.

        stereoOutputPrefix = os.path.join(tempFolder, 'stereoWorkDir/stereo')
        pointCloudPath     = stereoOutputPrefix + '-PC.tif'
        if (not os.path.exists(pointCloudPath)) or carry:
            cmd = ('parallel_stereo --corr-timeout 400 --alignment affineepipolar --subpixel-mode 1' +
                                  ' --disable-fill-holes ' +  options.leftPath + ' ' + options.rightPath + 
                                  ' ' + stereoOutputPrefix + ' --processes 8 --threads-multiprocess 4' +
                                  ' --threads-singleprocess 32 --compute-error-vector')
            print cmd
            os.system(cmd)
            #--nodes-list PBS_NODEFILE --processes 4 --threads-multiprocess 16 --threads-singleprocess 32
        else:
            print 'Stereo file ' + pointCloudPath + ' already exists, skipping stereo step.'

        stereoTime = time.time()
        logging.info('Stereo finished in %f seconds', stereoTime - startTime)

        # Compute percentage of good pixels
        percentGood = IsisTools.getStereoGoodPixelPercentage(stereoOutputPrefix)
        print 'Stereo completed with good pixel percentage: ' + str(percentGood)
        logging.info('Final stereo completed with good pixel percentage: %s', str(percentGood))


        # Find out the center latitude of the mosaic
        centerLat = IsisTools.getCubeCenterLatitude(options.leftPath, tempFolder)

        # Generate a DEM
        demPrefix = os.path.join(outputFolder, 'p2d')
        demPath   = demPrefix + '-DEM.tif'
        if (not os.path.exists(demPath)) or carry:
            cmd = ('point2dem --errorimage -o ' + demPrefix + ' ' + pointCloudPath + 
                            ' -r moon --tr 1 --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + 
                            ' +lat_0=0 +a=1737400 +b=1737400 +units=m" --nodata -32767')
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + pointCloudPath + ' already exists, skipping point2dem step.'

        # Create a hillshade image to visualize the output
        hillshadePath = os.path.join(outputFolder, 'outputHillshade.tif')
        if (not os.path.exists(hillshadePath)) or carry:
            cmd = 'hillshade ' + demPath + ' -o ' + hillshadePath
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + hillshadePath + ' already exists, skipping hillshade step.'


        # Call script to compare LOLA data with the DEM
        if options.lolaPath or carry:
            lolaDiffStatsPath  = os.path.join(outputFolder, 'LOLA_diff_stats.txt')
            lolaDiffPointsPath = os.path.join(outputFolder, 'LOLA_diff_points.csv')
            compareDemToLola(options.lolaPath, demPath, lolaDiffStatsPath, lolaDiffPointsPath, carry)

        # Call script to compare LOLA data with the ASU DEM
        if options.asuPath or carry:
            lolaAsuDiffStatsPath  = os.path.join(outputFolder, 'ASU_LOLA_diff_stats.txt')
            lolaAsuDiffPointsPath = os.path.join(outputFolder, 'ASU_LOLA_diff_points.csv')
            compareDemToLola(options.lolaPath, options.asuPath, lolaAsuDiffStatsPath, lolaAsuDiffPointsPath, carry)

        statsTime = time.time()
        logging.info('Final results finished in %f seconds', statsTime - stereoTime)

        # Clean up temporary files
        if not options.keep:
            print 'Removing temporary files'
            IsisTools.removeIfExists(mainMosaicCroppedPath)
            IsisTools.removeIfExists(stereoMosaicCroppedPath)
            IsisTools.removeIntermediateStereoFiles(stereoOutputPrefix)
            #if (hadToCreateTempFolder): Not done since stereo output needs to be retained
            #    IsisTools.removeFolderIfExists(tempFolder)


        endTime = time.time()

        logging.info('makeDemAndCompare.py finished in %f seconds', endTime - startTime)
        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
