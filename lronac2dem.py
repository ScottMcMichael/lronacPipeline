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

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two LRONAC pairs, trying to use LOLA data for increased accuracy.
Uses the StereoCalibrationProcess function to align the input images.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

# Generates a KML file to describe a set of GDC points on the moon
def generateKmlFromGdcPoints(inputFolder, outputFolder, filename, pointSkip, color, size, forceOperation):

    # TODO: Clean up the way this works!!!
    # Determine input and output paths
    inputPath      = os.path.join(inputFolder,  'SBA_check-outputGdcPoints.csv')
    if not os.path.exists(inputPath):
        inputPath = os.path.join(inputFolder,  'out-initialGdcPoints.csv') # Large GDC point set
    if not os.path.exists(inputPath):
        inputPath = os.path.join(inputFolder,  'dem-trans_reference.csv') # pc_align output


    outputFilename = os.path.splitext(filename)[0] + '.kml'
    outputPath     = os.path.join(outputFolder, outputFilename)
    kmlName        = os.path.splitext(filename)[0]


    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputPath)):
        print 'File ' + outputPath + ' already exists, skipping kml generation.'
        return True

    # Generate the new file
    cmd = ('calibrationReport.py --input ' + inputPath + ' --output ' + outputPath + 
          ' --name ' + kmlName +  ' --skip ' + str(pointSkip) + ' --color ' + color + ' --size ' + size)
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    # - Don't throw an exception here since this is only a debug step
    if not os.path.exists(outputPath):
        print ('=======> WARNING: Kml point generation failed to create output file ' + outputPath + 
               ' from input file ' + inputPath)
        return False

    return True


# Calls the ISIS noproj function on a cube
def noprojCubePair(inputCube, outputCube, matchCube, pvlPath, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCube)):
        print 'File ' + outputCube + ' already exists, skipping noproj.'
        return True

    # Generate the new file
    cmd = ('noproj from= '  + inputCube + ' to= '    + outputCube + 
                 ' match= ' + matchCube + ' specs= ' + pvlPath)
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCube):
        raise Exception('noproj failed to create output file ' + outputCube + 
                        ' from input file ' + inputCube)

    return True


# Creates a mosaic from two noproj'd input cubes.
def createMosaic(leftCube, rightCube, outputCube, workDir, forceOperation):


    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCube)):
        print 'File ' + outputCube + ' already exists, skipping mosaic creation.'
        return True


    # Call lronacjitreg to determine any remaining offset between the two input files
    jitRegOutputPath = os.path.join(workDir, 'jitregResults.txt')
    cmd = ('lronacjitreg --correlator-type 2 --kernel 15 15 --output-log ' + jitRegOutputPath + 
                        ' ' + leftCube + ' ' + rightCube)
    print cmd
    os.system(cmd)

    # Read in the output from lronacjitreg
    jitRegOffsets = IsisTools.readJitregFile(jitRegOutputPath)
    logging.info('For cubes %s and %s', leftCube, rightCube)
    logging.info('- jitreg offsets = %s', str(jitRegOffsets))

    # Set intermediate mosaic file path and start on mosaic
    mosaicCube = os.path.join(workDir, 'mosaic.cub')
    cmd = 'cp ' + leftCube + ' ' + mosaicCube
    print cmd
    os.system(cmd)
    
    # TODO: Try out more advanced ISIS mosaic merging functions!
    # Create the mosaic, applying offsets from jitreg (converting into handmos conventions)
    cmd = ('handmos from= ' + rightCube + ' mosaic= ' + mosaicCube + 
                  ' outsample= ' + str(int(round(1 - jitRegOffsets[0]))) + 
                  ' outline= ' + str(int(round(1 - jitRegOffsets[1]))) + 
                  ' matchbandbin=FALSE priority=ontop')
    print cmd
    os.system(cmd)

    # Call cubenorm to improve the mosaic appearance
    cmd = 'cubenorm from= ' + mosaicCube + ' to= ' + outputCube
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCube):
        raise Exception('Failed to create mosaic file ' + outputCube + 
                        ' from input files ' + leftCube + ' and ' + rightCube)

    return True


# Generate a comparison of a DEM to the LOLA point cloud
def compareDemToLola(lolaPath, demPath, outputPath, force):

    if force or not os.path.exists(outputPath):
        cmd = ('lola_compare --absolute 1 --limit-hist=2 -l "' + lolaPath + 
                             '" -d ' + demPath + ' -o ' + outputPath)
        print cmd
        os.system(cmd)
        
    return True

#--------------------------------------------------------------------------------

def main():

    print '#################################################################################'
    print "Running lronac2dem.py"

    try:
        try:
            usage = "usage: lronac2dem.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
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
            #parser.add_option("--keep", action="store_true", dest="keep",
            #                  help="Do not delete the temporary files.")
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

        print "Beginning processing....."

        startTime = time.time()
        
        # Set this to true to force steps after it to activate
        carry = False

        # Verify input files are present
        if not os.path.exists(options.leftPath):
            raise Exception('Input file ' + options.leftPath + ' not found!')
        if not os.path.exists(options.rightPath):
            raise Exception('Input file ' + options.rightPath + ' not found!')
        if not os.path.exists(options.stereoLeft):
            raise Exception('Input file ' + options.stereoLeft + ' not found!')
        if not os.path.exists(options.stereoRight):
            raise Exception('Input file ' + options.stereoRight + ' not found!')
        if not os.path.exists(options.lolaPath):
            raise Exception('Input file ' + options.lolaPath + ' not found!')

        # Set up the output folders
        outputFolder  = os.path.dirname(options.prefix)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder) 
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 
        
        # Set up logging
        logPath = options.prefix + '-Log.txt'
        logging.basicConfig(filename=logPath,level=logging.INFO)


        # Set up output paths for the stereo calibration call
        leftCorrectedPath        = os.path.join(tempFolder, 'leftFinalCorrected.cub')
        rightCorrectedPath       = os.path.join(tempFolder, 'rightFinalCorrected.cub')
        leftStereoCorrectedPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.cub')
        rightStereoCorrectedPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.cub')
        correct = os.path.join(tempFolder, 'rightStereoFinalCorrected.cub')

        # Generate a kml plot of the input LOLA data
        lolaKmlPath = os.path.join(tempFolder, 'lolaRdrPoints.kml')
        if not os.path.exists(lolaKmlPath):
            cmd = ('calibrationReport.py --input '  + options.lolaPath + 
                                       ' --output ' + lolaKmlPath + 
                                       ' --name '   + 'lolaRdrPoints' + 
                                       ' --skip '   + str(100) + ' --color ' + 'blue', + ' --size small')
            print cmd
            os.system(cmd)
                       
        carry = False
                       
        # Correct all four input images at once
        caughtException = False
        try:
            if (not os.path.exists(leftCorrectedPath)) or carry:
                print '\n=============================================================================\n'
                cmd = ('stereoDoubleCalibrationProcess.py --left '  + options.leftPath + 
                                                        ' --right ' +  options.rightPath + 
                                                        ' --stereo-left '  + options.stereoLeft + 
                                                        ' --stereo-right ' + options.stereoRight + 
                                                        ' --keep --lola '  + options.lolaPath + 
                                                        ' --outputL '  + leftCorrectedPath + 
                                                        ' --outputR '  + rightCorrectedPath + 
                                                        ' --outputSL ' + leftStereoCorrectedPath + 
                                                        ' --outputSR ' + rightStereoCorrectedPath + 
                                                        ' --workDir '  + options.workDir + 
                                                        ' --log-path ' + logPath)
                print cmd
                os.system(cmd)
                print '\n============================================================================\n'
                
                # Convert GDC output files into KML plots 
                # - This is just to help with debugging
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'initialGdcCheck'),            tempFolder, 'pairGdcCheckInitial.csv',           1,    'blue', 'normal', carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'posCorrectGdcCheck'),         tempFolder, 'pairGdcCheckPos.csv',               1,    'green', 'normal',  carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'posCorrectStereoGdcCheck'),   tempFolder, 'pairGdcStereoCheckPos.csv',         1,    'green', 'normal',  carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'stereoGlobalAdjustGdcCheck'), tempFolder, 'pairGdcCheckGlobalAdjustStero.csv', 1,    'blue', 'normal', carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'pcAlignStereoGdcCheck'),      tempFolder, 'pairGdcCheckPcAlign.csv',           1,    'red', 'normal',  carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'finalGdcCheck'),              tempFolder, 'pairGdcCheckFinal.csv',             1,    'white', 'normal', carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'finalStereoGdcCheck'),        tempFolder, 'pairGdcCheckFinalStereo.csv',       1,    'white', 'normal',  carry)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'gdcPointsLargeComp'),         tempFolder, 'inputGdcPoints.csv',                1000, 'blue', 'tiny', True)
                #generateKmlFromGdcPoints(os.path.join(tempFolder, 'initialGdcCheck'),            tempFolder, 'dem-trans_source.csv',              1000, 'blue', 'normal', , False)
                generateKmlFromGdcPoints(os.path.join(tempFolder, 'pcAlignOutput'),              tempFolder, 'dem-trans_reference.csv',           1000, 'red', 'tiny',  True)

                print 'Finished generating KML plots'
                
        except: 
            caughtException = True
            print 'Caught an exception!'

        # Delay check for left path to allow debug KML files to be generated
        if caughtException or not os.path.exists(leftCorrectedPath):
            raise Exception('Failed to run stereo calibration process!')

        print '\n-------------------------------------------------------------------------\n'

        # Generate a PVL file that we need for noproj
        pvlPath   = os.path.join(tempFolder, 'noprojInstruments_fullRes.pvl')
        isHalfRes = False # TODO: Check to see if this is true!
        if not os.path.exists(pvlPath):
            print 'Writing PVL'
            IsisTools.writeLronacPvlFile(pvlPath, isHalfRes)


        # Noproj the corrected data
        print 'Starting noproj calls'
        leftNoprojPath        = os.path.join(tempFolder, 'leftFinalCorrected.noproj.cub')
        rightNoprojPath       = os.path.join(tempFolder, 'rightFinalCorrected.noproj.cub')
        leftStereoNoprojPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.noproj.cub')
        rightStereoNoprojPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.noproj.cub')

        noprojCubePair(leftCorrectedPath,        leftNoprojPath,        leftCorrectedPath,       pvlPath, carry)
        noprojCubePair(rightCorrectedPath,       rightNoprojPath,       leftCorrectedPath,       pvlPath, carry)
        noprojCubePair(leftStereoCorrectedPath,  leftStereoNoprojPath,  leftStereoCorrectedPath, pvlPath, carry)
        noprojCubePair(rightStereoCorrectedPath, rightStereoNoprojPath, leftStereoCorrectedPath, pvlPath, carry)


        # Combine the noproj files to make a mosaic.
        # - This also takes care of the cubenorm step.
        # - This step takes a while.
        print 'Starting mosaic calls'
        mainMosaicPath      = os.path.join(tempFolder, 'mainMosaic.cub')
        stereoMosaicPath    = os.path.join(tempFolder, 'stereoMosaic.cub')
        mainMosaicWorkDir   = os.path.join(tempFolder, 'mainMosaicWorkDir/')
        stereoMosaicWorkDir = os.path.join(tempFolder, 'stereoMosaicWorkDir/')
        if not os.path.exists(mainMosaicWorkDir):
            os.mkdir(mainMosaicWorkDir)
        if not os.path.exists(stereoMosaicWorkDir):
            os.mkdir(stereoMosaicWorkDir)


        createMosaic(leftNoprojPath,       rightNoprojPath,       mainMosaicPath,   mainMosaicWorkDir,   carry)
        createMosaic(leftStereoNoprojPath, rightStereoNoprojPath, stereoMosaicPath, stereoMosaicWorkDir, carry)

        # TODO: Testing use of the handmos blend feature!
        #mainMosaicPath   = os.path.join(tempFolder, 'mosaicBlend.cub')
        #stereoMosaicPath = os.path.join(tempFolder, 'mosaicStereoBlend.cub')

        # For testing, crop the mosaic that will be passed into the stereo function to reduce processing time
        if (options.cropAmount > 0):
            mainMosaicCroppedPath   = os.path.join(tempFolder, 'mainMosaicCropped.cub')
            stereoMosaicCroppedPath = os.path.join(tempFolder, 'stereoMosaicCropped.cub')
            if (not os.path.exists(mainMosaicCroppedPath)) or carry:
                cmd = ('crop from= ' + mainMosaicPath   + ' to= ' + mainMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))
                print cmd
                os.system(cmd)
            if (not os.path.exists(stereoMosaicCroppedPath) or carry):
                cmd = ('crop from= ' + stereoMosaicPath + ' to= ' + stereoMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))
                print cmd
                os.system(cmd)
            stereoInputLeft  = mainMosaicCroppedPath
            stereoInputRight = stereoMosaicCroppedPath
        else:
            stereoInputLeft  = mainMosaicPath
            stereoInputRight = stereoMosaicPath

        print '\n-------------------------------------------------------------------------\n'

        # Call stereo to generate a point cloud from the two images
        # - This step takes a really long time.

        stereoOutputPrefix = os.path.join(tempFolder, 'stereoWorkDir/stereo')
        pointCloudPath     = stereoOutputPrefix + '-PC.tif'
        if (not os.path.exists(pointCloudPath)) or carry:
            cmd = ('parallel_stereo --corr-timeout 400 --alignment affineepipolar --subpixel-mode 1' +
                                  ' --disable-fill-holes ' +  stereoInputLeft + ' ' + stereoInputRight + 
                                  ' ' + stereoOutputPrefix + ' --processes 8 --threads-multiprocess 4' +
                                  ' --threads-singleprocess 32 --compute-error-vector')
            print cmd
            os.system(cmd)
            #--nodes-list PBS_NODEFILE --processes 4 --threads-multiprocess 16 --threads-singleprocess 32
        else:
            print 'Stereo file ' + pointCloudPath + ' already exists, skipping stereo step.'

        # Compute percentage of good pixels
        percentGood = IsisTools.getStereoGoodPixelPercentage(stereoOutputPrefix)
        print 'Stereo completed with good pixel percentage: ' + str(percentGood)
        logging.info('Final stereo completed with good pixel percentage: %s', str(percentGood))


        # Find out the center latitude of the mosaic
        centerLat = IsisTools.getCubeCenterLatitude(mainMosaicPath, tempFolder)

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

        # Create a hillshade image to check if the central errors are gone
        hillshadePath = os.path.join(outputFolder, 'outputHillshade.tif')
        if (not os.path.exists(hillshadePath)) or carry:
            cmd = 'hillshade ' + demPath + ' -o ' + hillshadePath
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + hillshadePath + ' already exists, skipping hillshade step.'



        # Call script to compare LOLA data with the DEM
        lolaDiffStatsPath = os.path.join(outputFolder, 'LOLA_diff_stats.txt')
        compareDemToLola(options.lolaPath, demPath, lolaDiffStatsPath, carry)
        

        # Call script to compare LOLA data with the ASU DEM
        lolaAsuDiffStatsPath = os.path.join(outputFolder, 'ASU_LOLA_diff_stats.txt')
        if (options.asuPath):
            compareDemToLola(options.lolaPath, options.asuPath, lolaAsuDiffStatsPath, carry)
        
        

        # Clean up temporary files
    #        if not options.keep:
    #            os.remove(tempTextPath)


        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
