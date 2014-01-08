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
From two pairs of LRONAC images generates two geometrically calibrated mosaic images.
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
    if not os.path.exists(inputPath):
        print ('=======> WARNING: Kml point generation could not find input in input folder ' + inputFolder)
        return False

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

    # Determine working directory
    outputFolder   = os.path.dirname(outputCube)
    outputBaseName = os.path.basename(outputCube)
    folderName     = outputBaseName + '_workDir' 
    tempDirectory  = os.path.join(outputFolder, folderName)

    # Execute command to go to temp directory, execute noproj, then clean up
    cmd = ('mkdir -p '           + tempDirectory + ' && '
           + ' cd '              + tempDirectory + ' && '
           + ' noproj from='     + inputCube
           +       ' match= '    + matchCube
           +       ' specs= '    + pvlPath
           +          ' to='     + outputCube    + ' && '
           + ' cd .. && rm -rf ' + tempDirectory)

    ## Generate the new file
    #cmd = ('noproj from= '  + inputCube + ' to= '    + outputCube + 
    #             ' match= ' + matchCube + ' specs= ' + pvlPath)
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


# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('calibrationReport.py')
    IsisTools.checkIfToolExists('noproj')
    IsisTools.checkIfToolExists('lronacjitreg')
    IsisTools.checkIfToolExists('handmos')
    IsisTools.checkIfToolExists('cubenorm')
    IsisTools.checkIfToolExists('stereoDoubleCalibrationProcess.py')

    return True

#--------------------------------------------------------------------------------

def main():

    print '#################################################################################'
    print "Running lronac2refinedMosaics.py"

    try:
        try:
            usage = "usage: lronac2refinedMosaics.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
            inputGroup = optparse.OptionGroup(parser, 'Input Paths')
            inputGroup.add_option("--left",  dest="leftPath",  help="Path to LE .IMG file")
            inputGroup.add_option("--right", dest="rightPath", help="Path to RE .IMG file")            
            inputGroup.add_option("--stereo-left",  dest="stereoLeft", 
                                  help="Path to LE .IMG file with overlapping view of --left file")
            inputGroup.add_option("--stereo-right", dest="stereoRight", 
                                  help="Path to RE .IMG file with overlapping view of --right file")

            inputGroup.add_option("--lola",    dest="lolaPath", help="Path to LOLA DEM")
            
            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--output-folder", dest="outputFolder",  
                              help="Output directory")
            parser.add_option("--workDir", dest="workDir",  
                              help="Folder to put temporary files in")
            parser.add_option("--log-path",  dest="logPath",        
                              help="Where to write the output log file.")

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
            if not options.outputFolder: 
                parser.error("Need output folder")

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
        if not os.path.exists(options.stereoLeft):
            raise Exception('Input file ' + options.stereoLeft + ' not found!')
        if not os.path.exists(options.stereoRight):
            raise Exception('Input file ' + options.stereoRight + ' not found!')
        if not os.path.exists(options.lolaPath):
            raise Exception('Input file ' + options.lolaPath + ' not found!')

        # Set up the output folders
        outputFolder  = options.outputFolder
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder) 
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

        # Set up logging
        if not options.logPath:
            options.logPath = os.path.join(options.outputFolder, '/lronac2refinedMosaicLog.txt')
        logging.basicConfig(filename=options.logPath,level=logging.INFO)    

        # Set up final output paths
        filename         = os.path.splitext(options.leftPath)[0] + '.correctedMosaic.cub'
        outputPathMain   = os.path.join(outputFolder, os.path.basename(filename))
        filename         = os.path.splitext(options.stereoLeft)[0] + '.correctedMosaic.cub'
        outputPathStereo = os.path.join(outputFolder, os.path.basename(filename))

        # Set up output paths for the stereo calibration call
        leftCorrectedPath        = os.path.join(tempFolder, 'leftFinalCorrected.cub')
        rightCorrectedPath       = os.path.join(tempFolder, 'rightFinalCorrected.cub')
        leftStereoCorrectedPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.cub')
        rightStereoCorrectedPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.cub')

        # Generate a kml plot of the input LOLA data
        lolaKmlPath = os.path.join(tempFolder, 'lolaRdrPoints.kml')
        if not os.path.exists(lolaKmlPath):
            cmd = ('calibrationReport.py --input '  + options.lolaPath + 
                                       ' --output ' + lolaKmlPath + 
                                       ' --name '   + 'lolaRdrPoints' + 
                                       ' --skip '   + str(100) + ' --color ' + 'blue' + ' --size small')
            print cmd
            os.system(cmd)

        carry = False

        # Set corrected image paths
        filename                 = os.path.splitext(options.leftPath)[0]    + '.geoCorrected.cub'
        leftCorrectedPath        = os.path.join(outputFolder, os.path.basename(filename))
        filename                 = os.path.splitext(options.rightPath)[0]   + '.geoCorrected.cub'
        rightCorrectedPath       = os.path.join(outputFolder, os.path.basename(filename))
        filename                 = os.path.splitext(options.stereoLeft)[0]  + '.geoCorrected.cub'
        leftStereoCorrectedPath  = os.path.join(outputFolder, os.path.basename(filename))
        filename                 = os.path.splitext(options.stereoRight)[0] + '.geoCorrected.cub'
        rightStereoCorrectedPath = os.path.join(outputFolder, os.path.basename(filename))

        # Correct all four input images at once
        caughtException = False
        doubleCalWorkFolder = os.path.join(tempFolder, 'doubleCal')
        try:
          if (not os.path.exists(leftCorrectedPath)) or carry:
              print '\n=============================================================================\n'

              cmd = ('stereoDoubleCalibrationProcess.py --left '          + options.leftPath + 
                                                      ' --right '         + options.rightPath + 
                                                      ' --stereo-left '   + options.stereoLeft + 
                                                      ' --stereo-right '  + options.stereoRight + 
                                                      ' --keep --lola '   + options.lolaPath + 
                                                      ' --output-folder ' + outputFolder + 
                                                      ' --workDir '       + doubleCalWorkFolder + 
                                                      ' --log-path '      + options.logPath)
              print cmd
              os.system(cmd)
              print '\n============================================================================\n'

        except: 
            caughtException = True
            print 'Caught an exception!'

        # Convert GDC output files into KML plots 
        # - This is just to help with debugging
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'initialGdcCheck'),            tempFolder, 'pairGdcCheckInitial.csv',           1,    'blue', 'normal', carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'posCorrectGdcCheck'),         tempFolder, 'pairGdcCheckPos.csv',               1,    'green', 'normal',  carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'posCorrectStereoGdcCheck'),   tempFolder, 'pairGdcStereoCheckPos.csv',         1,    'green', 'normal',  carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'stereoGlobalAdjustGdcCheck'), tempFolder, 'pairGdcCheckGlobalAdjustStero.csv', 1,    'blue', 'normal', carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'pcAlignStereoGdcCheck'),      tempFolder, 'pairGdcCheckPcAlign.csv',           1,    'red', 'normal',  carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'finalGdcCheck'),              tempFolder, 'pairGdcCheckFinal.csv',             1,    'white', 'normal', carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'finalStereoGdcCheck'),        tempFolder, 'pairGdcCheckFinalStereo.csv',       1,    'white', 'normal',  carry)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'gdcPointsLargeComp'),         tempFolder, 'inputGdcPoints.csv',                1000, 'blue', 'tiny', carry)
        #generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'initialGdcCheck'),            tempFolder, 'dem-trans_source.csv',              1000, 'blue', 'normal', , False)
        generateKmlFromGdcPoints(os.path.join(doubleCalWorkFolder, 'pcAlignOutput'),              tempFolder, 'dem-trans_reference.csv',           1000, 'red', 'tiny',  carry)

        print 'Finished generating KML plots'


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

        correctTime = time.time()
        logging.info('Nav correction finished in %f seconds', correctTime - startTime)

        # Noproj the corrected data
        leftNoprojPath        = os.path.join(tempFolder, 'leftFinalCorrected.noproj.cub')
        rightNoprojPath       = os.path.join(tempFolder, 'rightFinalCorrected.noproj.cub')
        leftStereoNoprojPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.noproj.cub')
        rightStereoNoprojPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.noproj.cub')

        # Set up thread objects
        leftThread        = threading.Thread(target=noprojCubePair, 
                                               args=(leftCorrectedPath, leftNoprojPath, 
                                                     leftCorrectedPath, pvlPath, carry))
        rightThread       = threading.Thread(target=noprojCubePair, 
                                               args=(rightCorrectedPath, rightNoprojPath, 
                                                     leftCorrectedPath, pvlPath, carry))
        leftStereoThread  = threading.Thread(target=noprojCubePair, 
                                             args  =(leftStereoCorrectedPath, leftStereoNoprojPath, 
                                                     leftStereoCorrectedPath, pvlPath, carry))
        rightStereoThread = threading.Thread(target=noprojCubePair, 
                                             args  =(rightStereoCorrectedPath, rightStereoNoprojPath, 
                                                     leftStereoCorrectedPath, pvlPath, carry))
        print 'Starting noproj call threads'
        leftThread.start()
        rightThread.start()
        leftStereoThread.start()
        rightStereoThread.start()

        print 'Waiting for noproj threads to complete...'
        leftThread.join()
        rightThread.join()
        leftStereoThread.join()
        rightStereoThread.join()

        print 'noproj threads finished.'

        noprojTime = time.time()
        logging.info('noproj finished in %f seconds', noprojTime - correctTime)

        # Combine the noproj files to make a mosaic.
        # - This also takes care of the cubenorm step.
        # - This step takes a while.
        print 'Starting mosaic calls'
        mainMosaicWorkDir   = os.path.join(tempFolder, 'mainMosaicWorkDir/')
        stereoMosaicWorkDir = os.path.join(tempFolder, 'stereoMosaicWorkDir/')
        if not os.path.exists(mainMosaicWorkDir):
            os.mkdir(mainMosaicWorkDir)
        if not os.path.exists(stereoMosaicWorkDir):
            os.mkdir(stereoMosaicWorkDir)

        # Set up thread objects
        leftThread  = threading.Thread(target=createMosaic, 
                                         args=(leftNoprojPath, rightNoprojPath, 
                                               outputPathMain, mainMosaicWorkDir, carry))
        rightThread = threading.Thread(target=createMosaic, 
                                         args=(leftStereoNoprojPath, rightStereoNoprojPath, 
                                               outputPathStereo,     stereoMosaicWorkDir,  carry))
        print 'Starting mosaic call threads'
        leftThread.start()
        rightThread.start()

        print 'Waiting for mosaic threads to complete...'
        leftThread.join()
        rightThread.join()

        mosaicTime = time.time()
        logging.info('Mosaics finished in %f seconds', mosaicTime - noprojTime)

        # Clean up temporary files
    #        if not options.keep:
    #            os.remove(tempTextPath)


        endTime = time.time()

        logging.info('lronac2refinedMosaics finished in %f seconds', endTime - startTime)
        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
