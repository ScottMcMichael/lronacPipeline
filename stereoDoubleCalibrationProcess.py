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

import os, glob, optparse, re, shutil, subprocess, string, time, math, logging, threading

import IsisTools

import positionCorrector, rotationCorrector, lronacCameraRotationCorrector

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Corrects the geometry of two pairs of input LRONAC images to minimize projection errors.  
Uses pixel pairs between all images and a LOLA point cloud which must be in the correct region.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#===============================================================

# TODO: Move this!
# Counts the number of lines in a file
def getFileLineCount(filePath):
    f = open(filePath)
    i = 0
    for line in f:
        i = i + 1
    return i

# TODO: Make this a standalone function!
# Get a single file ready to process
def prepareImgFile(inputPath, outputFolder, keep):

    # Generate the appropriate paths in the output folder
    filename       = os.path.splitext(inputPath)[0] + '.cub'
    initialCubFile = os.path.join(outputFolder, os.path.basename(filename))
    filename       = os.path.splitext(inputPath)[0] + '.lronaccal.cub'
    calCubFile     = os.path.join(outputFolder, os.path.basename(filename))
    filename       = os.path.splitext(inputPath)[0] + '.lronaccal.lronacecho.cub'
    echoCubFile    = os.path.join(outputFolder, os.path.basename(filename))

    # Quit immediately if the output file already exists
    if os.path.exists(echoCubFile):
        print 'File ' + echoCubFile + ' already exists, skipping init calls'
        return echoCubFile

    if not os.path.exists(initialCubFile):
        cmd = 'lronac2isis from= ' + inputPath + ' to= ' + initialCubFile
        print cmd
        os.system(cmd)

    if not os.path.exists(calCubFile):
        cmd = 'lronaccal from= ' + initialCubFile + ' to= ' + calCubFile
        print cmd
        os.system(cmd)

    if not os.path.exists(echoCubFile):
        cmd = 'lronacecho from= ' + calCubFile + ' to= ' + echoCubFile
        print cmd
        os.system(cmd)

    cmd = 'spiceinit from= ' + echoCubFile
    print cmd
    os.system(cmd)

    if not keep:
        IsisTools.removeIfExists(initialCubFile)
        IsisTools.removeIfExists(calCubFile)

    return echoCubFile


# Apply the position offset specified in the IK kernel to either an LE or RE camera.
def applyInterCameraPositionOffset(inputCubePath, outputCubePath, workingDirectory, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping position offset correction.'
        return True

    # Run the process
    cmdArgs = ['--keep', '--input', inputCubePath, 
               '--output',  outputCubePath, 
               '--workDir', workingDirectory]
    print cmdArgs
    positionCorrector.main(cmdArgs)


    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Position offset correction failed to create output file ' + outputCubePath + 
                        ' from input file ' + inputCubePath)

    return True



# Samples pixels pairs for bundle adjustment from the output of the stereo command.
# - Returns the output path.
def extractPixelPairsFromStereoResults(disparityImagePath, outputDirectory, outputFileName, 
                                       sampleInterval, forceOperation):

    # Quit immediately if the output file already exists
    outputPixelPath = os.path.join(outputDirectory, outputFileName)
    if (not forceOperation) and (os.path.exists(outputPixelPath)):
        print 'File ' + outputPixelPath + ' already exists, skipping pixel pair sampling.'
        return outputPixelPath

    # Run the process
    cmd = ('pixelPairsFromStereo ' + disparityImagePath + ' ' + outputPixelPath + 
          ' -p ' + str(sampleInterval))
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputPixelPath):
        raise Exception('Pixel sampling failed to create output file ' + outputPixelPath + 
                        ' from input file ' + disparityImagePath)

    numPairs = getFileLineCount(outputPixelPath)
    logging.info('In stereo file %s found %d point pairs', disparityImagePath, numPairs)

    return outputPixelPath


# Applies a planet-centered rotation and correction to the nav data of a cube
# - If the CK and SPK paths are not specified, paths are automatically generated.
# - Set forceOperation to run the operation even if the output file already exists
def applyNavTransform(inputCubePath, outputCubePath, transformMatrixPath, workDir, 
                      ckPath, spkPath, pcAlignTrans, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
        return True

    
    # Set up the transformation command
    cmdArgs = ['--keep', '--input', inputCubePath, '--output', outputCubePath, 
               '--transformPath', transformMatrixPath, '--workDir', workDir]
    
    # Set up CK and SPK manual paths if they were specified
    if (ckPath):
        cmdArgs.append('--ck')
        cmdArgs.append(ckPath)
    if (spkPath):
        cmdArgs.append('--spk')
        cmdArgs.append(spkPath)

    # Set flag if the transform comes from pc_align as opposed to the lronacAngleSolver
    if (pcAlignTrans):
        cmdArgs.append('--pcAlignTrans')
    
    # Execute the transformation command
    print cmdArgs
    rotationCorrector.main(cmdArgs)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Nav transform failed to create output file ' + outputCubePath + 
                        ' from input file ' + inputCubePath)

    return True


# Tries to compute the internal angle between an LE/RE image pair.
# - Output GDC points serve as a check to make sure the images are in roughly the correct place.
def checkAdjacentPairAlignment(leftInputPath, rightInputPath, outputDirectory,  surfaceElevation=0, forceOperation=False):

    # Figure out output paths
    if not os.path.exists(outputDirectory):
        os.mkdir(outputDirectory)
    sbaOutputPrefix = os.path.join(outputDirectory, 'SBA_check')
    defaultGdcPath  = sbaOutputPrefix + '-outputGdcPoints.csv'

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(defaultGdcPath)):
        return True

    # Run the process
    cmd = ['lronacAngleDoubleSolver',  '--outputPrefix',            sbaOutputPrefix, 
                                       '--leftCubePath',            leftInputPath, 
                                       '--rightCubePath',           rightInputPath, 
                                       '--elevation',               str(surfaceElevation)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outputText, err = p.communicate()

    # Check to make sure we actually created the file
    if not os.path.exists(defaultGdcPath):
        raise Exception('Adjacency check failed to create output file ' + defaultGdcPath + 
                        ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    return True


# Generate a modified IK kernel to adjust the rotation between an LE/RE camera pair.
def applyInterCameraPairRotation(leftInputPath, rightInputPath, newRotationPath, outputCubePath, 
                                 ckPath, spkPath, workDir, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
        return True

    # Generate the new file
    cmdArgs = ['--keep', '--output', outputCubePath, 
               '--rotation', newRotationPath, '--left', leftInputPath, 
               '--right', rightInputPath, '--ck', ckPath, '--spk', spkPath,
               '--workDir', workDir]
    print cmdArgs
    lronacCameraRotationCorrector.main(cmdArgs)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Inter-camera rotation failed to create output file ' + outputCubePath + 
                        ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    return True


# Calls stereo functions to generate a disparity image and returns the path to it.
def callStereoCorrelation(leftInputPath, rightInputPath, outputPrefix, correlationTimeout, forceOperation):

    # Quit immediately if the output file already exists
    disparityImagePath = outputPrefix + '-D.tif'
    if (not forceOperation) and (os.path.exists(disparityImagePath)):
        print 'File ' + disparityImagePath + ' already exists, skipping stereo computation.'
        return disparityImagePath

#    # Use parallel stereo call, steps 0 and 1 only.  Other options to try and increase speed.
#    cmd = ('parallel_stereo  --corr-max-levels 3 --compute-error-vector --entry-point 0 --stop-point 2' + 
#           ' --alignment none --subpixel-mode 1 --disable-fill-holes --processes 8 ' + 
#           ' --threads-multiprocess 4 --threads-singleprocess 32 --cost-mode 0 --corr-timeout ' + 
#           str(correlationTimeout) + ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix)
#    print cmd
#    os.system(cmd)

    # Call only the initial stereo stages to cut down on processing time

    # Stage 0 (fast)
    cmd = ('stereo_pprc --corr-max-levels 3 --compute-error-vector --cost-mode 0 --alignment none' + 
                      ' --subpixel-mode 1 --disable-fill-holes --corr-timeout ' + str(correlationTimeout) + 
                      ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix)
    print cmd
    os.system(cmd)
    # Stage 1 (slow!)
    cmd = ('stereo_corr --corr-max-levels 3 --compute-error-vector --cost-mode 0 --alignment none' + 
                      ' --subpixel-mode 1 --disable-fill-holes --corr-timeout ' + str(correlationTimeout) + 
                      ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix)
    print cmd
    os.system(cmd)



    # Generate disparity images for debugging
    #cmd = 'disparitydebug ' + disparityImagePath
    #print cmd
    #os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(disparityImagePath):
        raise Exception('Stereo processing failed to create output file ' + disparityImagePath + 
                        ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    # Compute percentage of good pixels
    percentGood = IsisTools.getStereoGoodPixelPercentage(outputPrefix)
    print 'Stereo completed with good pixel percentage: ' + str(percentGood)
    logging.info('For cubes %s and %s', leftInputPath, rightInputPath)
    logging.info('- Stereo completed with good pixel percentage: %s', str(percentGood))

    return disparityImagePath


# Quickly obtains a list of matched pixel pairs between two images
def getInterestPointPairs(leftInputPath, rightInputPath, outputPath, surfaceElevation=0, forceOperation=False):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputPath)):
        print 'File ' + outputPath + ' already exists, skipping point pair calculation.'
    else:

        # Generate intermediate binary file
        binaryPath = outputPath + '.bin'
        cmd = ('stereoIpFind '+ leftInputPath + ' ' + rightInputPath + ' ' + binaryPath +
                              ' --elevationGuess ' + str(surfaceElevation))
        print cmd
        os.system(cmd)

        if not os.path.exists(binaryPath):
            raise Exception('stereoIpFind failed to create output file ' + outputPath + 
                            ' from input files ' + leftInputPath + ' and ' + rightInputPath)

        # Convert from the binary file to an easy to read CSV final
        cmd = 'matchBinaryToCsv '+ binaryPath + ' -o ' + outputPath
        print cmd
        os.system(cmd)

        # Check to make sure we actually created the file
        if not os.path.exists(outputPath):
            raise Exception('matchBinaryToCsv call failed on file ' + outputPath)

        # Clean up temporary file
        IsisTools.removeIfExists(binaryPath)

    # Count the number of point pairs we found
    numPairs = getFileLineCount(outputPath)
    print('For cubes %s and %s found %d point pairs' % (leftInputPath, rightInputPath, numPairs))
    logging.info('For cubes %s and %s', leftInputPath, rightInputPath)
    logging.info('- Found %d point pairs', numPairs)

    return numPairs


# Compares the backprojected locations of matched pixel pairs
def evaluateAccuracy(leftCubePath, rightCubePath, ipFindOutputPath, workDir=''):

    print 'Evaluating accuracy between cubes ' + leftCubePath + ' and ' + rightCubePath

    # Run check every N lines in the file
    lineSkip = 15
    
    f = open(ipFindOutputPath)
    i = 0
    sumDistance = 0.0
    count       = 0.0
    for line in f:
        if i % lineSkip == 0: # Compare locations on a sampling of lines
            # Get the matching pixel coordinates from this line
            leftSample, leftLine, rightSample, rightLine = line.split(',')
            
            #print 'LeftPixel =  ' + leftSample  + ', ' + leftLine
            #print 'RightPixel = ' + rightSample + ', ' + rightLine
            
            # Obtain the backprojected GCC location for the pixels in the two images
            leftLoc  = IsisTools.getPixelLocInCube(leftCubePath,  leftSample,  leftLine,  workDir)
            rightLoc = IsisTools.getPixelLocInCube(rightCubePath, rightSample, rightLine, workDir)
            
            #print 'leftLoc  = ' + str(leftLoc)
            #print 'rightLoc = ' + str(rightLoc)
            
            # Determine the distance between the GCC points and accumulate
            distance = math.sqrt(math.pow(leftLoc[0] - rightLoc[0], 2.0) + 
                                 math.pow(leftLoc[1] - rightLoc[1], 2.0) + 
                                 math.pow(leftLoc[2] - rightLoc[2], 2.0))
            #print 'Distance = %.2f' % distance
            sumDistance = sumDistance + distance
            count = count + 1.0
        i = i + 1
   
    if (i == 0):
        raise Exception("Can't evaluate pixel pairs, match file " + ipFindOutputPath + ' is empty!')

    # Determine the mean distance
    meanDistance = sumDistance / count
    return meanDistance


# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('pc_align')
    IsisTools.checkIfToolExists('lronacAngleDoubleSolver')
    IsisTools.checkIfToolExists('lronac2isis')
    IsisTools.checkIfToolExists('lronaccal')
    IsisTools.checkIfToolExists('lronacecho')
    IsisTools.checkIfToolExists('spiceinit')
    IsisTools.checkIfToolExists('positionCorrector.py')
    IsisTools.checkIfToolExists('pixelPairsFromStereo')
    IsisTools.checkIfToolExists('rotationCorrector.py')
    IsisTools.checkIfToolExists('lronacCameraRotationCorrector.py')
    IsisTools.checkIfToolExists('stereo_corr')
    IsisTools.checkIfToolExists('stereo_pprc')
    IsisTools.checkIfToolExists('crop')

    return True

# Looks in the pc_align output folder for the variably-named log file
def findOutputLog(folder):

    for root, dirs, files in os.walk(folder):
        for f in files:
            if 'log' in f:
                return os.path.join(folder, f)
    return False

#==========================================================================================

def main(argsIn):

    print "Started stereoDoubleCalibrationProcess.py"

    try:
        try:
            usage = "usage: stereoDoubleCalibrationProcess.py TODO [--manual]\n  "
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
            parser.add_option_group(inputGroup)
  
            outputGroup = optparse.OptionGroup(parser, 'Output Paths')

            outputGroup.add_option("--output-folder", dest="outputFolder",  
                                   help="Output folder to store results in")
            outputGroup.add_option("--log-path",  dest="logPath",        
                                   help="Where to write the output log file.")
            parser.add_option_group(outputGroup)
  
            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")
          
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args(argsIn)

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.stereoLeft: 
                parser.error("Need stereo left input path")
            if not options.stereoRight: 
                parser.error("Need stereo right input path")
            if not options.outputFolder: 
                parser.error("Need output folder")
            if not options.lolaPath: 
                parser.error("Need LOLA DEM path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        # Make sure we have all the functions we need
        functionStartupCheck()

        # Set this to true to force steps after it
        carry = False

        outputFolder  = options.outputFolder
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
            options.logPath = options.workDir + '/stereoDoubleCalLog.txt'
        logging.basicConfig(filename=options.logPath,level=logging.INFO)      
        

        # Set output paths
        leftBaseName        = os.path.basename(os.path.splitext(options.leftPath)[0])
        rightBaseName       = os.path.basename(os.path.splitext(options.rightPath)[0])
        leftStereoBaseName  = os.path.basename(os.path.splitext(options.stereoLeft)[0])
        rightStereoBaseName = os.path.basename(os.path.splitext(options.stereoRight)[0])

        outputPathLeft        = os.path.join(outputFolder, leftBaseName        + '.geoCorrected.cub')
        outputPathRight       = os.path.join(outputFolder, rightBaseName       + '.geoCorrected.cub')
        outputPathStereoLeft  = os.path.join(outputFolder, leftStereoBaseName  + '.geoCorrected.cub')
        outputPathStereoRight = os.path.join(outputFolder, rightStereoBaseName + '.geoCorrected.cub')


        # Convert the input files from IMG files to spiceinit'ed cubes in the output folder
        leftThread        = threading.Thread(target=prepareImgFile, 
                                               args=(options.leftPath,    tempFolder, options.keep))
        rightThread       = threading.Thread(target=prepareImgFile, 
                                               args=(options.rightPath,   tempFolder, options.keep))
        leftStereoThread  = threading.Thread(target=prepareImgFile, 
                                               args=(options.stereoLeft,  tempFolder, options.keep))
        rightStereoThread = threading.Thread(target=prepareImgFile, 
                                               args=(options.stereoRight, tempFolder, options.keep))

        print 'Starting data init threads'
        leftThread.start()
        rightThread.start()
        leftStereoThread.start()
        rightStereoThread.start()

        print 'Waiting for data init threads to complete...'
        leftThread.join()
        rightThread.join()
        leftStereoThread.join()
        rightStereoThread.join()

        print 'data init threads finished.'

        # Get the file paths generated by prepareImgFile
        spiceInitLeftPath        = os.path.join(tempFolder, leftBaseName        + '.lronaccal.lronacecho.cub')
        spiceInitRightPath       = os.path.join(tempFolder, rightBaseName       + '.lronaccal.lronacecho.cub')
        spiceInitStereoLeftPath  = os.path.join(tempFolder, leftStereoBaseName  + '.lronaccal.lronacecho.cub')
        spiceInitStereoRightPath = os.path.join(tempFolder, rightStereoBaseName + '.lronaccal.lronacecho.cub')


        initTime = time.time()
        logging.info('Init complete in %f seconds', initTime - startTime)


        # Get the expected surface elevation in meters
        expectedSurfaceElevation = IsisTools.getCubeElevationEstimate(spiceInitLeftPath, tempFolder)

        # DEBUG: Check angle solver on input LE/RE images!
#        checkAdjacentPairAlignment(spiceInitLeftPath, spiceInitRightPath, 
#                                   os.path.join(tempFolder, 'initialGdcCheck'), 
#                                   expectedSurfaceElevation, carry)

        # Apply LE/RE LRONAC position offsets to each of the input files
        posOffsetCorrectedLeftPath = os.path.join(tempFolder, leftBaseName + '.posOffsetCorrected.cub')
        leftPosCorrectWorkDir      = os.path.join(tempFolder, 'leftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitLeftPath, posOffsetCorrectedLeftPath, leftPosCorrectWorkDir, carry)

        posOffsetCorrectedRightPath = os.path.join(tempFolder, rightBaseName + '.posOffsetCorrected.cub')
        rightPosCorrectWorkDir      = os.path.join(tempFolder, 'rightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitRightPath, posOffsetCorrectedRightPath, rightPosCorrectWorkDir, carry)

        posOffsetCorrectedStereoLeftPath = os.path.join(tempFolder, leftStereoBaseName + '.posOffsetCorrected.cub')
        stereoLeftPosCorrectWorkDir      = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoLeftPath, posOffsetCorrectedStereoLeftPath, 
                                       stereoLeftPosCorrectWorkDir, carry)

        posOffsetCorrectedStereoRightPath = os.path.join(tempFolder, rightStereoBaseName + '.posOffsetCorrected.cub')
        stereoRightPosCorrectWorkDir      = os.path.join(tempFolder, 'stereoRightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoRightPath, posOffsetCorrectedStereoRightPath,  
                                       stereoRightPosCorrectWorkDir, carry)

        # DEBUG: Check angle solver on input LE/RE images!
#        checkAdjacentPairAlignment(posOffsetCorrectedLeftPath, posOffsetCorrectedRightPath, 
#                                   os.path.join(tempFolder, 'posCorrectGdcCheck'), 
#                                   expectedSurfaceElevation, carry)
#        checkAdjacentPairAlignment(posOffsetCorrectedStereoLeftPath, posOffsetCorrectedStereoRightPath, 
#                                   os.path.join(tempFolder, 'posCorrectStereoGdcCheck'), 
#                                   expectedSurfaceElevation, carry)

        positionTime = time.time()
        logging.info('Position offset complete in %f seconds', positionTime - initTime)

        print '\n-------------------------------------------------------------------------\n'

        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefixLeft   = os.path.join(tempFolder, 'stereoOutputLeft/out')
        disparityImageLeft = callStereoCorrelation(posOffsetCorrectedLeftPath, 
                                                   posOffsetCorrectedStereoLeftPath, 
                                                   stereoPrefixLeft, 400, carry)

        # Extract a small number of matching pixel locations from the LE and RE disparity images ( < 300 pairs)
        pixelPairsLeftSmall  = extractPixelPairsFromStereoResults(disparityImageLeft, tempFolder, 
                                                                  'stereoPixelPairsLeftSmall.csv', 
                                                                  800, carry)


        # Get small number of matching pixels for the right side quickly
        pixelPairsRightSmall = os.path.join(tempFolder, 'stereoPixelPairsRightSmall.csv')
        getInterestPointPairs(posOffsetCorrectedRightPath, posOffsetCorrectedStereoRightPath, 
                              pixelPairsRightSmall, expectedSurfaceElevation, carry)

        print '\n-------------------------------------------------------------------------\n'

        # Perform cross-stereo matching of LE/RE cubes from opposite pairs

        # Crop each input cube to half its width
        # - Take the right half of the LE cubes and the left half of the RE cubes
        leftPosCorrectedCropped        = os.path.join(tempFolder, leftBaseName        + '.cropped.cub')
        rightPosCorrectedCropped       = os.path.join(tempFolder, rightBaseName       + '.cropped.cub')
        leftStereoPosCorrectedCropped  = os.path.join(tempFolder, leftStereoBaseName  + '.cropped.cub')
        rightStereoPosCorrectedCropped = os.path.join(tempFolder, rightStereoBaseName + '.cropped.cub')

        # Check the image size and determine if it is full resolution
        imageSize     = IsisTools.getCubeSize(posOffsetCorrectedLeftPath)
        print 'Input image size = ' + str(imageSize)
        logging.info('Image size = %s', str(imageSize))
        if (imageSize[0] < 5000):
            logging.info('Stopping processing on half-width image until they are tested!')
            raise Exception('Stopping processing on half-width image until they are tested!')
        
        
        cropWidth     = imageSize[0] / 2
        isisCropStart = cropWidth + 1 # One-based starting pixel for ISIS LE crop

        if not os.path.exists(leftPosCorrectedCropped):
            cmd = ('crop sample=' + str(isisCropStart) + ' from= ' + posOffsetCorrectedLeftPath + 
                       ' to= ' + leftPosCorrectedCropped)
            print cmd
            os.system(cmd)
        if not os.path.exists(rightPosCorrectedCropped):
            cmd = ('crop nsamples=' + str(cropWidth) + ' from= ' + posOffsetCorrectedRightPath + 
                       ' to= ' + rightPosCorrectedCropped)
            print cmd
            os.system(cmd)
        if not os.path.exists(leftStereoPosCorrectedCropped):
            cmd = ('crop sample=' + str(isisCropStart) + ' from= ' + posOffsetCorrectedStereoLeftPath + 
                       ' to= ' + leftStereoPosCorrectedCropped)
            print cmd
            os.system(cmd)
        if not os.path.exists(rightStereoPosCorrectedCropped):
            cmd = ('crop nsamples=' + str(cropWidth) + ' from= ' + posOffsetCorrectedStereoRightPath + 
                       ' to= ' + rightStereoPosCorrectedCropped)
            print cmd
            os.system(cmd)

        # Get small numbers of pixels for the two cross pairs quickly

#        # Use stereo command on two cross-pair cubes
#        # - Timeouts are shorter on these since there is probably less image overlap
        
        # First is left in main pair to right in the stereo pair
        numLeftCrossPairs = 0
        pixelPairsLeftCrossSmall = os.path.join(tempFolder, 'stereoPixelPairsLeftCrossSmall.csv')
        tempLeftCrossPixelPairs  = os.path.join(tempFolder, 'tempLeftCrossPixelPairs.csv')
        try:
            numLeftCrossPairs  = getInterestPointPairs(leftPosCorrectedCropped, rightStereoPosCorrectedCropped, 
                                                       tempLeftCrossPixelPairs, expectedSurfaceElevation, carry)
            # Add in offset to the LE image so that the pixel coordinates are in the full, not cropped, frame
            IsisTools.modifyPixelPairs(tempLeftCrossPixelPairs, pixelPairsLeftCrossSmall, cropWidth, 0, 0, 0)
            usingLeftCross = True
        except:
            print 'Failed to find left-cross match, ignoring this data source.'
            usingLeftCross = False

        # Next is left in the stereo pair to right in the main pair
        numRightCrossPairs = 0
        pixelPairsRightCrossSmall = os.path.join(tempFolder, 'stereoPixelPairsRightCrossSmall.csv')
        tempRightCrossPixelPairs  = os.path.join(tempFolder, 'tempRightCrossPixelPairs.csv')
        try:
            numRightCrossPairs = getInterestPointPairs(leftStereoPosCorrectedCropped, rightPosCorrectedCropped, 
                                                       tempRightCrossPixelPairs, expectedSurfaceElevation, carry)
            # Add in offset to the LE image so that the pixel coordinates are in the full, not cropped, frame
            IsisTools.modifyPixelPairs(tempRightCrossPixelPairs, pixelPairsRightCrossSmall, cropWidth, 0, 0, 0)
            usingRightCross = True
        except:
            print 'Failed to find right-cross match, ignoring this data source.'
            usingRightCross = False

        # Left and right cross pixels may not be used depending on the image overlap
        # - If a small number of pixel pairs was found that suggests we may have a high ratio of bad pixels
        leftCrossElems   = []
        rightCrossElems  = []
        if usingLeftCross and (numLeftCrossPairs > 50):
            leftCrossElems   = ['--matchingPixelsLeftCrossPath', pixelPairsLeftCrossSmall]
        if usingRightCross and (numRightCrossPairs > 50):
            rightCrossElems  = ['--matchingPixelsRightCrossPath', pixelPairsRightCrossSmall]

        pixelTime = time.time()
        logging.info('Pixel pair finding complete in %f seconds', pixelTime - positionTime)

        print '\n-------------------------------------------------------------------------\n'



 
#        meanLeftError   = evaluateAccuracy(posOffsetCorrectedLeftPath,       posOffsetCorrectedStereoLeftPath,
#                                           pixelPairsLeftSmall,  os.path.join(tempFolder, 'finalGdcCheck'))
#        meanRightError  = evaluateAccuracy(posOffsetCorrectedRightPath,      posOffsetCorrectedStereoRightPath, 
#                                           pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheck'))
 
#        print '-----> POS OFFSET'
#        print '=====> Mean left   pair error = %.4f meters' % meanLeftError
#        print '=====> Mean right  pair error = %.4f meters' % meanRightError


        # Compute all rotations and translations between the four cubes
        # - Computes the following transforms:
        #   = main   LE   to main   RE   (camera rotation offset in IK file)
        #   = main   pair to stereo pair (applied to ck and spk kernel files)
        #   = stereo LE   to stereo RE   (camera rotation offset in IK file)
        sbaOutputPrefix     = os.path.join(tempFolder, 'SBA_solution')
        smallGdcFile        = sbaOutputPrefix + "-outputGdcPoints.csv"
        localRotationPath   = sbaOutputPrefix + "-localRotationMatrix.csv"
        globalTransformPath = sbaOutputPrefix + "-globalTransformMatrix.csv"
        stereoRotationPath  = sbaOutputPrefix + "-stereoLocalRotationMatrix.csv"
        solvedParamsPath    = sbaOutputPrefix + "-finalParamState.csv"
        mainIpFindPath      = sbaOutputPrefix + "-mainIpFindPixels.csv"
        stereoIpFindPath    = sbaOutputPrefix + "-stereoIpFindPixels.csv"
        if (not os.path.exists(globalTransformPath)) or carry:
            cmd = ['lronacAngleDoubleSolver',  '--outputPrefix',            sbaOutputPrefix, 
                                               '--matchingPixelsLeftPath',  pixelPairsLeftSmall, 
                                               '--matchingPixelsRightPath', pixelPairsRightSmall, 
                                               '--leftCubePath',            posOffsetCorrectedLeftPath, 
                                               '--rightCubePath',           posOffsetCorrectedRightPath, 
                                               '--leftStereoCubePath',      posOffsetCorrectedStereoLeftPath, 
                                               '--rightStereoCubePath',     posOffsetCorrectedStereoRightPath, 
                                               '--elevation',               str(expectedSurfaceElevation)]
            cmd = cmd + leftCrossElems + rightCrossElems
            print cmd
            print '-------'
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            outputText, err = p.communicate()
            
            # Extract pertinent output information and log it
            initialErrorLine = outputText.find('>>>>')
            initialLineEnd   = outputText.find('\n', initialErrorLine)
            medianErrorLine  = outputText.find('>>>>', initialLineEnd)
            medianLineEnd    = outputText.find('\n', medianErrorLine)
            meanErrorLine    = outputText.find('>>>>', medianLineEnd)
            meanLineEnd      = outputText.find('\n', meanErrorLine)
            changeErrorLine  = outputText.find('>>>>', meanLineEnd)
            changeLineEnd    = outputText.find('\n', changeErrorLine)            
            
            logging.info(outputText[initialErrorLine:initialLineEnd])
            logging.info(outputText[medianErrorLine :medianLineEnd ])
            logging.info(outputText[meanErrorLine   :meanLineEnd   ])
            logging.info(outputText[changeErrorLine :changeLineEnd ])
            
            
            print '====='
            print outputText

        else:
            print 'Skipping stereo transform calculation step'

        sbaTime = time.time()
        logging.info('SBA complete in %f seconds', sbaTime - pixelTime)

        # Apply the planet-centered rotation/translation to both cameras in the stereo pair.
        # - This corrects the stereo pair relative to the main pair.
        # - The RE relative to LE corrections are performed later for convenience.
        leftStereoAdjustedPath  = os.path.join(tempFolder, leftStereoBaseName + '.stereoAdjusted.cub')
        leftSteroCorrectWorkDir = os.path.join(tempFolder, 'stereoLeftStereoCorrection/')
        applyNavTransform(posOffsetCorrectedStereoLeftPath, leftStereoAdjustedPath, 
                          globalTransformPath, leftSteroCorrectWorkDir, '', '', False, carry)


        rightStereoAdjustedPath  = os.path.join(tempFolder, rightStereoBaseName + '.stereoAdjusted.cub')
        rightSteroCorrectWorkDir = os.path.join(tempFolder, 'stereoRightStereoCorrection/')
        applyNavTransform(posOffsetCorrectedStereoRightPath, rightStereoAdjustedPath, 
                          globalTransformPath, rightSteroCorrectWorkDir, '', '', False, carry)

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(leftStereoAdjustedPath, rightStereoAdjustedPath, 
                                   os.path.join(tempFolder, 'stereoGlobalAdjustGdcCheck'), 
                                   expectedSurfaceElevation, carry)

#        meanMainError   = evaluateAccuracy(posOffsetCorrectedLeftPath,       posOffsetCorrectedRightPath, 
#                                           mainIpFindPath,       os.path.join(tempFolder, 'finalGdcCheck'))
#        meanStereoError = evaluateAccuracy(leftStereoAdjustedPath, rightStereoAdjustedPath, 
#                                           stereoIpFindPath,     os.path.join(tempFolder, 'finalGdcCheck'))
#        meanLeftError   = evaluateAccuracy(posOffsetCorrectedLeftPath,       leftStereoAdjustedPath,
#                                           pixelPairsLeftSmall,  os.path.join(tempFolder, 'finalGdcCheck'))
#        meanRightError  = evaluateAccuracy(posOffsetCorrectedRightPath,      rightStereoAdjustedPath, 
#                                           pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheck'))
 
#        print '-----> STEREO CORRECTED'
#        print '=====> Mean main   pair error = %.4f meters' % meanMainError
#        print '=====> Mean stereo pair error = %.4f meters' % meanStereoError
#        print '=====> Mean left   pair error = %.4f meters' % meanLeftError
#        print '=====> Mean right  pair error = %.4f meters' % meanRightError

        print '\n-------------------------------------------------------------------------\n'

        # Extract a large number of matching pixel locations (many thousands) from the LE/LE and RE/RE.
        # - The skip number is a row and column skip.
        pixelPairsLeftLarge = extractPixelPairsFromStereoResults(disparityImageLeft, tempFolder, 
                                                                 'stereoPixelPairsLeftLarge.csv', 8, carry)

        # TODO: Many changes needed before RE images can be used here!
        #pixelPairsRightLarge = extractPixelPairsFromStereoResults(disparityImageRight, tempFolder, 'stereoPixelPairsRightLarge.csv', 16, False)

        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        # - All this step does is use stereo intersection to determine a lat/lon/alt coordinate for each pixel pair in the large data set.  No optimization is performed.
        largeGdcFolder = os.path.join(tempFolder, 'gdcPointsLargeComp/')
        if not os.path.exists(largeGdcFolder):
            os.mkdir(largeGdcFolder)
        largeGdcPrefix = os.path.join(tempFolder, 'gdcPointsLargeComp/out')
        largeGdcFile   = largeGdcPrefix + '-initialGdcPoints.csv'
        if (not os.path.exists(largeGdcFile)) or carry:
            cmd = ('lronacAngleDoubleSolver --outputPrefix '           + largeGdcPrefix + 
                                          ' --matchingPixelsLeftPath ' + pixelPairsLeftLarge + 
                                          ' --leftCubePath '           + posOffsetCorrectedLeftPath + 
                                          ' --leftStereoCubePath '     + posOffsetCorrectedStereoLeftPath + 
                                          ' --initialOnly --initialValues ' + solvedParamsPath + 
                                          ' --elevation ' + str(expectedSurfaceElevation))
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        prepTime = time.time()
        logging.info('pc_align prep finished in %f seconds', prepTime - sbaTime)

        # TODO: Why does rotation always move the points somewhere else?
        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/dem')
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        #transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')
        pcAlignTransformPath  = pcAlignOutputPrefix + '-inverse-transform.txt'
        pcAlignFolder         = os.path.dirname(pcAlignOutputPrefix)
        if (not os.path.exists(pcAlignTransformPath)) or carry:
            # TODO: Confirm which input order works best
            #cmd = 'pc_align --highest-accuracy --max-displacement 1500 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            cmd = ('pc_align --highest-accuracy --max-displacement 200 --datum D_MOON ' + 
                   '--save-inv-transformed-reference-points ' + largeGdcFile + 
                   ' ' + options.lolaPath + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only')
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'

        if not os.path.exists(pcAlignTransformPath):
            raise Exception('pc_align call failed!')

        # Copy the pc_align log to the output folder
        pcAlignLogPath = findOutputLog(pcAlignFolder)
        shutil.copyfile(pcAlignLogPath, os.path.join(outputFolder, 'pcAlignLog.txt'))
        

        alignTime = time.time()
        logging.info('pc_align finished in %f seconds', alignTime - prepTime)


        print '\n-------------------------------------------------------------------------\n'

        # Now go back and apply the pc_align computed transform to all four cameras.
        # - This step corrects the four camera positions relative to the LOLA data.
        leftCkPath       = os.path.join(tempFolder, 'leftFinalCk.bc')
        leftSpkPath      = os.path.join(tempFolder, 'leftFinalSpk.bsp')
        leftFinalWorkDir = os.path.join(tempFolder, 'leftFullCorrection')
        applyNavTransform(posOffsetCorrectedLeftPath, outputPathLeft, 
                          pcAlignTransformPath, leftFinalWorkDir, leftCkPath, leftSpkPath, True, carry)

        partialCorrectedRightPath = os.path.join(tempFolder, rightBaseName + '.partial_corrected.cub')
        rightCkPath               = os.path.join(tempFolder, 'rightFinalCk.bc') 
        rightSpkPath              = os.path.join(tempFolder, 'rightFinalSpk.bsp')
        rightFinalWorkDir          = os.path.join(tempFolder, 'rightFullCorrection/')
        applyNavTransform(posOffsetCorrectedRightPath, partialCorrectedRightPath, 
                          pcAlignTransformPath, rightFinalWorkDir, rightCkPath, rightSpkPath, True, carry)


        leftStereoCkPath       = os.path.join(tempFolder, 'leftStereoFinalCk.bc') 
        leftStereoSpkPath      = os.path.join(tempFolder, 'leftStereoFinalSpk.bsp')
        leftStereoFinalWorkDir = os.path.join(tempFolder, 'leftStereoFullCorrection')
        applyNavTransform(leftStereoAdjustedPath, outputPathStereoLeft, 
                          pcAlignTransformPath, leftStereoFinalWorkDir, 
                          leftStereoCkPath, leftStereoSpkPath, True, carry)


        partialCorrectedStereoRightPath = os.path.join(tempFolder, rightStereoBaseName + '.partial_corrected.cub')
        rightStereoCkPath               = os.path.join(tempFolder, 'rightStereoFinalCk.bc')
        rightStereoSpkPath              = os.path.join(tempFolder, 'rightStereoFinalSpk.bsp')
        rightStereoFinalWorkDir         = os.path.join(tempFolder, 'rightStereoFullCorrection/')
        applyNavTransform(rightStereoAdjustedPath, partialCorrectedStereoRightPath, 
                          pcAlignTransformPath, rightStereoFinalWorkDir, 
                          rightStereoCkPath, rightStereoSpkPath, True, carry)

        # At this point the left images are hopefully in the correct position and we can apply the offset of the RE cameras


#        meanMainError   = evaluateAccuracy(outputPathLeft,       partialCorrectedRightPath, 
#                                           mainIpFindPath,       os.path.join(tempFolder, 'finalGdcCheck'))
#        meanStereoError = evaluateAccuracy(outputPathStereoLeft, partialCorrectedStereoRightPath, 
#                                           stereoIpFindPath,     os.path.join(tempFolder, 'finalGdcCheck'))
#        meanLeftError   = evaluateAccuracy(outputPathLeft,       outputPathStereoLeft,
#                                           pixelPairsLeftSmall,  os.path.join(tempFolder, 'finalGdcCheck'))
#        meanRightError  = evaluateAccuracy(partialCorrectedRightPath,      partialCorrectedStereoRightPath, 
#                                           pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheck'))
 

#        print '-----> ALL GLOBAL'
#        print '=====> Mean main   pair error = %.4f meters' % meanMainError
#        print '=====> Mean stereo pair error = %.4f meters' % meanStereoError
#        print '=====> Mean left   pair error = %.4f meters' % meanLeftError
#        print '=====> Mean right  pair error = %.4f meters' % meanRightError

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(outputPathStereoLeft, partialCorrectedStereoRightPath, 
                                   os.path.join(tempFolder, 'pcAlignStereoGdcCheck'), 
                                   expectedSurfaceElevation, carry)


        navTime = time.time()
        logging.info('Nav transforms finished in %f seconds', navTime - alignTime)


        print '\n-------------------------------------------------------------------------\n'

        # Apply local transforms to both pairs of images!

        # Apply the local rotation to the adjusted RE cube
        mainLocalWorkDir = os.path.join(tempFolder, 'mainLocalCorrection')
        applyInterCameraPairRotation(outputPathLeft, partialCorrectedRightPath, 
                                    localRotationPath, outputPathRight, 
                                    rightCkPath, rightSpkPath, mainLocalWorkDir, carry)

        # Apply the local rotation to the adjusted stereo RE cube
        stereoLocalWorkDir = os.path.join(tempFolder, 'stereoLocalCorrection')
        applyInterCameraPairRotation(outputPathStereoLeft, partialCorrectedStereoRightPath, 
                                     stereoRotationPath, outputPathStereoRight, 
                                     rightStereoCkPath, rightStereoSpkPath, stereoLocalWorkDir, carry)


        # DEBUG: Check angle solver on adjusted LE/RE images!
        checkAdjacentPairAlignment(outputPathLeft, outputPathRight, 
                                   os.path.join(tempFolder, 'finalGdcCheck'), 
                                   expectedSurfaceElevation, carry)

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(outputPathStereoLeft, outputPathStereoRight, 
                                   os.path.join(tempFolder, 'finalStereoGdcCheck'), 
                                   expectedSurfaceElevation, carry)

        localTime = time.time()
        logging.info('Local transforms finished in %f seconds', localTime - navTime)

        print '\n-------------------------------------------------------------------------\n'

        # One last pair of checks to compute accuracy
        # - This calls campt for pixel pairs and finds the GCC location difference
        meanMainError   = evaluateAccuracy(outputPathLeft,       outputPathRight, 
                                           mainIpFindPath,       os.path.join(tempFolder, 'finalGdcCheck'))
        meanStereoError = evaluateAccuracy(outputPathStereoLeft, outputPathStereoRight, 
                                           stereoIpFindPath,     os.path.join(tempFolder, 'finalGdcCheck'))
        meanLeftError   = evaluateAccuracy(outputPathLeft,       outputPathStereoLeft,
                                           pixelPairsLeftSmall,  os.path.join(tempFolder, 'finalGdcCheck'))
        meanRightError  = evaluateAccuracy(outputPathRight,      outputPathStereoRight, 
                                           pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheck'))
 
        print '=====> Mean main   pair error = %.4f meters' % meanMainError
        print '=====> Mean stereo pair error = %.4f meters' % meanStereoError
        print '=====> Mean left   pair error = %.4f meters' % meanLeftError
        print '=====> Mean right  pair error = %.4f meters' % meanRightError

        logging.info('=====> Mean main   pair error = %.4f meters' % meanMainError)
        logging.info('=====> Mean stereo pair error = %.4f meters' % meanStereoError)
        logging.info('=====> Mean left   pair error = %.4f meters' % meanLeftError)
        logging.info('=====> Mean right  pair error = %.4f meters' % meanRightError)


        # All finished!  We should have a fully calibrated version of each of the four input files.

        # Clean up temporary files
        if not options.keep:
            print 'Deleting temporary files'
            # Init files
            IsisTools.removeIfExists(spiceInitLeftPath)
            IsisTools.removeIfExists(spiceInitRightPath)
            IsisTools.removeIfExists(spiceInitStereoRightPath)
            IsisTools.removeIfExists(spiceInitStereoLeftPath)

            # Position correction files
            IsisTools.removeIfExists(posOffsetCorrectedLeftPath)
            IsisTools.removeIfExists(posOffsetCorrectedRightPath)
            IsisTools.removeIfExists(posOffsetCorrectedStereoLeftPath)
            IsisTools.removeIfExists(posOffsetCorrectedStereoRightPath)
            IsisTools.removeFolderIfExists(leftPosCorrectWorkDir)
            IsisTools.removeFolderIfExists(rightPosCorrectWorkDir)
            IsisTools.removeFolderIfExists(stereoLeftPosCorrectWorkDir)
            IsisTools.removeFolderIfExists(stereoRightPosCorrectWorkDir)
            # Stereo output
            stereoOutputLeftFolder = os.path.dirname(stereoPrefixLeft)
            IsisTools.removeFolderIfExists(stereoOutputLeftFolder) #TURN THIS ON AFTER TESTING COMPLETE!
            IsisTools.removeIfExists(pixelPairsLeftSmall)
            IsisTools.removeIfExists(pixelPairsRightSmall)
            IsisTools.removeIfExists(pixelPairsLeftCrossSmall)
            IsisTools.removeIfExists(pixelPairsRightCrossSmall)
            IsisTools.removeIfExists(leftPosCorrectedCropped)
            IsisTools.removeIfExists(rightPosCorrectedCropped)
            IsisTools.removeIfExists(leftStereoPosCorrectedCropped)
            IsisTools.removeIfExists(rightStereoPosCorrectedCropped)

            ## Remove all the SBA files
            #fileList = [ f for f in os.listdir(tempFolder) if f.startswith("SBA_solution") ]
            #for f in fileList:
            #    IsisTools.removeIfExists(os.path.join(tempFolder, f))

            # Remove stereo corrected files
            IsisTools.removeIfExists(leftStereoAdjustedPath)
            IsisTools.removeFolderIfExists(leftSteroCorrectWorkDir)
            IsisTools.removeIfExists(rightStereoAdjustedPath)
            IsisTools.removeFolderIfExists(rightSteroCorrectWorkDir)
            # Clean pc_align steps
            IsisTools.removeIfExists(pixelPairsLeftLarge)
            #IsisTools.removeFolderIfExists(largeGdcFolder)
            #IsisTools.removeFolderIfExists(pcAlignFolder)
            # Clean out final file correction
            IsisTools.removeIfExists(leftCkPath)
            IsisTools.removeIfExists(leftSpkPath)
            IsisTools.removeFolderIfExists(leftFinalWorkDir)
            IsisTools.removeIfExists(partialCorrectedRightPath)
            IsisTools.removeIfExists(rightCkPath)
            IsisTools.removeIfExists(rightSpkPath)
            IsisTools.removeFolderIfExists(rightFinalWorkDir)
            IsisTools.removeIfExists(leftStereoCkPath)
            IsisTools.removeIfExists(leftStereoSpkPath)
            IsisTools.removeFolderIfExists(leftStereoFinalWorkDir)
            IsisTools.removeIfExists(partialCorrectedStereoRightPath)
            IsisTools.removeIfExists(rightStereoCkPath)
            IsisTools.removeIfExists(rightStereoSpkPath)
            IsisTools.removeFolderIfExists(rightStereoFinalWorkDir)
            # Remove local transform folders
            IsisTools.removeFolderIfExists(mainLocalWorkDir)
            IsisTools.removeFolderIfExists(stereoLocalWorkDir)
            ## Remove check folders
            #IsisTools.removeFolderIfExists(os.path.join(tempFolder, 'pcAlignStereoGdcCheck'))
            #IsisTools.removeFolderIfExists(os.path.join(tempFolder, 'stereoGlobalAdjustGdcCheck'))
            #IsisTools.removeFolderIfExists(os.path.join(tempFolder, 'finalGdcCheck'))
            #IsisTools.removeFolderIfExists(os.path.join(tempFolder, 'finalStereoGdcCheck'))

            #if (hadToCreateTempFolder):
            #    IsisTools.removeFolderIfExists(tempFolder)

        endTime = time.time()

        print "Finished stereo calibration process in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
