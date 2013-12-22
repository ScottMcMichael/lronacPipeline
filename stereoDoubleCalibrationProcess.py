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

import os, glob, optparse, re, shutil, subprocess, string, time, math

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

#===============================================================


# TODO: Make this a standalone function!
# Get a single file ready to process
def prepareImgFile(inputPath, outputFolder):

    # Generate the appropriate paths in the output folder
    temp           = os.path.splitext(inputPath)[0] + '.cub'
    initialCubFile = os.path.join(outputFolder, os.path.basename(temp))
    temp           = os.path.splitext(inputPath)[0] + '.lronaccal.cub'
    calCubFile     = os.path.join(outputFolder, os.path.basename(temp))
    temp           = os.path.splitext(inputPath)[0] + '.lronaccal.lronacecho.cub'
    echoCubFile    = os.path.join(outputFolder, os.path.basename(temp))

    # Quit immediately if the output file already exists
    if os.path.exists(echoCubFile):
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


    return echoCubFile




# Apply the position offset specified in the IK kernel to either an LE or RE camera.
def applyInterCameraPositionOffset(inputCubePath, outputCubePath, workingDirectory, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping position offset correction.'
        return True

    # Run the process
    cmd = 'positionCorrector.py --keep --input ' + inputCubePath + ' --output ' + outputCubePath + ' --workDir ' + workingDirectory
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Position offset correction failed to create output file ' + outputCubePath + ' from input file ' + inputCubePath)

    return True



# Samples pixels pairs for bundle adjustment from the output of the stereo command.
# - Returns the output path.
def extractPixelPairsFromStereoResults(disparityImagePath, outputDirectory, outputFileName, sampleInterval, forceOperation):

    # Quit immediately if the output file already exists
    outputPixelPath = os.path.join(outputDirectory, outputFileName)
    if (not forceOperation) and (os.path.exists(outputPixelPath)):
        print 'File ' + outputPixelPath + ' already exists, skipping pixel pair sampling.'
        return outputPixelPath

    # Run the process
    cmd = 'pixelPairsFromStereo -i ' + disparityImagePath + ' -o ' + outputPixelPath + ' -p ' + str(sampleInterval)
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputPixelPath):
        raise Exception('Pixel sampling failed to create output file ' + outputPixelPath + ' from input file ' + disparityImagePath)

    return outputPixelPath




# Applies a planet-centered rotation and correction to the nav data of a cube
# - If the CK and SPK paths are not specified, paths are automatically generated.
# - Set forceOperation to run the operation even if the output file already exists
def applyNavTransform(inputCubePath, outputCubePath, transformMatrixPath, workDir, ckPath, spkPath, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
        return True

    # Set up CK and SPK manual paths if they were specified
    ckLine  = ''
    spkLine = ''
    if (ckPath):
        ckLine  = ' --ck '  + ckPath
    if (spkPath):
        spkLine = ' --spk ' + spkPath

    # Execute the transformation command
    cmd = 'rotationCorrector.py --keep --input ' + inputCubePath + ' --output ' + outputCubePath + ' --transformPath ' + transformMatrixPath + ' --workDir ' + workDir + ckLine + spkLine
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Nav transform failed to create output file ' + outputCubePath + ' from input file ' + inputCubePath)

    return True




# Tries to compute the internal angle between an LE/RE image pair.
# - Output GDC points serve as a check to make sure the images are in roughly the correct place.
def checkAdjacentPairAlignment(leftInputPath, rightInputPath, outputDirectory, forceOperation):

    # Figure out output paths
    if not os.path.exists(outputDirectory):
        os.mkdir(outputDirectory)
    sbaOutputPrefix = os.path.join(outputDirectory, 'SBA_check')
    defaultGdcPath  = sbaOutputPrefix + '-outputGdcPoints.csv'

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(defaultGdcPath)):
        return True

    # Run the process
    #cmd = 'lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + outputGdcPath + ' ' + leftInputPath + ' ' + rightInputPath
    cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefix + ' --leftCubePath ' + leftInputPath + ' --rightCubePath ' + rightInputPath
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(defaultGdcPath):
        raise Exception('Adjacency check failed to create output file ' + defaultGdcPath + ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    return True




# Generate a modified IK kernel to adjust the rotation between an LE/RE camera pair.
def applyInterCameraPairRotation(leftInputPath, rightInputPath, newRotationPath, outputCubePath, ckPath, spkPath, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
        return True

    # Generate the new file
    cmd = 'lronacCameraRotationCorrector.py --keep --output ' + outputCubePath + ' --rotation ' + newRotationPath + ' --left ' + leftInputPath + ' --right ' + rightInputPath + ' --ck ' + ckPath + ' --spk ' + spkPath 
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Inter-camera rotation failed to create output file ' + outputCubePath + ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    return True



# Calls stereo functions to generate a disparity image and returns the path to it.
def callStereoCorrelation(leftInputPath, rightInputPath, outputPrefix, correlationTimeout, forceOperation):

    # Quit immediately if the output file already exists
    disparityImagePath = outputPrefix + '-D.tif'
    if (not forceOperation) and (os.path.exists(disparityImagePath)):
        print 'File ' + disparityImagePath + ' already exists, skipping stereo computation.'
        return disparityImagePath

    # Use parallel stereo call, steps 0 and 1 only.  Other options to try and increase speed.
    cmd = 'parallel_stereo  --corr-max-levels 3 --compute-error-vector --entry-point 0 --stop-point 2 --alignment none --subpixel-mode 1 --disable-fill-holes --processes 8 --threads-multiprocess 4 --threads-singleprocess 32 --cost-mode 0 --corr-timeout ' + str(correlationTimeout) + ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix
    # TODO: Alignment methods mess up the disparity numbers without some further processing!
#    cmd = 'parallel_stereo --compute-error-vector --entry-point 0 --stop-point 2 --alignment affineepipolar --subpixel-mode 1 --disable-fill-holes --processes 8 --threads-multiprocess 4 --threads-singleprocess 32 --cost-mode 0 --corr-timeout ' + str(correlationTimeout) + ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix
    print cmd
    os.system(cmd)

    # Generate disparity images for debugging
    #cmd = 'disparitydebug ' + disparityImagePath
    #print cmd
    #os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(disparityImagePath):
        raise Exception('Stereo processing failed to create output file ' + disparityImagePath + ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    return disparityImagePath

# Counts the number of lines in a file
def getFileLineCount(filePath):
    f = open(filePath)
    i = 0
    for line in f:
        i = i + 1
    return i


# Compares the backprojected locations of matched pixel pairs
def evaluateAccuracy(leftCubePath, rightCubePath, ipFindOutputPath, workDir=''):

    # Run check every N lines in the file
    lineSkip = 25
    
    f = open(ipFindOutputPath)
    i = 0
    sumDistance = 0.0
    count       = 0.0
    for line in f:
        i = i + 1
        if i % lineSkip == 0: # Compare locations on a sampling of lines
            i = 0
            # Get the matching pixel coordinates from this line
            leftSample, leftLine, rightSample, rightLine = line.split(',')
            
            # Obtain the backprojected GCC location for the pixels in the two images
            leftLoc  = IsisTools.getPixelLocInCube(leftCubePath,  leftSample,  leftLine,  workDir)
            rightLoc = IsisTools.getPixelLocInCube(rightCubePath, rightSample, rightLine, workDir)
            
            # Determine the distance between the GCC points and accumulate
            distance = math.sqrt(math.pow(leftLoc[0] - rightLoc[0], 2.0) + math.pow(leftLoc[1] - rightLoc[1], 2.0) + math.pow(leftLoc[2] - rightLoc[2], 2.0))
#            print 'Distance = %.2f' % distance
            sumDistance = sumDistance + distance
            count = count + 1.0
    
    # Determine the mean distance
    meanDistance = sumDistance / count
    return meanDistance


#==========================================================================================

def main():

    print "Started stereoDoubleCalibrationProcess.py"

    try:
        try:
            usage = "usage: stereoDoubleCalibrationProcess.py TODO [--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--left",  dest="leftPath",  help="Path to LE .IMG file")
            parser.add_option("--right", dest="rightPath", help="Path to RE .IMG file")
            
            parser.add_option("--stereo-left",  dest="stereoLeft",  help="Path to LE .IMG file with overlapping view of --left file")
            parser.add_option("--stereo-right", dest="stereoRight", help="Path to RE .IMG file with overlapping view of --right file")
  
            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")
            parser.add_option("--lola",    dest="lolaPath", help="Path to LOLA DEM")

            parser.add_option("--outputL",  dest="outputPathLeft",        help="Where to write the output LE file.")
            parser.add_option("--outputR",  dest="outputPathRight",       help="Where to write the output RE file.")
            parser.add_option("--outputSL", dest="outputPathStereoLeft",  help="Where to write the output Stereo LE file.")
            parser.add_option("--outputSR", dest="outputPathStereoRight", help="Where to write the output Stereo RE file.")

            
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
                parser.error("Need LOLA DEM path")
            if not options.outputPathLeft: 
                parser.error("Need left output path")
            if not options.outputPathRight: 
                parser.error("Need right output path")
            if not options.outputPathStereoLeft: 
                parser.error("Need stereo left output path")
            if not options.outputPathStereoRight: 
                parser.error("Need stereo right output path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder  = os.path.dirname(options.outputPathLeft)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder) 
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder)

        # Convert the input files from IMG files to spiceinit'ed cubes in the output folder
        spiceInitLeftPath        = prepareImgFile(options.leftPath,    tempFolder)
        spiceInitRightPath       = prepareImgFile(options.rightPath,   tempFolder)
        spiceInitStereoLeftPath  = prepareImgFile(options.stereoLeft,  tempFolder)
        spiceInitStereoRightPath = prepareImgFile(options.stereoRight, tempFolder)

        # DEBUG: Check angle solver on input LE/RE images!
        checkAdjacentPairAlignment(spiceInitLeftPath, spiceInitRightPath, os.path.join(tempFolder, 'initialGdcCheck'), False)

        # Apply LE/RE LRONAC position offsets to each of the input files
        posOffsetCorrectedLeftPath = os.path.join(tempFolder, 'left.posOffsetCorrected.cub')
        thisWorkDir                = os.path.join(tempFolder, 'leftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitLeftPath, posOffsetCorrectedLeftPath, thisWorkDir, False)

        posOffsetCorrectedRightPath = os.path.join(tempFolder, 'right.posOffsetCorrected.cub')
        thisWorkDir                 = os.path.join(tempFolder, 'rightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitRightPath, posOffsetCorrectedRightPath, thisWorkDir, False)

        posOffsetCorrectedStereoLeftPath = os.path.join(tempFolder, 'stereoLeft.posOffsetCorrected.cub')
        thisWorkDir                      = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoLeftPath, posOffsetCorrectedStereoLeftPath, thisWorkDir, False)

        posOffsetCorrectedStereoRightPath = os.path.join(tempFolder, 'stereoRight.posOffsetCorrected.cub')
        thisWorkDir                       = os.path.join(tempFolder, 'stereoRightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoRightPath, posOffsetCorrectedStereoRightPath,  thisWorkDir, False)

        # DEBUG: Check angle solver on input LE/RE images!
        checkAdjacentPairAlignment(posOffsetCorrectedLeftPath, posOffsetCorrectedRightPath, os.path.join(tempFolder, 'posCorrectGdcCheck'), False)
        checkAdjacentPairAlignment(posOffsetCorrectedStereoLeftPath, posOffsetCorrectedStereoRightPath, os.path.join(tempFolder, 'posCorrectStereoGdcCheck'), False)

        print '\n-------------------------------------------------------------------------\n'

        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefixLeft   = os.path.join(tempFolder, 'stereoOutputLeft/out')
        disparityImageLeft = callStereoCorrelation(posOffsetCorrectedLeftPath, posOffsetCorrectedStereoLeftPath, stereoPrefixLeft, 400, False)

        #raise Exception('Done running left stereo!')

        # Perform initial stereo step on two RE cubes to generate a large number of point correspondences
        stereoPrefixRight   = os.path.join(tempFolder, 'stereoOutputRight/out')
        disparityImageRight = callStereoCorrelation(posOffsetCorrectedRightPath, posOffsetCorrectedStereoRightPath, stereoPrefixRight, 400, False)


        # Extract a small number of matching pixel locations from the LE and RE disparity images ( < 300 pairs)
        pixelPairsLeftSmall = extractPixelPairsFromStereoResults(disparityImageLeft, tempFolder, 'stereoPixelPairsLeftSmall.csv', 800, False)
        pixelPairsRightSmall = extractPixelPairsFromStereoResults(disparityImageRight, tempFolder, 'stereoPixelPairsRightSmall.csv', 800, False)

        print '\n-------------------------------------------------------------------------\n'

        # Perform cross-stereo matching of LE/RE cubes from opposite pairs

        # Crop each input cube to half its width
        # - Take the right half of the LE cubes and the left half of the RE cubes
        leftPosCorrectedCropped        = os.path.join(tempFolder, 'leftCropped.cub')
        rightPosCorrectedCropped       = os.path.join(tempFolder, 'rightCropped.cub')
        leftStereoPosCorrectedCropped  = os.path.join(tempFolder, 'leftStereoCropped.cub')
        rightStereoPosCorrectedCropped = os.path.join(tempFolder, 'rightStereoCropped.cub')

        if not os.path.exists(leftPosCorrectedCropped):
            cmd = 'crop sample=2531 from=' + posOffsetCorrectedLeftPath + ' to= ' + leftPosCorrectedCropped
            print cmd
            os.system(cmd)
        if not os.path.exists(rightPosCorrectedCropped):
            cmd = 'crop nsamples=2532 from=' + posOffsetCorrectedRightPath + ' to= ' + rightPosCorrectedCropped
            print cmd
            os.system(cmd)
        if not os.path.exists(leftStereoPosCorrectedCropped):
            cmd = 'crop sample=2531 from=' + posOffsetCorrectedStereoLeftPath + ' to= ' + leftStereoPosCorrectedCropped
            print cmd
            os.system(cmd)
        if not os.path.exists(rightStereoPosCorrectedCropped):
            cmd = 'crop nsamples=2532 from=' + posOffsetCorrectedStereoRightPath + ' to= ' + rightStereoPosCorrectedCropped
            print cmd
            os.system(cmd)


        # Use stereo command on two cross-pair cubes
        # - Timeouts are shorter on these since there is probably less image overlap
        
        # First is left in main pair to right in the stereo pair
        stereoPrefixLeftCross   = os.path.join(tempFolder, 'stereoOutputLeftCross/out')
        try:
            disparityImageLeftCross = callStereoCorrelation(leftPosCorrectedCropped, rightStereoPosCorrectedCropped, stereoPrefixLeftCross, 100, False)
            usingLeftCross = True
        except:
            print 'Failed to find left-cross match, ignoring this data source.'
            usingLeftCross = False

        # Next is left in the stereo pair to right in the main pair
        stereoPrefixRightCross   = os.path.join(tempFolder, 'stereoOutputRightCross/out')
        try:
            disparityImageRightCross = callStereoCorrelation(leftStereoPosCorrectedCropped, rightPosCorrectedCropped, stereoPrefixRightCross, 100, False)
            usingRightCross = True
        except:
            print 'Failed to find right-cross match, ignoring this data source.'
            usingRightCross = False

        # Extract a small number of matching pixel locations from the disparity images ( < 300 pairs)
        # - The pixels are extracted more densely because there is much less overlap area to work with.
        if usingLeftCross:
            pixelPairsLeftCrossSmall = extractPixelPairsFromStereoResults(disparityImageLeftCross, tempFolder, 'stereoPixelPairsLeftCrossSmall.csv', 400, False)
            numLeftCrossPairs = getFileLineCount(pixelPairsLeftCrossSmall)
        if usingRightCross:
            pixelPairsRightCrossSmall = extractPixelPairsFromStereoResults(disparityImageRightCross, tempFolder, 'stereoPixelPairsRightCrossSmall.csv', 400, False)
            numRightCrossPairs = getFileLineCount(pixelPairsRightCrossSmall)

        # Left and right cross pixels may not be used depending on the image overlap
        # - If a small number of pixel pairs was found that suggests we may have a high ratio of bad pixels
        leftCrossString  = ''
        rightCrossString = ''
        if usingLeftCross and (numLeftCrossPairs > 100):
            leftCrossString  = ' --matchingPixelsLeftCrossPath '  + pixelPairsLeftCrossSmall
        if usingRightCross and (numRightCrossPairs > 100):
            rightCrossString = ' --matchingPixelsRightCrossPath ' + pixelPairsRightCrossSmall

        print '\n-------------------------------------------------------------------------\n'

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
        if not os.path.exists(globalTransformPath):
            cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefix + ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --matchingPixelsRightPath ' + pixelPairsRightSmall + leftCrossString + rightCrossString + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --rightCubePath ' + posOffsetCorrectedRightPath + ' --leftStereoCubePath ' + posOffsetCorrectedStereoLeftPath + ' --rightStereoCubePath ' + posOffsetCorrectedStereoRightPath
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo transform calculation step'

#        raise Exception('buggin out')

#        # DEBUG - Confirm output results are roughly the same
        zeroParamsPath           = '/byss/moon/lronacPipeline_V2/zeroParamsFile.csv'
        scPosCube     = os.path.join(tempFolder, 'stereoLeft.sc.cub')
        sbaOutputPrefixDebug     = os.path.join(tempFolder, 'SBA_solutionDEBUG3')
#        posOffsetCorrectedStereoLeftPath
        cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefixDebug+ ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --leftStereoCubePath ' + scPosCube + " --initialOnly --initialValues " + zeroParamsPath 
#        print cmd
#        os.system(cmd)   
#        raise Exception('buggin out')

#        # DEBUG - Confirm output results are roughly the same
        sbaOutputPrefixDebug     = os.path.join(tempFolder, 'SBA_solutionDEBUG')
        cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefixDebug+ ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --leftStereoCubePath ' + posOffsetCorrectedStereoLeftPath + " --initialOnly --initialValues " + solvedParamsPath 
#        print cmd
#        os.system(cmd)   
#        raise Exception('buggin out')
        
        
        # Apply the planet-centered rotation/translation to both cameras in the stereo pair.
        # - This corrects the stereo pair relative to the main pair.
        # - The RE relative to LE corrections are performed later for convenience.
#        zeroRotParamsPath           = os.path.join(tempFolder, 'zeroRotGlobalTransformMatrix.csv')
        leftStereoAdjustedPath = os.path.join(tempFolder, 'leftStereoAdjusted.cub')
        thisWorkDir            = os.path.join(tempFolder, 'stereoLeftStereoCorrection/')
        applyNavTransform(posOffsetCorrectedStereoLeftPath, leftStereoAdjustedPath, globalTransformPath, thisWorkDir, '', '', False)
#        raise Exception('buggin out')
        
#        # DEBUG - Outputs here should be identical to the ones from the previous DEBUG call!
        sbaOutputPrefixDebug     = os.path.join(tempFolder, 'SBA_solutionDEBUG2')
        cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefixDebug+ ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --leftStereoCubePath ' + leftStereoAdjustedPath + " --initialOnly --initialValues " + zeroParamsPath 
#        print cmd
#        os.system(cmd)
        
#        raise Exception('buggin out')
        
        rightStereoAdjustedPath = os.path.join(tempFolder, 'rightStereoAdjusted.cub')
        thisWorkDir             = os.path.join(tempFolder, 'stereoRightStereoCorrection/')
        applyNavTransform(posOffsetCorrectedStereoRightPath, rightStereoAdjustedPath, globalTransformPath, thisWorkDir, '', '', False)

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(leftStereoAdjustedPath, rightStereoAdjustedPath, os.path.join(tempFolder, 'stereoGlobalAdjustGdcCheck'), False)


        # DEBUG: Re-run the SBA solver with the global adjustment applied to the stereo images.
        # - Since we just applied the solved for transform, we expect global transform parameters to be near zero.
#        TODO: Nail down why this is not happening!!!
        sbaGlobalCheckOutputPrefix   = os.path.join(tempFolder, 'globalSbaCheck/SBA_solution')
        globalSbaCheckTransformPath = sbaGlobalCheckOutputPrefix + "-globalTransformMatrix.csv"
        if not os.path.exists(globalSbaCheckTransformPath):
            cmd = 'lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefix + ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --matchingPixelsRightPath ' + pixelPairsRightSmall + leftCrossString + rightCrossString + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --rightCubePath ' + posOffsetCorrectedRightPath + ' --leftStereoCubePath ' + leftStereoAdjustedPath + ' --rightStereoCubePath ' + rightStereoAdjustedPath
            print cmd
            os.system(cmd)

        #raise Exception('done with SBA global check')

        print '\n-------------------------------------------------------------------------\n'

        # Extract a large number of matching pixel locations (many thousands) from the LE/LE and RE/RE.
        # - The skip number is a row and column skip.
        pixelPairsLeftLarge = extractPixelPairsFromStereoResults(disparityImageLeft, tempFolder, 'stereoPixelPairsLeftLarge.csv', 20, False)

        # TODO: Many changes needed before RE images can be used here!
        #pixelPairsRightLarge = extractPixelPairsFromStereoResults(disparityImageRight, tempFolder, 'stereoPixelPairsRightLarge.csv', 16, False)

        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        # - All this step does is use stereo intersection to determine a lat/lon/alt coordinate for each pixel pair in the large data set.  No optimization is performed.
        largeGdcFolder = os.path.join(tempFolder, 'gdcPointsLargeComp/')
        if not os.path.exists(largeGdcFolder):
            os.mkdir(largeGdcFolder)
        largeGdcPrefix = os.path.join(tempFolder, 'gdcPointsLargeComp/out')
        largeGdcFile   = largeGdcPrefix + '-initialGdcPoints.csv'
        if not os.path.exists(largeGdcFile):
            cmd = 'lronacAngleDoubleSolver --outputPrefix ' + largeGdcPrefix + ' --matchingPixelsLeftPath ' + pixelPairsLeftLarge + ' --leftCubePath ' + posOffsetCorrectedLeftPath + ' --leftStereoCubePath ' + posOffsetCorrectedStereoLeftPath + " --initialOnly --initialValues " + solvedParamsPath 
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        # TODO: Why does rotation always move the points somewhere else?
        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/dem')
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        #transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')
        pcAlignTransformPath  = pcAlignOutputPrefix + '-inverse-transform.txt'
        if not os.path.exists(pcAlignTransformPath):
            # TODO: Confirm which input order works best
            #cmd = 'pc_align --highest-accuracy --max-displacement 1500 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            cmd = 'pc_align --highest-accuracy --max-displacement 1500 --datum D_MOON --save-inv-transformed-reference-points ' + largeGdcFile + ' ' + options.lolaPath + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'
            
        if not os.path.exists(pcAlignTransformPath):
            raise Exception('pc_align call failed!')

#        raise Exception('pc align rot test')

        print '\n-------------------------------------------------------------------------\n'

        # Now go back and apply the pc_align computed transform to all four cameras.
        # - This step corrects the four camera positions relative to the LOLA data.
        leftCkPath  = os.path.join(tempFolder, 'leftFinalCk.bc')
        leftSpkPath = os.path.join(tempFolder, 'leftFinalSpk.bsp')
        thisWorkDir = os.path.join(tempFolder, 'leftFullCorrection')
        applyNavTransform(posOffsetCorrectedLeftPath, options.outputPathLeft, pcAlignTransformPath, thisWorkDir, leftCkPath, leftSpkPath, False)

        partialCorrectedRightPath = os.path.join(tempFolder, 'partial_corrected_RE.cub')
        rightCkPath               = os.path.join(tempFolder, 'rightFinalCk.bc') 
        rightSpkPath              = os.path.join(tempFolder, 'rightFinalSpk.bsp')
        thisWorkDir               = os.path.join(tempFolder, 'rightFullCorrection/')
        applyNavTransform(posOffsetCorrectedRightPath, partialCorrectedRightPath, pcAlignTransformPath, thisWorkDir, rightCkPath, rightSpkPath, False)


        leftStereoCkPath  = os.path.join(tempFolder, 'leftStereoFinalCk.bc') 
        leftStereoSpkPath = os.path.join(tempFolder, 'leftStereoFinalSpk.bsp')
        thisWorkDir       = os.path.join(tempFolder, 'leftStereoFullCorrection')
        applyNavTransform(leftStereoAdjustedPath, options.outputPathStereoLeft, pcAlignTransformPath, thisWorkDir, leftStereoCkPath, leftStereoSpkPath, False)

        partialCorrectedStereoRightPath = os.path.join(tempFolder, 'partial_corrected_stereo_RE.cub')
        rightStereoCkPath               = os.path.join(tempFolder, 'rightStereoFinalCk.bc')
        rightStereoSpkPath              = os.path.join(tempFolder, 'rightStereoFinalSpk.bsp')
        thisWorkDir                     = os.path.join(tempFolder, 'rightStereoFullCorrection/')
        applyNavTransform(rightStereoAdjustedPath, partialCorrectedStereoRightPath, pcAlignTransformPath, thisWorkDir, rightStereoCkPath, rightStereoSpkPath, False)

        # At this point the left images are hopefully in the correct position and we can apply the offset of the RE cameras

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(options.outputPathStereoLeft, partialCorrectedStereoRightPath, os.path.join(tempFolder, 'pcAlignStereoGdcCheck'), False)


        print '\n-------------------------------------------------------------------------\n'

        # Apply local transforms to both pairs of images!

        # Apply the local rotation to the adjusted RE cube
        applyInterCameraPairRotation(options.outputPathLeft, partialCorrectedRightPath, localRotationPath, options.outputPathRight, rightCkPath, rightSpkPath, False)

        # Apply the local rotation to the adjusted stereo RE cube
        applyInterCameraPairRotation(options.outputPathStereoLeft, partialCorrectedStereoRightPath, stereoRotationPath, options.outputPathStereoRight, rightStereoCkPath, rightStereoSpkPath, False)

        # DEBUG: Check angle solver on adjusted LE/RE images!
        checkAdjacentPairAlignment(options.outputPathLeft, options.outputPathRight, os.path.join(tempFolder, 'finalGdcCheck'), False)

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        checkAdjacentPairAlignment(options.outputPathStereoLeft, options.outputPathStereoRight, os.path.join(tempFolder, 'finalStereoGdcCheck'), False)
        
        print '\n-------------------------------------------------------------------------\n'

        # One last pair of checks to compute accuracy
        # - This calls campt for pixel pairs and finds the GCC location difference
        meanMainError = evaluateAccuracy(options.outputPathLeft, options.outputPathRight, mainIpFindPath, os.path.join(tempFolder, 'finalGdcCheck'))
        meanStereoError = evaluateAccuracy(options.outputPathStereoLeft, options.outputPathStereoRight, stereoIpFindPath, os.path.join(tempFolder, 'finalGdcCheck'))
        meanLeftError = evaluateAccuracy(options.outputPathLeft, options.outputPathStereoLeft, pixelPairsLeftSmall, os.path.join(tempFolder, 'finalGdcCheck'))
        meanRightError = evaluateAccuracy(options.outputPathRight, options.outputPathStereoRight, pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheck'))
        
        print '=====> Mean main   pair error = %.4f meters' % meanMainError
        print '=====> Mean stereo pair error = %.4f meters' % meanStereoError
        print '=====> Mean left   pair error = %.4f meters' % meanLeftError
        print '=====> Mean right  pair error = %.4f meters' % meanRightError

        # All finished!  We should have a fully calibrated version of each of the four input files.

        # Clean up temporary files
#        if not options.keep:
#            os.remove(tempTextPath)

        # --- Usability ---
        # TODO: Check for half-size pairs, handle accordingly
        # TODO: Clean up files and names

        # --- Accuracy improvement ---
        # TODO: Get Oleg's changes to look at error numbers
        # TODO: Experiment with more advanced mosaic merging tools

        endTime = time.time()

        print "Finished stereo calibration process in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
