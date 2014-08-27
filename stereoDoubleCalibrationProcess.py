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

import os, glob, optparse, re, shutil, subprocess, string, time, math, logging, threading, numpy

import IrgFileFunctions, IrgIsisFunctions, IrgAspFunctions

import positionCorrectorV2, rotationCorrectorV2

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
        IrgFileFunctions.removeIfExists(initialCubFile)
        IrgFileFunctions.removeIfExists(calCubFile)

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
    #print cmdArgs
    positionCorrectorV2.main(cmdArgs)


    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Position offset correction failed to create output file ' + outputCubePath + 
                        ' from input file ' + inputCubePath)

    return True


def extractAffinesFromStereoPprcLog(logPath):
    """Returns the affine transforms used in a stereo_pprc call in text format"""

    # Make sure the log file exists
    if not os.path.exists(logPath):
        raise Exception('Error: stereo_pprc log file ' + logPath + ' not found!')
    
    leftAffine  = '' # Default is no transforms used
    rightAffine = ''
    
    f = open(logPath, 'r') # Open the log file
    count = -1
    numberSets = [] # Need to read four sets of numbers
    for line in f:
            
        if (count > 0): # Reading the transform (may spread over variable number of lines)
            parenGroups = re.findall( r'\([ \de.,-]*\)', line) # Find all number sets on this line

            for p in parenGroups: # Record number sets as they are found
                count = count - 1
                numberSets.append(p[1:-1]) # Remove ()
 
            if count <= 0: # Found all four sets, pack up and quit
                leftAffine  = numberSets[0] + ', ' + numberSets[1]
                rightAffine = numberSets[2] + ', ' + numberSets[3]
                return (leftAffine, rightAffine)

            if not parenGroups:
            	continue # Did not find any transform elements on this line

        elif 'Aligning left and right images using affine matrices' in line:
            count = 4 # Need to read in four values

    raise Exception('Error: Failed to extract left and right affine transforms!')

# Samples pixels pairs for bundle adjustment from the output of the stereo command.
# - Returns the output path.
def extractPixelPairsFromStereoResults(disparityImagePath, outputPixelPath, 
                                       sampleInterval, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputPixelPath)):
        print 'File ' + outputPixelPath + ' already exists, skipping pixel pair sampling.'
        return outputPixelPath

    # TODO: Move this out to a function?
    # Determine the affine transforms which were used (stored in the log file)
    leftAffine  = '1 0 0 0 1 0 0 0 1' # Default both transforms to identity.
    rightAffine = '1 0 0 0 1 0 0 0 1'
    stereoFolder = os.path.dirname(disparityImagePath)
    stereoFiles  = os.listdir(stereoFolder)
    #print str(stereoFiles)
    logFile = None
    for f in stereoFiles: 
        if 'log-stereo_pprc' in f: 
            if f > logFile: # Get the log file with the highest number
                logFile = f
    #print logFile
    if (logFile):
        logFile     = os.path.join(stereoFolder, logFile)
        leftAffine, rightAffine = extractAffinesFromStereoPprcLog(logFile)
        
        leftAffine  = leftAffine.replace (',', ' ') # Strip out the commas
        rightAffine = rightAffine.replace(',', ' ')

    # Run the process
    cmd = ('pixelPairsFromStereo ' + disparityImagePath + ' ' + outputPixelPath + 
          ' --pointSpacing ' + str(sampleInterval) + ' --leftAffine ' + leftAffine + ' --rightAffine ' + rightAffine)
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputPixelPath):
        raise Exception('Pixel sampling failed to create output file ' + outputPixelPath + 
                        ' from input file ' + disparityImagePath)

    numPairs = IrgFileFunctions.getFileLineCount(outputPixelPath)
    logging.info('In stereo file %s found %d point pairs', disparityImagePath, numPairs)

    return numPairs


# Applies a planet-centered rotation and correction to the nav data of a cube
# - Set forceOperation to run the operation even if the output file already exists
#def applyNavTransform(inputCubePath, outputCubePath, transformMatrixPath, workDir, 
#                      ckPath, spkPath, pcAlignTrans, forceOperation):
def applyNavTransform(inputCubePath, outputCubePath, transformMatrixPath,
                      pcAlignTrans, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCubePath)):
        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
        return True

    
    # Set up the transformation command
    cmdArgs = ['--input', inputCubePath, '--output', outputCubePath, 
               '--transformPath', transformMatrixPath]
    
    # Set flag if the transform comes from pc_align as opposed to the lronacAngleSolver
    if (pcAlignTrans):
        cmdArgs.append('--pcAlignTrans')
    
    # Execute the transformation command
    #print cmdArgs
    rotationCorrectorV2.main(cmdArgs)

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

    # Run the process silently
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
    cmd = ('stereo_pprc --corr-max-levels 3 --compute-error-vector --cost-mode 0 --alignment AffineEpipolar' + 
                      ' --subpixel-mode 1 --disable-fill-holes --corr-timeout ' + str(correlationTimeout) + 
                      ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix)
    print cmd
    os.system(cmd)
    # Stage 1 (slow!)
    cmd = ('stereo_corr --corr-max-levels 3 --compute-error-vector --cost-mode 0 --alignment AffineEpipolar' + 
                      ' --subpixel-mode 1 --disable-fill-holes --corr-timeout ' + str(correlationTimeout) + 
                      ' ' + leftInputPath + ' ' + rightInputPath + ' ' + outputPrefix)
    print cmd
    os.system(cmd)

    print '*******'

    # Generate disparity images for debugging
    #cmd = 'disparitydebug ' + disparityImagePath
    #print cmd
    #os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(disparityImagePath):
        raise Exception('Stereo processing failed to create output file ' + disparityImagePath + 
                        ' from input files ' + leftInputPath + ' and ' + rightInputPath)

    # Compute percentage of good pixels
    percentGood = IrgAspFunctions.getStereoGoodPixelPercentage(outputPrefix)
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
        cmd = ['stereoIpFind', leftInputPath, rightInputPath, binaryPath,
                              '--elevationGuess', str(surfaceElevation)]
        print cmd
        # Call and record the output text
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        textOut, err = p.communicate()
        print textOut

        # Find the triangulation error in the output text
        basePos  = textOut.find('Triangulation Err:')
        startPos = textOut.find(':', basePos)
        splitPos = textOut.find('+-', basePos)
        endPos   = textOut.find('meters', startPos)
    
        # Quit if we didn't find the output error
        if (startPos < 0) or (splitPos < 0) or (endPos < 0):
            raise Exception('stereoIpFind call for file ' + outputPath + ' failed')
            
        baseError = float(textOut[startPos+1:splitPos-1])
        variance  = float(textOut[splitPos+2:endPos-1])

        # Check that the error is reasonable
        BASE_ERROR_LIMIT = 400.0
        VARIANCE_LIMIT   = 200.0
        if (baseError > BASE_ERROR_LIMIT) or (variance > VARIANCE_LIMIT):
           logging.info('Error: stereoIpFind file '+ outputPath  + ' has very high error: ' + textOut[basePos:endPos+6])
           raise Exception('Error: stereoIpFind returned very high error: ' + textOut[basePos:endPos+6])
           #print 'Warning: stereoIpFind returned very high error: ' + textOut[basePos:endPos+6]
           


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
        IrgFileFunctions.removeIfExists(binaryPath)

    # Count the number of point pairs we found
    numPairs = IrgFileFunctions.getFileLineCount(outputPath)
    print('For cubes %s and %s found %d point pairs' % (leftInputPath, rightInputPath, numPairs))
    logging.info('For cubes %s and %s', leftInputPath, rightInputPath)
    logging.info('- Found %d point pairs', numPairs)

    return numPairs


# Compares the backprojected locations of matched pixel pairs
def evaluateAccuracy(leftCubePath, rightCubePath, ipFindOutputPath, workDir='', savePoints=False):

    print 'Evaluating accuracy between cubes ' + leftCubePath + ' and ' + rightCubePath

    if savePoints:
        leftPointPath  = os.path.join(workDir, 'leftEvalPoints.csv')
        rightPointPath = os.path.join(workDir, 'rightEvalPoints.csv')
        IrgFileFunctions.createFolder(workDir)
        leftFile  = open(leftPointPath, 'w')
        rightFile = open(rightPointPath, 'w')

    # Run check every N lines in the file
    lineSkip = 15
    
    f = open(ipFindOutputPath)
    i = 0
    distanceList = []
    for line in f:
        if i % lineSkip == 0: # Compare locations on a sampling of lines
            # Get the matching pixel coordinates from this line
            leftSample, leftLine, rightSample, rightLine = line.split(',')
            
            #print 'LeftPixel =  ' + leftSample  + ', ' + leftLine
            #print 'RightPixel = ' + rightSample + ', ' + rightLine
            
            # Obtain the backprojected GCC location for the pixels in the two images
            leftLoc  = IrgIsisFunctions.getPixelLocInCube(leftCubePath,  leftSample,  leftLine)
            rightLoc = IrgIsisFunctions.getPixelLocInCube(rightCubePath, rightSample, rightLine)
            
            if savePoints:
                leftFile.write (str(leftLoc ['gdc'][1]) + ', ' +  str(leftLoc ['gdc'][0]) + ', ' + str(leftLoc ['gdc'][2]) + '\n');
                rightFile.write(str(rightLoc['gdc'][1]) + ', ' +  str(rightLoc['gdc'][0]) + ', ' + str(rightLoc['gdc'][2]) + '\n');
            
            #print 'leftLoc  = ' + str(leftLoc)
            #print 'rightLoc = ' + str(rightLoc)
            
            # Determine the distance between the GCC points and accumulate
            distance = math.sqrt(math.pow(leftLoc['gcc'][0] - rightLoc['gcc'][0], 2.0) + 
                                 math.pow(leftLoc['gcc'][1] - rightLoc['gcc'][1], 2.0) + 
                                 math.pow(leftLoc['gcc'][2] - rightLoc['gcc'][2], 2.0))
            #print 'Distance = %.2f' % distance
            distanceList.append(distance)
        i = i + 1
        
    if savePoints:
        leftFile.close()
        rightFile.close()
    
    if (i == 0):
        raise Exception("Can't evaluate pixel pairs, match file " + ipFindOutputPath + ' is empty!')
    
    # Return the median distance
    return numpy.median(distanceList)

def makeZeroParamsPath(outputPath):
    """Generates an all-zero parameter file for lronacAngleDoubleSolver to be used
       as an initial state for debugging checks"""
    
    NUM_PARAMS = 12   
    fakeLog = open(outputPath, 'w')
    for i in range(0, 12):
        fakeLog.write('0.0')
    fakeLog.close()
    
    return os.path.exists(outputPath)
    
## Generate a modified IK kernel to adjust the rotation between an LE/RE camera pair.
#def applyInterCameraPairRotation(leftInputPath, rightInputPath, newRotationPath, outputCubePath, 
#                                 ckPath, spkPath, workDir, forceOperation):
#
#    # Quit immediately if the output file already exists
#    if (not forceOperation) and (os.path.exists(outputCubePath)):
#        print 'File ' + outputCubePath + ' already exists, skipping nav transform.'
#        return True
#
#    # Generate the new file
#    cmdArgs = ['--keep', '--output', outputCubePath, 
#               '--rotation', newRotationPath, '--left', leftInputPath, 
#               '--right', rightInputPath, 
#               '--workDir', workDir]
#    if ckPath: # Add optional arguments
#        cmdArgs.append('--ck')
#        cmdArgs.append(ckPath)
#    if spkPath:
#        cmdArgs.append('--spk')
#        cmdArgs.append(spkPath)
#    #print cmdArgs
#    lronacCameraRotationCorrector.main(cmdArgs)
#
#    # Check to make sure we actually created the file
#    if not os.path.exists(outputCubePath):
#        raise Exception('Inter-camera rotation failed to create output file ' + outputCubePath + 
#                        ' from input files ' + leftInputPath + ' and ' + rightInputPath)
#
#    return True
    


# Looks in the pc_align output folder for the variably-named log file
def findLatestPcAlignOutputLog(folder):

    for root, dirs, files in os.walk(folder):
        for f in files:
            if 'log' in f:
                return os.path.join(folder, f)
    return False

# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IrgFileFunctions.checkIfToolExists('pc_align')
    IrgFileFunctions.checkIfToolExists('lronacAngleDoubleSolver')
    IrgFileFunctions.checkIfToolExists('lronac2isis')
    IrgFileFunctions.checkIfToolExists('lronaccal')
    IrgFileFunctions.checkIfToolExists('lronacecho')
    IrgFileFunctions.checkIfToolExists('spiceinit')
    IrgFileFunctions.checkIfToolExists('positionCorrector.py')
    IrgFileFunctions.checkIfToolExists('pixelPairsFromStereo')
    IrgFileFunctions.checkIfToolExists('rotationCorrector.py')
    IrgFileFunctions.checkIfToolExists('stereo_corr')
    IrgFileFunctions.checkIfToolExists('stereo_pprc')
    IrgFileFunctions.checkIfToolExists('stereoIpFind')
    IrgFileFunctions.checkIfToolExists('crop')

    return True

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
        expectedSurfaceElevation = IrgIsisFunctions.getCubeElevationEstimate(spiceInitLeftPath, tempFolder)


        # Apply LE/RE LRONAC position offsets to each of the input files
        print 'Applying known camera offset to main LE camera...'
        posOffsetCorrectedLeftPath = os.path.join(tempFolder, leftBaseName + '.posOffsetCorrected.cub')
        leftPosCorrectWorkDir      = os.path.join(tempFolder, 'leftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitLeftPath, posOffsetCorrectedLeftPath, leftPosCorrectWorkDir, carry)

        print 'Applying known camera offset to main RE camera...'
        posOffsetCorrectedRightPath = os.path.join(tempFolder, rightBaseName + '.posOffsetCorrected.cub')
        rightPosCorrectWorkDir      = os.path.join(tempFolder, 'rightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitRightPath, posOffsetCorrectedRightPath, rightPosCorrectWorkDir, carry)

        print 'Applying known camera offset to main LE camera...'
        posOffsetCorrectedStereoLeftPath = os.path.join(tempFolder, leftStereoBaseName + '.posOffsetCorrected.cub')
        stereoLeftPosCorrectWorkDir      = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoLeftPath, posOffsetCorrectedStereoLeftPath, 
                                       stereoLeftPosCorrectWorkDir, carry)

        print 'Applying known camera offset to stereo RE camera...'
        posOffsetCorrectedStereoRightPath = os.path.join(tempFolder, rightStereoBaseName + '.posOffsetCorrected.cub')
        stereoRightPosCorrectWorkDir      = os.path.join(tempFolder, 'stereoRightPosCorrectDir')
        applyInterCameraPositionOffset(spiceInitStereoRightPath, posOffsetCorrectedStereoRightPath,  
                                       stereoRightPosCorrectWorkDir, carry)
        positionTime = time.time()
        logging.info('Position offset complete in %f seconds', positionTime - initTime)

        print '\n-------------------------------------------------------------------------\n' # -------------------------------------------------


        # Use stereoIpFind to find a small (<500) set of paired pixels from each pair of images relatively quickly
        # - Image pairs that do not yield enough matched pixels will have all their matched pixels thrown out.

        MIN_NUM_SMALL_PIXELS = 75  # Must have at least this many matches between images, otherwise they are all probably junk
        
        numMainPairsSmall         = 0
        numStereoPairsSmall       = 0
        numLeftCrossPairsSmall    = 0
        numRightCrossPairsSmall   = 0
        leftPixelsCmdParams       = []
        rightPixelsCmdParams      = []
        leftCrossPixelsCmdParams  = []
        rightCrossPixelsCmdParams = []
        pixelPairsLeftSmall       = os.path.join(tempFolder, 'stereoPixelPairsLeftSmall.csv')
        pixelPairsRightSmall      = os.path.join(tempFolder, 'stereoPixelPairsRightSmall.csv')
        pixelPairsLeftCrossSmall  = os.path.join(tempFolder, 'stereoPixelPairsLeftCrossSmall.csv')
        pixelPairsRightCrossSmall = os.path.join(tempFolder, 'stereoPixelPairsRightCrossSmall.csv')
        pixelPairsLeftFail        = os.path.join(tempFolder, 'stereoPixelPairsLeftFail.txt') # Flag files
        pixelPairsRightFail       = os.path.join(tempFolder, 'stereoPixelPairsRightFail.txt')
        pixelPairsLeftCrossFail   = os.path.join(tempFolder, 'stereoPixelPairsLeftCrossFail.txt')
        pixelPairsRightCrossFail  = os.path.join(tempFolder, 'stereoPixelPairsRightCrossFail.txt')
        usingLeftPixels           = False
        usingRightPixels          = False
        usingLeftCrossPixels      = False
        usingRightCrossPixels     = False

        #TODO: Move this into a function
        if os.path.exists(pixelPairsLeftFail):
            print 'Already failed to find left pixel matches, skipping this step'
        else:
            try:  # Left in main pair to right in main pair
                print '--> Looking for Left-LeftStereo point matches...'
                numLeftPairsSmall  = getInterestPointPairs(posOffsetCorrectedLeftPath, posOffsetCorrectedStereoLeftPath, 
                                                           pixelPairsLeftSmall, expectedSurfaceElevation, carry)
                usingLeftPixels = (numLeftPairsSmall > MIN_NUM_SMALL_PIXELS)
            except Exception, e:
                print str(e)
                print 'Failed to find LE-LE match, ignoring this data source.'
                os.system('touch ' + pixelPairsLeftFail) # Set failure flag in case we need to re-run the script

        if os.path.exists(pixelPairsRightFail):
            print 'Already failed to find right pixel matches, skipping this step'
        else:
            try:  # Right in main pair to right in the stereo pair
                print '--> Looking for Right-RightStereo point matches...'
                numRightPairsSmall  = getInterestPointPairs(posOffsetCorrectedRightPath, posOffsetCorrectedStereoRightPath, 
                                                            pixelPairsRightSmall, expectedSurfaceElevation, carry)
                usingRightPixels = (numRightPairsSmall > MIN_NUM_SMALL_PIXELS)
            except Exception, e:
                print str(e)
                os.system('touch ' + pixelPairsRightFail) # Set failure flag in case we need to re-run the script
                print 'Failed to find RE-RE match, ignoring this data source.'
            
        if os.path.exists(pixelPairsLeftCrossFail):
            print 'Already failed to find left cross pixel matches, skipping this step'
        else:
            try:  # Left in main pair to right in the stereo pair
                print '--> Looking for Left-RightStereo point matches...'
                numLeftCrossPairsSmall  = getInterestPointPairs(posOffsetCorrectedLeftPath, posOffsetCorrectedStereoRightPath, 
                                                                pixelPairsLeftCrossSmall, expectedSurfaceElevation, carry)
                usingLeftCrossPixels = (numLeftCrossPairsSmall > MIN_NUM_SMALL_PIXELS)
            except Exception, e:
                print str(e)
                os.system('touch ' + pixelPairsLeftCrossFail) # Set failure flag in case we need to re-run the script
                print 'Failed to find left-cross match, ignoring this data source.'
            
        if os.path.exists(pixelPairsRightCrossFail):
            print 'Already failed to find right cross pixel matches, skipping this step'
        else:
            try:  # Left in the stereo pair to right in the main pair
                print '--> Looking for LeftStereo-Right point matches...'
                numRightCrossPairsSmall = getInterestPointPairs(posOffsetCorrectedStereoLeftPath, posOffsetCorrectedRightPath, 
                                                                pixelPairsRightCrossSmall, expectedSurfaceElevation, carry)
                usingRightCrossPixels = (numRightCrossPairsSmall > MIN_NUM_SMALL_PIXELS)
            except Exception, e:
                print str(e)
                os.system('touch ' + pixelPairsRightCrossFail) # Set failure flag in case we need to re-run the script
                print 'Failed to find right-cross match, ignoring this data source.'


        # Quit if we did not find enough matching pixels in any of the image pairs!
        if not (usingLeftPixels or usingRightPixels or usingLeftCrossPixels or usingRightCrossPixels):
            raise Exception('ERROR: No image pair had the minimum number of matching pixels found!')

        # If we are using each set of matching pixels, set up the command line parameters for it
        if usingLeftPixels:
            leftPixelsCmdParams        = ['--matchingPixelsLeftPath',       pixelPairsLeftSmall]
        if usingRightPixels:
            rightPixelsCmdParams       = ['--matchingPixelsRightPath',      pixelPairsRightSmall]
        if usingLeftCrossPixels:
            leftCrossPixelsCmdParams   = ['--matchingPixelsLeftCrossPath',  pixelPairsLeftCrossSmall]
        if usingRightCrossPixels:
            rightCrossPixelsCmdParams  = ['--matchingPixelsRightCrossPath', pixelPairsRightCrossSmall]

        #TODO: Log this to file
        print '\n<<< Pixel pair count summary >>>'
        print 'Left  to stereo left : ' + str(numLeftPairsSmall)
        print 'Right to stereo right: ' + str(numRightPairsSmall)
        print 'Left  to stereo right: ' + str(numLeftCrossPairsSmall)
        print 'Right to stereo left : ' + str(numRightCrossPairsSmall)

        pixelTime = time.time()
        logging.info('Pixel pair finding complete in %f seconds', pixelTime - positionTime)

        print '\n-------------------------------------------------------------------------\n' # -------------------------------------------------

        # Compute all rotations and translations between the four cubes
        # - Computes the following transforms:
        #   = main   LE   to main   RE   (camera rotation offset in IK file)
        #   = main   pair to stereo pair (applied to ck and spk kernel files)
        #   = stereo LE   to stereo RE   (camera rotation offset in IK file)
        sbaOutputPrefix     = os.path.join(tempFolder, 'SBA_solution')
        smallGdcFile        = sbaOutputPrefix + "-outputGdcPoints.csv"
        localRotationPath   = sbaOutputPrefix + "-localRotationMatrix.csv"       # RE to LE rotation
        globalTransformPath = sbaOutputPrefix + "-globalTransformMatrix.csv"     # Global transform of stereo cameras
        stereoRotationPath  = sbaOutputPrefix + "-stereoLocalRotationMatrix.csv" # Stereo RE to LE rotation
        solvedParamsPath    = sbaOutputPrefix + "-finalParamState.csv"
        mainIpFindPath      = sbaOutputPrefix + "-mainIpFindPixels.csv"
        stereoIpFindPath    = sbaOutputPrefix + "-stereoIpFindPixels.csv"
        if (not os.path.exists(globalTransformPath)) or carry:
            print 'Solving for inter-camera transforms...'
            cmd = ['lronacAngleDoubleSolver',  '--outputPrefix',            sbaOutputPrefix, 
                                               '--leftCubePath',            posOffsetCorrectedLeftPath, 
                                               '--rightCubePath',           posOffsetCorrectedRightPath, 
                                               '--leftStereoCubePath',      posOffsetCorrectedStereoLeftPath, 
                                               '--rightStereoCubePath',     posOffsetCorrectedStereoRightPath, 
                                               '--elevation',               str(expectedSurfaceElevation)]
            cmd = cmd + leftPixelsCmdParams + rightPixelsCmdParams + leftCrossPixelsCmdParams + rightCrossPixelsCmdParams
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

        # Apply the planet-centered rotation/translation to all but the main camera 
        #  (it is held constant while the others are moved around it)
        
# TODO: Change rotation variable names to reflect no local rotations!
        rightAdjustedPath = os.path.join(tempFolder, rightBaseName + '.stereoAdjusted.cub')
        applyNavTransform(posOffsetCorrectedRightPath, rightAdjustedPath,
                          localRotationPath, False, carry)

        leftStereoAdjustedPath = os.path.join(tempFolder, leftStereoBaseName + '.stereoAdjusted.cub')
        applyNavTransform(posOffsetCorrectedStereoLeftPath, leftStereoAdjustedPath,
                          globalTransformPath, False, carry)

        rightStereoAdjustedPath = os.path.join(tempFolder, rightStereoBaseName + '.stereoAdjusted.cub')
        applyNavTransform(posOffsetCorrectedStereoRightPath, rightStereoAdjustedPath,
                          stereoRotationPath, False, carry)



  
        #DEBUG: Now that the LE and RE images are corrected relative to each other, verify that our intersection error is still the same
        print 'Running check to look at moved small points'
        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        # - All this step does is use stereo intersection to determine a lat/lon/alt coordinate for each pixel pair in the large data set. No optimization is performed.
        smallGdcTestFolder = os.path.join(tempFolder, 'gdcPointsSmallTest/')
        if not os.path.exists(smallGdcTestFolder):
            os.mkdir(smallGdcTestFolder)
        smallGdcTestPrefix = os.path.join(tempFolder, 'gdcPointsSmallTest/out')
        smallGdcTestFile = smallGdcTestPrefix + '-initialGdcPoints.csv'
        zeroParamsPath = os.path.join(tempFolder, 'zeroParamsPath.csv')
        if (not os.path.exists(smallGdcTestFile)) or carry:
            makeZeroParamsPath(zeroParamsPath)
            # What are initial error computations using the files we just fixed?
            cmd = ['lronacAngleDoubleSolver', '--outputPrefix', smallGdcTestPrefix,
                   '--leftCubePath',        posOffsetCorrectedLeftPath,
                   '--rightCubePath',       rightAdjustedPath,
                   '--leftStereoCubePath',  leftStereoAdjustedPath,
                   '--rightStereoCubePath', rightStereoAdjustedPath,
                   '--initialOnly', '--initialValues', zeroParamsPath,
                   '--elevation', str(expectedSurfaceElevation)]
            cmd = cmd + leftPixelsCmdParams + rightPixelsCmdParams + leftCrossPixelsCmdParams + rightCrossPixelsCmdParams
          
            ## What are initial error computations using the parameters and files we just solved for?
            #cmd = ['lronacAngleDoubleSolver', '--outputPrefix', smallGdcTestPrefix,
            # '--leftCubePath', posOffsetCorrectedLeftPath,
            # '--rightCubePath', posOffsetCorrectedRightPath,
            # '--leftStereoCubePath', posOffsetCorrectedStereoLeftPath,
            # '--rightStereoCubePath', posOffsetCorrectedStereoRightPath,
            # '--initialOnly', '--initialValues', solvedParamsPath,
            # '--elevation', str(expectedSurfaceElevation)]
            #cmd = cmd + leftPixelsCmdParams + rightPixelsCmdParams + leftCrossPixelsCmdParams + rightCrossPixelsCmdParams
          
            print cmd
            p = subprocess.call(cmd)
        
        # Generate KML for the output!
        testKmlFile = os.path.join(tempFolder, 'testSmallGdc.kml')
        if (not os.path.exists(testKmlFile)) or carry:
           cmd = 'calibrationReport.py --size tiny --skip 10 --name smallGdcTestMovedPts --input ' + smallGdcTestFile + ' --output ' + testKmlFile
           print cmd
           os.system(cmd)
        
        else:
           print 'Skipping small GDC TEST file creation step'


        print '\n-------------------------------------------------------------------------\n' # ----------------------------------------------------

        # Now we need a large number of matching pixels to generate a point cloud to align with the LOLA points
        # - In order to do this we use the initial disparity generating steps of the ASP stereo process.

        # Choose the image pair with the largest number of matching pixels (in the small step) to do this with.
        # - TODO: Do this for all pairs above a threshold?
        print 'Running stereo correlation to obtain a set of dense stereo correspondence pairs...'
        CORRELATION_TIMEOUT = 400 # Max number of seconds to spend on each tile before giving up.
        stereoPrefixLeft       = ''
	stereoPrefixRight      = '' # Initialize these variables
        stereoPrefixLeftCross  = ''
        stereoPrefixRightCross = ''
        if ( (numLeftPairsSmall > numRightPairsSmall) and (numLeftPairsSmall > numLeftCrossPairsSmall) and (numLeftPairsSmall > numRightCrossPairsSmall) ): # Test left
            # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
            largePixelCmdParam = ' --matchingPixelsLeftPath '
            stereoPrefixLeft   = os.path.join(tempFolder, 'stereoOutputLeft/out')
            disparityImage     = callStereoCorrelation(posOffsetCorrectedLeftPath, posOffsetCorrectedStereoLeftPath, 
                                                       stereoPrefixLeft, CORRELATION_TIMEOUT, carry)
            
        elif ( (numRightPairsSmall > numLeftCrossPairsSmall) and (numRightPairsSmall > numRightCrossPairsSmall)): # Test right
            # Perform initial stereo step on two RE cubes to generate a large number of point correspondences
            largePixelCmdParam = ' --matchingPixelsRightPath '
            stereoPrefixRight  = os.path.join(tempFolder, 'stereoOutputRight/out')
            disparityImage     = callStereoCorrelation(posOffsetCorrectedRightPath, posOffsetCorrectedStereoRightPath, 
                                                       stereoPrefixRight, CORRELATION_TIMEOUT, False)#carry)
            
        elif (numLeftCrossPairsSmall > numRightCrossPairsSmall): # Test left cross
            # Perform initial stereo step on LE to stereo RE cubes to generate a large number of point correspondences
            largePixelCmdParam    = ' --matchingPixelsLeftCrossPath '
            stereoPrefixLeftCross = os.path.join(tempFolder, 'stereoOutputLeftCross/out')
            disparityImage        = callStereoCorrelation(posOffsetCorrectedLeftPath, posOffsetCorrectedStereoRightPath, 
                                                          stereoPrefixLeftCross, CORRELATION_TIMEOUT, carry)
            
        else: # Right cross is largest
            # Perform initial stereo step on stereo LE to RE cubes to generate a large number of point correspondences
            largePixelCmdParam     = ' --matchingPixelsRightCrossPath '
            stereoPrefixRightCross = os.path.join(tempFolder, 'stereoOutputRightCross/out')
            disparityImage         = callStereoCorrelation(posOffsetCorrectedStereoLeftPath, posOffsetCorrectedRightPath, 
                                                           stereoPrefixRightCross, CORRELATION_TIMEOUT, carry)

        # TODO: Allow this to happen from more than one input source?
        # Extract a large number of matching pixel locations (many thousands) from the stereo output.
        # - The skip number is a row and column skip.
        pixelPairsLarge    = os.path.join(tempFolder, 'stereoPixelPairsLarge.csv')
        # TODO: Determine the best value here
        numPixelPairsLarge = extractPixelPairsFromStereoResults(disparityImage, pixelPairsLarge, 10, carry)

        
        print '\n-------------------------------------------------------------------------\n' # --------------------------------------------------

        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        # - All this step does is use stereo intersection to determine a lat/lon/alt coordinate for each pixel pair in the large data set.  No optimization is performed.
        largeGdcFolder = os.path.join(tempFolder, 'gdcPointsLargeComp/')
        if not os.path.exists(largeGdcFolder):
            os.mkdir(largeGdcFolder)
        largeGdcPrefix = os.path.join(tempFolder, 'gdcPointsLargeComp/out')
        largeGdcFile   = largeGdcPrefix + '-initialGdcPoints.csv'
        if (not os.path.exists(largeGdcFile)) or carry:
            print 'Calculating GDC coordinates of large set of intersected points...'
            cmd = ('lronacAngleDoubleSolver --outputPrefix '            + largeGdcPrefix + 
                                          largePixelCmdParam            + pixelPairsLarge +
                                          ' --leftCubePath '            + posOffsetCorrectedLeftPath + 
                                          ' --leftStereoCubePath '      + posOffsetCorrectedStereoLeftPath +
                                          ' --rightCubePath '           + posOffsetCorrectedRightPath + 
                                          ' --rightStereoCubePath '     + posOffsetCorrectedStereoRightPath + 
                                          ' --initialOnly --initialValues ' + solvedParamsPath + 
                                          ' --elevation ' + str(expectedSurfaceElevation))
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        prepTime = time.time()
        logging.info('pc_align prep finished in %f seconds', prepTime - sbaTime)

        print '\n-------------------------------------------------------------------------\n' # ----------------------------------------------------

        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        
        # The max-displacement threshold will be adjusted until we are either using a certain number of LOLA points
        #  or until we are using a certain percentage of the input points.
        MIN_NUM_LOLA_POINTS    = 4000
        MIN_LOLA_PERCENTAGE    = 0.005
        STARTING_DISPLACEMENT  = 50
        MAX_MAX_DISPLACEMENT   = 1800
        DISPLACEMENT_INCREMENT = 50

        # Determine the number of points we want
        numLolaPoints          = IrgFileFunctions.getFileLineCount(options.lolaPath) - 1        
        currentMaxDisplacement = STARTING_DISPLACEMENT
        minNumPointsToUse      = min(MIN_NUM_LOLA_POINTS, MIN_LOLA_PERCENTAGE*float(numLolaPoints))
        
        print 'Using pc_align to compute a transform between intersected points and LOLA points...'
        print 'Starting pc_align, looking for ' + str(minNumPointsToUse) + ' lola point matches.' 
        logging.info('Starting pc_align, looking for %d lola point matches.', minNumPointsToUse)
        
        pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/dem')
        pcAlignTransformPath  = pcAlignOutputPrefix + '-inverse-transform.txt'
        pcAlignFolder         = os.path.dirname(pcAlignOutputPrefix)
        endErrorPath          = pcAlignOutputPrefix + '-end_errors.csv'
        if (not os.path.exists(pcAlignTransformPath)) or carry:
            
            while(True):
                cmd = ('pc_align --highest-accuracy --max-displacement ' + str(currentMaxDisplacement) + ' --datum D_MOON ' + 
                       '--save-inv-transformed-reference-points ' + largeGdcFile + 
                       ' ' + options.lolaPath + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only')
                print cmd
                os.system(cmd)
                
                if not os.path.exists(endErrorPath):
   	            numLolaPointsUsed = 0 # Failed to produce any output, maybe raising the error cutoff will help?
                else:
                    numLolaPointsUsed = IrgFileFunctions.getFileLineCount(endErrorPath) - 1
            
                if (numLolaPointsUsed >= minNumPointsToUse):
                    break # Success!
                elif (currentMaxDisplacement >= MAX_MAX_DISPLACEMENT): # Hit the maximum max limit!
                    raise Exception('Error! Unable to find a good value for max-displacement in pc_align.  Wanted '
                                    + str(minNumPointsToUse) + ' points, only found ' + str(numLolaPointsUsed))
                else: # Try again with a higher max limit
                    print ('Trying pc_align again, only got ' + str(numLolaPointsUsed)
                           + ' lola point matches with value ' + str(currentMaxDisplacement))
                    currentMaxDisplacement = currentMaxDisplacement + DISPLACEMENT_INCREMENT            
            
        else:
            print 'Skipping pc_align step'

        if not os.path.exists(pcAlignTransformPath):
            raise Exception('pc_align call failed!')

        # Copy the pc_align log to the output folder
        pcAlignLogPath = findLatestPcAlignOutputLog(pcAlignFolder)
        shutil.copyfile(pcAlignLogPath, os.path.join(outputFolder, 'pcAlignLog.txt'))
        

        alignTime = time.time()
        logging.info('pc_align finished in %f seconds', alignTime - prepTime)

        print '\n-------------------------------------------------------------------------\n'

        # Now go back and apply the pc_align computed transform to all four cameras.
        # - This step corrects the four camera positions relative to the LOLA data.
        print 'Applying pc_align transform to main LE image...'
        applyNavTransform(posOffsetCorrectedLeftPath, outputPathLeft, 
                          pcAlignTransformPath, True, carry)

        print 'Applying pc_align transform to main RE image...'
        applyNavTransform(posOffsetCorrectedRightPath, outputPathRight, 
                          pcAlignTransformPath, True, carry)

        print 'Applying pc_align transform to stereo LE image...'
        applyNavTransform(leftStereoAdjustedPath, outputPathStereoLeft, 
                          pcAlignTransformPath, True, carry)

        print 'Applying pc_align transform to stereo RE image...'
        applyNavTransform(rightStereoAdjustedPath, outputPathStereoRight, 
                          pcAlignTransformPath, True, carry)
        

        navTime = time.time()
        logging.info('Nav transforms finished in %f seconds', navTime - alignTime)


        print '\n-------------------------------------------------------------------------\n'

        #DEBUG: Now that the LE and RE images are fully corrected, see if our points match the pc_align transformed points
        print 'Running  check to look at moved large points'
        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        # - All this step does is use stereo intersection to determine a lat/lon/alt coordinate for each pixel pair in the large data set.  No optimization is performed.
        largeGdcTestFolder = os.path.join(tempFolder, 'gdcPointsLargeTest/')
        if not os.path.exists(largeGdcTestFolder):
            os.mkdir(largeGdcTestFolder)
        largeGdcTestPrefix = os.path.join(tempFolder, 'gdcPointsLargeTest/out')
        largeGdcTestFile   = largeGdcTestPrefix + '-initialGdcPoints.csv'
        zeroParamsPath     = os.path.join(tempFolder, 'zeroParamsPath.csv')
        if (not os.path.exists(largeGdcTestFile)) or carry:
            makeZeroParamsPath(zeroParamsPath)
            cmd = ('lronacAngleDoubleSolver --outputPrefix '            + largeGdcTestPrefix + 
                                          largePixelCmdParam            + pixelPairsLarge +
                                          ' --leftCubePath '            + outputPathLeft + 
                                          ' --leftStereoCubePath '      + outputPathStereoLeft +
                                          ' --rightCubePath '           + outputPathRight + 
                                          ' --rightStereoCubePath '     + outputPathStereoRight + 
                                          ' --initialOnly --initialValues ' + zeroParamsPath + 
                                          ' --elevation ' + str(expectedSurfaceElevation))
            print cmd
            os.system(cmd)

        # Generate KML for the output!
        testKmlFile = os.path.join(tempFolder, 'testLargeGdc.kml')
        if (not os.path.exists(testKmlFile)) or carry:
            cmd = 'calibrationReport.py --size tiny --skip 200 --name largeGdcTestMovedPts --input ' + largeGdcTestFile + ' --output ' + testKmlFile
            print cmd
            os.system(cmd)
            
        else:
            print 'Skipping large GDC TEST file creation step'
# ----------------

        print '\n-------------------------------------------------------------------------\n'
        print 'Starting last set of geocorrection accuracy checks...'
        
        # One last pair of checks to compute accuracy.
        # - This calls campt for pixel pairs and finds the GCC location difference.
        # - This is done for each image pairing that was used.
        medianMainError   = evaluateAccuracy(outputPathLeft,       outputPathRight, 
                                           mainIpFindPath,       os.path.join(tempFolder, 'finalGdcCheckMain'), True)
        medianStereoError = evaluateAccuracy(outputPathStereoLeft, outputPathStereoRight, 
                                           stereoIpFindPath,     os.path.join(tempFolder, 'finalGdcCheckStereo'), True)
        
        print        '=====> Median main   pair error = %.4f meters' % medianMainError
        print        '=====> Median stereo pair error = %.4f meters' % medianStereoError
        logging.info('=====> Median main   pair error = %.4f meters' % medianMainError)
        logging.info('=====> Median stereo pair error = %.4f meters' % medianStereoError)
        
        # DEBUG!!!!!!
        # Generate KML files to display the points
        leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckMain/leftEvalPoints.csv')
        leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckMain/leftEvalPoints.kml')
        cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftMainEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
        os.system(cmd)
        rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckMain/rightEvalPoints.csv')
        rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckMain/rightEvalPoints.kml')
        cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightMainEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
        os.system(cmd)
        
        leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckStereo/leftEvalPoints.csv')
        leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckStereo/leftEvalPoints.kml')
        cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftStereoEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
        os.system(cmd)
        rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckStereo/rightEvalPoints.csv')
        rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckStereo/rightEvalPoints.kml')
        cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightStereoEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
        os.system(cmd)
        
        
        
        # The next four checks are not used in every case
        if usingLeftPixels:
            medianLeftError   = evaluateAccuracy(outputPathLeft,       outputPathStereoLeft,
                                                pixelPairsLeftSmall,  os.path.join(tempFolder, 'finalGdcCheckLeft'), True)
            print        '=====> Median left  pair error = %.4f meters' % medianLeftError
            logging.info('=====> Median left  pair error = %.4f meters' % medianLeftError)            
        
            leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckLeft/leftEvalPoints.csv')
            leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckLeft/leftEvalPoints.kml')
            cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftLeftEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
            os.system(cmd)
            rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckLeft/rightEvalPoints.csv')
            rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckLeft/rightEvalPoints.kml')
            cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightLeftEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
            os.system(cmd)
        
        if usingRightPixels:
            medianRightError  = evaluateAccuracy(outputPathRight,      outputPathStereoRight, 
                                               pixelPairsRightSmall, os.path.join(tempFolder, 'finalGdcCheckRight'), True)
            print        '=====> Median right pair error = %.4f meters' % medianRightError
            logging.info('=====> Median right pair error = %.4f meters' % medianRightError)
        
            leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckRight/leftEvalPoints.csv')
            leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckRight/leftEvalPoints.kml')
            cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftRightEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
            os.system(cmd)
            rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckRight/rightEvalPoints.csv')
            rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckRight/rightEvalPoints.kml')
            cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightRightEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
            os.system(cmd)
        
        if usingLeftCrossPixels:
            medianLeftCrossError  = evaluateAccuracy(outputPathLeft,           outputPathStereoRight,
                                                    pixelPairsLeftCrossSmall,  os.path.join(tempFolder, 'finalGdcCheckLeftCross'), True)
            print        '=====> Median left  cross pair error = %.4f meters' % medianLeftCrossError
            logging.info('=====> Median left  cross pair error = %.4f meters' % medianLeftCrossError)
        
            leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckLeftCross/leftEvalPoints.csv')
            leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckLeftCross/leftEvalPoints.kml')
            cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftLeftCrossEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
            os.system(cmd)
            rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckLeftCross/rightEvalPoints.csv')
            rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckLeftCross/rightEvalPoints.kml')
            cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightLeftCrossEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
            os.system(cmd)
        
        
        if usingRightCrossPixels:
            medianRightCrossError = evaluateAccuracy(outputPathStereoLeft,     outputPathRight,
                                                    pixelPairsRightCrossSmall, os.path.join(tempFolder, 'finalGdcCheckRightCross'), True)
            print        '=====> Median right cross pair error = %.4f meters' % medianRightCrossError
            logging.info('=====> Median right cross pair error = %.4f meters' % medianRightCrossError)
        
            leftCsvFile = os.path.join(tempFolder, 'finalGdcCheckRightCross/leftEvalPoints.csv')
            leftKmlFile = os.path.join(tempFolder, 'finalGdcCheckRightCross/leftEvalPoints.kml')
            cmd = 'calibrationReport.py --color red --size tiny --skip 1 --name leftRightCrossEvalPoints --input ' + leftCsvFile + ' --output ' + leftKmlFile
            os.system(cmd)
            rightCsvFile = os.path.join(tempFolder, 'finalGdcCheckRightCross/rightEvalPoints.csv')
            rightKmlFile = os.path.join(tempFolder, 'finalGdcCheckRightCross/rightEvalPoints.kml')
            cmd = 'calibrationReport.py --color blue --size tiny --skip 1 --name rightRightCrossEvalPoints --input ' + rightCsvFile + ' --output ' + rightKmlFile
            os.system(cmd)
        
        # All finished!  We should have a fully calibrated version of each of the four input files.

        # Clean up temporary files
        if not options.keep:
            print 'Deleting temporary files'
            # Init files
            IrgFileFunctions.removeIfExists(spiceInitLeftPath)
            IrgFileFunctions.removeIfExists(spiceInitRightPath)
            IrgFileFunctions.removeIfExists(spiceInitStereoRightPath)
            IrgFileFunctions.removeIfExists(spiceInitStereoLeftPath)

            # Position correction files
            IrgFileFunctions.removeIfExists(posOffsetCorrectedLeftPath)
            IrgFileFunctions.removeIfExists(posOffsetCorrectedRightPath)
            IrgFileFunctions.removeIfExists(posOffsetCorrectedStereoLeftPath)
            IrgFileFunctions.removeIfExists(posOffsetCorrectedStereoRightPath)
            IrgFileFunctions.removeFolderIfExists(leftPosCorrectWorkDir)
            IrgFileFunctions.removeFolderIfExists(rightPosCorrectWorkDir)
            IrgFileFunctions.removeFolderIfExists(stereoLeftPosCorrectWorkDir)
            IrgFileFunctions.removeFolderIfExists(stereoRightPosCorrectWorkDir)
            # Stereo output
            IrgFileFunctions.removeFolderIfExists(os.path.dirname(stereoPrefixLeft))
            IrgFileFunctions.removeFolderIfExists(os.path.dirname(stereoPrefixRight))
            IrgFileFunctions.removeFolderIfExists(os.path.dirname(stereoPrefixLeftCross))
            IrgFileFunctions.removeFolderIfExists(os.path.dirname(stereoPrefixRightCross))
            #IrgFileFunctions.removeIfExists(pixelPairsLeftSmall)
            #IrgFileFunctions.removeIfExists(pixelPairsRightSmall)
            #IrgFileFunctions.removeIfExists(pixelPairsLeftCrossSmall)
            #IrgFileFunctions.removeIfExists(pixelPairsRightCrossSmall)

            ## Remove all the SBA files
            #fileList = [ f for f in os.listdir(tempFolder) if f.startswith("SBA_solution") ]
            #for f in fileList:
            #    IrgFileFunctions.removeIfExists(os.path.join(tempFolder, f))

            # Remove stereo corrected files
            IrgFileFunctions.removeIfExists(leftStereoAdjustedPath)
            #IrgFileFunctions.removeFolderIfExists(leftStereoCorrectWorkDir)
            IrgFileFunctions.removeIfExists(rightStereoAdjustedPath)
            #IrgFileFunctions.removeFolderIfExists(rightStereoCorrectWorkDir)
            # Clean pc_align steps
            IrgFileFunctions.removeIfExists(pixelPairsLarge)   # TODO: Update some of these for recent changes!
            #IrgFileFunctions.removeFolderIfExists(largeGdcFolder)
            #IrgFileFunctions.removeFolderIfExists(pcAlignFolder)
            # Clean out final file correction
            #IrgFileFunctions.removeIfExists(leftCkPath)
            #IrgFileFunctions.removeIfExists(leftSpkPath)
            #IrgFileFunctions.removeFolderIfExists(leftFinalWorkDir)
            #IrgFileFunctions.removeIfExists(rightCkPath)
            #IrgFileFunctions.removeIfExists(rightSpkPath)
            #IrgFileFunctions.removeFolderIfExists(rightFinalWorkDir)
            #IrgFileFunctions.removeIfExists(leftStereoCkPath)
            #IrgFileFunctions.removeIfExists(leftStereoSpkPath)
            #IrgFileFunctions.removeFolderIfExists(leftStereoFinalWorkDir)
            #IrgFileFunctions.removeIfExists(rightStereoCkPath)
            #IrgFileFunctions.removeIfExists(rightStereoSpkPath)
            #IrgFileFunctions.removeFolderIfExists(rightStereoFinalWorkDir)
            # Remove local transform folders
            #IrgFileFunctions.removeFolderIfExists(mainLocalWorkDir)
            #IrgFileFunctions.removeFolderIfExists(stereoLocalWorkDir)
            ## Remove check folders
            #IrgFileFunctions.removeFolderIfExists(os.path.join(tempFolder, 'pcAlignStereoGdcCheck'))
            #IrgFileFunctions.removeFolderIfExists(os.path.join(tempFolder, 'stereoGlobalAdjustGdcCheck'))
            #IrgFileFunctions.removeFolderIfExists(os.path.join(tempFolder, 'finalGdcCheck'))
            #IrgFileFunctions.removeFolderIfExists(os.path.join(tempFolder, 'finalStereoGdcCheck'))
            IrgFileFunctions.removeIfExists(zeroParamsPath)

            #if (hadToCreateTempFolder):
            #    IrgFileFunctions.removeFolderIfExists(tempFolder)

        endTime = time.time()

        print "Finished stereo calibration process in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
