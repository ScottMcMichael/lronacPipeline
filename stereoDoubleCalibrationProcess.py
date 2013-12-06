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

job_pool = [];

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Applies the LROC offset from spacecraft position to an LROC cube's spice data
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def add_job( cmd, num_working_threads=4 ):
    if ( len(job_pool) >= num_working_threads):
        job_pool[0].wait();
        job_pool.pop(0);
    print cmd;
    job_pool.append( subprocess.Popen(cmd, shell=True) );

def wait_on_all_jobs():
    print "Waiting for jobs to finish";
    while len(job_pool) > 0:
        job_pool[0].wait();
        job_pool.pop(0);


#--------------------------------------------------------------------------------

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
    cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + inputCubePath + ' --output ' + outputCubePath + ' --transformPath ' + transformMatrixPath + ' --workDir ' + workDir + ckLine + spkLine
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCubePath):
        raise Exception('Nav transform failed to create output file ' + outputCubePath + ' from input file ' + inputCubePath)

    return True

#--------------------------------------------------------------------------------

def main():

    print "Started stereoDoubleCalibrationProcess.py"

    try:
        try:
            usage = "usage: stereoDoubleCalibrationProcess.py TODO [--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--left",  dest="leftPath",  help="Path to LE .IMG file")
            parser.add_option("--right", dest="rightPath", help="Path to RE .IMG file")
            
            parser.add_option("--stereo-left", dest="stereoLeft", help="Path to LE .IMG file with overlapping view of --left file")
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

            #TODO: Use default and required tags above!
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
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckPre.csv')
        if not os.path.exists(pairGdcCheckPath):
            #cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + spiceInitLeftPath + ' ' + spiceInitRightPath
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + spiceInitLeftPath + ' ' + spiceInitRightPath
            print cmd
            os.system(cmd)
        #print 'QUITTING EARLY'
        #return 0 

#TODO: Do we need to specify all these SPK files?
        # Apply LE/RE LRONAC position offsets to each of the input files
        posOffsetCorrectedLeftPath    = os.path.join(tempFolder, 'left.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'leftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitLeftPath + ' --output ' + posOffsetCorrectedLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedRightPath    = os.path.join(tempFolder, 'right.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitRightPath + ' --output ' + posOffsetCorrectedRightPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedStereoLeftPath    = os.path.join(tempFolder, 'stereoLeft.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedStereoLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitStereoLeftPath + ' --output ' + posOffsetCorrectedStereoLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedStereoRightPath    = os.path.join(tempFolder, 'stereoRight.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedStereoRightPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoRightPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitStereoRightPath + ' --output ' + posOffsetCorrectedStereoRightPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        # Fail if required files are not present
        if (not os.path.exists(posOffsetCorrectedLeftPath)) or (not os.path.exists(posOffsetCorrectedRightPath)) or (not os.path.exists(posOffsetCorrectedStereoLeftPath)):
            print 'Missing required processed input files!  Processing stopped.'
            return 0


        # DEBUG: Check angle solver on input LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckInitial.csv')
        if not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedRightPath
            print cmd
            os.system(cmd)
        
        #print 'QUITTING EARLY'
        #return 0 


        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefixLeft   = os.path.join(tempFolder, 'stereoOutputLeft')
        disparityImageLeft = stereoPrefixLeft + '-D.tif'
        if not os.path.exists(disparityImageLeft):
            #cmd = 'stereo --entry-point 0 ' + options.leftPath + ' ' + options.rightPath + ' ' + stereoPrefix + ' --compute-low-res-disparity-only'
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefixLeft
            print cmd
            os.system(cmd)
            # Stage 1 (slow!)
            cmd = 'stereo_corr --cost-mode 0 --corr-timeout 400 ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefixLeft
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Left stereo step'
        if (not os.path.exists(disparityImageLeft)):
            print 'Left stereo processing failed!  Processing stopped.'
            return 0

        # Perform initial stereo step on two RE cubes to generate a large number of point correspondences
        stereoPrefixRight   = os.path.join(tempFolder, 'stereoOutputRight')
        disparityImageRight = stereoPrefixRight + '-D.tif'
        if not os.path.exists(disparityImageRight):
            #cmd = 'stereo --entry-point 0 ' + options.leftPath + ' ' + options.rightPath + ' ' + stereoPrefix + ' --compute-low-res-disparity-only'
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + posOffsetCorrectedRightPath + ' ' + posOffsetCorrectedStereoRightPath + ' ' + stereoPrefixRight
            print cmd
            os.system(cmd)
            # Stage 1 (slow!)
            cmd = 'stereo_corr --cost-mode 0 --corr-timeout 400 ' + posOffsetCorrectedRightPath + ' ' + posOffsetCorrectedStereoRightPath + ' ' + stereoPrefixRight
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Right stereo step'
        if (not os.path.exists(disparityImageRight)):
            print 'Right stereo processing failed!  Processing stopped.'
            return 0

        # Extract a small number of matching pixel locations from the disparity images ( < 300)
        pixelPairsLeftSmall = os.path.join(tempFolder, 'stereoPixelPairsLeftSmall.csv')
        if not os.path.exists(pixelPairsLeftSmall):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageLeft + ' -o ' + pixelPairsLeftSmall + ' -p 800'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping left small pair extraction step'

        pixelPairsRightSmall = os.path.join(tempFolder, 'stereoPixelPairsRightSmall.csv')
        if not os.path.exists(pixelPairsRightSmall):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageRight + ' -o ' + pixelPairsRightSmall + ' -p 800'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping right small pair extraction step'

        #print 'QUITTING EARLY'
        #return 0 



# TESTING - Perform cross-stereo matching of LE/RE cubes from opposite pairs

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
        stereoPrefixLeftCross   = os.path.join(tempFolder, 'stereoOutputLeftCross')
        disparityImageLeftCross = stereoPrefixLeftCross + '-D.tif'
        if not os.path.exists(disparityImageLeftCross):
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + leftPosCorrectedCropped + ' ' + rightStereoPosCorrectedCropped + ' ' + stereoPrefixLeftCross
            print cmd
            os.system(cmd)
            # Stage 1 (slow!)
            cmd = 'stereo_corr --cost-mode 0 --corr-timeout 100 ' + leftPosCorrectedCropped + ' ' + rightStereoPosCorrectedCropped + ' ' + stereoPrefixLeftCross
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Left to Right stereo step'
        if (not os.path.exists(disparityImageLeftCross)):
            print 'Left Cross stereo processing failed!  Processing stopped.'
            return 0

        # Next is left in the stereo pair to right in the main pair
        stereoPrefixRightCross   = os.path.join(tempFolder, 'stereoOutputRightCross')
        disparityImageRightCross = stereoPrefixRightCross + '-D.tif'
        if not os.path.exists(disparityImageRightCross):
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + leftStereoPosCorrectedCropped + ' ' + rightPosCorrectedCropped + ' ' + stereoPrefixRightCross
            print cmd
            os.system(cmd)
            # Stage 1 (slow!)
            cmd = 'stereo_corr --cost-mode 0 --corr-timeout 100 ' + leftStereoPosCorrectedCropped + ' ' + rightPosCorrectedCropped + ' ' + stereoPrefixRightCross
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Right to Left stereo step'
        if (not os.path.exists(disparityImageRightCross)):
            print 'Right to Left stereo processing failed!  Processing stopped.'
            return 0


        # Extract a small number of matching pixel locations from the disparity images ( < 300)
        pixelPairsLeftCrossSmall = os.path.join(tempFolder, 'stereoPixelPairsLeftCrossSmall.csv')
        if not os.path.exists(pixelPairsLeftCrossSmall):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageLeftCross + ' -o ' + pixelPairsLeftCrossSmall + ' -p 400'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Left to Right small pair extraction step'

        pixelPairsRightCrossSmall = os.path.join(tempFolder, 'stereoPixelPairsRightCrossSmall.csv')
        if not os.path.exists(pixelPairsRightCrossSmall):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageRightCross + ' -o ' + pixelPairsRightCrossSmall + ' -p 400'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping Right to Left small pair extraction step'


        #print 'QUITTING EARLY'
        #return 0 

# DONE TESTING




        # Compute all rotations and translations between the four cubes
        sbaOutputPrefix     = os.path.join(tempFolder, 'SBA_solution')
        smallGdcFile        = sbaOutputPrefix + "-outputGdcPoints.csv"
        localRotationPath   = sbaOutputPrefix + "-localRotationMatrix.csv"
        globalTransformPath = sbaOutputPrefix + "-globalTransformMatrix.csv"
        stereoRotationPath  = sbaOutputPrefix + "-stereoLocalRotationMatrix.csv"
        if not os.path.exists(globalTransformPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefix + ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --matchingPixelsRightPath ' + pixelPairsRightSmall + ' --matchingPixelsLeftCrossPath ' + pixelPairsLeftCrossSmall +  ' --matchingPixelsRightCrossPath ' + pixelPairsRightCrossSmall + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedRightPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + posOffsetCorrectedStereoRightPath
            #cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleDoubleSolver --outputPrefix ' + sbaOutputPrefix + ' --matchingPixelsLeftPath ' + pixelPairsLeftSmall + ' --matchingPixelsRightPath ' + pixelPairsRightSmall + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedRightPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + posOffsetCorrectedStereoRightPath
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo transform calculation step'

        #print 'QUITTING EARLY'
        #return 0 

        # Apply the transform to the second camera pair!  The transform is in moon coordinates.
        tempLeftStereoPath = os.path.join(tempFolder, 'leftStereoAdjusted.cub')
        if not os.path.exists(tempLeftStereoPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoLeftStereoCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedStereoLeftPath + ' --output ' + tempLeftStereoPath + ' --transformPath ' + globalTransformPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd) # TODO: Do we need to specify the SPK and CK file paths?

        tempRightStereoPath = os.path.join(tempFolder, 'rightStereoAdjusted.cub')
        if not os.path.exists(tempLeftStereoPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoRightStereoCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedStereoRightPath + ' --output ' + tempRightStereoPath + ' --transformPath ' + globalTransformPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd) # TODO: Do we need to specify the SPK and CK file paths?

        # TODO: Move all this stuff to the double solver?
        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckStereo.csv')
        if not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + tempLeftStereoPath + ' ' + tempRightStereoPath
            print cmd
            os.system(cmd)

        #print 'QUITTING EARLY'
        #return 0 


        # Extract a large number of matching pixel locations (thousands) from the LE/LE and RE/RE stereo.
        # - Right now the skip is rows and columns
        pixelPairsLeftLarge = os.path.join(tempFolder, 'stereoPixelPairsLeftLarge.csv')
        if not os.path.exists(pixelPairsLeftLarge):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageLeft + ' -o ' + pixelPairsLeftLarge + ' -p 16'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping left large pair extraction step'

        # TODO: Many changes needed before RE images can be used here!
#        pixelPairsRightLarge = os.path.join(tempFolder, 'stereoPixelPairsRightLarge.csv')
#        if not os.path.exists(pixelPairsRightLarge):
#            cmd = '/home/smcmich1/repo/lronacPipelineBuild/pixelPairsFromStereo -i ' + disparityImageRight + ' -o ' + pixelPairsRightLarge + ' -p 16'
#            print cmd
#            os.system(cmd)
#        else:
#            print 'Skipping right large pair extraction step'

        pixelPairsLarge = os.path.join(tempFolder, 'stereoPixelPairsLarge.csv')
        if not os.path.exists(pixelPairsLarge):
            #TODO: Can't use RE pixel pairs without applying their local correction earlier!
            #cmd = 'cat < ' + pixelPairsLeftLarge + ' < ' + pixelPairsRightLarge + ' > ' + pixelPairsLarge
            cmd = 'cp ' + pixelPairsLeftLarge + ' ' + pixelPairsLarge
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large pair merging extraction step'

        #print 'QUITTING EARLY'
        #return 0 # Don't forget to set the cropped file!


# TODO: Move this functionality to the double solver so it does not need to be hacked
        globalParamsFile  = sbaOutputPrefix + "-finalParamState-cropped.csv"
        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        largeGdcFile = os.path.join(tempFolder, 'gdcPointsLarge.csv')
        if True:#not os.path.exists(largeGdcFile):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + largeGdcFile + ' --matchingPixelsPath ' + pixelPairsLarge + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + " --worldTransform --includePosition --initialOnly --initialValues " + globalParamsFile 
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        # TODO: Why does rotation always move the points somewhere else?
        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/dem')
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        #transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')
        transformMatrixFile   = pcAlignOutputPrefix + '-inverse-transform.txt'
        if True:#not os.path.exists(transformMatrixFile):
            #cmd = 'pc_align --highest-accuracy --max-displacement 600 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            cmd = 'pc_align --highest-accuracy --max-displacement 600 --datum D_MOON --save-inv-transformed-reference-points ' + largeGdcFile + ' ' + options.lolaPath + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'

        #print 'QUITTING EARLY'
        #return 0 

  #NOTE: Currently testing with identify transformation................................

# TODO: Try to identify some kernels we don't need to load in C++ in order to speed position correction up
#       Also maybe skip the unload step.


        # Now go back and apply the pc_align computed transform to the left and right input image
        leftCkPath  = os.path.join(tempFolder, 'leftFinalCk.bc') #TODO: Do we need to specify these?
        leftSpkPath = os.path.join(tempFolder, 'leftFinalSpk.bsp')
        if True:#not os.path.exists(options.outputPathLeft): # Do this no matter what?
            thisWorkDir = os.path.join(tempFolder, 'leftFullCorrection')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedLeftPath + ' --output ' + options.outputPathLeft + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + leftCkPath + ' --spk ' + leftSpkPath
            print cmd
            os.system(cmd)

        tempRightPath = os.path.join(tempFolder, 'partial_corrected_RE.cub')
        rightCkPath   = os.path.join(tempFolder, 'rightFinalCk.bc') # Need these later to pass to internal angle correction function
        rightSpkPath  = os.path.join(tempFolder, 'rightFinalSpk.bsp')
        if True:#not os.path.exists(tempRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightFullCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedRightPath + ' --output ' + tempRightPath + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath
            print cmd
            os.system(cmd)

        # Do the same correction for the stereo-left and stereo-right images
        leftStereoCkPath  = os.path.join(tempFolder, 'leftStereoFinalCk.bc') #msopck can't handle the default kernel paths
        leftStereoSpkPath = os.path.join(tempFolder, 'leftStereoFinalSpk.bsp')
        if True:#not os.path.exists(options.outputPathStereoLeft):
            thisWorkDir = os.path.join(tempFolder, 'leftStereoFullCorrection')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + tempLeftStereoPath + ' --output ' + options.outputPathStereoLeft + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + leftStereoCkPath + ' --spk ' + leftStereoSpkPath
            print cmd
            os.system(cmd)

        tempRightStereoPath2 = os.path.join(tempFolder, 'partial_corrected_stereo_RE.cub')
        rightStereoCkPath   = os.path.join(tempFolder, 'rightStereoFinalCk.bc') # Need these later to pass to internal angle correction function
        rightStereoSpkPath  = os.path.join(tempFolder, 'rightStereoFinalSpk.bsp')
        if True:#not os.path.exists(tempRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightStereoFullCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + tempRightStereoPath + ' --output ' + tempRightStereoPath2 + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + rightStereoCkPath + ' --spk ' + rightStereoSpkPath
            print cmd
            os.system(cmd)

        ## TODO: Re-run the point computations to see if they match the pc_align output!

        # At this point the left images are hopefully in the correct position and we can apply the offset of the RE cameras

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckMidStereo.csv')
        if not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + options.outputPathStereoLeft + ' ' + tempRightStereoPath2
            print cmd
            os.system(cmd)


        #print 'QUITTING EARLY'
        #return 0 

        # Apply local transforms to both pairs of images!

        # Apply the local rotation to the adjusted RE cube
        if True:#not os.path.exists(options.outputPathRight):
            cmd = '/home/smcmich1/repo/lronacPipeline/lronacCameraRotationCorrector.py --keep --output ' + options.outputPathRight + ' --rotation ' + localRotationPath + ' --left ' + options.outputPathLeft + ' --right ' + tempRightPath + ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath 
            print cmd
            os.system(cmd)

        # Apply the local rotation to the adjusted stereo RE cube
        if True:#not os.path.exists(options.outputPathStereoRight):
            cmd = '/home/smcmich1/repo/lronacPipeline/lronacCameraRotationCorrector.py --keep --output ' + options.outputPathStereoRight + ' --rotation ' + stereoRotationPath + ' --left ' + options.outputPathStereoLeft + ' --right ' + tempRightStereoPath2 + ' --ck ' + rightStereoCkPath + ' --spk ' + rightStereoSpkPath 
            print cmd
            os.system(cmd)

        # DEBUG: Check angle solver on adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckFinal.csv')
        if True:#not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + options.outputPathLeft + ' ' + options.outputPathRight
            print cmd
            os.system(cmd)

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckFinalStereo.csv')
        if True:#not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipelineBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + options.outputPathStereoLeft + ' ' + options.outputPathStereoRight
            print cmd
            os.system(cmd)


        # Clean up temporary files
#        if not options.keep:
#            os.remove(tempTextPath)

      

        endTime = time.time()

        print "Finished stereo calibration process in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
