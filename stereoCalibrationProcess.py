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

    # TODO: Is there a way to tell if this action has already been performed?
    cmd = 'spiceinit from= ' + echoCubFile
    print cmd
    os.system(cmd)


    return echoCubFile

#--------------------------------------------------------------------------------

def main():

    print "Started stereoCalibrationProcess.py"

    try:
        try:
            usage = "usage: rotationCorrector.py [--output <path>][--manual]\n  "
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
            cmd = '/home/smcmich1/repo/lronacPipeline/cmakeBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + spiceInitLeftPath + ' ' + spiceInitRightPath
            print cmd
            os.system(cmd)
        #print 'QUITTING EARLY'
        #return 0 

        # Apply LE/RE LRONAC position offsets to each of the input files
        posOffsetCorrectedLeftPath = os.path.join(tempFolder, 'left.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'leftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitLeftPath + ' --output ' + posOffsetCorrectedLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedRightPath = os.path.join(tempFolder, 'right.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitRightPath + ' --output ' + posOffsetCorrectedRightPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedStereoLeftPath = os.path.join(tempFolder, 'stereoLeft.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedStereoLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --keep --input ' + spiceInitStereoLeftPath + ' --output ' + posOffsetCorrectedStereoLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedStereoRightPath = os.path.join(tempFolder, 'stereoRight.posOffsetCorrected.cub')
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
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedRightPath
            print cmd
            os.system(cmd)
        
        #print 'QUITTING EARLY'
        #return 0 


        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefix   = os.path.join(tempFolder, 'stereoOutput')
        disparityImage = stereoPrefix + '-D.tif'
        if not os.path.exists(disparityImage):
            #cmd = 'stereo --entry-point 0 ' + options.leftPath + ' ' + options.rightPath + ' ' + stereoPrefix + ' --compute-low-res-disparity-only'
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefix
            print cmd
            os.system(cmd)
            # Stage 1 (slow!)
            cmd = 'stereo_corr --cost-mode 0 --corr-timeout 400 ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefix
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo step'

        if (not os.path.exists(disparityImage)):
            print 'Stereo processing failed!  Processing stopped.'
            return 0

        #print 'QUITTING EARLY'
        #return 0 

        # Extract a small number of matching pixel locations ( < 100)
        pixelPairsSmall = os.path.join(tempFolder, 'stereoPixelPairsSmall.csv')
        if not os.path.exists(pixelPairsSmall):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/stereoPixelPairExtractor -i ' + disparityImage + ' -o ' + pixelPairsSmall + ' -p 1700'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping small pair extraction step'



        # Compute the global rotation and offet between the two LE cubes
        globalParamsFile = os.path.join(tempFolder, 'solvedParams.csv')
        globalParamsFileMatrix = os.path.join(tempFolder, 'solvedParams.csv.matrix.csv') # TODO: Improve how this is handled!
        smallGdcFile     = os.path.join(tempFolder, 'gdcPointsSmall.csv')
        if True: #not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath ' + globalParamsFile + ' --gdcPointsOutPath ' + smallGdcFile + ' --matchingPixelsPath ' + pixelPairsSmall + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' --worldTransform --includePosition'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo transform calculation step'
        # Mean projection error here is 2.525159 with p=1700
        # Mean projection error here is 2.640027 with p=1699

        #print 'QUITTING EARLY'
        #return 0 

        # Apply the transform to the second camera pair!  The transform is in moon coordinates.
        tempLeftStereoPath = os.path.join(tempFolder, 'leftStereoAdjusted.cub')
        if True: #not os.path.exists(tempLeftStereoPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoLeftStereoCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedStereoLeftPath + ' --output ' + tempLeftStereoPath + ' --transformPath ' + globalParamsFileMatrix + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd) # TODO: Do we need to specify the SPK and CK file paths?

        tempRightStereoPath = os.path.join(tempFolder, 'rightStereoAdjusted.cub')
        if True: #not os.path.exists(tempLeftStereoPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoRightStereoCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedStereoRightPath + ' --output ' + tempRightStereoPath + ' --transformPath ' + globalParamsFileMatrix + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd) # TODO: Do we need to specify the SPK and CK file paths?

        # DEBUG: Check angle solver on stereo adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckStereo.csv')
        if True: #not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipeline/cmakeBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + tempLeftStereoPath + ' ' + tempRightStereoPath
            print cmd
            os.system(cmd)

        #print 'QUITTING EARLY'
        #return 0 


        # Extract a large number of matching pixel locations (thousands)
        # - Right now the skip is rows and columns
        pixelPairsLarge = os.path.join(tempFolder, 'stereoPixelPairsLarge.csv')
        if True: #not os.path.exists(pixelPairsLarge):
            cmd = '/home/smcmich1/repo/lronacPipeline/cmakeBuild/pixelPairsFromStereo -i ' + disparityImage + ' -o ' + pixelPairsLarge + ' -p 10'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large pair extraction step'

        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        largeGdcFile = os.path.join(tempFolder, 'gdcPointsLarge.csv')
        if True: #not os.path.exists(largeGdcFile):
            cmd = '/home/smcmich1/repo/lronacPipeline/cmakeBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + largeGdcFile + ' --matchingPixelsPath ' + pixelPairsLarge + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + " --worldTransform --includePosition --initialOnly --initialValues " + globalParamsFile 
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/dem')
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        #transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')
        transformMatrixFile   = pcAlignOutputPrefix + '-inverse-transform.txt'
        if True: #not os.path.exists(transformMatrixFile):
            #cmd = 'pc_align --highest-accuracy --max-displacement 600 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix + ' --compute-translation-only'
            cmd = 'pc_align --highest-accuracy --max-displacement 600 --datum D_MOON --save-inv-transformed-reference-points ' + largeGdcFile + ' ' + options.lolaPath + ' -o ' + pcAlignOutputPrefix# + ' --compute-translation-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'

        print 'QUITTING EARLY'
        return 0 

  #NOTE: Currently testing with identify transformation................................

# TODO: Try to identify some kernels we don't need to load in C++ in order to speed position correction up
#       Also maybe skip the unload step.


        # Now go back and apply the pc_align computed transform to the left and right input image
        leftCkPath  = os.path.join(tempFolder, 'leftFinalCk.bc') #TODO: Do we need to specify these?
        leftSpkPath = os.path.join(tempFolder, 'leftFinalSpk.bsp')
        if True: #not os.path.exists(options.outputPathLeft): # Do this no matter what?
            thisWorkDir = os.path.join(tempFolder, 'leftFullCorrection')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedLeftPath + ' --output ' + options.outputPathLeft + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + leftCkPath + ' --spk ' + leftSpkPath
            print cmd
            os.system(cmd)

        tempRightPath = os.path.join(tempFolder, 'partial_corrected_RE.cub')
        rightCkPath   = os.path.join(tempFolder, 'rightFinalCk.bc') # Need these later to pass to internal angle correction function
        rightSpkPath  = os.path.join(tempFolder, 'rightFinalSpk.bsp')
        if True: #not os.path.exists(tempRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightFullCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedRightPath + ' --output ' + tempRightPath + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath
            print cmd
            os.system(cmd)

        # Do the same correction for the stereo-left and stereo-right images
        leftStereoCkPath  = os.path.join(tempFolder, 'leftStereoFinalCk.bc') #msopck can't handle the default kernel paths
        leftStereoSpkPath = os.path.join(tempFolder, 'leftStereoFinalSpk.bsp')
        if True: #not os.path.exists(options.outputPathStereoLeft):
            thisWorkDir = os.path.join(tempFolder, 'leftStereoFullCorrection')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + tempLeftStereoPath + ' --output ' + options.outputPathStereoLeft + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + leftStereoCkPath + ' --spk ' + leftStereoSpkPath
            print cmd
            os.system(cmd)

        tempRightStereoPath2 = os.path.join(tempFolder, 'partial_corrected_stereo_RE.cub')
        rightStereoCkPath   = os.path.join(tempFolder, 'rightStereoFinalCk.bc') # Need these later to pass to internal angle correction function
        rightStereoSpkPath  = os.path.join(tempFolder, 'rightStereoFinalSpk.bsp')
        if True: #not os.path.exists(tempRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightStereoFullCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + tempRightStereoPath + ' --output ' + tempRightStereoPath2 + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + rightStereoCkPath + ' --spk ' + rightStereoSpkPath
            print cmd
            os.system(cmd)

        ## TODO: Re-run the point computations to see if they match the pc_align output!


        # At this point the left images are hopefully in the correct position and we can compute the offset of the RE camera


        # DEBUG: Check angle solver on adjusted LE/RE images!
        pairGdcCheckPath = os.path.join(tempFolder, 'pairGdcCheckFinal.csv')
        if True: #not os.path.exists(pairGdcCheckPath):
            cmd = '/home/smcmich1/repo/lronacPipeline/cmakeBuild/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + pairGdcCheckPath + ' ' + options.outputPathLeft + ' ' + tempRightPath
            print cmd
            os.system(cmd)
        
        #print 'QUITTING EARLY'
        #return 0 


        # Compute the local rotation between the adjusted LE and RE cubes
        checkGdcFile = os.path.join(tempFolder, 'gdcPointsCheckFinalRot.csv') # Record GDC points for debugging
        if True: #not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/lronacPipeline/lronacCameraRotationCorrector.py --keep --output ' + options.outputPathRight + ' --left ' + options.outputPathLeft + ' --right ' + tempRightPath + ' --gdcLogPath ' + checkGdcFile +  ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath 
            print cmd
            os.system(cmd)

        # Do the same thing for the two stereo cubes
        checkGdcFileStereo = os.path.join(tempFolder, 'gdcPointsCheckFinalRotStereo.csv') # Record GDC points for debugging
        if True: #not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/lronacPipeline/lronacCameraRotationCorrector.py --keep --output ' + options.outputPathStereoRight + ' --left ' + options.outputPathStereoLeft + ' --right ' + tempRightStereoPath2 + ' --gdcLogPath ' + checkGdcFileStereo +  ' --ck ' + rightStereoCkPath + ' --spk ' + rightStereoSpkPath 
            print cmd
            os.system(cmd)


        # Clean up temporary files
#        if not options.keep:
#            os.remove(tempTextPath)

      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
