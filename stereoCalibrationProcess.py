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

# Reads in a lronacAngleSolver state file and produces a string of the outputs
def readSolvedState(filePath):
    
    if not os.path.exists(filePath):
        print 'Error: file ' + filePath + ' does not exist!'
        return ''

    # Read through file one line at a time
    f = open(filePath)
    #params = []
    paramString = ''
    limit = 6
    for line in f:
        #params.append(float(line)) # Each line just contains a floating point value
        newNum       = float(line.strip())
        newNumString = '%f' % (newNum) # Remove any scientific notation
        paramString  = paramString + newNumString + ' '
        limit        = limit - 1
        if limit == 0: # Only read in the first six parameters
            break
    f.close()

    return paramString


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

#TODO: Support for file based logging of results

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

            parser.add_option("--outputL", dest="outputPathLeft",
                              help="Where to write the output LE file.")
            parser.add_option("--outputR", dest="outputPathRight",
                              help="Where to write the output RE file.")
            
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
            if not options.lolaPath: 
                parser.error("Need LOLA DEM path")
            if not options.outputPathLeft: 
                parser.error("Need left output path")
            if not options.outputPathRight: 
                parser.error("Need right output path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder  = os.path.dirname(options.outputPathLeft)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp3/'
        if (options.workDir):
            tempFolder = options.workDir
        if not os.path.exists(tempFolder):
            os.mkdir(tempFolder) 

  
        # TODO: Wrap the following functions into one call!
  
        # Convert the input files from IMG files to spiceinit'ed cubes in the output folder
        spiceInitLeftPath       = prepareImgFile(options.leftPath,   tempFolder)
        spiceInitRightPath      = prepareImgFile(options.rightPath,  tempFolder)
        spiceInitStereoLeftPath = prepareImgFile(options.stereoLeft, tempFolder)

        # Apply LE/RE LRONAC position offsets to each of the input files
        posOffsetCorrectedLeftPath = os.path.join(tempFolder, 'left.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'leftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ' + spiceInitLeftPath + ' --output ' + posOffsetCorrectedLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedRightPath = os.path.join(tempFolder, 'right.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ' + spiceInitRightPath + ' --output ' + posOffsetCorrectedRightPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        posOffsetCorrectedStereoLeftPath = os.path.join(tempFolder, 'stereoLeft.posOffsetCorrected.cub')
        if not os.path.exists(posOffsetCorrectedStereoLeftPath):
            thisWorkDir = os.path.join(tempFolder, 'stereoLeftPosCorrectDir')
            cmd = '/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ' + spiceInitStereoLeftPath + ' --output ' + posOffsetCorrectedStereoLeftPath + ' --workDir ' + thisWorkDir
            print cmd
            os.system(cmd)

        # Fail if required files are not present
        if (not os.path.exists(posOffsetCorrectedLeftPath)) or (not os.path.exists(posOffsetCorrectedRightPath)) or (not os.path.exists(posOffsetCorrectedStereoLeftPath)):
            print 'Missing required processed input files!  Processing stopped.'
            return 0

        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefix   = os.path.join(tempFolder, 'stereoOutput')
        disparityImage = stereoPrefix + '-D.tif'
        if not os.path.exists(disparityImage):
            #cmd = 'stereo --entry-point 0 ' + options.leftPath + ' ' + options.rightPath + ' ' + stereoPrefix + ' --compute-low-res-disparity-only'
            # Stage 0 (fast)
            cmd = 'stereo_pprc ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefix
            print cmd
            os.system(cmd)
            # Stage 1 (sloooow)
            cmd = 'stereo_corr ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' ' + stereoPrefix
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo step'

        if (not os.path.exists(disparityImage)):
            print 'Stereo processing failed!  Processing stopped.'
            return 0

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
        smallGdcFile     = os.path.join(tempFolder, 'gdcPointsSmall.csv')
        if not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath ' + globalParamsFile + ' --gdcPointsOutPath ' + smallGdcFile + ' --matchingPixelsPath ' + pixelPairsSmall + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + ' --worldTransform --includePosition'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo transform calculation step'


        ## Read in the final state parameters (global rotation and shift)
        #initialParamString = readSolvedState(globalParamsFile)
        #if not initialParamString:
        #    return 1

        # Extract a large number of matching pixel locations (several thousand)
        pixelPairsLarge = os.path.join(tempFolder, 'stereoPixelPairsLarge.csv')
        if not os.path.exists(pixelPairsLarge):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/stereoPixelPairExtractor -i ' + disparityImage + ' -o ' + pixelPairsLarge + ' -p 100'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large pair extraction step'

        # Compute the 3d coordinates for each pixel pair using the rotation and offset computed earlier
        largeGdcFile = os.path.join(tempFolder, 'gdcPointsLarge.csv')
        if not os.path.exists(largeGdcFile):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + largeGdcFile + ' --matchingPixelsPath ' + pixelPairsLarge + ' ' + posOffsetCorrectedLeftPath + ' ' + posOffsetCorrectedStereoLeftPath + " --worldTransform --includePosition --initialOnly --initialValues " + globalParamsFile # + initialParamString + "'"
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        # TODO: Apply rotation and offset to second pair (need absolute rot changer?)

#        print 'QUITTING EARLY'
#        return 0 

        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        #transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')
        transformMatrixFile   = os.path.join(tempFolder, 'pcAlignOutput-transform.txt')
        if not os.path.exists(transformMatrixFile):
            pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput/')
            cmd = 'pc_align --highest-accuracy --max-displacement 600 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix# + ' --compute-translation-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'


        # Now go back and apply the pc_align computed transform to the left and right input image
        leftCkPath  = os.path.join(tempFolder, 'leftFinalCk.bc') #TODO: Do we need to specify these?
        leftSpkPath = os.path.join(tempFolder, 'leftFinalSpk.bsp')
        if not os.path.exists(options.outputPathLeft): # Do this no matter what?
            thisWorkDir = os.path.join(tempFolder, 'leftFullCorrection')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedLeftPath + ' --output ' + options.outputPathLeft + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + leftCkPath + ' --spk ' + leftSpkPath
            print cmd
            os.system(cmd)

        ## DEBUG: Correct the stereo-left image and re-run the point computations to see if they match the pc_align output
        #stereoLeftOutPath = os.path.join(tempFolder, 'stereo-left_final.cub') # TODO: Get from output path
        #if not os.path.exists(stereoLeftOutPath):
        #    cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + options.stereoLeftPath + ' --output ' + stereoLeftOutPath + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir
        #    print cmd
        #    os.system(cmd)

        tempRightPath = os.path.join(tempFolder, 'partial_corrected_RE.cub')
        rightCkPath  = os.path.join(tempFolder, 'rightFinalCk.bc') # Need these later to pass to internal angle correction function
        rightSpkPath = os.path.join(tempFolder, 'rightFinalSpk.bsp')
        if not os.path.exists(tempRightPath):
            thisWorkDir = os.path.join(tempFolder, 'rightFullCorrection/')
            cmd = '/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --keep --input ' + posOffsetCorrectedRightPath + ' --output ' + tempRightPath + ' --transformPath ' + transformMatrixFile + ' --workDir ' + thisWorkDir + ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath
            print cmd
            os.system(cmd)

        # At this point the left image is hopefully in the correct position and we can compute the offset of the RE camera

        # Compute the local rotation between the adjusted LE and RE cubes                                                                                                                                                                                                                                        
        checkGdcFile = os.path.join(tempFolder, 'gdcPointsCheck.csv') # Record GDC points for debugging
        if True: #not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/lronacPipeline/lronacCameraRotationCorrector.py --keep --output ' + options.outputPathRight + ' --left ' + options.outputPathLeft + ' --right ' + tempRightPath + ' --gdcLogPath ' + checkGdcFile +  ' --ck ' + rightCkPath + ' --spk ' + rightSpkPath 
            print cmd
            os.system(cmd)

        # If we are feeling generous we could also correct the "stereo-left" input cube and its RE pair
        # - This would require that we first apply the rotation correction found earlier, then the one found using pc_align



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
