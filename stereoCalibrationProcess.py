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
        paramString = paramString + line.strip() + ' '
        limit = limit - 1
        if limit == 0: # Only read in the first six parameters
            break
    f.close()

    return paramString

#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    print "Started stereoCalibrationProcess.py"

    try:
        try:
            usage = "usage: rotationCorrector.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--left",  dest="leftPath",  help="Path to LE .cub file")
            parser.add_option("--right", dest="rightPath", help="Path to RE .cub file")
            parser.add_option("--lola",  dest="lolaPath",  help="Path to LOLA DEM")

            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output (RE) file.")
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.lolaPath: 
                parser.error("Need LOLA DEM path")
#            if not options.outputPath: 
#                parser.error("Need output path")


            #TODO: Find LE and RE paths in args

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        outputFolder  = os.path.dirname(options.outputPath)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
#        if not os.path.exists(tempFolder):
#            os.mkdir(tempFolder) 

        # File must already have had spiceinit called

        #TODO: Extract all file names!

        # Perform initial stereo step on two LE cubes to generate a large number of point correspondences
        stereoPrefix   = os.path.join(tempFolder, 'stereoOutput')
        disparityImage = stereoPrefix + '-D_sub.tif'
        if not os.path.exists(disparityImage):
            cmd = 'stereo --entry-point 0 ' + options.leftPath + ' ' + options.rightPath + ' ' + stereoPrefix + ' --compute-low-res-disparity-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo step'
        #TODO: Make this only perform the requested stages!

        # Extract a small number of matching pixel locations ( < 100)
        
        pixelPairsSmall = os.path.join(tempFolder, 'stereoPixelPairsSmall.csv')
        if not os.path.exists(pixelPairsSmall):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/stereoPixelPairExtractor -i ' + disparityImage + ' -o ' + pixelPairsSmall + ' -p 100'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping small pair extraction step'

        # Compute the global rotation and offet between the two LE cubes
        globalParamsFile = os.path.join(tempFolder, 'solvedParams.csv')
        smallGdcFile     = os.path.join(tempFolder, 'gdcPointsSmall.csv')
        if not os.path.exists(globalParamsFile):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath ' + globalParamsFile + ' --gdcPointsOutPath ' + smallGdcFile + ' --matchingPixelsPath ' + pixelPairsSmall + ' ' + options.leftPath + ' ' + options.rightPath + ' --worldTransform --includePosition'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping stereo transform calculation step'

        # Read in the final state parameters (global rotation and shift)
        initialParamString = readSolvedState(globalParamsFile)
        if not initialParamString:
            return 1

        # Extract a large number of matching pixel locations (several thousand)
        pixelPairsLarge = os.path.join(tempFolder, 'stereoPixelPairsLarge.csv')
        if not os.path.exists(pixelPairsLarge):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/stereoPixelPairExtractor -i ' + disparityImage + ' -o ' + pixelPairsLarge + ' -p 5'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large pair extraction step'

        # Compute the 3d coordinates for each pixel pair
        largeGdcFile = os.path.join(tempFolder, 'gdcPointsLarge.csv')
        if not os.path.exists(largeGdcFile):
            cmd = '/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath ' + largeGdcFile + ' --matchingPixelsPath ' + pixelPairsLarge + ' ' + options.leftPath + ' ' + options.rightPath + ' --worldTransform --includePosition --initialOnly --initialValues ' + initialParamString
            print cmd
            os.system(cmd)
        else:
            print 'Skipping large GDC file creation step'

        # TODO: Apply rotation and offset to second pair (need absolute rot changer?)

        # Use pc-align to compare points to LOLA DEM, compute rotation and offset
        #largeGdcTransformedFile = os.path.joint(tempFolder, 'gdcPointsTransformedLarge.csv')
        transformedPointsFile = os.path.join(tempFolder, 'pcAlignOutput-trans_source.csv')        
        if not os.path.exists(transformedPointsFile):
            pcAlignOutputPrefix   = os.path.join(tempFolder, 'pcAlignOutput')
            cmd = 'pc_align --max-displacement 2000 --datum D_MOON --max-num-reference-points 25000000 --save-transformed-source-points ' + options.lolaPath + ' ' + largeGdcFile + ' -o ' + pcAlignOutputPrefix
            #--compute-translation-only'
            print cmd
            os.system(cmd)
        else:
            print 'Skipping pc_align step'

        #TODO: Generate final output file!

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
