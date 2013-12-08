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

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two LRONAC pairs, trying to use LOLA data for increased accuracy.
'''
    sys.exit()



#--------------------------------------------------------------------------------

# Generates a KML file to describe a set of GDC points on the moon
def generateKmlFromGdcPoints(inputFolder, outputFolder, filename, pointSkip, color, forceOperation):

    # Determine input and output paths
    inputPath      = os.path.join(inputFolder,  filename)
    outputFilename = os.path.splitext(filename)[0] + '.kml'
    outputPath     = os.path.join(outputFolder, outputFilename)
    kmlName        = os.path.splitext(filename)[0]


    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputPath)):
        print 'File ' + outputPath + ' already exists, skipping kml generation.'
        return True

    # Generate the new file
    cmd = 'calibrationReport.py --input ' + inputPath + ' --output ' + outputFilename + ' --name ' + kmlName +  ' --skip ' + str(pointSkip) + ' --color ' + color
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputPath):
        raise Exception('Kml point generation failed to create output file ' + outputPath + ' from input file ' + inputPath)

    return True


# Calls the ISIS noproj function on a cube
def noprojCubePair(inputCube, outputCube, matchCube, pvlPath, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCube)):
        print 'File ' + outputCube + ' already exists, skipping noproj.'
        return True

    # Generate the new file
    cmd = 'noproj from= ' + inputCube + ' to= ' + outputCube + ' match= ' + matchCube + ' specs= ' + pvlPath
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCube):
        raise Exception('noproj failed to create output file ' + outputCube + ' from input file ' + inputCube)

    return True


# Creates a mosaic from two noproj'd input cubes.
def createMosaic(leftCube, rightCube, outputCube, workDir, forceOperation):

    # Quit immediately if the output file already exists
    if (not forceOperation) and (os.path.exists(outputCube)):
        print 'File ' + outputCube + ' already exists, skipping mosaic creation.'
        return True

    # Call lronacjitreg to determine any remaining offset between the two input files
    jitRegOutputPath = os.path.join(workDir, 'jitregResults.txt')
    cmd = 'lronacjitreg --correlator-type 2 --kernel 15 15 --output-log ' + jitRegOutputPath + ' ' + leftCube + ' ' + rightCube;
    print cmd
    os.system(cmd)

    # Read in the output from lronacjitreg
    jitRegOffsets = IsisTools.readJitregFile(jitRegOutputPath)

    # Set intermediate mosaic file path and start on mosaic
    mosaicCube = os.path.join(workDir, 'mosaic.cub')
    cmd = 'cp ' + leftCube + ' ' + mosaicCube
    print cmd
    os.system(cmd)

    # TODO: Try out more advanced ISIS mosaic merging functions!
    # Create the mosaic, applying offsets from jitreg (converting into handmos conventions)
    cmd = 'handmos from= ' + rightCube + ' mosaic= ' + mosaicCube + ' outsample= ' + str(1 - jitRegOffsets(0)) + ' outline= ' + str(1 - jitRegOffsets(1)) + ' matchbandbin=FALSE priority=ontop'
    print cmd
    os.system(cmd)

    # Call cubenorm to improve the mosaic appearance
    cmd = 'cubenorm from= ' + mosaicCube + ' to= ' + outputCube
    print cmd
    os.system(cmd)

    # Check to make sure we actually created the file
    if not os.path.exists(outputCube):
        raise Exception('Failed to create mosaic file ' + outputCube + ' from input files ' + leftCube + ' and ' + rightCube)

    return True

#--------------------------------------------------------------------------------

def main():

    print "lronac2dem.py"

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

            parser.add_option("--prefix",  dest="prefix",        help="Output prefix.")

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

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

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


        # Set up output paths for the stereo calibration call
        leftCorrectedPath        = os.path.join(tempFolder, 'leftFinalCorrected.cub')
        rightCorrectedPath       = os.path.join(tempFolder, 'rightFinalCorrected.cub')
        leftStereoCorrectedPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.cub')
        rightStereoCorrectedPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.cub')
        correct = os.path.join(tempFolder, 'rightStereoFinalCorrected.cub')

        # Correct all four input images at once
        if not os.path.exists(leftCorrectedPath):
            cmd = '~/repo/lronacPipeline/stereoDoubleCalibrationProcess.py --left ' + options.leftPath + ' --right ' +  options.rightPath + ' --stereo-left ' + options.stereoLeft + ' --stereo-right ' + options.stereoRight + ' --lola ' + options.lolaPath + ' --keep --outputL ' + leftCorrectedPath + ' --outputR ' + rightCorrectedPath + ' --outputSL ' + leftStereoCorrectedPath + ' --outputSR ' + rightStereoCorrectedPath + ' --workDir ' + options.workDir
            print cmd
            os.system(cmd)

        # Convert GDC output files into KML plots 
        # - This is just to help with debugging
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'pairGdcCheckInitial.csv',          1,      'blue', False)
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'pairGdcCheckFinal.csv',            1,      'red',  False)
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'gdcPointsCheckFinalRot.csv',       1,      'blue', False)
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'gdcPointsCheckFinalRotStereo.csv', 1,      'blue', False)
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'gdcPointsLarge.csv',               100000, 'blue', False)
        #generateKmlFromGdcPoints(options.workDir, tempFolder, 'dem-trans_source.csv',             'blue', False)
        generateKmlFromGdcPoints(options.workDir, tempFolder, 'dem-trans_reference.csv',          100000, 'red',  False)

        # Generate a PVL file that we need for noproj
        pvlPath   = os.path.join(tempFolder, 'noprojInstruments_fullRes.pvl')
        isHalfRes = False # TODO: Check to see if this is true!
        IsisTools.writeLronacPvlFile(pvlPath, isHalfRes)


        # Noproj the corrected data

        leftNoprojPath        = os.path.join(tempFolder, 'leftFinalCorrected.noproj.cub')
        rightNoprojPath       = os.path.join(tempFolder, 'rightFinalCorrected.noproj.cub')
        leftStereoNoprojPath  = os.path.join(tempFolder, 'leftStereoFinalCorrected.noproj.cub')
        rightStereoNoprojPath = os.path.join(tempFolder, 'rightStereoFinalCorrected.noproj.cub')

        noprojCubePair(leftCorrectedPath,        leftNoprojPath,        leftCorrectedPath,       pvlPath, False)
        noprojCubePair(rightCorrectedPath,       rightNoprojPath,       leftCorrectedPath,       pvlPath, False)
        noprojCubePair(rightCorrectedPath,       leftStereoNoprojPath,  leftStereoCorrectedPath, pvlPath, False)
        noprojCubePair(rightStereoCorrectedPath, rightStereoNoprojPath, leftStereoCorrectedPath, pvlPath, False)


        # Combine the noproj files to make a mosaic.
        # - This also takes care of the cubenorm step.
        # - This step takes a while.

        mainMosaicPath      = os.path.join(tempFolder, 'mainMosaic.cub')
        stereoMosaicPath    = os.path.join(tempFolder, 'stereoMosaic.cub')
        mainMosaicWorkDir   = os.path.join(tempFolder, 'mainMosaicWorkDir/')
        stereoMosaicWorkDir = os.path.join(tempFolder, 'stereoMosaicWorkDir/')

        createMosaic(leftNoprojPath,       rightNoprojPath,       mainMosaicPath,   mainMosaicWorkDir,   forceOperation)
        createMosaic(leftStereoNoprojPath, rightStereoNoprojPath, stereoMosaicPath, stereoMosaicWorkDir, forceOperation)


        # Call stereo to generate a point cloud from the two images
        # - This step takes a really long time.

        stereoOutputPrefix = os.path.join(tempFolder, 'stereoWorkDir/stereo')
        pointCloudPath     = stereoOutputPrefix + '-PC.tif'
        if not os.path.exists(pointCloudPath):
            cmd = 'parallel_stereo --corr-timeout 400 --alignment affineepipolar --subpixel-mode 1 --disable-fill-holes ' +  mainMosaicPath + ' ' + stereoMosaicPath + ' ' + stereoWorkDir + ' --processes 8 --threads-multiprocess 4 --threads-singleprocess 32 --compute-error-vector'
            print cmd
            os.system(cmd)
            #--nodes-list PBS_NODEFILE --processes 4 --threads-multiprocess 16 --threads-singleprocess 32
        else:
            print 'Stereo file ' + pointCloudPath + ' already exists, skipping stereo step.'

        # Find out the center latitude of the mosaic
        centerLat = IsisTools.getCubeCenterLatitude(mainMosaicPath, tempFolder)

        # Generate a DEM
        demPath = os.path.join(outputFolder, 'outputDEM.tif')
        if not os.path.exists(demPath):
            cmd = 'point2dem --errorimage -o ' + demPath + ' ' + pointCloudPath + ' -r moon --tr 1 --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + ' +lat_0=0 +a=1737400 +b=1737400 +units=m" --nodata -32767'
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + pointCloudPath + ' already exists, skipping point2dem step.'

        # Create a hillshade image to check if the central errors are gone
        hillshadePath = os.path.join(outputFolder, 'outputHillshade.tif')
        if not os.path.exists(demPath):
            cmd = 'hillshade ' + demPath + ' -o ' + hillshadePath
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + hillshadePath + ' already exists, skipping hillshade step.'






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





