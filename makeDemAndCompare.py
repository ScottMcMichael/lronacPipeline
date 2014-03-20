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

import sys, os, glob, optparse, re, subprocess, string, time, logging, threading

import IsisTools, pointErrorToKml

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two input mosaics and optionally compares to ASU and LOLA DEMs.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

MAP_PROJECT_METERS_PER_PIXEL  = 0.5 # Map resolution in meters
MAX_VALID_TRIANGULATION_ERROR = 15 # Meters
DEM_METERS_PER_PIXEL          = 1.0

DEM_NODATA    = -32767
MOON_RADIUS   = 1737400
SUBPIXEL_MODE = 1 # 1 = fast, 2 = accurate

ACCURATE_PIXEL_LIMIT = 5



# Generate a comparison of a DEM to the LOLA point cloud
def compareDemToLola(lolaPath, demPath, outputPath, csvPath, force):

    # Generate LOLA comparison files
    if force or not os.path.exists(outputPath):
        cmd = ('lola_compare --absolute --limit-hist=2 ' + demPath +  ' "' + lolaPath + 
                             '" -o ' + outputPath + ' -c ' + csvPath)
        print cmd
        os.system(cmd)

    # Make KML plot of the point errors
    kmlPath = os.path.splitext(csvPath)[0] + '.kml'
    if not os.path.exists(kmlPath):
        demName = os.path.dirname(demPath)
        cmdArgs = [csvPath, kmlPath, '--name', demName, '--skip', 4, '--errorLimit', 10]
        pointErrorToKml.main(cmdArgs)
    
    return True


# Generate a map projected version of an input image using the output DEM
def mapProjectImage(inputImage, demPath, outputPath, resolution, centerLat, force):
    
    if (not os.path.exists(outputPath)) or force:
        cmd = ('parallel_mapproject.py ' + demPath         + ' ' + inputImage + ' ' + outputPath +
                                       ' --tr ' + str(resolution) + ' --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + 
                                       ' +lat_0=0 +a='+str(MOON_RADIUS)+' +b='+str(MOON_RADIUS)+' +units=m" --nodata ' + str(DEM_NODATA) +
                                       ' --suppress-output --num-threads 24')
        
        print cmd
        os.system(cmd)
    else:
        print 'Map projected file  ' + outputPath + ' already exists, skipping map project step.'

    return True

# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('lola_compare')
    IsisTools.checkIfToolExists('parallel_stereo')
    IsisTools.checkIfToolExists('parallel_mapproject.py')
    IsisTools.checkIfToolExists('point2dem')
    IsisTools.checkIfToolExists('hillshade')
    IsisTools.checkIfToolExists('crop')

    return True

#--------------------------------------------------------------------------------

def main(argsIn):

    print '#################################################################################'
    print "Running makeDemAndCompare.py"

    try:
        try:
            usage = "usage: makeDemAndCompare.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
            inputGroup = optparse.OptionGroup(parser, 'Input Paths')
            inputGroup.add_option("--left",  dest="leftPath",  help="Path to left  cube file")
            inputGroup.add_option("--right", dest="rightPath", help="Path to right cube file")            


            inputGroup.add_option("--lola",    dest="lolaPath", help="Path to LOLA DEM")
            inputGroup.add_option("--asu",     dest="asuPath",  help="Path to ASU DEM")
            
            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--prefix",  dest="prefix",   help="Output prefix.")

            parser.add_option("--log-path",  dest="logPath",        
                              help="Where to write the output log file.")

            parser.add_option("--crop",  dest="cropAmount", 
                              help="Crops the output image to reduce processing time.")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args(argsIn)

            if not options.leftPath: 
                parser.error("Need left input path")
            if not options.rightPath: 
                parser.error("Need right input path")
            if not options.prefix: 
                parser.error("Need output prefix")            

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

        # Set up the output folders
        outputFolder  = os.path.dirname(options.prefix)
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
            options.logPath = options.prefix + '-Log.txt'
        logging.basicConfig(filename=options.logPath,level=logging.INFO)


        # If specified, crop the inputs that will be passed into the stereo function to reduce processing time
        mainMosaicCroppedPath   = os.path.join(tempFolder, 'mainMosaicCropped.cub')
        stereoMosaicCroppedPath = os.path.join(tempFolder, 'stereoMosaicCropped.cub')
        if options.cropAmount and (options.cropAmount > 0):
            if (not os.path.exists(mainMosaicCroppedPath)) or carry:
                cmd = ('crop from= ' + options.leftPath   + ' to= ' + mainMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount) )#+ ' line=24200')
                print cmd
                os.system(cmd)
            if (not os.path.exists(stereoMosaicCroppedPath) or carry):
                cmd = ('crop from= ' + options.rightPath + ' to= ' + stereoMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount) )# + ' line=24200')
                print cmd
                os.system(cmd)
            options.leftPath  = mainMosaicCroppedPath
            options.rightPath = stereoMosaicCroppedPath

        print '\n-------------------------------------------------------------------------\n'

        # Call stereo to generate a point cloud from the two images
        # - This step takes a really long time.

        stereoOutputPrefix = os.path.join(tempFolder, 'stereoWorkDir/stereo')
        stereoOutputFolder = os.path.join(tempFolder, 'stereoWorkDir')
        pointCloudPath     = stereoOutputPrefix + '-PC.tif'
        
        stereoOptionString = ('--corr-timeout 400 --alignment-method AffineEpipolar --subpixel-mode ' + str(SUBPIXEL_MODE) + 
                              ' ' + options.leftPath + ' ' + options.rightPath + 
                              ' ' + stereoOutputPrefix + ' --processes 8 --threads-multiprocess 4' +
                              ' --threads-singleprocess 32 --compute-error-vector' + ' --filter-mode 1' +
                              ' --erode-max-size 5000 --max-valid-triangulation-error ' + str(MAX_VALID_TRIANGULATION_ERROR))
  
        if (not os.path.exists(pointCloudPath)) or carry:
            cmd = ('parallel_stereo ' + stereoOptionString)
            print cmd
            os.system(cmd)
            
            ## Preprocessing, correlation, and subpixel refinement
            #cmd = ('parallel_stereo --entry-point 0 --stop-point 3' + stereoOptionString)
            #print cmd
            #os.system(cmd)
            #
            ## TODO: Delete -D files
            #
            ## Filtering, triangulation, and post-processing
            #cmd = ('parallel_stereo --entry-point 3 --stop-point 6' + stereoOptionString)
            #print cmd
            #os.system(cmd)
            #
            ## TODO: Delete other stereo files!
            #if not options.keep:
            #    IsisTools.removeIntermediateStereoFiles(stereoOutputPrefix) # Limited clear
            #    #IsisTools.removeFolderIfExists(stereoOutputFolder) # Larger clear
            
            
        else:
            print 'Stereo file ' + pointCloudPath + ' already exists, skipping stereo step.'

# TODO: Split this up and do piecemeal deletions (stages 0-5)


        stereoTime = time.time()
        logging.info('Stereo finished in %f seconds', stereoTime - startTime)

        # Compute percentage of good pixels
        percentGood = IsisTools.getStereoGoodPixelPercentage(stereoOutputPrefix)
        print 'Stereo completed with good pixel percentage: ' + str(percentGood)
        logging.info('Final stereo completed with good pixel percentage: %s', str(percentGood))


        # Find out the center latitude of the mosaic
        centerLat = IsisTools.getCubeCenterLatitude(options.leftPath, tempFolder)

        # Go ahead and set up all the output paths
        # -- Deliverables
        demPath               = options.prefix + '-DEM.tif'
        intersectionErrorPath = options.prefix + '-IntersectionErr.tif'
        hillshadePath         = options.prefix + '-Hillshade.tif'
        colormapPath          = options.prefix + '-Colormap.tif'
        mapProjectLeftPath    = options.prefix + '-MapProjLeft.tif'
        mapProjectRightPath   = options.prefix + '-MapProjRight.tif'
        # -- Diagnostic
        accuratePixelMask     = options.prefix + '-AccuratePixelMask.tif'
        intersectionViewPathX = options.prefix + '-IntersectionErrorX.tif'
        intersectionViewPathY = options.prefix + '-IntersectionErrorY.tif'
        intersectionViewPathZ = options.prefix + '-IntersectionErrorZ.tif'
        lolaDiffStatsPath     = options.prefix + '-LOLA_diff_stats.txt'
        lolaDiffPointsPath    = options.prefix + '-LOLA_diff_points.csv'
        lolaAsuDiffStatsPath  = options.prefix + '-ASU_LOLA_diff_stats.txt'
        lolaAsuDiffPointsPath = options.prefix + '-ASU_LOLA_diff_points.csv'


        # Generate a DEM
        if (not os.path.exists(demPath)) or carry:
            cmd = ('point2dem --errorimage -o ' + options.prefix + ' ' + pointCloudPath + 
                            ' -r moon --tr ' + str(DEM_METERS_PER_PIXEL) + ' --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + 
                            ' +lat_0=0 +a='+str(MOON_RADIUS)+' +b='+str(MOON_RADIUS)+' +units=m" --nodata ' + str(DEM_NODATA))
            print cmd
            os.system(cmd)
        else:
            print 'DEM file ' + pointCloudPath + ' already exists, skipping point2dem step.'

        # Create a hillshade image to visualize the output
        if (not os.path.exists(hillshadePath)) or carry:
            cmd = 'hillshade ' + demPath + ' -o ' + hillshadePath
            print cmd
            os.system(cmd)
        else:
            print 'Output file ' + hillshadePath + ' already exists, skipping hillshade step.'

        # Create a colorized version of the hillshade
        # - Uses a blue-red color map from here: http://www.sandia.gov/~kmorel/documents/ColorMaps/
        # - TODO: Modify colormap so we can save the legend file to the output directory!
        if (not os.path.exists(colormapPath)) or carry:
            lutFilePath = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'colorProfileBlueRed.csv')
            cmd = 'colormap ' + demPath + ' -o ' + colormapPath + ' -s ' + hillshadePath + ' --lut-file ' + lutFilePath
            print cmd
            os.system(cmd)
        else:
            print 'Output file ' + colormapPath + ' already exists, skipping colormap step.'


        ## Create a 3d mesh of the point cloud
        #meshPath   = os.path.join(outputFolder, 'mesh.ive')
        #meshPrefix = os.path.join(outputFolder, 'mesh')
        #cmd = 'point2mesh ' + pointCloudPath + ' ' + options.leftPath + ' -o ' + meshPrefix
        #if not os.path.exists(meshPath):
        #    print cmd
        #    os.system(cmd)

        # Convert the intersection error to a viewable format
        cmdX = 'gdal_translate -ot byte -scale 0 10 0 255 -outsize 50% 50% -b 1 ' + intersectionErrorPath + ' ' + intersectionViewPathX
        cmdY = 'gdal_translate -ot byte -scale 0 10 0 255 -outsize 50% 50% -b 2 ' + intersectionErrorPath + ' ' + intersectionViewPathY
        cmdZ = 'gdal_translate -ot byte -scale 0 10 0 255 -outsize 50% 50% -b 3 ' + intersectionErrorPath + ' ' + intersectionViewPathZ
        if not os.path.exists(intersectionViewPathX) or carry:
            print cmdX
            os.system(cmdX)
        if not os.path.exists(intersectionViewPathY) or carry:
            print cmdY
            os.system(cmdY)
        if not os.path.exists(intersectionViewPathZ) or carry:
            print cmdZ
            os.system(cmdZ)

        # Generate a good pixel mask from the intersection error (all pixels below limit in white)
        if not os.path.exists(accuratePixelMask):
            cmd = ('maskFromIntersectError ' + intersectionErrorPath + ' ' + accuratePixelMask +
                                         ' --scaleOutput --thresholds ' + str(ACCURATE_PIXEL_LIMIT))
            print cmd
            os.system(cmd)

        hillshadeTime = time.time()
        logging.info('DEM and hillshade finished in %f seconds', hillshadeTime - stereoTime)


        # Call script to compare LOLA data with the DEM
        if options.lolaPath:
            compareDemToLola(options.lolaPath, demPath, lolaDiffStatsPath, lolaDiffPointsPath, carry)

        # Call script to compare LOLA data with the ASU DEM
        if options.asuPath:
            compareDemToLola(options.lolaPath, options.asuPath, lolaAsuDiffStatsPath, lolaAsuDiffPointsPath, carry)


        # Generate a map projected version of the left and right images
        # - This step is done last since it is so slow!
        #mapProjectImage(options.leftPath,  demPath, mapProjectLeftPath,  MAP_PROJECT_METERS_PER_PIXEL, centerLat, carry)
        #mapProjectImage(options.rightPath, demPath, mapProjectRightPath, MAP_PROJECT_METERS_PER_PIXEL, centerLat, carry)

        mapProjectTime = time.time()
        logging.info('Map project finished in %f seconds', mapProjectTime - hillshadeTime)

        # Clean up temporary files
        if not options.keep:
            print 'Removing temporary files'
            IsisTools.removeIfExists(mainMosaicCroppedPath)
            IsisTools.removeIfExists(stereoMosaicCroppedPath)
            #IsisTools.removeIntermediateStereoFiles(stereoOutputPrefix) # Limited clear
            IsisTools.removeFolderIfExists(stereoOutputFolder) # Larger clear
            #if (hadToCreateTempFolder): Not done since stereo output needs to be retained
            #    IsisTools.removeFolderIfExists(tempFolder)


        endTime = time.time()

        logging.info('makeDemAndCompare.py finished in %f seconds', endTime - startTime)
        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
