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

import IrgIsisFunctions, IrgAspFunctions, IrgFileFunctions, IrgMathFunctions
import IrgGeoFunctions, pointErrorToKml

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

# For convenience some of the constants used in the file are set up here

MAP_PROJECT_METERS_PER_PIXEL  = 0.5 # Map resolution in meters
DEM_METERS_PER_PIXEL          = 1.0

DEM_NODATA    = -32767
MOON_RADIUS   = 1737400
SUBPIXEL_MODE = 3 # 1 = fast, 2 = accurate, 3 = compromise

PIXEL_ACCURACY_THRESHOLDS = [4, 16, 64]



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
def mapProjectImage(inputImage, demPath, outputPath, resolution, centerLat, nodeFilePath, force):

    if (not os.path.exists(outputPath)) or force:
        cmd = ('parallel_mapproject.py ' + demPath         + ' ' + inputImage + ' ' + outputPath +
                                       ' --tr ' + str(resolution) + ' --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + 
                                       ' +lat_0=0 +a='+str(MOON_RADIUS)+' +b='+str(MOON_RADIUS)+' +units=m" --nodata ' + str(DEM_NODATA) +
                                       ' --suppress-output --tile-size 4096')
        
	if nodeFilePath: # Needed to operate across multiple computers
	    cmd = cmd + ' --node-file ' + nodeFilePath
        print cmd
        os.system(cmd)

    else:
        print 'Map projected file  ' + outputPath + ' already exists, skipping map project step.'

    if not os.path.exists(outputPath):
        raise Exception('Failed to create map projected image ' + outputPath)

    return True


def writeColorMapInfo(colormapPath, lutFilePath, demPath, outputPath):
    """Writes a file containing the color map information"""
    
    colormapPercentiles = []
    numColorSteps = IrgFileFunctions.getFileLineCount(lutFilePath)
    for i in range(0,numColorSteps): # This loop generates the percentage values from the color profile we are using
        colormapPercentiles.append(i / float(numColorSteps-1))
    # Get the min and max elevation of the DEM
    elevationBounds = IrgGeoFunctions.getImageStats(demPath) 
    # Get elevation values for each percentile
    colormapPercentileValues = IrgMathFunctions.getPercentileValues(elevationBounds[0][0], elevationBounds[0][1], colormapPercentiles)
    
    # Now write out a version of the LUT file with values added
    inputFile  = open(lutFilePath, 'r')
    outputFile = open(outputPath,  'w')

    # Add a header line to the output file
    outputFile.write('Percent of range, R, G, B, Elevation (meters above datum)\n')

    # Write a copy of the input file with the elevation values appended to each line    
    for (line, value) in zip(inputFile, colormapPercentileValues):
        newLine = line[:-1] + ', ' + str(round(value,2)) + '\n'
        outputFile.write(newLine) 
        
    inputFile.close()
    outputFile.close()
    



# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IrgFileFunctions.checkIfToolExists('lola_compare')
    IrgFileFunctions.checkIfToolExists('parallel_stereo')
    IrgFileFunctions.checkIfToolExists('parallel_mapproject.py')
    IrgFileFunctions.checkIfToolExists('point2dem')
    IrgFileFunctions.checkIfToolExists('hillshade')
    IrgFileFunctions.checkIfToolExists('colormap')
    IrgFileFunctions.checkIfToolExists('crop')
    IrgFileFunctions.checkIfToolExists('tar')
    IrgFileFunctions.checkIfToolExists('gdal_translate')
    IrgFileFunctions.checkIfToolExists('maskFromIntersectError')

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
            
            inputGroup.add_option("--node-file", dest="nodeFilePath", 
                                  help="Path to file containing list of available nodes")

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


        # Go ahead and set up all the output paths
        # -- Deliverables
        demPath               = options.prefix + '-DEM.tif'
        intersectionErrorPath = options.prefix + '-IntersectionErr.tif'
        hillshadePath         = options.prefix + '-Hillshade.tif'
        colormapPath          = options.prefix + '-Colormap.tif'
        colormapLegendPath    = options.prefix + '-ColormapLegend.csv'
        mapProjectLeftPath    = options.prefix + '-MapProjLeft.tif'
        mapProjectRightPath   = options.prefix + '-MapProjRight.tif'
        confidenceLevelPath   = options.prefix + '-Confidence.tif'
        confidenceLegendPath  = options.prefix + '-ConfidenceLegend.csv'
        # -- Diagnostic
        intersectionViewPathX    = options.prefix + '-IntersectionErrorX.tif'
        intersectionViewPathY    = options.prefix + '-IntersectionErrorY.tif'
        intersectionViewPathZ    = options.prefix + '-IntersectionErrorZ.tif'
        lolaDiffStatsPath        = options.prefix + '-LOLA_diff_stats.txt'
        lolaDiffPointsPath       = options.prefix + '-LOLA_diff_points.csv'
        lolaAsuDiffStatsPath     = options.prefix + '-ASU_LOLA_diff_stats.txt'
        lolaAsuDiffPointsPath    = options.prefix + '-ASU_LOLA_diff_points.csv'
        mapProjectLeftUint8Path  = options.prefix + '-MapProjLeftUint8.tif'
        mapProjectRightUint8Path = options.prefix + '-MapProjRightUint8.tif'


        # If specified, crop the inputs that will be passed into the stereo function to reduce processing time
        mainMosaicCroppedPath   = os.path.join(tempFolder, 'mainMosaicCropped.cub')
        stereoMosaicCroppedPath = os.path.join(tempFolder, 'stereoMosaicCropped.cub')
        if options.cropAmount and (options.cropAmount > 0):

            # Verify input files are present
            if not os.path.exists(options.leftPath):
                raise Exception('Input file ' + options.leftPath + ' not found!')
            if not os.path.exists(options.rightPath):
                raise Exception('Input file ' + options.rightPath + ' not found!')

            if (not os.path.exists(mainMosaicCroppedPath)) or carry:
                cmd = ('crop from= ' + options.leftPath   + ' to= ' + mainMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))# + ' line=24200')
                print cmd
                os.system(cmd)
            if (not os.path.exists(stereoMosaicCroppedPath) or carry):
                cmd = ('crop from= ' + options.rightPath + ' to= ' + stereoMosaicCroppedPath + 
                           ' nlines= ' + str(options.cropAmount))#  + ' line=24200')
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
                              ' --job-size-w 4096 --job-size-h 4096 ' + # Reduce number of tile files created
                              ' ' + stereoOutputPrefix + ' --processes 10 --threads-multiprocess 4' +
                             ' --threads-singleprocess 32 --compute-error-vector' + ' --filter-mode 1' +
                              ' --erode-max-size 5000 --subpixel-kernel 35 35 --subpixel-max-levels 0')
  
        if (not os.path.exists(pointCloudPath) and not os.path.exists(demPath)) or carry:

            # Verify input files are present
            if not os.path.exists(options.leftPath):
                raise Exception('Input file ' + options.leftPath + ' not found!')
            if not os.path.exists(options.rightPath):
                raise Exception('Input file ' + options.rightPath + ' not found!')

            cmd = ('parallel_stereo ' + stereoOptionString)
            print cmd
            os.system(cmd)
 

            # Compute percentage of good pixels
            percentGood = IrgAspFunctions.getStereoGoodPixelPercentage(stereoOutputPrefix)
            print 'Stereo completed with good pixel percentage: ' + str(percentGood)
            logging.info('Final stereo completed with good pixel percentage: %s', str(percentGood))
                      
        else:
            print 'Stereo file ' + pointCloudPath + ' already exists, skipping stereo step.'


        stereoTime = time.time()
        logging.info('Stereo finished in %f seconds', stereoTime - startTime)


        # Find out the center latitude of the mosaic
        if os.path.exists(options.leftPath):
            centerLat = IrgIsisFunctions.getCubeCenterLatitude(options.leftPath, tempFolder)
        elif os.path.exists(demPath): # Input file has been deleted but we still have the info
            demInfo = IrgGeoFunctions.getImageGeoInfo(demPath, False)
            centerLat = demInfo['standard_parallel_1']
        else:
            raise Exception("Can't delete the input files before creating the DEM!")

        # Generate a DEM
        if (not os.path.exists(demPath)) or carry:
            
            #TODO: Add a +lon_0 entry here to set the central meridian?
            
            # Equirectangular style projection
            # - Latitude of true scale = center latitude = lat_ts
            # - Latitude of origin = 0 = lat+0
            # - Longitude of projection center = Central meridian = lon+0
            cmd = ('point2dem --dem-hole-fill-len 15  --remove-outliers --errorimage -o ' + options.prefix + ' ' + pointCloudPath + 
                            ' -r moon --tr ' + str(DEM_METERS_PER_PIXEL) + ' --t_srs "+proj=eqc +lat_ts=' + str(centerLat) + 
                            ' +lat_0=0 +a='+str(MOON_RADIUS)+' +b='+str(MOON_RADIUS)+' +units=m" --nodata ' + str(DEM_NODATA))
            os.system(cmd)
        else:
            print 'DEM file ' + demPath + ' already exists, skipping point2dem step.'

        # Create a hillshade image to visualize the output
        if (not os.path.exists(hillshadePath)) or carry:
            cmd = 'hillshade ' + demPath + ' -o ' + hillshadePath
            print cmd
            os.system(cmd)
        else:
            print 'Output file ' + hillshadePath + ' already exists, skipping hillshade step.'

        # Create a colorized version of the hillshade
        # - Uses a blue-red color map from here: http://www.sandia.gov/~kmorel/documents/ColorMaps/
        if (not os.path.exists(colormapPath)) or (not os.path.exists(colormapLegendPath)) or carry:
            # The color LUT is kept with the source code
            lutFilePath = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'colorProfileBlueRed.csv')
            
            
            # Generate the initial version of colormap
            colormapTempPath = options.prefix + '-ColormapTemp.tif'
            cmd = 'colormap ' + demPath + ' -o ' + colormapTempPath + ' -s ' + hillshadePath + ' --lut-file ' + lutFilePath
            print cmd
            os.system(cmd)
            
            # Now convert to the final output version (remove transparency layer) and remove the temp file
            IrgFileFunctions.stripRgbImageAlphaChannel(colormapTempPath, colormapPath)            
            os.remove(colormapTempPath)
            
            # Generate another file storing the colormap info
            writeColorMapInfo(colormapPath, lutFilePath, demPath, colormapLegendPath)
        else:
            print 'Output file ' + colormapPath + ' already exists, skipping colormap step.'


        ## Create a 3d mesh of the point cloud
        #meshPath   = os.path.join(outputFolder, 'mesh.ive')
        #meshPrefix = os.path.join(outputFolder, 'mesh')
        #cmd = 'point2mesh ' + pointCloudPath + ' ' + options.leftPath + ' -o ' + meshPrefix
        #if not os.path.exists(meshPath):
        #    print cmd
        #    os.system(cmd)

        if not options.keep: # Remove stereo folder here to cut down on file count before mapproject calls
             IrgFileFunctions.removeFolderIfExists(stereoOutputFolder)


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

        # Generate a confidence plot from the intersection error
        if not os.path.exists(confidenceLevelPath):
            thresholdString = str(PIXEL_ACCURACY_THRESHOLDS)[1:-1].replace(',', '')
            cmd = ('maskFromIntersectError ' + intersectionErrorPath + ' ' + confidenceLevelPath
                                             + ' --legend ' + confidenceLegendPath + 
                                        ' --scaleOutput --thresholds ' + thresholdString)
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
        mapProjectImage(options.leftPath,  demPath, mapProjectLeftPath,  MAP_PROJECT_METERS_PER_PIXEL, centerLat, options.nodeFilePath, carry)
        mapProjectImage(options.rightPath, demPath, mapProjectRightPath, MAP_PROJECT_METERS_PER_PIXEL, centerLat, options.nodeFilePath, carry)

        # Generate 8 bit versions of the mapproject files for debugging
        cmdLeft  = 'gdal_translate -scale -ot byte ' + mapProjectLeftPath  + ' ' + mapProjectLeftUint8Path
        cmdRight = 'gdal_translate -scale -ot byte ' + mapProjectRightPath + ' ' + mapProjectRightUint8Path
        if not os.path.exists(mapProjectLeftUint8Path) or carry:
            print cmdLeft
            os.system(cmdLeft)
        if not os.path.exists(mapProjectRightUint8Path) or carry:
            print cmdRight
            os.system(cmdRight)

        mapProjectTime = time.time()
        logging.info('Map project finished in %f seconds', mapProjectTime - hillshadeTime)

        # Clean up temporary files
        if not options.keep:
            print 'Removing temporary files'
            IrgFileFunctions.removeIfExists(mainMosaicCroppedPath)
            IrgFileFunctions.removeIfExists(stereoMosaicCroppedPath)
            #IrgFileFunctions.removeIntermediateStereoFiles(stereoOutputPrefix) # Limited clear
            IrgFileFunctions.removeFolderIfExists(stereoOutputFolder) # Larger clear
            #if (hadToCreateTempFolder): Not done since stereo output needs to be retained
            #    IrgFileFunctions.removeFolderIfExists(tempFolder)


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
