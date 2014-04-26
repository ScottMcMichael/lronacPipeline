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

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging, threading
import simplekml, matplotlib

import IsisTools, IrgIsisFunctions, IrgGeoFunctions, IrgFileFunctions, IrgMathFunctions, IrgAspFunctions

import lronac2refinedMosaics
import makeDemAndCompare


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a stereo DEM from two LRONAC pairs using SBA and LOLA for increased accuracy.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------

def writeLabelFile(imagePath, outputPath, dataSetName, versionId, description, extraData=None):
    """Write out a .LBL file formatted for the PDS"""
    
    # Call functions to automatically obtain some data from the referenced image
    imageSize   = IrgGeoFunctions.getImageSize(imagePath)
    boundingBox = IrgGeoFunctions.getImageBoundingBox(imagePath)    

    # Obtain the ASP version string
    aspVersionString = IrgAspFunctions.getAspVersionStrings()
    
    imageGeoInfo = IrgGeoFunctions.getImageGeoInfo(imagePath)

    projCenterLatitude  = imageGeoInfo['standard_parallel_1']
    projCenterLongitude = imageGeoInfo['central_meridian']
        
    # Currently assuming pixels are the same size
    metersPerPixel = abs(imageGeoInfo['pixel size'][0])
    
    # Compute pixels per degree
    lonSpanDegrees     = boundingBox[1] - boundingBox[0]
    latSpanDegrees     = boundingBox[3] - boundingBox[2]
    pixelsPerDegreeLon = imageSize[0] / lonSpanDegrees
    pixelsPerDegreeLat = imageSize[1] / latSpanDegrees
    pixelsPerDegree    = (pixelsPerDegreeLat + pixelsPerDegreeLon) / 2.0

    # Computed by dividing 'Origin' by 'Pixel Size'
    lineProjOffset   = imageGeoInfo['origin'][0] / imageGeoInfo['pixel size'][1]
    sampleProjOffset = imageGeoInfo['origin'][1] / imageGeoInfo['pixel size'][0]
        
    
    labelFile = open(outputPath, 'w')
    
    labelFile.write('PDS_VERSION_ID            = PDS3\n')
    
    labelFile.write('/* The source image data definition. */\n')
    labelFile.write('^IMAGE        = ' + os.path.basename(outputPath) +'\n')
    
    labelFile.write('/* Identification Information  */\n')
    labelFile.write('DATA_SET_ID               = ""\n') # Someone will tell us what to put here
    labelFile.write('DATA_SET_NAME             = ""\n') # Someone will tell us what to put here
    labelFile.write("VOLUME_ID                 = ''\n") # Someone will tell us what to put here
    labelFile.write('PRODUCER_INSTITUTION_NAME = "NASA AMES RESEARCH CENTER"\n')
    labelFile.write('PRODUCER_ID               = NASA IRG\n')
    labelFile.write('PRODUCER_FULL_NAME        = "ZACHARY MORATTO"\n')
    labelFile.write("PRODUCT_ID                = " + dataSetName + "\n")
    labelFile.write("PRODUCT_VERSION_ID        = " + versionId + "\n")
    labelFile.write('PRODUCT_TYPE              = "RDR"\n')
    labelFile.write('INSTRUMENT_HOST_NAME      = "LUNAR RECONNAISSANCE ORBITER"\n')
    labelFile.write('INSTRUMENT_HOST_ID        = "LRO"\n')
    labelFile.write('INSTRUMENT_NAME           = "LUNAR RECONNAISSANCE ORBITER CAMERA"\n')
    labelFile.write('INSTRUMENT_ID             = "LROC"\n')
    labelFile.write('TARGET_NAME               = MOON\n')
    labelFile.write('MISSION_PHASE_NAME        = "NOMINAL MISSION"\n')
    labelFile.write("""RATIONALE_DESC            = "Created at the request of NASA's Exploration\n""") 
    labelFile.write('                            Systems Mission Directorate to support future\n')
    labelFile.write('                            human exploration"\n')
    labelFile.write('SOFTWARE_NAME             = "'+ aspVersionString[0] +' | '+ aspVersionString[2] +'"\n')
    labelFile.write('DESCRIPTION               = "' + description + '"\n')
    labelFile.write('\n')
    labelFile.write('/* Time Parameters */\n')
    labelFile.write('PRODUCT_CREATION_TIME        = ' + time.strftime("%Y-%m-%dT%H:%M:%S") + '\n')
    labelFile.write('\n')
    labelFile.write('/* NOTE:                                                                   */\n')
    labelFile.write('/* This raster image is composed of a set of pixels that represent finite  */\n')
    labelFile.write('/* areas, and not discrete points.  The center of the upper left pixel is  */\n')
    labelFile.write('/* defined as line and sample (1.0,1.0). The                               */\n')
    labelFile.write('/* [LINE,SAMPLE]_PROJECTION_OFFSET elements are the pixel offset from line */\n')
    labelFile.write('/* and sample (1.0,1.0) to the map projection origin (defined by the       */\n')
    labelFile.write('/* CENTER_LATITUDE and CENTER_LONGITUDE elements).  These offset values    */\n')
    labelFile.write('/* are positive when the map projection origin is to the right or below    */\n')
    labelFile.write('/* the center of the upper left pixel.                                     */\n')
    if extraData: # Location for additional notes
        labelFile.write(extraData)
    labelFile.write('\n')
    labelFile.write('OBJECT = IMAGE_MAP_PROJECTION\n')
    labelFile.write('    MAP_PROJECTION_TYPE          = EQUIRECTANGULAR\n') # Specified by +proj=eqc
    labelFile.write('    PROJECTION_LATITUDE_TYPE     = PLANETOCENTRIC\n')  #From gdalinfo?
    labelFile.write('    A_AXIS_RADIUS                = 1737.4 <KM>\n') # Fixed lunar radius
    labelFile.write('    B_AXIS_RADIUS                = 1737.4 <KM>\n')
    labelFile.write('    C_AXIS_RADIUS                = 1737.4 <KM>\n')
    labelFile.write('    COORDINATE_SYSTEM_NAME       = PLANETOCENTRIC\n') #From gdalinfo?
    labelFile.write('    POSITIVE_LONGITUDE_DIRECTION = EAST\n') #From gdalinfo?
    labelFile.write('    KEYWORD_LATITUDE_TYPE        = PLANETOCENTRIC\n') #From gdalinfo?
    labelFile.write('    /* NOTE:  CENTER_LATITUDE and CENTER_LONGITUDE describe the location   */\n')
    labelFile.write('    /* of the center of projection, which is not necessarily equal to the  */\n')
    labelFile.write('    /* location of the center point of the image.                          */\n')
    labelFile.write('    CENTER_LATITUDE              = ' + str(projCenterLatitude)  + ' <DEG>\n')
    labelFile.write('    CENTER_LONGITUDE             = ' + str(projCenterLongitude) + ' <DEG>\n')
    labelFile.write('    LINE_FIRST_PIXEL             = 1\n')
    labelFile.write('    LINE_LAST_PIXEL              = ' + str(imageSize[1] + 1) + '\n')
    labelFile.write('    SAMPLE_FIRST_PIXEL           = 1\n')
    labelFile.write('    SAMPLE_LAST_PIXEL            = ' + str(imageSize[0] + 1) + '\n')
    labelFile.write('    MAP_PROJECTION_ROTATION      = 0.0 <DEG>\n') #From gdalinfo (probably always zero)
    labelFile.write('    MAP_RESOLUTION               = ' + str(round(pixelsPerDegree,2)) +' <PIX/DEG>\n')
    labelFile.write('    MAP_SCALE                    = ' + str(round(metersPerPixel,4)) + ' <METERS/PIXEL>\n')
    labelFile.write('    MAXIMUM_LATITUDE             = ' + str(boundingBox[3]) + ' <DEG>\n') 
    labelFile.write('    MINIMUM_LATITUDE             = ' + str(boundingBox[2]) + ' <DEG>\n')
    labelFile.write('    EASTERNMOST_LONGITUDE        = ' + str(boundingBox[0]) + ' <DEG>\n')
    labelFile.write('    WESTERNMOST_LONGITUDE        = ' + str(boundingBox[1]) + ' <DEG>\n')
    labelFile.write('    LINE_PROJECTION_OFFSET       = ' + str(round(lineProjOffset,2))  +' <PIXEL>\n')
    labelFile.write('    SAMPLE_PROJECTION_OFFSET     = ' + str(round(sampleProjOffset,2)) +' <PIXEL>\n')
    labelFile.write('END_OBJECT = IMAGE_MAP_PROJECTION\n')
    labelFile.write('\n')
    labelFile.write('END\n')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    labelFile.close()
    return True

# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IrgFileFunctions.checkIfToolExists('lronac2refinedMosaics.py')
    IrgFileFunctions.checkIfToolExists('makeDemAndCompare.py')
    return True

#--------------------------------------------------------------------------------

def main():

    print '#################################################################################'
    print "Running lronac2dem.py"

    try:
        try:
            usage = "usage: lronac2dem.py [--output <path>][--manual]\n  "
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
            inputGroup.add_option("--asu",     dest="asuPath",  help="Path to ASU DEM")
            
            inputGroup.add_option("--node-file", dest="nodeFilePath", 
                                  help="Path to list of available computing nodes")

            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--outputPrefix",  dest="outputPrefix",   help="Output prefix.")

            parser.add_option("--crop",  dest="cropAmount", 
                              help="Crops the output image to reduce processing time.")

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
                parser.error("Need LOLA data path")
            if not options.outputPrefix: 
                parser.error("Need output prefix")

        except optparse.OptionError, msg:
            raise Usage(msg)

        startTime = time.time()

        # Set up the output folders
        outputFolder  = os.path.dirname(options.outputPrefix)
        inputBaseName = os.path.basename(options.leftPath)
        tempFolder    = outputFolder + '/' + inputBaseName + '_stereoCalibrationTemp/'
        if (options.workDir):
            tempFolder = options.workDir
        IrgFileFunctions.createFolder(outputFolder)
        hadToCreateTempFolder = not os.path.exists(tempFolder)
        IrgFileFunctions.createFolder(tempFolder)


        # Set up logging
        logPath = options.outputPrefix + '-Log.txt'
        logging.basicConfig(filename=logPath,level=logging.INFO)

        # These are the output mosaic paths
        filename         = os.path.splitext(options.leftPath)[0] + '.correctedMosaic.cub'
        mainMosaicPath   = os.path.join(tempFolder, os.path.basename(filename))
        filename         = os.path.splitext(options.stereoLeft)[0] + '.correctedMosaic.cub'
        stereoMosaicPath = os.path.join(tempFolder, os.path.basename(filename))


        # List of all the output files that will be created by makeDemAndCompare.py
        demPath               = options.outputPrefix + '-DEM.tif'
        intersectionErrorPath = options.outputPrefix + '-IntersectionErr.tif'
        hillshadePath         = options.outputPrefix + '-Hillshade.tif'
        colormapPath          = options.outputPrefix + '-Colormap.tif'
        colormapLegendPath    = options.outputPrefix + '-ColormapLegend.csv'
        confidencePath        = options.outputPrefix + '-Confidence.tif'
        confidenceLegendPath  = options.outputPrefix + '-ConfidenceLegend.csv'
        mapProjectLeftPath    = options.outputPrefix + '-MapProjLeft.tif'
        mapProjectRightPath   = options.outputPrefix + '-MapProjRight.tif'
        
        compressedInputPath  = options.outputPrefix + '-CompressedInputs.tar.bz2'
        compressedOutputPath = options.outputPrefix + '-CompressedOutputs.tar'
        compressedDebugPath  = options.outputPrefix + '-CompressedDiagnostics.tar.bz2'
        
        # List of .LBL text files to be created
        demLabelPath               = options.outputPrefix + '-DEM.LBL' 
        intersectionErrorLabelPath = options.outputPrefix + '-IntersectionErr.LBL'
        hillshadeLabelPath         = options.outputPrefix + '-Hillshade.LBL'
        colormapLabelPath          = options.outputPrefix + '-Colormap.LBL'
        confidenceLabelPath        = options.outputPrefix + '-Confidence.LBL'
        mapProjectLeftLabelPath    = options.outputPrefix + '-MapProjLeft.LBL'
        mapProjectRightLabelPath   = options.outputPrefix + '-MapProjRight.LBL'

        # List of diagnostic files we want to archive
        #  Currently specified below.  Should be expanded.



        # Call lronac2refinedMosaics.py if outputs do not already exist
        if ( (not os.path.exists(mainMosaicPath) or not os.path.exists(stereoMosaicPath)) 
             and not os.path.exists(demPath) ):
            refineTemp = os.path.join(tempFolder, 'refinement')
            cmdArgs = ['--left',          options.leftPath,
                       '--right',         options.rightPath,
                       '--stereo-left',   options.stereoLeft,
                       '--stereo-right',  options.stereoRight,
                       '--lola',          options.lolaPath,
                       '--workDir',       refineTemp,
                       '--output-folder', tempFolder,
                       '--log-path',      logPath]
            if options.keep:
                cmdArgs.append('--keep')
            print cmdArgs
            lronac2refinedMosaics.main(cmdArgs)
        
            # Copy the pc_align log file to the output folder
            pcAlignLogPath = os.path.join(tempFolder, 'pcAlignLog.txt')
            shutil.copyfile(pcAlignLogPath, options.outputPrefix + '-PcAlignLog.txt')

            # Check that we successfully created the output files
            if (not os.path.exists(mainMosaicPath) or not os.path.exists(stereoMosaicPath)):
                raise Exception('lronac2refinedMosaics failed to produce mosaics!')
 
        # Call makeDemAndCompare.py (This won't do much if outputs already exist)
        cmdArgs = ['--left',     mainMosaicPath, 
                   '--right',    stereoMosaicPath, 
                   '--lola',     options.lolaPath, 
                   '--workDir',  tempFolder, 
                   '--prefix',   options.outputPrefix, 
                   '--log-path', logPath]
        if options.keep:
            cmdArgs.append('--keep')
        if options.asuPath:
            cmdArgs.append('--asu')
            cmdArgs.append(options.asuPath)
        if options.cropAmount:
            cmdArgs.append('--crop')
            cmdArgs.append(str(options.cropAmount))
        if options.nodeFilePath:
            cmdArgs.append('--node-file')
            cmdArgs.append(options.nodeFilePath)
        print cmdArgs
        makeDemAndCompare.main(cmdArgs)
        
        print 'Finished image processing, setting up final output files...'
        logging.info('Started label file generation')

        # Get the data set name
        if (options.asuPath):  # If ASU has a name for this set, use it
            startPos    = options.asuPath.rfind('_') + 1
            stopPos     = options.asuPath.rfind('.') - 1
            dataSetName = options.asuPath[startPos:endPos]
        else: # Generate a data set in format NAC_DTM_MXXXXXX_MXXXXXXX
            dataSetName = IsisTools.makeDataSetName(options.leftPath, options.stereoLeft)

        # Get some data for the label file
        thisFilePath = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'lronac2dem.py')
        versionId   = IrgFileFunctions.getLastGitTag(thisFilePath)

        demDescription        = 'High-resolution NAC digital terrain models in GeoTIFF format. DEM is IEEE floating point TIFF, 32 bits/sample, 1 samples/pixel in single image plane configuration. The NoDATA value is -3.40282266e+038.'
        hillshadeDescription  = 'Shaded-relief derived from NAC digital terrain models. Image is a TIFF image, 8 bits/sample, 1 samples/pixel in single image plane configuration.'
        colormapDescription   = 'Color shaded-relief derived from NAC digital terrain models. Image is 3-channel RGB TIFF, 8 bits/sample, 3 samples/pixel in single image plane configuration.'
        confidenceDescription = 'Discretized visualization of intersection error in NAC digital terrain model. Image is a TIFF image, 8 bits/sample, 1 samples/pixel in single image plane configuration.'
        mapProjectDescription = 'High-resolution NAC orthoprojected image in GeoTIFF format. Image is IEEE floating point TIFF, 32 bits/sample, 1 samples/pixel in single image plane configuration. The NoDATA value is -3.40282266e+038.'

        # Copy the colormap legend file in to the colormap label file!
        colormapLegendText = ('\n/* Colormap legend information: */\n')
        colormapLegendFile = open(colormapLegendPath, 'r')
        for line in colormapLegendFile:
            colormapLegendText = colormapLegendText + '/*  ' + line.strip() + '     */\n'
            
        # Copy the confidence legend file in to the confidence label file!
        confidenceLegendText = ('\n/* Confidence legend information: */\n')
        confidenceLegendFile = open(confidenceLegendPath, 'r')
        for line in confidenceLegendFile:
            confidenceLegendText = confidenceLegendText + '/*  ' + line.strip() + '     */\n'

        # Generate label files for each of the output files
        writeLabelFile(demPath,               demLabelPath,               dataSetName, versionId, demDescription)
        writeLabelFile(hillshadePath,         hillshadeLabelPath,         dataSetName, versionId, hillshadeDescription)
        writeLabelFile(colormapPath,          colormapLabelPath,          dataSetName, versionId, colormapDescription,   colormapLegendText)
        writeLabelFile(confidencePath,        confidenceLabelPath,        dataSetName, versionId, confidenceDescription, confidenceLegendText)
        writeLabelFile(mapProjectLeftPath,    mapProjectLeftLabelPath,    dataSetName, versionId, mapProjectDescription)
        writeLabelFile(mapProjectRightPath,   mapProjectRightLabelPath,   dataSetName, versionId, mapProjectDescription)
        
        logging.info('Done with label files, compressing results.')

        # Compress the input files to save disk space
        if not os.path.exists(compressedInputPath):
            
            listOfInputFiles = [options.leftPath, options.rightPath, options.stereoLeft, options.stereoRight, options.lolaPath]
            if options.asuPath: # Handle optional input
                listOfInputFiles.append(options.asuPath)
            IrgFileFunctions.tarFileList(listOfInputFiles, compressedInputPath, True)
            

        # Tar the output files to save file count (they are already compressed)
        if not os.path.exists(compressedOutputPath):
            
            listOfDeliveryFiles = [demPath,      hillshadePath,      colormapPath,      confidencePath,      mapProjectLeftPath,      mapProjectRightPath,
                                   demLabelPath, hillshadeLabelPath, colormapLabelPath, confidenceLabelPath, mapProjectLeftLabelPath, mapProjectRightLabelPath]

            # The file names as they should appear in the tar file
            baseTarPath  = os.path.join(outputFolder, dataSetName + '_')
            demTarBase   = baseTarPath + 'DEM'
            hillTarBase  = baseTarPath + 'Hillshade'
            colorTarBase = baseTarPath + 'Colormap'
            confTarBase  = baseTarPath + 'Confidence'
            leftTarBase  = baseTarPath + 'MapProjLeft'
            rightTarBase = baseTarPath + 'MapProjRight'

            listOfTarDeliveryFiles = [demTarBase+'.TIF', hillTarBase+'.TIF', colorTarBase+'.TIF', confTarBase+'.TIF', leftTarBase+'.TIF', rightTarBase+'.TIF',
                                      demTarBase+'.LBL', hillTarBase+'.LBL', colorTarBase+'.LBL', confTarBase+'.LBL', leftTarBase+'.LBL', rightTarBase+'.LBL']

            # This function can handle renaming the files as they get put in the tar.
            IrgFileFunctions.tarFileList(listOfDeliveryFiles, compressedOutputPath, False, listOfTarDeliveryFiles)


        # Compress the diagnostic files to save disk space
        if not os.path.exists(compressedDebugPath):
            
            listOfDebugFiles = [options.outputPrefix + '-Log.txt',
                                options.outputPrefix + '-PcAlignLog.txt',
                                options.outputPrefix + '-LOLA_diff_points.csv',
                                options.outputPrefix + '-LOLA_diff_points.kml',]
            IrgFileFunctions.tarFileList(listOfDebugFiles, compressedDebugPath, True)


        if not options.keep:
            print 'Deleting temporary files'
            IrgFileFunctions.removeIfExists(mainMosaicPath)
            IrgFileFunctions.removeIfExists(stereoMosaicPath)
        
        endTime = time.time()

        logging.info('lronac2dem finished in %f seconds', endTime - startTime)
        print "Finished in " + str(endTime - startTime) + " seconds."
        print '#################################################################################'
        return 0

    except Usage, err:
        print err
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
