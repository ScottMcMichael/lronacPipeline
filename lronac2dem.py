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

import IsisTools, simplekml, matplotlib

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

def writeLabelFile(outputPath, dataSetName, versionId, imageSize, boundingBox):
    """Write out a .LBL file formatted for the PDS"""
    
    labelPath = os.path.join(outputFolder, 'labelFile.txt')
    labelFile = open(labelPath, 'w')
    
    labelFile.write('PDS_VERSION_ID            = PDS3\n') #TODO
    
    labelFile.write('/* The source image data definition. */\n')
    labelFile.write('^IMAGE        = ' + os.path.basename(outputPath) +'\n')
    
    labelFile.write('/* Identification Information  */\n')
    labelFile.write('DATA_SET_ID               = "LRO-L-LROC-5-RDR-V1.0"\n') #TODO
    labelFile.write('DATA_SET_NAME             = "LRO MOON LROC 5 RDR V1.0"\n')
    labelFile.write("VOLUME_ID                 = 'LROLRC_2001'\n")
    labelFile.write('PRODUCER_INSTITUTION_NAME = "NASA Ames"\n')
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
    labelFile.write("""RATIONALE_DESC            = "Created at the request of NASA's Exploration\n""") #TODO!!!!!!!!!!!!!
    labelFile.write('                            Systems Mission Directorate to support future\n')
    labelFile.write('                            human exploration"\n')
    labelFile.write('SOFTWARE_NAME             = "ISIS 3.2.1 with SER enhancements | SOCET SET v5.5\n') #TODO!!!!!!!!!!!!!!
    labelFile.write('DESCRIPTION               = "High-resolution NAC ditital terrain models in GeoTIFF format. DEM is IEEE floating point TIFF, 32 bits/sample, 1 samples/pixel in single image plane configuration. The NoDATA value is -3.40282266e+038."\n')
    labelFile.write('\n')
    labelFile.write('/* Time Parameters */\n')
    labelFile.write('START_TIME                   = "N/A"\n')
    labelFile.write('STOP_TIME                    = "N/A"\n')
    labelFile.write('SPACECRAFT_CLOCK_START_COUNT = "N/A"\n')
    labelFile.write('SPACECRAFT_CLOCK_STOP_COUNT  = "N/A"\n')
    labelFile.write('PRODUCT_CREATION_TIME        = ' + time.strftime("%Y-%m-%dT%H:%M:%S") + '\n')
    labelFile.write('\n')
    labelFile.write('/* NOTE:                                                                   */\n')
    labelFile.write('/* This raster image is composed of a set of pixels that represent finite  */\n') #TODO!!!!!!!!!!!!!!!!!!!!!!!!
    labelFile.write('/* areas, and not discrete points.  The center of the upper left pixel is  */\n')
    labelFile.write('/* defined as line and sample (1.0,1.0). The                               */\n')
    labelFile.write('/* [LINE,SAMPLE]_PROJECTION_OFFSET elements are the pixel offset from line */\n')
    labelFile.write('/* and sample (1.0,1.0) to the map projection origin (defined by the       */\n')
    labelFile.write('/* CENTER_LATITUDE and CENTER_LONGITUDE elements).  These offset values    */\n')
    labelFile.write('/* are positive when the map projection origin is to the right or below    */\n')
    labelFile.write('/* the center of the upper left pixel. This definition was adopted in      */\n')
    labelFile.write('/* November 2011 by the LROC team.                                         */\n')
    labelFile.write('\n')
    labelFile.write('OBJECT = IMAGE_MAP_PROJECTION\n')
    labelFile.write('    ^DATA_SET_MAP_PROJECTION     = "DSMAP.CAT"\n') #TODO
    labelFile.write('    MAP_PROJECTION_TYPE          = EQUIRECTANGULAR\n')
    labelFile.write('    PROJECTION_LATITUDE_TYPE     = PLANETOCENTRIC\n') #TODO
    labelFile.write('    A_AXIS_RADIUS                = 1737.4 <KM>\n')
    labelFile.write('    B_AXIS_RADIUS                = 1737.4 <KM>\n')
    labelFile.write('    C_AXIS_RADIUS                = 1737.4 <KM>\n')
    labelFile.write('    COORDINATE_SYSTEM_NAME       = PLANETOCENTRIC\n') #TODO
    labelFile.write('    POSITIVE_LONGITUDE_DIRECTION = EAST\n') #TODO
    labelFile.write('    KEYWORD_LATITUDE_TYPE        = PLANETOCENTRIC\n') #TODO
    labelFile.write('    /* NOTE:  CENTER_LATITUDE and CENTER_LONGITUDE describe the location   */\n') #TODO
    labelFile.write('    /* of the center of projection, which is not necessarily equal to the  */\n')
    labelFile.write('    /* location of the center point of the image.                          */\n')
    labelFile.write('    CENTER_LATITUDE              = 25.0 <DEG>\n') #TODO
    labelFile.write('    CENTER_LONGITUDE             = 180.0 <DEG>\n')
    labelFile.write('    LINE_FIRST_PIXEL             = 1\n')
    labelFile.write('    LINE_LAST_PIXEL              = ' + str(imageSize[1] + 1) + '\n')
    labelFile.write('    SAMPLE_FIRST_PIXEL           = 1\n')
    labelFile.write('    SAMPLE_LAST_PIXEL            = ' + str(imageSize[0] + 1) + '\n')
    labelFile.write('    MAP_PROJECTION_ROTATION      = 0.0 <DEG>\n') #TODO
    labelFile.write('    MAP_RESOLUTION               = 15161.675 <PIX/DEG>\n') #TODO
    labelFile.write('    MAP_SCALE                    = 2.00 <METERS/PIXEL>\n') #TODO
    labelFile.write('    MAXIMUM_LATITUDE             = ' + boundingBox[0] + ' <DEG>\n') #TODO!!!!!!!!!!!!!!!!!1
    labelFile.write('    MINIMUM_LATITUDE             = ' + boundingBox[0] + ' <DEG>\n')
    labelFile.write('    EASTERNMOST_LONGITUDE        = ' + boundingBox[0] + ' <DEG>\n')
    labelFile.write('    WESTERNMOST_LONGITUDE        = ' + boundingBox[0] + ' <DEG>\n')
    labelFile.write('    LINE_PROJECTION_OFFSET       = 379713.5 <PIXEL>\n') #TOOO
    labelFile.write('    SAMPLE_PROJECTION_OFFSET     = -1802158.5 <PIXEL>\n') #TODO
    labelFile.write('END_OBJECT = IMAGE_MAP_PROJECTION\n')
    labelFile.write('\n')
    labelFile.write('END\n')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    labelFile.close()
    return True

# Makes sure all needed functions are found in the PATH
def functionStartupCheck():

    # These calls will raise an exception if the tool is not found
    IsisTools.checkIfToolExists('lronac2refinedMosaics.py')
    IsisTools.checkIfToolExists('makeDemAndCompare.py')
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
            
            parser.add_option_group(inputGroup)

            # The default working directory path is kind of ugly...
            parser.add_option("--workDir", dest="workDir",  help="Folder to store temporary files in")

            parser.add_option("--output-folder",  dest="outputFolder",   help="Output folder.")

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
            if not options.outputFolder: 
                parser.error("Need output folder")

        except optparse.OptionError, msg:
            raise Usage(msg)

        startTime = time.time()

        # Set up the output folders
        outputPrefix  = options.outputFolder + '/output'
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
        logPath = outputPrefix + '-Log.txt'
        logging.basicConfig(filename=logPath,level=logging.INFO)


        # These are the output mosaic paths
        filename         = os.path.splitext(options.leftPath)[0] + '.correctedMosaic.cub'
        mainMosaicPath   = os.path.join(tempFolder, os.path.basename(filename))
        filename         = os.path.splitext(options.stereoLeft)[0] + '.correctedMosaic.cub'
        stereoMosaicPath = os.path.join(tempFolder, os.path.basename(filename))

        # Call lronac2refinedMosaics.py
        if (not os.path.exists(mainMosaicPath) or not os.path.exists(stereoMosaicPath)):    
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
        shutil.copyfile(pcAlignLogPath, os.path.join(outputFolder, 'pcAlignLog.txt'))

        # Check that we successfully created the output files
        if (not os.path.exists(mainMosaicPath) or not os.path.exists(stereoMosaicPath)):
            raise Exception('lronac2refinedMosaics failed to produce mosaics!')

        # Call makeDemAndCompare.py
        outputDemPath = outputPrefix + '-DEM.tif'
        cmdArgs = ['--left',     mainMosaicPath, 
                   '--right',    stereoMosaicPath, 
                   '--lola',     options.lolaPath, 
                   '--asu',      options.asuPath, 
                   '--workDir',  tempFolder, 
                   '--prefix',   outputPrefix, 
                   '--log-path', logPath]
        if options.keep:
            cmdArgs.append('--keep')
        if options.cropAmount:
            cmdArgs.append('--crop')
            cmdArgs.append(str(options.cropAmount))
        print cmdArgs
        makeDemAndCompare.main(cmdArgs)

        # Get the data set name
        if (options.asuPath):  # If ASU has a name for this set, use it
            startPos    = options.asuPath.rfind('_') + 1
            stopPos     = options.asuPath.rfind('.') - 1
            dataSetName = options.asuPath[startPos:endPos]
        else: # Generate a data set in format NAC_DTM_MXXXXXX_MXXXXXXX
            dataSetName = IsisTools.getDataSetName(options.leftPath, options.leftStereoPath)

        # Get some data for the label file
        thisFilePath = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'lronac2dem.py')
        versionId   = IsisTool.getLastGitTag(thisFilePath)
        imageSize   = IsisTools.getCubeSize(outputDemPath)
        boundingBox = IsisTools.getCubeBoundingBox(outputDemPath, tempFolder)

        # Generate label file
        labelPath = os.path.join(outputFolder, 'labelFile.txt')
        writeLabelFile(outputPath, dataSetName, versionId, imageSize, boundingBox)
        

        # Compress the input files to save disk space
        compressedPath = os.path.join(options.outputFolder, 'compressedInputs.tar.bz2')
        if not os.path.exists(compressedPath):
            cmd = ('tar -jcvf ' + compressedPath + ' ' + options.leftPath   + ' ' + options.rightPath
                                                 + ' ' + options.stereoLeft + ' ' + options.stereoRight
                                                 + ' ' + options.lolaPath   + ' ' + options.asuPath)
            print cmd
            os.system(cmd)


        if not options.keep:
            print 'Deleting temporary files'
            IsisTools.removeIfExists(mainMosaicPath)
            IsisTools.removeIfExists(stereoMosaicPath)
        
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
