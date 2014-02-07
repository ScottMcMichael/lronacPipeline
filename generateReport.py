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

import os, glob, optparse, re, shutil, subprocess, sys, string, time, urllib, urllib2, simplekml

import IsisTools
import matplotlib.pyplot as plt
import numpy as np

job_pool = [];

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tools for compiling statistics for lronacPipeline
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
    #print "Waiting for jobs to finish";
    while len(job_pool) > 0:
        job_pool[0].wait();
        job_pool.pop(0);


#--------------------------------------------------------------------------------

# Reads in the results of a DEM alignment operation
def readPcAlignResults(outputFolder, prefix):

    # Output transform is in this format (excepting whitespace):
    #  0.9999968710678928    0.002478390929168969 0.0003397540659725261 -1024.133917521918
    # -0.002477949646550771  0.9999960937499494  -0.001293155207595954   4739.357603680983
    # -0.0003429576829447021 0.001292309267933415 0.9999991061579924      124.7476123821689
    #  0                     0                    0                         1

    # File locations (PREFIX = align_ASU and align_LOLA):
    # outputFolder/PREFIX-transform.txt
    # outputFolder/PREFIX-beg_errors.txt
    # outputFolder/PREFIX-end_errors.txt

    # Get the file paths
    transformPath    = outputFolder + '/' + prefix + '-transform.txt'
    initialErrorPath = outputFolder + '/' + prefix + '-beg_errors.txt'
    finalErrorPath   = outputFolder + '/' + prefix + '-end_errors.txt'
    
#TODO: Make sure to pull in ASP update which inverts output matrix if needed!
    # Read in the transform matrix
    transformMatrix = []
    transformFile = open(transformPath, 'r')
    for line in transformFile:
        transformMatrix.append(float, line.split(' '))

    # For now just return the matrix
    return transformMatrix


#def backupData(localFolder):

#    #TODO: Don't duplicate these!
#    supercomputerSourceFolder = '/u/smcmich1/data/lronacPipeline'

#    demFolderList = [ \
#                      'ARISTARCHU2',  'FEOKTISTOV',   'HORTENSIUS1',  'KINGCRATER2',   'LICHTENBER7',  'MRECRISIUM2',  'ORIENTALE2',    'SEISMCLAND',\
#                      'ARISTARCHU3',  'FRESH1',       'HORTENSIUS2',  'KINGCRATER3',   'LICHTENBER8',  'MRECRISIUM3',  'ORIENTALE3',    'SLIPHER1',\
#                      'ARISTARCHU4',  'FRESH3',       'HORTENSIUS3',  'KINGCRATER4',   'LICHTENBER9',  'MRINGENII1',   'PLANCKFLOOR',   'SOSIRILLE',\
#                      'ARISTARCHU5',  'FRESHMELT1',   'HORTENSIUS5',  'KUGLERRIDGE',   'LINNECRATER',  'MRINGENII2',   'PYTHAGORAS1',   'SPARIM1',\
#                      'BHABHAPLAIN',  'FRESHMELT2',   'HORTENSIUS6',  'LICHTENBER1',   'LUNA16',       'MRINGENII3',   'RANGER',        'SPARIM2',\
#                      'CMPTNBELK2',   'FRSHCRATER4',  'HORTENSIUS7',  'LICHTENBER10',  'LUNA20',       'MRINGENII4',   'REINER1',       'SPARIM3',\
#                      'CMPTNBELK3',   'FRSHCRATER5',  'IMBRIUM',      'LICHTENBER11',  'LUNA24',       'NEARLENTS1',   'REINER2',       'SPARIM4',\
#                      'COMPTONBELK',  'FRSHCRATER8',  'IMBRIUM2',     'LICHTENBER12',  'MARIUS1',      'NRTHCRTRI1',   'REINER3',       'SULPICIUS1',\
#                      'EIMMARTA',     'GRUITHUISE1',  'IMPACTMELT1',  'LICHTENBER13',  'MARIUS2',      'NRTHCRTRI2',   'REINER4',       'SULPICIUS2',\
#                      'ENDYMION',     'GRUITHUISE2',  'IMPACTMELT2',  'LICHTENBER2',   'MARIUS3',      'NRTHCRTRII',   'RUMKERDOME1',   'SULPICIUS3',\
#                      'ERATOSTHNS1',  'GRUITHUISE3',  'IMPACTMELT3',  'LICHTENBER3',   'MARIUS4',      'NRTHCRTRIII',  'RUMKERDOME2',   'VIRTANEN1',\
#                      'ESALL_CR1',    'GRUITHUISE4',  'INACALDERA1',  'LICHTENBER4',   'MOOREF1',      'OBLIVIONIS1',  'RUMKERDOME3',   'VIRTANEN2',\
#                      'ESALL_MP1',    'GRUITHUISE5',  'INACALDERA3',  'LICHTENBER5',   'MOOREF2',      'OPENHIMERF',   'RUMKERDOME4',   'VIRTANEN3',\
#                      'ESALL_SR12',   'HIGHESTPOIN',  'ISISOSIRIS',   'LICHTENBER6',   'MRECRISIUM1',  'ORIENTALE1',   'RUMKERDOMES5',  'VITELLO']

#    try:
  
#        if not os.path.exists(localFolder):
#            os.makedirs(localFolder)

#        # Operate on each folder
#        for f in demFolderList:
	
#            # Get folders
#            folderPath              = os.path.join(localFolder,               f)
#            supercomputerFolderPath = os.path.join(supercomputerSourceFolder, f)

#            if not os.path.exists(folderPath):
#                os.makedirs(folderPath)

#            # Build file request
#            inputFile = '*.lronaccal.lronacecho.noproj.mosaic.norm.cub'
#            supercomputerFilePath = '"' + os.path.join(supercomputerFolderPath, inputFile) + '"'

#            # Check if we have both input files
#            cubCounter = len(glob.glob1(folderPath,"*.cub"))
#            if cubCounter < 2: # If we don't, request file copy
#                IsisTools.grabSupercomputerFile(supercomputerFilePath, folderPath)            
    
#        wait_on_all_jobs() # Wait on all file grab requests
        
#    except Exception,e: # Catch any errors, the program will move on to the next folder
#        print "Caught: ", e


# Grabs the output files from the supercomputer
def grabResultFiles(localFolder):

    supercomputerSourceFolder = '/u/smcmich1/data/lronacPipeline'

    demFolderList = [ \
                      'ARISTARCHU2',  'FEOKTISTOV',   'HORTENSIUS1',  'KINGCRATER2',   'LICHTENBER7',  'MRECRISIUM2',  'ORIENTALE2',    'SEISMCLAND',\
                      'ARISTARCHU3',  'FRESH1',       'HORTENSIUS2',  'KINGCRATER3',   'LICHTENBER8',  'MRECRISIUM3',  'ORIENTALE3',    'SLIPHER1',\
                      'ARISTARCHU4',  'FRESH3',       'HORTENSIUS3',  'KINGCRATER4',   'LICHTENBER9',  'MRINGENII1',   'PLANCKFLOOR',   'SOSIRILLE',\
                      'ARISTARCHU5',  'FRESHMELT1',   'HORTENSIUS5',  'KUGLERRIDGE',   'LINNECRATER',  'MRINGENII2',   'PYTHAGORAS1',   'SPARIM1',\
                      'BHABHAPLAIN',  'FRESHMELT2',   'HORTENSIUS6',  'LICHTENBER1',   'LUNA16',       'MRINGENII3',   'RANGER',        'SPARIM2',\
                      'CMPTNBELK2',   'FRSHCRATER4',  'HORTENSIUS7',  'LICHTENBER10',  'LUNA20',       'MRINGENII4',   'REINER1',       'SPARIM3',\
                      'CMPTNBELK3',   'FRSHCRATER5',  'IMBRIUM',      'LICHTENBER11',  'LUNA24',       'NEARLENTS1',   'REINER2',       'SPARIM4',\
                      'COMPTONBELK',  'FRSHCRATER8',  'IMBRIUM2',     'LICHTENBER12',  'MARIUS1',      'NRTHCRTRI1',   'REINER3',       'SULPICIUS1',\
                      'EIMMARTA',     'GRUITHUISE1',  'IMPACTMELT1',  'LICHTENBER13',  'MARIUS2',      'NRTHCRTRI2',   'REINER4',       'SULPICIUS2',\
                      'ENDYMION',     'GRUITHUISE2',  'IMPACTMELT2',  'LICHTENBER2',   'MARIUS3',      'NRTHCRTRII',   'RUMKERDOME1',   'SULPICIUS3',\
                      'ERATOSTHNS1',  'GRUITHUISE3',  'IMPACTMELT3',  'LICHTENBER3',   'MARIUS4',      'NRTHCRTRIII',  'RUMKERDOME2',   'VIRTANEN1',\
                      'ESALL_CR1',    'GRUITHUISE4',  'INACALDERA1',  'LICHTENBER4',   'MOOREF1',      'OBLIVIONIS1',  'RUMKERDOME3',   'VIRTANEN2',\
                      'ESALL_MP1',    'GRUITHUISE5',  'INACALDERA3',  'LICHTENBER5',   'MOOREF2',      'OPENHIMERF',   'RUMKERDOME4',   'VIRTANEN3',\
                      'ESALL_SR12',   'HIGHESTPOIN',  'ISISOSIRIS',   'LICHTENBER6',   'MRECRISIUM1',  'ORIENTALE1',   'RUMKERDOMES5',  'VITELLO']


    # List of all the temporary files we want copied
    fileList = [  'ASU_diff_stats.txt', \
                  'LOLA_ASU_diff_stats.txt' ,\
                  'LOLA_diff_stats.txt' ,\
                  'align_ASU_PC-transform.txt' ,\
                  'align_LOLA_PC-transform.txt' ,\
                  'align_ASU-DEM.tif' ,\
                  'align_LOLA-DEM.tif' ,\
                  'stereo-DEM.tif' ,\
                  'geodiff_ASU-diff.tif']
    fileList2 = []

    try:
  
        if not os.path.exists(localFolder):
            os.makedirs(localFolder)

        # Operate on each folder
        for f in demFolderList:
	
            # Get folders
            folderPath              = os.path.join(localFolder ,              f)
            supercomputerFolderPath = os.path.join(supercomputerSourceFolder, f)

            if not os.path.exists(folderPath):
                os.makedirs(folderPath)

            # Add the ASU input DEM to the file copy list
            fileList.append('NAC_DTM_' + f + '.TIF')

            # For each folder, copy each result
            for c in fileList:
		            # Get file paths
                outputFilePath        = os.path.join(folderPath,              c)
                supercomputerFilePath = os.path.join(supercomputerFolderPath, c)

                if not os.path.exists(outputFilePath): # Request file copy
                    IsisTools.grabSupercomputerFile(supercomputerFilePath, outputFilePath)

            fileList.pop() # Remove the last ASU input DEM from the list
		            
        
        wait_on_all_jobs() # Wait on all file grab requests
        
    except Exception,e: # Catch any errors, the program will move on to the next folder
        print "Caught: ", e


def accumulateStatistics(prefix, filePath, meanList, stdDevList, meanHistogram, dataStorage):

    try:
        meanValue, stdDev, percentile, histogram = IsisTools.readLolaCompareFile(filePath)

    except Exception,e: # Catch any errors, the program will move on to the next folder
        #print "Caught: ", e
        #print "Unable to process data in file " + filePath
        return False

#    print 'Read data from ' + filePath

    meanList.append(meanValue)
    stdDevList.append(stdDev)
  
    # Record the raw percentile data
    dataStorage[prefix] = percentile
  
    if not meanHistogram:
        meanHistogram.extend(histogram)
    else:
        for o, n in zip(meanHistogram, histogram):
            o = o + n

    return True # We successfully added data

def readShiftAmounts(path, prefix, shiftDict):

    # Make sure the file exists
    if not os.path.exists(path):
        return False

    # Read the three numbers from the file
    shiftFile = open(path, 'r')
    line      = shiftFile.readline()
    parts     = line.split(' ')
    transform = (float(parts[0]), float(parts[1]), float(parts[2]))
    shiftDict[prefix] = transform
    
    return True

# Obtains the bounding box (as specified by ASU) for each DTM
def parseBoundingBoxes(filePath, storage):
        
    logFile = open(filePath, 'r')
    currentDtm = ''
    for line in logFile:
        if (line.find('ASU DTM') >= 0):         # Started a new DTM
            eqPos       = line.find('=')
            asuUrl	    = line[eqPos+2:]
            asuFileName = os.path.basename(asuUrl)
            noExt       = os.path.splitext(asuFileName)[0]
            currentDtm  = noExt[8:]
            print currentDtm
        elif (line.find('Min lat') >= 0):
            eqPos          = line.find('=')
            numberText     = line[eqPos+2:]
            value          = float(numberText)
            index          = currentDtm + '_minLat'
            storage[index] = value
        elif (line.find('Min lon') >= 0):
            eqPos          = line.find('=')
            numberText     = line[eqPos+2:]
            value          = float(numberText)
            index          = currentDtm + '_minLon'
            storage[index] = value
        elif (line.find('Max lat') >= 0):
            eqPos          = line.find('=')
            numberText     = line[eqPos+2:]
            value          = float(numberText)
            index          = currentDtm + '_maxLat'
            storage[index] = value
        elif (line.find('Max lon') >= 0):
            eqPos          = line.find('=')
            numberText     = line[eqPos+2:]
            value          = float(numberText)
            index          = currentDtm + '_maxLon'
            storage[index] = value

def generateKml(dataFolder, usedFolderList, bbStorage, asuMeanList):
    
    # Initialize kml document
    kml = simplekml.Kml()
    kml.document.name = 'test'
    kml.hint = 'target=moon'
    
    minError = min(asuMeanList)
    maxError = max(asuMeanList)
    errorRange = maxError - minError
    
    # Try to set up bounding box for each used folder
    i = 0
    for f in usedFolderList:
        minLat = bbStorage[f+'_minLat']
        maxLat = bbStorage[f+'_maxLat']
        minLon = bbStorage[f+'_minLon']
        maxLon = bbStorage[f+'_maxLon']
        
        poly = kml.newpolygon(name=f, outerboundaryis=[(maxLat,minLon), (maxLat,maxLon), (minLat, maxLon), (minLat,minLon), (maxLat,minLon)])
#        poly.style.polystyle.color   = simplekml.Color.red # Why is polygon not working?
#        poly.style.polystyle.fill    = 1
#        poly.style.polystyle.outline = 1
        poly.style.linestyle.width   = 5
        
        # Generate a color based on the error value:  white (low error) <--> red (high error)
        colorVal = 255 - 255*(asuMeanList[i] - minError)/errorRange
        poly.style.linestyle.color   = simplekml.Color.rgb(255,colorVal,colorVal,255)
        
        i = i + 1

      
    
    # Save kml document
    outputPath = dataFolder + '/bbList.kml'
    kml.save(outputPath)
    


# Generates a set of plots describing the results
def generatePlots(dataFolder):

    # Search through all the output files and pull out results

    asuMeanList       = []
    asuStdDevList     = []
    asuMeanPercentile = []
    asuMeanHistogram  = []

    lolaMeanList       = []
    lolaStdDevList     = []
    lolaMeanPercentile = []
    lolaMeanHistogram  = []

    lolaAsuMeanList       = []
    lolaAsuStdDevList     = []
    lolaAsuMeanPercentile = []
    lolaAsuMeanHistogram  = []
    
    dataStorage      = dict()
    shiftAmountsAsu  = dict()
    shiftAmountsLola = dict()

    # Parse bounding boxes
    logPath = dataFolder + '/logFile.txt'
    bbStorage = dict()
    parseBoundingBoxes(logPath, bbStorage)

    usedFolderList = []
    for f in os.listdir(dataFolder):

        folderPath = os.path.join(dataFolder, f)
	
        # Get the file paths for this folder
        asuDiffPath         = folderPath + '/ASU_diff_stats.txt'
        lolaDiffPath        = folderPath + '/LOLA_diff_stats.txt'
        lolaAsuDiffPath     = folderPath + '/LOLA_ASU_diff_stats.txt'
        asuTransformPath    = folderPath + '/align_ASU_PC-transform.txt'
        lolaTransformPath   = folderPath + '/align_LOLA_PC-transform.txt'
        stereoDemPath       = folderPath + '/stereo-DEM.tif'
        
        # Create output file paths for this folder
        asuTransformLogPath  = folderPath + '/asuLocalShift.txt'
        lolaTransformLogPath = folderPath + '/lolaLocalShift.txt'
               
        # Compute the shift amount in local coordinates for ASU and LOLA
        if (not os.path.exists(asuTransformLogPath)  and os.path.exists(asuTransformPath) ):
            cmd = '~/repot/StereoPipeline/src/asp/Tools/transformConvert --output-log ' + asuTransformLogPath  + ' ' + asuTransformPath  + ' ' + stereoDemPath
            add_job(cmd, 2)
        if (not os.path.exists(lolaTransformLogPath) and os.path.exists(lolaTransformPath) ):
            cmd = '~/repot/StereoPipeline/src/asp/Tools/transformConvert --output-log ' + lolaTransformLogPath + ' ' + lolaTransformPath + ' ' + stereoDemPath
            add_job(cmd, 2)
        wait_on_all_jobs()
        
        readAsuTransform  = readShiftAmounts(asuTransformLogPath,  f, shiftAmountsAsu)
        readLolaTransform = readShiftAmounts(lolaTransformLogPath, f, shiftAmountsLola)
         
               
        # Accumulate the statistics
        readAsu     = accumulateStatistics(f+'_asu',  asuDiffPath,     asuMeanList,     asuStdDevList,     asuMeanHistogram,     dataStorage)
        readLola    = accumulateStatistics(f+'_lola', lolaDiffPath,    lolaMeanList,    lolaStdDevList,    lolaMeanHistogram,    dataStorage )
        readLolaAsu = accumulateStatistics(f+'_comp', lolaAsuDiffPath, lolaAsuMeanList, lolaAsuStdDevList, lolaAsuMeanHistogram, dataStorage )

        if readAsu and readLola and readLolaAsu and readAsuTransform and readLolaTransform:
            usedFolderList.append(f) # Keep track of the folders we read data from
        else:
            # Remove list entries from partial successes
            if readAsu:
                asuMeanList.pop()
            if readLola:
                lolaMeanList.pop()
            if readLolaAsu:
                lolaAsuMeanList.pop()
                
                
            if readAsu or readLola or readLolaAsu:
                print 'Partial success for folder ' + f

    # Generate some google earth KML output
    generateKml(dataFolder, usedFolderList, bbStorage, asuMeanList)
    
    return True


    # Scatter plot of the local coordinate x/y offsets
    sScaling = 1.0
    xAsu  = []
    yAsu  = []
    sAsu  = []
    xLola = []
    yLola = []
    sLola = []
    meanAsuX  = 0
    meanAsuY  = 0
    meanAsuS  = 0
    meanLolaX = 0
    meanLolaY = 0
    meanLolaS = 0
    for f in usedFolderList:
        xAsu.append(shiftAmountsAsu[f][0])
        yAsu.append(shiftAmountsAsu[f][1])       
        meanAsuX = meanAsuX + shiftAmountsAsu[f][0]
        meanAsuY = meanAsuY + shiftAmountsAsu[f][1]
        
        s = abs(shiftAmountsAsu[f][2]) * sScaling
        sAsu.append(s)
        meanAsuS = meanAsuS + s

        xLola.append(shiftAmountsLola[f][0])
        yLola.append(shiftAmountsLola[f][1])
        sLola.append(abs(shiftAmountsLola[f][2]))
        meanLolaX = meanLolaX + shiftAmountsLola[f][0]
        meanLolaY = meanLolaY + shiftAmountsLola[f][1]        

        s = abs(shiftAmountsLola[f][2]) * sScaling
        sLola.append(s)
        meanLolaS = meanLolaS + s

    meanAsuX  = meanAsuX  / len(usedFolderList)
    meanAsuY  = meanAsuY  / len(usedFolderList)
    meanAsuS  = meanAsuS  / len(usedFolderList)                
    meanLolaX = meanLolaX / len(usedFolderList)
    meanLolaY = meanLolaY / len(usedFolderList)
    meanLolaS = meanLolaS / len(usedFolderList)

    plt.scatter(xAsu,      yAsu,      sAsu,      c='r', marker='o', label='ASU offsets')
    plt.scatter(xLola,     yLola,     sLola,     c='b', marker='o', label='LOLA offsets')
    plt.scatter(meanAsuX,  meanAsuY,  meanAsuS,  c='g', marker=(5,1), label='mean ASU offset')
    plt.scatter(meanLolaX, meanLolaY, meanLolaS, c='y', marker=(5,1), label='mean LOLA offset')
    plt.grid(color='gray', linestyle='dashed')
    plt.xlabel('X shift in meters')
    plt.ylabel('Y shift in meters')
    plt.title('Alignment offsets in NED frame')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/localShift.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()


    # We now have a list of means/stds and sums of percentiles and histograms
    # - Get means for all of them
    numElements = float(len(asuMeanList))
    asuMeanOfMeans       = sum(asuMeanList      ) / numElements
    asuMeanOfStdDevs     = sum(asuStdDevList    ) / numElements
    lolaMeanOfMeans      = sum(lolaMeanList     ) / numElements
    lolaMeanOfStdDevs    = sum(lolaStdDevList   ) / numElements
    lolaAsuMeanOfMeans   = sum(lolaAsuMeanList  ) / numElements
    lolaAsuMeanOfStdDevs = sum(lolaAsuStdDevList) / numElements

    for i in asuMeanHistogram:
      i = i / numElements

    for i in lolaMeanHistogram:
      i = i / numElements

    for i in lolaAsuMeanHistogram:
      i = i / numElements

    # Now generate plots

    # Three more plots with all data points

    numEls = 20 # Number of elements in percentile plots (leaving off last one containing error values)
    yMin   = 0
    yMax   = 10
    xAxis  = range(0, numEls)
    xAxis  = [x * 0.05 for x in xAxis]
    asuMeanPercentile     = [0] * numEls
    lolaMeanPercentile    = [0] * numEls
    lolaAsuMeanPercentile = [0] * numEls

    # us vs ASU
    for f in usedFolderList:
        # Plot data points for this folder
        dataStorage[f+'_asu'].pop() # Strip of error values in last bin
        plt.plot(xAxis, dataStorage[f+'_asu'], 'o', label=f)
        
        # Accumulate mean value
        for i in range(0,numEls):
            asuMeanPercentile[i] = asuMeanPercentile[i] + dataStorage[f+'_asu'][i]
            
    # Finish and plot mean value
    for i in range(0,numEls): 
        asuMeanPercentile[i] = asuMeanPercentile[i] / len(usedFolderList)
    plt.plot(xAxis, asuMeanPercentile, '-', label='mean percentiles')
    
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('Us vs ASU pixel percentages for all samples')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesASU.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()



    # us vs LOLA
    for f in usedFolderList:
        # Plot data points for this folder
        dataStorage[f+'_lola'].pop() # Strip of error values in last bin
        plt.plot(xAxis, dataStorage[f+'_lola'], 'o', label=f)
        
        # Accumulate mean value
        for i in range(0,numEls):
            lolaMeanPercentile[i] = lolaMeanPercentile[i] + dataStorage[f+'_lola'][i]
            
    # Finish and plot mean value
    for i in range(0,numEls): 
        lolaMeanPercentile[i] = lolaMeanPercentile[i] / len(usedFolderList)
    plt.plot(xAxis, lolaMeanPercentile, '-', label='mean percentiles')
    
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('Us vs LOLA pixel percentages for all samples')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesLOLA.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()
    
    
    
    # LOLA vs ASU
    for f in usedFolderList:
        # Plot data points for this folder
        dataStorage[f+'_comp'].pop() # Strip of error values in last bin
        plt.plot(xAxis, dataStorage[f+'_comp'], 'o', label=f)
        
        # Accumulate mean value
        for i in range(0,numEls):
            lolaAsuMeanPercentile[i] = lolaAsuMeanPercentile[i] + dataStorage[f+'_comp'][i]
            
    # Finish and plot mean value
    for i in range(0,numEls): 
        lolaAsuMeanPercentile[i] = lolaAsuMeanPercentile[i] / len(usedFolderList)
    plt.plot(xAxis, lolaAsuMeanPercentile, '-', label='mean percentiles')
    
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('LOLA vs ASU pixel percentages for all samples')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesLolaASU.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()    
    
    
    # Plot the three mean percentiles on one chart
    plt.plot(xAxis, asuMeanPercentile    , label='us vs ASU')
    plt.plot(xAxis, lolaMeanPercentile   , label='us vs LOLA')
    plt.plot(xAxis, lolaAsuMeanPercentile, label='ASU vs LOLA')
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('Cumulative pixel percentiles averaged across DEMs')
    plt.legend(loc='upper left')
    plt.savefig(dataFolder + '/meanPercentiles.png')
    plt.clf()


    # Now plot the mean error for each folder
    yMin = 0
    yMax = 6
    xAxis = np.arange(len(asuMeanList))
    barwidth = 0.2
    plt.bar(xAxis,            asuMeanList,     barwidth, color='r', label='us vs ASU')
    plt.bar(xAxis+1*barwidth, lolaMeanList,    barwidth, color='g', label='us vs LOLA')
    plt.bar(xAxis+2*barwidth, lolaAsuMeanList, barwidth, color='b', label='ASU vs LOLA')
    plt.xticks(xAxis+.5, usedFolderList, size='small', rotation='vertical')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Mean difference between DEMs')
    plt.savefig(dataFolder + '/means.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()



#TODO: Multi-plot of this in addition to mean

#    --> Not using this plot anymore
#    # Mean histograms
#    yMin = 0
#    yMax = 20
#    xAxis = range(0,len(asuMeanHistogram)) #TODO: Label with number of std's
#    plt.plot(xAxis, asuMeanHistogram,     label='us vs ASU')
#    plt.plot(xAxis, lolaMeanHistogram,    label='us vs LOLA')
#    plt.plot(xAxis, lolaAsuMeanHistogram, label='ASU vs LOLA')
#    plt.ylim(yMin, yMax)
#    plt.ylabel('Percent pixels in bin')
#    plt.title('Difference histograms averaged across DEMs')
#    plt.savefig(dataFolder + '/meanHistograms.png')
#    plt.clf()

#    # Simple plot example
#    x = np.arange(0, 10, 0.2)
#    y = np.sin(x)
#    plt.plot(x, y)
#    plt.show()


#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    try:
        try:
            usage = "usage: lronacPipeline.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.set_defaults(delete =True)
            parser.set_defaults(threads=4)
            parser.set_defaults(fakePvl=True)
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

#            if not args: parser.error("need .IMG files")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        dataFolder = '/byss/moon/lronacPipeline'
        grabResultFiles(dataFolder)
#        backupData(dataFolder)
#        generatePlots(dataFolder)



        print "Finished"
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

	# To more easily debug this program, comment out this catch block.
    # except Exception, err:
    #     sys.stderr.write( str(err) + '\n' )
    #     return 1


if __name__ == "__main__":
    sys.exit(main())
