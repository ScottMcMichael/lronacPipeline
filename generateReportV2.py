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

import matplotlib.pyplot as plt
import numpy as np

import IsisTools

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tools for compiling statistics for lronacPipeline
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------


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
#                grabFile(supercomputerFilePath, folderPath)            
    
#        wait_on_all_jobs() # Wait on all file grab requests
        
#    except Exception,e: # Catch any errors, the program will move on to the next folder
#        print "Caught: ", e


# Grabs the output files from the supercomputer
def grabResultFiles(localFolder, force):

    #supercomputerSourceFolder = '/u/smcmich1/data/lronacPipeline'
    supercomputerSourceFolder = '/nobackupnfs2/oalexan1/scott'
    #supercomputerSourceFolder = '/u/smcmich1/data/lronacProduction'

    # List of processed production data
    #demFolderList = ['NAC_DTM_M151318807_M181974094']
    
    
    # Entire ASU folder list
    demFolderList = [ \
                        'ARISTARCHU2',  'FECUNPIT',	     'HORTENSIUS1',	 'LARMORQ3',       'LICHTENBER9',   'MRECRISIUM2',	 'PRINZVENT',     'SPARIM4',    \
                        'ARISTARCHU3',  'FEOKTISTOV',    'HORTENSIUS2',	 'LARMORQ4',       'LINNECRATER',   'MRECRISIUM3',	 'PYTHAGORAS1',   'SULPICIUS1', \
                        'ARISTARCHU4',  'FRESH1',	     'HORTENSIUS3',	 'LARMORQ5',       'LOWELL1',	    'MRINGENII1',	 'RANGER',	      'SULPICIUS2', \
                        'ARISTARCHU5',  'FRESH3',	     'HORTENSIUS5',	 'LASELMASIF1',    'LOWELL2',	    'MRINGENII2',	 'REINER1',       'SULPICIUS3', \
                        'ATLAS1',	    'FRESHMELT1',    'HORTENSIUS6',	 'LASELMASIF2',    'LUNA16',	    'MRINGENII3',	 'REINER2',       'VIRTANEN1',  \
                        'ATLAS2',	    'FRESHMELT2',    'HORTENSIUS7',	 'LASELMASIF3',    'LUNA20',	    'MRINGENII4',	 'REINER3',       'VIRTANEN2',  \
                        'ATLAS5',	    'FRSHCRATER10',  'IMBRIUM',	     'LASELMASIF5',    'LUNA24',	    'NEARLENTS1',	 'REINER4',       'VIRTANEN3',  \
                        'ATLAS6',	    'FRSHCRATER14',  'IMBRIUM2',	 'LICHTENBER1',    'LUNOKHOD2_1',   'NECHO2',	     'RUMKERDOME1',   'VITELLO',    \
                        'BHABHAPLAIN',  'FRSHCRATER15',  'IMPACTMELT1',	 'LICHTENBER10',   'LUNOKHOD2_2',   'NEWCRATER1',	 'RUMKERDOME2',   'WEIRDCRTR',  \
                        'CAUCHY',	    'FRSHCRATER4',   'IMPACTMELT2',	 'LICHTENBER11',   'MANILUS',	    'NRTHCRTRI1',	 'RUMKERDOME3', 
                        'CMPLXSCRP',    'FRSHCRATER5',   'IMPACTMELT3',	 'LICHTENBER12',   'MARIUS1',	    'NRTHCRTRI2',	 'RUMKERDOME4', 
                        'CMPTNBELK2',   'FRSHCRATER8',   'INACALDERA1',	 'LICHTENBER13',   'MARIUS2',	    'NRTHCRTRII',	 'RUMKERDOMES5',
                        'CMPTNBELK3',   'GRUITHUISE1',   'INACALDERA3',	 'LICHTENBER2',    'MARIUS3',	    'NRTHCRTRIII',	 'SEARESGRBN',  
                        'COMPTONBELK',  'GRUITHUISE2',   'ISISOSIRIS',	 'LICHTENBER3',    'MARIUS4',	    'OBLIVIONIS1',	 'SEISMCLAND',
                        'COPERNGRAB',   'GRUITHUISE3',   'KINGCRATER2',	 'LICHTENBER4',    'MOOREF1',	    'OPENHIMERF',	 'SLIPHER1',
                        'DALEMBERT',    'GRUITHUISE4',   'KINGCRATER3',	 'LICHTENBER5',    'MOOREF2',	    'ORIENTALE1',	 'SOSIRILLE',
                        'EIMMARTA',     'GRUITHUISE5',   'KINGCRATER4',	 'LICHTENBER6',    'MOROZOVE',	    'ORIENTALE2',	 'SPARIM1',
                        'ENDYMION',     'HANSTEENAL1',   'KUGLERRIDGE',	 'LICHTENBER7',    'MOSCOVNSE1',    'ORIENTALE3',	 'SPARIM2',
                        'ERATOSTHNS1',  'HIGHESTPOIN',   'LARMORQ1',	 'LICHTENBER8',    'MRECRISIUM1',   'PLANCKFLOOR',	 'SPARIM3' ]

    demFolderList.sort() # Sort these names for convenience

    # List of all the output files we want copied
    fileList = [  'results/ASU_LOLA_diff_points.csv', \
                  'results/ASU_LOLA_diff_points.kml', \
                  'results/ASU_LOLA_diff_stats.txt' ,\
                  'results/LOLA_diff_points.csv', \
                  'results/LOLA_diff_points.kml', \
                  'results/LOLA_diff_stats.txt' ,\
                  'results/output-Log.txt' ,\
                  'results/pcAlignLog.txt' ,\
                  'stdOutLog.txt' ,\
# TODO: Option to only grab the small output files
#                  'results/outputHillshade.tif' ,\
#                  'results/p2d-DEM.tif' ,\
#                  'results/p2d-IntersectionErr.tif',\
                  'workDir/refinement/lolaRdrPoints.kml' ,\
                  'workDir/refinement/inputGdcPoints.kml' ,\
                  'workDir/refinement/transformedGdcPoints.kml' ,\
                  'workDir/refinement/beg-errors.kml' ,\
                  'workDir/refinement/end-errors.kml' ,\
                  'workDir/refinement/pairGdcCheckFinal.kml' ,\
                  'workDir/refinement/pairGdcCheckFinalStereo.kml']

    try:
  
        if not os.path.exists(localFolder):
            os.makedirs(localFolder)

        dirListFile = os.path.join(localFolder, 'desiredFileList.txt')
        dirFile     = open(dirListFile, 'w')

        # Operate on each folder
        for f in demFolderList:

            # Get folders
            folderPath              = os.path.join(localFolder ,              f)
            supercomputerFolderPath = os.path.join(supercomputerSourceFolder, f)

            if not os.path.exists(folderPath):
                os.makedirs(folderPath)

            # For each folder, copy each result
            for c in fileList:
	            # Get file paths
                outputFilePath        = os.path.join(folderPath,              c)
                supercomputerFilePath = os.path.join(supercomputerFolderPath, c)

                #if not os.path.exists(outputFilePath) or force: # Request file copy
                #    IsisTools.grabSupercomputerFile(supercomputerFilePath, outputFilePath)    
                
                # Write the SC file path to the temp file
                relativeFilePath = os.path.join(f, c)
                dirFile.write(relativeFilePath + '\n')
        
        dirFile.close()
    
        # Use rsync to grab all the specified files at once
        #cmd = 'rsync -av --files-from=' + dirListFile + ' smcmich1@pfe22.nas.nasa.gov:/u/smcmich1/data/lronacPipeline/ ' + localFolder
        cmd = 'rsync -av --files-from=' + dirListFile + ' smcmich1@pfe22.nas.nasa.gov:/nobackupnfs2/oalexan1/scott/ ' + localFolder
        #cmd = 'rsync -av --files-from=' + dirListFile + ' smcmich1@pfe22.nas.nasa.gov:/u/smcmich1/data/lronacProduction/ ' + localFolder
        print cmd
        os.system(cmd)
    
    except Exception,e: # Catch any errors, the program will move on to the next folder
        print "Caught: ", e


# Read statistics from a file and append to lists
def accumulateStatistics(prefix, filePath, meanList, stdDevList, meanHistogram, dataStorage):

    try:
        meanValue, stdDev, percentile, histogram = IsisTools.readLolaCompareFile(filePath)

    except Exception,e: # Catch any errors, the program will move on to the next folder
        #print "Caught: ", e
        #print "Unable to process data in file " + filePath
        return False

    #print 'Read data from ' + filePath

    meanList.append(meanValue)
    stdDevList.append(stdDev)
  
    # Record the raw percentile data
    dataStorage[prefix] = percentile
  
    if not meanHistogram: # First histogram
        meanHistogram.extend(histogram)
    else: # Accumulate each histogram bin
        for o, n in zip(meanHistogram, histogram):
            o = o + n 

    return True # We successfully added data

#def readShiftAmounts(path, prefix, shiftDict):
#
#    # Make sure the file exists
#    if not os.path.exists(path):
#        return False
#
#    # Read the three numbers from the file
#    shiftFile = open(path, 'r')
#    line      = shiftFile.readline()
#    parts     = line.split(' ')
#    transform = (float(parts[0]), float(parts[1]), float(parts[2]))
#    shiftDict[prefix] = transform
#    
#    return True

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
    #bbStorage = dict()
    #parseBoundingBoxes(logPath, bbStorage)
    # ---> TODO: Get log file and restore this!

    usedFolderList = []
    foldersInDirectory = os.listdir(dataFolder)
    foldersInDirectory.sort()
    for f in foldersInDirectory:

        folderPath = os.path.join(dataFolder, f)

        # Get the file paths for this folder
        lolaDiffPath        = folderPath + '/results/LOLA_diff_stats.txt'
        lolaAsuDiffPath     = folderPath + '/results/ASU_LOLA_diff_stats.txt'
        #stereoDemPath       = folderPath + '/stereo-DEM.tif'
        
        if not os.path.exists(lolaDiffPath): # Ignore directories if the LOLA data is not there
            continue
        
        # Accumulate the statistics
        readLola    = accumulateStatistics(f+'_lola', lolaDiffPath,    lolaMeanList,    lolaStdDevList,    lolaMeanHistogram,    dataStorage )
        readLolaAsu = accumulateStatistics(f+'_comp', lolaAsuDiffPath, lolaAsuMeanList, lolaAsuStdDevList, lolaAsuMeanHistogram, dataStorage )

        # Check the mean LOLA error from the last result
        LOLA_ERROR_LIMIT = 50 # Max LOLA error before we consider the run a failure
        lastLolaError = lolaMeanList[-1]
        if lastLolaError >= LOLA_ERROR_LIMIT:
            print 'Folder ' + f + ' considered failure due to high mean LOLA error of ' + str(lastLolaError)
        
        if readLola and readLolaAsu and (lastLolaError < LOLA_ERROR_LIMIT):
            usedFolderList.append(f) # Keep track of the folders we read data from
            print 'Read folder ' + f
        else:
            # Remove list entries from partial successes
            if readLola:
                lolaMeanList.pop()
            if readLolaAsu:
                lolaAsuMeanList.pop()
                   
            if readLola or readLolaAsu:
                print 'Partial success for folder ' + f
            #else:
            #    print 'Failed to read in folder ' + f

    # TODO: Restore this once BB data read from the log file!
#    # Generate some google earth KML output
#    generateKml(dataFolder, usedFolderList, bbStorage, asuMeanList)
    

    # We now have a list of means/stds and sums of percentiles and histograms
    # - Get means for all of them
    numElements = float(len(lolaMeanList))
    lolaMeanOfMeans      = sum(lolaMeanList     ) / numElements
    lolaMeanOfStdDevs    = sum(lolaStdDevList   ) / numElements
    lolaAsuMeanOfMeans   = sum(lolaAsuMeanList  ) / numElements
    lolaAsuMeanOfStdDevs = sum(lolaAsuStdDevList) / numElements

    for i in lolaMeanHistogram:
      i = i / numElements

    for i in lolaAsuMeanHistogram:
      i = i / numElements

    # Compute standard deviation of the LOLA and ASU means
    lolaStdDev    = np.std(lolaMeanList)
    lolaAsuStdDev = np.std(lolaAsuMeanList)


    # Write results into a condensed file
    condensedDataPath = os.path.join(dataFolder, 'resultsSummary.csv')
    condensedFile     = open(condensedDataPath, 'w')
    print 'Writing condensed result file ' + condensedDataPath
    i = 0
    condensedFile.write('Data_set,  Us_vs_Lola,  ASU_vs_Lola,  Change\n')
    for f in usedFolderList:
        condensedFile.write(f + ', ' + str(lolaMeanList[i]) + ', ' + str(lolaAsuMeanList[i]) + ', ' + str(lolaMeanList[i] - lolaAsuMeanList[i]) + '\n')
        i = i + 1
    condensedFile.close()
    
    
    # Now generate plots

    numEls = 20 # Number of elements in percentile plots (leaving off last one containing error values)
    yMin   = 0
    yMax   = 10
    xAxis  = range(0, numEls)
    xAxis  = [x * 0.05 for x in xAxis]
    lolaMeanPercentile    = [0] * numEls
    lolaAsuMeanPercentile = [0] * numEls


    # us vs LOLA
    lolaComparisonPath = os.path.join(dataFolder, 'allPercentilesLOLA.png')
    for f in usedFolderList:
        # Plot data points for this folder
        dataStorage[f+'_lola'].pop() # Strip of error values in last bin
        plt.plot(xAxis, dataStorage[f+'_lola'], 'o', label=f)
        
        # Accumulate mean value for each percentile across data sets
        for i in range(0,numEls):
            lolaMeanPercentile[i] = lolaMeanPercentile[i] + dataStorage[f+'_lola'][i]
            
    # Finish and plot mean value
    for i in range(0,numEls): 
        lolaMeanPercentile[i] = lolaMeanPercentile[i] / len(usedFolderList)
    plt.plot(xAxis, lolaMeanPercentile, '-', color='r', linewidth=3, label='mean percentiles')
    
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('Us vs LOLA pixel percentages for all samples')
    #lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print 'Writing plot ' + lolaComparisonPath
    plt.savefig(lolaComparisonPath)#, bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()
    
    
    
    # LOLA vs ASU
    lolaAsuComparisonPath = os.path.join(dataFolder, 'allPercentilesLolaASU.png')
    for f in usedFolderList:
        # Plot data points for this folder
        dataStorage[f+'_comp'].pop() # Strip of error values in last bin
        plt.plot(xAxis, dataStorage[f+'_comp'], 'o', label=f)
        
        # Accumulate mean value for each percentile across data sets
        for i in range(0,numEls):
            lolaAsuMeanPercentile[i] = lolaAsuMeanPercentile[i] + dataStorage[f+'_comp'][i]
            
    # Finish and plot mean value
    for i in range(0,numEls): 
        lolaAsuMeanPercentile[i] = lolaAsuMeanPercentile[i] / len(usedFolderList)
    plt.plot(xAxis, lolaAsuMeanPercentile, '-', color='b', linewidth=3, label='mean percentiles')
    
    plt.grid(color='gray', linestyle='dashed')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.xlabel('Percent of pixels less than difference')
    plt.title('LOLA vs ASU pixel percentages for all samples')
    #lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print 'Writing plot ' + lolaAsuComparisonPath
    plt.savefig(lolaAsuComparisonPath)#, bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()    
    
    # Compute the standard deviations for LOLA and ASU in each percentile
    lolaPercentileStdList    = []
    lolaAsuPercentileStdList = []
    for i in range(0,numEls): # For each percentile
        lolaDiffSum    = 0
        lolaAsuDiffSum = 0
        for f in usedFolderList: # For each data set
            lolaDiff       = dataStorage[f+'_lola'][i] - lolaMeanPercentile[i]
            lolaAsuDiff    = dataStorage[f+'_comp'][i] - lolaAsuMeanPercentile[i]
            lolaDiffSum    = lolaDiffSum    + lolaDiff*lolaDiff
            lolaAsuDiffSum = lolaAsuDiffSum + lolaAsuDiff*lolaAsuDiff
        
        # This is the standard deviation for this percentile
        lolaPercentileStd    = lolaDiffSum       / len(usedFolderList)
        lolaAsuPercentileStd = lolaAsuDiffSum    / len(usedFolderList)
        print lolaPercentileStd
        
        lolaPercentileStdList.append(lolaPercentileStd)
        lolaAsuPercentileStdList.append(lolaAsuPercentileStd)
    
    lolaPercentileErrBars    = np.row_stack((lolaPercentileStdList,    lolaPercentileStdList))
    lolaAsuPercentileErrBars = np.row_stack((lolaAsuPercentileStdList, lolaAsuPercentileStdList))
    
    #print lolaPercentileStdList
    #print '---'
    #print lolaAsuPercentileStdList
    
    # Plot the mean percentiles on one chart
    (_, capsLola, _) = plt.errorbar(xAxis, lolaMeanPercentile   , color='r', label='us vs LOLA',  yerr=lolaPercentileErrBars,    linewidth=3, elinewidth=4, capsize=6)
    (_, capsAsu,  _) = plt.errorbar(xAxis, lolaAsuMeanPercentile, color='b', label='ASU vs LOLA', yerr=lolaAsuPercentileErrBars, linewidth=2, elinewidth=2, capsize=6)
    for cap in capsLola:
        cap.set_markeredgewidth(4)
    for cap in capsAsu:
        cap.set_markeredgewidth(4)
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
    yMax = 10
    xAxis = np.arange(len(lolaMeanList))
    barwidth = 0.4
    # Sort the values according to the ASU error for this plot
    zippedValues = zip(lolaAsuMeanList, lolaMeanList, usedFolderList)
    zippedValues.sort()
    lolaAsuMeanList, lolaMeanList, usedFolderList = zip(*zippedValues)
    plt.bar(xAxis+1*barwidth, lolaMeanList,    barwidth, linewidth=0, color='r', label='us vs LOLA')
    plt.bar(xAxis+2*barwidth, lolaAsuMeanList, barwidth, linewidth=0, color='b', label='ASU vs LOLA')
    #plt.xticks(xAxis+.5, usedFolderList, size='small', rotation='vertical')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Mean difference between DEMs')
    plt.savefig(dataFolder + '/means.png', bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()


#--------------------------------------------------------------------------------

def main():

    try:
        try:
            usage = "usage: lronacPipeline.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--folder", dest="folder",
                              help="Data storage folder")
            parser.add_option("--force", action="store_true", dest="force",
                              help="Force overwrite of existing files.")
            parser.add_option("--get-results", action="store_true", dest="getResults",
                              help="Retrieve remote result files.")
            parser.add_option("--make-plots", action="store_true", dest="makePlots",
                              help="Generates plots.")
                
            (options, args) = parser.parse_args()

            if not options.folder:
                raise Exception('Must pass in the data folder!')

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        # All data is stored in this directory
        dataFolder = options.folder


        if options.getResults:
            grabResultFiles(dataFolder, options.force)
#        backupData(dataFolder)
        if options.makePlots:
            generatePlots(dataFolder)


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
