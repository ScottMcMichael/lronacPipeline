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

import os, glob, optparse, re, shutil, subprocess, sys, string, time, urllib, urllib2

import matplotlib.pyplot as plt
import numpy as np

job_pool = [];

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
This program operates on LRO (.IMG) files, and performs the

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

# Copies one file from the supercomputer - Caller needs to wait for the job to finish!
def grabFile(supercomputerPath, localPath):

    numThreads = 1 #TODO

    cmd = 'sup scp smcmich1@bridge3.nas.nasa.gov:' + supercomputerPath  + ' ' + localPath
#    print cmd #TESTING
    add_job(cmd, numThreads)


# Grabs the output files from the supercomputer
def grabResultFiles(outputFolder):

    supercomputerSourceFolder = '/nobackupp1/smcmich1/data/lronacPipeline'

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


    # List of all the temporary files we want deleted
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
  
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)

        # Operate on each folder
        for f in demFolderList:
	
            # Get folders
            folderPath              = os.path.join(outputFolder,              f)
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
                    grabFile(supercomputerFilePath, outputFilePath)

            fileList.pop() # Remove the last ASU input DEM from the list
		            
        
        wait_on_all_jobs() # Wait on all file grab requests
        
    except Exception,e: # Catch any errors, the program will move on to the next folder
        print "Caught: ", e


# Reads the results from a single diff_stats.txt file
def readStatsFile(filePath):

    # Read in the output file to extract the CenterLatitude value
    meanValue  = -32768
    stdDev     = -32768
    percentile = []
    histogram  = []

    statsFile = open(filePath, 'r')
    for line in statsFile:
        if (line.find('Mean') >= 0):         # Parse out the mean
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            meanValue = float(numStr)
        elif (line.find('deviation') >= 0):  # Parse out the standard deviation
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            stdDev = float(numStr)
        elif (line.find('Percentile') >= 0): # Add a line to the percentile distribtion
            eqPt   = line.find('=')
            numStr = line[eqPt+2:]
            percentile.append(float(numStr))
        elif (line.find('<-->') >= 0): # Add a line to the histogram
            eqPt   = line.rfind('=') # Get the second equal sign
            perPt  = line.find('%')
            numStr = line[eqPt+1:perPt-1]
            histogram.append(float(numStr))

    # Make sure we found the desired values
    if (meanValue == -32768) or (stdDev == -32768):
        raise Exception("Unable to find statistics in file " + filePath)
    
    return (meanValue, stdDev, percentile, histogram)


def accumulateStatistics(prefix, filePath, meanList, stdDevList, meanHistogram, dataStorage):

    try:
        meanValue, stdDev, percentile, histogram = readStatsFile(filePath)

    except Exception,e: # Catch any errors, the program will move on to the next folder
        #print "Caught: ", e
        #print "Unable to process data in file " + filePath
        return False

    print 'Read data from ' + filePath

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
    
    dataStorage = dict()

    usedFolderList = []
    for f in os.listdir(dataFolder):

        folderPath = os.path.join(dataFolder, f)
	
        # Get the file paths for this folder
        asuDiffPath     = folderPath + '/ASU_diff_stats.txt'
        lolaDiffPath    = folderPath + '/LOLA_diff_stats.txt'
        lolaAsuDiffPath = folderPath + '/LOLA_ASU_diff_stats.txt'
               
        # Accumulate the statistics
        readAsu     = accumulateStatistics(f+'_asu',  asuDiffPath,     asuMeanList,     asuStdDevList,     asuMeanHistogram,     dataStorage)
        readLola    = accumulateStatistics(f+'_lola', lolaDiffPath,    lolaMeanList,    lolaStdDevList,    lolaMeanHistogram,    dataStorage )
        readLolaAsu = accumulateStatistics(f+'_comp', lolaAsuDiffPath, lolaAsuMeanList, lolaAsuStdDevList, lolaAsuMeanHistogram, dataStorage )

        if readAsu and readLola and readLolaAsu:
            usedFolderList.append(f) # Keep track of the folders we read data from
        else:
            if readAsu or readLola or readLolaAsu:
                print 'Partial success for folder ' + f
    

    #TODO: Back out local vertical/horizontal offset, put on scatter plot

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

#TODO: Different version of this plot for large data sets
   
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
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesASU.png')
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
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesLOLA.png')
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
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(dataFolder + '/allPercentilesLolaASU.png')
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
    barwidth = 0.3
    plt.bar(xAxis,            asuMeanList,     barwidth, color='r', label='us vs ASU')
    plt.bar(xAxis+1*barwidth, lolaMeanList,    barwidth, color='g', label='us vs LOLA')
    plt.bar(xAxis+2*barwidth, lolaAsuMeanList, barwidth, color='b', label='ASU vs LOLA')
    plt.xticks(xAxis+.5, usedFolderList, size='small', rotation='vertical')
    plt.ylim(yMin, yMax)
    plt.ylabel('Difference in meters')
    plt.legend(loc='upper right')
    plt.title('Mean difference between DEMs')
    plt.savefig(dataFolder + '/means.png')
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

        dataFolder = '/home/smcmich1/data/lronacResults'
        grabResultFiles(dataFolder)
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
