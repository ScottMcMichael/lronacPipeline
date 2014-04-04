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

import IsisTools, IrgFileFunctions, IrgIsisFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''Tool for managing lronac pipeline production data'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



#--------------------------------------------------------------------------------

# All status information is saved in this folder
PROCESSING_FOLDER      = '/u/smcmich1/data/lronacProduction' 
STATUS_FOLDER          = '/u/smcmich1/projects/lronacPipeline/dataStatus'
INCOMPLETE_STATUS_PATH = os.path.join(STATUS_FOLDER, 'incompleteStatus.txt')
COMPLETE_STATUS_PATH   = os.path.join(STATUS_FOLDER, 'completeStatus.txt')
ARCHIVE_STATUS_PATH    = os.path.join(STATUS_FOLDER, 'archiveStatus.txt')
LOU_STORAGE_PATH       = 'lfe4:/lou/s2i/smcmich1/LRONAC_DTM'


def folderToDataSetName(folder):
        """Given a folder, get the data set name"""
        
        justFolder = os.path.basename(folder)
        if not justFolder.startswith('NAC_DTM_'):
            raise Exception('Folder does not contain a data set!')
        
        dataSetName = justFolder[8:]
        return dataSetName

class DataSet:
    """Class describing the status of a data set"""
    # Constants
    INVALID    = 0
    INCOMPLETE = 1
    COMPLETE   = 2
    ARCHIVED   = 3
    
    # Variables
    status          = INVALID
    meanLolaError   = 0.0
    stdDevLolaError = 0.0
    name            = ''
    
    def __init__(self, name, evaluateStatus=True):
        """Constructor"""
        self.name = name
        if evaluateStatus:
            self.evaluateStatus()

    def evaluateStatus(self):
        """Fills in the class information and returns true the data folder was found"""
        
        # Get folder name and make sure it exists    
        folder = self.getLocalFolder()
        if not os.path.exists(folder):
            self.status = self.INVALID
            return False
        
        if self._isFolderArchived(folder):
            self.status = self.ARCHIVED
            
        else: # Data is still local, see if it is completed.
            
            # List of the files that need to be present to be considered complete
            requiredOutputList = [ 'output-CompressedInputs.tar.bz2',
                                   'output-CompressedOutputs.tar.bz2']
        
            # If all the required files are present the status is complete.
            self.status = self.COMPLETE
            for f in requiredOutputList:
                fullPath = os.path.join(folder, 'results/' + f)
                # If any file is missing return incomplete status
                if not os.path.exists(fullPath):
                    status.status = DataSetStatus.INCOMPLETE
                
        # Fill in some additional information if available
        statsPath = os.path.join(folder, 'results/output-LOLA_diff_stats.txt')
        if os.path.exists(statsPath):
            meanValue, stdDev, percentile, histogram = IsisTools.readLolaCompareFile(filePath)
            status.meanLolaError   = meanValue
            status.stdDevLolaError = stdDev

        return True

    def _isFolderArchived(self, folder):
        """Returns true if the folder has been archived and cleared.\n
           - This is indicated by a folder that is almost empty except for 'downloadLog.txt'"""
           
        # TODO: Verify the number of entries on this list!
        itemList = os.listdir(folder)
        return (len(itemList) == 3) and (itemList[0] == 'downloadLog.txt')

    def getLocalFolder(self):
        """Get the local storage location"""
    
        folderName = 'NAC_DTM_' + self.name       
        fullPath   = os.path.join(PROCESSING_FOLDER, folderName)
        return fullPath
    
    def getPbsJobName(self):
        """Returns the abbreviated data set name used with PBS"""
        
        # The short name is the last 7 characters from each name
        spacerPos = name.find('_')
        shortName = name[spacerPos-7:spacerPos+1] + name[-7:]
        return shortName
        
    
    def getStatusString():
        """Returns a string describing the data set"""
        
        s = (self.name + ' ' +
              str(self.meanLolaError) + ' ' + str(self.stdDevLolaError) )
        return s
    
    def isOnLocalDisk(self):
        return (self.status == self.INCOMPLETE) or (self.status == self.COMPLETE)
    
    def isComplete(self):
        return self.status == self.COMPLETE
    
def findPbsMatch(pbsName):
    """Return a list of possible paths for the PBS job name"""

    results = []

    # Loop through all local data folders
    foldersInDirectory = os.listdir(PROCESSING_FOLDER)
    foldersInDirectory.sort()
    for f in foldersInDirectory:
        try:
            dataSetName = folderToDataSetName(f)
        except: # Skip folders that do not contain a data set
            continue
        
        # Get the status of this data set
        ds = DataSet(dataSetName)
        
        if ds.getPbsJobName() == pbsName:
            results.append(ds.getLocalFolder())
    
    return results

def recordLocalDataStatus(dryRun=False):
    """Records the status of all locally stored data sets\n
       - Currently this generates a list of incomplete and complete files"""
    
    if not dryRun:
        # Replace the existing local status files
        incompleteFile = open(INCOMPLETE_STATUS_PATH, 'w')
        completeFile   = open(COMPLETE_STATUS_PATH,   'w')
        archiveFile    = open(ARCHIVE_STATUS_PATH,    'w')
    
    # Loop through all local data folders
    foldersInDirectory = os.listdir(PROCESSING_FOLDER)
    foldersInDirectory.sort()
    for f in foldersInDirectory:
        try:
            dataSetName = folderToDataSetName(f)
        except: # Skip folders that do not contain a data set
            continue
        
        # Get the status of this data set
        ds = DataSet(dataSetName)
        
        if ds.status == DataSet.INVALID:
            continue # Skip invalid data sets

        dataSetStatusString = ds.getStatusString()

        if ds.status == DataSet.ARCHIVED:
            if dryRun:
                print ' - ARCHIVED file - ' + dataSetStatusString
            else: # Add to log file file
                archiveFile.write(dataSetStatusString + '\n') 
        
        if ds.status == DataSet.COMPLETE:
            if dryRun:
                print ' - COMPLETE file - ' + dataSetStatusString
            else: # Add to log file file
                completeFile.write(dataSetStatusString + '\n')

        if ds.status == DataSet.INCOMPLETE:
            if dryRun:
                print ' - INCOMPLETE file - ' + dataSetStatusString
            else: # Add to log file file
                incompleteFile.write(dataSetStatusString + '\n')
        
    if not dryRun:
        # Finished writing the output files
        completeFile.close()
        incompleteFile.close()
        archive.close()
    
    

def flagIncompleteDataSets(dryRun=False):
    """Flags incomplete data sets for repeat processing"""
    
    # Loop through all data sets in the incomplete folder
    incompleteFile = open(INCOMPLETE_STATUS_PATH, 'r')
    for line in incompleteFile:
        entries     = line.split(' ')
        dataSetName = entries[0]
    
        # Get the status of this data set
        ds = DataSet(dataSetName)
        
        # Check that it is actually incomplete
        if ds.status != DataSet.INCOMPLETE:
            print 'WARNING: non-incomplete set found in incomplete file:'
            print line
            continue
        
        # For incomplete data sets, remove the downloadLog.txt file
        # - This will cause the auto-processor code to run the data set again.
        flagFile = os.path.join(f, 'downloadLog.txt')
        if dryRun:
            print '- Remove file ' + flagFile
        else:  # Actually removet the file
            IrgFileFunctions.removeIfExists(flagFile)
    
    # Need to ru-run the update function to update this file
    incompleteFile.close()

    
    
def archiveDataSet(dataSetName, deleteLocalFiles=False, dryRun=False):
    """Archives a data set to long term storage on Lou"""
    
    # The local files to be archived are already in tar files
    
    # Make sure the data set is complete before archiving it
    ds = DataSet(dataSetName)
    if not ds.isComplete():
        raise Exception('Attempting to archive incomplete data set: ' + dataSetName)
    localFolder = ds.getLocalFolder()
    
    # In order to archive the file it just has to be moved over to the Lou filesystem.
    # - Tape archiving is performed automatically when needed.
    # - All tar files are stored in the same directory
    filesToArchive = ['output-CompressedOutputs.tar.bz2']
    
    for f in filesToArchive:
        fullName  = 'NAC_DTM_' + dataSetName + '_'
        inputPath = os.path.join(localFolder, fullName)
        outputPath
        cmd = 'shiftc ' + fullPath + ' ' + LOU_STORAGE_PATH
        if dryRun:
            print '- ' + cmd
        else: # Actually transfer the file
            print cmd
            os.system(cmd)   
    
    if deleteLocalFiles: # Remove almost all files
        
        # Keep some result files and move them to the main data set directory
        filesToKeep   = ['output-LOLA_diff_stats.txt']
        resultsFolder = os.path.join(localFolder, 'results/')
        for f in filesToKeep:
            inputPath  = os.path.join(resultsFolder, f)
            outputPath = os.path.join(localFolder,   f)
            if dryRun:
                print 'mv ' + inputPath + ' ' + outputPath
            else: # Actually move the files
                shutil.move(inputPath, outputPath)
        
        # Now that we have saved a few files, delete everything else
        
        cmdRmImg   = 'rm ' + localFolder + '/*.IMG'
        rdrPath    = os.path.join(localFolder + 'lolaRdrPoints.csv')
        prtPath    = os.path.join(localFolder + 'print.prt')
        logPath    = os.path.join(localFolder + 'stdOutputLog.txt')
        workDir    = os.path.join(localFolder + '/workDir')
        resultsDir = os.path.join(localFolder + '/results')
        
        if dryRun:
            print '- ' + cmdRmImg
            print '- rm ' + rdrPath
            print '- rm ' + ptrPath
            print '- rm ' + logPath
            print '- rm ' + workDir
            print '- rm ' + resultsDir
        else: # Actually delete the files
            os.system(cmdRmImg)
            IrgFileFunctions.removeIfExists(rdrPath)
            IrgFileFunctions.removeIfExists(ptrPath)
            IrgFileFunctions.removeIfExists(logPath)
            IrgFileFunctions.removeFolderIfExists(workDir)
            IrgFileFunctions.removeFolderIfExists(resultsDir)
        
    
    
    

def archiveAllCompletedResults(deleteLocalFiles=False, dryRun=False):
    """Archives all the folders in the completed results archive folder"""

    # Loop through all data sets in the incomplete folder
    completeFile = open(COMPLETE_STATUS_PATH, 'r')
    if not dryRun:
        archiveFile  = open(ARCHIVE_STATUS_PATH,  'a')
    for line in completeFile:
        entries     = line.split(' ')
        dataSetName = entries[0]
        
        # Get the status of this data set
        ds = DataSet(dataSetName)
        
        # Check that it is actually incomplete
        if ds.status != DataSet.INCOMPLETE:
            print 'WARNING: non-incomplete set found in incomplete file:'
            print line
            continue
        
        archiveDataSet(dataSetName, deleteLocalFiles, dryRun)
        if not dryRun:
            archiveFile.write(line)

    completeFile.close()
    if not dryRun:
        archiveFile.close()
    
    if not dryRun:
        # Clear the completed results file
        completeFile = open(COMPLETE_STATUS_PATH, 'w')
        completeFile.close()
    
    
    

def restoreArchivedDataSet(dataSetName):
    """Restores an archived data set from Lou"""
    print 'TODO'
    
#def lookupArchivedDataSets():
#    """Make a list of all the data sets currently stored on Lou"""
#    print 'TODO'
    

def main():


    try:
        try:
            usage = "usage: archiveManager.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
            parser.add_option("--find-pbs", dest="pbsName",
                              help="Get the full path to the data set for a PBS job.",
                              default=None)

            parser.add_option("--check-local", action="store_true",
                              dest="checkLocal", default=False,
                              help="Update local data set listings")
            
            parser.add_option("--archive-results", action="store_true",
                              dest="archive", default=False,
                              help="Archive completed files to Lou.")

            # TODO: Archive restore function

            parser.add_option("--flag-incomplete", action="store_true",
                              dest="flagIncomplete", default=False,
                              help="Flag incomplete data sets for repeat processing.")
            
            parser.add_option("--clear", action="store_true", dest="clear", default=False,
                              help="Clear files after archiving them.")
            
            parser.add_option("--dry-run", action="store_true", dest="dryRun", default=True, #TODO: Update
                              help="Don't touch any data, just display actions to be taken.")
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

        except optparse.OptionError, msg:
            raise Usage(msg)

        if options.pbsName: # Search for all folders this PBS name may apply to
            outputList = findPbsMatch(options.pbsName)
            print outputList
            return 0 # Don't combine with other commands
            

        if options.checkLocal:
            recordLocalDataStatus(options.dryRun)

        if options.archive:
            print 'TODO!'
        
        if options.flagIncomplete:
            flagIncompleteDataSets(options.dryRun)

        
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())