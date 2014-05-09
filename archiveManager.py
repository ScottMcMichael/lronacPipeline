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

import os, glob, optparse, re, shutil, subprocess, string, time, pipes

import workerpool

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
LOU_LOCAL_STORAGE_PATH = '/lou/s2i/smcmich1/LRONAC_DTM'
LOU_STORAGE_PATH       = 'lfe4:' + LOU_LOCAL_STORAGE_PATH
LOU_SSH_TARGET         = 'smcmich1@lfe4'

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
    
    def __init__(self, name, evaluateStatus=True, verbose=False):
        """Constructor"""
        self.name = name
        if evaluateStatus:
            self.evaluateStatus(verbose)

    def evaluateStatus(self, verbose=False):
        """Fills in the class information and returns true the data folder was found"""
        
        # Get folder name and make sure it exists    
        folder = self.getLocalFolder()
        if not os.path.exists(folder):
            self.status = self.INVALID
            return False
        
        if self._isFolderArchived(folder):
            self.status = self.ARCHIVED

        elif self._hasErrors(folder, verbose): # Data is local, but may have errors
            self.status = self.INCOMPLETE
            # This check is important since not all errors are fatal!
            
        else: # Data is still local, see if it is completed.
        
            # List of the files that need to be present to be considered complete
            requiredOutputList = [ 'output-CompressedInputs.tar.bz2',
                                   'output-CompressedOutputs.tar']
        
            # If all the required files are present the status is complete.
            self.status = self.COMPLETE
            for f in requiredOutputList:
                fullPath = os.path.join(folder, 'results/' + f)
                # If any file is missing return incomplete status
                if not os.path.exists(fullPath):
                    self.status = DataSet.INCOMPLETE
                
        # Fill in some additional information if available
        statsPath = os.path.join(folder, 'results/output-LOLA_diff_stats.txt')
        if os.path.exists(statsPath):
            meanValue, stdDev, percentile, histogram = IsisTools.readLolaCompareFile(statsPath)
            self.meanLolaError   = meanValue
            self.stdDevLolaError = stdDev

        return True

    def _isFolderArchived(self, localFolder):
        """Returns true if the folder has been archived and cleared.\n
           - This is indicated by a folder that has certain files moved into it"""
           
        itemList = os.listdir(localFolder)
        return ('output-CompressedInputs.tar.bz2' in itemList) and ('downloadLog.txt' in itemList)

    def _hasErrors(self, localFolder, verbose):
        """Searches for errors in the log file"""
        
        # Open the log file
        logPath = os.path.join(localFolder, 'stdOutLog.txt')
        if not os.path.exists(logPath): # If the log isn't there just return true to mark as incomplete
            return True
        logFile = open(logPath, 'r')
        
        # List of errors we know to look for
        errorList = ['Traceback (most recent call last):',
                     'caught an exception',
                     ' ERROR**',
                     'Disk quota exceeded',
                     'Failed to read',
                     'Unable to open',
                     'file size exceeded',
                     'unrecognised option']
        
        # Search for each of these errors in the log file
        for line in logFile:
            for e in errorList:
                if (e in line): # Error was found, print a big message and return
                    if verbose:
                        print '+++++ Found an error in the following line: +++++'
                        print line
                        print ' in file: ' + logPath
                        print '+++++++++++++++++++++++++++++++++++++++++++++++++\n'
                    logFile.close()
                    return True
        
        logFile.close()
        return False # No errors found

    def getLocalFolder(self):
        """Get the local storage location"""
    
        folderName = 'NAC_DTM_' + self.name       
        fullPath   = os.path.join(PROCESSING_FOLDER, folderName)
        return fullPath
    
    def getPbsJobName(self):
        """Returns the abbreviated data set name used with PBS"""
        
        # The short name is the last 7 characters from each name
        spacerPos = self.name.find('_')
        shortName = self.name[spacerPos-7:spacerPos+1] + self.name[-7:]
        return shortName
        
    
    def getStatusString(self):
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
        print 'Writing new status files based on folders in ' + PROCESSING_FOLDER
    
    # Loop through all local data folders
    foldersInDirectory = os.listdir(PROCESSING_FOLDER)
    foldersInDirectory.sort()
    completeCount   = 0
    incompleteCount = 0
    archivedCount   = 0
    for f in foldersInDirectory:
        try:
            dataSetName = folderToDataSetName(f)
        except: # Skip folders that do not contain a data set
            continue
        
        # Get the status of this data set
        ds = DataSet(dataSetName, verbose=True)
        
        if ds.status == DataSet.INVALID:
            continue # Skip invalid data sets

        dataSetStatusString = ds.getStatusString()

        if ds.status == DataSet.ARCHIVED:
            archivedCount = archivedCount + 1
            if dryRun:
                print ' - ARCHIVED file   - ' + dataSetStatusString
            else: # Add to log file file
                archiveFile.write(dataSetStatusString + '\n') 
        
        if ds.status == DataSet.COMPLETE:
            completeCount = completeCount + 1
            if dryRun:
                print ' - COMPLETE file   - ' + dataSetStatusString
            else: # Add to log file file
                completeFile.write(dataSetStatusString + '\n')

        if ds.status == DataSet.INCOMPLETE:
            incompleteCount = incompleteCount + 1
            if dryRun:
                print ' - INCOMPLETE file - ' + dataSetStatusString
            else: # Add to log file file
                incompleteFile.write(dataSetStatusString + '\n')
        
    print 'Found ' + str(completeCount)   + ' complete data sets'
    print 'Found ' + str(incompleteCount) + ' incomplete data sets'
    print 'Found ' + str(archivedCount)   + ' archived data sets'

    if not dryRun:
        # Finished writing the output files
        completeFile.close()
        incompleteFile.close()
        archiveFile.close()
    
    

def flagIncompleteDataSets(dryRun=False):
    """Flags incomplete data sets for repeat processing"""
    
    if not os.path.exists(INCOMPLETE_STATUS_PATH):
        raise Exception('Incomplete file list not present, run --check-local to generate it.')

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
        flagFile = os.path.join(ds.getLocalFolder(), 'downloadLog.txt')
        if dryRun:
            print '- Remove file ' + flagFile
        else:  # Actually removet the file
            IrgFileFunctions.removeIfExists(flagFile)
    
    # Need to ru-run the update function to update this file
    incompleteFile.close()

    
    
def removeCompressedOutputs(dataSetName, dryRun=False):
    """Removes the large local output files in a completed run"""
       
    # Make sure the data set is complete before removing any files
    ds = DataSet(dataSetName)
    if not ds.isComplete():
        raise Exception('Attempting to clear compressed on  incomplete data set: ' + dataSetName)
    localFolder = ds.getLocalFolder()
    compressedFilePath = os.path.join(localFolder, 'results/output-CompressedOutputs.tar')
    if not os.path.exists(compressedFilePath):
        raise Exception('Attempting to clear compressed but output TAR file is missing!')
    
    # Each of these files will be deleted
    filesToDelete = ['output-Colormap.tif',
                     'output-Confidence.tif',
                     'output-DEM.tif',
                     'output-Hillshade.tif',
                     'output-MapProjLeft.tif',
                     'output-MapProjRight.tif',
                     'output-MapProjLeftUint8.tif',  # Files after her are not compressed but can
                     'output-MapProjRightUint8.tif', #  quickly be regenerated from compressed files.
                     'output-IntersectionErrorX.tif',
                     'output-IntersectionErrorY.tif',
                     'output-IntersectionErrorZ.tif']
    
    for f in filesToDelete: # For each file in the list above
        inputPath = os.path.join(localFolder, 'results/' + f)
        
        if dryRun: # Print output action
            print '- rm ' + inputPath
        else: # Actually delete the file
            IrgFileFunctions.removeIfExists(inputPath)
   
   

def removeAllCompletedCompressedOutputs(dryRun=False):
    """Removes the large local output files in all completed runs"""

    if not os.path.exists(COMPLETE_STATUS_PATH):
        raise Exception('Complete file list not present, run --check-local to generate it.')

    # Loop through all data sets in the complete folder
    completeFile = open(COMPLETE_STATUS_PATH, 'r')
    for line in completeFile:
        entries     = line.split(' ')
        dataSetName = entries[0]
        
        removeCompressedOutputs(dataSetName, dryRun)

    completeFile.close()

   
def archiveDataSet(dataSetName, deleteLocalFiles=False, dryRun=False):
    """Archives a data set to long term storage on Lou"""
    
    # The local files to be archived are already in tar files
    
    # Make sure the data set is complete before archiving it
    ds = DataSet(dataSetName)
    if not ds.isComplete():
        raise Exception('Attempting to archive incomplete data set: ' + dataSetName)
    localFolder   = ds.getLocalFolder()
    resultsFolder = os.path.join(localFolder, 'results/')
    
    # In order to archive the file it just has to be moved over to the Lou filesystem.
    # - Tape archiving is performed automatically when needed.
    # - All tar files are stored in the same directory
    filesToArchive = ['output-CompressedOutputs.tar']
    
    for f in filesToArchive:
        # Get storage name to use
        archiveName  = 'NAC_DTM_' + dataSetName + '_CompressedOutputs.tar'
        inputPath    = os.path.join(resultsFolder, f)
        outputPath   = os.path.join(LOU_STORAGE_PATH, archiveName)
        
        # Move the file over to Lou
        #cmd = ['shiftc', '--wait', '--verify', inputPath,  outputPath]
        cmd = 'shiftc --wait --verify ' + inputPath + ' ' + outputPath
        print cmd
        if not dryRun: # Actually transfer the file
            os.system(cmd)
            #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            #textOutput, err = p.communicate()
            #print textOutput
            #print 'ppp'
            #print err
            #print textOutput.find('done')
            #if not ('done' in textOutput):
            #    raise Exception('Error transferring file!')
    
            # Verify that the file was actually transferred
            resp = subprocess.call(['ssh', LOU_SSH_TARGET, 'test -e ' + pipes.quote(LOU_LOCAL_STORAGE_PATH)])
            if resp == 0:
                print 'shiftc file transfer confirmed.'
            else:
                raise Exception('Failed to transfer file using shiftc!')

    if deleteLocalFiles: # Remove almost all files
        
        # Keep some result files and move them to the main data set directory
        filesToKeep   = ['output-LOLA_diff_stats.txt',
                         'output-Log.txt',
                         'output-PcAlignLog.txt',
                         'output-CompressedInputs.tar.bz2',
                         'output-CompressedDiagnostics.tar.bz2']
        for f in filesToKeep:
            inputPath  = os.path.join(resultsFolder, f)
            outputPath = os.path.join(localFolder,   f)
            print 'mv ' + inputPath + ' ' + outputPath
            if not dryRun: # Actually move the files
                shutil.move(inputPath, outputPath)
        
        # Now that we have saved a few files, delete everything else
        
        cmdRmImg   = 'rm ' + localFolder + '/*.IMG'
        rdrPath    = os.path.join(localFolder, 'lolaRdrPoints.csv')
        prtPath    = os.path.join(localFolder, 'print.prt')
        #logPath    = os.path.join(localFolder, 'stdOutputLog.txt')
        workDir    = os.path.join(localFolder, 'workDir')
        resultsDir = os.path.join(localFolder, 'results')
        
        print cmdRmImg
        print 'rm ' + rdrPath
        print 'rm ' + prtPath
        #print 'rm ' + logPath
        print 'rm -rf ' + workDir
        print 'rm -rf ' + resultsDir
        if not dryRun: # Actually delete the files
            os.system(cmdRmImg)
            IrgFileFunctions.removeIfExists(rdrPath)
            IrgFileFunctions.removeIfExists(prtPath)
            #IrgFileFunctions.removeIfExists(logPath)
            IrgFileFunctions.removeFolderIfExists(workDir)
            IrgFileFunctions.removeFolderIfExists(resultsDir)
    
    

def archiveAllCompletedResults(deleteLocalFiles=False, dryRun=False):
    """Archives all the folders in the completed results archive folder"""

    if not os.path.exists(COMPLETE_STATUS_PATH):
        raise Exception('Complete file list not present, run --check-local to generate it.')

    # Loop through all data sets in the complete folder
    completeFile = open(COMPLETE_STATUS_PATH, 'r')
    if not dryRun:
        archiveFile  = open(ARCHIVE_STATUS_PATH,  'a')
    for line in completeFile:
        entries     = line.split(' ')
        dataSetName = entries[0]
        
        # Get the status of this data set
        ds = DataSet(dataSetName)
        
        # Check that it is actually complete
        if not ds.isComplete():
            print 'WARNING: incomplete set found in complete file:'
            print dataSetName
            continue
        
        archiveDataSet(dataSetName, deleteLocalFiles, dryRun)
        if not dryRun:
            archiveFile.write(line)
    
    # Finished archiving data, now update files.

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


def getDataSetListFromFile(filePath):
    """Gets a list of data sets from one of the status files.
    
       Call one of these: getDataSetListFromFile(COMPLETE_STATUS_PATH)
                          getDataSetListFromFile(INCOMPLETE_STATUS_PATH)
                          getDataSetListFromFile(ARCHIVE_STATUS_PATH)
    """
    
    if not os.path.exists(filePath):
        raise Exception('File list ' + filePath + ' not present, run --check-local to generate it.')

    # Loop through all data sets in the log file
    inputFile = open(filePath, 'r')
    outputList = []
    for line in inputFile:
        entries     = line.split(' ')
        dataSetName = entries[0]
        outputList.append(dataSetName)
    inputFile.close
    return outputList


# TODO: replace with a generic function!
def fixCompletedDataSets(numThreads):
    """Applies a specific fix to completed local data sets after an output format change"""
    
    # Get a list of all the completed data sets
    dataSetList = getDataSetListFromFile(COMPLETE_STATUS_PATH)

    
    # Make a pool with multiple worker threads
    pool = workerpool.WorkerPool(size=numThreads)
    
    # Perform the mapping
    pool.map(fixCompleteDataSet, dataSetList)
    
    # TODO: How are exceptions handled?
    
    # Send shutdown jobs to all threads, and wait until all the jobs have been completed
    pool.shutdown()
    pool.wait()
    
    
def fixCompleteDataSet(dataSetName):
    """Applies a specific fix to completed local data sets after an output format change"""

    # Make sure the data set is complete before removing any files
    ds = DataSet(dataSetName)
    if not ds.isComplete():
        #raise Exception('Attempting to fix incomplete data set: ' + dataSetName)
        print 'Cannot fix incomplete data set: ' + dataSetName
        return
    
    localFolder   = ds.getLocalFolder()
    resultsFolder = os.path.join(localFolder, 'results/')
    outputTarPath = os.path.join(resultsFolder, 'output-CompressedOutputs.tar.bz2')
    
    print 'Fixing ' + localFolder
    if not os.path.exists(outputTarPath):
        print 'Skipping folder because bz2 tar file is not present'
        return

    # Check if the local files are still there (two should be a safe check)
    confPath = os.path.join(resultsFolder, 'output-Confidence.tif')
    demPath  = os.path.join(resultsFolder, 'output-DEM.tif')
    localFilesPresent = os.path.exists(confPath) and os.path.exists(demPath)
    
    if not localFilesPresent: # Unpack them from the output TAR file!
        print 'unpack ' + outputTarPath
        os.chdir(resultsFolder) # First move to the output directory
        cmd = 'tar -xjvf ' + outputTarPath
        print cmd
        os.system(cmd)
          
    # Each of these files will be deleted
    filesToDelete = ['output-Colormap.tif',
                     'output-Colormap.LBL',
                     'output-Confidence.tif',
                     'output-Confidence.LBL']
    
    # Delete all the files that need to be regenerated
    for f in filesToDelete: # For each file in the list above
        path = os.path.join(resultsFolder, f)
        print 'Delete file: ' + path
        os.remove(path)
        
    # Now call the main processing script on the file
    # - For fast jobs we don't need to use PBS
    cmd = 'jobWrapperV2.sh ' + localFolder
    print cmd
    os.system(cmd)
    
    # The main processing script will take care of all the work, including
    #  packing everything back up in to a TAR file.
   
    # Now clean up the large files we decompressed.
    print 'Removing large backed up files...'
    removeCompressedOutputs(dataSetName, dryRun=False)
    
    # Remove the old TAR file    
    print 'Deleting ' + outputTarPath
    os.remove(outputTarPath)

    print 'Done fixing ' + localFolder



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

            parser.add_option("--clear-compressed", action="store_true",
                              dest="clearCompressed", default=False,
                              help="Clear all local files in completed folders which have been compressed.")

            # TODO: Archive restore function

            parser.add_option("--flag-incomplete", action="store_true",
                              dest="flagIncomplete", default=False,
                              help="Flag incomplete data sets for repeat processing.")
            
            parser.add_option("--fix-complete", action="store_true",
                              dest="fixComplete", default=False,
                              help="Apply a fix to completed local data sets.")
            
            parser.add_option("--clear", action="store_true", dest="clear", default=False,
                              help="Clear files after archiving them.")
            
            parser.add_option("--threads", type="int", dest="numThreads", default=8,
                              help="Number of threads to use.")
            
            parser.add_option("--dry-run", action="store_true", dest="dryRun", default=False,
                              help="Don't touch any data, just display actions to be taken.")
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

        except optparse.OptionError, msg:
            raise Usage(msg)

        if options.pbsName: # Search for all folders this PBS name may apply to
            outputList = findPbsMatch(options.pbsName)
            for i in outputList:
                print i
            return 0 # Don't combine with other commands
            

        if options.checkLocal:
            recordLocalDataStatus(options.dryRun)

        if options.archive:
            archiveAllCompletedResults(options.clear, options.dryRun)
        
        if options.clearCompressed:
            removeAllCompletedCompressedOutputs(options.dryRun)
        
        if options.flagIncomplete:
            flagIncompleteDataSets(options.dryRun)

        if options.fixComplete:
            fixCompletedDataSets(options.numThreads)
        
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main())
