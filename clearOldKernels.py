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

import sys, os, glob, optparse, re, shutil, subprocess, string, time, logging, threading, math


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generate a DERT landscape directory from orthoimage and DEM inputs
'''
    sys.exit()


# Sample kernel name:
# moc42_2013029_2013030_v02.bc

def isNewerKernel(k1, k2):
    """Returns True if the second file is a newer kernel than the first"""
    
    # All date numbers must match and the version must be higher   
    return (k1[6:21] == k2[6:21]) and (k2[23:25] > k1[23:25])

def clearOldKernels(kernelFolder, dryRun=False):
    
    KERNEL_PREFIX = 'moc42_20'
    MAX_LOOKAHEAD = 4
    
    filesToDelete = []
    
    # Loop through each file in the directory
    print 'Finding files to delete...'
    fileList = os.listdir(kernelFolder)
    fileList.sort()
    for i, f in enumerate(fileList):
        #print f
        if not f.startswith(KERNEL_PREFIX):  # Skip entries that are not what we are looking for.
            continue
        
        # Look at the next several elements and see if any are newer versions of this kernel.
        try:
            for c in range(i+1, i+MAX_LOOKAHEAD+1):
                #print ' - ' + fileList[c]
                if isNewerKernel(f, fileList[c]): # Check if this file is newer.
                    filesToDelete.append(i)       # If so, add the current file to the delete list
                    break                         #  and break out of this loop.
        except:
            pass # Use this try-catch code to take care of falling off the end of the list
        
    # Now go through and delete all the files
    print 'Deleting files...'
    for i in filesToDelete:
        fullPath = os.path.join(kernelFolder, fileList[i])
        if dryRun: # Just print the files to be deleted
            print fullPath
        else: # Delete the file
            os.remove(fullPath)
        
    return len(filesToDelete)



# Set up a pair of files to be loaded into DERT
def main(argsIn):
    
    try:
        usage = "usage: clearOldKernels.py [--help][--manual]\n  "
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("--manual", action="callback", callback=man,
                          help="Read the manual.")
        
        parser.add_option("--inputFolder", dest="inputFolder", help="Folder to clear out kernels in")
        
        parser.add_option("--dryRun", action="store_true", dest="dryRun",
                          help="Just print the names of the files that would be deleted.")
            
        (options, args) = parser.parse_args()

        if not options.inputFolder:
            parser.error('Missing required inputs!')

    except optparse.OptionError, msg:
        raise Usage(msg)

    print 'Operating on folder ' + options.inputFolder
    
    clearOldKernels(options.inputFolder, options.dryRun)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    