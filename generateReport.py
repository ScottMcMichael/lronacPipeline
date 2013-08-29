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

            # For each folder, copy each result
            for c in fileList:

		            # Get file paths
                outputFilePath        = os.path.join(folderPath,              c)
                supercomputerFilePath = os.path.join(supercomputerFolderPath, c)

                if not os.path.exists(outputFilePath): # Request file copy
                    grabFile(supercomputerFilePath, outputFilePath)
		            
        
        wait_on_all_jobs() # Wait on all file grab requests
        
    except Exception,e: # Catch any errors, the program will move on to the next folder
        print "Caught: ", e


# Generates a set of plots describing the results
def generatePlots(dataFolder):

    # Search through all the output files and pull out results


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
#        grabResultFiles(dataFolder)
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
