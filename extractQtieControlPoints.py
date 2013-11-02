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

job_pool = [];

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Applies the LROC offset from spacecraft position to an LROC cube's spice data
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


#--------------------------------------------------------------------------------


def main():

    print "Started qtie point extractor"

    try:
        try:
            usage = "usage: extractQtieControlPoints.py [--output <path>][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--cnetPath",  dest="cnetPath",  help="Path to qtie .net file")

            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output text file.")
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--keep", action="store_true", dest="keep",
                              help="Do not delete the temporary files.")
            (options, args) = parser.parse_args()

            if not options.cnetPath: 
                parser.error("Need cnet path")
            if not options.outputPath: 
                parser.error("Need output path")


        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()


        # Convert from binary control net file to PVL file
        workingFolder = os.path.dirname(options.outputPath)
        inputFile     = os.path.basename(options.cnetPath)
        pvlPath       = os.path.join(workingFolder, inputFile) + '.pvl'
        cmd = "cnetbin2pvl from=" + options.cnetPath + " to=" + pvlPath;
        print cmd
        os.system(cmd)

        # Make sure the files are ready

        # If the file already exists, delete it and rewrite it.
        if os.path.exists(options.outputPath):
            os.remove(options.outputPath)
        if not os.path.exists(pvlPath):
            print 'Error, file ' + pvlPath + ' not found!'
            return 0

        # Find and modify this line: TKFRAME_-85610_ANGLES    = ( 0.0, 0.0, 0.0 )

        # Now convert the pvl file into a convenient csv format
        # - Each line in the output file is one point: sample1, line1, sample2, line2,

        print 'reading file ' + pvlPath
        i = open(pvlPath,    'r')
        f = open(options.outputPath, 'w')
        for line in i:
            #print 'input --> ' + line
            if (line.find('End_Object') >= 0): # Done with info for a point, move to next line
                f.write('\n')
            elif ( (line.find(' Sample ') >= 0) or (line.find(' Line ') >= 0) ): # Copy point info
                pos   = line.find('=')
                value = line[pos+2:-1]
                f.write(value + ', ')
        i.close()
        f.close()

        # Clean up temporary files
        if not options.keep:
            os.remove(pvlPath)


        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main())
