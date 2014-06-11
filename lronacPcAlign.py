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

import IrgFileFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Call pc_align in a robust manner for matching LRONAC DTMs to LOLA points.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------


def robustLronacPcAlignCall(demPointsPath, lolaPath, outputPrefix):

    # Use pc-align to compare points to LOLA DEM, compute rotation and offset
    
    # The max-displacement threshold will be adjusted until we are either using a certain number of LOLA points
    #  or until we are using a certain percentage of the input points.
    MIN_NUM_LOLA_POINTS    = 4000
    MIN_LOLA_PERCENTAGE    = 0.005
    STARTING_DISPLACEMENT  = 50    # --max-displacement starts at this value..
    DISPLACEMENT_INCREMENT = 50    #                    Increments by this value...
    MAX_MAX_DISPLACEMENT   = 1800  #                    And we quit when it hits this value.


    # Determine the number of points we want
    numLolaPoints          = IrgFileFunctions.getFileLineCount(lolaPath) - 1
    currentMaxDisplacement = STARTING_DISPLACEMENT
    minNumPointsToUse      = min(MIN_NUM_LOLA_POINTS, int(MIN_LOLA_PERCENTAGE*float(numLolaPoints)))
    
    print 'Using pc_align to compute a transform between intersected points and LOLA points...'
    print 'Starting pc_align, looking for ' + str(minNumPointsToUse) + ' lola point matches.' 
   
    transformPath  = outputPrefix + '-inverse-transform.txt'
    endErrorPath   = outputPrefix + '-end_errors.csv'
        
    while(True):
        cmd = ('pc_align --highest-accuracy --max-displacement ' + str(currentMaxDisplacement) + ' --datum D_MOON ' + 
               '--save-inv-transformed-reference-points ' + demPointsPath + 
               ' ' + lolaPath + ' -o ' + outputPrefix + ' --compute-translation-only')
        print cmd
        os.system(cmd)
        
        if not os.path.exists(endErrorPath):
            numLolaPointsUsed = 0 # Failed to produce any output, maybe raising the error cutoff will help?
        else:
            numLolaPointsUsed = IrgFileFunctions.getFileLineCount(endErrorPath) - 1
    
        if (numLolaPointsUsed >= minNumPointsToUse):
            break # Success!
        elif (currentMaxDisplacement >= MAX_MAX_DISPLACEMENT): # Hit the maximum max limit!
            raise Exception('Error! Unable to find a good value for max-displacement in pc_align.  Wanted '
                            + str(minNumPointsToUse) + ' points, only found ' + str(numLolaPointsUsed))
        else: # Try again with a higher max limit
            print ('Trying pc_align again, only got ' + str(numLolaPointsUsed)
                   + ' lola point matches with value ' + str(currentMaxDisplacement))
            currentMaxDisplacement = currentMaxDisplacement + DISPLACEMENT_INCREMENT            


    if not os.path.exists(transformPath):
        raise Exception('pc_align call failed!')

    return True    


def main(argsIn):

    print "Started lronacPcAlign.py"

    try:
        try:
            usage = "usage: lronacPcAlign.py <lronac DTM> <lola data> <output prefix> [--manual]\n  "
            parser = optparse.OptionParser(usage=usage)

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args(argsIn)

            if len(args) < 3:
                raise Exception('Not enough arguments provided.')
              
            options.demPath      = args[0]
            options.lolaPath     = args[1]
            options.outputPrefix = args[2]

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        robustLronacPcAlignCall(options.demPath, options.lolaPath, options.outputPrefix)
      

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
