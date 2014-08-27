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

import IrgIsisFunctions

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
This program generates plots of pc_align output

'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



def readErrorFile(filePath):
    """Reads in a pc_align error file and stores the result in a dict"""

    contents = dict()

# The file looks like this:
# # latitude,longitude,height above datum (meters),error(meters)
# 24.0979890758,-48.8327835539,602.184303282174,317.7317728955579
#...

    latVals   = []
    lonVals   = []
    elevation = []
    errVals   = []

    print 'Reading data file ' + filePath
    errFile = open(filePath, 'r')
    for line in errFile:
        # Skip the header line
        if line[0] == '#':
            continue

        # Otherwise read the elements
        elements = line.strip().split(',')
        latVals.append  (float(elements[0]))
        lonVals.append  (float(elements[1]))
        elevation.append(float(elements[2]))
        errVals.append  (float(elements[3]))

        #print elements

    errFile.close()

    contents['latVals']   = latVals
    contents['lonVals']   = lonVals
    contents['elevation'] = elevation
    contents['errVals']   = errVals

    return contents

def plotError(fileContents, outputPath, skip):
    """Generates a plot of error by location"""
  

    # Scatter plot of the input pointst scaled to error
    sScaling = 100.0 # This is the maximum point size
    #for lon, lat, alt, err in zip(fileContents['latVals'], fileContents['lonVals'], fileContents['elevations'], fileContents['errVals']):
    
    lats   = fileContents['latVals']
    lons   = fileContents['lonVals']
    alts   = fileContents['elevation']
    errors = fileContents['errVals']

    minError = min(errors)
    maxError = max(errors)
    errorRange = maxError - minError

    sizes = []
    for e in errors:
        newSize = (e - minError) / errorRange
        sizes.append(newSize * sScaling)

    # Subsample points for the plot
    lons  = lons [::skip]
    lats  = lats [::skip]
    sizes = sizes[::skip]
    
    #TODO: Use color instead of size to plot error!
    plt.scatter(lons, lats, sizes, c='r', marker='o', label='Size proportional to error score')
    plt.grid(color='gray', linestyle='dashed')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Plot of errors in pc_align output file from ' + str(minError) + ' to ' + str(maxError))
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    print 'Writing output image ' + outputPath
    plt.savefig(outputPath, bbox_extra_artists=[lgd], bbox_inches='tight')
    plt.clf()




#--------------------------------------------------------------------------------

def main():

    val = IrgIsisFunctions.getCubeCenterLatitude('/u/smcmich1/data/lronacProduction/NAC_DTM_M115108088_M115114873/workDir/M115114873LE.correctedMosaic.cub')
    print val
    return 0


    try:
        try:
            usage = "usage: makePcAlignPlots.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-i", "--input", dest="inputPath",
                              help="Path to the input file.")
            parser.add_option("-o", "--output", dest="outputPath",
                              help="Where to write the output file.")
            parser.add_option("-s", "--skip", dest="skip",
                              help="Only plot every N points.")

            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

            if not options.inputPath:  parser.error("Need input file!")
            if not options.outputPath: parser.error("Need output file!")
            
            if not options.skip:  options.skip = 1
  
        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        #filePath   = '/home/smcmich1/data/stereoCorrectionTest/M120168714LE.corrected.cub_stereoCalibrationTemp2/pcAlignOutput-beg_errors.csv'
        #outputPath = '/home/smcmich1/data/stereoCorrectionTest/begErrors.png'

        fileContents = readErrorFile(options.inputPath)
    
        plotError(fileContents, options.outputPath, int(options.skip))

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
