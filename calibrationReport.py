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
Tool for converting lat, lon, alt .csv files into a KML visualization of the points.
'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

#--------------------------------------------------------------------------------


def plotPoints(pointList, kml, color, size, prefix, lineSkip):
    
    numPoints = len(pointList) / 3
    
#    minError = min(asuMeanList)
#    maxError = max(asuMeanList)
#    errorRange = maxError - minError
    
    style = simplekml.Style()
    if color=='blue':
        style.labelstyle.color = simplekml.Color.blue
    elif color=='red':
        style.labelstyle.color = simplekml.Color.red
    elif color=='green':
        style.labelstyle.color = simplekml.Color.green
    elif color=='yellow':
        style.labelstyle.color = simplekml.Color.yellow
    else: 
        style.labelstyle.color = simplekml.Color.white

    if size=='small':
        style.labelstyle.scale    = 0
        style.iconstyle.scale     = 0.7
        style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/open-diamond.png'
        style.iconstyle.color = style.labelstyle.color
    elif size=='tiny':
        style.labelstyle.scale    = 0
        style.iconstyle.scale     = 0.5
        style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
        style.iconstyle.color = style.labelstyle.color
    # All other words trigger the default

    # Plot each point
    counter = 0
    for i in range (0,numPoints, lineSkip):
    
        point = kml.newpoint(name=prefix+str(i), coords=[(pointList[i*3], pointList[i*3+1], pointList[i*3+2])], \
                              gxaltitudemode= simplekml.AltitudeMode.absolute)
        point.style   = style
        point.extrude = 1
        counter = counter + 1
        
#        # Generate a color based on the error value:  white (low error) <--> red (high error)
#        colorVal = 255 - 255*(asuMeanList[i] - minError)/errorRange
#        poly.style.linestyle.color   = simplekml.Color.rgb(255,colorVal,colorVal,255)
 
    print 'Added ' + str(counter) + ' KML points'

    return kml

#--------------------------------------------------------------------------------

def main():

    try:
        try:
            usage = "usage: calibrationReport.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--input",  dest="inputPath",  help="Path to GDC point file")
            parser.add_option("--output", dest="outputPath", help="Path of output kml file to write")
            parser.add_option("--name",   dest="name",       help="Name used to identify KML data")
            parser.add_option("--color",  dest="color",      help="Color used to plot points")
            parser.add_option("--size",   dest="size",       help="Size shortcut: (normal / small / tiny)")
            parser.add_option("--skip",   dest="skip",       help="Only sample every N points")
            
            (options, args) = parser.parse_args()

            if not options.inputPath: 
                parser.error("Missing input path")

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        # Set options defaults
        if not options.outputPath:
            options.outputPath = options.inputPath + ".kml"
        if not options.name:
            options.name  = 'gdcToKml'
        if not options.color:
            options.color = 'red'
        if not options.size:
            options.color = 'normal'
        if not options.skip:
            options.skip = 1

        # Read point locations
        points = IsisTools.readPositions(options.inputPath)
        if not points:
            return 0

        # Initialize kml document
        kml = simplekml.Kml()
        kml.document.name = options.name
        kml.hint = 'target=moon'

        kml = plotPoints(points, kml, options.color, options.size, 'I_', int(options.skip))

        # Save kml document
        kml.save(options.outputPath)

        print "Finished calibrationReport"
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
