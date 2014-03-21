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

import sys, os, optparse, sys, string, simplekml

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Generates a KML plot of lat/lon/alt/error csv files.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg




def readPositions(positionFilePath):
    """Reads in a list of GDC coordinates with errors"""

    if not os.path.exists(positionFilePath):
        print 'File ' + positionFilePath + ' is missing!'
        return []

    #TODO: Read this from the file?
    MEAN_MOON_RADIUS = 1737400 # from D_MOON

    pointList = []

    f = open(positionFilePath, 'r')
    for line in f:

        if line.find('#') < 0: # Skip lines containing the comment symbol
            strings = line.split(',')
            
            #print strings
            if (len(strings) < 5): # PC_align format
                lat       = float(strings[1])
                lon       = float(strings[0])
                height    = float(strings[2]) - MEAN_MOON_RADIUS
                diff      = float(strings[3])
                thisPoint = (lon, lat, height, diff)
            else: # lolaCompare format
                lat        = float(strings[0])
                lon        = float(strings[1])
                lolaHeight = float(strings[2])
                ourHeight  = float(strings[3])
                diff       = float(strings[4])
                #thisPoint  = (lon, lat, lolaHeight, diff, ourHeight)
                thisPoint  = (lon, lat, ourHeight, diff, lolaHeight)
            pointList.append(thisPoint)
    f.close()

    return pointList




def generateKml(pointList, outputPath, pointSkip, maxErrorLimit, name, color):
    """Generates a KML plot for the input points"""
        
    # Initialize kml document
    kml = simplekml.Kml()
    kml.document.name = name
    kml.hint = 'target=moon'
    
    # Compute the min and max point error
    minError = pointList[0][3]
    maxError = minError
    for p in pointList:
        if p[3] < minError:
            minError = p[3]
        if p[3] > maxError:
            maxError = p[3]
    if (maxErrorLimit > 0): # Apply max error limit if used
        maxError = maxErrorLimit
    errorRange = maxError - minError
    
    print 'min = ' + str(minError)
    print 'max = ' + str(maxError)
    
    # Plot each point
    #style = simplekml.Style()
    counter = 0
    for i in range (0, len(pointList), int(pointSkip)):
    
        lon    = pointList[i][0]
        lat    = pointList[i][1]
        height = pointList[i][2]
        point = kml.newpoint(name=str(counter), coords=[(lon, lat, height)],
                              gxaltitudemode= simplekml.AltitudeMode.absolute)
        
        #print point
        
        #point.style   = style
        point.extrude = 0
        counter       = counter + 1
        
        point.style.labelstyle.scale    = 0
        point.style.iconstyle.scale     = 0.7
        point.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/open-diamond.png'
        
        # Generate a color based on the error value:  white (low error) <--> red (high error)
        thisError = pointList[i][3]
        colorVal = int(255.0 - 255.0*(thisError - minError)/errorRange)
        if (colorVal < 0):
            colorVal = 0
        #print str(pointList[i][3]) + ' --> ' + str(colorVal)
        if color == 'blue':
            point.style.iconstyle.color   = simplekml.Color.rgb(colorVal,colorVal,255,255)
        elif color == 'green':
            point.style.iconstyle.color   = simplekml.Color.rgb(colorVal,255,colorVal,255)
        else: # red
            point.style.iconstyle.color   = simplekml.Color.rgb(255,colorVal,colorVal,255)
    
    
    # Save kml document
    kml.save(outputPath)
    return counter
    


#--------------------------------------------------------------------------------------------
    
def main(argsIn):
    
    try:
        try:
            usage = "usage: pointErrorToKml.py [options] inputPath outputPath\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.set_defaults(skip=1)
            parser.set_defaults(errorLimit=0)
            parser.set_defaults(name='errorPoints')
            parser.set_defaults(color='red')
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("--skip",       dest="skip",       help="Point skip (default none)")
            parser.add_option("--errorLimit", dest="errorLimit", help="Max displayed error")
            parser.add_option("--name",       dest="name",       help="KML name")
            parser.add_option("--color",      dest="color",      help="Set to red, blue, or green")
            (options, args) = parser.parse_args(argsIn)
    
        except optparse.OptionError, msg:
            raise Usage(msg)
    
        if len(args) != 2:
            print usage
            raise Usage('Must pass in input and output paths!')
    
        print "Beginning processing....."
    
        pointList = readPositions(args[0])
        
        generateKml(pointList, args[1], options.skip, float(options.errorLimit), options.name, options.color)
    
    
        print "Finished"
        return 0
    
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

    