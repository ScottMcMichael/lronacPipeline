#!/usr/bin/env python

import os, glob, optparse, re, shutil, subprocess, sys, string, time

import simplekml



def parseBoundingBoxes(folder):
    '''Gets the bounding box of one of our files'''
        
    downloadLogPath = os.path.join(folder, 'downloadLog.txt')

    # Search the download log for the bounding box    
    f = open(downloadLogPath, 'r')
    imgFiles = []
    minLon = None
    minLat = None
    maxLon = None
    maxLat = None
    for line in f:
        if 'minLon' in line:
            numStart = line.find(':') + 1
            minLon   = line[numStart:]
        if 'maxLon' in line:
            numStart = line.find(':') + 1
            maxLon   = line[numStart:]
        if 'minLat' in line:
            numStart = line.find(':') + 1
            minLat   = line[numStart:]
        if 'maxLat' in line:
            numStart = line.find(':') + 1
            maxLat   = line[numStart:]

    f.close()

    # Make sure we got all four values
    if (not minLon) or (not minLat) or (not maxLon) or (not maxLat):
        raise Exception('Could not find BB in folder ' + folder)

    return (minLon, minLat, maxLon, maxLat)


def generateKml(outputPath, bbDict, colors):
    '''Plots all of the bounding boxes with the specified colors'''
    
    # Initialize kml document
    kml = simplekml.Kml()
    kml.document.name = 'DTM Locations'
    kml.hint = 'target=moon'
    
    #minError = min(asuMeanList)
    #maxError = max(asuMeanList)
    #errorRange = maxError - minError
    
    # Try to set up bounding box for each used folder
    i = 0
    for k, v in bbDict:
        (minLon, minLat, maxLon, maxLat) = v
        
        poly = kml.newpolygon(name=f, outerboundaryis=[(maxLat,minLon), (maxLat,maxLon), (minLat, maxLon), (minLat,minLon), (maxLat,minLon)])
#        poly.style.polystyle.color   = simplekml.Color.red # Why is polygon not working?
#        poly.style.polystyle.fill    = 1
#        poly.style.polystyle.outline = 1
        poly.style.linestyle.width   = 5
        
        # For now just draw white lines
        if colors[i] == 'red':
            poly.style.linestyle.color   = simplekml.Color.rgb(255,0,0,255)
        else: # Default is white
            poly.style.linestyle.color   = simplekml.Color.rgb(255,255,255,255)
        
        ## Generate a color based on the error value:  white (low error) <--> red (high error)
        #colorVal = 255 - 255*(asuMeanList[i] - minError)/errorRange
        #poly.style.linestyle.color   = simplekml.Color.rgb(255,colorVal,colorVal,255)
        
        i = i + 1
    
    # Save kml document
    kml.save(outputPath)


def getInputImages(folder):
    '''Find out the input images used in a given folder'''

    downloadLogPath = os.path.join(folder, 'downloadLog.txt')
    
    f = open(downloadLogPath, 'r')
    imgFiles = []
    for line in f:
        if '.IMG' in line:
            nameStart = line.rfind('/') + 1
            imgName = line[nameStart:]
            imgFiles.append(imgName)
    f.close()

    if len(imgFiles) != 4:
        raise Exception('Failed to find four IMG files in folder ' + folder)

    return imgFiles

def getAsuImageList():
    '''Returns dictionary containing input images for each ASU data set'''
    
    sourceFile = '/u/smcmich1/projects/lronacPipeline/logFile.txt'
    
    asuDict = {}
    setName = None
    f = open(sourceFile, 'r')
    for line in f:
        
        if 'ASU DTM' in line:
            
            # Grab the name from the line
            nameStart = line.rfind('/') + 1
            setName = line[nameStart:]
            asuDict[setName] = []
        
        if '.IMG' in line:
            # Get the name and append to the list for this set
            nameStart = line.rfind('/') + 1
            imgName = line[nameStart:]
            asuDict[setName].append(imgName)
    f.close()
    
    print asuDict
    
    return(asuDict)
    

def findAsuMatches():
    
    # Location of processed data
    dataFolder = '/u/smcmich1/data/lronacProduction'
    bbKmlPath  = '/u/smcmich1/data/lronacProduction/boundingBoxes.kml'
    
    print 'Building ASU dictionary...'
    asuImageList = getAsuImageList()
    
    print 'Looping through pipeline output...'
    
    # Loop through all folders in the directory
    foldersInDirectory = os.listdir(dataFolder)
    foldersInDirectory.sort()
    bbDict = {}
    colors = []
    i = 0
    for f in foldersInDirectory:

        folderPath = os.path.join(dataFolder, f)

        # Get the bounding box from the folder
        bb = parseBoundingBoxes(folderPath)
        bbDict[f] = bb

        # Get the list of input images used by that folder
        imageList = getInputImages(folderPath)
        
        # Search for those images in the ASU list
        allContained = True
        for k, v in asuImageList: # Loop through all ASU data sets
        
            for image in imageList: # Loop through all image for current folder
                if image not in v:
                    allContained = False
                    break
                
            if allContained:
                print 'Folder ' + f + ' matched ASU data set ' + k
                break
            
        # Handle unmatched DTMs and assign color for KML plot
        if not allContained:
            colors.append('white')
            print 'Folder ' + f + ' did not match any data set'
        else:
            colors.append('red')
            
        i = i + 1
    
    # Create a kml file containing the bounding box of each DTM    
    generateKml(bbKmlPath, bbDict, colors)
        
    return 0
    
    


def main():
    
    findAsuMatches()
    
    pass



if __name__ == "__main__":
    sys.exit(main())