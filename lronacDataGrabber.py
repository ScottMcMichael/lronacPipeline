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

#sys.path.append('/home/smcmich1/.local/lib/python2.7/site-packages/')
from bs4 import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

#TODO: Clean this up!
#sys.path.append('/home/smcmich1/programs/mechanize-0.2.5/')

import mechanize

job_pool = [];

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing file for LRONAC batch processing '''

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


# Gets the download links to the LE and RE parts of a given LRONAC ID
def getLinksForImgFile(productId):

  pdsUrl = 'http://wms.lroc.asu.edu/lroc/search'

  # Open browser object to ASU data search page 
  br = mechanize.Browser()
  br.open(pdsUrl)

  # Get unnamed form handle, set product ID filter
  br.form = list(br.forms())[0]
  control = br.form.find_control("filter[product_id]")
  control.value = productId #'M112646261'

  # Submit the form, then parse the response of the form submission
  response       = br.submit()
  parsedResponse = BeautifulSoup(response.read()); 

  # Get the links to the LE and RE file pages
  resultsNode = parsedResponse.find(id="table")

  leftRegex  = "M[0-9]*LE$"
  rightRegex = "M[0-9]*RE$"
  leftLink  = 'NOT_FOUND'
  rightLink = 'NOT_FOUND'
  for line in resultsNode.findAll('a'):
    if re.search( leftRegex, line.get('href')):
      leftLink = 'http://wms.lroc.asu.edu' + line.get('href')
    if re.search( rightRegex, line.get('href')):
      rightLink = 'http://wms.lroc.asu.edu' + line.get('href')

  # Extract the left and right EDR paths
  leftEdrPath  = 'NOT_FOUND'
  rightEdrPath = 'NOT_FOUND'
  if (leftLink != 'NOT_FOUND'):
    leftPage  = BeautifulSoup(urllib2.urlopen(leftLink).read())
    for link in leftPage.findAll('a'):
      if link.string == 'Download EDR':
        leftEdrPath = link.get('href')
  if (rightLink != 'NOT_FOUND'):
    rightPage = BeautifulSoup(urllib2.urlopen(rightLink).read())
    for link in rightPage.findAll('a'):
      if link.string == 'Download EDR':
        rightEdrPath = link.get('href')

  # Return the output
  results = (leftEdrPath, rightEdrPath)
#  print results
  return results

# Retrieves LOLA data from the WUSTL REST web interface
def retrieveLolaFile(minLat, maxLat, minLon, maxLon, outputFolder):

    # We add this much padding on each side of the LOLA data to give the DEM some room to shift
    LOLA_EXTRA_SPACE = 0.25

    # Build a query to the WUSTL REST interface for LOLA data
    lolaUrl        = 'http://oderest.rsl.wustl.edu/test/'
    baseQuery      = '?query=lolardr&results=v&'
    locationParams = ('maxlat='      + str(maxLat+LOLA_EXTRA_SPACE) +
                      '&minlat='     + str(minLat-LOLA_EXTRA_SPACE) +
                      '&westernlon=' + str(minLon-LOLA_EXTRA_SPACE) +
                      '&easternlon=' + str(maxLon+LOLA_EXTRA_SPACE))

    outputPath = os.path.join(outputFolder, 'lolaRdrPoints.csv')
    if os.path.exists(outputPath):
        return True

    queryUrl = lolaUrl + baseQuery + locationParams
    print queryUrl
    
    # Parse the response
    parsedPage = BeautifulSoup(urllib2.urlopen((queryUrl)).read())
    #print parsedPage.prettify()

    # Find the link containing '_pts_csv.csv' and download it
    found = False
        
    for url in parsedPage.findAll('url'):
        if (url.string.find('_pts_csv.csv') >= 0):
            os.system("wget --output-document=" + outputPath + "  " + url.string)
            found = True
            
    return found


# Obtains the full list of files required to replicate ASU's DEMs from their webpage
def getDataList(outputFilePath):
    baseUrl     = "http://wms.lroc.asu.edu/lroc/rdr_product_select?page="
    currentPage = 1

    # Get URL to current page
    currentPageUrl = baseUrl + str(currentPage)

    # Parse the current page
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((currentPageUrl)).read())

    #print parsedIndexPage.prettify()

    # Figure out how many pages in total
    largestPage = 1
    #pageNavNode = parsedIndexPage.find(id="dtm_select_pagenav")
    for line in parsedIndexPage.findAll('a'):
        index = line.get('href').find("page=")
        if ( index >= 0 ):  
            if (line.string.find("Next") < 0):
                page = int(line.string)
                if (page > largestPage):
                    largestPage = page
    print "Found " + str(largestPage) + " pages of DEMs"

    # Loop through all index pages and collect DEM pages
    dtmPageList = []
    for currentPage in range(1, largestPage+1):
        currentPageUrl  = baseUrl + str(currentPage) + '&sort=time_reverse'
        parsedIndexPage = BeautifulSoup(urllib2.urlopen(currentPageUrl).read())
        #tableNode       = parsedIndexPage.find(id="dtm_select_selectiontable")
        for line in parsedIndexPage.findAll('a'):	
            if (line.get('href').find("view_rdr")	> 0):
                dtmPageList.append('http://wms.lroc.asu.edu' + line.get('href'))

    print "Found " + str(len(dtmPageList)) + " DEM pages"

    outputFile = open(outputFilePath, 'w')

    # Loop through all individual pages and get download links
    for p in dtmPageList:

        print p
        thisPage = BeautifulSoup(urllib2.urlopen(p).read())

        downloadSections = thisPage.findAll(attrs={"class": "download_container"})
        if len(downloadSections) > 1:
            downloadSection  = downloadSections[1] # Want the second of two instances of this
        else:
            downloadSection = thisPage # Probably going to fail later

        # Find the two input files (the links are not here but we can get the names)
        firstImgFile  = "NOT_FOUND"
        secondImgFile = "NOT_FOUND"
        #imgFileRegex  = "_M[0-9]*_[a-zA-Z0-9]*$"
        imgFileRegex  = "_M[0-9]*_[0-9C]*M$" 
        for line in downloadSection.findAll('a'):
            matchObj = re.search( imgFileRegex, line.get('href'))
            if matchObj:
                #startIndex = line.string.rfind("_M")  + 1
                #stopIndex  = line.string.find("_", startIndex) 
                #imgFile    = line.string[startIndex:stopIndex]
                if (firstImgFile  == "NOT_FOUND"):
                    firstImgFile = line.get('href')
                else:
                    if line.get('href') != firstImgFile:
                        secondImgFile = line.get('href')
                        break

#        print firstImgFile
#        print secondImgFile

	# Find the ASU DEM	
        demLink = "NOT_FOUND"
        for line in downloadSection.findAll('a'):
            #if  line.get('href').find(".TIF") >= 0:
            if line.string and line.string.find("(32-bit GeoTIFF)") >= 0:
                demLink = line.get('href')
                break

        #print demLink
        if (demLink == 'NOT_FOUND'):
            print 'Failed to find link to DEM page'

        # Init to flag values
        minLat = -999
        maxLat = -999
        minLon = -999
        maxLon = -999

        if (demLink!= "NOT_FOUND"): # Get real DEM link and boundaries
            url     = 'http://wms.lroc.asu.edu' + demLink
            demPage = BeautifulSoup(urllib2.urlopen(url).read())

            # Get the full download link
            pos          = demLink.rfind('/')
            demName      = demLink[pos+1:]
            demLinkFound = False
            searchText   = demName+'.TIF'
            print searchText
            for line in demPage.findAll('a'):
                if line.get('href').find(searchText) > 0:
                    demLink = line.get('href')
                    demLinkFound = True
                    break
            if not demLinkFound: # Note if we failed to find the entire DEM link
                for line in demPage.findAll('a'):
                    print line
                demLink = 'NOT_FOUND'

            # Get the lat/lon boundaries
            # - The exact position within the HTML is semi-random so we have to check the text
            tableTop     = demPage.find(attrs={"class": "presentable_data"})           
            positionNode = tableTop.contents[1]   
            rows         = positionNode.findAll('tr')
         
            #print tableTop.prettify()
   
            for i in range(0,3): # For each of three rows of interest
                cols = rows[i].findAll('td')
                for j in range(0, 3, 2): # Hit indices 0 and 2
                    if   (cols[j].string.find('Minimum Latitude') >= 0):
                        minLat = cols[j+1].string
                    elif (cols[j].string.find('Maximum Latitude') >= 0):
                        maxLat = cols[j+1].string
                    elif (cols[j].string.find('Western-most Longitude') >= 0):
                        minLon = cols[j+1].string
                    elif (cols[j].string.find('Eastern-most Longitude') >= 0):
                        maxLon = cols[j+1].string

        else:
            demLink = p # Record the full URL so it is easier to check why we failed to find the DEM

        # Track down the links to the input files
        if (firstImgFile != "NOT_FOUND"):
          startIndex = firstImgFile.rfind("_M")  + 1
          stopIndex  = firstImgFile.find("_", startIndex) 
          productId  = firstImgFile[startIndex:stopIndex]          
          print 'Finding links for image ' + productId
          firstImgDownloadPaths  = getLinksForImgFile(productId)
        if (secondImgFile != "NOT_FOUND"):
          startIndex = secondImgFile.rfind("_M")  + 1
          stopIndex  = secondImgFile.find("_", startIndex) 
          productId  = secondImgFile[startIndex:stopIndex]          
          print 'Finding links for image ' + productId
          secondImgDownloadPaths  = getLinksForImgFile(productId)


        # Log the results
        outputFile.write('----------------------------------------------\n')
        if (demLink != "NOT_FOUND"):
           outputFile.write('ASU DTM = ' + demLink + '\n')
        if ( (minLat != -999) and (maxLat != -999) and (minLon != -999) and (maxLon != -999) ):
            outputFile.write('Min lat = ' + minLat + '\n')
            outputFile.write('Max lat = ' + maxLat + '\n')
            outputFile.write('Min lon = ' + minLon + '\n')
            outputFile.write('Max lon = ' + maxLon + '\n')
        if (firstImgFile != "NOT_FOUND"):
            outputFile.write(firstImgDownloadPaths[0] + '\n')
            outputFile.write(firstImgDownloadPaths[1] + '\n')
        if (secondImgFile != "NOT_FOUND"):
            outputFile.write(secondImgDownloadPaths[0] + '\n')
            outputFile.write(secondImgDownloadPaths[1] + '\n')

    # Done obtaining data, close the file.
    outputFile.close()

# Fetches all of the files listed in the log folder and puts them in different directories
def retrieveDataFiles(logPath, outputDir, name=''):
	
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    currentOutputFolder = outputDir
    minLat = 0
    maxLat = 0
    minLon = 0
    maxLon = 0
    nameBlock = False
    for line in open(logPath, 'r'):
        if (line.find('ASU') >= 0):
            nameBlock = False # Stop blocking on a name when we get to a new one
            
            # Create new folder for this DEM
            eqPos       = line.find('=')
            asuUrl	    = line[eqPos+2:]
            asuFileName = os.path.basename(asuUrl)
            noExt       = os.path.splitext(asuFileName)[0]
            demName     = noExt[8:]
            if (name != '') and (demName != name): # Not the DEM we are looking for   
                nameBlock = True # Ignore data until we hit the next name
                continue         # Keep looking through the file
                
            currentOutputFolder = outputDir           + demName + '/'
            asuDemPath          = currentOutputFolder + asuFileName
            if not os.path.exists(currentOutputFolder):
                os.makedirs(currentOutputFolder)

            # wget DEM
            #			print "wget --directory-prefix=" + currentOutputFolder + "  " + asuUrl
            if not os.path.exists(asuDemPath.strip()):
                os.system("wget -P " + currentOutputFolder + "  " + asuUrl)
        elif nameBlock:
            pass # Do nothing until the block is cleared
        elif (line.find('.IMG') >= 0):
            # wget image if we don't already have it
            imgUrl	    = line
            imgFileName = os.path.basename(imgUrl).strip() # Need to clear whitespace
            imgCopyPath = os.path.join(currentOutputFolder + imgFileName)
            if not os.path.exists(imgCopyPath):
		os.system("wget -P " + currentOutputFolder + "  " + imgUrl)

        # Read bounding box
        elif (line.find('Min lat') >= 0):
            minLat = float(line[10:])
        elif (line.find('Max lat') >= 0):
            maxLat = float(line[10:])
        elif (line.find('Min lon') >= 0):
            minLon = float(line[10:])
        elif (line.find('Max lon') >= 0): # The last BB entry, download the LOLA data file
            maxLon = float(line[10:])
            
            retrieveLolaFile(minLat, maxLat, minLon, maxLon, currentOutputFolder)


#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    print "Started lronacDataGrabber.py"

    try:
        try:
            usage = "usage: lronacDataGrabber.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-f", "--fetch", action="store_true",
                              dest="fetch",
                              help="Fetch the list of data locations (slow).")
            parser.add_option("-i", "--input-file", dest="inputFile",
                              help="Specifies the data list file to read from.",
                              default='lronacDataSourceList.txt')
            parser.add_option("-o", "--output-folder", dest="outputFolder",
                              help="Specifies the folder to copy the data to.",
                              default='./')
            parser.add_option("-n", "--name", dest="name",
                              help="Only get the data for the DTM with this name.",
                              default='')
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

            if (not options.fetch) and (not options.inputFile):
                print 'Error: Need to either fetch or specify an input file!'
                return 1

        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

#        retrieveLolaFile(12.0, 12.1, 10.0, 10.1, '~/repot/lronacPipeline')

        # If option was passed in, build the list of data locations.
        if options.fetch:
            print 'Searching for data sources'
            getDataList(options.inputFile)
            
        else: # Download all of the data we need 
            if options.name:
                print 'Retrieving data for location ' + options.name
            else:
                print 'Retrieving ALL data files'
            retrieveDataFiles(options.inputFile, options.outputFolder, options.name)

        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
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
