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

sys.path.append('/home/smcmich1/.local/lib/python2.7/site-packages/')
from BeautifulSoup import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

#TODO: Clean this up!
sys.path.append('/home/smcmich1/programs/mechanize-0.2.5/')

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



# Retrieves LOLA data from the WUSTL REST web interface
def retrieveLolaFile(minLat, maxLat, minLon, maxLon, outputFolder):

    # Build a query to the WUSTL REST interface for LOLA data
    lolaUrl        = 'http://oderest.rsl.wustl.edu/test/'
    baseQuery      = '?query=lolardr&results=v&'
    locationParams = 'maxlat='+str(maxLat)+'&minlat='+str(minLat)+'&westernlon='+str(minLon)+'&easternlon='+str(maxLon)

    queryUrl = lolaUrl + baseQuery + locationParams

    # Parse the response
    parsedPage = BeautifulSoup(urllib2.urlopen((queryUrl)).read())
    
    # Find the link containing '_pts_csv.csv' and download it
    found = False
    for url in parsedPage.findAll('url'):
        if (url.string.find('_pts_csv.csv') >= 0):
            os.system("wget -P " + outputFolder + "  " + url.string)
            found = True
            
    return found


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
	response = br.submit()
	parsedResponse = BeautifulSoup(response.read())


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
	return results


# Obtains the full list of files required to replicate ASU's DEMs from their webpage
def getDataList():
    baseUrl     = "http://wms.lroc.asu.edu/lroc/dtm_select?page="
    currentPage = 1

    # Get URL to current page
    currentPageUrl = baseUrl + str(currentPage)

    # Parse the current page
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((currentPageUrl)).read())


    # Figure out how many pages in total
    largestPage = 1
    pageNavNode = parsedIndexPage.find(id="dtm_select_pagenav")
    for line in pageNavNode.findAll('a'):
        index = line.get('href').find("page=")
        if ( index >= 0 ):  
            if (line.string.find("Next") < 0):
                page = int(line.string)
                if (page > largestPage):
                    largestPage = page
    print "Found " + str(largestPage) + " pages of DEMs"

	# Loop through all index pages and collect DEM pages
    dtmPageList = []
    for currentPage in range(1,largestPage+1):
        currentPageUrl  = baseUrl + str(currentPage)
        parsedIndexPage = BeautifulSoup(urllib2.urlopen((currentPageUrl)).read())
        tableNode       = parsedIndexPage.find(id="dtm_select_selectiontable")
        for line in tableNode.findAll('a'):	
            if (line.get('href').find("dtm_detail")	> 0):
                dtmPageList.append('http://wms.lroc.asu.edu/' + line.get('href'))
		
    print "Found " + str(len(dtmPageList)) + " DEM pages"

    outputFilePath = 'logFile.txt'
    outputFile = open(outputFilePath, 'w')

	# Loop through all individual pages and get download links
    for p in dtmPageList:

#		print p
        thisPage        = BeautifulSoup(urllib2.urlopen((p)).read())
        downloadSection = thisPage.find(id="dtm_downloads")	

        # Find the two input files (the links are not here but we can get the names)
        firstImgFile  = "NOT_FOUND"
        secondImgFile = "NOT_FOUND"
        imgFileRegex  = "_M[0-9]*_[a-zA-Z0-9]*.IMG$"
        for line in downloadSection.findAll('a'):

            matchObj = re.search( imgFileRegex, line.get('href'))
            if matchObj:
                    startIndex = line.string.rfind("_M")  + 1
                    stopIndex  = line.string.find("_", startIndex) 
                    imgFile    = line.string[startIndex:stopIndex]
                    if (firstImgFile  == "NOT_FOUND"):
                        firstImgFile = imgFile
                    else:
                        if imgFile != firstImgFile:
                            secondImgFile = imgFile
                            break

		# Find the ASU DEM	
        demLink = "NOT_FOUND"
        for line in downloadSection.findAll('a'):
            if  line.get('href').find(".TIF") >= 0:
                demLink = line.get('href')
                break

        # Get the lat/lon boundaries
        positionNode = thisPage.find(id="detailtable")
        rows         = positionNode.findAll('tr')
   		
        cols   = rows[0].findAll('td')
        minLat = cols[1].string
        maxLat = cols[3].string
   		
        cols   = rows[1].findAll('td')
        minLon = cols[1].string
        maxLon = cols[3].string

		# Track down the links to the input files
        if (firstImgFile != "NOT_FOUND"):
        	print 'Finding links for image ' + firstImgFile
        	firstImgDownloadPaths  = getLinksForImgFile(firstImgFile)
        if (secondImgFile != "NOT_FOUND"):
        	print 'Finding links for image ' + secondImgFile
        	secondImgDownloadPaths = getLinksForImgFile(secondImgFile)

        # Log the results
        outputFile.write('----------------------------------------------\n')
        outputFile.write('ASU DTM = ' + demLink + '\n')
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
def retrieveDataFiles(logPath, outputDir):
	
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    currentOutputFolder = outputDir
    minLat = 0
    maxLat = 0
    minLon = 0
    maxLon = 0
    for line in open(logPath, 'r'):
        if (line.find('ASU') >= 0):
            # Create new folder for this DEM
            eqPos       = line.find('=')
            asuUrl	    = line[eqPos+2:]
            asuFileName = os.path.basename(asuUrl)
            noExt       = os.path.splitext(asuFileName)[0]
            demName     = noExt[8:]
            currentOutputFolder = outputDir           + demName + '/'
            asuDemPath          = currentOutputFolder + asuFileName
            if not os.path.exists(currentOutputFolder):
                os.makedirs(currentOutputFolder)

            # wget DEM
            #			print "wget --directory-prefix=" + currentOutputFolder + "  " + asuUrl
            if not os.path.exists(asuDemPath.strip()):
                os.system("wget -P " + currentOutputFolder + "  " + asuUrl)
 
        elif (line.find('.IMG') >= 0):
            # wget image
            imgUrl	    = line
            imgFileName = os.path.basename(imgUrl)
            imgCopyPath = currentOutputFolder + '/' + imgFileName
            #			print "wget -P " + currentOutputFolder + "  " + imgUrl
            print imgCopyPath
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
            TODO: Verify this!
            retrieveLolaFile(minLat, maxLat, minLon, maxLon, currentOutputFolder):




#--------------------------------------------------------------------------------

#TODO: Support for file based logging of results

def main():

    print "Started lronacPipeline.py"

    try:
        try:
            usage = "usage: lronacPipeline.py [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.set_defaults(lowMem =False)
#            parser.set_defaults(threads=4)
            parser.add_option("-i", "--input-folder", dest="inputFolder",
                              help="Specifies the folder to operate on.")
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

#            if not args: parser.error("need .IMG files")

        except optparse.OptionError, msg:
            raise Usage(msg)

#TODO: Verify input folder is present!

        print "Beginning processing....."

        startTime = time.time()

        retrieveLolaFile(12.0, 12.1, 10.0, 10.1, '~/repot/lronacPipeline')


#        getDataList()
	
        # Download all of the data we need 
        print 'Retrieving data files'
        #retrieveDataFiles('logFile.txt', options.inputFolder)

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
