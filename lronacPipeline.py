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


import os, glob, optparse, re, shutil, subprocess, sys, string, time, urllib2
from BeautifulSoup import BeautifulSoup

#TODO: Clean this up!
sys.path.append('/home/smcmich1/mechanize-0.2.5/')

import mechanize

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


def getLinksForImgFile(productId):

	pdsUrl = 'http://wms.lroc.asu.edu/lroc/search'

	# Open browser object to ASU data search page 
	br = mechanize.Browser()
	response = br.open(pdsUrl)
 
	# Get unnamed form handle, set product ID filter
	br.form = list(br.forms())[0]
	control = br.form.find_control("filter[product_id]")
	control.value = (productId #'M112646261'

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

	# Parse the left and right side result pages
	leftPage  = BeautifulSoup(urllib2.urlopen(leftLink).read())
	rightPage = BeautifulSoup(urllib2.urlopen(rightLink).read())

	# Extract the left and right EDR paths
	leftEdrPath = 'NOT_FOUND'
	for link in leftPage.findAll('a'):
		if link.string == 'Download EDR':
			leftEdrPath = link.get('href')
#			print 'left = ' + leftEdrPath

	rightEdrPath = 'NOT_FOUND'
	for link in rightPage.findAll('a'):
		if link.string == 'Download EDR':
			rightEdrPath = link.get('href')
#			print 'right = ' + rightEdrPath

	# Return the output
	results = (leftEdrPath, rightEdrPath)
	return results



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

		print p
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
				print imgFile
				if (firstImgFile  == "NOT_FOUND"):
					firstImgFile = imgFile
					print "first!"
				else:
					if imgFile != firstImgFile:
						secondImgFile = imgFile
						print "second!"
						break

		# Find the ASU DEM	
		demLink = "NOT_FOUND"
		for line in downloadSection.findAll('a'):
			if  line.get('href').find(".TIF") >= 0:
				demLink = line.get('href')
				break

		outputFile.write('----------------------------------------------\n')
		outputFile.write('First  Input = ' + firstImgFile + '\n')
		outputFile.write('Second Input = ' + secondImgFile + '\n')
		outputFile.write('ASU DTM      = ' + demLink + '\n')

	outputFile.close()



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



		

#TODO: Add GIT login name/email


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
