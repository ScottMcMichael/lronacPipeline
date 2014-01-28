#!/bin/bash

## This script calls lronacPipeline with the passed in arguments and pipes output to a text file

# Grab the only input argument
INPUT_FOLDER=$1
echo Input folder = $INPUT_FOLDER

## Make sure the environment is correct
##/home1/smcmich1/bashStartup.sh
#source /home1/smcmich1/.profile

# Change to the correct directory
#cd ~/projects/lronacPipeline

#echo "ISISROOT  = $ISISROOT"
#echo "ISISDATA  = $ISISDATA"
#echo "ISIS3DATA = $ISIS3DATA"

# Get the path to the ASU DEM
ASU_FILE=`ls $INPUT_FOLDER/*.TIF`


# Get the paths to the IMG files
IMG_FILES=`ls $INPUT_FOLDER/*.IMG`
set -- $IMG_FILES # Assigns the components to $1, $2, $3, $4
LEFT_IMG=$1
RIGHT_IMG=$2
LEFT_S_IMG=$3
RIGHT_S_IMG=$4

# Get the path to the LOLA file, to make it cleaner
LOLA_FILE_IN=`ls $INPUT_FOLDER/*.csv`
LOLA_FILE_OUT=$INPUT_FOLDER/lolaRdrPoints.csv
if [ ! -f $LOLA_FILE_OUT ]; then
    mv $LOLA_FILE_IN $LOLA_FILE_OUT
fi

# Set output directory options
WORKDIR=$INPUT_FOLDER/workDir
OUTPUT_FOLDER=$INPUT_FOLDER/results

#echo $LEFT_IMG
#echo $RIGHT_IMG
#echo $LEFT_S_IMG
#echo $RIGHT_S_IMG
#echo $ASU_FILE
#echo $LOLA_FILE_OUT

#echo $PATH

# Now call the lronac2dem.py function with the correct inputs
lronac2dem.py --keep --left $LEFT_IMG --right $RIGHT_IMG --stereo-left $LEFT_S_IMG --stereo-right $RIGHT_S_IMG --lola $LOLA_FILE_OUT --asu $ASU_FILE --workDir $WORKDIR --output-folder $OUTPUT_FOLDER

exit 0


