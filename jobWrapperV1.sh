#!/bin/bash

## This script calls lronacPipeline with the passed in arguments and pipes output to a text file

# Grab the only input argument
echo "Input folder = $1"

## Make sure the environment is correct
##/home1/smcmich1/bashStartup.sh
source /home1/smcmich1/.profile

# Change to the correct directory
cd ~/projects/lronacPipeline

echo "ISISROOT  = $ISISROOT"
echo "ISISDATA  = $ISISDATA"
echo "ISIS3DATA = $ISIS3DATA"

## Call the working script with correct parameters
if [ $# -lt 2 ]; then # 1 argument = high-memory
    echo "./lronacPipeline.py --input-folder $1 >& $1/outputLog_p2.txt"
    ./lronacPipeline.py --input-folder $1 >& $1/outputLog_p2.txt
else # 2 argument = low-memory
    echo "./lronacPipeline.py --input-folder $1 $2 >& $1/outputLog_p1.txt"
    ./lronacPipeline.py --input-folder $1 $2 >& $1/outputLog_p1.txt
fi
