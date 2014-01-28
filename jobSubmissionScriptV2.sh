#!/bin/bash

## This script requests jobs for each of the desired folders

## Each job is one compute node on one input folder

# Move to data folder
cd /u/smcmich1/data/lronacPipeline

# Find all the files with the .lola data points downloaded
FILES=$(find . -name '*.csv')

#echo $FILES

# Move back to execution folder
cd /u/smcmich1/projects/lronacPipeline

for f in $FILES
do

    # Get the name of the data set    
    LOCAL_FOLDER=$(dirname $f)
    PRETTY_NAME=${LOCAL_FOLDER:2}
    
    # This is the full path to the data set
    FULL_DIRECTORY=/u/smcmich1/data/lronacPipeline/$PRETTY_NAME

    STD_OUT_PATH=/u/smcmich1/data/lronacPipeline/$PRETTY_NAME/stdOutLog.txt
    ERR_OUT_PATH=/u/smcmich1/data/lronacPipeline/$PRETTY_NAME/errorLog.txt

    # Only run if the last output file is not present
    ASU_STATS_FILE=$FULL_DIRECTORY/results/ASU_LOLA_diff_stats.txt
    if [ -f $ASU_STATS_FILE ]; then
        echo "Running script for $PRETTY_NAME"
        
        # Submit the job using a westmere (cheap) CPU
        qsub -q normal -N ${PRETTY_NAME} -l walltime="8:00:00" -W group_list=s1219 -j oe -e $ERR_OUT_PATH -o $STD_OUT_PATH-S /bin/bash -V -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY

    else
        echo $PRETTY_NAME = Already finished!
    fi


#done


# Sample individual submissions

# Submit the job using a westmere (cheap) CPU
qsub -q normal -N marius3 -e /u/smcmich1/data/lronacPipeline/MARIUS3/errorOut.txt -o /u/smcmich1/data/lronacPipeline/MARIUS3/outOut.txt -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -V -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh /u/smcmich1/data/lronacPipeline/MARIUS3












