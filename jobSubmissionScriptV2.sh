#!/bin/bash

## This script requests jobs for each of the desired folders

## Each job is one compute node on one input folder

# Move to data folder
cd /nobackupp1/smcmich1/data/lronacPipeline

# Find all the files with the .lola data points downloaded
FILES=$(find . -name '*.csv')

# Move back to execution folder
cd /nobackupp1/smcmich1/projects/lronacPipeline

#for f in $FILES
#do

#    # Get the name of the data set    
#    LOCAL_FOLDER=$(dirname $f)
#    PRETTY_NAME=${LOCAL_FOLDER:2}
    
#    # This is the full path to the data set
#    FULL_DIRECTORY=/nobackupp1/smcmich1/data/lronacPipeline/$PRETTY_NAME

#    # Only run if the last output file is not present
#    ASU_STATS_FILE=$FULL_DIRECTORY/output/ASU_LOLA_diff_stats.txt
#    if [ -f $ASU_STATS_FILE ]; then
#        echo "Running script for $PRETTY_NAME"
        
#        # Submit the job using a westmere (cheap) CPU
#        qsub -q normal -N ${PRETTY_NAME} -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY

#    else
#        echo $PRETTY_NAME = Already finished!
#    fi


#done


# Sample individual submissions

# Submit the job using a westmere (cheap) CPU
qsub -q normal -N rumkerdome2 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapperV2.sh /nobackupp1/smcmich1/data/lronacPipeline/RUMKERDOME2











