#!/bin/bash

## This script requests jobs for each of the desired folders

## Each job is one compute node on one input folder

# Return the current number of active PBS jobs
function getNumActiveJobs(){
  NUM_LINES="qstat -u smcmich1 | wc -l"
  NUM_TASKS=$((NUM_LINES - 3))
  return $NUM_TASKS
}


#---------------------------------------------

DATA_FOLDER=/nobackupnfs2/oalexan1/scott
#DATA_FOLDER=/u/smcmich1/data/lronacPipeline


# Move to data folder
cd $DATA_FOLDER

# Find all the files with the .lola data points downloaded
FILES=$(find . -name 'NAC_DTM*')

#echo $FILES

# Move back to execution folder
cd /u/smcmich1/projects/lronacPipeline

# Limit number of batch jobs to submit
SIMULTANEOUS_JOB_LIMIT=8
SLEEP_TIME=300 # Check number of active processes every five minutes
declare -i LIMIT=8
declare -i ONE=1

for f in $FILES
do

    # Get the name of the data set    
    LOCAL_FOLDER=$(dirname $f)
    PRETTY_NAME=${LOCAL_FOLDER:2}
    
    # This is the full path to the data set
    FULL_DIRECTORY=$DATA_FOLDER/$PRETTY_NAME

    STD_OUT_PATH=$FULL_DIRECTORY/stdOutLog.txt
    ERR_OUT_PATH=$FULL_DIRECTORY/errorLog.txt

    # Only run if the last output file is not present
    ASU_STATS_FILE=$FULL_DIRECTORY/results/ASU_LOLA_diff_stats.txt
    #echo $ASU_STATS_FILE
    if [ ! -e "$ASU_STATS_FILE" ]; then
        echo "Running script for $PRETTY_NAME"


        # Check the number of currently running jobs
        NUM_ACTIVE_JOBS=getNumActiveJobs()
        do
            if [ $NUM_ACTIVE_JOBS -lt $SIMULTANEOUS_JOB_LIMIT]
       
                # Submit the job using a westmere (cheap) CPU
                qsub -q normal -N ${PRETTY_NAME} -l walltime="8:00:00" -W group_list=s1219 -j oe -e $ERR_OUT_PATH -o $STD_OUT_PATH -S /bin/bash -V -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY       
                    
            else # Wait for a while
                CURRENT_TIME=date +"%T"
                echo "$CURRENT_TIME - Waiting for a job to finish."
                sleep $SLEEP_TIME    
                NUM_ACTIVE_JOBS=getNumActiveJobs()
            fi
        done # End waiting loop
        
#        LIMIT=$(($LIMIT-$ONE))        
    fi

    #if [ $LIMIT -eq  0 ]; then
    #    exit 0
    #fi

done


# Sample individual submissions

# Submit the job using a westmere (cheap) CPU
#qsub -q normal -N marius3 -e /u/smcmich1/data/lronacPipeline/MARIUS3/errorOut.txt -o /u/smcmich1/data/lronacPipeline/MARIUS3/outOut.txt -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -V -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh /u/smcmich1/data/lronacPipeline/MARIUS3


#exit









