#!/bin/bash

## This script requests jobs for each of the desired folders

## Each job is one compute node on one input folder


# Return the current number of active PBS jobs
function getNumActiveJobs(){
  NUM_LINES=$(qstat -u smcmich1 | wc -l)
  if [ $NUM_LINES -eq 0 ]; then
      echo 0
  else # Lines returned means at least one active job
      NUM_TASKS=$((NUM_LINES - 3)) # num tasks = num lines - header lines
      echo $NUM_TASKS
  fi
}

# Check if the given job name is already in the queue
function checkForJobName(){

 NUM_LINES=$(qstat -u smcmich1 | grep $1 | wc -l)
 if [ $NUM_LINES -gt 0 ]; then # Job found
     echo 1 
 else # Job not found in queue
     echo 0
 fi
}

#---------------------------------------------

DATA_FOLDER=/nobackupnfs2/oalexan1/scott
#DATA_FOLDER=/u/smcmich1/data/lronacPipeline

# Move to data folder
cd $DATA_FOLDER

## Find all the files with the ASU data downloaded
#FILES=$(find . -name 'NAC_DTM*')

# Find all the files with the LOLA data data downloaded
FILES=$(find . -name 'lolaRdrPoints.csv')


#echo $FILES
echo Obtained file list

# Move back to execution folder
cd /u/smcmich1/projects/lronacPipeline

# Limit number of batch jobs to submit
SIMULTANEOUS_JOB_LIMIT=6
SLEEP_TIME=120 # Check number of active processes every five minutes

for f in $FILES
do

    # Get the name of the data set    
    LOCAL_FOLDER=$(dirname $f)
    PRETTY_NAME=${LOCAL_FOLDER:2}
    
    #echo $LOCAL_FOLDER
    #echo $PRETTY_NAME

    # This is the full path to the data set
    FULL_DIRECTORY=$DATA_FOLDER/$PRETTY_NAME

    STD_OUT_PATH=$FULL_DIRECTORY/stdOutLog.txt
    ERR_OUT_PATH=$FULL_DIRECTORY/errorLog.txt

    # Only run if the last output file is not present and this job is not in the queue
    #ASU_STATS_FILE=$FULL_DIRECTORY/results/ASU_LOLA_diff_stats.txt
    LOLA_STATS_FILE=$FULL_DIRECTORY/results/LOLA_diff_stats.txt
    #echo $ASU_STATS_FILE
    if [ ! -e "$LOLA_STATS_FILE" ] && [ $(checkForJobName $PRETTY_NAME) -eq 0 ]; then

        # Check the number of currently running jobs
        NUM_ACTIVE_JOBS=$(getNumActiveJobs)
        while : # Loop until we break out
        do
            if [ "$NUM_ACTIVE_JOBS" -lt "$SIMULTANEOUS_JOB_LIMIT" ]; then
 
                echo Num active jobs = $NUM_ACTIVE_JOBS
                echo "Running script for $PRETTY_NAME"
      
                # Submit the job using a westmere (cheap) CPU
                qsub -q normal -N ${PRETTY_NAME} -l walltime="8:00:00" -W group_list=s1219 -j oe -e $ERR_OUT_PATH -o $STD_OUT_PATH -S /bin/bash -V -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY       
               break # Move on to the next data set
            else # Wait for a while
                CURRENT_TIME=$(date +"%T")
                echo "$CURRENT_TIME - Waiting for a job to finish."
                sleep $SLEEP_TIME    
                NUM_ACTIVE_JOBS=$(getNumActiveJobs)
            fi
        done # End waiting loop
        
    fi

done

echo jobSubmissionScript completed!

# Sample individual submissions

# Submit the job using a westmere (cheap) CPU
#qsub -q normal -N aristarchu2 -e $DATA_FOLDER/ARISTARCHU2/errorOut.txt -o $DATA_FOLDER/ARISTARCHU2/outOut.txt -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -V -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $DATA_FOLDER/ARISTARCHU2


#exit


