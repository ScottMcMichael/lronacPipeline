#!/bin/bash

# This script processes LRONAC pairs from a list.
# - Each processing job is started after its data is downloaded.  
# - Only one download will be going at once but processing jobs
#    will be started up to a maximum number of simultaneous jobs.

# Each job is one compute node on one input folder


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

# Cuts the last 7 characters from each name and
#  appends them together.
# - Needed to get under the 15 character QSUB limit.
function getShortName(){

  DESIRED_SIZE=7

  FULL_NAME=$1 # Input (long) name
  NAME1=$(cut -d "_" -f 1 <<< "$FULL_NAME")
  NAME2=$(cut -d "_" -f 2 <<< "$FULL_NAME")

  SIZE1=${#NAME1}
  SIZE2=${#NAME2}

  DIFF1=$[$SIZE1-$DESIRED_SIZE]
  DIFF2=$[$SIZE2-$DESIRED_SIZE]

  # This removes characters from the beginning of the names
  #  since the ending numbers are likely to change more.
  SHORT1=${NAME1:DIFF1:${#NAME1}}
  SHORT2=${NAME2:DIFF2:${#NAME2}}
 
  SPACER='_'
  OUTPUT=$SHORT1$SPACER$SHORT2
  echo $OUTPUT
}

#---------------------------------------------

INPUT_FILE=/u/smcmich1/data/lronacProduction/lronacPairList60.csv

echo "Starting jobSubmissionScriptProd.sh with input file $INPUT_FILE"

#OUTPUT_FOLDER=/nobackupnfs2/oalexan1/scott/production
OUTPUT_FOLDER=/u/smcmich1/data/lronacProduction

# Limit number of batch jobs to submit
SIMULTANEOUS_JOB_LIMIT=2
SLEEP_TIME=120 # Check number of active processes every five minutes

# Loop until we have exhausted our source file
LAST_LINE=0
while :
do

    # Call the data grabber script
    # - The index of the last line we got is used to start the next search on the correct line.
    # - Unfortunately this script must be completely silent for this to work.
    echo "Downloading next data set..."
    GRABBER_OUTPUT=$(productionDataGrabber.py -i $INPUT_FILE --starting-line $LAST_LINE -o $OUTPUT_FOLDER)
    echo $GRABBER_OUTPUT

    # Check if we actually got data
    LAST_LINE=$(cut -d " " -f 1 <<< "$GRABBER_OUTPUT")
    #echo $LAST_LINE
    if [ "$LAST_LINE" -eq "-1" ]; then
        echo ">>>>> Ran out of data lines in input file, stopping script <<<<<"
        return 0
    fi 

    # Parse results
    FULL_DIRECTORY=$(cut -d " " -f 2 <<< "$GRABBER_OUTPUT")
    POS=`expr "$GRABBER_OUTPUT" : '.*NAC_DTM_'`
    PRETTY_NAME=${GRABBER_OUTPUT:POS}
    PRETTY_NAME=${PRETTY_NAME%?} # Remove trailing slash

    # Verify that the last input file was downloaded before processing this data set
    LOLA_POINTS_PATH=$FULL_DIRECTORY/lolaRdrPoints.csv
    if [ ! -e "$LOLA_POINTS_PATH" ]; then
        echo "LOLA point path $LOLA_POINTS_PATH expected but not found, skipping this data set"
        continue
    fi

    echo "Downloaded data set $PRETTY_NAME, moving to process..."
 
    # Unfortunately the job name is limited to 15 characters!
    JOB_NAME=$(getShortName $PRETTY_NAME)
    echo "Job name = $JOB_NAME"

    STD_OUT_PATH=$FULL_DIRECTORY/stdOutLog.txt
    ERR_OUT_PATH=$FULL_DIRECTORY/errorLog.txt

    # Only run if the last output file is not present and this job is not in the queue
    ASU_STATS_FILE=$FULL_DIRECTORY/results/ASU_LOLA_diff_stats.txt
    LOLA_STATS_FILE=$FULL_DIRECTORY/results/LOLA_diff_stats.txt    
    #echo $ASU_STATS_FILE
    if [ ! -e "$LOLA_STATS_FILE" ] && [ $(checkForJobName $JOB_NAME) -eq 0 ]; then

        # Check the number of currently running jobs
        NUM_ACTIVE_JOBS=$(getNumActiveJobs)
        while : # Loop until we break out
        do
            if [ "$NUM_ACTIVE_JOBS" -lt "$SIMULTANEOUS_JOB_LIMIT" ]; then
 
                echo Num active jobs = $NUM_ACTIVE_JOBS
                echo "Running script for $PRETTY_NAME"
      
                # Submit the job using a westmere (cheap) CPU
                echo qsub -q long -N ${JOB_NAME} -l walltime="80:00:00" -W group_list=s1219 -j oe -e $ERR_OUT_PATH -o $STD_OUT_PATH -S /bin/bash -V -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY       

                qsub -q long -N ${JOB_NAME} -l walltime="80:00:00" -W group_list=s1219 -j oe -e $ERR_OUT_PATH -o $STD_OUT_PATH -S /bin/bash -V -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $FULL_DIRECTORY       

               break # Move on to the next data seta

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

#TODO: Will need to re-run modified v2 script to catch "stragglers" where the download finished but the processing stopped.

# Sample individual submissions

# Submit the job using a westmere (cheap) CPU
#qsub -q normal -N aristarchu2 -e $DATA_FOLDER/ARISTARCHU2/errorOut.txt -o $DATA_FOLDER/ARISTARCHU2/outOut.txt -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -V -l select=1:ncpus=12:model=wes -m eb -- /u/smcmich1/projects/lronacPipeline/jobWrapperV2.sh $DATA_FOLDER/ARISTARCHU2


#exit


