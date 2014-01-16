#!/bin/bash

## This script requests jobs for each of the desired folders

## Each job is one compute node on one input folder

#set tooldir=/nobackupp1/smcmich1/projects/lronacPipeline
#set datadir=/nobackupp1/smcmich1/data/lronacPipeline

# Move to data folder
cd /nobackupp1/smcmich1/data/lronacPipeline

# Find all the files with the .lola data points downloaded
FILES=$(find . -name '*.csv')

# Move back to execution folder
cd /nobackupp1/smcmich1/projects/lronacPipeline

for f in $FILES
do
    
    FOLDER=$(dirname $f)

    PRETTY_NAME=${FOLDER:2}

    # Count up the number of diff_stats.txt files which are present
    NUM_STATS_FILES=$(ls ~/data/lronacPipeline/$PRETTY_NAME | grep diff_stats.txt | wc -l)

    # Only run if one or more of the diff_stats.txt file are missing    
    if [ $NUM_STATS_FILES -lt 3 ]; then
        echo $PRETTY_NAME = $NUM_STATS_FILES
        
         subjob=$(qsub -q normal -N ${PRETTY_NAME}_1 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapperV1.sh /nobackupp1/smcmich1/data/lronacPipeline/$PRETTY_NAME --low-mem)

         qsub -q normal -N ${PRETTY_NAME}_2 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:bigmem=true:model=wes -m eb -W depend=afterok:$subjob -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapperV1.sh /nobackupp1/smcmich1/data/lronacPipeline/$PRETTY_NAME

    else
        echo $PRETTY_NAME = Already finished!
    fi


#done


# Sample individual submissions

# First part (up to pc_align)
#subjob=$(qsub -q normal -N rumkerdome2_1 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapper.sh /nobackupp1/smcmich1/data/lronacPipeline/RUMKERDOME2 --low-mem)

# Second part (pc_align to finish)
#qsub -q normal -N rumkerdome2_2 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:bigmem=true:model=wes -m eb -W depend=afterok:$subjob -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapper.sh /nobackupp1/smcmich1/data/lronacPipeline/RUMKERDOME2


# First part (up to pc_align)
#subjob=$(qsub -q normal -N kuglerridge_1 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:model=wes -m eb -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapper.sh /nobackupp1/smcmich1/data/lronacPipeline/KUGLERRIDGE --low-mem)

# Second part (pc_align to finish)
#qsub -q normal -N kuglerridge_2 -l walltime="8:00:00" -W group_list=s1219 -j oe -S /bin/bash -C $PWD -l select=1:ncpus=12:bigmem=true:model=wes -m eb -W depend=afterok:$subjob -- /nobackupp1/smcmich1/projects/lronacPipeline/jobWrapper.sh /nobackupp1/smcmich1/data/lronacPipeline/KUGLERRIDGE












