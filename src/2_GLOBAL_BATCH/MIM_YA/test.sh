#!/bin/bash
#SBATCH --job-name=TEST
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
$parent_path
$SLURM_WORKING_DIR
$SLURM_SUBMIT_DIR
$SLURM_REMOTE_CWD

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    SCRIPT_PATH=$(realpath $0)
fi

# getting location of software_name 
SHARED_PATH=$(dirname $(dirname $(SCRIPT_PATH)))
# separating the software_name from path
SOFTWARE_NAME=$(basename $SHARED_PATH)
# target location to copy project
LOCAL_SOFTWARE_FOLDER='/src'
# corrected path for target
LOCAL_PATH=$LOCAL_SOFTWARE_FOLDER/$SOFTWARE_NAME
$LOCAL_PATH