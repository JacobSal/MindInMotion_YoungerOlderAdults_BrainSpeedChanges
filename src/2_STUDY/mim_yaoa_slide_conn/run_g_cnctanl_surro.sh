#!/bin/bash
#SBATCH --job-name=G_CNCTANL_PROCESS # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --mem-per-cpu=30000mb# Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_slide_conn/_slurm_logs/%j_g_cnctanl_surro.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on, use slurm command 'sinfo -s'
# sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_slide_conn/run_g_cnctanl_surro.sh
module load matlab/2023a

# set linux workspace
# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    TMP_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    TMP_PATH=$(realpath $0)
fi
export SCRIPT_DIR=$(dirname $TMP_PATH)
export STUDY_DIR=$SCRIPT_DIR
export SRC_DIR=$(dirname $(dirname $STUDY_DIR))
cd $STUDY_DIR


echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"


# Create a temporary directory on scratch
mkdir -p $STUDY_DIR/_slurm_scratch/$SLURM_JOB_ID

# Kick off matlab
matlab -nodisplay < $STUDY_DIR/g_cnctanl_surro.m

# Cleanup local work directory
rm -rf $STUDY_DIR/_slurm_scratch/$SLURM_JOB_ID