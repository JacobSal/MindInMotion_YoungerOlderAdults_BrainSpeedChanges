#!/bin/bash
#SBATCH --job-name=MCC_DIPFIT # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=48 # Number of CPU cores per task
#SBATCH --mem-per-cpu=12000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=05:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_slurm_logs/%j_test_c_run_mim_mcc_dipfit.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"; bigmem
# sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/test_c_run_mim_mcc_dipfit.sh

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
export SRC_DIR=$(dirname $SCRIPT_DIR)
export STUDY_DIR=$SCRIPT_DIR
cd $STUDY_DIR

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module load mcr/2020a
echo $MCRROOT
echo $LD_LIBRARY_PATH

# FUNCTION PARAMETERS
export DONOT_RECREATE=false;
cond_vals=$(echo "[1.65,0.33,0.33,0.01,0.126,2.5*10^-14]"); # csf, gray, scalp, skull, white, air
# export cond_vals=$(echo "[1.65,0.33,0.33,0.0042,0.126,2.5*10^-14]"); # csf, gray, scalp, skull, white, air
export MCC_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_Functions/src/_mcc_funcs/mcc_dipfit"

# SET SUBJECT DIRECTORIES
# export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams"
export SUBJ_HEADMOD="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset"
export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej"

# %% PART 2
 export s=("H1004")
 # JACOB SAL(08/23/2023); JS Updated (02/21/2025)

# %% LOOP through a particular cohort of subjects
export curr_f=$SUBJ_EEG/$s/head_model/dipfit_struct.mat
export mri_f=$SUBJ_HEADMOD/$s/MRI
# export mri_f=$SUBJ_HEADMOD/$s/MRI/m2m_$s
export set_f=$SUBJ_EEG/$s/clean/*.set
export out_f=$SUBJ_EEG/$s/head_model/
# %% printouts
echo "Processing Subject $s"
echo "MRI folder: $SUBJ_HEADMOD/$s/MRI"
echo "ICA .set file path: $SUBJ_EEG/$s/clean/*.set"
if test -f "$curr_f" && $DONOT_RECREATE;
then
	echo "$s headmodel file already generated."
else
	echo "Calculating Headmodel..."
	echo $mri_f
	echo $set_f
	echo $out_f
	# %% create output folder for source.mat
	mkdir $SUBJ_EEG/"$s"/head_model/
	if $DONOT_RECREATE
	then
		force_recreate=0;
	else
		force_recreate=1;
	fi
	# %% run program
	eval $MCC_DIR/_out/run_mcc_dipfit.sh $MCRROOT "$mri_f" "$set_f" "$out_f" "VOL_CONDUCTIVITIES"  "$cond_vals" "FORCE_RECREATE" "$force_recreate"
	wait
	echo "done: $s"
fi
exit

# $MCC_PATH/_out/run_mim_mcc_dipfit.sh 
