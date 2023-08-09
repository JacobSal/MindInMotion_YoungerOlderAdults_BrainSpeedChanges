#!/bin/bash
#SBATCH --job-name=EXE_DIP # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=48 # Number of CPU cores per task
#SBATCH --mem-per-cpu=15000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/_hpg_logs/%j_mcc_dipole_fit_exe.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"; bigmem

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module purge
module load mcr/2020b
echo $MCRROOT
echo $LD_LIBRARY_PATH
# SET SUBJECT DIRECTORIES
export MCC_PATH="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit"
export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/05192023_YAN33_OAN79_prep_verified"
export SUBJ_HEADMOD="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset"
# SET ENVIORNMENT VARIABLES
export XAPPLRESDIR="$MCRROOT/X11/app-defaults"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$MCRROOT/runtime/glnxa64:$MCRROOT/bin/glnxa64:$MCRROOT/sys/os/glnxa64:$MCRROOT/sys/opengl/lib/glnxa64"
export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}$MCRROOT/bin/glnxa64/glibc-2.17_shim.so"

cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# LOOP through a particular cohort of subjects
for s in ${SUBJ_RUN[@]};
do
	#%% printouts
	echo Processing file $s
	echo $SUBJ_HEADMOD/$s/MRI
	echo $SUBJ_EEG/$s/clean/*.set
	#%% create output folder for source.mat
	mkdir $SUBJ_EEG/"$s"/head_model/
	#%% run program
	# $MCC_PATH/mim_mcc_dipfit $SUBJ_HEADMOD/"$s"/MRI $SUBJ_EEG/"$s"/clean/*.set
	$MCC_PATH/mim_mcc_dipfit $SUBJ_HEADMOD/"$s"/MRI $SUBJ_EEG/"$s"/clean/*.set $SUBJ_EEG/"$s"/head_model/
	echo "Dipole fit for $s is done."
done
