#!/bin/bash
#SBATCH --job-name=EXE_DIP # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=64 # Number of CPU cores per task
#SBATCH --mem-per-cpu=15000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/%j_mcc_dipole_fit_exe.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on, use slurm command 'sinfo -s'; bigmem

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
export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/07042023_OA_prep_verified"
export SUBJ_HEADMOD="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset"
# SET ENVIORNMENT VARIABLES
export XAPPLRESDIR="$MCRROOT/X11/app-defaults"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$MCRROOT/runtime/glnxa64:$MCRROOT/bin/glnxa64:$MCRROOT/sys/os/glnxa64:$MCRROOT/sys/opengl/lib/glnxa64"
export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}$MCRROOT/bin/glnxa64/glibc-2.17_shim.so"

# cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# ./mim_mcc_dipfit $SUBJ_HEADMOD/H2020/MRI $SUBJ_EEG/H2020/clean/H2020_cleanEEG_EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set
SUBJ_H2=("H2017" "H2002" "H2007" "H2008" "H2013" "H2015" "H2020" "H2021" "H2022" "H2023" "H2025" "H2026" "H2027" "H2033" "H2034" "H2036" "H2037" "H2038" "H2039" "H2041" "H2042" "H2052" "H2059" "H2062" "H2072" "H2082" "H2090" "H2095" "H2111" "H2117") # JACOB SAL(02/23/2023)
cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# LOOP through a particular cohort of subjects
for s in ${SUBJ_H2[@]};
do
	echo "Processing file $s"
	echo "$SUBJ_HEADMOD/$s/MRI"
	echo "$SUBJ_EEG/$s/clean/$s_cleanEEG_EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set"
	$MCC_PATH/mim_mcc_dipfit $SUBJ_HEADMOD/"$s"/MRI $SUBJ_EEG/"$s"/clean/"$s"_cleanEEG_EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set
	echo "Dipole fit for $s is done."
done
