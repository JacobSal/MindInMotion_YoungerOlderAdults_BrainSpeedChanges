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
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/_hpg_logs/%j_mcc_dipole_fit_exe.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"; bigmem
# sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/source_mim_mcc_dipfit_exe.sh

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
# (EDIT)!
# SET SUBJECT DIRECTORIES
export MCC_PATH="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/"

#%% SET ENVIORNMENT VARIABLES
echo Setting up environment variables
echo ---
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;
echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
# Preload glibc_shim in case of RHEL7 variants
test -e /usr/bin/ldd &&  ldd --version |  grep -q "(GNU libc) 2\.17"  \
		&& export LD_PRELOAD="${MCRROOT}/bin/glnxa64/glibc-2.17_shim.so"
export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}${MCRROOT}/bin/glnxa64/libmwlaunchermain.so"
eval $MCC_PATH/_out/mim_mcc_dipfit "$json_file"
