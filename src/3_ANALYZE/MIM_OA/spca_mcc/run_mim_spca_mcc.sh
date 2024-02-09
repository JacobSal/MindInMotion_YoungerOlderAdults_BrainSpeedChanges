#!/bin/bash
#SBATCH --job-name=EXE_SPCA # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=20 # Number of CPU cores per task
#SBATCH --mem-per-cpu=20000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc/_hpg_logs/%j_mim_spca_mcc.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"; bigmem
# sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc/run_mim_spca_mcc.sh

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
# (EDIT)! SET SUBJECT DIRECTORIES
export ICA_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams"
export MCC_PATH="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc"
# export JSON_FPATH="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc/spca_mcc_params"
export JSON_FPATH=""
export STUDY_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/test_spca"
# export STUDY_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/01122024_spca_analysis"

# SET ENVIORNMENT VARIABLES
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
# export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}/lib64/libgfortran.so.3"

$MCC_PATH/_out/mim_spca_mcc "$ICA_DIR" "$STUDY_DIR"  "$JSON_FPATH"