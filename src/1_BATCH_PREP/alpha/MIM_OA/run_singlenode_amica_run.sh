#!/bin/bash
#SBATCH --job-name=ALL_AMICA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsalminen@ufl.edu
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/%j_amica_out.log
#SBATCH --nodes=64
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000mb
#SBATCH --distribution=cyclic:cyclic
#SBATCH --time=48:00:00
#SBATCH --account=dferris
#SBATCH --qos=dferris-b
#SBATCH --partition=hpg2-compute
# NOTE: (04/22/2023) SalminenJ, Seems to time out after ~15 subject runs (6hr time limit). Moving to a 48hr cycle ( I think this is max for hpg2-compute).
SUBJ_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/05192023_YAN33_OAN79_prep_verified"
# REGEXP_SUBJ_DIR="$SUBJ_DIR/*/clean/*.sh"
REGEXP_SUBJ_DIR="$SUBJ_DIR/*/amica/*.param"
cd SUBJ_DIR
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "slurm_mem_per_cpu $SLURM_MEM_PER_CPU"
echo "slurm_mem_per_gpu $SLURM_MEM_PER_GPU"
echo "slurm_mem_per_node $SLURM_MEM_PER_NODE"
module load ufrc
module load intel/2020 openmpi/4.0.3

for f in $REGEXP_SUBJ_DIR
do
# FAILSAFE
# Check if "$f" FILE exists and is a regular file and then only copy it #
  if [ -f "$f" ]
  then
    echo "Processing $f file..."
    srun --mpi=pmix_v3 /blue/dferris/share/s.peterson/test/AMICA_15/amica15ub $f
  else
    echo "Warning: Some problem with \"$f\""
  fi
done

