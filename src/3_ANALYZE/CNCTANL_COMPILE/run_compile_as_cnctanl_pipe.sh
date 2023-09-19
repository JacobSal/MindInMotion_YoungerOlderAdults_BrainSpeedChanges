#!/bin/sh
#SBATCH --job-name=COMPILE_EXE # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --mem-per-cpu=4000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/CNCTANL_COMPILE/_hpg_logs/%j_compile.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on, use slurm command 'sinfo -s'; bigmem

# The location for the output from the compiler
export FILE_OUT="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/CNCTANL_COMPILE/_out"

# Check to make sure directories exist
mkdir $FILE_OUT

# purge modules currently loaded and load matlab compiler resources
module purge
module load matlab/2020b

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"

cd $FILE_IN
matlab -nodisplay < "/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/CNCTANL_COMPILE/compile_as_cnctanl_pipe.m"