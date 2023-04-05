#!/bin/sh
#SBATCH --job-name=RUN_EXE # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem-per-cpu=8000mb# Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=72:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_YA_proc/%j_MIM_YA_conn_batch.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=bigmem # cluster to run on, use slurm command 'sinfo -s'; bigmem

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

# purge modules currently loaded and load matlab compiler resources
module purge
module load mcr/2020b 

# LINUX dynamic link libraries
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}\
/apps/matlab/mcr/2020b/v99/runtime/glnxa64:\
/apps/matlab/mcr/2020b/v99/bin/glnxa64:\
/apps/matlab/mcr/2020b/v99/sys/os/glnxa64:\
/apps/matlab/mcr/2020b/v99/extern/bin/glnxa64"

# LINUX dynamic link libraries
export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}\
/apps/matlab/mcr/2020b/v99/bin/glnxa64/glibc-2.17_shim.so"

cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_YA_proc/_out
./main_func var1 var2
