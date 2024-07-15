#!/bin/bash
#SBATCH --job-name=vscode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=08:00:00 # example 8 hrs

hostname;date;pwd
export XDG_RUNTIME_DIR=${SLURM_TMPDIR}

# HPG install of vscode. Comment the next line out with '#' if using a personal install
module load vscode.

# For either the HPG or personal install of vscode
code tunnel