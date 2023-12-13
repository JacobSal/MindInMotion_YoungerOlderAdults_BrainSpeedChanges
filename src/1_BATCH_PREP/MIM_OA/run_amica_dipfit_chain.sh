#!/bin/bash
cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/

# sbatch run_amica_dipfit_chain.sh
# sentence=$(sbatch --time=21-00:00:00 run_amica_dipfit_chain.sh) # get the output from sbatch
# stringarray=($sentence)                            # separate the output in words
# jobid=(${stringarray[3]})                          # isolate the job ID
# sentence="$(squeue -j $jobid)"            # read job's slurm status
# stringarray=($sentence) 
# jobstatus=(${stringarray[12]})            # isolate the status of job number jobid
# if [ "$jobstatus" = "R" ];then
	# sbatch source_mim_mcc_dipfit_exe.sh
  # insert here relevant code to run next job
# fi

srun run_amica_dipfit_chain.sh
srun source_mim_mcc_dipfit_exe.sh