#!/bin/bash
#SBATCH --job-name=ANTS_MRI_NORM # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem-per-cpu=10000mb# Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=16:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_BrainSpeedChanges/src/_slurm_logs/%j_d_run_ants_norm.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on, use slurm command 'sinfo -s'
#%% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/d_run_ants_norm.bash
ml gcc/5.2.0
ml ants

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
export SRC_DIR=$(dirname $(dirname $TMP_PATH))
export SCRIPT_DIR=$SRC_DIR/prep_scripts
export STUDY_DIR=$SRC_DIR
cd $STUDY_DIR

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

export DONOT_RECREATE=true;
# SET SUBJECT DIRECTORIES
export SUBJ_HEADMOD="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset"
export MNI_TEMPLATE="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/_resources/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii"

export SUBJ_RUN=("H1002" "H1004" "H1007" "H1009"
 "H1010" "H1011" "H1012" "H1013" "H1017" "H1018" "H1019"
 "H1020" "H1022" "H1024" "H1025" "H1026" "H1027" "H1029" "H1030"
 "H1031" "H1032" "H1033" "H1034" "H1035"
 "H1036" "H1037" "H1038" "H1039" "H1041"
 "H1042" "H1044" "H1045" "H1046" "H1047" "H1048"
 "H2002" "H2007" "H2008" "H2012_FU"
 "H2013" "H2015" "H2017" "H2018_FU" "H2020" "H2021"
 "H2022" "H2023" "H2025" "H2026" "H2027"
 "H2033" "H2034" "H2037" "H2038" "H2039"
 "H2042" "H2052" "H2059" "H2062" "H2082"
 "H2090" "H2095" "H2111" "H2117"
 "H3029" "H3034" "H3039" "H3053"
 "H3063" "H3072" "H3092" "H3103" "H3107" "H3120"
 "NH3006" "NH3007" "NH3008" "NH3010" "NH3021"
 "NH3026" "NH3030" "NH3036" "NH3040"
 "NH3041" "NH3043" "NH3054"
 "NH3055" "NH3058" "NH3059" "NH3066"
 "NH3068" "NH3069" "NH3070" "NH3074"
 "NH3076" "NH3086" "NH3090" "NH3102"
 "NH3104" "NH3105" "NH3106" "NH3108" "NH3110"
 "NH3112" "NH3113" "NH3114" "NH3123" "NH3128" "NH3129") # JACOB SAL(08/23/2023)
# export SUBJ_RUN=("H2117")
for s in ${SUBJ_RUN[@]};
do
	# export mri_f=$SUBJ_HEADMOD/$s/MRI/mri_acpc_rs.nii
    export mri_f=$SUBJ_HEADMOD/$s/MRI/"$s"_MRI_acpc_rs.nii
	# %% printouts
	echo "Processing Subject $s"
	echo "MRI nifti: $SUBJ_HEADMOD/$s/MRI"
	if test -f "$curr_f" && $DONOT_RECREATE;
	then
		echo "$s headmodel file already generated."
	else
		echo "Calculating Headmodel..."
		echo $mri_f
		cd "$SUBJ_HEADMOD/$s/MRI"
		antsRegistrationSyNQuick.sh -d 3 -f "$MNI_TEMPLATE" -m $mri_f -t s -o ants | tee myRegOutput.txt
		echo "done: $s"
	fi
done
exit
