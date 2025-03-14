#!/bin/bash
#SBATCH --job-name=HEADRECO # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --mem-per-cpu=8000mb# Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_slurm_logs/run_headreco_batch-%j.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"
# sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/run_prep_headreco_batch.sh

module purge
module load matlab/2020a
#(02/09/2025) JS, need to use 2020a due to compatability issues with the headreco functions. (https://github.com/simnibs/simnibs/issues/343)
export DONOT_RECREATE=true;
cd /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/prep_scripts/
export PATH="/blue/dferris/jsalminen/GitHub/simnibs_install/bin:$PATH"
export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/"
# export SUBJ_RUN=("H1002" "H1004" "H1007" "H1009"
#  "H1010" "H1011" "H1012" "H1013" "H1017" "H1018" "H1019"
#  "H1020" "H1022" "H1024" "H1025" "H1026" "H1027" "H1029" "H1030"
#  "H1031" "H1032" "H1033" "H1034" "H1035"
#  "H1036" "H1037" "H1038" "H1039" "H1041"
#  "H1042" "H1044" "H1045" "H1046" "H1047" "H1048"
#  "H2002" "H2007" "H2008" "H2012_FU"
#  "H2013" "H2015" "H2017" "H2018_FU" "H2020" "H2021"
#  "H2022" "H2023" "H2025" "H2026" "H2027"
#  "H2033" "H2034" "H2037" "H2038" "H2039"
#  "H2042" "H2052" "H2059" "H2062" "H2082"
#  "H2090" "H2095" "H2111" "H2117"
#  "H3029" "H3034" "H3039" "H3053"
#  "H3063" "H3072" "H3092" "H3103" "H3107" "H3120"
#  "NH3006" "NH3007" "NH3008" "NH3010" "NH3021"
#  "NH3026" "NH3030" "NH3036" "NH3040"
#  "NH3041" "NH3043" "NH3054"
#  "NH3055" "NH3058" "NH3059" "NH3066"
#  "NH3068" "NH3069" "NH3070" "NH3074"
#  "NH3076" "NH3086" "NH3090" "NH3102"
#  "NH3104" "NH3105" "NH3106" "NH3108" "NH3110"
#  "NH3112" "NH3113" "NH3114" "NH3123" "NH3128" "NH3129") # JACOB SAL(08/23/2023)
# export SUBJ_RUN=("H2012_FU" "H2018_FU" "H3120" "NH3129")
# export SUBJ_RUN=("NH3023" "NH3028")
export SUBJ_RUN=("NH3028")

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

for s in ${SUBJ_RUN[@]}; # LOOP through a particular cohort of subjects
do
    cd $SUBJ_DIR/$s/MRI #change directory to where the ACPC aligned mri is
    export curr_f=./m2m_"$s"/"$s"_masks_contr.nii
    #%% printouts
	echo "Processing Subject $s"
	echo "MRI nifti: $SUBJ_HEADMOD/$s/MRI"
	if test -f "$curr_f" && $DONOT_RECREATE;
	then
		echo "$s segmentation file already generated."
	else
		echo "Calculating segmentation..."
		headreco preparevols $s "$s"_MRI_acpc_rs.nii && headreco preparecat $s && headreco cleanvols $s
        echo "Segmentation for $s participant is done. Please manually transfer the file generated."
	fi
	# cd $SUBJ_DIR/$s/MRI/Processed_fiducials #change directory to where the ACPC aligned mri is

done



