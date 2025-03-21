#!/bin/bash
#SBATCH --job-name=ALL_AMICA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsalminen@ufl.edu
#SBATCH --output=/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_slurm_logs/%j_amica_out.log
#SBATCH --nodes=64
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000mb
#SBATCH --distribution=cyclic:cyclic
#SBATCH --time=48:00:00
#SBATCH --account=dferris
#SBATCH --qos=dferris-b
#SBATCH --partition=hpg-default
# NOTE: (04/22/2023) SalminenJ, Seems to time out after ~15 subaject runs (6hr time limit). Moving to a 48hr cycle ( I think this is max for hpg2-compute).
# sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/b_run_singlenode_amica.sh

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
module load intel/2020 openmpi/4.1.5

#%% PARAMS
export DONOT_RECREATE=false;

# export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/01132024_antsnorm_iccREEG0p65_iccREMG0p4_skull0p0042"
# export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams"
export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/_data/MIM_dataset/_studies/02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej"

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
 "NH3112" "NH3113" "NH3114" "NH3123" "NH3128" "NH3129"
 "H3046" "H3047" "H3073" "H3077"
 "H3092" "NH3023" "NH3025" "NH3027"
 "NH3028" "NH3051" "NH3056" "NH3071"
 "NH3082" "NH3123") # JACOB SAL(08/23/2023); JS Updated (02/21/2025)

# export SUBJ_RUN=("H3046" "H3047" "H3073" "H3077" 
# "H3092" "NH3023" "NH3025" "NH3027" 
# "NH3028" "NH3051" "NH3056" "NH3071" 
# "NH3082" "NH3123")
#%% LOOP
for s in ${SUBJ_RUN[@]}
do
	export param_f=$SUBJ_DIR/$s/clean/*.param
	export chk_f=$SUBJ_DIR/$s/clean/W
	echo "Processing $s eeg file..."
	if test -f "$chk_f" && $DONOT_RECREATE;
	then
		echo "ICA weight file is already generated: $chk_f"
	else
		if test -f $param_f;
		then
			echo "Calculating ICA weights..."
			#srun --mpi=pmix_v3 /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/AMICA_15/amica15ub $param_f
			srun --mpi=pmix_v3 /blue/dferris/share/s.peterson/test/AMICA_15/amica15ub $param_f
			#%% This "wait" may be needed to ensure shared libraries aren't accessed multiple times?
			wait
			echo "done: $s"
		else
			echo "$param_f does not exist. Preprocess EEG for $s."
		fi
	fi
done

