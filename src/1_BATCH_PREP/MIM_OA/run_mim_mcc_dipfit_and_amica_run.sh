#!/bin/bash
#SBATCH --job-name=EXE_DIP # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=64 # Use one node
#SBATCH --ntasks=64 # Run a single tasks
#SBATCH --cpus-per-task=1 # Number of CPU cores per task
#SBATCH --mem-per-cpu=1000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/_hpg_logs/%j_mcc_dipole_fit_exe.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg2-compute # cluster to run on  use slurm command "sinfo -s"; bigmem
# sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/run_mim_mcc_dipfit_and_amica_run.sh

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module purge
module load ufrc
module load intel/2020 openmpi/4.0.3

export SUBJ_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/07142023_OAN79_iccRX0p55_iccREMG0p3_changparams"
export BASHSCRIPT_HDM="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/"

# cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# ./mim_mcc_dipfit $SUBJ_HEADMOD/H2020/MRI $SUBJ_EEG/H2020/clean/H2020_cleanEEG_EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set
# SUBJ_DONE=("H2017" "H2002" "H2007" "H2008" "H2013" "H2015" "H2020" "H2021" "H2022" "H2023" "H2025" "H2026" "H2027" "H2033" "H2034" "H2036" "H2037" "H2038" "H2039" "H2041" "H2042" "H2052" "H2059" "H2062")
# SUBJ_RUN=("H2072" "H2082" "H2090" "H2095" "H2111" "H2117""H3018" "H3029" "H3034" "H3039" "H3042" "H3046" "H3047" "H3053" "H3063" "H3072" "H3073" "H3077" "H3092" "H3103" "H3107" "H3120" "NH3006"  "NH3007"  "NH3008"  "NH3010" "NH3021"  "NH3025"  "NH3026" "NH3030"  "NH3036" "NH3041"  "NH3043"  "NH3051"  "NH3054"  "NH3055"  "NH3056"  "NH3058" "NH3059"  "NH3066"  "NH3068"  "NH3069"  "NH3070"  "NH3071"  "NH3074" "NH3076"  "NH3082"  "NH3086"  "NH3090"  "NH3102"  "NH3104"  "NH3105"  "NH3106" "NH3108"  "NH3110"  "NH3112"  "NH3113"  "NH3114" "NH3123" "NH3128") # JACOB SAL(02/23/2023)
# "H1002" "H1004" "H1007" "H1009" "H1010" "H1011" "H1012" "H1013" "H1017" "H1018" "H1019"
# "H1020" "H1022" "H1024" "H1026" "H1027" "H1029" "H1030" "H1031" "H1032" "H1033" "H1034" "H1035"
# "H1036" "H1037" "H1038" "H1039" "H1041" "H1042" "H1044" "H1045" "H1047" "H1047" 
export SUBJ_RUN=("H2002" "H2007"
 "H2008" "H2013" "H2015" "H2017" "H2020" "H2021" "H2022" "H2023" "H2025" "H2026" "H2027" "H2033" "H2034"
 "H2037" "H2038" "H2039" "H2042" "H2052" "H2059" "H2062" "H2072" "H2082" "H2090"
 "H2095" "H2111" "H2117" "H3018" "H3029" "H3034" "H3039" "H3042" "H3046" "H3047" "H3053" "H3063"
 "H3072" "H3073" "H3077" "H3092" "H3103" "H3107" "H3120" "NH3006" "NH3007" "NH3008" "NH3010"
 "NH3021" "NH3025" "NH3026" "NH3030" "NH3036" "NH3041" "NH3043" "NH3051" "NH3054" "NH3055"
 "NH3056" "NH3058" "NH3059" "NH3066" "NH3068" "NH3069" "NH3070" "NH3071" "NH3074" "NH3076"
 "NH3082" "NH3086" "NH3090" "NH3102" "NH3104" "NH3105" "NH3106" "NH3108" "NH3110" "NH3112"
 "NH3113" "NH3114" "NH3123" "NH3128") # JACOB SAL(05/19/2023)
# export SUBJ_RUN=("NH3025" "NH3076" "NH3082" "NH3086" "NH3090" "NH3102" "NH3104" "NH3105" "NH3106" "NH3108"
 # "NH3110" "NH3112" "NH3113" "NH3114" "NH3123" "NH3128")

cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# LOOP through a particular cohort of subjects
for s in ${SUBJ_RUN[@]};
do
	#%% create output folder for source.mat
	# mkdir $SUBJ_EEG/"$s"/head_model/
	#%% run program
	export shf=$SUBJ_DIR/$s/clean/*.param
	echo "Processing AMICA for $shf file..."
	srun --mpi=pmix_v3 /blue/dferris/share/s.peterson/test/AMICA_15/amica15ub $shf
	#%% printouts
	echo "$s AMICA is done"
done

echo "Running headmodel execution"
source $BASHSCRIPT_HDM/source_mim_mcc_dipfit_exe.sh
echo "Dipole fit and AMICA for $s is done."
