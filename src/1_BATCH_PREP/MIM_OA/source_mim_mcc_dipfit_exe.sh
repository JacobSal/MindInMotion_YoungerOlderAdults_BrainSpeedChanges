#!/bin/bash
#SBATCH --job-name=EXE_DIP # Job name
#SBATCH --mail-type=ALL # Mail events (NONE  BEGIN  END  FAIL  ALL)
#SBATCH --mail-user=jsalminen@ufl.edu # Where to send mail
#SBATCH --nodes=1 # Use one node
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=48 # Number of CPU cores per task
#SBATCH --mem-per-cpu=20000mb # Total memory limit
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/_hpg_logs/%j_mcc_dipole_fit_exe.log # Standard output
#SBATCH --account=dferris # Account name
#SBATCH --qos=dferris-b # Quality of service name
#SBATCH --partition=hpg-default # cluster to run on  use slurm command "sinfo -s"; bigmem
# sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/source_mim_mcc_dipfit_exe.sh

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
# (EDIT)!
export DONOT_RECREATE=true;
# cond_vals=$(echo "[1.65,0.33,0.33,0.01,0.126,2.5*10^-14]"); # csf, gray, scalp, skull, white, air
export cond_vals=$(echo "[1.65,0.33,0.33,0.0042,0.126,2.5*10^-14]"); # csf, gray, scalp, skull, white, air
# export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/08202023_OAN82_iccRX0p65_iccREMG0p4_changparams"
# export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/08202023_OAN82_iccRX0p65_iccREMG0p3_newparams"
export SUBJ_EEG="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/08202023_OAN82_iccRX0p60_iccREMG0p3_newparams"
# SET SUBJECT DIRECTORIES
export MCC_PATH="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/"
export SUBJ_HEADMOD="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset"

export SUBJ_RUN=("H1002" "H1004" "H1007" "H1009"
 "H1010" "H1011" "H1012" "H1013" "H1017" "H1018" "H1019"
 "H1020" "H1022" "H1024" "H1025" "H1026" "H1027" "H1029" "H1030"
 "H1031" "H1032" "H1033" "H1034" "H1035"
 "H1036" "H1037" "H1038" "H1039" "H1041"
 "H1042" "H1044" "H1045" "H1046" "H1047" "H1048"
 "H2002" "H2007" "H2008"
 "H2013" "H2015" "H2017" "H2020" "H2021"
 "H2022" "H2023" "H2025" "H2026" "H2027"
 "H2033" "H2034" "H2037" "H2038" "H2039"
 "H2042" "H2052" "H2059" "H2062" "H2082"
 "H2090" "H2095" "H2111" "H2117"
 "H3029" "H3034" "H3039" "H3053"
 "H3063" "H3072" "H3092" "H3103" "H3107"
 "NH3006" "NH3007" "NH3008" "NH3010" "NH3021"
 "NH3026" "NH3030" "NH3036" "NH3040"
 "NH3041" "NH3043" "NH3054"
 "NH3055" "NH3058" "NH3059" "NH3066"
 "NH3068" "NH3069" "NH3070" "NH3074"
 "NH3076" "NH3086" "NH3090" "NH3102"
 "NH3104" "NH3105" "NH3106" "NH3108" "NH3110"
 "NH3112" "NH3113" "NH3114" "NH3123" "NH3128") # JACOB SAL(08/23/2023)
# cd /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit
# %% SET ENVIORNMENT VARIABLES
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

# %% LOOP through a particular cohort of subjects
# CSF_VALS=(1.65);
# GRAY_VALS=(0.33);
# SCALP_VALS=(0.33);
# SKULL_VALS=(0.01,0.0042);
# WHITE_VALS=(0.126);
# AIR_VALS=(2.5*10^-14);
for s in ${SUBJ_RUN[@]};
do
	export curr_f=$SUBJ_EEG/$s/head_model/dipfit_struct.mat
	export mri_f=$SUBJ_HEADMOD/$s/MRI
	export set_f=$SUBJ_EEG/$s/clean/*.set
	export out_f=$SUBJ_EEG/$s/head_model/
	# %% printouts
	echo "Processing Subject $s"
	echo "MRI folder: $SUBJ_HEADMOD/$s/MRI"
	echo "ICA .set file path: $SUBJ_EEG/$s/clean/*.set"
	if test -f "$curr_f" && $DONOT_RECREATE;
	then
		echo "$s headmodel file already generated."
	else
		echo "Calculating Headmodel..."
		echo $mri_f
		echo $set_f
		echo $out_f
		# %% create output folder for source.mat
		mkdir $SUBJ_EEG/"$s"/head_model/
		if $DONOT_RECREATE
		then
			force_recreate=0;
		else
			force_recreate=1;
		fi
		# %% run program
		eval $MCC_PATH/_out/mim_mcc_dipfit "$mri_f" "$set_f" "$out_f" "$cond_vals" "FORCE_RECREATE" "$force_recreate"
		wait
		echo "done: $s"
	fi
done
exit

# $MCC_PATH/_out/run_mim_mcc_dipfit.sh 
