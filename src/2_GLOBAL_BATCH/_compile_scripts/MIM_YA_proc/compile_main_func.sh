#!/bin/sh
# Submodules that you want to add to compile your function (adds all directories within)
A_SUBMODS_IN=("/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/Granger_Geweke_Causality"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/_external"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/cnctanl"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/eeglab"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/ES"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/MIM"
"/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/NJ")

for i in ${A_SUBMODS_IN[@]}; do
    echo "$i"
done

# Submodules that you want to add to compile your function (adds only immediate directory)
I_SUBMODS_IN=("/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0")

for i in ${I_SUBMODS_IN[@]}; do
    echo "$i"
done

# The location of the function
FILE_IN=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/iota/MIM_YA_proc/

# The location for the output from the compiler
FILE_OUT=/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/iota/MIM_YA_proc/_out

# Check to make sure directories exist
mkdir $FILE_OUT

# purge modules currently loaded and load matlab compiler resources
module purge
module load matlab/2020a

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"

cd $FILE_IN

mcc -m main_func_v2.m -d $FILE_OUT -a ${A_SUBMODS_IN} -I ${I_SUBMODS_IN} -R -singleCompThread 
