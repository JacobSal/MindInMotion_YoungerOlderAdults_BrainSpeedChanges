#!/bin/sh
SUBJ_DIR="/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/07042023_OA_prep_verified"
# REGEXP_SUBJ_DIR="$SUBJ_DIR/*/clean/*.sh"
REGEXP_SUBJ_DIR="$SUBJ_DIR/*/amica/*.sh"

for f in $REGEXP_SUBJ_DIR
do
# FAILSAFE
# Check if "$f" FILE exists and is a regular file and then only copy it #
  if [ -f "$f" ]
  then
    echo "Processing $f file..."
    sbatch $f
  else
    echo "Warning: Some problem with \"$f\""
  fi
done