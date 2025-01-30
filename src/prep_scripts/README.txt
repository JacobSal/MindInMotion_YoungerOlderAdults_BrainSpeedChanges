STEPS:
1) batch_preprocess.m
	* This is a preprocessing script that 
	cleans the EEG data and generates a .set/.fdt file pair that contains all the conditions in the study.
2) run_singlenode_amica_run.sh
	* This generates ICA decompositions from step (1)
3) source_mim_mcc_dipfit_exe.sh
	* This generates dipole fits for independent components in (2)
4*) ants_norm.bash (note: only needs to be ran once)
	* This normalizes each subject's MRI to MNI-space
4) ants_write_dip_csv.m
	* This generates a csv file of all the dipoles fits generated in (3)
5) ants_norm_apply_trans.bash
	* this uses the .csv file in (4) and the MRI file in (4*) to transform the dipole positions to MNI space.
6) ants_norm_dips.m (run_ants_norm_dips.sh)
	* This updates the .set file in (1) to include the Finite Element Dipole fits.


NOTES:
Email from Sumire:
"""
ml gcc/5.2.0 
ml ants

should load ants. Check if loading was successful by entering "antsRegistation"  if it displays you a manual it was correctly loaded.

# APPLICATION
MNI_Template=/Users/sumiresato/Documents/MR_Templates/ANTs_c0Template_T1_IXI555_MNI152_GS_brain.nii 
antsRegistrationSyNQuick.sh -d 3 -f $MNI_template -m ./cat_brain.nii -t s -o ants | tee myRegOutput.txt

Briefly, -f is your MNI template (fixed), -m is your image in native space (moving) -o is the prefix for your output images.

"""