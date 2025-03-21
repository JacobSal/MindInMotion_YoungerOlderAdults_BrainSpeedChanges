%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
ADD_ALL_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
SUBJ_PICS = {{'H3046','H3047','H3073','H3077','H3092', ...
    'NH3023','NH3025','NH3028', ...
    'NH3051','NH3056','NH3071','NH3082','NH3123'}};
%% ===================================================================== %%
%## PARAMS
%- datetime override
study_fName = sprintf('copy_study');
%- soft define
save_dir_bluedrive = [PATHS.data_dir filesep 'MIM_dataset', filesep '_studies' filesep study_fName];
load_dir_bluedrive = [PATHS.data_dir filesep 'MIM_dataset'];
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = load_dir_bluedrive;
%- create new study directory
if ~exist(save_dir_bluedrive,'dir')
    mkdir(save_dir_bluedrive);
end
%% (CONVERT MRI .MAT TO .NII) ========================================== %%
%## MRI .NII ACPC_RS
subj_chars = [SUBJ_PICS{:}];
for subj_i = 1:length(subj_chars)
    mri_mat_fpath = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI'];    
    %-- transfer if exist
    mri_mat = ft_read_mri([mri_mat_fpath filesep 'mri_acpc_rs.mat']);
    ft_write_mri([mri_mat_fpath filesep sprintf('%s_MRI_acpc_rs.nii',subj_chars{subj_i})],mri_mat, ...
        'dataformat','nifti')
end
