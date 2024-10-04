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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = fileparts(STUDY_DIR);
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_mri_study');
%% ===================================================================== %%
%## PARAMS
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (TRANSFER IMU/LS CSV DATA TO _STUDIES) ============================== %%
save_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset', filesep '_studies' filesep 'arkaprava_data'];
%## IMU DATA
subj_chars = [SUBJ_PICS{:}];
for subj_i = 45:length(subj_chars)
    %## TRANSFER ALL IMU CSV'S
    try
        folder_to = [save_dir filesep subj_chars{subj_i}];
        folder_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'IMU' filesep 'Raw'];
        transfer_folder(folder_from,folder_to,'*.csv');
    catch e
        fprintf('%s',getReport(e));
    end
    %## TRANSFER ALL LS CSV'S
    % folder_to = [save_dir filesep subj_chars{subj_i}];
    % folder_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'Loadsol' filesep 'Raw'];
    % transfer_folder(folder_from,folder_to,'*.txt');
end
%% (TRANSFER IMU/LS CSV DATA TO MIM_DATASET) =========================== %%
save_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset'];
%## IMU DATA
subj_chars = [SUBJ_PICS{:}];
for subj_i = 1:length(subj_chars)
    %## TRANSFER ALL IMU CSV'S
    folder_to = [save_dir filesep subj_chars{subj_i} filesep 'IMU' filesep 'Raw'];
    folder_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'IMU' filesep 'Raw'];
    transfer_folder(folder_from,folder_to,'*.csv');
    %## TRANSFER ALL LS CSV'S
    % folder_to = [save_dir filesep subj_chars{subj_i} filesep 'Loadsol' filesep 'Raw'];
    % folder_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'Loadsol' filesep 'Raw'];
    % transfer_folder(folder_from,folder_to,'*.txt');
end