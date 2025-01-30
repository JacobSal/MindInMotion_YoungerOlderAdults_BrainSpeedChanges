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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
%% ===================================================================== %%
%## PARAMS
%- datetime override 
% dt = '03232023_AS_Bishoy';
dt = 'test_ls_alg';
study_fName = sprintf('copy_study');
%- soft define
save_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset' filesep '_studies' filesep sprintf('%s',dt)];
load_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ASSIGN FILE PATHS
M_MIND_IN_MOTION_DIR    = load_dir; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
files_loadsol           = cell(1,length([SUBJ_ITERS{:}]));
subj_save_ls            = cell(1,length([SUBJ_ITERS{:}]));
files_imu               = cell(1,length([SUBJ_ITERS{:}]));
subj_save_imu           = cell(1,length([SUBJ_ITERS{:}]));
subj_name               = cell(1,length([SUBJ_ITERS{:}]));
stack_iter              = 0;
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    cnt = stack_iter + 1;
    sub_idx = SUBJ_ITERS{group_i};
    for subj_i = sub_idx % 1:length(SUBJ_PICS{group_i})
        %- assign filepaths for imu and loadsol
        folder_ls_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Raw'];
        folder_imu_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Raw'];
        folder_ls_save = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Imported'];
        folder_imu_save = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Imported'];
        %- store filepaths in cells
        files_loadsol{cnt} = folder_ls_to; %[files_loadsol{group_i}; {folder_ls_to}];
        files_imu{cnt} = folder_imu_to; %[files_imu{group_i}; {folder_imu_to}];
        subj_save_ls{cnt} = folder_ls_save; %[subj_save_ls{group_i}; {folder_ls_save}];
        subj_save_imu{cnt} = folder_imu_save; %[subj_save_imu{group_i}; {folder_imu_save}];
        subj_name{cnt} = SUBJ_PICS{group_i}{subj_i};
        %- create new study directory
        if ~exist(folder_ls_to,'dir')
            mkdir(folder_ls_to);
        end
        %- create new study directory
        if ~exist(folder_ls_save,'dir')
            mkdir(folder_ls_save);
        end
        %- create new study directory
        if ~exist(folder_imu_to,'dir')
            mkdir(folder_imu_to);
        end
        %- create new study directory
        if ~exist(folder_imu_save,'dir')
            mkdir(folder_imu_save);
        end
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%% MERGE LIVEAMPS TO EEG STRUCT
subj_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\H1007';
liveamp_fpaths = {[subj_fpath filesep 'EEG' filesep 'Raw' filesep 'Cort_GY'],...
    [subj_fpath filesep 'EEG' filesep 'Raw' filesep 'Cort_WR'],...
    [subj_fpath filesep 'EEG' filesep 'Raw' filesep 'Noise_GY'],...
    [subj_fpath filesep 'EEG' filesep 'Raw' filesep 'Noise_WR']};
[EEG] = mim_merge_liveamps(liveamp_fpaths);

% datetime('20190911-152201','InputFormat','yyyyMMdd-HHmmss')
outputFileName = [fileNameNoExt,'.set']; %(ex: H1002_EEG.set)
disp(['Saving: ', fullfile(EEGoutputFolder,outputFileName)]);
EEG = pop_saveset( EEG, 'filepath', EEGoutputFolder, 'filename', fileNameNoExt);
disp('Done saving.');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab redraw;

%% CONVERT LOADSOL TXT TO EVENTS
% subj_chars = [SUBJ_PICS{:}];
ls_save_dir = [save_dir filesep '_figs' filesep 'LS'];
if ~exist(ls_save_dir,'dir')
    mkdir(ls_save_dir);
end
ALLEEG_LS = [];
for cnt = 1:length([SUBJ_ITERS{:}])
    [loadsol_table, rec_start_struct] = convert_loadsol_to_table(files_loadsol{cnt});
    for f_i = 1:length(loadsol_table)
        if ~exist([ls_save_dir filesep subj_name{cnt}],'dir')
            mkdir([ls_save_dir filesep subj_name{cnt}])
        end
        EEG_LS = create_loadsol_struct(loadsol_table{f_i},rec_start_struct(f_i),[ls_save_dir filesep subj_name{cnt}]);
        close('all');
        %- manually assigning subject designation
        EEG_LS.subject = subj_name{cnt};
        %- save set
        [EEG_LS] = pop_saveset( EEG_LS, 'filepath', subj_save_ls{cnt}, 'filename', sprintf('%s_allTrials_LS.set',subj_name{cnt}));
        %- 
        ALLEEG_LS = [ALLEEG_LS; EEG_LS];
    end
end
%% MERGE EEG LIVEAMPS TO LOADSOL
mim_ls_merge_eeg()

%% CONVERT IMU TXT TO EVENT
%- change these variables
% imu_save_dir = ['/path/where/to/save/outputfiles'];
% subj_name = {'H1030','H2039','H3107'};
% files_imu = {'/path/to/H1030','/path/to/H2039','/path/to/H3107'};
imu_save_dir = [save_dir filesep 'IMU'];
if ~exist(imu_save_dir,'dir')
    mkdir(imu_save_dir);
end
%- don't change these variables (unless you want to).
ALLEEG_IMU = [];
SENSORS_TO_ANALYZE = {'Back'};
TRIALS_CHANGE = {'rest','TM_rest_1','TM_rest_2'};
TRIALS_NEW_NAMES = {'Rest','Rest','Rest'};
for cnt = 1:length(subj_name)
    [imu_table,rec_start_struct] =  convert_imu_to_table(files_imu{cnt});
    for f_i = 1:length(imu_table)
        %- required fields for rec_start_struct(f_i)
        rec_start_struct(f_i).subjectName = subj_name{cnt};
        [~,idx_end] = regexp(rec_start_struct(f_i).trialName,sprintf('%s_',subj_name{cnt}));
        trial_name = rec_start_struct(f_i).trialName; trial_name = trial_name(idx_end+1:end);
        if any(strcmp(trial_name,TRIALS_CHANGE))
            rec_start_struct(f_i).trialName = sprintf('%s_%s',subj_name{cnt},TRIALS_NEW_NAMES{strcmp(trial_name,TRIALS_CHANGE)});
        end
        % rec_start_struct(f_i).datetime = ['example_datetime'];
        rec_start_struct(f_i).filename = [subj_name{cnt} '_' rec_start_struct(f_i).trialName];
        rec_start_struct(f_i).filepath = subj_save_imu{cnt};
    end
    %- turn imu_table into EEGLAB compatible structure
    EEG_IMU = create_imu_struct(imu_table,rec_start_struct,SENSORS_TO_ANALYZE);
    %- save set
    % [EEG_IMU] = pop_saveset( EEG_IMU, 'filepath', rec_start_struct(1).filepath, 'filename', sprintf('%s_allTrials_IMU.set',subj_name{cnt}));
%     par_save(loadsol_events,subj_save_imu{cnt},sprintf('%s_LS_%i',subj_name{cnt},f_i));
    %## Convert Accelerometer Frame to Body Frame using Quanternions
    %-  NOTE: You'll need this toolbox to run this function: https://github.com/xioTechnologies/Gait-Tracking-With-x-IMU.git
    if ~exist([imu_save_dir filesep EEG_IMU.subject],'dir')
        mkdir([imu_save_dir filesep EEG_IMU.subject])
    end
    EEG_IMU = imu_get_body_frame(EEG_IMU, [imu_save_dir filesep EEG_IMU.subject],...
        'TRIAL_CROPPING_XLSX',[load_dir filesep '_studies' filesep 'subject_mgmt' filesep 'Trial_Cropping_V2_test.xlsx']);
    %-
    [EEG_IMU] = pop_saveset(EEG_IMU, 'filepath', subj_save_ls{cnt}, 'filename', sprintf('%s_allTrials_LS.set',subj_name{cnt}));
    %- 
    % ALLEEG_IMU = [ALLEEG_IMU; EEG_IMU];
end
%% USE IMU & Loadsol to get time locked gait events

%% HELPFUL CODE CELLS
%## TRANSFER LOADSOL DATA FROM R:\ TO M:\
%{
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data\';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Raw'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Raw'];
        transfer_folder(folder_from,folder_to,'*.txt');
    end
end
%}
%% TRANSFER IMU DATA FROM R:\ TO M:\
%{
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data\';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Raw'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Raw'];
        transfer_folder(folder_from,folder_to,'*.csv');
    end
end
%}
%% TRANSFER MERGED EEG DATA FROM R:\ TO M:\
%{
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data\';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Raw'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Raw'];
        transfer_folder(folder_from,folder_to,'*.csv');
    end
end
%}