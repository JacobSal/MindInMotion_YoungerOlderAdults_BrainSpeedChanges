%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/0_DATA_PREP/2_RAW_DATA_VIS/MIM_OA/run_process_ls_imu.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave','software') % might be needed to plot dipole plots?
%## TIME
tic
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '0_DATA_PREP' filesep '2_RAW_DATA_VIS' filesep 'MIM_OA'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
addpath([PATH_ROOT filesep REPO_NAME filesep 'submodules' filesep 'Gait Tracking With x-IMU' filesep 'Quaternions']);
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
colormap(linspecer);
close all;
% dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
dt = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%- Subject Directory information
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [DATA_DIR filesep DATA_SET]; %[STUDIES_DIR filesep sprintf('%s',dt) filesep 'ls_imu_proc'];
TRIAL_CROPPING_XLSX = [DATA_DIR filesep DATA_SET filesep '_studies' filesep 'subject_mgmt' filesep 'Trial_Cropping_V2_test.xlsx'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ASSIGN FILE PATHS
M_MIND_IN_MOTION_DIR    = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
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
        fprintf('loading subject %s...\n',SUBJ_PICS{group_i}{subj_i})
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
%% CONVERT LOADSOL TXT TO EVENTS
ls_save_dir = [save_dir filesep '_figs' filesep 'LS'];
if ~exist(ls_save_dir,'dir')
    mkdir(ls_save_dir);
end
EEG_LS = [];
ALLEEG_LS = cell(length(subj_name),1);
% all_subj_inds = 1:length(subj_name);
parfor (cnt = 1:length(subj_name),length(subj_name))
    [loadsol_table, rec_start_struct] = convert_loadsol_to_table(files_loadsol{cnt});
    if ~exist([ls_save_dir filesep subj_name{cnt}],'dir')
        mkdir([ls_save_dir filesep subj_name{cnt}]);
    end
    for f_i = 1:length(loadsol_table)
        TMP = create_loadsol_struct(loadsol_table{f_i},rec_start_struct(f_i),[ls_save_dir filesep subj_name{cnt}]);
        if ~exist([ls_save_dir filesep subj_name{cnt}],'dir')
            mkdir([ls_save_dir filesep subj_name{cnt}])
        end
        %- manually assigning subject designation
        TMP.subject = subj_name{cnt};
        %- save set
        [TMP] = pop_saveset( TMP, 'filepath', subj_save_ls{cnt}, 'filename', sprintf('%s_allTrials_LS.set',subj_name{cnt}));
        %- 
        EEG_LS = [EEG_LS; TMP];
    end
%     ALLEEG_LS{cnt} = EEG_LS;
end
%% CONVERT IMU TXT TO EVENT
imu_save_dir = [save_dir filesep '_figs' filesep 'IMU'];
if ~exist(imu_save_dir,'dir')
    mkdir(imu_save_dir);
end
EEG_IMU = [];
ALLEEG_IMU = cell(length(subj_name),1);
SENSORS_TO_ANALYZE = {'Back'};
TRIALS_CHANGE = {'rest','SP_1p0_1'};
TRIALS_NEW_NAMES = {'Rest','SP_1p0_1'};
parfor (cnt = 1:length(subj_name),length(subj_name)/2)
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
%     [EEG_IMU] = pop_saveset( EEG_IMU, 'filepath', rec_start_struct(1).filepath, 'filename', sprintf('%s_allTrials_IMU.set',subj_name{cnt}));
%     par_save(loadsol_events,subj_save_imu{cnt},sprintf('%s_LS_%i',subj_name{cnt},f_i));
    %## Convert Accelerometer Frame to Body Frame using Quanternions
    %-  NOTE: You'll need this toolbox to run this function: https://github.com/xioTechnologies/Gait-Tracking-With-x-IMU.git
    if ~exist([imu_save_dir filesep EEG_IMU.subject],'dir')
        mkdir([imu_save_dir filesep EEG_IMU.subject])
    end
    EEG_IMU = imu_get_body_frame(EEG_IMU, [imu_save_dir filesep EEG_IMU.subject],...
        'TRIAL_CROPPING_XLSX',TRIAL_CROPPING_XLSX);
    [EEG_IMU] = pop_saveset( EEG_IMU, 'filepath', rec_start_struct(1).filepath, 'filename', sprintf('%s_allTrials_IMU.set',subj_name{cnt}));
    %-
%     ALLEEG_IMU{cnt} = EEG_IMU;
end
%% USE IMU & Loadsol to get time locked gait events
% TRIAL_IND = 1;
% ACC_INDS = [1,2,3]; %[2,3,4];
% for cnt = 1:length([SUBJ_ITERS{:}])
%     EEG_IMU = pop_loadset('filepath', subj_save_imu{cnt}, 'filename', sprintf('%s_allTrials_IMU.set',subj_name{cnt}));
%     EEG_LS = pop_loadset('filepath', subj_save_ls{cnt}, 'filename', sprintf('%s_allTrials_LS.set',subj_name{cnt}));
% %     timewarp_imu_ls();
% %     EEG_IMU=eeg_checkset(EEG_IMU,'loaddata');
%     for TRIAL_IND = 1:length(EEG_IMU.event)-1
%         trial_name = EEG_IMU.event(TRIAL_IND).trialName;
%         trial_idxs = [EEG_IMU.event(TRIAL_IND).latency, EEG_IMU.event(TRIAL_IND+1).latency];
%         trial_idxs = [trial_idxs(1)+500,trial_idxs(1)+1000];
%         pos_in = EEG_IMU.data(ACC_INDS,trial_idxs(1):trial_idxs(2))';
%         quat_in = quatern2rotMat(quaternConj(EEG_IMU.etc.sensor_frame_ori(trial_idxs(1):trial_idxs(2),:)));
%         fprintf('making video for subject %s and trial %s\n',subj_name{cnt},trial_name);
%         disp_imu_movie(pos_in,quat_in,[imu_save_dir filesep EEG_IMU.subject filesep sprintf('%i_6DOFmovie',TRIAL_IND)]);
%     end
% end
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
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
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