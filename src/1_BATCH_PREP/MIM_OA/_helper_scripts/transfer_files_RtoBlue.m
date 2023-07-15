%   Project Title: Transfer major project files from long term storage
%   drive to hypercomputer storage
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230223.0
%   Previous Version: n/a
%   Summary: 
%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%## TIMEd
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
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
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '1_BATCH_PREP' filesep 'MIM_OA_YA'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = true;
setWorkspace
%% PARPOOL SETUP
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
    'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## (MIND IN MOTION) DATASET SPECIFIC (05/24/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H1008','H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
SUBJ_NO_MRI = {'H2010', 'H2036', 'H2041', 'H2072', 'H3018','H3120'};
SUBJ_1YA = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1032','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047','H1047'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2017', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2037', 'H2038',...
    'H2039', 'H2042', 'H2052', 'H2059', 'H2062', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- (OY) Subject Picks 
% SUBJ_PICS = {SUBJ_1YA}; 
% GROUP_NAMES = {'H1000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_1YA)}; 
%- (OA) Subject Picks 
SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)}; 
%% (PROCESSING PARAMS) ================================================= %%
%## hard define
%- dataset name
DATA_SET = 'MIM_dataset';
%- datetime override
OA_PREP_FNAME = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
%## soft define
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FNAME]; % JACOB,SAL(02/23/2023)
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
save_dir = [STUDIES_DIR filesep sprintf('%s',OA_PREP_FNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (FIND FILES) ======================================================== %%
%## (HELPER SCRIPT) TRANSFER ELEC_ALIGNED.MAT & VOL.MAT & VALIDATION FIGURES FROM M:\ TO R:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_to = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01'];
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'vol.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_mesh.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_sens_1.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_sens_2.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_sourceplot.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
    end
end
%% (HELPER SCRIPT) TRANSFER HEADMODEL & ELECTRODE_ALIGNED DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'HeadModel'];
        if exist(folder_to,'dir')
            rmdir(folder_to,'s')
        end
        %- custom electrode locations & headmodel for SKULL_0.01
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'HeadModel' filesep 'skull_0p01'];
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned_init.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        %- custom electrode locations and headmodels for SKULL_0.0042
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'HeadModel' filesep 'skull_0p0042'];
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned_init.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'vol.mat'];
        if exist(file_from,'file')
%             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
    end
end
%% (HELPER SCRIPT) TRANSFER MRI & Fiducial DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'MRI'];
%         if exist(folder_to,'dir')
%             rmdir(folder_to,'s')
%         end
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' ];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'HeadModel'];
%         delete(folder_to)
%         if exist(folder_to,'dir')
%             rmdir(folder_to)
%         end
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        %- custom electrode locations for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',SUBJ_PICS{group_i}{subj_i},file_from)
        end
        %- custom headmodel for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',SUBJ_PICS{group_i}{subj_i},file_from)
        end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Raw' filesep ''];
%         copyfile(file_from,folder_to);
        %- custom electrode locations from digital 3D scan
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.txt'];    
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc rostral-caudal reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc_rs.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to ctf reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'ctf_fiducials.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri segmentation using SimNIBS aligned to acpc/s coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep sprintf('%s_MRI_acpc_rs.nii',SUBJ_PICS{group_i}{subj_i})];
%         copyfile(tmp,folder_to);
        %-----------------------------------------------------------------%
        %{
        %- electrode coordinates from HEADSCAN aligned to MRI mesh
        tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'elec_aligned.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
                continue;
            end
        else
            fprintf('copying file %s.\n',tmp);
            copyfile(tmp,folder_to);
        end
        
        %- subject specific headmodel volume 
        tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'vol.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
                continue;
            end
        else
            fprintf('copying file %s.\n',tmp);
            copyfile(tmp,folder_to);
        end
        %}
    end
end
%% (HELPER SCRIPT) TRANSFER CHANLOCS DATA FROM R:\ TO M:\
BAD_SUBJS = {}; %{'NH3004','NH3009'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan'];
        delete(folder_to)
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        copyfile(file_from,folder_to);
    end
end
%% (HELPER SCRIPT) TRANSFER FILES FROM TEH AMICA FOLDER TO THE CLEAN FOLDER FOR ALLS
%SUBJECTS
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];
subjs = [SUBJ_PICS{:}];
% dt = '07042023_OA_prep_verified';
% dt = '04182023_YA_N37_prep_verified';
dt = '05192023_YAN33_OAN79_prep_verified';
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
for subj_i = 1:length(subjs)
    folder_from = [save_dir filesep subjs{subj_i} filesep 'amica'];
    folder_to = [save_dir filesep subjs{subj_i} filesep 'clean'];
    if ~exist([folder_to filesep 'W'],'file')
        fprintf('Transfering folder for subject %s...\n',subjs{subj_i});
        transfer_folder(folder_from,folder_to,'.');
    end
end
%% (HELPER SCRIPT) TRANSFER MERGED EEG DATA FROM R:\ TO M:\
BAD_SUBJS = {'NH3004','NH3009'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data\';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        % R:\Ferris-Lab\share\MindInMotion\Data\H3039\EEG\Trials
        folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
%         folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Merged'];
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Merged'];
%         transfer_folder(folder_from,folder_to,'*EEGandLSandIMU.set');
%         transfer_folder(folder_from,folder_to,'*EEGandLSandIMU.fdt');
        tmp = dir([folder_to filesep '*.set']);
        if (isempty(tmp) && exist(folder_from,'dir')) && ~any(strcmp(SUBJ_PICS{group_i}{subj_i},BAD_SUBJS))
            transfer_folder(folder_from,folder_to,'*.set');
            transfer_folder(folder_from,folder_to,'*.fdt');
        end
    end
end
%% HELPER SCRIPT
EMAIL_CHAR = 'jsalminen@ufl.edu';
LOOP_VAR = 1:length(fPaths);
avg_ref_pca_reduction = 1; % length({'EEG'})
for subj_i = LOOP_VAR
    %## RUN MAIN_FUNC
    amica_out_fPath = [save_dir filesep subjectNames{subj_i} filesep 'amica'];
    float_fPath = [save_dir filesep subjectNames{subj_i} filesep 'clean'];
    set_fName = sprintf('%s_cleanEEG.set',subjectNames{subj_i});
    if ~exist(amica_out_fPath,'dir')
        mkdir(amica_out_fPath);
    end
    if exist(float_fPath,'dir')
        %- find .set file
        tmp = dir([float_fPath filesep '*.set']);
        %- load EEG
        EEG = pop_loadset('filepath',float_fPath,'filename',tmp.name);
        %- set the .fdt (float) file
        tmp = split(tmp.name,'.');
        float_fName = [tmp{1} '.fdt'];
        %- create bash and param files
        [EEG,cmd_out] = mim_prep_hpg_amica(EEG,[float_fPath filesep float_fName],float_fPath,EMAIL_CHAR,avg_ref_pca_reduction);
        fprintf('%s\n',cmd_out{2});
    end
    if ~ispc
        system([cmd_out,''])
    else
        fprintf('run this on unix\n');
    end
end