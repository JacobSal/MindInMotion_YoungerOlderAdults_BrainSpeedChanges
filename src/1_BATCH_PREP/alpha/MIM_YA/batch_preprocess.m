%   Project Title: Run BATCH processing for EEG filtering and artifact cleaning
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230223.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/alpha/MIM/run_batch_preprocess.sh

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
    DO_UNIX = false;
    PATH_EXT = 'M';
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '1_BATCH_PREP' filesep 'alpha' filesep 'MIM'];
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
%     eeg_options;
    % see. eeg_optionsbackup.m for all eeglab options.
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
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
%     pop_editoptions('option_parallel',0,'option_storedisk',1,...
%         'option_saveversion6',0,'option_cachesize',1000,'option_savetwofiles',1);
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
%## DATASET SPECIFIC
%- MIND IN MOTION
% SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
%     'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
%     'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H1008','H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1032','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047','H1048'}; % JACOB,SAL (04/18/2023)
SUBJ_2HMA = {'H2017', 'H2010', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2036', 'H2037', 'H2038',...
    'H2039', 'H2041', 'H2042', 'H2052', 'H2059', 'H2062', 'H2072', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3HMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120'}; % JACOB,SAL(02/23/2023)
SUBJ_3NHMA = {'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)

TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
%- Subject Picks
% SUBJ_PICS = {SUBJ_2HMA,SUBJ_3HMA,SUBJ_3NHMA}; % JACOB,SAL(02/23/2023)
% GROUP_NAMES = {'H2000''s','H3000''s','NH3000''s'}; % JACOB,SAL(02/23/2023)
% SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3HMA),1:length(SUBJ_3NHMA)}; % JACOB,SAL(02/23/2023)
SUBJ_PICS = {SUBJ_YNG}; % JACOB,SAL(04/18/2023)
GROUP_NAMES = {'H1000''s'}; % JACOB,SAL(04/18/2023)
SUBJ_ITERS = {1:length(SUBJ_YNG)}; % JACOB,SAL(04/18/2023)
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET]; 
%% ===================================================================== %%
%## PROCESSING PARAMS
%- datetime override
% dt = '24022023_OA_prep'
% dt = '07042023_OA_prep_verified';
dt = '04182023_YA_N37_prep_verified';
%## PATH & TEMP STUDY NAME
%- hard define
RUN_AMICA = false;
%- soft define
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## Store fNames and fPaths
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
cnt = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        if ~isempty(tmp)
            fNames{cnt} = tmp.name;
        else
            fNames{cnt} = [];
        end
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file')>1)
        stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
        cnt = cnt + 1;
    end
end
tmp = ~cellfun(@isempty,fNames);
fNames = fNames(tmp);
tmp = ~cellfun(@isempty,fPaths);
fPaths = fPaths(tmp);
tmp = ~cellfun(@isempty,subjectNames);
subjectNames = subjectNames(tmp);
%% SET POOLSIZE
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(fPaths)]);
else
    POOL_SIZE = 1;
end
%% (PARFOR) GENERATE CONNECTIVITY METRICS 
LOOP_VAR = 1:length(fPaths);
amica_cmd = cell(length(fPaths),1);
params = cell(length(fPaths),1);
parfor (subj_i = LOOP_VAR,POOL_SIZE)
% for subj_i = 1
    fprintf('Running subject %s...',subjectNames{subj_i})
    %## PREP for MAIN_FUNC
    if ~exist([save_dir filesep subjectNames{subj_i}],'dir')
        mkdir([save_dir filesep subjectNames{subj_i}]);
    end
    %## RUN MAIN_FUNC
    if ~exist([save_dir filesep subjectNames{subj_i} filesep 'clean' filesep sprintf('%s_cleanEEG.set',subjectNames{subj_i})],'file') || true
        [EEG,amica_cmd{subj_i},params{subj_i}] = main_func(subjectNames{subj_i},fPaths{subj_i},...
                [save_dir filesep subjectNames{subj_i}],STUDIES_DIR);
        fprintf('%s\n',amica_cmd{subj_i}{2})
    end
end
% if RUN_AMICA
%     for subj_i = LOOP_VAR
%         system([amica_cmd{subj_i}{2},''])
%     end
% end