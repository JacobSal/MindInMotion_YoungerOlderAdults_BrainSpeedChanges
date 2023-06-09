%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_test/1_paper_MIM/run_rerun_conn_mats.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
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
else  % isunix=
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep '4_conn_analysis' filesep 'AS'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
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
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ===================================================================== %%
%## Data PATHS
fprintf('\nData PATHS\n');
%## DATASET SPECIFIC
SUBJ_PICS = {{'02','03','04','05','09','11','15','16','18','19','21','22',...
            '23','24','25','27','28','29','30','31','32','33','35','36','38'}};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
%% (PARAMETERS) ======================================================== %%
fprintf('\nData Processing Parameters\n');
%## Hard Defines
%- dataset specific
DATA_SET = 'AS_dataset';
%- study group and saving
COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
%- connectivity process
CONN_FREQS = (1:100);
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
CNCTANL_TOOLBOX = 'sift'; %'bsmart'
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
DO_BOOTSTRAP = true;
DO_PHASE_RND = true;
%- datetime override
dt = '05252023_bounces_1h2h2bm_JS';
%## Soft Define
%- combinations of events and conditions
EVENT_COND_COMBOS = cell(length(COND_CHARS)*length(EVENT_CHARS),1);
cnt = 1;
for cond_i = 1:length(COND_CHARS)
    for event_i = 1:length(EVENT_CHARS)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',COND_CHARS{cond_i},EVENT_CHARS{event_i});
        cnt = cnt + 1;
    end
end
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
study_fName_2 = sprintf('%s_CONN_study',[EVENT_COND_COMBOS{:}]);
cluster_info_fpath = [STUDIES_DIR filesep 'as_cluster_info' filesep 'cluster_info.mat'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    %## (AS) ATTACH CLUSTER
%     MAIN_STUDY.cluster = MAIN_STUDY.urcluster;
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(MAIN_STUDY);
end
%% INITIALIZE PARFOR LOOP VARS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(MAIN_ALLEEG)]);
else
    POOL_SIZE = 1;
end
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
%- HARD DEFINES
%% PERFORM PERMUTATION STATISTICS
%{
STUDIES = cell(1,length(load_trials));
ALLEEGS = cell(1,length(load_trials));
%- Create STUDY & ALLEEG structs
for cond_i = 1:length(load_trials)
    if DO_UNIX
        study_fName = sprintf('%s_MIM_study_UNIX',load_trials{cond_i});
    else
        study_fName = sprintf('%s_MIM_study',load_trials{cond_i});
    end
    if ~exist([load_dir filesep study_fName '.study'],'file')
        error('ERROR. study file does not exist');
    else
        [STUDIES{cond_i},ALLEEGS{cond_i}] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);
    end
end
%## CONNECTIVITY STATISTICS RERUN 
%- Reconfigure
tmpALLEEGS = cell(1,length(ALLEEGS{1}));
for subj_i = 1:length(ALLEEGS{1})
    for cond_i = 1:length(load_trials)
        tmpALLEEGS{subj_i} = [tmpALLEEGS{subj_i} ALLEEGS{cond_i}(subj_i)];
    end
end
clear ALLEEGS
%- Generate Connectivity Statistics 
parfor subj_i = 1:length(tmpALLEEGS)
    %- Rerun statistics 
    fprintf('%s) Running wrapper_cnctanl_stats.m ...',tmpALLEEGS{subj_i}(1).subject);
    [tmpALLEEGS{subj_i}] = wrapper_cnctanl_stats(tmpALLEEGS{subj_i},...
        'DO_PHASERAND',DO_PHASERAND,...
        'DO_BOOTSTRAP',DO_BOOTSTRAP);
    %{
    for cond_i = 1:length(tmpALLEEGS{subj_i})
         %- load statistics
         EEG = tmpALLEEGS{subj_i}(cond_i);
         EEG = cnctanl_loadCAT(EEG,'NonzeroTest');
        for conn_i = 1:length(CONN_METHODS)
            %- calculate average connnectivity for each component pair across time
            [statMat,extract_sig] = gen_connMatrix(EEG,CONN_METHODS{conn_i},(1:length(EEG.CAT.Conn.freqs)),STAT_ALPHA);
            %- 
            EEG.etc.js_processing(cond_i).META.(CONN_METHODS{conn_i}).connExtract(1).freqs = EEG.CAT.Conn.freqs;
            EEG.etc.js_processing(cond_i).META.(CONN_METHODS{conn_i}).connExtract(1).mats = statMat;
            EEG.etc.js_processing(cond_i).META.(CONN_METHODS{conn_i}).connExtract(1).sigs = extract_sig;
            %##
        end
        [~,~] = pop_saveset(EEG,...
            'filepath',EEG.filepath,'filename',EEG.filename,...
            'savemode','twofiles');
    end
    %}
end
%}
%% HELPERS
%{
for cond_i = 1:length(load_trials)
    study_fName = sprintf('%s_MIM_study',load_trials{cond_i});
    [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);

    %## ROBUST SAVE STUDY
    %- fix paths from dferris to M
%     for subj_i = 1:length(MAIN_STUDY.datasetinfo)
%         MAIN_STUDY.datasetinfo(subj_i).filepath = convertPath2Drive(MAIN_ALLEEG(subj_i).filepath,'M');
%     end
%     [MAIN_STUDY,MAIN_ALLEEG] = pop_saavestudy( MAIN_STUDY, MAIN_ALLEEG,...
%                         'filename',study_fName,'filepath',save_dir);

    %- fix paths from M to dferris
    for subj_i = 1:length(MAIN_STUDY.datasetinfo)
        MAIN_STUDY.datasetinfo(subj_i).filepath = convertPath2UNIX(fullfile(MAIN_ALLEEG(subj_i).filepath),'dferris');
%         MAIN_ALLEEG(subj_i).filepath = convertPath2UNIX(MAIN_ALLEEG(subj_i).filepath,'dferris');
    end
    [MAIN_STUDY,MAIN_ALLEEG] = pop_savestudy( MAIN_STUDY, MAIN_ALLEEG,...
                        'filename',[study_fName '_UNIX'],...
                        'filepath',load_dir);
end
%}
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}