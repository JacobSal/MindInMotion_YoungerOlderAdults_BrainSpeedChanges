%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_OA_proc/run_epoch_process.sh

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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'gamma' filesep 'NJ_proc'];
%- addpath for _functions folder
addpath(run_dir)
addpath(source_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP
if ~ispc
%     eeg_options;
    % see. eeg_optionsbackup.m for all eeglab options.
%     pop_editoptions('option_parallel',0,'option_storedisk',1,...
%         'option_saveversion6',0,'option_cachesize',1000,'option_savetwofiles',1);
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
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS && VARIABLES
% cnctanl_config = struct('var_name',[],...
%                         'var_value',[])
%- hardcode data_dir
DATA_SET = 'jacobsenN_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
SUBJ_RISERS = {'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34',...
                     'S35','S36','S37','S38','S39','S40','S41','S42','S43','S44',...
                     'S46','S47','S48'};
SUBJ_FALLERS = {};
SUBJ_PICS = {SUBJ_RISERS,SUBJ_FALLERS};
SUBJ_ITERS = {1:length(SUBJ_RISERS),1:length(SUBJ_FALLERS)};
GROUP_NAMES = {'risers','fallers'};
%% ===================================================================== %%
%## PROCESSING PARAMS
% pop_editoptions('option_parallel',1);
params = [];
%- study group and saving
params.OVERRIDE_DIPFIT =  true;
% SAVE_EEG = false; % saves EEG structures throughout processing
%- component rejection crit
params.THRESHOLD_DIPFIT_RV = 0.15;
params.THRESHOLD_BRAIN_ICLABEL = 0.50;
%- MIM specific epoching
params.EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
%- epoching params
params.TRIAL_TYPES = {'pre','post'};
params.EPOCH_TIME_LIMITS = [-3,3]; % [-1,3] captures gait events well
params.TRIAL_LENGTH = 0; % trial length in seconds
params.PER_OVERLAP = 0.0; % percent overlap between epochs
params.TRIAL_BEGIN_STR = '';
params.TRIAL_END_STR = '';
params.EVENT_TRIAL_PARSER = '';
params.EVENT_COND_PARSER = '';
%- connectivity process
params.CONN_FREQS = (1:100);
params.CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
%- subj_i_genConnMeas
params.CNCTANL_TOOLBOX = 'sift'; %'bsmart'
params.WINDOW_LENGTH = 0.5;
params.WINDOW_STEP_SIZE = 0.025;
params.NEW_SAMPLE_RATE = [];
params.DO_BOOTSTRAP = true;
% ASSIGN_BOOTSTRAP_MEAN = false;
% SAVE_CONN_BOOTSTRAP = false;
%- subj_i_genConnStats
params.DO_PHASE_RND = true;
%- eeglab_cluster.m spectral params
params.FREQ_LIMITS = [1,100];
params.CYCLE_LIMITS = [3,0.8];
params.SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
params.FREQ_FAC = 4;
params.PAD_RATIO = 2;
%% global script chain (VERSION 1)
%- datetime override
dt = '05082023_NJ_Standing';
%## PATH & TEMP STUDY NAME
%- hard define
DO_CONN_ANL = false;
SAVE_EEG = true;
SESSION_NUMBER = '1';
study_fName_3 = sprintf('%s_EPOCH_study',[params.TRIAL_TYPES{:}]);
study_fName_4 = sprintf('%s_CONN_study',[params.TRIAL_TYPES{:}]);
%- soft define
path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
mniMRI = fullfile(path2BEM, 'standard_mri.mat');
mniVol = fullfile(path2BEM, 'standard_vol.mat');
mniChan1005 = fullfile(path2BEM,'elec','standard_1005.elc');
TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_3 '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_3 '_UNIX.study'],'filepath',load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_3 '.study'],'filepath',load_dir);
    end
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
SUFFIX_PATH_EPOCHED = 'EPOCHED';
%- clear vars for memory
% clear MAIN_ALLEEG
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
%     tmp = cell(length(MAIN_ALLEEG),1);
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),POOL_SIZE)
%     for subj_i = LOOP_VAR
    if ~isempty(fPaths_out{subj_i})
        %- Parse out components
        components = comps_out(:,subj_i);
        components = sort(components(components ~= 0));
        %## LOAD EEG DATA
        EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        try
            %## RUN MAIN_FUNC
            [TMP_EEG] = main_func_v2(EEG,save_dir,components,...
                'SAVE_EEG',SAVE_EEG,...
                'CONN_METHODS',params.CONN_METHODS,... %
                'FREQS',params.CONN_FREQS,...
                'CNCTANL_TOOLBOX',params.CNCTANL_TOOLBOX,... 
                'DO_BOOTSTRAP',params.DO_BOOTSTRAP,...
                'DO_PHASE_RND',params.DO_PHASE_RND,...
                'WINDOW_LENGTH',params.WINDOW_LENGTH,...
                'WINDOW_STEP_SIZE',params.WINDOW_STEP_SIZE);
%                 pop_mergeset();
            EEG.CAT = TMP_EEG.CAT;
            [EEG] = pop_saveset(EEG,...
                'filepath',EEG.filepath,'filename',EEG.filename,...
                'savemode','twofiles');
            tmp{subj_i} = EEG;
        catch e
            rmv_subj(subj_i) = 1;
            EEG.CAT = struct([]);
            tmp{subj_i} = EEG;
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
        end
    end
end
pop_editoptions('option_computeica', 0);
%% SAVE BIG STUDY
MAIN_STUDY.etc.pipe_params = params;
[ALLEEG,MAIN_STUDY] = parfunc_rmv_subjs(tmp,MAIN_STUDY,rmv_subj);
%- Save
[MAIN_STUDY,ALLEEG] = parfunc_save_study(MAIN_STUDY,ALLEEG,...
                                        study_fName_4,save_dir,...
                                        'STUDY_COND',[]);
%% Version History
%{
v2.0; (04/28/2023) JS: Splitting up the epoching, plotting, and
connectivity calculation process to avoid bugs with STUDY files.
v1.0; (11/11/2022), JS: really need to consider updating bootstrap
    algorithm with parallel computing. Taking ~ 1 day per
    condition for all subjects and the bottle neck is entirely the
    bootstrap.

    Note: validateattributes and assert functions may be helpful
    in more clearly defining function inputs.
        e.g.  DO_PHASE_RND = true;
          errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
          validationFcn = @(x) assert(islogical(x),errorMsg);
v1.0; (12/5/2022) Need to adapt this to include all conditions
    within each SUBJ structure so connectivity can be calculated
    for the ALLEEG structure rather than the EEG structure.
    *** maybe try to ditch the SUBJ strucutre entirely for this
    round?
v1.0.01132023.0 : Initializing versioning for future iterations.
%}

