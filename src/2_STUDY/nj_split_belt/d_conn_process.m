%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/NJ/run_d_conn_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;s
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'NJ'];
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
%% (PARAMETERS) ======================================================== %%
%## Hard Defines
%- datset name
DATA_SET = 'jacobsenN_dataset';
%- datetime override
% dt = '07162023_NJ_standing_customheadmods';
dt = '09082023_NJ_gait_customheadmods';
%- epoching params
EPOCH_METHOD = 'sliding_window_gait';
%* sliding window
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
%* gait
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-0.25,3]; %[-1,3]; %[-0.5,5]; % [-1,3] captures gait events well , [-0.5,5] captures gait events poorly
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
switch EPOCH_METHOD
    case 'sliding_window_standing'
        SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
        TRIAL_TYPES = {'pre','post'};
    case 'sliding_window_gait'
        SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED_GAIT';
        TRIAL_TYPES = {'B1','B2','B3','P1','P2','SB1','SB2','none'};
    case 'gait'
        SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
        TRIAL_TYPES = {'B1','B2','B3','P1','P2','SB1','SB2','none'};
    otherwise
        error('error. choose a valid EPOCH_METHOD');        
end
%- connecitivty modeling
CONN_FREQS = (1:200);
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
CNCTANL_TOOLBOX = 'sift'; %'bsmart'
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
%- connectivity statistics 
% NEW_SAMPLE_RATE = [];
DO_BOOTSTRAP = true;
DO_PHASE_RND = true;
%## Soft Defines
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
study_save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
study_load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
conn_save_dir = [study_save_dir filesep '_figs' filesep 'conn'];
%- create new study directory
if ~exist(study_save_dir,'dir')
    mkdir(study_save_dir);
end
if ~exist(conn_save_dir,'dir')
    mkdir(conn_save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([study_load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file %s does not exist',[study_load_dir filesep study_fName_1 '.study']);
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',study_load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',study_load_dir);
    end
    %## load chang's algorithmic clustering
    %* cluster parameters
%     pick_cluster = 8;
%     clustering_weights.dipoles = 1;
%     clustering_weights.scalp = 0;
%     clustering_weights.ersp = 0;
%     clustering_weights.spec = 0;
%     cluster_alg = 'kmeans';
%     do_multivariate_data = 1;
%     evaluate_method = 'min_rv';
%     clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
%         '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
%         '_spec_',num2str(clustering_weights.spec)];
%     %* load cluster information
%     cluster_load_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
%     outputdir = [cluster_load_dir filesep clustering_method,...
%         filesep num2str(pick_cluster)];
%     tmp = load([outputdir filesep sprintf('cluster_update_%i.mat',pick_cluster)]);
%     cluster_update = tmp.cluster_update;
%     MAIN_STUDY.cluster = cluster_update;
    %- get inds
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
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),POOL_SIZE)
% for subj_i = LOOP_VAR
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
        EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
    end
    fprintf('%s) Processing componets:\n',EEG.subject)
    fprintf('%i,',components'); fprintf('\n');
    %- re-epoch
    ALLEEG = cell(1,length(TRIAL_TYPES));
    for i = 1:length(EEG.etc.cond_files)
        if ~ispc
            ALLEEG{i} = pop_loadset('filepath',convertPath2UNIX(EEG.etc.cond_files(i).fPath),'filename',EEG.etc.cond_files(i).fName);
        else
            ALLEEG{i} = pop_loadset('filepath',convertPath2Drive(EEG.etc.cond_files(i).fPath),'filename',EEG.etc.cond_files(i).fName);
        end
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    try
        %## RUN MAIN_FUNC
        [TMP,t_out] = main_func(ALLEEG,conn_save_dir,...
            'CONN_COMPONENTS',components,...
            'CONN_METHODS',CONN_METHODS,... 
            'FREQS',CONN_FREQS,...
            'CNCTANL_TOOLBOX',CNCTANL_TOOLBOX,... 
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE);
%         for i=1:length(TMP)
%             par_save(TMP(i).CAT)
%         end
        BIG_CAT = cat(1,TMP(:).CAT);
        EEG.etc.COND_CAT = BIG_CAT;
        EEG.etc.conn_table = t_out;
        fName = strsplit(EEG.filename,'.'); fName = [fName{1} '.mat'];
        par_save(t_out,EEG.filepath,fName,'_conntable');
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
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
        disp(e);
    end
end
pop_editoptions('option_computeica',0);
%% SAVE BIG STUDY
% [ALLEEG,MAIN_STUDY] = parfunc_rmv_subjs(tmp,MAIN_STUDY,rmv_subj);
%- Save
[MAIN_STUDY,tmp] = parfunc_save_study(MAIN_STUDY,tmp,...
                                        study_fName_1,study_save_dir,...
                                        'STUDY_COND',[]);
%% Version History
%{
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
