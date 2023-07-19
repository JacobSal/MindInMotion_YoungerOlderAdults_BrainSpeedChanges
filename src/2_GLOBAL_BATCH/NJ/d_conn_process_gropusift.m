%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: s

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/NJ/run_d_conn_plotting.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'NJ'];
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
% addpath([submodules_dir filesep 'groupSIFT'])
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
%## PATHS
%- hardcode data_dir
DATA_SET = 'jacobsenN_dataset';
%- datetime override
% dt = '06292023_NJ_Standing';
dt = '07162023_NJ_standing_customheadmods';
%- connectiviy specific
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % AS (06/22/2023)
conn_meas = 'dDTF08';
%## soft define
%- combinations of events and conditions
TRIAL_TYPES = {'pre','post'};
%- path for local data
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'conn'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(TRIAL_TYPES)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    %## load chang's algorithmic clustering
    %* cluster parameters
    pick_cluster = 8;
    clustering_weights.dipoles = 1;
    clustering_weights.scalp = 0;
    clustering_weights.ersp = 0;
    clustering_weights.spec = 0;
    cluster_alg = 'kmeans';
    do_multivariate_data = 1;
    evaluate_method = 'min_rv';
    clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
        '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
        '_spec_',num2str(clustering_weights.spec)];
    %* load cluster information
    cluster_load_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
    outputdir = [cluster_load_dir filesep clustering_method,...
        filesep num2str(pick_cluster)];
    tmp = load([outputdir filesep sprintf('cluster_update_%i.mat',pick_cluster)]);
    cluster_update = tmp.cluster_update;
    STUDY.cluster = cluster_update;
    %- get inds
    [comps_out,main_cl_inds,outlier_cl_inds,valid_clusters] = eeglab_get_cluster_comps(STUDY);
end
%% MAKE CONDITIONS SUBJ STACKS
parfor subj_i = 1:length(ALLEEG)
    EEG = ALLEEG(subj_i);
    [EEG,stats_bootstrap,stats_nonzero] = cnctanl_bootstrap_test(EEG);
    for cond_i = 1:length(EEG)
        pop_saveset(EEG(cond_i),'filename',EEG(cond_i).filename,'filepath',[save_dir filesep 'groupsift']);
    end
end
%% GROUPSIFT STEPS
%## eegplugin_groupSIFT.m
% uimenu( submenu, 'label', '1.Run SIFT batch',                     'callback', 'pop_groupSIFT_runSiftBatch');
% uimenu( submenu, 'label', '2.Validate AR models',                 'callback', 'pop_groupSIFT_validateArModels');
% uimenu( submenu, 'label', '3.Convert to group anatomical ROIs',   'callback', 'pop_groupSIFT_convertToGroupAnatomicalRois');
% uimenu( submenu, 'label', '4.Compute t-stats & p-values',         'callback', 'pop_groupSIFT_computeTstatsAndPvalues');
% uimenu( submenu, 'label', '5 Show pre-selected ROIs',             'callback', 'pop_groupSIFT_showPreselectedRois');
% uimenu( submenu, 'label', '6.View results & Export for movie',    'callback', 'pop_groupSIFT_viewResultsAndExportForMovie');
%## 1) generate connectivity values
%## 2) validate AR models
pop_groupSIFT_validateArModels()
%## 3) cluster and fit to anatomical ROIs
pop_groupSIFT_convertToGroupAnatomicalRois()
%## 4) Compute t-stats & p-values (bootstrapped?)
pop_groupSIFT_computeTstatsAndPvalues()
%## 5) Display ROIS
% pop_groupSIFT_showPreselectedRois
%## 6) View Connectivity Results
% pop_groupSIFT_viewResultsAndExportForMovie

