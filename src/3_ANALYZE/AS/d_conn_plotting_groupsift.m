%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: s

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/AS/run_d_conn_plotting_groupsift.sh

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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'AS'];
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
addpath([submodules_dir filesep 'groupSIFT'])
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
%## HARDCODE
CLUSTER_ITERS = [2,3,4,5,6,7,8,9,10,11];
CLUSTER_ASSIGNMENTS = {'RPPa','Cuneus','Precuneus','RSuppMotor','LPPa','LSM','RSM','LTemp','Cing','LSuppMotor'}; 
% (06/27/2023) JS, unsure on these as of yet.
CONN_MEAS_ANLYZ = 'dDTF08';
%- hardcode data_dir
DATA_SET = 'AS_dataset';
COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
ALPHA = 0.05;
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
dt = '08182023_bounces_1h_2h_2bm_JS';
%- connectiviy specific
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % AS (06/22/2023)
%## soft define
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
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
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
    POOL_SIZE = min([SLURM_POOL_SIZE,length(EVENT_COND_COMBOS)*4]);
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
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
%## VARIABLE INITIALIZATION
tmp = ALLEEG(1).etc.conn_table;
% CONN_METHODS = unique(tmp.t_conn_meas);
conn_conds = unique(tmp.t_fNames);
conn_comps = [tmp.t_conn_comps{1}];
idx = find(strcmp(CONN_MEAS_ANLYZ,tmp.t_conn_meas));
%- get unique frequencies tested
conn_freqs = tmp.t_conn_freqs(idx);
tmpsz = cellfun(@size,conn_freqs,'UniformOutput',false); tmpsz = cellfun(@max,conn_freqs); tmpsz = max(tmpsz);
store_freqs = zeros(length(conn_freqs),tmpsz);
for i = 1:length(conn_freqs)
    store_freqs(i,1:length(conn_freqs{i})) = conn_freqs{i};
end
uniq_freqs = unique(store_freqs,'rows');
%% REMOVE COMPS
%- select components from EEG
groupsift_fpath = [load_dir filesep '_groupsift'];
eeg_dir = dir([groupsift_fpath filesep '*.set']);
for f_i = 1:length(eeg_dir)
    EEG = pop_loadset('filename',eeg_dir(f_i).name,'filepath',eeg_dir(f_i).folder);
    EEG.etc.urreject.ic_keep = 1:size(EEG.icaweights);
    if length(EEG.etc.urreject.ic_keep) < 2 || isempty(EEG.etc.urreject)
        fprintf('** Subject %s rejected.\n',EEG.subject);
%         tmp_rmv_subjs(subj_i) = 1;
    else
        EEG.icachansind = EEG.etc.urreject.ic_keep;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
%             EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
        end
        ica_weights = (EEG.icaact/EEG.data(EEG.icachansind,:))/EEG.icasphere;
        ica_winv = (EEG.data(EEG.icachansind,:)/EEG.icaact);
        tmp = sum(sqrt((EEG.icaweights-ica_weights).^2),[1,2]);
        fprintf('%7ssum(sqrt((icaweights_new-icaweights_old).^2)) = %0.3f\n','',tmp);
        %## Update ALLEEG & comps_out
        EEG.icaweights = ica_weights;
        EEG.icawinv = ica_winv;
    end
    EEG = pop_saveset(EEG,'filename',eeg_dir(f_i).name,'filepath',eeg_dir(f_i).folder);
end

%% GROUPSIFT STEPS
%## eegplugin_groupSIFT.m
pop_editoptions('option_single', 1);
% uimenu( submenu, 'label', '1.Run SIFT batch',                     'callback', 'pop_groupSIFT_runSiftBatch');
% uimenu( submenu, 'label', '2.Validate AR models',                 'callback', 'pop_groupSIFT_validateArModels');
% uimenu( submenu, 'label', '3.Convert to group anatomical ROIs',   'callback', 'pop_groupSIFT_convertToGroupAnatomicalRois');
% uimenu( submenu, 'label', '4.Compute t-stats & p-values',         'callback', 'pop_groupSIFT_computeTstatsAndPvalues');
% uimenu( submenu, 'label', '5 Show pre-selected ROIs',             'callback', 'pop_groupSIFT_showPreselectedRois');
% uimenu( submenu, 'label', '6.View results & Export for movie',    'callback', 'pop_groupSIFT_viewResultsAndExportForMovie');
%## 1) generate connectivity values
pop_groupSIFT_runSiftBatch()
%## 2) validate AR models
pop_groupSIFT_validateArModels()
%## 3) cluster and fit to anatomical ROIs
pop_groupSIFT_convertToGroupAnatomicalRois()
%## 4) Compute t-stats & p-values (bootstrapped?)
pop_groupSIFT_computeTstatsAndPvalues()
%## 5) Display ROIS
% pop_groupSIFT_showPreselectedRois()
%## 6) View Connectivity Results
% pop_groupSIFT_viewResultsAndExportForMovie()
