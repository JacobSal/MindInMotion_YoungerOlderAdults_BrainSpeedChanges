%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: This code was modified from

%- run sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/AS_proc/run_epoch_process.sh

% sbatch /blue/dferris/bishoy.pramanik/GitHub/Jacob/src/AS_PingPong_Connectivity/AS_proc_v2/run_conn_process.sh

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
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'gamma' filesep 'AS_proc'];
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
%- Hardcodes directory (data_dir) to .SET files for raw EEG data
DATA_SET = 'AS_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
SUBJ_PICS = {{'02','03','04','05','09','11','15','16','18','19','21','22',...
            '23','24','25','27','28','29','30','31','32','33','35','36','38'}};
% SUBJ_ITERS = {[1,5]};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
%COND_CHARS = {'cooperative','competitive'};
FNAME_POSTFIX = 'iCanNoiseEMG_postAMICA_postDIPFIT_cat12_AutoSelectComponents_ManuallySelectComponents_GoProCheckedEvents_Epoched_ALL_NoWarp_ChangeEpochRejection_Downsample_v2';
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET];

%% Data Processing Parameters
fprintf('\nData Processing Parameters\n');
% pop_editoptions('option_parallel',1);
%- study group and saving
SAVE_EEG = false; % saves EEG structures throughout processing
%-
COND_FIELD_PARSER = 'bounces'; %'condlabel';
EVENT_FIELD_PARSER = 'type';
% COND_CHARS = {'human','BM'};
COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
EVENT_COND_COMBOS = cell(length(COND_CHARS)*length(EVENT_CHARS),1);
cnt = 1;
for cond_i = 1:length(COND_CHARS)
    for event_i = 1:length(EVENT_CHARS)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',COND_CHARS{cond_i},EVENT_CHARS{event_i});
        cnt = cnt + 1;
    end
end
EPOCH_TIME_LIMITS = [-1,1];
PARSE_TYPE = 'Custom';
%- eeglab_cluster.m spectral params
FREQ_LIMITS = [1,100];
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
PAD_RATIO = 2;
%% global script chain (VERSION 1)
%- datetime override
%dt = '16032023_OA_subset'; % OA
% dt = '03292023_JS_test_v2'; 
% dt = '05082023_ASdataset';
dt = '05222023_run_JS';
%## PATH & TEMP STUDY NAME
fprintf('\nPATH and TEMP .study Names\n');
%## PARAMS
%- hard define
SAVE_EEG = true;
SESSION_NUMBER = '1';
SUFFIX_PATH_EPOCHED = 'EPOCHED';
study_fName_1 = sprintf('%s_all_comps_study',[EVENT_COND_COMBOS{:}]);
study_fName_2 = sprintf('%s_reduced_comps_study',[EVENT_COND_COMBOS{:}]);
study_fName_3 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
%- soft define
cluster_info_fpath = [STUDIES_DIR filesep 'as_cluster_info' filesep 'cluster_info.mat'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Determining paths and file names
fprintf('\nDetermining paths and file names\n');
% - Storing fNames and fPaths for saving data
conditions      = cell(1,length([SUBJ_ITERS{:}]));
groups          = cell(1,length([SUBJ_ITERS{:}]));
sessions        = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fPath   = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fName   = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;

for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}];
        fNames{cnt} = ['Pilot'  SUBJ_PICS{group_i}{subj_i} '_' FNAME_POSTFIX '.set'];
        fprintf('Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file'))
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;

    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = ['Pilot' SUBJ_PICS{group_i}{subj_i}];
        tmp = strjoin(COND_CHARS,'_'); 
        conditions{cnt} = tmp; %tmp{:};
        groups{cnt} = num2str(group_i);
        sessions{cnt} = '1';
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%% CREATE STUDY
fprintf('\nCreating .study files\n');
%- NOTE: need a dipole fit before this step can be done.
%- params
%## Create STUDY & ALLEEG structs
if ~exist([save_dir filesep study_fName_1 '.study'],'file') %|| true
    fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
    [MAIN_ALLEEG] = as_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                        conditions,groups,sessions);
%     [MAIN_STUDY,MAIN_ALLEEG] = mim_create_study(MAIN_ALLEEG,study_fName_1,save_dir);
    [MAIN_STUDY,MAIN_ALLEEG] = as_create_study(MAIN_ALLEEG,cluster_info_fpath,study_fName_1,save_dir);
    [MAIN_STUDY,MAIN_ALLEEG] = std_checkset(MAIN_STUDY,MAIN_ALLEEG);
    [MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                                            study_fName_1,save_dir,...
                                            'RESAVE_DATASETS','on');
    fprintf(1,'\n==== DONE: CLUSTERING SUBJECT DATA ====\n');
else
    fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',save_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',save_dir);
    end
    fprintf(1,'\n==== DONE: LOADING CLUSTER STUDY DATA ====\n');
end
%% INITIALIZE PARFOR LOOP VARS
fprintf('\nInitialize PARFOR Loop Vars\n');
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
%- clear vars for memory
% clear MAIN_ALLEEG
%% GENERATE EPOCH MAIN FUNC
%## PARFOR LOOP
parfor (subj_i = LOOP_VAR,POOL_SIZE)
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    fprintf('Running subject %s\n',EEG.subject)
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end
    EEG.icaact    = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    %## PARSE TRIALS
    fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED filesep [EVENT_COND_COMBOS{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[EVENT_COND_COMBOS{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    if exist([fPath filesep fName],'file') && false
        try
            ALLEEG = pop_loadset('filename',fName,'filepath',fPath);
            tmp{subj_i} = ALLEEG;
        catch e
            rmv_subj(subj_i) = 1;
            EEG.timewarp = struct([]);
            EEG.urevent = [];
            tmp{subj_i} = []; %EEG;
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
        end
    else
        try
            %## EPOCH
            ALLEEG = as_parse_trials(EEG,COND_CHARS,EVENT_CHARS,EPOCH_TIME_LIMITS,...
                    'EVENT_FIELD_TRIAL',COND_FIELD_PARSER,...
                    'EVENT_FIELD_CONDITION',EVENT_FIELD_PARSER);
            %## SAVE ONE BIG EEG FILE
            ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
            %{
            %- timewarp for all conditions
            events = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
            timewarp = make_timewarp(ALLEEG,events,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = median(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select( ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %}
            %- STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
            %- 
            ALLEEG = eeg_checkset(ALLEEG,'eventconsistency');
            ALLEEG = eeg_checkset(ALLEEG);
            ALLEEG = eeg_checkamica(ALLEEG);
            %- save
            [ALLEEG] = pop_saveset(ALLEEG,'savemode','twofiles',...
                    'filename',fName,...
                    'filepath',fPath,...
                    'version','6');
            tmp{subj_i} = ALLEEG;
        catch e
            rmv_subj(subj_i) = 1;
            EEG.timewarp = struct([]);
            EEG.urevent = [];
            tmp{subj_i} = []; %EEG;
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
        end
    end
end
%% SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
tmp = tmp(~cellfun(@isempty,tmp));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding %s %s\n',EEG.subject,addFs{j})
        end
    end 
    tmp{subj_i} = EEG;
end
tmp = cellfun(@(x) [[]; x], tmp);
%##
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName_3,...
                                'filename',study_fName_3,...
                                'filepath',save_dir);
%## (AS) ATTACH CLUSTER STRUCT
STUDY.cluster = MAIN_STUDY.cluster;
STUDY.urcluster = MAIN_STUDY.urcluster;
%## ADD ANATOMICAL LABELS
[~, atlas_cell] = add_anatomical_labels(STUDY,ALLEEG);
STUDY.etc.add_anatomical_labels = atlas_cell;
%## SAVE
% [ALLEEG,MAIN_STUDY] = parfunc_rmv_subjs(tmp,MAIN_STUDY,rmv_subj);
%- Save
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            study_fName_3,save_dir,...
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