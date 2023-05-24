%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_YA_proc/run_epoch_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'gamma' filesep 'MIM_YA_proc'];
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
    % see. eeg_options.m 
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
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
%- MIND IN MOTION
SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)
%- Subject Picks
SUBJ_PICS = {SUBJ_YNG};
SUBJ_ITERS = {1:length(SUBJ_YNG)}; % CHANG,LIU(12/26/2023)
GROUP_NAMES = {'H1000s'};
%- Subject Directory Information
OA_PREP_FPATH = '04182023_YA_N37_prep_verified'; % JACOB,SAL(04/10/2023)
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
%% ===================================================================== %%
%## PROCESSING PARAMS
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
% params.TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
params.TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
params.EPOCH_TIME_LIMITS = [-0.5,4.5]; % [-1,3] captures gait events well
params.TRIAL_LENGTH = 3*60; % trial length in seconds
params.PER_OVERLAP = 0.0; % percent overlap between epochs
params.TRIAL_BEGIN_STR = 'TrialStart';
params.TRIAL_END_STR = 'TrialEnd';
params.EVENT_TRIAL_PARSER = 'type';
params.EVENT_COND_PARSER = 'cond';
%- eeglab_cluster.m spectral params
params.FREQ_LIMITS = [1,100];
params.CYCLE_LIMITS = [3,0.8];
params.SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
params.FREQ_FAC = 4;
params.PAD_RATIO = 2;
%% global script chain (VERSION 1)
%- datetime override
dt = '05082023_YA_N33_fem_customelec';
%## PATH & TEMP STUDY NAME
%- hard define
SAVE_EEG = true;
SESSION_NUMBER = '1';
SUFFIX_PATH_EPOCHED = 'EPOCHED';
study_fName_1 = sprintf('%s_all_comps_study',[params.TRIAL_TYPES{:}]);
study_fName_2 = sprintf('%s_reduced_comps_study',[params.TRIAL_TYPES{:}]);
study_fName_3 = sprintf('%s_EPOCH_study',[params.TRIAL_TYPES{:}]);
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
%% Store fNames and fPaths
conditions      = cell(1,length([SUBJ_ITERS{:}]));
groups          = cell(1,length([SUBJ_ITERS{:}]));
sessions        = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
chanlocs_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
vol_fPaths = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        fNames{cnt} = tmp.name;
        %- Chanlocs fPaths
%         chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
        chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep '' filesep 'CustomElectrodeLocations.mat'];
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(params.TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = GROUP_NAMES{group_i};
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%% CREATE STUDY
%## Create STUDY & ALLEEG structs
if ~exist([save_dir filesep study_fName_1 '.study'],'file') %|| true
    fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
    [MAIN_ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                        conditions,groups,sessions,...
                        'SAVE_EEG',SAVE_EEG,...
                        'CHANLOCS_FPATHS',chanlocs_fPaths);
    [MAIN_STUDY,MAIN_ALLEEG] = mim_create_study(MAIN_ALLEEG,study_fName_1,save_dir);
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
    fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED filesep [params.TRIAL_TYPES{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[params.TRIAL_TYPES{:}]);
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
            [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,params.TRIAL_TYPES,params.EPOCH_TIME_LIMITS);

            %## SAVE ONE BIG EEG FILE
            ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
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
            %{
            new_warp = [];
            new_warp.latencies = cat(1,timewarp_struct.latencies);
            new_warp.epochs = cat(2,timewarp_struct.epochs);
            new_warp.eventSequence = timewarp_struct(1).eventSequence;
            new_warp.warpto = timewarp_struct(1).warpto;
            %}
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %- STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = params.EPOCH_TIME_LIMITS;
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
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
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
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
%- NOTES: spec_mode, 'psd' does not work.
[STUDY,ALLEEG] = eeglab_cluster(STUDY,ALLEEG,...
                    'SPEC_MODE',params.SPEC_MODE,...
                    'FREQ_LIMITS',params.FREQ_LIMITS,...
                    'CYCLE_LIMITS',params.CYCLE_LIMITS,...
                    'FREQ_FAC',params.FREQ_FAC,...
                    'PAD_RATIO',params.PAD_RATIO);
%## PCA reduction algorithm
%- RECALCULATE ICAACT MATRICES
%     ALLEEG = eeg_checkset(ALLEEG,'loaddata');
%     for subj_i = 1:length(ALLEEG)
%         if isempty(ALLEEG(subj_i).icaact)
%             fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
%             ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
%         end
%     end
% [STUDY,ALLEEG,~,~] = cluster_pca_reduce(STUDY,ALLEEG);
%## ICA reduction algorithm
[STUDY,~,~] = cluster_ica_reduce(STUDY);
%## ADD ANATOMICAL LABELS
[~, atlas_cell] = add_anatomical_labels(STUDY,ALLEEG);
STUDY.etc.add_anatomical_labels = atlas_cell;
STUDY.etc.pipe_params = params;
% [ALLEEG,MAIN_STUDY] = parfunc_rmv_subjs(tmp,MAIN_STUDY,rmv_subj);
%- Save
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            study_fName_3,save_dir,...
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