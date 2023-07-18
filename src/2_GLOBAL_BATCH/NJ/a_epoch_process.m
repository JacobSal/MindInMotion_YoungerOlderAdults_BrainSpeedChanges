%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/NJ/run_a_epoch_process.sh

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
%% (DATASET INFORMATION) =============================================== %%
%## DATASET SPECIFIC
SUBJ_STANDING_OLD = {'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34',...
                     'S35','S36','S37','S38','S39','S40','S41','S42','S43','S44',...
                     'S46','S47','S48'};
SUBJ_STANDING = {'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34',...
                 'S35','S36','S37','S38','S39','S41','S42','S43','S44',...
                 'S46','S47','S48'};
SUBJ_GAIT = {'S18','S19','S20','S21','S22','S23','S24','S25','S26','S27',...
             'S28','S29','S30','S31','S32','S33','S34',...
             'S35','S36','S37','S38','S39','S41','S42','S43','S44',...
             'S46','S47','S48'};
%- STANDING OLD GROUP      
% SUBJ_PICS = {SUBJ_STANDING_OLD};
% SUBJ_ITERS = {1:length(SUBJ_STANDING_OLD)};
% GROUP_NAMES = {'young adults'};
%- STANDING GROUP
SUBJ_PICS = {SUBJ_STANDING};
SUBJ_ITERS = {1:length(SUBJ_STANDING)};
GROUP_NAMES = {'young adults'};
%- GAIT GROUP
% SUBJ_PICS = {SUBJ_GAIT};
% SUBJ_ITERS = {1:length(SUBJ_GAIT)};
% GROUP_NAMES = {'young adults'};
%- TEST
% SUBJ_PICS = {SUBJ_STANDING};
% SUBJ_ITERS = {[1,2]};
% GROUP_NAMES = {'young adults'};

%% (PARAMETERS) ======================================================== %%
%## hard define
%- save epoch'd eegs as seperate files (for connectivity)
SAVE_ALLEEG = true;
%- epoching params
DO_SLIDING_WINDOW = true;
%* sliding window
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
%* gait
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-1,3]; %[-1,3]; %[-0.5,5]; % [-1,3] captures gait events well , [-0.5,5] captures gait events poorly
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
if DO_SLIDING_WINDOW
    SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
    TRIAL_TYPES = {'pre','post'};
else
    SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
    TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
end
%- Dataset Params
%* datetime override
% dt = '06292023_NJ_Standing';
dt = '07162023_NJ_standing_customheadmods';
% dt = 'test';
SUBFOLDER_CHAR = 'standing';
SESSION_NUMBER = '1';
%* hardcode data_dir
DATA_SET = 'jacobsenN_dataset';
%## soft define
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET]; 
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
study_fName_1 = sprintf('%s_all_comps_study',[TRIAL_TYPES{:}]);
study_fName_2 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
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
% dipfit_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = cell(1,length([SUBJ_ITERS{:}]));
% vol_fPaths = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep SUBFOLDER_CHAR];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        fNames{cnt} = tmp.name;
        %- Chanlocs fPaths
%         chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
%         chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
%         dipfit_fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'head_model' filesep 'dipfit_struct.mat'];
%         dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm.mat'];
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file'))
%         fprintf('DIPFIT Exists: %i\n',exist(dipfit_fPaths{cnt},'file'));
%         fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{cnt},'file'));
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = GROUP_NAMES{group_i};
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
% inds = logical(cellfun(@(x) exist(x,'file'),dipfit_norm_fPaths));
% chanlocs_fPaths = chanlocs_fPaths(inds);
% dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
% fPaths = fPaths(inds);
% fNames = fNames(inds);
% sessions = sessions(inds);
% groups = groups(inds);
% conditions = conditions(inds);
% subjectNames = subjectNames(inds);
%% CREATE STUDY
%## Create STUDY & ALLEEG structs
if ~exist([save_dir filesep study_fName_1 '.study'],'file') %|| true
    fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
    [MAIN_ALLEEG] = nj_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                        conditions,groups,sessions); %,...
%                         'SAVE_EEG',SAVE_EEG); %,...
%                         'CHANLOCS_FPATHS',chanlocs_fPaths);
    [MAIN_STUDY,MAIN_ALLEEG] = nj_create_study(MAIN_ALLEEG,study_fName_1,save_dir);
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
% parfor (subj_i = LOOP_VAR,POOL_SIZE)
for subj_i = LOOP_VAR
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    fprintf('Running subject %s\n',EEG.subject)
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    
    %## PARSE TRIALS
    epoched_fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED];
    fPath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[TRIAL_TYPES{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    try
        %## EPOCH
        [ALLEEG,timewarp_struct] = nj_parse_trials(EEG,DO_SLIDING_WINDOW,...
            'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,...
            'STD_TIMEWARP',STD_TIMEWARP,...
            'COND_CHARS',TRIAL_TYPES);
        %## REMOVE USELESS EVENT FIELDS
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_ALLEEG
            for i = 1:length(ALLEEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fName = sprintf(REGEX_FNAME,ALLEEG(i).condition);
                tmp_fPath = [epoched_fPath filesep tmp_fName];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(ALLEEG(i),...
                    'filepath',tmp_fPath,'filename',[tmp_fName '.set']);
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = [tmp_fName '.set'];
            end
        end
        %## MERGE SETS
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## TIMEWARPING?
        if ~DO_SLIDING_WINDOW
            ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
            ALLEEG.etc.cond_files = cond_files;
            %- timewarp for across condition
            timewarp = make_timewarp(ALLEEG,TIMEWARP_EVENTS,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select(ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
        end
        %- checks
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
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
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
                                'name',study_fName_2,...
                                'filename',study_fName_2,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
                                        
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

