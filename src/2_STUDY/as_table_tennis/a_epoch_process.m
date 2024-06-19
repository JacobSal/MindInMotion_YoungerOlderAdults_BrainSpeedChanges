%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: This code was modified from

%- run sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/AS/run_a_epoch_process.sh

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
%% DEFINE SOURCE DIRECTORY & CD ======================================== %%
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'AS'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
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
% SUBJ_ITERS={[1,2,3]};
%% (IMPORT) ============================================================ %%
% import_dir = 'R:\Ferris-Lab\astudnicki\Project_PingPong\STUDY\220620_iCanNoiseEMG\AutoSelectComponents\ManuallySelectComponents\Epoched_checked\ALL_EpochSubjectHit_NoWarp_ChangeEpochRejection_Downsample';
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset';
% group_i = 1;
% for subj_i = 1:length(SUBJ_PICS{group_i})
%     tmp_dir = [import_dir filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}];
%     tmpset = dir([tmp_dir filesep '*.set']);
%     tmpfdt = dir([tmp_dir filesep '*.fdt']);
%     fprintf('Copying Data...\n')
%     mkdir([save_dir filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}])
%     copyfile([tmpset.folder filesep tmpset.name],[save_dir filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}] );
%     copyfile([tmpfdt.folder filesep tmpfdt.name],[save_dir filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}] );
% end
%##
% import_dir = 'R:\Ferris-Lab\astudnicki\Project_PingPong\STUDY\220620_iCanNoiseEMG\AutoSelectComponents';
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset';
% tmpset = dir([import_dir filesep '*.set']);
% tmpfdt = dir([import_dir filesep '*.fdt']);
% for i = 1:length(SUBJ_PICS{1})
%     inds = cellfun(@(x) contains(x,SUBJ_PICS{1}(i)),{tmpset.name});
%     mkdir([save_dir filesep 'Pilot' SUBJ_PICS{1}{i} filesep 'full_set'])
%     copyfile([tmpset(inds).folder filesep tmpset(inds).name],[save_dir filesep 'Pilot' SUBJ_PICS{1}{i} filesep 'full_set']);
%     copyfile([tmpfdt(inds).folder filesep tmpfdt(inds).name],[save_dir filesep 'Pilot' SUBJ_PICS{1}{i} filesep 'full_set']);
% end
%% (PARAMETERS) ======================================================== %%
fprintf('\nData Processing Parameters\n');
%## Hard Defines
%- dataset specific
DATA_SET = 'AS_dataset';
%- study group and saving
% SESSION_NUMBER = '1';
SUFFIX_PATH_EPOCHED = 'EPOCHED';
SAVE_ALLEEG = true;
% SAVE_EEG = false; % saves EEG structures throughout processing
%- epoching parameters
% CONDLABEL_CHARS = {'competitive','cooperative','moving_serve','stationary_serve'};
% (08/18/2023) JS, not using this
% EVENTS_TIMEWARP = {'Subject_hit','Subject_receive','Subject_hit'};
% (08/18/2023) JS, not using this anymore
% COND_CHARS = {'2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
% (08/18/2023) JS, no longer epoching out '1Bounce_BM"
% EVENT_CHARS = {'Subject_hit','Subject_receive'};
% (08/18/2023) JS, no comment but marking that this is still in use.
% EPOCH_METHOD = 'event_locked';
% (08/18/2023) JS, switch from 'timewarp' to 'event_locked'
% EPOCH_TIME_LIMITS = [-0.150,0.5]; 
% (08/18/2023) JS, wasn't using this for timewarp epoching but now am.
% CONDLABEL_CHARS = {};
% (07/27/2023) this doesn't do anything
% (08/18/2023) it now does something;
%- epoching params
epoch_method = 'sub_epoch';
switch epoch_method
    case 'sliding_window'
        EPOCH_PARAMS = struct('epoch_method','sliding_window',...
            'event_chars',{{'Subject_hit','Subject_receive'}},...
            'events_timewarp',{{'Subject_hit','Subject_receive','Subject_hit'}},...
            'cond_chars',{{'2Bounce_Human','1Bounce_Human','Serve_Human'}},...
            'std_timewarp',3,...
            'struct_field_cond','bounces',...
            'struct_field_event','type',...
            'regexp_slidingwindow',{{'Standing_Baseline'}},...
            'window_length',1.5+1.5,...
            'percent_overlap',0,...
            'epoch_length_timelim',[-1.5,1.5]);
    case 'timewarp'
        EPOCH_PARAMS = struct('epoch_method','timewarp',...
            'event_chars',{{'Subject_hit','Subject_receive'}},...
            'events_timewarp',{{'Subject_hit','Subject_receive','Subject_hit'}},...
            'cond_chars',{{'2Bounce_Human','1Bounce_Human','Serve_Human'}},...
            'std_timewarp',3,...
            'struct_field_cond','bounces',...
            'struct_field_event','type',...
            'regexp_slidingwindow',{{}},...
            'epoch_length_timelim',[-1.5,1.5]);
    case 'event_locked'
        EPOCH_PARAMS = struct('epoch_method','event_locked',...
            'event_chars',{{'Subject_hit'}},...
            'events_timewarp',{{}},...
            'cond_chars',{{'2Bounce_Human','1Bounce_Human','Serve_Human'}},...
            'std_timewarp',0,...
            'struct_field_cond','bounces',...
            'struct_field_event','type',...
            'regexp_slidingwindow',{{}},...
            'epoch_length_timelim',[-1.5,1.5]);
        %- (12/28/2023) JS, may be interesting to understand what is going
        %on around 'Subject_receive' events... sticking with 'Subject_hit'
        %for now but changing epoching mechanism.
        %- (12/18/2023) JS, trying [-0.15,1.5] from [-0.25,0.5];
        %- (12/28/2023) JS, trying to understand brain act. before ball hit
        % so using [-1,0.5]
    case 'sub_epoch'
        EPOCH_PARAMS = struct('epoch_method','sub_epoch',...
            'event_chars',{{'Subject_hit'}},...
            'events_timewarp',{{}},...
            'cond_chars',{{'2Bounce_Human','1Bounce_Human','Serve_Human'}},...
            'std_timewarp',[],...
            'struct_field_cond','bounces',...
            'struct_field_event',[],...
            'regexp_slidingwindow',{{}},...
            'epoch_length_timelim',[]);
        EPOCH_PARAMS_stand = struct('epoch_method','sliding_window',...
            'event_chars',{{'Subject_hit'}},...
            'events_timewarp',{{}},...
            'cond_chars',{{}},...
            'std_timewarp',0,...
            'struct_field_cond','condlabel',...
            'struct_field_event','type',...
            'regexp_slidingwindow',{{'StandingBaseline'}},...
            'window_length',1.5+1.4960,...
            'percent_overlap',0,...
            'epoch_length_timelim',[-1.5,1.4960]);
    otherwise
        error('error. choose a valid EPOCH_METHOD');
end
%- datetime override
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
% dt = '08182023_bounces_1h_2h_2bm_JS';
% dt = '12182023_bounces_1h_2h_2bm_JS_0p25-1';
% dt = '12282023_bounces_1h_2bm_JS_n1-0p5';
% dt = '01182023_subjrec_2bounces_1h_2bm_JS_n5-1p5';
% dt = '01252023_subjrec_2bounces_rally_serve_human_JS_n5-1p5';
% dt = '01292023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75';
% dt = '01312023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75';
% dt = '02012023_subjrec_2bounces_rally_serve_human_epochfix_JS_n1p5-1p5';
dt = '02012023_subjrec_2bounces_rally_serve_human_epochfixfix_JS_n1p5-1p5';
%## Soft Define
%- combinations of events and conditions
EVENT_COND_COMBOS = cell(length(EPOCH_PARAMS.cond_chars)*length(EPOCH_PARAMS.event_chars),1);
cnt = 1;
for cond_i = 1:length(EPOCH_PARAMS.cond_chars)
    for event_i = 1:length(EPOCH_PARAMS.event_chars)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',EPOCH_PARAMS.cond_chars{cond_i},EPOCH_PARAMS.event_chars{event_i});
        cnt = cnt + 1;
    end
end
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = 'all_comps_study';
study_fName_2 = 'epoch_study';
cluster_info_fpath = [STUDIES_DIR filesep 'as_cluster_info' filesep 'cluster_info.mat'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
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
subject_chars    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
fNames_stand = cell(1,length([SUBJ_ITERS{:}]));
fPaths_stand = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep 'Pilot' SUBJ_PICS{group_i}{subj_i}];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        fNames{cnt} = tmp.name;
%         fNames{cnt} = sprintf('Pilot%s_iCanNoiseEMG_postAMICA_postDIPFIT_cat12_AutoSelectComponents_ManuallySelectComponents_GoProCheckedEvents_Epoched_ALL_NoWarp_ChangeEpochRejection_Downsample_v2.set',SUBJ_PICS{group_i}{subj_i});
        fPaths_stand{cnt} = [OUTSIDE_DATA_DIR filesep 'Pilot' SUBJ_PICS{group_i}{subj_i} filesep 'full_set'];
        tmp = dir([fPaths_stand{cnt} filesep '*.set']);
        fNames_stand{cnt} = tmp.name;
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('.set/.fdt Exists: %s\n',[fPaths{cnt} filesep fNames{cnt}]);
%         fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;

    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subject_chars{cnt} = ['Pilot' SUBJ_PICS{group_i}{subj_i}];
        tmp = strjoin(EPOCH_PARAMS.cond_chars,'_'); 
        conditions{cnt} = tmp; %tmp{:};
        groups{cnt} = num2str(group_i);
        sessions{cnt} = '1';
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%% CREATE STUDY
fprintf('\nCreating .study files\n');
%- params
%- NOTE: need a dipole fit before this step 
%## Create STUDY & ALLEEG structs
% if ~exist([save_dir filesep study_fName_1 '.study'],'file') %|| true
%     fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
%     [MAIN_ALLEEG] = as_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
%                         conditions,groups,sessions);
%     [MAIN_STUDY,MAIN_ALLEEG] = as_create_study(MAIN_ALLEEG,cluster_info_fpath,study_fName_1,save_dir);
%     [MAIN_STUDY,MAIN_ALLEEG] = std_checkset(MAIN_STUDY,MAIN_ALLEEG);
%     [MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
%                                             study_fName_1,save_dir,...
%                                             'RESAVE_DATASETS','on');
%     fprintf(1,'\n==== DONE: CLUSTERING SUBJECT DATA ====\n');
% else
%     fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
%     if ~ispc
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',save_dir);
%     else
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',save_dir);
%     end
%     fprintf(1,'\n==== DONE: LOADING CLUSTER STUDY DATA ====\n');
% end
[MAIN_ALLEEG] = as_create_alleeg(fNames,fPaths,subject_chars,save_dir,...
                        conditions,groups,sessions);
[MAIN_STUDY,MAIN_ALLEEG] = as_create_study(MAIN_ALLEEG,cluster_info_fpath,study_fName_1,save_dir);
[MAIN_STUDY,MAIN_ALLEEG] = std_checkset(MAIN_STUDY,MAIN_ALLEEG);
%% INITIALIZE PARFOR LOOP VARS
fprintf('\nInitialize PARFOR Loop Vars\n');
% fPaths = {MAIN_ALLEEG.filepath};
% fNames = {MAIN_ALLEEG.filename};
% LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(subject_chars));
rmv_subj = zeros(1,length(subject_chars));
alleeg_fpaths = cell(length(subject_chars),1);
%- clear vars for memory
% clear MAIN_ALLEEG
%% GENERATE EPOCH MAIN FUNC
fprintf('Running Epoching Process...\n');
%## PARFOR LOOP
% for subj_i = LOOP_VAR
parfor (subj_i = 1:length(subject_chars),SLURM_POOL_SIZE)
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
    epoched_fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED];
    fPath = [epoched_fPath filesep [EVENT_COND_COMBOS{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[EVENT_COND_COMBOS{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    EEG_stand = pop_loadset('filepath',fPaths_stand{subj_i},'filename',fNames_stand{subj_i});
    %- Recalculate ICA Matrices && Book Keeping
    EEG_stand = eeg_checkset(EEG_stand,'loaddata');
    if isempty(EEG_stand.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG_stand.subject);
        EEG_stand.icaact = (EEG_stand.icaweights*EEG_stand.icasphere)*EEG_stand.data(EEG_stand.icachansind,:);
    end
    EEG_stand.icaact    = reshape( EEG_stand.icaact, size(EEG_stand.icaact,1), EEG_stand.pnts, EEG_stand.trials);
    EEG_stand.subject = EEG.subject;
    %## MAKE SURE SETS HAVE SAME ICS
    [v1,i1] = setdiff(EEG.etc.ic_classification.ICLabel.classifications,EEG_stand.etc.ic_classification.ICLabel.classifications,'rows','stable')
    [v2,i2] = setdiff(EEG_stand.etc.ic_classification.ICLabel.classifications,EEG.etc.ic_classification.ICLabel.classifications,'rows','stable')
    if isempty(v1)&&~isempty(v2)
        fprintf('%s) removing component(s):',EEG_stand.subject); fprintf('%i,\n',i2);
        EEG_stand = pop_subcomp(EEG_stand,i2,0,0);
    elseif ~isempty(v1)&&isempty(v2)
        fprintf('%s) removing component(s):',EEG.subject); fprintf('%i,\n',i1);
        EEG = pop_subcomp(EEG,i1,0,0);
    elseif ~isempty(v1)&&~isempty(v2)
        fprintf('%s) removing component(s):',EEG.subject); fprintf('%i,\n',i1);
        EEG = pop_subcomp(EEG,i1,0,0);
        fprintf('%s) removing component(s):',EEG_stand.subject); fprintf('%i,\n',i2);
        EEG_stand = pop_subcomp(EEG_stand,i2,0,0);
    else
        fprintf('%s) .set have the same number of ic''s\n',EEG.subject);
    end
    
    if EEG.nbchan ~= EEG_stand.nbchan
        [vals,inds] = setdiff(1:length(EEG.chanlocs),1:length(EEG_stand.chanlocs))
        EEG = pop_select(EEG,'nochannel',inds);
    end
    if EEG.srate ~= EEG_stand.srate
        error('The two datasets must have the same sampling rate');
    end
    if EEG.trials > 1 || EEG_stand.trials > 1
        if EEG.pnts ~= EEG_stand.pnts
%                 error('The two epoched datasets must have the same number of points');
        end
        if EEG.xmin ~= EEG_stand.xmin
            EEG_stand.xmin = EEG.xmin;
%                 fprintf('Warning: the two epoched datasets do not have the same time onset, adjusted');
        end
        if EEG.xmax ~= EEG_stand.xmax
            EEG_stand.xmax = EEG.xmax;
%                 fprintf('Warning: the two epoched datasets do not have the same time offset, adjusted');
        end
    end
    [cell_of_struct] = chk_struct_fields({EEG,EEG_stand});
    EEG = cell_of_struct{1};
    EEG_stand = cell_of_struct{2};
%     cell_of_struct = cellfun(@(x) [[]; x], cell_of_struct);
%     tmp = pop_mergeset(EEG,EEG_stand)
%     tmp = pop_mergeset([ALLEEG;EEG_stande],1:length(ALLEEG)+1,1);
    try
        %## (FUNCTION) EPOCH
        %- extract hit data
        [ALLEEG,timewarp_struct] = as_parse_trials(EEG,'epoch_params',EPOCH_PARAMS);
        %- extract standing
        [EEG_stande,timewarp_struct_stand] = as_parse_trials(EEG_stand,'epoch_params',EPOCH_PARAMS_stand);      

        %## REMOVE USELESS FIELD EEG_stand
        if isfield(EEG_stande.event,'trialName')
            EEG_stande.event = rmfield(EEG_stande.event,'trialName');
        end
        if isfield(EEG_stande.event,'channel')
            EEG_stande.event = rmfield(EEG_stande.event,'channel');
        end
        if isfield(EEG_stande.event,'code')
            EEG_stande.event = rmfield(EEG_stande.event,'code');
        end
        if isfield(EEG_stande.event,'bvtime')
            EEG_stande.event = rmfield(EEG_stande.event,'bvtime');
        end
        if isfield(EEG_stande.event,'bvmknum')
            EEG_stande.event = rmfield(EEG_stande.event,'bvmknum');
        end
        if isfield(EEG_stande.event,'datetime')
            EEG_stande.event = rmfield(EEG_stande.event,'datetime');
        end
        if isfield(EEG_stande.event,'visible')
            EEG_stande.event = rmfield(EEG_stande.event,'visible');
        end
        if isfield(EEG_stande.event,'visible')
            EEG_stande.event = rmfield(EEG_stande.event,'visible');
        end
        %## REMOVE USELESS EVENT FIELDS
        for i = 1:length(ALLEEG)
            if isfield(ALLEEG(i).event,'trialName')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'trialName');
            end
            if isfield(ALLEEG(i).event,'channel')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'channel');
            end
            if isfield(ALLEEG(i).event,'code')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'code');
            end
            if isfield(ALLEEG(i).event,'bvtime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvtime');
            end
            if isfield(ALLEEG(i).event,'bvmknum')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvmknum');
            end
            if isfield(ALLEEG(i).event,'datetime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'datetime');
            end
            if isfield(ALLEEG(i).event,'visible')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'visible');
            end
            if isfield(ALLEEG(i).event,'visible')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'visible');
            end
        end
        %## SAVE ONE BIG EEG FILE
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_ALLEEG
            for i = 1:length(ALLEEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,ALLEEG(i).condition)];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(ALLEEG(i),...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition);
            end
            alleeg_fpaths{subj_i} = cond_files;
            %- append EEG_stande
%             REGEX_FNAME = 'cond_%s';
%             tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,EEG_stande.condition)];
%             if ~exist(tmp_fPath,'dir')
%                 mkdir(tmp_fPath)
%             end
%             [~] = pop_saveset(EEG_stande,...
%                 'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],EEG_stande.condition));
%             cond_files(i+1).fPath = tmp_fPath;
%             cond_files(i+1).fName = sprintf([REGEX_FNAME '.set'],EEG_stande.condition);
        end
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## AGGREGATE TIMEWARP
        if strcmp(EPOCH_PARAMS.epoch_method,'timewarp')
            %- timewarp for across condition
            timewarp = make_timewarp(ALLEEG,EPOCH_PARAMS.events_timewarp,'baselineLatency',0, ...
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
    %         ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
            %- 
        end
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
try
    fprintf('Bugged Subjects:\n');
    fprintf('%s\n',MAIN_ALLEEG(cellfun(@isempty,tmp)).subject);
    bugged_subjs = MAIN_ALLEEG(cellfun(@isempty,tmp)).subject;
    tmp = tmp(~cellfun(@isempty,tmp));
catch
    bugged_subjs = [];
    tmp = tmp(~cellfun(@isempty,tmp));
    fprintf('No Bugged Subjects.\n');
end
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
tmp = eeg_checkset(tmp,'eventconsistency');
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName_2,...
                                'filename',study_fName_2,...
                                'filepath',save_dir);
%## (AS) ATTACH CLUSTER STRUCT
STUDY.cluster = MAIN_STUDY.cluster;
STUDY.urcluster = MAIN_STUDY.urcluster;
%## ASSIGN PARAMETERS
STUDY.etc.a_epoch_process.epoch_params = EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_chars = EVENT_COND_COMBOS;
STUDY.etc.a_epoch_process.subjs_bugged = bugged_subjs;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            study_fName_2,save_dir,...
                                            'STUDY_COND',[]);
%%
SUBJ_PICS = {{MAIN_STUDY.datasetinfo.subject}};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
cluster_struct = STUDY.urcluster;
cluster_struct_orig = STUDY.urcluster;
% subj_chars_orig = SUBJ_PICS{1};
% subj_chars_orig = cellfun(@(x) [{} sprintf('Pilot%s',x)],subj_chars_orig);
subj_chars_orig = SUBJ_PICS{1};
subj_chars = {STUDY.datasetinfo.subject};
subj_keep = zeros(length(subj_chars),1);
for subj_i = 1:length(subj_chars_orig)
    disp(any(strcmp(subj_chars_orig{subj_i},subj_chars)));
    if any(strcmp(subj_chars_orig{subj_i},subj_chars))
        subj_keep(subj_i) = 1;
    end
end
subjs_rmv = find(~subj_keep);
subj_keep = find(subj_keep);
[val,ord] = sort(subj_keep);
for cli = 2:length(cluster_struct)
    si = cluster_struct(cli).sets;
    ci = cluster_struct(cli).comps;
    keep_si = setdiff(si,subjs_rmv);
    tmp = cluster_struct(cli).preclust.preclustdata;
    tmp_preclust = [];
    tmp_si = [];
    tmp_ci = [];
    for i = 1:length(keep_si)
        tmp_si = [tmp_si, repmat(ord(keep_si(i) == val),1,sum(keep_si(i) == si))];
        tmp_ci = [tmp_ci, ci(keep_si(i) == si)];
        
        tmp_preclust = [tmp_preclust; tmp(keep_si(i) == si,:)];
    end
    cluster_struct(cli).sets = tmp_si;
    cluster_struct(cli).comps = tmp_ci;
    cluster_struct(cli).preclust.preclustdata = tmp_preclust;
end
cluster_struct(1).sets = [cluster_struct(2:end).sets];
cluster_struct(1).comps = [cluster_struct(2:end).comps];
%-
STUDY.cluster = cluster_struct;
STUDY.etc.a_epoch_process.epoch_params = EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_chars = EVENT_COND_COMBOS;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                        STUDY.filename,STUDY.filepath,...
                        'RESAVE_DATASETS','off');
[MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                        'recovered',MAIN_STUDY.filepath,...
                        'RESAVE_DATASETS','off');
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