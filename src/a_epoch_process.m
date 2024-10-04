%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_a_epoch_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- study group and saving
SESSION_NUMBER = '1';
SAVE_ALLEEG = true;
SAVE_EEG = true; %true;
OVERRIDE_DIPFIT = true;
%- epoching params
DO_SLIDING_WINDOW = false;
RECALC_ICA_STUDY = false;
% %* sliding window
% WINDOW_LENGTH = 6; % sliding window length in seconds
% PERCENT_OVERLAP = 0.0; % percent overlap between epochs
% %* gait
% EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
% STD_TIMEWARP = 3;
% EPOCH_TIME_LIMITS = [-0.5,4.5]; %[-1,3]; %[-0.5,5]; % [-1,3] captures gait events well , [-0.5,5] captures gait events poorly
% % (10/13/2023) changing from [-1,4.25] to [-0.5,4.5] to match chang's
% % (10/25/2023) changing from [-0.5,4.5] to [-1,4.25] as it seems to help
% % with frequency decomposition artifact during ERSP creation
% % paper
% % (01/23/2024) changing from [-1,4.25] to [-0.5,4.5] to match chang
% TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
% if DO_SLIDING_WINDOW
%     SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
%     DEF_EPOCH_PARAMS.gait_trial_chars = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
% else
%     SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
%     DEF_EPOCH_PARAMS.gait_trial_chars = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% end
%## EPOCH PARAMS
SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
DEF_EPOCH_PARAMS = struct('epoch_method','timewarp',...
    'percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'baseline_time_lims',[-0.5,4.5-2],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'cond_field','cond',...
    'appx_cond_len',3*60,...
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%- Study Name
STUDY_DIR_FNAME = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- Subject Directory information
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_FNAME_ALLCOMP = 'all_comps_study';
STUDY_FNAME_EPOCH = 'epoch_study';
%## soft define
studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
ica_data_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies' filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
save_dir = [studies_dir filesep sprintf('%s',STUDY_DIR_FNAME)];
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
        fPaths{cnt} = [ica_data_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        try
            fNames{cnt} = tmp.name;
            %- Chanlocs fPaths
    %         chanlocs_fPaths{cnt} = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
            chanlocs_fPaths{cnt} = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
    %         dipfit_fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'head_model' filesep 'dipfit_struct.mat'];
%             dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm.mat'];
            dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm_ants.mat'];
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
    %         fprintf('DIPFIT Exists: %i\n',exist(dipfit_fPaths{cnt},'file'));
            fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{cnt},'file'));
        catch e
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('%s\n',getReport(e))
            dipfit_norm_fPaths{cnt} = [];
        end
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(DEF_EPOCH_PARAMS.gait_trial_chars,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = GROUP_NAMES{group_i};
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
inds = logical(cellfun(@(x) exist(x,'file'),dipfit_norm_fPaths));
chanlocs_fPaths = chanlocs_fPaths(inds);
dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
fPaths = fPaths(inds);
fNames = fNames(inds);
sessions = sessions(inds);
groups = groups(inds);
conditions = conditions(inds);
subjectNames = subjectNames(inds);
%% Create STUDY & ALLEEG structs
if ~exist([save_dir filesep STUDY_FNAME_ALLCOMP '.study'],'file') || RECALC_ICA_STUDY
    try
        fprintf('Creating ALLEEG...\n');
        [MAIN_ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                            conditions,groups,sessions);
        fprintf('Generating STUDY...\n');
        [MAIN_STUDY,MAIN_ALLEEG] = mim_create_study(MAIN_ALLEEG,STUDY_FNAME_ALLCOMP,save_dir);
        [MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                                        MAIN_STUDY.filename,MAIN_STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
    catch e
        fprintf('\n%s\n',getReport(e));
        exit();
    end
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_ALLCOMP '_UNIX.study'],'filepath',save_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_ALLCOMP '.study'],'filepath',save_dir);
    end
end

%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
alleeg_fpaths = cell(length(MAIN_ALLEEG),1);
%- clear vars for memory
% clear MAIN_ALLEEG
%% UNIT TEST =========================================================== %%
EEG = pop_loadset('filepath',[fPaths{1}],'filename',[fNames{1}]);
%- Recalculate ICA Matrices && Book Keeping
EEG = eeg_checkset(EEG,'loaddata');
if isempty(EEG.icaact)
    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end
%## PARSE TRIALS
epoched_fPath = [save_dir filesep EEG.subject];
fPath = [epoched_fPath filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}]];
fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[DEF_EPOCH_PARAMS.gait_trial_chars{:}]);
if ~exist(fPath,'dir')
    mkdir(fPath)
end
%- parse
%## REMOVE USELESS EVENT FIELDS (Improve Load Time)
for i = 1:length(EEG)
    if isfield(EEG.event,'trialName')
        EEG.event = rmfield(EEG.event,'trialName');
    end
    if isfield(EEG.event,'channel')
        EEG.event = rmfield(EEG.event,'channel');
    end
    if isfield(EEG.event,'code')
        EEG.event = rmfield(EEG.event,'code');
    end
    if isfield(EEG.event,'bvtime')
        EEG.event = rmfield(EEG.event,'bvtime');
    end
    if isfield(EEG.event,'bvmknum')
        EEG.event = rmfield(EEG.event,'bvmknum');
    end
    if isfield(EEG.event,'datetime')
        EEG.event = rmfield(EEG.event,'datetime');
    end
end
%## EPOCH
[ALLEEG,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',DEF_EPOCH_PARAMS);

%% GENERATE EPOCH MAIN FUNC
%## PARFOR LOOP
parfor (subj_i = LOOP_VAR,SLURM_POOL_SIZE)
    %## LOAD EEG DATA
    EEG = MAIN_ALLEEG(subj_i);
%     EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    fprintf('Running subject %s\n',EEG.subject)
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    
    %## PARSE TRIALS
    epoched_fPath = [save_dir filesep EEG.subject];
    fPath = [epoched_fPath filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[DEF_EPOCH_PARAMS.gait_trial_chars{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    try
         %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
        for i = 1:length(EEG)
            if isfield(EEG.event,'trialName')
                EEG.event = rmfield(EEG.event,'trialName');
            end
            if isfield(EEG.event,'channel')
                EEG.event = rmfield(EEG.event,'channel');
            end
            if isfield(EEG.event,'code')
                EEG.event = rmfield(EEG.event,'code');
            end
            if isfield(EEG.event,'bvtime')
                EEG.event = rmfield(EEG.event,'bvtime');
            end
            if isfield(EEG.event,'bvmknum')
                EEG.event = rmfield(EEG.event,'bvmknum');
            end
            if isfield(EEG.event,'datetime')
                EEG.event = rmfield(EEG.event,'datetime');
            end
        end
        %## EPOCH
        [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',DEF_EPOCH_PARAMS);
        
       
        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
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
                [~] = pop_saveset(ALLEEG(i),'savemode','twofiles',...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition);
            end
            alleeg_fpaths{subj_i} = cond_files;
        end
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## timewarp for across condition
        if strcmp(DEF_EPOCH_PARAMS.epoch_method,'timewarp')
            timewarp = make_timewarp(ALLEEG,DEF_EPOCH_PARAMS.tw_events,'baselineLatency',0, ...
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
            ALLEEG.etc.epoch.epoch_limits = DEF_EPOCH_PARAMS.epoch_time_lims;
        end
        %## STRUCT EDITS
        ALLEEG.urevent = []; % might be needed
        ALLEEG.etc.epoch.epoch_limits = DEF_EPOCH_PARAMS.epoch_time_lims;
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
fprintf('Bugged Subjects: %s',MAIN_ALLEEG(cellfun(@isempty,tmp)).subject);
tmp = tmp(~cellfun(@isempty,tmp));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
tmp = util_resolve_struct(tmp,{tmp.subject});
%##
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_EPOCH,...
                                'filename',STUDY_FNAME_EPOCH,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
STUDY.etc.a_epoch_process.epoch_params = DEF_EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
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

