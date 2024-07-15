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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- study group and saving
SESSION_NUMBER = '1';
SAVE_ALLEEG = true;
SAVE_EEG = true; 
OVERRIDE_DIPFIT = true;
%- datetime override
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- Subject Directory information
preproc_dir_name = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%## soft define
DATA_DIR = [PATHS.src_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep preproc_dir_name]; % JACOB,SAL(02/23/2023)
study_fName_1 = 'all_comps_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name)];
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
dipfit_norm_fPaths = cell(1,length([SUBJ_ITERS{:}]));
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
        try
            fNames{cnt} = tmp.name;
            %- Dipfit
            dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm_ants.mat'];
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')));
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
fprintf('Creating ALLEEG...\n');
%## CREATE ALLEEG STRUCTURE & ADD ICLABEL
[ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                    conditions,groups,sessions);
fprintf('Generating STUDY...\n');
%## AUTOMATICALLY REMOVE NON-BRAIN COMPONENTS
% [MAIN_STUDY,MAIN_ALLEEG] = mim_create_study(ALLEEG,study_fName_1,save_dir);
%## SAVE STUDY
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName_1,...
                                'filename',study_fName_1,...
                                'filepath',save_dir);
[MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                                MAIN_STUDY.filename,MAIN_STUDY.filepath,...
                                'RESAVE_DATASETS','on');

%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
alleeg_fpaths = cell(length(MAIN_ALLEEG),1);
%- clear vars for memory
% clear MAIN_ALLEEG
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
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i})';
    disp(size(fields(tmp{subj_i})'));
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false); 
    out = [out{:}];
    addFs = fsPrev(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
        end
    end 
%     tmp{subj_i} = EEG;
    tmp{subj_i} = orderfields(EEG);
end
%- CONCATENATE tmp
tmp = cellfun(@(x) [[]; x], tmp);
%##
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName_2,...
                                'filename',study_fName_2,...
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

