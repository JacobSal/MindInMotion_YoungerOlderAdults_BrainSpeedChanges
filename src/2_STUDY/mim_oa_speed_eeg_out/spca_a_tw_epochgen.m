%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/run_spca_a_tw_epochgen.sh

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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name & 
DATA_SET = 'MIM_dataset';
%- Study Name
STUDY_DIR_FNAME = '03232024_spca_analysis_OA';
%- Subject Directory information
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_FNAME_REST = 'rest_epoch_study_yaoa';
STUDY_FNAME_GAIT = 'gait_epoch_study_yaoa';
%-
studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
ica_data_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies' filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
save_dir = [studies_dir filesep sprintf('%s',STUDY_DIR_FNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%- study group and saving
% SAVE_ALLEEG = false;
%## EPOCH PARAMS
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
    'slide_cond_chars',{{'rest'}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'do_recalc_epoch',true);
%- compute measures for spectrum and ersp
ICLABEL_EYE_CUTOFF = 0.75;
% FORCE_RECALC_ERSP = false;
% DO_BASELINE_CORRECTION = false;
% DO_RECREATE_SPCA = false;
%- statistics & conditions
% speed_trials = {'0p25','0p5','0p75','1p0'};
% terrain_trials = {'flat','low','med','high'};
% COND_DESIGNS = {speed_trials,terrain_trials};
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% SPCA_PARAMS = struct('analysis_type','component',...
%     'event_char','RHS',...
%     'epoch_min_max',[1,4.25],...
%     'n_resamples',100,...
%     'timewarp_events',{{'RHS','LHS','LTO','RTO'}},...
%     'condition_base','rest',...
%     'condition_gait',{{'flat','low','med','high','0p25','0p5','0p75','1p0'}});
%% GRAB SUBJECT .SET & DIPFIT FILES ==================================== %%
subj_chars    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = zeros(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; 
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [ica_data_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
        tmp_gait = dir([fPaths{cnt} filesep '*.set']);
        try
            fNames{cnt} = tmp_gait.name;
            %- Chanlocs fPaths
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            dipfit_norm_fPaths(cnt) = 1;
        catch e
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('%s\n',getReport(e))
        end
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subj_chars{cnt} = SUBJ_PICS{group_i}{subj_i};
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
inds = logical(dipfit_norm_fPaths);
fPaths = fPaths(inds);
fNames = fNames(inds);
subj_chars = subj_chars(inds);
%% ===================================================================== %%
%## GENERATE EPOCH MAIN FUNC
tmp_gait = cell(1,length(fPaths));
tmp_rest = cell(1,length(fPaths));
% rmv_subj = zeros(1,length(fPaths));
%## PARFOR LOOP
% parfor (subj_i = 1:length(fPaths),floor(length(fPaths)/3))
parfor subj_i = 1:length(fPaths)
% for subj_i = 1:length(subjectNames)
    %## PARAM COPIES
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    %## CHECK DATA EXISTENCE
    epoched_fPath = [save_dir filesep subj_chars{subj_i} filesep 'GAIT_EPOCHED'];
    %- gait
    gait_fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
    gait_fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subj_chars{subj_i},[tmp_epoch_params.gait_trial_chars{:}]);
    if ~exist(gait_fpath,'dir')
        mkdir(gait_fpath)
    end
    %- rest
    rest_fpath = [epoched_fPath filesep tmp_epoch_params.slide_cond_chars{:}];
    rest_fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subj_chars{subj_i},'rest');
    if ~exist(rest_fpath,'dir')
        mkdir(rest_fpath)
    end
    if exist([gait_fpath filesep gait_fname],'file') && exist([rest_fpath filesep rest_fname],'file')
        fprintf('Epoch data exists... Loading subject %s\n',subj_chars{subj_i})
        EEG = pop_loadset('filepath',gait_fpath,'filename',gait_fname);
        tmp_gait{subj_i} = EEG;
        EEG = pop_loadset('filepath',rest_fpath,'filename',rest_fname);
        tmp_rest{subj_i} = EEG;
        EEG = struct.empty;
    %## CREATE EPOCH EEG
    else
        %## LOAD EEG DATA
        [~,EEG,~] = eeglab_loadICA(fNames{subj_i},fPaths{subj_i});
        fprintf('Running subject %s\n',EEG.subject)
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        EEG = iclabel(EEG);
        clssts = EEG.etc.ic_classification.ICLabel.classifications;
        bad_eye_ics = find(clssts(:,3) > ICLABEL_EYE_CUTOFF);
        EEG = pop_subcomp(EEG,bad_eye_ics,0,0);
        EEG = eeg_checkset(EEG,'loaddata');
        EEG.etc.spca.eye_ic_rej = bad_eye_ics;
        ics_orig = 1:size(EEG.icaweights,2);
        tmp_cut = ics_orig;
        tmp_cut(bad_eye_ics) = [];
        [valc,ordc] = sort(tmp_cut);
        unmix = [valc; ordc];
        EEG.etc.spca.unmix_mat = unmix;
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %## PARSE GAIT TRIALS
        %- REMOVE USELESS EVENT FIELDS (Improve Load Time)
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
        [tmp_eeg,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',tmp_epoch_params);
        tmp_eeg = pop_mergeset(tmp_eeg,1:length(tmp_eeg),1);
        %## timewarp for across condition
        if strcmp(tmp_epoch_params.epoch_method,'timewarp')
            timewarp = make_timewarp(tmp_eeg,tmp_epoch_params.tw_events,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(tmp_eeg.epoch),goodepochs);
            %- reject outlier strides
            tmp_eeg = pop_select(tmp_eeg,'notrial',sedi);
            %- store timewarp structure in EEG
            tmp_eeg.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            tmp_eeg.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            tmp_eeg.urevent = []; % might be needed
            tmp_eeg.etc.epoch.epoch_limits = tmp_epoch_params.epoch_time_lims;
        end
        %## STRUCT EDITS
        tmp_eeg.urevent = []; % might be needed
        tmp_eeg.etc.epoch.epoch_limits = tmp_epoch_params.epoch_time_lims;
        %- checks
        tmp_eeg = eeg_checkset(tmp_eeg,'eventconsistency');
        tmp_eeg = eeg_checkset(tmp_eeg);
        tmp_eeg = eeg_checkamica(tmp_eeg);
        %- save
        [tmp_eeg] = pop_saveset(tmp_eeg,'savemode','twofiles',...
                'filename',gait_fname,...
                'filepath',gait_fpath,...
                'version','6');
        tmp_gait{subj_i} = tmp_eeg;
        %## RESTING STATE
        % tmp_eeg = {};
        tmp_epoch_params.slide_cond_chars = {'rest'};
        tmp_epoch_params.epoch_method = 'sliding_window';
        tmp_epoch_params.percent_overlap = 0;
        [tmp_eeg,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',tmp_epoch_params);
        %- save
        [tmp_eeg] = pop_saveset(tmp_eeg,'savemode','twofiles',...
                'filename',rest_fname,...
                'filepath',rest_fpath,...
                'version','6');
        tmp_rest{subj_i} = tmp_eeg;
        EEG = struct.empty;
    end
end
%% ===================================================================== %%
%## SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s\n',subj_chars{cellfun(@isempty,tmp_rest)});
tmp_rest = tmp_rest(~cellfun(@isempty,tmp_rest));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
tmp_rest = util_resolve_struct(tmp_rest);
%##
[STUDY, ALLEEG] = std_editset([],tmp_rest,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_REST,...
                                'filename',STUDY_FNAME_REST,...
                                'filepath',save_dir);
[STUDY_REST,ALLEEG_REST] = std_checkset(STUDY,ALLEEG);
% [STUDY_REST,ALLEEG_REST] = parfunc_save_study(STUDY_REST,ALLEEG_REST,...
%                                         STUDY_REST.filename,STUDY_REST.filepath,...
%                                         'RESAVE_DATASETS','off');
STUDY = struct.empty;
ALLEEG = struct.empty;
%## ersp plot per cluster per condition
STUDY_REST = pop_statparams(STUDY_REST,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY_REST = pop_erspparams(STUDY_REST,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',[]);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_REST = pop_specparams(STUDY_REST,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
[~,~] = parfunc_save_study(STUDY_REST,ALLEEG_REST,...
                                        STUDY_REST.filename,STUDY_REST.filepath,...
                                        'RESAVE_DATASETS','off');
STUDY_REST = struct.empty;
ALLEEG_REST = struct.empty;
%% ===================================================================== %%
%## SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s\n',subj_chars{cellfun(@isempty,tmp_gait)});
tmp_gait = tmp_gait(~cellfun(@isempty,tmp_gait));
tmp_gait = util_resolve_struct(tmp_gait);
%##
[STUDY, ALLEEG] = std_editset([],tmp_gait,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_GAIT,...
                                'filename',STUDY_FNAME_GAIT,...
                                'filepath',save_dir);
[STUDY_GAIT,ALLEEG_GAIT] = std_checkset(STUDY,ALLEEG);
STUDY = struct.empty;
ALLEEG = struct.empty;
% [STUDY_GAIT,ALLEEG_GAIT] = parfunc_save_study(STUDY,ALLEEG,...
%                                         STUDY.filename,STUDY.filepath,...
%                                         'RESAVE_DATASETS','off');
%## CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG_GAIT)
    %- assign percondition timewarping
    ALLEEG_GAIT(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG_GAIT(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG_GAIT),size(ALLEEG_GAIT(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG_GAIT)
    allWarpTo(subj_i,:) = ALLEEG_GAIT(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])

%## (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG_GAIT(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
STUDY_GAIT.etc.averaged_warpto_events = averaged_warpto_events;
STUDY_GAIT.etc.timewarp_ntimes = TIMEWARP_NTIMES;
STUDY_GAIT.etc.ersp_crop_times = ERSP_CROP_TIMES;
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
disp(averaged_warpto_events);
%## ersp plot per cluster per condition
STUDY_GAIT = pop_statparams(STUDY_GAIT,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY_GAIT = pop_erspparams(STUDY_GAIT,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_CROP_TIMES);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_GAIT = pop_specparams(STUDY_GAIT,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
[~,~] = parfunc_save_study(STUDY_GAIT,ALLEEG_GAIT,...
                                        STUDY_GAIT.filename,STUDY_GAIT.filepath,...
                                        'RESAVE_DATASETS','off');
