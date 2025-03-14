%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/run_a_epoch_process.sh

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
ADD_ALL_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = STUDY_DIR;
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
    SRC_DIR = STUDY_DIR;
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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('ya');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- study group and saving
SESSION_CHAR = '1';
SAVE_ALLEEG = false;
%- epoching params
RECALC_ICA_STUDY = true;
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
    'do_recalc_epoch',true,...
    'do_combine_trials',true);
DO_RELOAD_ICA = true;
DEF_LOAD_STRUCT = struct('do_bem_dipfit',false,...
    'do_load_dipfit',true,...
    'dipfit_fpath',{''},...
    'chanloc_fpath',{''},...
    'mni_coord_transform',[0 0 0; 0 0 -1.5708; 1 1 1],...
    'dip_num',1,...
    'dip_plot','off');
REJ_STRUCT = struct( ...
    'powpow_params',struct('do_calc',true,...
        'upper_freq_lim',100, ...
        'input_data_int',2, ...
        'method_int',2, ...
        'num_iters',100, ...
        'do_plot',true), ...
    'psd_params',struct('fit_freqs',(2:40), ...
        'calc_freqs',[2,100], ...
        'slope_thresh',-0.2, ...
        'spec_perc',80), ...
    'iclabel_params',struct('iclabel_ver','lite', ...
        'class_ints_keep',[1,1,1], ...
        'class_threshs_keep',[0.5,0.75,0.9], ...
        'class_ints_rmv',[2,3], ...
        'class_threshs_rmv',[0.5,0.5]),...
    'do_valid_plots',true, ...
    'ica_rv_thresh',0.15, ...
    'plot_save_dir',{''});
CHK_STRUCT = struct( ...
    'brain_thresh_int',8, ...
    'brain_exempt_int',9, ...
    'cond_chk',struct('do_reject',true, ...
        'conds',{{'0p25','0p5','0p75','1p0','rest'}}), ... % ,'flat','high','low','med'
    'ic_cut_off',5);

%% (PATHS) ========================================================== %%
%- datset name
DATA_SET = 'MIM_dataset';
%- Study Name
% STUDY_DNAME = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
STUDY_DNAME = '01192025_mim_yaoa_nopowpow_crit_speed';
% STUDY_DNAME = 'dummy_study';
%- Subject Directory information
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_FNAME_CONT = 'contin_study';
STUDY_FNAME_EPOCH = 'epoch_study';
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
ica_data_dir = [studies_dir filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
group_cats = {'H1000','H2000','H3000'};
ICA_REGEXP = '%s_cleanEEG_EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set';
%-
subj_chars = [SUBJ_PICS{:}];
cond_chars  = cell(1,length(subj_chars));
group_chars = cell(1,length(subj_chars));
sess_chars = cell(1,length(subj_chars));
eeg_fnames = cell(1,length(subj_chars));
eeg_fpaths = cell(1,length(subj_chars));
chanlocs_fpaths = cell(1,length(subj_chars));
%(01/08/2025) JS, chanlocs_fpaths no longer necessary with new dipfit pipeline
dipfit_norm_fpaths = cell(1,length(subj_chars));
keep_inds = zeros(1,length(subj_chars));
%- loop
for subj_i = 1:length(subj_chars)
    %## PARAMS
    chk1 = false;
    chk2 = false;
    chk3 = false;
    fprintf('==== (%i) Loading Subject %s ====\n',subj_i,subj_chars{subj_i})
    %## GET GROUP
    vals = regexp(subj_chars{subj_i},'[nN]?[hH](\d)\d*','tokens');
    group_chars(subj_i) = group_cats(double(string(vals{1})));
    %## GET DIPFIT
    %- eeg fpath
    eeg_fpaths{subj_i} = [ica_data_dir filesep subj_chars{subj_i} filesep 'clean'];
    eeg_fnames{subj_i} = sprintf(ICA_REGEXP,subj_chars{subj_i});
    %- test existence of variables
    if exist([eeg_fpaths{subj_i} filesep eeg_fnames{subj_i}],'file')
        %- Chanlocs fPaths
        chanlocs_fpaths{subj_i} = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep subj_chars{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
        dipfit_norm_fpaths{subj_i} = [eeg_fpaths{subj_i} filesep 'dipfit_fem_norm_ants.mat'];
        %- Prints
        chk1 = exist([eeg_fpaths{subj_i} filesep eeg_fnames{subj_i}],'file') && exist([eeg_fpaths{subj_i} filesep 'W'],'file');
        chk2 = exist(dipfit_norm_fpaths{subj_i},'file');
        chk3 = exist(chanlocs_fpaths{subj_i},'file');
        fprintf('ICA exists: %i\n',chk1);
        fprintf('Normalized dipfit exists: %i\n',chk2);
        fprintf('Chanlocs exists: %i\n',chk3);
        %## ADD ADDITIONAL IDENTIFIERS
        ALLEEG = join(DEF_EPOCH_PARAMS.gait_trial_chars,'_'); 
        cond_chars{subj_i} = ALLEEG{:};
        sess_chars{subj_i} = SESSION_CHAR;
    end
    %-
    if ~(chk1 && chk2)   %~(chk1 && chk2 && chk3)       
        %(01/08/2025) JS, removing chk3 because of new dipfit pipeline 
        fprintf(2,'WARNING: Rejecting subject %s...\n',subj_chars{subj_i});
    end
    %-
    keep_inds(subj_i) = (chk1 && chk2); %(chk1 && chk2 && chk3);
    %(01/08/2025) JS, removing chk3 because of new dipfit pipeline 
end
%- remove subjects without a dipole fit
keep_inds = logical(keep_inds);
chanlocs_fpaths = chanlocs_fpaths(keep_inds);
dipfit_norm_fpaths = dipfit_norm_fpaths(keep_inds);
eeg_fpaths = eeg_fpaths(keep_inds);
eeg_fnames = eeg_fnames(keep_inds);
sess_chars = sess_chars(keep_inds);
group_chars = group_chars(keep_inds);
cond_chars = cond_chars(keep_inds);
subj_chars = subj_chars(keep_inds);

%% Create STUDY & ALLEEG structs
if ~exist([save_dir filesep STUDY_FNAME_CONT '.study'],'file') || RECALC_ICA_STUDY
    tmp_rmv_subjs = zeros(1,length(subj_chars));
    tmp_rej_crit_out = cell(1,length(subj_chars));    
    %## DIPOLE REJECTION
    fprintf('Creating ALLEEG...\n');
    ALLEEG = cell(1,length(subj_chars));
    parfor subj_i = 1:length(subj_chars)  
    % parfor subj_i = 1:5 %length(subj_chars) 
        try
            EEG = struct.empty;
            tmp_save_dir = [save_dir filesep subj_chars{subj_i} filesep 'ICA'];
            if ~exist(tmp_save_dir,'dir')
                mkdir(tmp_save_dir);
            end
            tmp_rej_struct = REJ_STRUCT;
            tmp_load_struct = DEF_LOAD_STRUCT;
            tmp_chk_struct = CHK_STRUCT;
            %## LOAD EEG & REJECT ICS
            fprintf('%s) Loading EEG data...\n',subj_chars{subj_i});
            tmp_load_struct.dipfit_fpath = dipfit_norm_fpaths{subj_i};
            % tmp_struct.chanloc_fpath = chanlocs_fpaths{subj_i}; 
            %(01/05/2025) JS, technically don't need anymore
            if DO_RELOAD_ICA || ~exist([tmp_save_dir filesep eeg_fnames{subj_i}],'file')
                %## LOAD ICA DATA
                EEG = mim_load_one_subj(eeg_fnames{subj_i},eeg_fpaths{subj_i},...
                    subj_chars{subj_i},cond_chars{subj_i},group_chars{subj_i},sess_chars{subj_i},...
                    'LOAD_STRUCT',tmp_load_struct)

                %## REJECT ICS
                fprintf('%s) Rejecting EEG independent components...\n',subj_chars{subj_i});
                tmp_rej_struct.plot_save_dir = tmp_save_dir;
                %- get ic component criteria
                [ic_rej_out] = mim_get_reject_subj_ics(EEG, ...
                    'REJ_STRUCT',tmp_rej_struct)
                par_save(ic_rej_out,tmp_save_dir,sprintf('%s_ICRej.mat',EEG.subject));
                %- perform rejections and log them
                [EEG,rej_struct_out,rmv_subj_flag] = mim_reject_subj_ics(EEG,ic_rej_out, ...
                    'CHK_STRUCT',tmp_chk_struct);
                if ~rmv_subj_flag
                    tmp_rej_crit_out{subj_i} = rej_struct_out;
                    par_save(rej_struct_out,[EEG.filepath filesep 'rej_struct_out.mat']);
                    % ALLEEG{subj_i} = pop_saveset(EEG, ...
                    %         'filepath',tmp_save_dir, ...
                    %         'filename',eeg_fnames{subj_i}, ...
                    %         'savemode','twofiles')
                    ALLEEG{subj_i} = EEG;
                else
                    tmp_rej_crit_out{subj_i} = rej_struct_out;
                    fprintf('%s) Subject removed...\n',subj_chars{subj_i});
                end
            else
                ALLEEG{subj_i} = pop_loadset('filepath',tmp_save_dir, ...
                    'filename',eeg_fnames{subj_i})
            end
        catch e
            fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
            % exit();
        end
    end
    %- save reject table
    tmp_rej_crit_out = tmp_rej_crit_out(~cellfun(@isempty,tmp_rej_crit_out));
    tmp_rej_crit_out = util_resolve_struct(tmp_rej_crit_out);
    tmp_rej_crit_out = struct2table(tmp_rej_crit_out);
    writetable(tmp_rej_crit_out,[save_dir filesep 'rejection_crit_table.xlsx']);

    %##
    ALLEEG = ALLEEG(~cellfun(@isempty,ALLEEG));
    ALLEEG = util_resolve_struct(ALLEEG);
    %##
    [STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                    'updatedat','off',...
                                    'savedat','off',...
                                    'name',STUDY_FNAME_CONT,...
                                    'filename',STUDY_FNAME_CONT,...
                                    'filepath',save_dir);
    [STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
    [STUDY,~] = parfunc_save_study(STUDY,ALLEEG,...
                                    STUDY.filename,STUDY.filepath,...
                                    'RESAVE_DATASETS','on');  
    
else
    if ~ispc
        [STUDY,~] = pop_loadstudy('filename',[STUDY_FNAME_CONT '_UNIX.study'],'filepath',save_dir);
    else
        [STUDY,~] = pop_loadstudy('filename',[STUDY_FNAME_CONT '.study'],'filepath',save_dir);
    end
end
%## UNIT TESTING
%{
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME) filesep 'iclabel_cluster_kmeansalt_rb3'];
cluster_study_dir = [cluster_dir filesep 'icrej_5'];
if ~ispc
    tmp = load('-mat',[cluster_study_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_dir filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
%}
%## RECOVER REJ STRUCT XSLX
tmp_rej_crit_out = cell(1,length(subj_chars)); 
for subj_i = 1:length(subj_chars)  
    
end
%- save reject table
tmp_rej_crit_out = tmp_rej_crit_out(~cellfun(@isempty,tmp_rej_crit_out));
tmp_rej_crit_out = util_resolve_struct(tmp_rej_crit_out);
tmp_rej_crit_out = struct2table(tmp_rej_crit_out);
writetable(tmp_rej_crit_out,[save_dir filesep 'rejection_crit_table.xlsx']);
%% INITIALIZE PARFOR LOOP VARS
%##
subj_chars = {STUDY.datasetinfo.subject};
eeg_fpaths = {STUDY.datasetinfo.filepath};
eeg_fnames = {STUDY.datasetinfo.filename};
ALLEEG = cell(1,length(STUDY.datasetinfo));
alleeg_fpaths = cell(length(STUDY.datasetinfo),1);
%- clear vars for memory
% clear MAIN_ALLEEG
%% UNIT TEST ======================================================== %%
%{
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
fPath = [epoched_fPath filesep [DEF_EPOCH_PA
rmv_subj = zeros(1,length(STUDY.datasetinfo));RAMS.gait_trial_chars{:}]];
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
%}
%% GENERATE EPOCH MAIN FUNC
%## PARFOR LOOP
parfor subj_i = 1:length(subj_chars)
    %## LOAD EEG DATA
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    % EEG = ALLEEG(subj_i);
    EEG = pop_loadset('filepath',eeg_fpaths{subj_i},'filename',eeg_fnames{subj_i});
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
    epoch_fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
    epoch_fname = sprintf('%s_%s_epochd.set',EEG.subject,[tmp_epoch_params.gait_trial_chars{:}]);
    if ~exist(epoch_fpath,'dir')
        mkdir(epoch_fpath)
    end
    %- parse
    try
        %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
        for i = 1:length(EEG)
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
        %## EPOCHsubj_chars
        [EEG,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',tmp_epoch_params);
        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_ALLEEG
            for i = 1:length(EEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,EEG(i).condition)];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(EEG(i),'savemode','twofiles',...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],EEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],EEG(i).condition);
            end
            alleeg_fpaths{subj_i} = cond_files;
        end
        EEG = pop_mergeset(EEG,1:length(EEG),1);
        EEG.etc.cond_files = cond_files;
        %## timewarp for across condition
        if strcmp(tmp_epoch_params.epoch_method,'timewarp')
            timewarp = make_timewarp(EEG,tmp_epoch_params.tw_events,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(EEG.epoch),goodepochs);
            %- reject outlier strides
            EEG = pop_select(EEG,'notrial',sedi);
            %- store timewarp structure in EEG
            EEG.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            EEG.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            EEG.urevent = []; % might be needed
            EEG.etc.epoch.epoch_limits = tmp_epoch_params.epoch_time_lims;
        end
        %## STRUCT EDITS
        EEG.urevent = []; % might be needed
        EEG.etc.epoch.epoch_limits = tmp_epoch_params.epoch_time_lims;
        %- checks
        EEG = eeg_checkset(EEG,'eventconsistency');
        EEG = eeg_checkset(EEG);
        EEG = eeg_checkamica(EEG);
        %- save
        [EEG] = pop_saveset(EEG,'savemode','twofiles',...
                'filename',epoch_fname,...
                'filepath',epoch_fpath,...
                'version','6');
        ALLEEG{subj_i} = EEG;
    catch e
        ALLEEG{subj_i} = []; %EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
    end
end
%% SAVE BIG STUDY
fprintf('==== Saving Epoch''d Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s',strjoin(subj_chars(cellfun(@isempty,ALLEEG)),','));
ALLEEG = ALLEEG(~cellfun(@isempty,ALLEEG));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
ALLEEG = util_resolve_struct(ALLEEG);
%##
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_EPOCH,...
                                'filename',STUDY_FNAME_EPOCH,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
%- epoch params
STUDY.etc.a_epoch_process.epoch_params = DEF_EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
%## TIMEWARP TIMINGS
warptos = nan(length(ALLEEG),length(ALLEEG(1).etc.timewarp_by_cond(1).warpto));
%- subject based timings
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    warptos(subj_i,:) = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
avg_warpto_events = floor(nanmean(warptos));
STUDY.etc.a_epoch_process.avg_warpto_events = avg_warpto_events;
%-
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

