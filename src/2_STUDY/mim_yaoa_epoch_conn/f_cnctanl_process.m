%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_epoch_conn/run_f_cnctanl_process.sh

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
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- connecitivty modeling
SRATE = 500;
CONN_FREQS = (4:50);
CONN_METHODS = {'dDTF08'}; %{'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
WINDOW_LENGTH = 0.4;
VERBOSITY_LEVEL = 1;
MORDER = 70;
WINDOW_STEP_SIZE = (1/SRATE)*MORDER;
%-
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
RESAMLE_FREQ = 225;
DO_BOOTSTRAP = false;
DO_PHASE_RND = false;
DO_STANDARD_TRIALS = false;
MIN_CONN_TRIALS = 100;
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
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%## PREPARE DATA PARAMS
DEF_PREPDATA = struct('VerbosityLevel',VERBOSITY_LEVEL,...
             'SignalType',{{'Components'}},...
             'VariableNames',[],...
             'Detrend',{{'verb',VERBOSITY_LEVEL,'method',{'linear'},...
                    'piecewise',{'seglength',0.33,'stepsize',0.0825},...
                    'plot',false}},...
             'NormalizeData',{{'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}}},...
             'resetConfigs',true,...
             'badsegments',[],...
             'newtrials',[],...
             'equalizetrials',false);
%## ESTIAMTE MODEL ORDER PARAMS
% DEF_ESTSELMOD = struct('modelingApproach',{{'Segmentation VAR',...
%                         'algorithm',{{'Vieira-Morf'}},...
%                         'winStartIdx',[],...
%                         'winlen',WINDOW_LENGTH,...
%                         'winstep',WINDOW_STEP_SIZE,...
%                         'taperfcn','rectwin',...
%                         'epochTimeLims',[],...
%                         'prctWinToSample',100,...
%                         'normalize',[],...
%                         'detrend',{'method','linear'},...
%                         'verb',VERBOSITY_LEVEL}},...
%     'morderRange',[40, 70],...
%     'downdate',true,...
%     'runPll',[],...
%     'icselector',{{'aic','hq'}},...
%     'winStartIdx',[],...
%     'epochTimeLims',[],...
%     'prctWinToSample',80,...
%     'plot',[],...
%     'verb',VERBOSITY_LEVEL);
DEF_ESTSELMOD = struct('modelingApproach',{{'Segmentation VAR',...
                        'algorithm',{{'Vieira-Morf'}},...
                        'winStartIdx',[],...
                        'winlen',WINDOW_LENGTH,...
                        'winstep',WINDOW_STEP_SIZE,...
                        'taperfcn','rectwin',...
                        'epochTimeLims',[],...
                        'prctWinToSample',100,...
                        'normalize',[],...
                        'detrend',{'method','linear'},...
                        'verb',VERBOSITY_LEVEL}},...
    'morderRange',[1,floor(SRATE*(WINDOW_LENGTH/2)-1)],...
    'downdate',true,...
    'RunInParallel',{{'profile','local',...
        'numWorkers',SLURM_POOL_SIZE}},...
    'icselector',{{'aic','hq'}},...
    'winStartIdx',[],...
    'epochTimeLims',[],...
    'prctWinToSample',80,...
    'plot',[],...
    'verb',VERBOSITY_LEVEL);
%##
DEF_PLOTORDERCRIT = struct('conditions',{DEF_EPOCH_PARAMS.gait_trial_chars},    ...
                            'icselector',{DEF_ESTSELMOD.icselector},  ...
                            'minimizer',{{'min'}}, ...
                            'prclim', 90);
%##
DEF_ESTDISPMVAR_CHK = struct('morder',MORDER,...
        'winlen',WINDOW_LENGTH,'winstep',WINDOW_STEP_SIZE,'verb',1);
%##
DEF_ESTFITMVAR = struct('connmethods',{CONN_METHODS}, ...
            'absvalsq',true,           ...
            'spectraldecibels',true,   ...
            'freqs',CONN_FREQS,        ...
            'verb',1);
%## PATHING
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
%- study dir & conn figure dir
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
conn_fig_dir = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'conn_valid_epoch'];
%- load study file
STUDY_FNAME_LOAD = 'all_comps_study';
STUDY_FNAME_SAVE = 'epoch_conn_study';
study_fpath = [studies_fpath filesep sprintf('%s',study_dir_name)];
%- load cluster
CLUSTER_K = 12;
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%- create new study directory
if ~exist(conn_fig_dir,'dir')
    mkdir(conn_fig_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
% if ~exist([study_fpath filesep STUDY_FNAME_LOAD '.study'],'file')
%     error('ERROR. study file does not exist');
%     exit(); %#ok<UNRCH>
% else
%     %## LOAD STUDY
%     if ~ispc
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_LOAD '_UNIX.study'],'filepath',study_fpath);
%     else
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_LOAD '.study'],'filepath',study_fpath);
%     end
%     cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
%     MAIN_STUDY.cluster = cl_struct;
%     [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[study_fpath filesep sprintf('%s_UNIX.study',STUDY_FNAME_LOAD)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[study_fpath filesep sprintf('%s.study',STUDY_FNAME_LOAD)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,~,~,valid_cls] = eeglab_get_cluster_comps(STUDY);
% %% MIM STUDY STATS
% %## FIND MINIMUM TRIALS ACROSS COHORT
% if DO_STANDARD_TRIALS
%     subj_chars = {MAIN_STUDY.datasetinfo.subject}';
%     trial_counts = cell(length(subj_chars),length(TRIAL_TYPES));
%     channel_rejs = cell(length(subj_chars),1);
%     sample_rej_perc = cell(length(subj_chars),1);
%     for subj_i=1:length(subj_chars)
%         EEG = MAIN_ALLEEG(subj_i);
%         for cond_i = 1:length(TRIAL_TYPES)
%             tmp = sum(strcmp({MAIN_STUDY.datasetinfo(subj_i).trialinfo.cond},TRIAL_TYPES{cond_i}));
%             trial_counts{subj_i,cond_i} = tmp;
%         end
%         channel_rejs{subj_i,1} = sum(~EEG.etc.channel_mask(strcmp({EEG.urchanlocs.type},'EEG')));
%         sample_rej_perc{subj_i,1} = (length(EEG.etc.clean_sample_mask)-sum(EEG.etc.clean_sample_mask))/length(EEG.etc.clean_sample_mask);
%     end
%     tbl_out = table(subj_chars,trial_counts(:,1),trial_counts(:,2),trial_counts(:,3),...
%         trial_counts(:,4),trial_counts(:,5),trial_counts(:,6),trial_counts(:,7),trial_counts(:,8),channel_rejs,sample_rej_perc,'VariableNames',{'subj_chars','0p25','0p5','0p75','1p0','flat','low','med','high','channel_rejs','sample_rej_perc'});
%     vals = tbl_out{:,2:9};
%     trial_mins = zeros(size(vals,2),1);
%     trial_maxs = zeros(size(vals,2),1);
%     for i = 1:size(vals,2)
%         trial_mins(i) = min([vals{:,i}]);
%         trial_maxs(i) = max([vals{:,i}]);
%     end
%     MIN_CONN_TRIALS = min(trial_mins);
% end
%% INITIALIZE PARFOR LOOP VARS
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
LOOP_VAR = 1:length(STUDY.datasetinfo);
tmp = cell(1,length(STUDY.datasetinfo));
rmv_subj = zeros(1,length(STUDY.datasetinfo));
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica',1);
%## PARFOR LOOP
EEG = [];
%## CASE DEBUG
%{
fPaths = {'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01\H1004\GAIT_EPOCHED\cond_0p5'};
fNames = {'cond_0p5.set'};
%}
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% for subj_i = LOOP_VAR
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %## LOAD EEG DATA
    eeg_savefpath = [fPaths{subj_i} filesep 'conn_epoch'];
    if ~exist(eeg_savefpath,'dir')
        mkdir(eeg_savefpath)
    end
    if ~exist([eeg_savefpath filesep fNames{subj_i}],'file')
        EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); 
            EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
        end
        %## RESAMPLE
        %- (04/24/2024) JS, the rationale for resampling is majorly based
        %on computation times, however, a argument could be made that by
        %removing higher frequency (i.e., downsampling) data the
        %connectivity model may fit more accurately to lower frequency
        %content. The latter point has not been proven, merely suggested in
        %the literature.
        %(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10044923/)
        % EEG = pop_resample(EEG,RESAMLE_FREQ);
        %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
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
        %## EPOCH
        [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,'EPOCH_PARAMS',DEF_EPOCH_PARAMS);
        for i = 1:length(ALLEEG)
            if DO_STANDARD_TRIALS
                %##
                inds = randi(length(ALLEEG(i).epoch),MIN_CONN_TRIALS,1);
                %##
                tmp_all = ALLEEG(i);
                tmp_all = pop_selectevent(tmp_all,'event',inds);
                while length(tmp_all.epoch)~=MIN_CONN_TRIALS
                    tmp_all = ALLEEG(i);
                    inds = randi(length(tmp_all.epoch),MIN_CONN_TRIALS,1)
                    tmp_all = pop_selectevent(tmp_all,'event',inds);
                end
                ALLEEG(i) = tmp_all;
            end
        end
        EEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        %- checks
        % EEG = eeg_checkset(EEG,'eventconsistency');
        % EEG = eeg_checkset(EEG);
        % EEG = eeg_checkamica(EEG);
        % %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); 
            EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
        end
        % pop_editoptions('option_computeica',1);
        try
            %## RUN MAIN_FUNC
            [TMP,t_out] = cnctanl_sift_pipe(EEG,components,conn_fig_dir,...
                'DO_PHASE_RND',DO_PHASE_RND,...
                'DO_BOOTSTRAP',DO_BOOTSTRAP,...
                'PREPDATA',DEF_PREPDATA,...
                'ESTSELMOD',DEF_ESTSELMOD,...
                'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
                'ESTFITMVAR',DEF_ESTFITMVAR,...
                'PLOTORDERCRIT',DEF_PLOTORDERCRIT);
    %         for i=1:length(TMP)
    %             par_save(TMP(i).CAT)
    %         end
            BIG_CAT = cat(1,TMP(:).CAT);
            EEG.etc.COND_CAT = BIG_CAT;
            EEG.etc.conn_table = t_out;
            EEG.etc.conn_meta.comps_out = comps_out;
            fName = strsplit(EEG.filename,'.'); fName = [fName{1} '.mat'];
            par_save(t_out,eeg_savefpath,fName,'_conntable');
            [EEG] = pop_saveset(EEG,...
                'filepath',eeg_savefpath,'filename',EEG.filename,...
                'savemode','twofiles');
            tmp{subj_i} = EEG;
        catch e
            rmv_subj(subj_i) = 1;
            EEG.CAT = struct([]);
            tmp{subj_i} = EEG;
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n',...
                     'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
            disp(e);
        end
    end
end
pop_editoptions('option_computeica',0);
%% SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s',MAIN_ALLEEG(cellfun(@isempty,tmp)).subject);
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
tmp = eeg_checkset(tmp,'eventconsistency');
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_SAVE,...
                                'filename',STUDY_FNAME_SAVE,...
                                'filepath',study_fpath);
%## ASSIGN PARAMETERS
STUDY.etc.a_epoch_process.epoch_chars = TRIAL_TYPES;
STUDY.etc.d_cnctanl_process.params = struct('DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'PREPDATA',DEF_PREPDATA,...
            'ESTSELMOD',DEF_ESTSELMOD,...
            'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
            'ESTFITMVAR',DEF_ESTFITMVAR,...
            'PLOTORDERCRIT',DEF_PLOTORDERCRIT,...
            'comps',comps_out);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY_FNAME_SAVE,study_fpath,...
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

