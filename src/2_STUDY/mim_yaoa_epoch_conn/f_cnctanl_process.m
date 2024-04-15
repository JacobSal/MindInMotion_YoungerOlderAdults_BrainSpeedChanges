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
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = SCRIPT_DIR;
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(SCRIPT_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- connecitivty modeling
CONN_FREQS = (4:50);
% (01/19/2024) JS, trying 4:50 from 1:100(see. Steven Peterson 2019 NeuroImage)
FREQ_BANDS = {CONN_FREQS;4:8;8:13;13:28};
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
CONN_METHODS = {'dDTF08','S'}; %{'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
WINDOW_LENGTH = 0.4;
RESAMLE_FREQ = 225;
% (01/19/2024) JS, trying 0.4 from 0.5 (see. Steven Peterson 2019 NeuroImage)
WINDOW_STEP_SIZE = 0.02;
% (01/19/2024) JS, trying 0.02 from 0.025 (see. Steven Peterson 2019 NeuroImage)
DO_BOOTSTRAP = false;
% (01/19/2024) JS, unsure to turn to false for quicker estimates?
DO_PHASE_RND = false;
DO_STANDARD_TRIALS = false;
MORDER = 32;
%- datetime override
study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
%## soft define
studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
conn_save_dir = [studies_dir filesep sprintf('%s',study_dir_name) filesep 'conn_valid'];
%- load study file
STUDY_FNAME = 'all_comps_study';
STUDY_FPATH = [studies_dir filesep sprintf('%s',study_dir_name)];
STUDY_FNAME_SAVE = 'epoch_conn_study';
%- load cluster
CLUSTER_DIR = [studies_dir filesep sprintf('%s',study_dir_name) filesep 'cluster'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
CLUSTER_STUDY_DIR = [CLUSTER_DIR filesep 'icrej_5'];
CLUSTER_K = 12;
%- create new study directory
if ~exist(conn_save_dir,'dir')
    mkdir(conn_save_dir);
end
%% LOAD EPOCH STUDY

%- Create STUDY & ALLEEG structs
if ~exist([STUDY_FPATH filesep STUDY_FNAME '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
    %## LOAD STUDY
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',STUDY_FPATH);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',STUDY_FPATH);
    end
    cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
    MAIN_STUDY.cluster = cl_struct;
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
end
%% MIM STUDY STATS
%## FIND MINIMUM TRIALS ACROSS COHORT
if DO_STANDARD_TRIALS
    subj_chars = {MAIN_STUDY.datasetinfo.subject}';
    trial_counts = cell(length(subj_chars),length(TRIAL_TYPES));
    channel_rejs = cell(length(subj_chars),1);
    sample_rej_perc = cell(length(subj_chars),1);
    for subj_i=1:length(subj_chars)
        EEG = MAIN_ALLEEG(subj_i);
        for cond_i = 1:length(TRIAL_TYPES)
            tmp = sum(strcmp({MAIN_STUDY.datasetinfo(subj_i).trialinfo.cond},TRIAL_TYPES{cond_i}));
            trial_counts{subj_i,cond_i} = tmp;
        end
        channel_rejs{subj_i,1} = sum(~EEG.etc.channel_mask(strcmp({EEG.urchanlocs.type},'EEG')));
        sample_rej_perc{subj_i,1} = (length(EEG.etc.clean_sample_mask)-sum(EEG.etc.clean_sample_mask))/length(EEG.etc.clean_sample_mask);
    end
    tbl_out = table(subj_chars,trial_counts(:,1),trial_counts(:,2),trial_counts(:,3),...
        trial_counts(:,4),trial_counts(:,5),trial_counts(:,6),trial_counts(:,7),trial_counts(:,8),channel_rejs,sample_rej_perc,'VariableNames',{'subj_chars','0p25','0p5','0p75','1p0','flat','low','med','high','channel_rejs','sample_rej_perc'});
    vals = tbl_out{:,2:9};
    trial_mins = zeros(size(vals,2),1);
    trial_maxs = zeros(size(vals,2),1);
    for i = 1:size(vals,2)
        trial_mins(i) = min([vals{:,i}]);
        trial_maxs(i) = max([vals{:,i}]);
    end
    MIN_CONN_TRIALS = min(trial_mins);
end
%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% for subj_i = LOOP_VAR
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %## LOAD EEG DATA
    %{
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    eeg_savefpath = [fPaths{subj_i} filesep 'conn'];
    if ~exist(eeg_savefpath,'dir')
        mkdir(eeg_savefpath)
    end
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
    end
    fprintf('%s) Processing componets:\n',EEG.subject)
    fprintf('%i,',components'); fprintf('\n');
    %- re-epoch
    ALLEEG = cell(1,length(TRIAL_TYPES));
    for i = 1:length(EEG.etc.cond_files)
        if ~ispc
            ALLEEG{i} = pop_loadset('filepath',convertPath2UNIX(EEG.etc.cond_files(i).fPath),'filename',EEG.etc.cond_files(i).fName);
        else
            ALLEEG{i} = pop_loadset('filepath',convertPath2Drive(EEG.etc.cond_files(i).fPath),'filename',EEG.etc.cond_files(i).fName);
        end
        if DO_STANDARD_TRIALS
            %##
            inds = randi(length(EEG.epoch),MIN_CONN_TRIALS,1);
            %##
            tmp_all = ALLEEG{i};
            tmp_all = pop_selectevent(tmp_all,'event',inds);
            while length(tmp_all.epoch)~=MIN_CONN_TRIALS
                tmp_all = ALLEEG{i};
                inds = randi(length(tmp_all.epoch),MIN_CONN_TRIALS,1)
                tmp_all = pop_selectevent(tmp_all,'event',inds);
            end
            ALLEEG{i} = tmp_all;
        end
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1)
    %}
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    eeg_savefpath = [fPaths{subj_i} filesep 'conn'];
    if ~exist(eeg_savefpath,'dir')
        mkdir(eeg_savefpath)
    end
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
    end
    EEG = pop_resample(EEG,RESAMLE_FREQ);
    try
        %## RUN MAIN_FUNC
        % [TMP,t_out] = cnctanl_sift_pipe(ALLEEG,components,CONN_METHODS,conn_save_dir,...
        %     'MORDER',MORDER,...
        %     'DO_PHASE_RND',DO_PHASE_RND,...
        %     'DO_BOOTSTRAP',DO_BOOTSTRAP,...
        %     'FREQS',CONN_FREQS,...
        %     'WINDOW_LENGTH',WINDOW_LENGTH,... 
        %     'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,...
        %     'GUI_MODE','nogui',...
        %     'VERBOSITY_LEVEL',1,...
        %     'ESTSELMOD_CFG',[]);
        %## RUN MAIN_FUNC
        [TMP,t_out] = cnctanl_sift_pipe(EEG,components,CONN_METHODS,conn_save_dir,...
            'MORDER',MORDER,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'FREQS',CONN_FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,... 
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,...
            'GUI_MODE','nogui',...
            'VERBOSITY_LEVEL',1,...
            'ESTSELMOD_CFG',[]);
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
                                'filepath',STUDY_FPATH);
%## ASSIGN PARAMETERS
STUDY.etc.a_epoch_process.epoch_chars = TRIAL_TYPES;
STUDY.etc.d_cnctanl_process.params = comps_out;
STUDY.etc.d_cnctanl_process.params = struct('CONN_METHODS',{CONN_METHODS},...
            'MORDER',MORDER,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'FREQS',CONN_FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,... 
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,...
            'GUI_MODE','nogui',...
            'VERBOSITY_LEVEL',1,...
            'ESTSELMOD_CFG',[]);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY_FNAME_SAVE,STUDY_FPATH,...
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

