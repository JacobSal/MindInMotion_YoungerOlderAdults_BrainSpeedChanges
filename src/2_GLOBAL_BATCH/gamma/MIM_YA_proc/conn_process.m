%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_YA_proc/run_conn_process.sh

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
if DO_UNIX
%     eeg_options;
    pop_editoptions('option_parallel',1,'option_storedisk',1);
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
    fprintf('Number of workers: %i',pp.NumWorkers);
    fprintf('Number of threads: %i',pp.NumThreads);
    %- make meta data directory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions('option_parallel',0,'option_storedisk',1);
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
% SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
%     'H1022','H1024','H1026','H1027','H1033','H1034'}; % CHANG,LIU(12/26/2022)
SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)
% SUBJ_MISSING_EVENT = {'H1002','H1048'};
TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
%- Subject Picks
SUBJ_PICS = {SUBJ_YNG};
SUBJ_ITERS = {1:length(SUBJ_YNG)}; % CHANG,LIU(12/26/2023)
% SUBJ_ITERS = {[],1:length(SUBJ_HMA),1:length(SUBJ_NMA)}; % CHANG,LIU(02/15/2023)
GROUP_NAMES = {'H1000s'};
%- Subject Directory Information
PREPROCESS_NAME = 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10'; % CHANG,LIU(02/15/2023)
OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-HY_202212';  % CHANG,LIU(02/15/2023)
% PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % CHANG,LIU(12/26/2022)
% OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % CHANG,LIU(12/26/2022)
if DO_UNIX
    OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR);
else
    OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR);
end
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping
%% ===================================================================== %%
%## PROCESSING PARAMS
% pop_editoptions('option_parallel',1);
%- study group and saving
SAVE_EEG = false; % saves EEG structures throughout processing
%- component rejection crit
THRESHOLD_DIPFIT_RV = 0.15;
THRESHOLD_BRAIN_ICLABEL = 0.60;
%- MIM specific epoching
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
%- subj_i_epoch
PARSE_TYPE = 'Constant'; %
EPOCH_TIME_LIMITS = [-1,3]; %[-2,2];
TRIAL_LENGTH = 3*60; % trial length in seconds
PER_OVERLAP = 0.0; % percent overlap between epochs
TRIAL_BEGIN_STR = 'TrialStart';
TRIAL_END_STR = 'TrialEnd';
EVENT_TRIAL_PARSER = 'type';
EVENT_COND_PARSER = 'cond';
%- connectivity process
CONN_FREQS = (1:100);
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
%- subj_i_genConnMeas
CNCTANL_TOOLBOX = 'sift'; %'bsmart'
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
DO_BOOTSTRAP = true;
% ASSIGN_BOOTSTRAP_MEAN = false;
% SAVE_CONN_BOOTSTRAP = false;
%- subj_i_genConnStats
DO_PHASE_RND = true;
%- eeglab_cluster.m spectral params
FREQ_LIMITS = [1,100];
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
PAD_RATIO = 2;
%% global script chain (VERSION 1)
%- datetime override
% dt = '26122023'; % YA 
% dt = '24022023_YA'; % YA (MORE SUBJECTS)
dt = '23032023_YA_fem'; % YA (03/23/2023)
% dt = '02022023';
%## PATH & TEMP STUDY NAME
%- hard define
SESSION_NUMBER = '1';
study_fName_1 = sprintf('copy_study');
study_fName_2 = 'reduced_comps_study';
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
dip_chn_fPath   = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fName   = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep PREPROCESS_NAME filesep SUBJ_PICS{group_i}{subj_i}];
        fNames{cnt} = [SUBJ_PICS{group_i}{subj_i} '_cleanEEG_' PREPROCESS_NAME '_' OUTSIDE_DATA_SUFFIX '.set'];
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file')>1)
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
%% CREATE STUDY
%- NOTE: need a dipole fit before this step can be done.
%- params
%## Create STUDY & ALLEEG structs
if ~exist([save_dir filesep study_fName_1 '.study'],'file') %|| true
    fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
    [ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,...
                        conditions,groups,sessions);
    %- NOTES: spec_mode, 'psd' does not work.
    [MAIN_ALLEEG,MAIN_STUDY] = eeglab_cluster(ALLEEG,...
                        study_fName_1,save_dir,...
                        'SPEC_MODE',SPEC_MODE,...
                        'FREQ_LIMITS',FREQ_LIMITS,...
                        'CYCLE_LIMITS',CYCLE_LIMITS,...
                        'FREQ_FAC',FREQ_FAC,...
                        'PAD_RATIO',PAD_RATIO);
    fprintf(1,'\n==== DONE: CLUSTERING SUBJECT DATA ====\n');
end
%% REDUCE CLUSTERS TO ONE COMPONENT PER SUBJECT
if ~exist([save_dir filesep study_fName_2 '.study'],'file')
    %## LOAD PREVIOUS STUDY
    if ~exist('MAIN_STUDY','var') || ~exist('MAIN_ALLEEG','var')
        fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
        if ~ispc
            [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',save_dir);
        else
            [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',save_dir);
        end
        fprintf(1,'\n==== DONE: LOADING CLUSTER STUDY DATA ====\n');
    end
    %## RECALCULATE ICAACT MATRICES
    MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
    for subj_i = 1:length(MAIN_ALLEEG)
        if isempty(MAIN_ALLEEG(subj_i).icaact)
            fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
            MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
        end
    end
    %## ADD ANATOMICAL LABELS
    [~, atlas_cell] = add_anatomical_labels(MAIN_STUDY,MAIN_ALLEEG);
    MAIN_STUDY.etc.add_anatomical_labels = atlas_cell;
    %## PCA reduction algorithm
%     [MAIN_STUDY,MAIN_ALLEEG,comps_out,outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
    %## ICA reduction algorithm
    [MAIN_STUDY, comps_out] = cluster_ica_reduce(MAIN_STUDY);
    %## SAVE STUDY
    [MAIN_STUDY,MAIN_ALLEEG] = eeglab_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                                            study_fName_2,save_dir);
else
    fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_2 '_UNIX.study'],'filepath',save_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_2 '.study'],'filepath',save_dir);
    end
    %## Extract components for each cluster & subject
    comps_out = zeros(length(MAIN_STUDY.cluster),length(MAIN_ALLEEG));
    for clus_i = 2:length(MAIN_STUDY.cluster)
        sets_i = MAIN_STUDY.cluster(clus_i).sets;
        for j = 1:length(sets_i)
            comps_out(clus_i,sets_i(j)) = MAIN_STUDY.cluster(clus_i).comps(j);
        end
    end
end

%% PLOT DIPOLES, PSD'S, && TOPOPLOTS
plot_fNames = {'allDipPlot_0','allDipPlot_1','allDipPlot_2','allSpecPlot','allTopoPlot'};
plot_chk = cellfun(@(x) ~exist([save_dir filesep sprintf('%s.jpg',x)],'file'),plot_fNames);
if any(plot_chk) && false
    fprintf('==== Making Dipole Plots ====\n');
    std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),...
                'mode','multicolor');
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_0.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_0.jpg')]);
    view([45,0,0])
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.jpg')]);
    view([0,-45,0])
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.jpg')]);
    close fig_i
    %- Spec plot
    fprintf('==== Making Spectogram Plots ====\n');
    [MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                        'components',...                               
                                        'recompute','on',...
                                        'spec','on',...
                                        'specparams',...
                                        {'specmode',SPEC_MODE,'freqfac',4,...
                                        'freqrange',FREQ_LIMITS});
    specMin = 10;
    specMax = 45;
    std_specplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.jpg')]);
    close fig_i
    %- Topo plot
    fprintf('==== Making Topograph Plots ====\n');
    [MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                        'components',...                               
                                        'recompute','on',...
                                        'scalp','on');
    std_topoplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster));
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.jpg')]);
    close fig_i
    %## POP VIEW PROPS
    if ~exist([save_dir filesep 'component_props'],'dir')
        mkdir([save_dir filesep 'component_props']);
        for cluster_i = 2:length(MAIN_STUDY.cluster)
            sets_clust = MAIN_STUDY.cluster(cluster_i).sets;
            for i = 1:length(sets_clust)
                subj_i = sets_clust(i);
                comps_clust = MAIN_STUDY.cluster(cluster_i).comps(i);
                hold on;
                pop_prop_extended(MAIN_ALLEEG(subj_i),0,comps_clust,NaN,...
                {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
                fig = gcf;
                fig.Position = [500 300 920 480]; 
                hold off;
                savefig(fig,[save_dir filesep 'component_props' filesep sprintf('%i_%s_viewprops_co%i_cl%i.fig',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust,cluster_i)]);
                close fig
            end
        end
    end
end
%% INITIALIZE PARFOR LOOP VARS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(MAIN_ALLEEG)]);
else
    POOL_SIZE = 1;
end
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
fPaths_out = cell(1,length(MAIN_ALLEEG));
fNames_out = cell(1,length(MAIN_ALLEEG));
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
%- HARD DEFINES
SUFFIX_PATH_EPOCHED = 'EPOCHED';
%- clear vars for memory
% clear MAIN_ALLEEG
%% GENERATE EPOCH MAIN FUNC
%## PARFOR LOOP
parfor (subj_i = LOOP_VAR,POOL_SIZE)
% for subj_i = LOOP_VAR
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    fprintf('Running subject %s\n',EEG.subject)
    %## INITIATLE CHECK
    chk = zeros(length(TRIAL_TYPES),1);
    for trial_i = 1:length(TRIAL_TYPES)
        out_fPath = [save_dir filesep EEG.subject filesep,...
                        SUFFIX_PATH_EPOCHED filesep sprintf('cond_%s',TRIAL_TYPES{trial_i})];
        out_fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,TRIAL_TYPES{trial_i});
        if exist([out_fPath filesep out_fName],'file')
            chk(trial_i) = 1;
        end
    end
    if all(chk)
        fprintf('Skipping epoching for subject %s...',EEG.subject)
        continue
    end
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end
    %## PARSE TRIALS
    for trial_i = 1:length(TRIAL_TYPES)
        events = mim_trial_event_override(EEG.subject,TRIAL_TYPES{trial_i},TRIAL_OVERRIDE_FPATH);
        if ~isempty(events)
            f = fields(EEG.event);
            for e = 1:length(events)
                for f_i = 1:length(f)
                    if ~isfield(events{e},f{f_i})
                        events{e}.(f{f_i}) = [];
                    end
                end
            end
            events = [events{:}];
            fprintf('appending events...\n');
            EEG.event = [EEG.event, events];
        end
    end
    %- let EEGLAB rearrange the event order
    EEG = eeg_checkset(EEG,'eventconsistency');
    %- parse
    try
        ALLEEG = mim_parse_trials(EEG,TRIAL_TYPES,EPOCH_TIME_LIMITS);
        %## Save EEG
        for trial_i = 1:length(ALLEEG)
            out_fPath = [save_dir filesep ALLEEG(trial_i).subject filesep SUFFIX_PATH_EPOCHED filesep sprintf('cond_%s',TRIAL_TYPES{trial_i})];
            if ~exist(out_fPath,'dir')
                mkdir(out_fPath);
            end
            fprintf(1,'Saving Subject %s\n',ALLEEG(trial_i).subject); 
            [ALLEEG(trial_i)] = pop_saveset(ALLEEG(trial_i),'savemode','twofiles',...
                'filename',ALLEEG(trial_i).filename,...
                'filepath',out_fPath);
        end
        fNames_out{subj_i} = {ALLEEG.filename};
        fPaths_out{subj_i} = {ALLEEG.filepath};
        tmp{subj_i} = ALLEEG;
%         fields_EEG = fields(ALLEEG);
    catch e
        rmv_subj(subj_i) = 1;
        EEG.timewarp = struct([]);
        tmp{subj_i} = EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
    end
end
%% CONNECTIVITY MAIN FUNC
%## PARFOR LOOP
EEG = [];
LOOP_VAR = find(~rmv_subj);
parfor (subj_i = LOOP_VAR,POOL_SIZE)
% for subj_i = 31
    %## INITIATE VARS
    ALLEEG = cell(1,length(TRIAL_TYPES));
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %## LOAD EEG DATA
    for trial_i = 1:length(TRIAL_TYPES)
        EEG = pop_loadset('filepath',fPaths_out{subj_i}{trial_i},'filename',fNames_out{subj_i}{trial_i});
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        ALLEEG{trial_i} = EEG;
    end
    ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
    fprintf('Running subject %s\n',ALLEEG(1).subject)
    %## RUN MAIN_FUNC
%     fPath_chk = [save_dir filesep EEG.subject filesep,...
%                 SUFFIX_PATH_EPOCHED filesep];
%     chk = cnctanl_chk_CAT(EEG,TRIAL_TYPES,fPaths_out{subj_i}{trial_i},...
%             'CHK_BOOTSTRAP',DO_BOOTSTRAP,...
%             'CHK_PHASE_RND',DO_PHASE_RND);
    if true %any(chk)
        [ALLEEG] = main_func_v2(ALLEEG,TRIAL_TYPES,save_dir,components,...
            'SAVE_EEG',SAVE_EEG,...
            'CONN_METHODS',CONN_METHODS,... %
            'FREQS',CONN_FREQS,...
            'CNCTANL_TOOLBOX',CNCTANL_TOOLBOX,... 
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE);
    end
    tmp{subj_i} = ALLEEG;
end

%% SAVE CONDITION STUDIES
fprintf('==== SAVING STUDIES ====\n');
for cond_i = 1:length(TRIAL_TYPES)
    tmp_alleeg = cell(1,length(tmp));
    study_fName = sprintf('%s_MIM_study',TRIAL_TYPES{cond_i});
    
    %- extract each condition for each subject
    for subj_i = 1:length(tmp)
        if length(tmp{subj_i}) == 1
            EEG = tmp{subj_i};
        else
            EEG = tmp{subj_i}(cond_i);
        end
        tmp_alleeg{subj_i} = EEG;
    end
    tmp_alleeg = cellfun(@(x) [[]; x], tmp_alleeg);
    %## CREATE NEW STUDY STRUCTURED
    [tmp_study, tmp_alleeg] = std_editset(MAIN_STUDY,tmp_alleeg,...
                                'rmclust','off',...
                                'addchannellabels','on',...
                                'name',study_fName,...
                                'commands',{'remove',find(rmv_subj)});
    %- study modifications
    tmp_study.urcluster = MAIN_STUDY.cluster;
    tmp_study.rmvd_subj_inds = find(rmv_subj);
    tmp_study.filename = [study_fName '.study'];
    tmp_study.name = study_fName;
    tmp_study.condition = TRIAL_TYPES{cond_i};
    %## ROBUST SAVE
    [tmp_study,tmp_alleeg] = eeglab_save_study(tmp_study,tmp_alleeg,...
                                        study_fName,save_dir,...
                                        'STUDY_COND',TRIAL_TYPES{cond_i});
end
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