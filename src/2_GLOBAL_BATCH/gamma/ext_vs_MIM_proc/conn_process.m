%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/ext_vs_MIM_proc/run_conn_process.sh

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
dt.Format = 'ddMMyyyy';
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
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
end
%## DEBUG: PATHROOT OVERRIDE
if DO_UNIX
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'gamma' filesep 'ext_vs_MIM_proc'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
setWorkspace
%% PARPOOL SETUP
if DO_UNIX
%     eeg_options;
    pop_editoptions('option_parallel',1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
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
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', 1440);
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'ext_vs_MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
%- MIND IN MOTION
SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
    'H1022','H1024','H1026','H1027','H1033','H1034'};
SUBJ_HMA = {'H2002', 'H2010', 'H2015', 'H2017', 'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2034', 'H2059', 'H2062', 'H2082', 'H2095'};
SUBJ_NMA = {'NH3008', 'NH3043', 'NH3055', 'NH3059', 'NH3069', ...
    'NH3070', 'NH3074', 'NH3086', 'NH3090', 'NH3104', 'NH3105', 'NH3106', 'NH3112', 'NH3114'};
DATASET_PROCESSNUM = 1;
TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
%- Subject Picks
SUBJ_PICS = {SUBJ_YNG,SUBJ_HMA,SUBJ_NMA};
% SUBJ_ITERS = {1:length(SUBJ_YNG),1:length(SUBJ_HMA),1:length(SUBJ_NMA)};
% SUBJ_ITERS = {[1,2],[],[]};
SUBJ_ITERS = {1:length(SUBJ_YNG),[],[]};
%- Subject Directory Information
PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % value to help with file looping
OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % value to help with file looping
if DO_UNIX
    OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR,'dferris');
else
    OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR,'M');
end
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping
%% ===================================================================== %%
%## PROCESSING PARAMS
% pop_editoptions('option_parallel',1);
%- study group and saving
GROUP_INT = 1;
SAVE_EEG = false; % saves EEG structures throughout processing
%- component rejection crit
THRESHOLD_DIPFIT_RV = 0.15;
THRESHOLD_BRAIN_ICLABEL = 0.50;
%- MIM specific epoching
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
%- subj_i_epoch
PARSE_TYPE = 'Constant'; %
EPOCH_TIME_LIMITS = [-1,1];
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
CNCTANL_TOOLBOX = 'bsmart'; %'sift'; %'bsmart'
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
DO_BOOTSTRAP = false;
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
% dt = '16012023';
%## PATH & TEMP STUDY NAME
%- hard define
study_fName = sprintf('copy_study');
%- soft define
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
conditions = cell(1,length([SUBJ_ITERS{:}]));
groups = cell(1,length([SUBJ_ITERS{:}]));
sessions = cell(1,length([SUBJ_ITERS{:}]));
subjectNames = cell(1,length([SUBJ_ITERS{:}]));
fNames = cell(1,length([SUBJ_ITERS{:}]));
fPaths = cell(1,length([SUBJ_ITERS{:}]));
dip_hdm_fPath = cell(1,length([SUBJ_ITERS{:}]));
dip_hdm_fName = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fPath = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fName = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        cnt = subj_i + stack_iter;
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep PREPROCESS_NAME filesep SUBJ_PICS{group_i}{subj_i}];
        fNames{cnt} = [SUBJ_PICS{group_i}{subj_i} '_cleanEEG_' PREPROCESS_NAME '_' OUTSIDE_DATA_SUFFIX '.set'];
        dip_hdm_fPath{cnt} = [fPaths{subj_i} filesep sprintf('%s_comp',SUBJ_PICS{group_i}{subj_i})];
        dip_hdm_fName{cnt} = 'dipfit_fem.mat';
        dip_chn_fPath{cnt} = fPaths{subj_i};
        dip_chn_fName{cnt} = 'dipfit_templateElec.mat';
        fprintf('Exists: %i\n',exist([fPaths{subj_i} filesep fNames{subj_i}],'file'))
    end
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        cnt = subj_i + stack_iter;
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = num2str(group_i);
        sessions{cnt} = '1';
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%% CREATE STUDY
%- NOTE: need a dipole fit before this step can be done.
%- params
%## Create STUDY & ALLEEG structs
if ~exist([save_dir filesep study_fName '.study'],'file') %|| true
    fprintf(1,'\n==== CLUSTERING SUBJECT DATA ====\n');
    %- NOTES: spec_mode, 'psd' does not work.
    [MAIN_ALLEEG,MAIN_STUDY] = feval(@eeglab_cluster,fNames,fPaths,...
                        subjectNames,study_fName,save_dir,...
                        conditions,groups,sessions,...
                        'SPEC_MODE',SPEC_MODE,...
                        'FREQ_LIMITS',FREQ_LIMITS,...
                        'CYCLE_LIMITS',CYCLE_LIMITS,...
                        'FREQ_FAC', FREQ_FAC,...
                        'PAD_RATIO', PAD_RATIO);
    %- Recalculate ICA Matrices && Book Keeping
    fprintf('==== Calculating ICA Matrices ====\n');
    MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
    for subj_i = 1:length(MAIN_ALLEEG)
        if isempty(MAIN_ALLEEG(subj_i).icaact)
            fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
            MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
        end
    end
    %- assiging important dipfit model information for later recall
    fprintf('==== Reassigning MRI for MNI plotting ====\n');
    % pop_clustedit(STUDY, ALLEEG, clusts); %## DEBUG
    path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
    mniMRI = fullfile(path2BEM, 'standard_mri.mat');
    mniVol = fullfile(path2BEM, 'standard_vol.mat');
    for subj_i = 1:length(MAIN_ALLEEG)
         MAIN_ALLEEG(subj_i).dipfit.mrifile = mniMRI;
         MAIN_ALLEEG(subj_i).dipfit.hdmfile = mniVol; %[dip_hdm_fName{subj_i} filesep dip_hdm_fName{subj_i}];
         MAIN_ALLEEG(subj_i).dipfit.coordformat = 'MNI';
         MAIN_ALLEEG(subj_i).dipfit.chanfile = [dip_chn_fPath{subj_i} filesep dip_chn_fName{subj_i}];
    end
    
    %## ROBUST SAVE STUDY
    if DO_UNIX
        %- Save STUDY for UNIX
        [MAIN_STUDY,MAIN_ALLEEG] = pop_savestudy( MAIN_STUDY, MAIN_ALLEEG,...
                            'filename',[study_fName '_UNIX'],'filepath',save_dir);
        %- Change datasetinfo filepaths and resave for accessing via DRIVE
        for subj_i = 1:length(MAIN_STUDY.datasetinfo)
            MAIN_STUDY.datasetinfo(subj_i).filepath = convertPath2Drive(MAIN_ALLEEG(subj_i).filepath,'M');
        end
        [MAIN_STUDY,MAIN_ALLEEG] = pop_savestudy( MAIN_STUDY, MAIN_ALLEEG,...
                            'filename',study_fName,'filepath',save_dir);
    else
        %- save STUDY for DRIVE
        for subj_i = 1:length(MAIN_STUDY.datasetinfo)
            MAIN_STUDY.datasetinfo(subj_i).filepath = convertPath2UNIX(MAIN_ALLEEG(subj_i).filepath,'dferris');
        end
        %- Change datasetinfo filepaths and resave for accessing via UNIX
        [MAIN_STUDY,MAIN_ALLEEG] = pop_savestudy( MAIN_STUDY, MAIN_ALLEEG,...
                            'filename',[study_fName '_UNIX'],'filepath',save_dir);
    end
    fprintf(1,'\n==== DONE: CLUSTERING SUBJECT DATA ====\n');
else
    
    %## Load STUDY
    if DO_UNIX
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',save_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',save_dir);
    end
    %## Recalculate ICA Matrices && Book Keeping
    fprintf('==== Recalculating ICA Matrices && checking .set file ====\n');
    MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
    for subj_i = 1:length(MAIN_ALLEEG)
        if isempty(MAIN_ALLEEG(subj_i).icaact)
            fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
            MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
        end
    end
end
%## PCA reduction algorithm
[MAIN_STUDY, MAIN_ALLEEG, comps_out, outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
%## ICA reduction algorithm
% [MAIN_STUDY, comps_out] = cluster_ica_reduce(MAIN_STUDY);
%## ADD ANATOMICAL LABELS
[~, atlas_cell] = add_anatomical_labels(MAIN_STUDY,MAIN_ALLEEG);
MAIN_STUDY.etc.add_anatomical_labels = atlas_cell;
%## PLOTS
std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),'mode','multicolor');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_0.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_0.jpg')]);
view([45,0,0])
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.jpg')]);
view([0,-45,0])
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.jpg')]);
%- Spec plot
[MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                    'components',...                               
                                    'recompute','on',....
                                    'spec','on',...
                                    'specparams',...
                                    {'specmode',SPEC_MODE,'freqfac',4,...
                                    'freqrange',FREQ_LIMITS});
specMin = 10;
specMax = 45;
std_specplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.jpg')]);
%- Topo plot
[MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                    'components',...                               
                                    'recompute','on',...
                                    'scalp','on');
std_topoplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster));
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.jpg')]);
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
        end
    end
end
%% GENERATE CONNECTIVITY METRICS (PARFOR)
ALLEEG = [];
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
parfor (subj_i = LOOP_VAR,POOL_SIZE)
%         for subj_i = 1
    %- check to see if ICA is available from main struct
    if isempty(MAIN_ALLEEG(subj_i).icaact)
        MAIN_ALLEEG(subj_i) = eeg_checkset(MAIN_ALLEEG(subj_i),'loaddata');
        MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
    end
    
    %## RUN MAIN_FUNC
    %(11/11/2022), JS: really need to consider updating bootstrap
    %algorithm with parallel computing. Taking ~ 1 day per
    %condition for all subjects and the bottle neck is entirely the
    %bootstrap.
    %
    %Note: validateattributes and assert functions may be helpful
    %in more clearly defining function inputs.
    % e.g.  DO_PHASE_RND = true;
    %       errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
    %       validationFcn = @(x) assert(islogical(x),errorMsg);
    %       

    % (12/5/2022) Need to adapt this to include all conditions
    % within each SUBJ structure so connectivity can be calculated
    % for the ALLEEG structure rather than the EEG structure.
    % *** maybe try to ditch the SUBJ strucutre entirely for this
    % round?
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %-
    [ALLEEG,pathsout] = main_func(MAIN_ALLEEG(subj_i),TRIAL_TYPES,save_dir,components,...
        'SAVE_EEG',SAVE_EEG,...
        'PATH_EXT',PATH_EXT,...
        'PARSE_TYPE',PARSE_TYPE,... %
        'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,...
        'TRIAL_LENGTH',TRIAL_LENGTH,...
        'PER_OVERLAP',PER_OVERLAP,...
        'EVENT_CHAR',EVENT_CHAR,...
        'TRIAL_BEGIN_STR',TRIAL_BEGIN_STR,...
        'TRIAL_END_STR',TRIAL_END_STR,...
        'EVENT_TRIAL_PARSER',EVENT_TRIAL_PARSER,...
        'EVENT_COND_PARSER',EVENT_COND_PARSER,...
        'CONN_METHODS',CONN_METHODS,... %
        'FREQS',CONN_FREQS,...
        'CNCTANL_TOOLBOX',CNCTANL_TOOLBOX,... 
        'DO_BOOTSTRAP',DO_BOOTSTRAP,...
        'WINDOW_LENGTH',WINDOW_LENGTH,...
        'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,...
        'DO_PHASE_RND',DO_PHASE_RND);
    disp(pathsout)
    tmp{subj_i} = ALLEEG;
end
%% SAVE CONDITION STUDIES
for cond_i = 1:length(TRIAL_TYPES)
    tmp_alleeg = cell(1,length(tmp));
    study_name = sprintf('%s_MIM_study',TRIAL_TYPES{cond_i});
    %- study modifications
    tmp_study = MAIN_STUDY;
    tmp_study.name = study_name;
    tmp_study.condition = TRIAL_TYPES{cond_i};
    %- extract each condition for each subject
    for subj_i = 1:length(tmp)
        tmp_alleeg{subj_i} = tmp{subj_i}(cond_i);
    end
    val = [];
    tmp_alleeg = cellfun(@(x) [val; x], tmp_alleeg, 'UnfiormOutput', false);
    %## ROBUST SAVE
    %- assign new filepaths for STUDY file 
    if DO_UNIX
        %- save STUDY for UNIX
        [~] = pop_savestudy(tmp_studys, tmp_alleeg,...
                            'filename',[study_name '_UNIX'],'filepath',save_dir);
        %- Change datasetinfo filepaths and resave for accessing via DRIVE
        for subj_i = 1:length(tmp_study.datasetinfo)
            tmp_study.datasetinfo(subj_i).filepath = convertPath2Drive(tmp_alleeg(subj_i).filepath,'M');
            tmp_study.datasetinfo(subj_i).filename = tmp_alleeg(subj_i).filename;
        end
        [~] = pop_savestudy( tmp_study, tmp_alleeg,...
                            'filename',study_name,'filepath',save_dir);
    else
        %- save STUDY for DRIVE
        [~] = pop_savestudy(tmp_study, tmp_alleeg,...
                            'filename',study_name,'filepath',save_dir);
        %- Change datasetinfo filepaths and resave for accessing via UNIX
        for subj_i = 1:length(tmp_study.datasetinfo)
            tmp_study.datasetinfo(subj_i).filepath = convertPath2UNIX(tmp_alleeg(subj_i).filepath,'dferris');
            tmp_study.datasetinfo(subj_i).filename = tmp_alleeg(subj_i).filename;
        end
        [~] = pop_savestudy( tmp_study, tmp_alleeg,...
                            'filename',[study_name '_UNIX'],'filepath',save_dir);
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:
dt = '26122022';
%## PROCESSING PARAMS
% pop_editoptions('option_parallel',1);
%- subjstruct and saving
SAVE_EEG = false; % saves EEG structures throughout processing
%- component rejection crit
THRESHOLD_DIPFIT_RV = 0.15;
THRESHOLD_BRAIN_ICLABEL = 0.50;
%- subj_i_epoch
PARSE_TYPE = 'Constant'; %
EPOCH_TIME_LIMITS = [-2,2];
TRIAL_LENGTH = 3*60; % trial length in seconds
PER_OVERLAP = 0.0; % percent overlap between epochs
TRIAL_BEGIN_STR = 'TrialStart';
TRIAL_END_STR = 'TrialEnd';
EVENT_TRIAL_PARSER = 'type';
EVENT_COND_PARSER = 'cond';
%- connectivity process
FREQS = (3:60);
CONN_METHODS = {'dDTF','GGC','ffDTF'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
% ANAT_THRESH = 30; % distance to anatomical location in (mm)
%- subj_i_genConnMeas
CNCTANL_TOOLBOX = 'sift'; %'bsmart'
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
DO_BOOTSTRAP = false;
ASSIGN_BOOTSTRAP_MEAN = false;
SAVE_CONN_BOOTSTRAP = false;
%- subj_i_genConnStats
DO_PHASE_RND = true;
%- eeglab_cluster.m spectral params
FREQ_LIMITS = [1,100];
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
PAD_RATIO = 2;
%}

%% TESTING ICA vs. PCA CLUSTERING
%{
[tmp_study,tmp_alleeg,comps_out] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
[tmp_study2] = cluster_ica_reduce(MAIN_STUDY);
%##COMPARISON PLOTS
std_dipplot(tmp_study2,MAIN_ALLEEG,'clusters',2:length(tmp_study2.cluster),'mode','multicolor');
std_dipplot(tmp_study,tmp_alleeg,'clusters',2:length(tmp_study.cluster),'mode','multicolor');
% view([45,-45,0])
% view([0,-45,0])
% view([0,0,45])
%- Spec plot ICA
specMin = 15;
specMax = 40;
std_specplot(tmp_study2,MAIN_ALLEEG,'clusters',1:length(tmp_study2.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Spec plot PCA
std_specplot(tmp_study,tmp_alleeg,'clusters',1:length(tmp_study.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Topo plot ICA
std_topoplot(tmp_study2,MAIN_ALLEEG,'clusters', 'all');
%- Topo plot PCA
std_topoplot(tmp_study,tmp_alleeg,'clusters', 'all');
%## POP VIEW PROPS (ICA)
if ~exist([save_dir filesep 'component_props'],'dir')
    mkdir([save_dir filesep 'component_props']);
end
for cluster_i = 2:length(tmp_study2.cluster)
    sets_clust = tmp_study2.cluster(cluster_i).sets;
    for i = 1:length(sets_clust)
        subj_i = sets_clust(i);
        comps_clust = tmp_study2.cluster(cluster_i).comps(i);
        hold on;
        pop_prop_extended(MAIN_ALLEEG(subj_i),0,comps_clust,NaN,...
        {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
        fig = gcf;
        fig.Position = [500 300 920 480]; 
        hold off;
%         savefig(fig,[save_dir filesep 'component_props' filesep sprintf('%i_%s_viewprops_co%i_cl%i.fig',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust,cluster_i)]);
    end
end
%## POP VIEW PROPS (PCA)
if ~exist([save_dir filesep 'component_props'],'dir')
    mkdir([save_dir filesep 'component_props']);
end
for cluster_i = 2:length(tmp_study.cluster)
    sets_clust = tmp_study.cluster(cluster_i).sets;
    for i = 1:length(sets_clust)
        subj_i = sets_clust(i);
        comps_clust = tmp_study.cluster(cluster_i).comps(i);
        hold on;
        pop_prop_extended(tmp_alleeg(subj_i),0,comps_clust,NaN,...
        {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
        fig = gcf;
        fig.Position = [500 300 920 480]; 
        hold off;
%         savefig(fig,[save_dir filesep 'component_props' filesep sprintf('%i_%s_viewprops_co%i_cl%i.fig',cluster_i,tmp_alleeg(subj_i).subject,comps_clust,cluster_i)]);
    end
end
%}