%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_proc/run_conn_process.sh

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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'gamma' filesep 'MIM_proc'];
%- cd to source directory
cd(source_dir)
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
TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
%- Subject Picks
SUBJ_PICS = {SUBJ_YNG};
SUBJ_ITERS = {1:length(SUBJ_YNG)}; % CHANG,LIU(12/26/2023)
% SUBJ_ITERS = {[],1:length(SUBJ_HMA),1:length(SUBJ_NMA)}; % CHANG,LIU(02/15/2023)
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
% dt = '21012023';
% dt = '02022023';
%## PATH & TEMP STUDY NAME
%- hard define
SESSION_NUMBER = '1';
study_fName = sprintf('copy_study');
%- soft define
path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
mniMRI = fullfile(path2BEM, 'standard_mri.mat');
mniVol = fullfile(path2BEM, 'standard_vol.mat');
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
dip_hdm_fPath   = cell(1,length([SUBJ_ITERS{:}]));
dip_hdm_fName   = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fPath   = cell(1,length([SUBJ_ITERS{:}]));
dip_chn_fName   = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(fdSUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep PREPROCESS_NAME filesep SUBJ_PICS{group_i}{subj_i}];
        fNames{cnt} = [SUBJ_PICS{group_i}{subj_i} '_cleanEEG_' PREPROCESS_NAME '_' OUTSIDE_DATA_SUFFIX '.set'];
        dip_hdm_fPath{cnt} = [fPaths{cnt} filesep sprintf('%s_comp',SUBJ_PICS{group_i}{subj_i})];
        dip_hdm_fName{cnt} = 'dipfit_fem.mat';
        dip_chn_fPath{cnt} = fPaths{cnt};
        dip_chn_fName{cnt} = 'dipfit_templateElec.mat';
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file')>1)
        fprintf('Headmodel Exists: %i\n',exist([dip_hdm_fPath{cnt} filesep dip_hdm_fName{cnt}],'file')>1)
        fprintf('Channel File Exists: %i\n',exist([dip_chn_fPath{cnt} filesep dip_chn_fName{cnt}],'file')>1)
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = num2str(group_i);
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
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
    fprintf('==== Reassigning MRI for MNI plotting ====\n');
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
    for subj_i = 1:length(MAIN_ALLEEG)
        if isempty(MAIN_ALLEEG(subj_i).dipfit.mrifile)
            MAIN_ALLEEG(subj_i).dipfit.mrifile = mniMRI;
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.mrifile)
            MAIN_ALLEEG(subj_i).dipfit.hdmfile = mniVol; %[dip_hdm_fName{subj_i} filesep dip_hdm_fName{subj_i}];
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.coordformat) 
            MAIN_ALLEEG(subj_i).dipfit.coordformat = 'MNI';
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.chanfile)
            MAIN_ALLEEG(subj_i).dipfit.chanfile = [dip_chn_fPath{subj_i} filesep dip_chn_fName{subj_i}];
        end
    end
    %## ROBUST SAVE STUDY
    [MAIN_STUDY,MAIN_ALLEEG] = eeglab_save_study(MAIN_STUDY,MAIN_ALLEEG);
    fprintf(1,'\n==== DONE: CLUSTERING SUBJECT DATA ====\n');
else
    fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
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
    %- assiging important dipfit model information for later recall
    fprintf('==== Reassigning MRI for MNI plotting ====\n');
    % pop_clustedit(STUDY, ALLEEG, clusts); %## DEBUG
    for subj_i = 1:length(MAIN_ALLEEG)
        if isempty(MAIN_ALLEEG(subj_i).dipfit.mrifile)
            MAIN_ALLEEG(subj_i).dipfit.mrifile = mniMRI;
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.mrifile)
            MAIN_ALLEEG(subj_i).dipfit.hdmfile = mniVol; %[dip_hdm_fName{subj_i} filesep dip_hdm_fName{subj_i}];
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.coordformat) 
            MAIN_ALLEEG(subj_i).dipfit.coordformat = 'MNI';
        end
        if isempty(MAIN_ALLEEG(subj_i).dipfit.chanfile)
            MAIN_ALLEEG(subj_i).dipfit.chanfile = [dip_chn_fPath{subj_i} filesep dip_chn_fName{subj_i}];
        end
    end
    fprintf(1,'\n==== DONE: LOADING CLUSTER STUDY DATA ====\n');
end
%% POST-STUDY PROCESSING
%{
save_ALLEEG = MAIN_ALLEEG;
save_STUDY = MAIN_STUDY;
MAIN_ALLEEG = save_ALLEEG;
MAIN_STUDY = save_STUDY;
%- reset dipfit values using Chang's FEM fit
for i = 1:length(MAIN_ALLEEG)
    MAIN_ALLEEG(i).dipfit.model = MAIN_ALLEEG(i).dipfit_fem.model;
end
%}
%## PCA reduction algorithm
% [MAIN_STUDY,MAIN_ALLEEG,comps_out,outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
%## ICA reduction algorithm
[MAIN_STUDY, comps_out] = cluster_ica_reduce(MAIN_STUDY);
%## ADD ANATOMICAL LABELS
[~, atlas_cell] = add_anatomical_labels(MAIN_STUDY,MAIN_ALLEEG);
MAIN_STUDY.etc.add_anatomical_labels = atlas_cell;
%## PLOTS
std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),...
            'mode','multicolor');
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
                                    'recompute','on',...
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
    %- check to see if ICA is available from main struct
    if isempty(MAIN_ALLEEG(subj_i).icaact)
        error('Subject %s needs their EEG.icact calculated',MAIN_ALLEEG(subj_i).subject);
    end
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %## RUN MAIN_FUNC
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
    tmp_alleeg = cellfun(@(x) [val; x], tmp_alleeg);
    %## ROBUST SAVE
    [STUDY,tmp_alleeg] = eeglab_save_study(STUDY,tmp_alleeg);
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
%% RECOOPING STUDY FILES AFTER ERRROR
%{
dt = '12022023';
study_fName = sprintf('copy_study');
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
PARSE_TYPE = 'Constant';
SUFFIX_PATH_EPOCHED = 'Epoched';

%## 
fprintf(1,'\n==== LOADING CLUSTER STUDY DATA ====\n');
if DO_UNIX
    [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',save_dir);
else
    [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',save_dir);
end
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
end
[MAIN_STUDY,MAIN_ALLEEG,comps_out,outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);

%##
for cond_i = 1:length(TRIAL_TYPES)
    tmp_alleeg = cell(1,length(MAIN_ALLEEG));
    study_name = sprintf('%s_MIM_study',TRIAL_TYPES{cond_i});
    %- study modifications
    tmp_study = MAIN_STUDY;
    tmp_study.name = study_name;
    tmp_study.condition = TRIAL_TYPES{cond_i};
    tmp_study.filename = [study_name '.study'];
    %- extract each condition for each subject
    for subj_i = 1:length(tmp_alleeg)
        epoch_save_dir = [save_dir filesep MAIN_ALLEEG(subj_i).subject filesep SUFFIX_PATH_EPOCHED filesep PARSE_TYPE filesep TRIAL_TYPES{cond_i}];
        EEG = pop_loadset('filename',sprintf('%s_%s_EPOCH_TMPEEG.set',MAIN_ALLEEG(subj_i).subject,TRIAL_TYPES{cond_i}),'filepath',epoch_save_dir);
        tmp_alleeg{subj_i} = EEG;
        if DO_UNIX
            tmp_study.datasetinfo(subj_i).condition = TRIAL_TYPES{cond_i};
            tmp_study.datasetinfo(subj_i).filepath = convertPath2UNIX(EEG.filepath);
            tmp_study.datasetinfo(subj_i).filename = EEG.filename;
        else
            tmp_study.datasetinfo(subj_i).condition = TRIAL_TYPES{cond_i};
            tmp_study.datasetinfo(subj_i).filepath = convertPath2Drive(EEG.filepath);
            tmp_study.datasetinfo(subj_i).filename = EEG.filename;
        end
    end
    val = [];
    tmp_alleeg = cellfun(@(x) [val; x], tmp_alleeg);
    disp(tmp_alleeg)
    %## ROBUST SAVE
    %- assign new filepaths for STUDY file 
    if DO_UNIX
        %- save STUDY for UNIX
        [~] = pop_savestudy( tmp_study, tmp_alleeg,...
                            'filename',[study_name '_UNIX'],'filepath',save_dir);
        %- Change datasetinfo filepaths and resave for accessing via DRIVE
        for subj_i = 1:length(tmp_study.datasetinfo)
            tmp_study.datasetinfo(subj_i).condition = TRIAL_TYPES{cond_i};
            tmp_study.datasetinfo(subj_i).filepath = convertPath2Drive(tmp_alleeg(subj_i).filepath);
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
            tmp_study.datasetinfo(subj_i).condition = TRIAL_TYPES{cond_i};
            tmp_study.datasetinfo(subj_i).filepath = convertPath2UNIX(tmp_alleeg(subj_i).filepath);
            tmp_study.datasetinfo(subj_i).filename = tmp_alleeg(subj_i).filename;
        end
        [~] = pop_savestudy( tmp_study, tmp_alleeg,...
                            'filename',[study_name '_UNIX'],'filepath',save_dir);
    end
end
%}