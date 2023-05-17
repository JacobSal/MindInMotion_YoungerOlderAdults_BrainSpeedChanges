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
if ~ispc
    % see. eeg_options.m 
    % see. eeg_optionsbackup.m for all eeglab options.
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
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
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
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
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
% SUBJ_MISSING_TRIAL_DATA = {};
SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)
%- Subject Picks
SUBJ_PICS = {SUBJ_YNG};
SUBJ_ITERS = {1:length(SUBJ_YNG)}; % CHANG,LIU(12/26/2023)
% SUBJ_ITERS = {[],1:length(SUBJ_HMA),1:length(SUBJ_NMA)}; % CHANG,LIU(02/15/2023)
GROUP_NAMES = {'H1000s'};
%- Subject Directory Information
% PREPROCESS_NAME = 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10'; % CHANG,LIU(02/15/2023)
% OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-HY_202212';  % CHANG,LIU(02/15/2023)
% PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % CHANG,LIU(12/26/2022)
% OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % CHANG,LIU(12/26/2022)
OA_PREP_FPATH = '04182023_YA_N37_prep_verified'; % JACOB,SAL(04/10/2023)
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping % CHANG,LIU(12/26/2022)
%% ===================================================================== %%
%## PROCESSING PARAMS
% pop_editoptions('option_parallel',1);
params = [];
%- study group and saving
params.OVERRIDE_DIPFIT =  true;
% SAVE_EEG = false; % saves EEG structures throughout processing
%- component rejection crit
params.THRESHOLD_DIPFIT_RV = 0.15;
params.THRESHOLD_BRAIN_ICLABEL = 0.50;
%- MIM specific epoching
params.EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
%- epoching params
% params.TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
params.TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% params.TRIAL_TYPES = {'0p25','0p5','0p75','1p0'};
% params.TRIAL_TYPES = {'flat','low','med','high'};
params.EPOCH_TIME_LIMITS = [-0.5,5]; % [-1,3] captures gait events well
params.TRIAL_LENGTH = 3*60; % trial length in seconds
params.PER_OVERLAP = 0.0; % percent overlap between epochs
params.TRIAL_BEGIN_STR = 'TrialStart';
params.TRIAL_END_STR = 'TrialEnd';
params.EVENT_TRIAL_PARSER = 'type';
params.EVENT_COND_PARSER = 'cond';
%- connectivity process params
params.CONN_FREQS = (1:100);
params.CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
params.CNCTANL_TOOLBOX = 'sift'; %'bsmart'
params.WINDOW_LENGTH = 0.5;
params.WINDOW_STEP_SIZE = 0.025;
params.NEW_SAMPLE_RATE = [];
params.DO_BOOTSTRAP = true;
% ASSIGN_BOOTSTRAP_MEAN = false;
% SAVE_CONN_BOOTSTRAP = false;
%- subj_i_genConnStats
params.DO_PHASE_RND = true;
%- eeglab_cluster.m spectral params
params.FREQ_LIMITS = [1,100];
params.CYCLE_LIMITS = [3,0.8];
params.SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
params.FREQ_FAC = 4;
params.PAD_RATIO = 2;
%% global script chain (VERSION 1)
%- datetime override
% dt = '26122023'; % YA 
% dt = '24022023_YA'; % YA (MORE SUBJECTS)
% dt = '23032023_YA_fem'; % YA (03/23/2023)
dt = '04172023_MIM_YA_fullset_N_speed_terrain_merge';
% dt = '02022023';
%## PATH & TEMP STUDY NAME
%- hard define
DO_CONN_ANL = false;
SAVE_EEG = true;
SESSION_NUMBER = '1';
study_fName_1 = sprintf('%s_all_comps_study',[params.TRIAL_TYPES{:}]);
study_fName_2 = sprintf('%s_reduced_comps_study',[params.TRIAL_TYPES{:}]);
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
chanlocs_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
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
        fNames{cnt} = tmp.name;
        %- Chanlocs fPaths
        chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(params.TRIAL_TYPES,'_'); 
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
                        conditions,groups,sessions,...
                        'SAVE_EEG',SAVE_EEG,...
                        'CHANLOCS_FPATHS',chanlocs_fPaths);
    %- NOTES: spec_mode, 'psd' does not work.
    [MAIN_ALLEEG,MAIN_STUDY] = eeglab_cluster(ALLEEG,...
                        study_fName_1,save_dir,...
                        'SPEC_MODE',params.SPEC_MODE,...
                        'FREQ_LIMITS',params.FREQ_LIMITS,...
                        'CYCLE_LIMITS',params.CYCLE_LIMITS,...
                        'FREQ_FAC',params.FREQ_FAC,...
                        'PAD_RATIO',params.PAD_RATIO);
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
    %## PCA reduction algorithm
%     [tmps,tmpa,~,outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG); %## DEBUG
    [MAIN_STUDY,MAIN_ALLEEG,~,~] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
    %## ICA reduction algorithm
%     [tmpS,comps_out,comps_rej] = cluster_ica_reduce(MAIN_STUDY); %## DEBUG
%     [MAIN_STUDY,~,~] = cluster_ica_reduce(MAIN_STUDY);
    %## ADD ANATOMICAL LABELS
    [~, atlas_cell] = add_anatomical_labels(MAIN_STUDY,MAIN_ALLEEG);
    MAIN_STUDY.etc.add_anatomical_labels = atlas_cell;
    %## SAVE STUDY
%     [MAIN_STUDY,MAIN_ALLEEG] = eeglab_save_study(MAIN_STUDY,MAIN_ALLEEG,...
%                                             'reduced_comps_pca_study',save_dir);
    [MAIN_STUDY,MAIN_ALLEEG] = eeglab_save_study(MAIN_STUDY,MAIN_ALLEEG,...
                                            study_fName_2,save_dir);
    %## Extract components for each cluster & subject
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(MAIN_STUDY);
else
    fprintf(1,'\n==== LOADING REDUCED COMPONENTS STUDY DATA ====\n');
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_2 '_UNIX.study'],'filepath',save_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_2 '.study'],'filepath',save_dir);
    end
    fprintf(1,'\n==== DONE: LOADING REDUCED COMPONENTS STUDY DATA ====\n');
    %## Extract components for each cluster & subject
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(MAIN_STUDY);
end

%% PLOT DIPOLES, PSD'S, && TOPOPLOTS
plot_fNames = {'allDipPlot_top','allDipPlot_sagittal','allDipPlot_coronal','allSpecPlot','allTopoPlot'};
plot_chk = cellfun(@(x) ~exist([save_dir filesep sprintf('%s.jpg',x)],'file'),plot_fNames);
if any(plot_chk) && false
    fprintf('==== Making Dipole Plots ====\n');
%     [~] = std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),...
%             'mode','multicolor','figure','on');
    %- main plot
    [~] = std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',main_cl_inds(2:end),...
                'mode','multicolor','figure','on');
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_top.jpg')]);
    view([45,0,0])
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_sagittal.jpg')]);
    view([0,-45,0])
    saveas(fig_i,[save_dir filesep sprintf('allDipPlot_coronal.jpg')]);
    %- outlier plot
    [~] = std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',outlier_cl_inds,...
                'mode','multicolor','figure','on');
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_top.jpg')]);
    view([45,0,0])
    saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_sagittal.jpg')]);
    view([0,-45,0])
    saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_coronal.jpg')]);
    %{
    for cluster_i = 2:length(MAIN_STUDY.cluster)
        std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',cluster_i);
        fig_i = get(groot,'CurrentFigure');
    %     saveas(fig_i,[save_dir filesep sprintf('DipPlot_0.fig')]);
        saveas(fig_i,[save_dir filesep sprintf('cl%i_DipPlot_top.jpg',cluster_i)]);
        view([45,0,0])
    %     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.fig')]);
        saveas(fig_i,[save_dir filesep sprintf('cl%i_allDipPlot_sagittal.jpg',cluster_i)]);
        view([0,-45,0])
    %     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.fig')]);
        saveas(fig_i,[save_dir filesep sprintf('cl%i_allDipPlot_coronal.jpg',cluster_i)]);
    end
    %}
    %- Spec plot
    fprintf('==== Making Spectogram Plots ====\n');
    [MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                        'components',...                               
                                        'recompute','on',...
                                        'spec','on',...
                                        'specparams',...
                                        {'specmode',params.SPEC_MODE,'freqfac',params.FREQ_FAC,...
                                        'freqrange',params.FREQ_LIMITS});
    specMin = 10;
    specMax = 45;
    std_specplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),...
                    'ylim',[specMin,specMax],'freqrange',[1,65]);
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.jpg')]);
    %- Topo plot
    fprintf('==== Making Topograph Plots ====\n');
    [MAIN_STUDY, MAIN_ALLEEG] = std_precomp(MAIN_STUDY, MAIN_ALLEEG,...
                                        'components',...                               
                                        'recompute','on',...
                                        'scalp','on');
    std_topoplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster));
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.jpg')]);
    close all
    %## POP VIEW PROPS
    if ~exist([save_dir filesep 'component_props'],'dir') 
        mkdir([save_dir filesep 'component_props']);
        %## RECALCULATE ICAACT MATRICES
        MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
        for subj_i = 1:length(MAIN_ALLEEG)
            if isempty(MAIN_ALLEEG(subj_i).icaact)
                fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
                MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
                MAIN_ALLEEG(subj_i).icaact    = reshape( MAIN_ALLEEG(subj_i).icaact, size(MAIN_ALLEEG(subj_i).icaact,1), MAIN_ALLEEG(subj_i).pnts, MAIN_ALLEEG(subj_i).trials);
            end
        end
        for cluster_i = 2:length(MAIN_STUDY.cluster)
            sets_clust = MAIN_STUDY.cluster(cluster_i).sets;
            for i = 1:length(sets_clust)
                subj_i = sets_clust(i);
                comps_clust = MAIN_STUDY.cluster(cluster_i).comps(i);
                hold on;
                pop_prop_extended(MAIN_ALLEEG(subj_i),0,comps_clust,NaN,...
                {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
                fig = get(groot,'CurrentFigure');
                fig.Name = sprintf('%s_cluster%i_ic%i',...
                                    cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust);
                fig.Position = [500 300 1280 720]; 
                hold off;
                
%                 saveas(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.fig',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust)]);
                saveas(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.png',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust)]);
%                 savefig(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.fig',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust)]);
                close all
            end
        end
    end
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