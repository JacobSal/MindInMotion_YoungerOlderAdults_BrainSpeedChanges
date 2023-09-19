% Create study for cluster ICs. This code only works for cluster without
% using ERSP. Precompute ERSP needed to be done on Hipergator
% Chang Liu - 2021-11-23 - V1

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%   NJacobsen notes
%   When timewarping data, save values as EEG.timewarp = timewarp;
%   EEG.timewarp.medianlatency = median(timewarp.latencies(:,:));%Warping to the median latency of my 5 events
%   By default, std_ersp will use the median of all subject's
%   timewarp.latencies(:,:) as 'timewarpms' unless individual subject 
%   warpto is indiciated using 'timewarpms', 'subject tw matrix'
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230417.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_b_ic_clustering_refine.sh

%{
%## RESTORE MATLABs
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
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
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_OA'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP ======================================================= %%
if ~ispc
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
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result? 
SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
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
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
% NOTE: see. statcondfieldtrip.m or std_stat.m
COND_EVENT_CHAR = 'cond';
%- clustering parameters
MIN_ICS_SUBJ = [2,3,4,5,6,7,8]; % iterative clustering
DISTS_FROM_MEDIAN = [10,20,30,40,50,60,70]; % distance in units (mm)
K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
cluster_ks = K_RANGE(1):K_RANGE(2);
% (08/21/2023) JS, this currenlty doesn't do anything but take up more
% memory.
REPEATED_CLUSTERING_STD = 3;
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',[],...
    'filepath',[],...
    'outliers',inf(),...
    'maxiter',200);
%- 
%- custom params
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
STD_PRECLUST_COMMAND = {'dipoles','weight',clustering_weights.dipoles};
%- iterative clustering parameters
n_iterations = 50;
outlier_sigma = 3;
%- datetime override
% dt = '05182023_MIM_OA_subset_N85_oldpipe';
% dt = '05192023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06122023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06282023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07152023_MIM_OAN79_subset_prep_verified_gait';
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## ADMIN SET

SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\dipole_1_scalp_0_ersp_0_spec_0';
%- convert SUB_DIR
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true,true,true,true};
CLUSTER_K_PICKS = [14,14,14,14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics2','temp_study_rejics3','temp_study_rejics7'};
CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics2' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics3' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics7' filesep '14']};
CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat'};
CLUSTER_STUDY_DIRS = {[CLUSTER_DIRS{1} filesep 'spec_data'],...
    [SUB_DIR filesep 'subjrejs_minics2'],...
    [SUB_DIR filesep 'subjrejs_minics3'],...
    [SUB_DIR filesep 'subjrejs_minics7']};

POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
%% ===================================================================== %%
for k_i = 3:length(CLUSTER_DIRS)
    %## CREATE & GRAB DIRECTORIES
    %- convert cluster directory
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    %- convert study directory
    if ~ispc
        cluster_study_dir = convertPath2UNIX(CLUSTER_STUDY_DIRS{k_i});
    else
        cluster_study_dir = convertPath2Drive(CLUSTER_STUDY_DIRS{k_i});
    end
    %## LOAD STUDY
    if LOAD_DIFFERENT_STUDY{k_i}
        %- Create STUDY & ALLEEG structs
        if ~exist([cluster_study_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
            error('error. study file does not exist');
        else
            if ~ispc
                tmp = load('-mat',[cluster_study_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
                TMP_STUDY = tmp.STUDY;
            else
                tmp = load('-mat',[cluster_study_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
                TMP_STUDY = tmp.STUDY;
            end
        end
    else
        error('error. assign a study file.')
    end 
    %- add cluster information
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    TMP_STUDY.cluster = cluster_update;
    clust_i = CLUSTER_K_PICKS(k_i);
    %## Create PowerPoint
    %- start powerpoint
    ppt = actxserver('powerpoint.application');
%         ppt.Visible = 1;% Open PowerPoint and make it visible
    ppt.Presentation.invoke;
    ppt.Presentation.Add();
    layout = ppt.ActivePresentation.SlideMaster.CustomLayouts.Item(1);
    ppt.ActivePresentation.Slides.AddSlide(1, layout);
    layout = ppt.ActiveWindow.Selection.SlideRange(1).CustomLayout;
    slides = ppt.ActivePresentation.Slides;
    %##
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cluster,~,nonzero_cluster] = eeglab_get_cluster_comps(TMP_STUDY);
    CLUSTER_PICKS = nonzero_cluster; %valid_cluster; %main_cl_inds(2:end);
    fprintf('Clusters with more than 50%% of subjects:'); fprintf('%i,',valid_cluster(1:end-1)); fprintf('%i',valid_cluster(end)); fprintf('\n');
    fprintf('Main cluster numbers:'); fprintf('%i,',main_cl_inds(1:end-1)); fprintf('%i',main_cl_inds(end)); fprintf('\n');
    %## CLUSTER SPECIFIC FIGURES
    cl_dir = dir([cluster_dir filesep '*.jpg']);
    inds = cellfun(@(x) contains(x,'topo_plot_averaged'),{cl_dir.name});
    topo_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    inds = cellfun(@(x) contains(x,'alldips_per_clust_top'),{cl_dir.name});
    diptop_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    inds = cellfun(@(x) contains(x,'alldips_per_clust_sagittal'),{cl_dir.name});
    dipsag_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    inds = cellfun(@(x) contains(x,'alldips_per_clust_coronal'),{cl_dir.name});
    dipcor_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    inds = cellfun(@(x) contains(x,'allSpecPlot_des1'),{cl_dir.name});
    psd_des1_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    inds = cellfun(@(x) contains(x,'allSpecPlot_des2'),{cl_dir.name});
    psd_des2_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
    %## PPT MAKE
    %## (SLIDES 4) Subject Specific ERSPs & Criteria
%             im_scale = 0.225;
    subj_inds = 1:length(TMP_STUDY.datasetinfo);
    for i = 1:length(subj_inds)
        subj_char = TMP_STUDY.datasetinfo(i).subject;
        %-
%                 subj_icscores = [load_dir filesep subj_char filesep 'ICA' filesep 'IC_score.jpg'];
        subj_topos = [load_dir filesep subj_char filesep 'ICA' filesep 'potential_brain_components_allcritera_customElectrode.jpg'];
        subj_keeppsds = [load_dir filesep subj_char filesep 'ICA' filesep 'KEEP_IC_PSD.jpg'];
        subj_powpow = [load_dir filesep subj_char filesep 'ICA' filesep 'powpowcat.jpg'];
        subj_rejpsds = [load_dir filesep subj_char filesep 'ICA' filesep 'reject_potential_brain_components_spectral.jpg'];
        subj_icstem = [load_dir filesep subj_char filesep 'ICA' filesep 'spectral_stem.jpg'];
        subj_icscore = [load_dir filesep subj_char filesep 'ICA' filesep 'IC_score.jpg'];
        %- Subject Specific IC Criteria
        TOP_DIST = 30;
        newSlide = slides.AddSlide(1,layout);
        Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
        Title1.TextFrame.TextRange.Text = sprintf('Subject %s',subj_char);
        im_scale = 0.35;
        if exist(subj_topos,'file')
            newSlide.Shapes.AddPicture(subj_rejpsds, 'msoFalse', 'msoTrue',0+100,0+TOP_DIST,2250*im_scale,1500*im_scale); %Left, top, width, height
        end
%         im_scale = 0.3;
%         if exist(subj_keeppsds,'file')
%             newSlide.Shapes.AddPicture(subj_icstem, 'msoFalse', 'msoTrue',875*im_scale+70,0+TOP_DIST,2250*im_scale,1500*im_scale); %Left, top, width, height
%         end
        %- Subject Specific IC Criteria
        TOP_DIST = 30;
        newSlide = slides.AddSlide(1,layout);
        Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
        Title1.TextFrame.TextRange.Text = sprintf('Subject %s',subj_char);
        im_scale = 0.35;
        if exist(subj_powpow,'file')
            newSlide.Shapes.AddPicture(subj_icscore, 'msoFalse', 'msoTrue',0+100,0+TOP_DIST,2250*im_scale,1500*im_scale); %Left, top, width, height
        end
        %- Subject Specific IC Criteria
        TOP_DIST = 70;
        im_scale = 0.275;
        newSlide = slides.AddSlide(1,layout);
        Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
        Title1.TextFrame.TextRange.Text = sprintf('Subject %s',subj_char);
        if exist(subj_topos,'file')
            newSlide.Shapes.AddPicture(subj_topos, 'msoFalse', 'msoTrue',0,0+TOP_DIST,1167*im_scale,875*im_scale); %Left, top, width, height
        end
        im_scale = 0.3;
        if exist(subj_keeppsds,'file')
            newSlide.Shapes.AddPicture(subj_keeppsds, 'msoFalse', 'msoTrue',875*im_scale+70,0+TOP_DIST,2250*im_scale,1500*im_scale); %Left, top, width, height
        end
        im_scale = 0.275;
        if exist(subj_powpow,'file')
            newSlide.Shapes.AddPicture(subj_powpow, 'msoFalse', 'msoTrue',0,0+TOP_DIST+875*im_scale,1167*im_scale,875*im_scale); %Left, top, width, height
        end
    end
    %## (SLIDE TITLE)
    dt = datetime;
    dt.Format = 'MMddyyyy';
    %- title slide
    newSlide = slides.AddSlide(1,layout);
    set(newSlide.Shapes.Title.Textframe.Textrange, 'Text', 'Subject Brain IC Criteria');
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    tmp ={TMP_STUDY.datasetinfo.subject};
    tmp = strjoin(tmp,', ');
    Title1.TextFrame.TextRange.Text = sprintf('Subjects: %s\n',tmp);
    ppt.ActivePresentation.SaveAs([cluster_dir filesep sprintf('%s_summary_subject_rejection_report.pptx',dt)]);
    % Close Powerpoint and delete the object
    ppt.ActivePresentation.Close;
    ppt.Quit;
    ppt.delete;
    fprintf('\nResults powerpoint saved here:\n%s\n',cluster_dir);
end