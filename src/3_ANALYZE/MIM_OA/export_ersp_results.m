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
% LOAD_DIFFERENT_STUDY = {true,true,true,true};
% CLUSTER_K_PICKS = [14,14,14,14];
% CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics2','temp_study_rejics3','temp_study_rejics7'};
% CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14'],...
%     [SUB_DIR filesep 'subjrejs_minics2' filesep '14'],...
%     [SUB_DIR filesep 'subjrejs_minics3' filesep '14'],...
%     [SUB_DIR filesep 'subjrejs_minics7' filesep '14']};
% CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat'};
% CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6'],...
%     [SUB_DIR filesep 'subjrejs_minics2'],...
%     [SUB_DIR filesep 'subjrejs_minics3'],...
%     [SUB_DIR filesep 'subjrejs_minics7']};
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6'};
CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14']};
CLUSTER_FILES = {'cluster_update_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6']};
%% ===================================================================== %%
for k_i = 1:length(CLUSTER_DIRS)
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
    %- assign plots_out directory
    plot_store_dir = [cluster_dir filesep 'plots_out'];
    if ~exist(plot_store_dir,'dir')
        error('error. directory does not exist: %s\n',plot_store_dir);
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
%     inds = cellfun(@(x) contains(x,'avgdips_per_clust_top'),{cluster_dir.name});
%     avg_diptop_files = cluster_dir(inds);
%     inds = cellfun(@(x) contains(x,'avgdips_per_clust_sagittal'),{cluster_dir.name});
%     avg_dipsag_files = cluster_dir(inds);
%     inds = cellfun(@(x) contains(x,'avgdips_per_clust_coronal'),{cluster_dir.name});
%     avg_dipcor_files = cluster_dir(inds);%- variables
    %## SUBJECT SPECIFIC IC REJECTION FIGURES
    %- assign subject specific IC rejection directory
%     sub_icrej_dir = dir([load_dir filesep '*' filesep 'ICA' filesep '*.jpg']);
%     %*
%     inds = cellfun(@(x) contains(x,'IC_score'),{sub_icrej_dir.name});
%     icscore_files = sub_icrej_dir(inds);
%     inds = cellfun(@(x) contains(x,'KEEP_IC_PSD'),{sub_icrej_dir.name});
%     ickeeppsd_files = sub_icrej_dir(inds);
%     inds = cellfun(@(x) contains(x,'powpowcat'),{sub_icrej_dir.name});
%     powpow_files = sub_icrej_dir(inds);
%     inds = cellfun(@(x) contains(x,'potential_brain_components_allcritera_customElectrode'),{sub_icrej_dir.name});
%     ickeeptopo_files = sub_icrej_dir(inds);
    %## PPT MAKE
    Image = [];
%         cnt = 1;
    cl_names = {TMP_STUDY.cluster.name};
    for k = main_cl_inds %length(cl_names):-1:1 %1:length(cl_names)
        ind1 = cellfun(@(x) strcmp(x,sprintf('Cluster_topo_plot_averaged_%i.jpg',k)),{topo_files.name});
        ind2 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldips_per_clust_top_%i.jpg',k)),{diptop_files.name});
        ind3 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldips_per_clust_sagittal_%i.jpg',k)),{dipsag_files.name});
        ind4 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldips_per_clust_coronal_%i.jpg',k)),{dipcor_files.name});
        ind5 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des1_cl%i.jpg',k)),{psd_des1_files.name});
        ind6 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des2_cl%i.jpg',k)),{psd_des2_files.name});
        %*
        % opts: (design 1, cluster 3) erspplot1_des1_cl3_stats_lim60,
        % erspplottype2_notfull_stats_compare-flat_des1_cl3,
        % erspplottype2_stats_compare-flat_des1_cl3,
        % erspplottype2_stats_des1_cl3,
        % erspplottype2_stats_notfull_des1_cl3,
        % erspplottype3_orig_des1_cl3_lim60, && subject specific e.g., H2008_within_des1_cl3_stats
        ersp_fpath = [cluster_dir filesep 'plots_out' filesep sprintf('%i',k)];
        %- checks
        chk1 = sum(ind1)>=1 && sum(ind2)>=1 && sum(ind3)>=1 && sum(ind4)>=1 && sum(ind5)>=1 && sum(ind6)>=1;
        chk2 = exist(ersp_fpath,'dir');
        if chk1 && chk2
            %-
            D1 = [diptop_files(ind2).folder filesep diptop_files(ind2).name];
            D2 = [dipsag_files(ind3).folder filesep dipsag_files(ind3).name];
            D3 = [dipcor_files(ind4).folder filesep dipcor_files(ind4).name];
            topo = [topo_files(ind1).folder filesep topo_files(ind1).name];
            psd1 = [psd_des1_files(ind5).folder filesep psd_des1_files(ind5).name];
            psd2 = [psd_des2_files(ind6).folder filesep psd_des2_files(ind6).name];
            %## (SLIDES 4) Subject Specific ERSPs & Criteria
            %{
%             im_scale = 0.225;
            subj_inds = unique(TMP_STUDY.cluster(k).sets);
            comp_inds = TMP_STUDY.cluster(k).comps;
            for i = 1:length(subj_inds)
                subj_i = subj_inds(i);
                comp_i = comp_inds(i);
                subj_char = TMP_STUDY.datasetinfo(subj_i).subject;
                %-
%                 subj_icscores = [load_dir filesep subj_char filesep 'ICA' filesep 'IC_score.jpg'];
                subj_topos = [load_dir filesep subj_char filesep 'ICA' filesep 'potential_brain_components_allcritera_customElectrode.jpg'];
                subj_psds = [load_dir filesep subj_char filesep 'ICA' filesep 'KEEP_IC_PSD.jpg'];
                subj_powpow = [load_dir filesep subj_char filesep 'ICA' filesep 'powpowcat.jpg'];
                subj_ersp = cell(length(TMP_STUDY.design),1);
                for des_i = 1:length(TMP_STUDY.design)
                    subj_ersp{des_i} = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('%s_within_des%i_cl%i_stats.jpg',subj_char,des_i,k)];
                    subj_psds{des_i} = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_psd_des%i_ic%i',subj_char,des_i,comp_i)];
                end
                %-
                newSlide = slides.AddSlide(1,layout);
                Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,40,400,70);
                Title1.TextFrame.TextRange.Text = sprintf('Subject %s\n(K=%i) Cluster %i',subj_char,clust_i,k);
                %- Subject Specific ERSP Plots
                im_scale = 0.35;
                vertical_move = 0;
                for des_i = 1:length(TMP_STUDY.design)
                    if exist(subj_ersp{des_i},'file')
                        newSlide.Shapes.AddPicture(subj_ersp{des_i}, 'msoFalse', 'msoTrue',0,0+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
                    end
                    vertical_move = vertical_move + 729*im_scale; 
                end
                %- Subject Specific IC Criteria
                TOP_DIST = 70;
                im_scale = 0.275;
                newSlide = slides.AddSlide(1,layout);
                Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
                Title1.TextFrame.TextRange.Text = sprintf('Subject %s\n(K=%i) Cluster %i',subj_char,clust_i,k);
                if exist(subj_topos,'file')
                    newSlide.Shapes.AddPicture(subj_topos, 'msoFalse', 'msoTrue',0,0+TOP_DIST,1167*im_scale,875*im_scale); %Left, top, width, height
                end
                if exist(subj_psds,'file')
                    newSlide.Shapes.AddPicture(subj_psds, 'msoFalse', 'msoTrue',875*im_scale+70,0+TOP_DIST,2250*im_scale,1500*im_scale); %Left, top, width, height
                end
                if exist(subj_powpow,'file')
                    newSlide.Shapes.AddPicture(subj_powpow, 'msoFalse', 'msoTrue',0,0+TOP_DIST+875*im_scale,1167*im_scale,875*im_scale); %Left, top, width, height
                end
            end 
            %}
            %## (SLIDES 4) Subject Specific ERSPs & Criteria
%             im_scale = 0.225;
            subj_inds = TMP_STUDY.cluster(k).sets;
            comp_inds = TMP_STUDY.cluster(k).comps;
            for i = 1:length(subj_inds)
                subj_i = subj_inds(i);
                comp_i = comp_inds(i);
                subj_char = TMP_STUDY.datasetinfo(subj_i).subject; %TMP_STUDY.subject{subj_i};
                fprintf('Making Slide for subject %s\n',subj_char);
                %-
%                 subj_icscores = [load_dir filesep subj_char filesep 'ICA' filesep 'IC_score.jpg'];
%                 subj_topos = [load_dir filesep subj_char filesep 'ICA' filesep 'potential_brain_components_allcritera_customElectrode.jpg'];
%                 subj_psds = [load_dir filesep subj_char filesep 'ICA' filesep 'KEEP_IC_PSD.jpg'];
%                 subj_powpow = [load_dir filesep subj_char filesep 'ICA' filesep 'powpowcat.jpg'];
                subj_topos = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_topo_ic%i.jpg',subj_char,comp_i)];
%                 subj_psds = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_topo_ic%i',subj_char,comp_i)];
                subj_dip_top = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_dip_top_ic%i.jpg',subj_char,comp_i)];
                subj_dip_cor = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_dip_coronal_ic%i.jpg',subj_char,comp_i)];
                subj_dip_sag = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_dip_sagittal_ic%i.jpg',subj_char,comp_i)];
                subj_psds = cell(length(TMP_STUDY.design),1);
                subj_ersp = cell(length(TMP_STUDY.design),1);
                for des_i = 1:length(TMP_STUDY.design)
                    subj_ersp{des_i} = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('%s_within_des%i_cl%i_stats.jpg',subj_char,des_i,k)];
                    subj_psds{des_i} = [cluster_dir filesep sprintf('%i',k) filesep sprintf('%s_psd_des%i_ic%i.jpg',subj_char,des_i,comp_i)];
                end
                %- (SUBSLIDE ERSPS)
                TOP_DIST = 70;
                newSlide = slides.AddSlide(1,layout);
                %- Subject Specific ERSP Plots
                im_scale = 0.35;
                vertical_move = 0;
                for des_i = 1:length(TMP_STUDY.design)
                    if exist(subj_ersp{des_i},'file')
                        newSlide.Shapes.AddPicture(subj_ersp{des_i}, 'msoFalse', 'msoTrue',100,0+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
                    end
                    vertical_move = vertical_move + 729*im_scale; 
                end
                %- topo
                im_scale = 0.275;
                if exist(subj_topos,'file')
                    newSlide.Shapes.AddPicture(subj_topos, 'msoFalse', 'msoTrue',0-50,0+TOP_DIST+20,781*im_scale,547*im_scale); %Left, top, width, height
                end
                Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,40,400,70);
                Title1.TextFrame.TextRange.Text = sprintf('Subject %s\n(K=%i) Cluster %i',subj_char,clust_i,k);
                %- (SUBSLIDE PSDS,TOPOS,DIPS) Subject Specific IC Criteria
                shift = 0;
                TOP_DIST = 70;
%                 im_scale = 0.275;
                newSlide = slides.AddSlide(1,layout);
                Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
                Title1.TextFrame.TextRange.Text = sprintf('Subject %s\n(K=%i) Cluster %i',subj_char,clust_i,k);
                %- topo
                im_scale = 0.3;
                if exist(subj_topos,'file')
                    newSlide.Shapes.AddPicture(subj_topos, 'msoFalse', 'msoTrue',0-30,0+TOP_DIST+10,781*im_scale,547*im_scale); %Left, top, width, height
                end
                %- dipoles
                im_scale = 0.35;
                if exist(subj_dip_top,'file')
                    %## crop the image
                    newSlide.Shapes.AddPicture(subj_dip_top, 'msoFalse', 'msoTrue',300,0+TOP_DIST-70,469*im_scale,547*im_scale); %Left, top, width, height
                end
                if exist(subj_dip_cor,'file')
                    newSlide.Shapes.AddPicture(subj_dip_cor, 'msoFalse', 'msoTrue',300+547*im_scale,0+TOP_DIST-70,469*im_scale,547*im_scale); %Left, top, width, height
                end
                if exist(subj_dip_sag,'file')
                    newSlide.Shapes.AddPicture(subj_dip_sag, 'msoFalse', 'msoTrue',300+547*im_scale*2,0+TOP_DIST-70,469*im_scale,547*im_scale); %Left, top, width, height
                end
                %- psds
                im_scale = 0.6;
                for des_i = 1:length(TMP_STUDY.design)
                    if exist(subj_psds{des_i},'file')
                        newSlide.Shapes.AddPicture(subj_psds{des_i}, 'msoFalse', 'msoTrue',150+shift,0+TOP_DIST+150,656*im_scale,583*im_scale); %Left, top, width, height
                    end
                    shift = shift + 656*im_scale;
                end
            end
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.325;
            TOP_DIST = 40;
            LEFT_DIST = -30;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
%                 if des_i == 1
%                     LEFT_DIST = -200;
%                 else
%                     LEFT_DIST = 0;
%                 end
%                 ersp1 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot1_des%i_cl%i_stats_lim60.jpg',des_i,k)];
%                 ersp2 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_notfull_stats_compare-flat_des%i_cl%i.jpg',des_i,k)];
%                 ersp3 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_stats_compare*.jpg']);
%                 ersp3 = [ersp3.folder filesep ersp3.name];
                ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp4,'file')
                    newSlide.Shapes.AddPicture(ersp4, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 729*im_scale; 
            end
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.33;
            TOP_DIST = 40;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
                if des_i == 1
                    LEFT_DIST = -220;
                else
                    LEFT_DIST = 0;
                end
%                 ersp1 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot1_des%i_cl%i_stats_lim60.jpg',des_i,k)];
                ersp2 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_notfull_stats_compare*.jpg']);
                ersp2 = [ersp2.folder filesep ersp2.name];
%                 ersp3 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_stats_compare*.jpg']);
%                 ersp3 = [ersp3.folder filesep ersp3.name];
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp2,'file')
                    newSlide.Shapes.AddPicture(ersp2, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
                end
%                 if exist(ersp3,'file')
%                     newSlide.Shapes.AddPicture(ersp3, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
%                 end
                vertical_move = vertical_move + 729*im_scale; 
            end
            %## (SLIDE 3a) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.33;
            TOP_DIST = 40;
            LEFT_DIST = -90;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
                ersp1 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot1_des%i_cl%i_stats_lim60.jpg',des_i,k)];
%                 ersp2 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_notfull_stats_compare-flat_des%i_cl%i.jpg',des_i,k)];
%                 ersp3 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_compare-flat_des%i_cl%i.jpg',des_i,k)];
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp1,'file')
                    newSlide.Shapes.AddPicture(ersp1, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 729*im_scale; 
            end
            %## (SLIDE 2) Add power spectra across conditions
            im_scale = 0.55;
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',50,10,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            if exist(topo,'file')
                newSlide.Shapes.AddPicture(psd1, 'msoFalse', 'msoTrue',0,50,875*im_scale,750*im_scale); %Left, top, width, height
            end
            if exist(topo,'file')
                newSlide.Shapes.AddPicture(psd2, 'msoFalse', 'msoTrue',875*im_scale,50,875*im_scale,750*im_scale); %Left, top, width, height
            end
            %## (SLIDE 1) Add topo plots and dips to slide
            im_scale = 0.225;
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',1200*im_scale*2,30,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            if exist(D1,'file')
                newSlide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',0,0-40,1200*im_scale,1575*im_scale); %Left, top, width, height
            end
            if exist(D2,'file')
                newSlide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',1200*im_scale,0-40,1200*im_scale,1575*im_scale); %Left, top, width, height
            end
            if exist(D3,'file')
                newSlide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',0,1575*im_scale-100,1200*im_scale,1575*im_scale); %Left, top, width, height
            end
            if exist(topo,'file')
                newSlide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',1200*im_scale*2,100,1200*im_scale*1.2,1575*im_scale*1.2); %Left, top, width, height
            end
            %## (SLIDE TITLE)
            %- title slide
            newSlide = slides.AddSlide(1,layout);
            set(newSlide.Shapes.Title.Textframe.Textrange, 'Text', sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel));
            subj_inds = TMP_STUDY.cluster(k).sets;
            subj_char = {TMP_STUDY.datasetinfo(subj_inds).subject};
            subj_char = strjoin(subj_char,', ');
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('Cluster Label: %s\nSubjects in cluster: %s\n',TMP_STUDY.cluster(k).name,subj_char);
        end
    end
    
    ppt.ActivePresentation.SaveAs([cluster_dir filesep sprintf('summary_cluster_report.pptx')]);
%     [Status,Message] = saveppt([ppt_savefPath filesep sprintf('summary_cluster_report.pptx')],'converttopdf');
%     fprintf('%s\n',Status);
    % Close Powerpoint and delete the object
    ppt.ActivePresentation.Close;
    ppt.Quit;
    ppt.delete;
    fprintf('\nResults powerpoint saved here:\n%s\n',cluster_dir);
end