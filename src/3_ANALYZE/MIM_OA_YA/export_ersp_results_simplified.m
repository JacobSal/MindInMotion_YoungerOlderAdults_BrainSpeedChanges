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
%- datetime override
% dt = '10252023_MIM_OAN70_noslowwalkers_gait_powpow0p25';
% dt = '10302023_MIM_OAN70_noslowwalkers_gait_iccREMG0p4_powpow0p1';
% dt = '10302023_MIM_OAN70_newnormalize_iccREMG0p4_powpow0p1';
% dt = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
dt = '11302023_MIM_OAN70_antsnormalize_iccREMG0p3_powpow0p1';

%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%- (EDIT!) convert SUB_DIR
% SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\10052023_MIM_OAN70_noslowwalkers_gait\cluster';
% SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\10252023_MIM_OAN70_noslowwalkers_gait_powpow0p25\cluster';
% SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\10252023_MIM_OAN70_noslowwalkers_gait_powpow0p20\cluster';
SUB_DIR = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [12];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_5' filesep '12']};
CLUSTER_FILES = {'cl_inf_12.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
SUB_GROUP_FNAME = []; %'H3000'; %[]; %'H2000';
SUB_GROUP_FNAME_REGEX = []; %'H3000''s'; %[]; %'H2000''s';
CLUSTER_CLIM_MATCH = [];
STUDY_GROUP_DESI = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable2','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
%% ===================================================================== %%
for k_i = 1:length(CLUSTER_DIRS)
    %## CREATE & GRAB DIRECTORIES
    %- convert cluster directory
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    if ~isempty(SUB_GROUP_FNAME_REGEX)
        spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
        plot_store_dir = [cluster_dir filesep 'plots_out' filesep SUB_GROUP_FNAME];
    else
        spec_data_dir = [cluster_dir filesep 'spec_data'];
        plot_store_dir = [cluster_dir filesep 'plots_out'];
    end
    if ~exist(spec_data_dir,'dir')
        error('spec_data dir does not exist');
    end
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    %- convert study directory
    if ~ispc
        cluster_study_dir = convertPath2UNIX(CLUSTER_STUDY_DIRS{k_i});
    else
        cluster_study_dir = convertPath2Drive(CLUSTER_STUDY_DIRS{k_i});
    end
    %## LOAD STUDY
    if ~ispc
        tmp = load('-mat',[cluster_study_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
        TMP_STUDY = tmp.STUDY;
    else
        tmp = load('-mat',[cluster_study_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
        TMP_STUDY = tmp.STUDY;
    end
    %* assign studies
    TMP_STUDY.design = STUDY_GROUP_DESI;
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
%     cl_dir = dir([cluster_dir filesep '*.jpg']);
%     inds = cellfun(@(x) contains(x,'Cluster_topo_avg'),{cl_dir.name});
%     topo_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
%     inds = cellfun(@(x) contains(x,'dipplot_alldipspc_top'),{cl_dir.name});
%     diptop_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
%     inds = cellfun(@(x) contains(x,'dipplot_alldipspc_sagittal'),{cl_dir.name});
%     dipsag_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
%     inds = cellfun(@(x) contains(x,'dipplot_alldipspc_coronal'),{cl_dir.name});
%     dipcor_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
%     inds = cellfun(@(x) contains(x,'allSpecPlot_des1'),{cl_dir.name});
%     psd_des1_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
%     inds = cellfun(@(x) contains(x,'allSpecPlot_des2'),{cl_dir.name});
%     psd_des2_files = cl_dir(inds); %{cluster_dir(inds).folder filesep cluster_dir(inds).name};
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
%         ind1 = cellfun(@(x) strcmp(x,sprintf('Cluster_topo_avg%i.jpg',k)),{topo_files.name});
%         ind2 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldipspc_top%i.jpg',k)),{diptop_files.name});
%         ind3 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldipspc_sagittal%i.jpg',k)),{dipsag_files.name});
%         ind4 = cellfun(@(x) strcmp(x,sprintf('dipplot_alldipspc_coronal%i.jpg',k)),{dipcor_files.name});
% %         ind5 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des1_cl%i.jpg',k)),{psd_des1_files.name});
% %         ind6 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des2_cl%i.jpg',k)),{psd_des2_files.name});
%         ind5 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des1_cl%i.pdf',k)),{psd_des1_files.name});
%         ind6 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des2_cl%i.pdf',k)),{psd_des2_files.name});
        topo = [cluster_dir filesep sprintf('%i_Cluster_topo_avg.jpg',k)];
        D1 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.jpg',k)];
        D2 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.jpg',k)];
        D3 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.jpg',k)];
%         ind5 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des1_cl%i.jpg',k)),{psd_des1_files.name});
%         ind6 = cellfun(@(x) strcmp(x,sprintf('allSpecPlot_des2_cl%i.jpg',k)),{psd_des2_files.name});
%         psd1 = [cluster_dir filesep sprintf('allSpecPlot_des1_cl%i.jpg',k)];
%         psd2 = [cluster_dir filesep sprintf('allSpecPlot_des2_cl%i.jpg',k)];

        %*
        % opts: (design 1, cluster 3) erspplot1_des1_cl3_stats_lim60,
        % erspplottype2_notfull_stats_compare-flat_des1_cl3,
        % erspplottype2_stats_compare-flat_des1_cl3,
        % erspplottype2_stats_des1_cl3,
        % erspplottype2_stats_notfull_des1_cl3,
        % erspplottype3_orig_des1_cl3_lim60, && subject specific e.g., H2008_within_des1_cl3_stats
        ersp_fpath = [cluster_dir filesep 'plots_out' filesep sprintf('%i',k)];
        ersp_f = [];
        %- checks
%         chk1 = sum(ind1)>=1 && sum(ind2)>=1 && sum(ind3)>=1 && sum(ind4)>=1 && sum(ind5)>=1 && sum(ind6)>=1;
        chk2 = exist(ersp_fpath,'dir');
%         disp(chk1 && chk2)
        if chk2 %chk1 %&& chk2
            %-
%             D1 = [diptop_files(ind2).folder filesep diptop_files(ind2).name];
%             D2 = [dipsag_files(ind3).folder filesep dipsag_files(ind3).name];
%             D3 = [dipcor_files(ind4).folder filesep dipcor_files(ind4).name];
%             topo = [topo_files(ind1).folder filesep topo_files(ind1).name];
%             psd1 = [psd_des1_files(ind5).folder filesep psd_des1_files(ind5).name];
%             psd2 = [psd_des2_files(ind6).folder filesep psd_des2_files(ind6).name];
            
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 50;
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
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_stats_notfull_des%i_cl%i*.jpg',des_i,k)]);
%                 ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_stats_notfull_des%i_cl%i.pdf',des_i,k)]);
                ersp_f = [ersp_f(1).folder filesep ersp_f(1).name];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,4075*im_scale,1312*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 1312*im_scale+10; 
            end
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
%                 ersp1 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot1_des%i_cl%i_stats_lim60.jpg',des_i,k)];
%                 ersp2 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_notfull_stats_compare*.jpg']);
%                 ersp2 = [ersp2.folder filesep ersp2.name];
%                 ersp3 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_stats_compare*.jpg']);
%                 ersp3 = [ersp3.folder filesep ersp3.name];
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
                ersp_f = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp_f = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.pdf',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,4075*im_scale,1312*im_scale); %Left, top, width, height
                end
%                 if exist(ersp3,'file')
%                     newSlide.Shapes.AddPicture(ersp3, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
%                 end
                vertical_move = vertical_move + 1312*im_scale+10; 
            end
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 100;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
%                 if des_i == 1
%                     LEFT_DIST = -220;
%                 else
%                     LEFT_DIST = 0;
%                 end
%                 ersp1 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot1_des%i_cl%i_stats_lim60.jpg',des_i,k)];
                ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_notfull_stats_compare*.jpg']);
                [val,ind] = max([ersp_f.datenum]);
                ersp_f = [ersp_f(ind).folder filesep ersp_f(ind).name];
%                 ersp3 = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_stats_compare*.jpg']);
%                 ersp3 = [ersp3.folder filesep ersp3.name];
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,2631*im_scale,1311*im_scale); %Left, top, width, height
                end
%                 if exist(ersp3,'file')
%                     newSlide.Shapes.AddPicture(ersp3, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3083*im_scale,729*im_scale); %Left, top, width, height
%                 end
                vertical_move = vertical_move + 1311*im_scale+10; 
            end
            %## (SLIDE 3b) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 100;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
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
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype3_notfull_stats_compare*.jpg']);
                [val,ind] = max([ersp_f.datenum]);
                ersp_f = [ersp_f(ind).folder filesep ersp_f(ind).name];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,2631*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 1311*im_scale+10; 
            end
            %## (SLIDE 3a) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot_des%i_cl%i_bootstats.jpg',des_i,k)];
%                 ersp2 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_notfull_stats_compare-flat_des%i_cl%i.jpg',des_i,k)];
%                 ersp3 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_compare-flat_des%i_cl%i.jpg',des_i,k)];
%                 ersp4 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',des_i,k)];
%                 ersp5 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
%                 ersp6 = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_orig_des%i_cl%i_lim60.jpg',des_i,k)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 1311*im_scale+10; 
            end
            %## (SLIDE 2) Add power spectra across conditions
            vertical_move = 0;
            im_scale = 0.25;
            TOP_DIST = 50;
            LEFT_DIST = 50;
            newSlide = slides.AddSlide(1,layout);
            fooof_fpath = [spec_data_dir filesep 'psd_calcs'];
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',50,10,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
            for des_i = 1:length(TMP_STUDY.design)
                psd_plot = dir([fooof_fpath filesep sprintf('TopoPSD_Violins_cl%i_des%i.jpg',k,des_i)]);
                psd_plot = [psd_plot.folder filesep psd_plot.name];
                if exist(psd_plot,'file')
                    newSlide.Shapes.AddPicture(psd_plot, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3411*im_scale,768*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 768*im_scale+10; 
            end
            %## (SLIDE 1) Add topo plots and dips to slide
            im_scale = 0.225;
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',1200*im_scale*2,30,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel);
            if exist(D1,'file')
                newSlide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',0,0+20,812*im_scale,974*im_scale); %Left, top, width, height
            end
            if exist(D2,'file')
                newSlide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',812*im_scale,0+20,812*im_scale,813*im_scale); %Left, top, width, height
            end
            if exist(D3,'file')
                newSlide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',812*im_scale*2,0+20,974*im_scale,813*im_scale); %Left, top, width, height
            end
            if exist(topo,'file')
                newSlide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',0,813*im_scale+30,834*im_scale,946*im_scale); %Left, top, width, height
            end
            %## (SLIDE TITLE)
            %- title slide
            newSlide = slides.AddSlide(1,layout);
            Txt1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Txt1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,TMP_STUDY.cluster(k).analabel);
%             set(newSlide.Shapes.Title.Textframe.Textrange, 'Text', sprintf('(K=%i) Cluster %i, %s',clust_i,k,TMP_STUDY.cluster(k).analabel));
            subj_inds = TMP_STUDY.cluster(k).sets;
            subj_char = {TMP_STUDY.datasetinfo(subj_inds).subject};
            subj_char = strjoin(subj_char,', ');
            Txt2 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,80,400,70);
            Txt2.TextFrame.TextRange.Text = sprintf('Cluster Label: %s\nSubjects in cluster: %s\n',TMP_STUDY.cluster(k).name,subj_char);
            %## ATLAS INFORMATION
%             SUBJ_ATLAS_INT = 1;
%             atlas_dir = dir([cluster_dir filesep sprintf('%i',k) filesep '*atlasinf*']);
%             atlas_fPath = [atlas_dir(SUBJ_ATLAS_INT).folder filesep atlas_dir(SUBJ_ATLAS_INT).name];
%             atlas_inf = readtable(atlas_fPath);
%             while size(atlas_inf,2)~=3 && SUBJ_ATLAS_INT<length(atlas_dir)
%                 SUBJ_ATLAS_INT = SUBJ_ATLAS_INT+1;
%                 atlas_fPath = [atlas_dir(SUBJ_ATLAS_INT).folder filesep atlas_dir(SUBJ_ATLAS_INT).name];
%                 tmp = readtable(atlas_fPath);
%                 atlas_inf = tmp;
%                 atlases = tmp(:,3);
%             end
            %- get atlas counts
%             atlas_inf = cell(length(subj_inds),1);
%             for i = 1:length(subj_inds)
%                 atlas_fPath = [atlas_dir(i).folder filesep atlas_dir(i).name];
%                 tmp = readtable(atlas_fPath);
%                 if size(tmp,2)==3
%                     atlas_inf{i} = {tmp.subject_dipole{:}};
%                  	atlases = tmp;
%                 end
%             end
%             atlas_inf = atlas_inf(~cellfun(@isempty,atlas_inf));
%             word_vec = cat(2,atlas_inf{:});
%             words = unique(word_vec);
%             ahist_fpath = [cluster_dir filesep sprintf('cluster%i_histogram_atlas_counts.jpg',k)];
%             if ~exist(ahist_fpath,'file')
%                 %## HISTOGRAM
%                 fig = figure('Position',[10,100,1200,900]);
%                 histogram(categorical(word_vec),words)
%                 ax = gca;
%                 set(gca,'TickLabelInterpreter','none')
%                 ylabel('Counts of Atlas Labels');
%                 xlabel('Atlas Determination')
%                 %## ADD HISTOGRAM
%                 im_scale = 0.25;
%                 saveas(fig,ahist_fpath);
%             end
%             if exist(ahist_fpath,'file')
%                 newSlide.Shapes.AddPicture(ahist_fpath, 'msoFalse', 'msoTrue',400,100,1875*im_scale*1.2,1406*im_scale*1.2); %Left, top, width, height
%             end
            %- add label
%             Txt2 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',300,00,650,100);
%             Txt2.TextFrame.TextRange.Text = sprintf('Atlases Used: %s',strjoin({atlases.atlas{:}},', '));
            
            %- table
            %{
            table = invoke(newSlide.Shapes, 'AddTable', size(atlas_inf,1)+1, size(atlas_inf,2), 400, 0, 400, 70);    % create Table
            table.Table.Cell(1,1).Shape.TextFrame.TextRange.Text = 'Centroid'; 
            table.Table.Cell(1,2).Shape.TextFrame.TextRange.Text = sprintf('Subject %s',TMP_STUDY.datasetinfo(subj_inds(SUBJ_ATLAS_INT)).subject);
            table.Table.Cell(1,3).Shape.TextFrame.TextRange.Text = 'Atlas';
            for i = 1:size(atlas_inf,1)
                centr = string(atlas_inf{i,1});
                subj =  string(atlas_inf{i,2});
                atlas = string(atlas_inf{i,3});
                table.Table.Cell(i+1,1).Shape.TextFrame.TextRange.Text = centr;
                table.Table.Cell(i+1,2).Shape.TextFrame.TextRange.Text = subj;
                table.Table.Cell(i+1,3).Shape.TextFrame.TextRange.Text = atlas;
            end
            for i = 1:size(atlas_inf,2)
                table.Table.Columns.Item(i).Width = 180;
            end
            %}
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