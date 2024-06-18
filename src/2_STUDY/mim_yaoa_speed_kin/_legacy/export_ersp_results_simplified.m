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
% dt = '11302023_MIM_OAN70_antsnormalize_iccREMG0p3_powpow0p1';
% dt = '12012023_OAYA104_icc0p65-0p4_changparams';
% dt = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
dt = '01232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3';
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
    %## PPT MAKE
    Image = [];
    cl_names = {TMP_STUDY.cluster.name};
    main_cl_inds = main_cl_inds(2:end);
   %%
    for k = flip(main_cl_inds,2) %length(cl_names):-1:1 %1:length(cl_names)
        %## atlas calc
        ahist_fpath = [cluster_dir filesep sprintf('cluster%i_histogram_atlas_counts.jpg',k)];
        if ~exist(ahist_fpath,'file') %|| true
            SUBJ_ATLAS_INT = 1;
            subj_inds = TMP_STUDY.cluster(k).sets;
            atlas_dir = dir([cluster_dir filesep sprintf('%i',k) filesep '*atlasinf*']);
            atlas_fPath = [atlas_dir(SUBJ_ATLAS_INT).folder filesep atlas_dir(SUBJ_ATLAS_INT).name];
            atlas_inf = readtable(atlas_fPath,'Delimiter',',');
            %- get atlas counts
            at_subj_char = cell(length(subj_inds)*size(atlas_inf,1),1);
            at_dip_char = cell(length(subj_inds)*size(atlas_inf,1),1);
            at_atlas_char = cell(length(subj_inds)*size(atlas_inf,1),1);
            cnta = 1;
            for i = 1:length(subj_inds)
                atlas_fPath = [atlas_dir(i).folder filesep atlas_dir(i).name];
                tmp = readtable(atlas_fPath,'Delimiter',','); %readtable(atlas_fPath);
                if size(tmp,2)==3
                    for a_i = 1:size(tmp,1)
    %                         atlas_inf(i,a_i,1) = tmp.subject_dipole(a_i);
    %                         atlas_inf(i,a_i,2) = tmp.atlas(a_i);
                        at_subj_char{cnta} = TMP_STUDY.datasetinfo(subj_inds(i)).subject;
                        at_dip_char(cnta) = tmp.subject_dipole(a_i);
                        at_atlas_char(cnta) = tmp.atlas(a_i);
                        cnta = cnta + 1;
                    end
                end
            end
            at_subj_char = categorical(at_subj_char(~cellfun(@isempty,at_subj_char)));
            at_dip_char = categorical(at_dip_char(~cellfun(@isempty,at_dip_char)));
            at_atlas_char = categorical(at_atlas_char(~cellfun(@isempty,at_atlas_char)));
            atlas_table = table(at_subj_char,at_dip_char,at_atlas_char);
            par_save(atlas_table,cluster_dir,sprintf('cluster%i_atlas_table.mat',k));
            %- histo setup
%             ahist_fpath = [cluster_dir filesep sprintf('cluster%i_histogram_atlas_counts.jpg',k)];
            %## HISTOGRAM
            at_chars = unique(atlas_table.at_atlas_char);
            colors = linspecer(length(at_chars));
            fig = figure('Position',[10,100,1620,720]);
            hold on;
            bx_iter = [];
            bx_hold = [];
            bin_hold = {};
            h1_hold = [];
            for i = 1:length(at_chars)
                at_ind = at_chars(i)==atlas_table.at_atlas_char;
                tmp_in = atlas_table.at_dip_char(at_ind);
                if i == 1
                    [h1,bin] = histcounts(tmp_in);
                    chk = h1~=0;
                    h1 = h1(chk);
                    bin = bin(chk);
                    bx = 1:length(h1);
                    bb = bar(bx,h1);
                    ax = gca;cfname_path
                    set(gca,'TickLabelInterpreter','none')
                    ylabel('Counts of Atlas Labels');
                    xlabel('Atlas Determination')
                    bb.CData = repmat(colors(i,:),[length(h1),1]);
                    bx_iter = length(h1);
                    bin_hold = [bin_hold, bin];
                    bx_hold = [bx_hold, bx];
                    h1_hold = [h1_hold, h1];
                else
                    [h1,bin] = histcounts(tmp_in);
                    chk = h1~=0;
                    h1 = h1(chk);
                    bin = bin(chk);
                    bx = (1:length(h1))+bx_iter;
                    bb = bar(ax,bx,h1);
                    bb.CData = repmat(colors(i,:),[length(h1),1]);
                    bx_iter = bx_iter + length(h1);
                    bin_hold = [bin_hold, bin];
                    bx_hold = [bx_hold, bx];
                    h1_hold = [h1_hold, h1];
                end
            end
            set(gca,'xtick',bx_hold,'xticklabel',bin_hold);
            xtickangle(60);
            legend(at_chars,'location', 'bestoutside','interpreter','none')
            hold off;
            %## ADD HISTOGRAM
            saveas(fig,ahist_fpath);
            close(fig);
            chk = ~strcmp('no_label_found',bin_hold);
            bin_hold = bin_hold(chk);
            h1_hold = h1_hold(chk);
            [val] = max(h1_hold);
            ind = find(h1_hold == val);
            if length(ind) > 1
                analabel_new = [sprintf('%s',bin_hold{ind(1)}),sprintf('-%s',bin_hold{ind(2:end)}),sprintf('-counts%i',val)];
            else
                analabel_new = sprintf('%s-N%i',bin_hold{ind},val);
            end
        else
            analabel_new = TMP_STUDY.cluster(k).analabel;
        end
        
        %##
        topo = [cluster_dir filesep sprintf('%i_Cluster_topo_avg.jpg',k)];
        D1 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.jpg',k)];
        D2 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.jpg',k)];
        D3 = [cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.jpg',k)];

        %*
        spca_fpath = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'spca'];
        ersp_fpath = [cluster_dir filesep 'plots_out' filesep sprintf('%i',k)];
        ersp_f = [];
        %- checks
        chk2 = exist(ersp_fpath,'dir');
        if chk2
            %## (SLIDE 3f_spca) Cluster Level GPM (common base)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpm_com.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpm_com.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move - 1311*im_scale+10; 
            end
            %## (SLIDE 3e_spca) Cluster Level ERSPs (common base)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersp_com.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersp_com.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move - 1311*im_scale+10; 
            end
            %## (SLIDE 3d_spca) Cluster Level GPM (clean)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpm.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpm.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move- 1311*im_scale+10; 
            end
            %## (SLIDE 3c_spca) Cluster Level ERSPs (clean)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersp.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersp.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move - 1311*im_scale+10; 
            end
            %## (SLIDE 3b_spca) Cluster Level GPM (original)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpmorig.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_gpmorig.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move - 1311*im_scale+10; 
            end
            %## (SLIDE 3a_spca) Cluster Level ERSPs (original)
            im_scale = 0.19;
            vertical_move = 1311*im_scale+10-20;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersporig.jpg',k,des_i)];
%                 ersp_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersporig.tiff',k,des_i)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3365*im_scale,1311*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move - 1311*im_scale+10; 
            end
            %## (SLIDE 3e) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype3_stats_notfull_des%i_cl%i*.jpg',des_i,k)]);
                ersp_f = [ersp_f(1).folder filesep ersp_f(1).name];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,4075*im_scale,1312*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 1312*im_scale+10; 
            end
            %## (SLIDE 3d) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 50;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',des_i,k)];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,4075*im_scale,1312*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 1312*im_scale+10; 
            end
            %## (SLIDE 3c) Cluster Level ERSPs
            vertical_move = 0;
            im_scale = 0.19;
            TOP_DIST = 30;
            LEFT_DIST = 100;
            %-
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = dir([ersp_fpath filesep sprintf('%i',des_i) filesep 'erspplottype2_notfull_stats_compare*.jpg']);
                [val,ind] = max([ersp_f.datenum]);
                ersp_f = [ersp_f(ind).folder filesep ersp_f(ind).name];
                if exist(ersp_f,'file')
                    newSlide.Shapes.AddPicture(ersp_f, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,2631*im_scale,1311*im_scale); %Left, top, width, height
                end
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
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
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
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                ersp_f = [ersp_fpath filesep sprintf('%i',des_i) filesep sprintf('erspplot_des%i_cl%i_bootstats.jpg',des_i,k)];
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
            Title1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
            for des_i = 1:length(TMP_STUDY.design)
                psd_plot = dir([fooof_fpath filesep sprintf('TopoPSD_Violins_cl%i_des%i.jpg',k,des_i)]);
                psd_plot = [psd_plot.folder filesep psd_plot.name];
                if exist(psd_plot,'file')
                    newSlide.Shapes.AddPicture(psd_plot, 'msoFalse', 'msoTrue',0+LEFT_DIST,0+TOP_DIST+vertical_move,3411*im_scale,768*im_scale); %Left, top, width, height
                end
                vertical_move = vertical_move + 768*im_scale+10; 
            end
            %## (SLIDE 2) ATLAS INFORMATION
            im_scale = 0.28;
            newSlide = slides.AddSlide(1,layout);
            if exist(ahist_fpath,'file')
                newSlide.Shapes.AddPicture(ahist_fpath, 'msoFalse', 'msoTrue',0,100,2531*im_scale*1.2,1125*im_scale*1.2); %Left, top, width, height
            end
            %## (SLIDE TITLE)
            %## (SLIDE 1) Add topo plots and dips to slide
            im_scale = 0.225;
            newSlide = slides.AddSlide(1,layout);
            Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',1200*im_scale*2+50,30,400,70);
            Title1.TextFrame.TextRange.Text = sprintf('(K=%i) Cluster %i, %s',clust_i,k,analabel_new);
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
            %- title slide
%             newSlide = slides.AddSlide(1,layout);
%             Txt1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',400,200,400,70);
%             Txt1.TextFrame.TextRange.Text = sprintf('(N=%i) Cluster %i, %s',length(TMP_STUDY.cluster(k).sets),k,analabel_new);
%             set(newSlide.Shapes.Title.Textframe.Textrange, 'Text', sprintf('(K=%i) Cluster %i, %s',clust_i,k,analabel_new));
            subj_inds = TMP_STUDY.cluster(k).sets;
            subj_char = {TMP_STUDY.datasetinfo(subj_inds).subject};
            subj_char = strjoin(subj_char,', ');
            Txt2 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',400,270,400,70);
            Txt2.TextFrame.TextRange.Text = sprintf('Cluster Label: %s\nSubjects in cluster: %s\n',TMP_STUDY.cluster(k).name,subj_char);
            %- add label
%             Txt2 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',300,00,650,100);
%             Txt2.TextFrame.TextRange.Text = sprintf('Atlases Used: %s',strjoin({atlases.atlas{:}},', '));
            
            %- table
%             table = invoke(newSlide.Shapes, 'AddTable', size(atlas_inf,1)+1, size(atlas_inf,2), 400, 0, 400, 70);    % create Table
%             table.Table.Cell(1,1).Shape.TextFrame.TextRange.Text = 'Centroid'; 
%             table.Table.Cell(1,2).Shape.TextFrame.TextRange.Text = sprintf('Subject %s',TMP_STUDY.datasetinfo(subj_inds(SUBJ_ATLAS_INT)).subject);
%             table.Table.Cell(1,3).Shape.TextFrame.TextRange.Text = 'Atlas';
%             for i = 1:size(atlas_inf,1)
%                 centr = string(atlas_inf{i,1});
%                 subj =  string(atlas_inf{i,2});
%                 atlas = string(atlas_inf{i,3});
%                 table.Table.Cell(i+1,1).Shape.TextFrame.TextRange.Text = centr;
%                 table.Table.Cell(i+1,2).Shape.TextFrame.TextRange.Text = subj;
%                 table.Table.Cell(i+1,3).Shape.TextFrame.TextRange.Text = atlas;
%             end
%             for i = 1:size(atlas_inf,2)
%                 table.Table.Columns.Item(i).Width = 180;
%             end
            
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