%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_ersp_clusters.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
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
%- define the directory to the src folderd
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
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('oa');
% fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
% fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
dt = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
%- study group and saving
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
study_fName_1 = 'all_comps_study';
study_fName_2 = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
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
k_i = 1;
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
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
    TMP_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
    TMP_STUDY = tmp.STUDY;
end
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(TMP_STUDY);
%% ===================================================================== %%

%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
%## ersp plot per cluster per condition
TMP_STUDY = pop_statparams(TMP_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
TMP_STUDY = pop_erspparams(TMP_STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
TMP_STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [TMP_STUDY] = std_makedesign(TMP_STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%%
%##
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%##
condition_gait = {'flat','low','med','high','0p25','0p5','0p75','1p0'};
% condition_pairs = {{'flat','low','med','high'},...
%     {'0p25','0p5','0p75','1p0'}};
allersp_clusters = par_load(spec_data_dir,'spca_ersp_clusters.mat');
allgpm_clusters = par_load(spec_data_dir,'spca_gpm_clusters.mat');
allersp = cell(4,length(GROUP_NAMES));
allgpm = cell(4,length(GROUP_NAMES));
%##
cond_iters = {1:4,5:8};

%%
for cl_i = main_cl_inds(2:13)
    %##
    dip1 = TMP_STUDY.cluster(cl_i).dipole;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    cfg              = [];
    cfg.roi        = dip1.posxyz;
    cfg.output     = 'multiple';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        % disp(labels.count(labels.count ~= 0))
        [val, indx] = max(label_i.count);
        if strcmp(label_i.name(indx),'no_label_found')
            sub_indx = find(label_i.count ~= 0 & label_i.count < val);
            if isempty(sub_indx)
                atlas_name = label_i.name{indx};
            end
        else
            atlas_name = label_i.name{indx};
        end
    end
    fprintf('Cluster %i) Anatomy Label: %s\n',cl_i,atlas_name);
    for des_i = 1:length(cond_iters)
        con_con = cond_iters{des_i};
        %##
    %     for g_i = 1:length(GROUP_NAMES)
    %         for cond_i = cond_iters
    %             tmp = squeeze(allersp_clusters(cl_i,cond_i,g_i,:,:,:));
    %             tmp(tmp == 0) = nan();
    %             tmp = squeeze(nanmean(tmp,1));
    %             allersp{cond_i,g_i} = tmp;
    %             tmp = squeeze(allgpm_clusters(cl_i,cond_i,g_i,:,:,:));
    %             tmp(tmp == 0) = nan();
    %             tmp = squeeze(nanmean(tmp,1));
    %             allgpm{cond_i,g_i} = tmp;
    %         end
    %     end
        allersp = cell(length(con_con),1); %cell(1,4);
        allgpm = cell(length(con_con),1); %cell(1,4);
        subjs = zeros(length(con_con),1);
        cnt = 1;
        for cond_i = con_con
%             tmp = reshape(squeeze(allersp_clusters(cl_i,cond_i,:,:,:,:)),size(allersp_clusters,4)*size(allersp_clusters,3),size(allersp_clusters,5),size(allersp_clusters,6));
            tmp = squeeze(allersp_clusters(cl_i,cond_i,:,:,:,:));
            tt = [];
            for i = 1:size(allersp_clusters,3)
                for j = 1:size(allersp_clusters,4)
                    tt = cat(3,tt,squeeze(tmp(i,j,:,:)));
                end
            end
            tt(tt == 0) = nan();
            inds = ~all(isnan(tt),[1,2]);
            tt = squeeze(tt(:,:,inds));
            subjs(cnt) = size(tt,3);
    %         tmp = squeeze(nanmean(tmp,1));
    %         allersp{cond_i} = tmp;
    %         allersp{cond_i} = reshape(tmp,size(tmp,2),size(tmp,3),size(tmp,1));
%             stack = zeros(size(tmp,2),size(tmp,3),size(tmp,1));
%             for subj_i = 1:size(tmp,1)
%                 stack(:,:,subj_i) = squeeze(tmp(subj_i,:,:));
%             end
            allersp{cnt} = tt;
%             tmp = reshape(squeeze(allgpm_clusters(cl_i,cond_i,:,:,:,:)),size(allgpm_clusters,4)*size(allgpm_clusters,3),size(allgpm_clusters,5),size(allgpm_clusters,6));
            tmp = squeeze(allgpm_clusters(cl_i,cond_i,:,:,:,:));
            tt = [];
            for i = 1:size(allgpm_clusters,3)
                for j = 1:size(allgpm_clusters,4)
                    tt = cat(3,tt,squeeze(tmp(i,j,:,:)));
                end
            end
            tt(tt == 0) = nan();
            inds = ~all(isnan(tt),[1,2]);
            tt = squeeze(tt(:,:,inds));
    %         tmp = squeeze(nanmean(tmp,1));
    %         allgpm{cond_i} = tmp;
    %         allgpm{cond_i} = reshape(tmp,size(tmp,2),size(tmp,3),size(tmp,1));
%             stack = zeros(size(tmp,2),size(tmp,3),size(tmp,1));
%             for subj_i = 1:size(tmp,1)
%                 stack(:,:,subj_i) = squeeze(tmp(subj_i,:,:));
%             end
            allgpm{cnt} = tt;
            cnt = cnt + 1;
        end
        %-
        for cond_i = 1:length(con_con)
            allersp{cond_i} = allersp{cond_i}(:,:,1:min(subjs));
            allgpm{cond_i} = allgpm{cond_i}(:,:,1:min(subjs));
        end
        %##
        allersp_com = cell(length(con_con),1);
        allersp_sb = cell(length(con_con),1);
        allgpm_com = cell(length(con_con),1);
        allgpm_sb  = cell(length(con_con),1);
        stackersp_com = [];
        stackersp_sb = [];
        stackgpm_com = [];
        stackgpm_sb = [];
        for cond_i = 1:length(con_con)
            stackersp_com = cat(3,stackersp_com,allersp{cond_i});
            stackgpm_com = cat(3,stackgpm_com,allgpm{cond_i});
            %-
            stackersp_sb = mean(allersp{cond_i},3);
            stackgpm_sb = mean(allgpm{cond_i},3);
            tmp1 = [];
            tmp2 = [];
            for subj_i = 1:size(allersp{cond_i},3)
               tmp1 = cat(3,tmp1,(allersp{cond_i}(:,:,subj_i) - stackersp_sb));
               tmp2 = cat(3,tmp2,(allgpm{cond_i}(:,:,subj_i) - stackgpm_sb));
            end
            allersp_sb{cond_i} = tmp1;
            allgpm_sb{cond_i} = tmp2;
    %         stackersp_com = cat(3,stackersp_com,allersp_sb{cond_i});
    %         stackgpm_com = cat(3,stackgpm_com,allgpm_sb{cond_i});
        end
        stackersp_com = mean(stackersp_com,3);
        stackgpm_com = mean(stackgpm_com,3);
        for cond_i = 1:length(con_con)
    %         allersp_com{cond_i} = allersp_sb{cond_i} - stackersp_com;
    %         allgpm_com{cond_i} = allgpm_sb{cond_i} - stackgpm_com;
            tmp1 = [];
            tmp2 = [];
            for subj_i = 1:size(allersp{cond_i},3)
                tmp1 = cat(3,tmp1,allersp{cond_i}(:,:,subj_i) - stackersp_com);
                tmp2 = cat(3,tmp2,allgpm{cond_i}(:,:,subj_i) - stackgpm_com);
            end
            allersp_com{cond_i} = tmp1;
            allgpm_com{cond_i} = tmp2;
        end
        %##
        allfreqs = 4:100;
        alltimes = 1:100;
        alltitles = std_figtitle('condnames',condition_gait(cond_iters{des_i}), ...
            'clustname', sprintf('%s',atlas_name)); %sprintf('CL%i',cl_i));
        %% NO BASELINE
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm,allfreqs,alltimes);
        %- subject specific plots
        %{
        for subj_i = 1:size(allersp{1},3)
            %- ERSP
            tmp_ersp = cell(4,1);
            for cond_i = 1:4
                tmp_ersp{cond_i} = squeeze(allersp{cond_i}(:,:,subj_i));
            end
            plot_txf_conds_tftopo(tmp_ersp,alltimes,allfreqs,alltitles,[],...
            2,linspecer,[4,100],[100,100,350,350])
            %- GPM
            tmp_gpm = cell(4,1);
            for cond_i = 1:4
                tmp_ersp{cond_i} = squeeze(allgpm{cond_i}(:,:,subj_i));
            end
            plot_txf_conds_tftopo(tmp_ersp,alltimes,allfreqs,alltitles,[],...
            2,linspecer,[4,100],[100,100,350,350])
        end
        %}
        CLIM = 1.5;
        %- calculate per condition means
        for cond_i = 1:length(con_con)
            allersp{cond_i} = squeeze(mean(allersp{cond_i},3));
            allgpm{cond_i} = squeeze(mean(allgpm{cond_i},3));
        end
        %##
        fig = plot_txf_conds_tftopo(allersp,alltimes,allfreqs,alltitles,pcond_ersp_crop{1},...
            CLIM,linspecer,[4,100],[100,100,350,350]);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_ersp.jpg',cl_i,des_i)]);
        %##
        fig = plot_txf_conds_tftopo(allgpm,alltimes,allfreqs,alltitles,pcond_gpm_crop{1},...
            CLIM,linspecer,[4,100],[100,100,350,350]);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_gpm.jpg',cl_i,des_i)]);
        %% WITHIN CONDITION BASELINED ERSPs
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp_sb,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm_sb,allfreqs,alltimes);
        %- subject specific plots
        %{
        for subj_i = 1:size(allersp_sb{1},3)
            tmp_ersp = cell(4,1);
            for cond_i = 1:4
                tmp_ersp{cond_i} = squeeze(allersp_sb{cond_i}(:,:,subj_i));
            end
            plot_txf_conds_tftopo(tmp_ersp,alltimes,allfreqs,alltitles,[],...
            2,linspecer,[4,100],[100,100,350,350])
        end
        %}
        %-
        CLIM = 1.5;
        FIG_POS = [100,100,350,350];
        FREQ_LIMS = [4,100];
        for cond_i = 1:length(con_con)
            allersp_sb{cond_i} = squeeze(mean(allersp_sb{cond_i},3));
            allgpm_sb{cond_i} = squeeze(mean(allgpm_sb{cond_i},3));
        end
        %##
        fig = plot_txf_conds_tftopo(allersp_sb,alltimes,allfreqs,alltitles,pcond_ersp_crop{1},...
            CLIM,linspecer,FREQ_LIMS,FIG_POS);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_ersp_sb.jpg',cl_i,des_i)]);
        %##
        fig = plot_txf_conds_tftopo(allgpm_sb,alltimes,allfreqs,alltitles,pcond_gpm_crop{1},...
            CLIM,linspecer,FREQ_LIMS,FIG_POS);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_gpm_sb.jpg',cl_i,des_i)]);
        %% ACROSS CONDITIONS BASELINED ERSPS
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp_com,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm_com,allfreqs,alltimes);
        %- subject specific plots
        %{
        for subj_i = 1:size(allersp_com{1},3)
            tmp_ersp = cell(4,1);
            for cond_i = 1:4
                tmp_ersp{cond_i} = squeeze(allersp_com{cond_i}(:,:,subj_i));
            end
            plot_txf_conds_tftopo(tmp_ersp,alltimes,allfreqs,alltitles,[],...
            2,linspecer,[4,100],[100,100,350,350])
        end
        %}
        %-
        CLIM = 1.5;
        FIG_POS = [100,100,350,350];
        FREQ_LIMS = [4,100];
        for cond_i = 1:length(con_con)
            allersp_com{cond_i} = squeeze(mean(allersp_com{cond_i},3));
            allgpm_com{cond_i} = squeeze(mean(allgpm_com{cond_i},3));
        end
        %##
        fig = plot_txf_conds_tftopo(allersp_com,alltimes,allfreqs,alltitles,pcond_ersp_crop{1},...
            CLIM,linspecer,FREQ_LIMS,FIG_POS);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_ersp_com.jpg',cl_i,des_i)]);
        %##
        fig = plot_txf_conds_tftopo(allgpm_com,alltimes,allfreqs,alltitles,pcond_gpm_crop{1},...
            CLIM,linspecer,FREQ_LIMS,FIG_POS);
        exportgraphics(fig,[cluster_dir filesep sprintf('cl%i_des%i_spca_gpm_com.jpg',cl_i,des_i)]);
    end
end

%##
% fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
% exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
%##
