%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh

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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%% Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
%- statistics & conditions
SPEED_VALS = {'0.25','0.5','0.75','1.0';
              '0p25','0p5','0p75','1p0'};
TERRAIN_VALS = {'flat','low','med','high'};
COLORS_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
%## hard define
%- FOOOF
settings = struct('peak_width_limits',[1,8],...
    'min_peak_height',0.05,...
    'max_n_peaks',3);
f_range = [3, 40];
theta_band = [4, 8];
alpha_band = [8 12];
beta_band  = [12 30];
%- datset name
DATA_SET = 'MIM_dataset';
%- cluster directory for study
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';

%- study info
SUB_GROUP_FNAME = 'group_spec';
% SUB_GROUP_FNAME = 'all_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## SET STUDY PATHS
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDY
%{
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
%}
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',spec_data_dir);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',spec_data_dir);
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
%% RE-POP PARAMS
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
    'groupstats',ERSP_STAT_PARAMS.groupstats,...
    'method',ERSP_STAT_PARAMS.method,...
    'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
    'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
    'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
    'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
%% (LOAD EXISTING TALBES && FORMAT STUDY)
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
tmp = load([save_dir filesep 'psd_band_power_stats.mat']);
psd_feature_stats = tmp.psd_feature_stats;
tmp = load([save_dir filesep 'fooof_pcond.mat']);
pcond = tmp.pcond;
tmp = load([save_dir filesep 'org_pcond.mat']);
pcond_org = tmp.pcond_org;
% tmp = load([save_dir filesep 'fooof_results_summary.mat']);
% fooof_group_results_org = tmp.fooof_group_results_org;
tmp = load([save_dir filesep 'fooof_diff_store.mat']);
fooof_diff_store = tmp.fooof_diff_store;
tmp = load([save_dir filesep 'fooof_apfit_store.mat']);
fooof_apfit_store = tmp.fooof_apfit_store;
tmp = load([save_dir filesep 'spec_data_original.mat']);
spec_data_original = tmp.spec_data_original;
tmp = load([save_dir filesep 'fooof_results.mat']);
fooof_results = tmp.fooof_results;
fooof_freq = fooof_results{1}{1,1}{1}.freqs;
%## TOPO PLOTS
% tmp_study = STUDY;
% RE_CALC = true;
% if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') || isfield(tmp_study.cluster,'topopol') 
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
% end
% if ~isfield(tmp_study.cluster,'topo'), tmp_study.cluster(1).topo = [];end
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% for i = 1:length(designs)
%     des_i = string(designs(i));
%     for j = 1:length(clusters) % For each cluster requested
%         cl_i = double(string(clusters(j)));
%         if isempty(tmp_study.cluster(cl_i).topo) || RE_CALC
%             inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == string(cl_i));
%             sets_i = unique([FOOOF_TABLE.subj_cl_ind(inds)]);
%             tmp_study.cluster(cl_i).sets = tmp_study.cluster(cl_i).sets(sets_i);
%             tmp_study.cluster(cl_i).comps = tmp_study.cluster(cl_i).comps(sets_i);
%             tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
%             STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
%             STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
%             STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
%             STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
%             STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
%         end
%     end
% end
tmp_study = STUDY;
RE_CALC = true;
if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') || isfield(tmp_study.cluster,'topopol') 
    tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
    tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
    tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
    tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
    tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
end
if ~isfield(tmp_study.cluster,'topo'), tmp_study.cluster(1).topo = [];end
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
for j = 1:length(CLUSTER_PICKS) % For each cluster requested
    cluster_i = CLUSTER_PICKS(j);
    if isempty(tmp_study.cluster(cluster_i).topo) || RE_CALC
        tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cluster_i);% Using this custom modified code to allow taking average within participant for each cluster
        STUDY.cluster(cluster_i).topox = tmp_study.cluster(cluster_i).topox;
        STUDY.cluster(cluster_i).topoy = tmp_study.cluster(cluster_i).topoy;
        STUDY.cluster(cluster_i).topoall = tmp_study.cluster(cluster_i).topoall;
        STUDY.cluster(cluster_i).topo = tmp_study.cluster(cluster_i).topo;
        STUDY.cluster(cluster_i).topopol = tmp_study.cluster(cluster_i).topopol;
    end
end
%## STATS
iter = 200; % in eeglab, the fdr stats will automatically *20
try
    STUDY.etc = rmfield(STUDY.etc,'statistics');
end
% STUDY = pop_statparams(STUDY,'groupstats','on','condstats','on','statistics','perm',...
%     'singletrials','off','mode','eeglab','effect','main','alpha',NaN,'mcorrect','fdr','naccu',iter);% If not using mcorrect, use none, Not sure why, if using fdr correction, none of these are significant
% 
STUDY = pop_statparams(STUDY, 'groupstats','off','condstats', 'on',...
            'method','perm',...
            'singletrials','off','mode','fieldtrip','fieldtripalpha',NaN,...
            'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',iter*20);

stats = STUDY.etc.statistics;
stats.paired{1} = 'on'; % Condition stats
stats.paired{2} = 'off'; % Group stats
%% ===================================================================== %%
%## PARAMS
%-
ATLAS_PATH = [PATHS.submods_dir,...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%-
SUBPLOT_BOTTOM = 0.2;
SUBPLOT_INIT_SHIFT = 0.05;
PG_SIZE = [6.5,9];
PLOT_STRUCT = struct('figure_position_inch',[3,3,6.5,9],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[],...
    'time_lims',[],...
    'subplot_width',0.15,...
    'subplot_height',0.65,...
    'shift_amnt',0.175,...
    'stats_title','F Statistic Mask',...
    'figure_title','');
%% ===================================================================== %%
clusters = unique(FOOOF_TABLE.cluster_id);
% [STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
txt_store = cell(length(clusters),1);
atlas_name_store = cell(length(clusters),1);
for k_i = 1:length(clusters)
    k = double(string(clusters(k_i)));
    %## ANATOMY
    dip1 = STUDY.cluster(k).all_diplocs;
    STUDY.cluster(k).centroid.dipole.posxyz = mean(dip1);
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'single';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    cfg.verbose = 0;
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if ~isempty(sub_indx)
                atlas_name = names{sub_indx};
            end
        else
            atlas_name = names{indx};
        end
    end
    %## ANATOMY
    
    dip1 = STUDY.cluster(k).centroid.dipole.posxyz;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name_ct = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'multiple';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    cfg.verbose = 0;
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if ~isempty(sub_indx)
                atlas_name_ct = names{sub_indx};
            end
        else
            atlas_name_ct = names{indx};
        end
    end
    txt_store{k} = [sprintf('CL%i: N=%i\n',k,length(STUDY.cluster(k).sets)),...
    sprintf('CL%i: %s\n',k,atlas_name),...
    sprintf('Dip Center: [%0.1f,%0.1f,%0.1f]\n',STUDY.cluster(k).dipole.posxyz),...
    sprintf('CENTROID: CL%i: %s\n',k,atlas_name_ct),...
    sprintf('CENTROID: Dip %0.1f,%0.1f,%0.1f]\n\n',STUDY.cluster(k).centroid.dipole.posxyz)];
    atlas_name_store{k_i} = sprintf('CL%i: %s\n',k,atlas_name);
    % atlas_name_store{k} = sprintf('CL%i: %s\n',k,atlas_name_ct);
end
cellfun(@(x) disp(x),txt_store);
%% ==================================================================== %%
%## MIM KINEMATICS
% meas_names_imu = {'nanmean_APexc_mean','nanmean_MLexc_mean','nanmean_APexc_COV','nanmean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
% meas_names_ls = {'nanmean_StepDur','nanmean_StepDur_cov','nanmean_GaitCycleDur_cov','nanmean_GaitCycleDur',};
%##
SCATTER_BOTTOM = 0.65;
im_resize = 0.7;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
subjects = unique(FOOOF_TABLE.subj_id);
design_chars = {'terrain','speed'};
group_chars = unique(FOOOF_TABLE.group_char);
conditions = unique(FOOOF_TABLE.cond_char);
PLOT_PARAMS = struct('color_map',linspecer(4),...
                'cond_labels',unique(FOOOF_TABLE.cond_char),'group_labels',unique(FOOOF_TABLE.group_char),...
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label','10*log_{10}(Flattened PSD)',...
                'title','','font_size',9,'y_lim',[-1,15],...
                'font_name','Arial','x_label','');
measure_name_plot = {'theta_avg_power','alpha_avg_power','beta_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \alpha','Mean \beta'};
% measure_name_plot = {'med_sub_flat','low_sub_flat','high_sub_flat'};
%% ================================================================= %%
% % measure_name_plot = {'theta_sub','alpha_sub','beta_sub'};
% T_FOOOF_TABLE = FOOOF_TABLE;
% T_FOOOF_TABLE.med_sub_flat = zeros(size(T_FOOOF_TABLE,1),1);
% T_FOOOF_TABLE.low_sub_flat = zeros(size(T_FOOOF_TABLE,1),1);
% T_FOOOF_TABLE.high_sub_flat = zeros(size(T_FOOOF_TABLE,1),1);
% 
% for des_i = 1
%     for subj_i = 1:length(subjects)
%         inds = T_FOOOF_TABLE.design_id == num2str(des_i) & T_FOOOF_TABLE.subj_id == string(subj_i);
%         T_plot = T_FOOOF_TABLE(inds,:);
%         % inds = psd_feature_stats.study 
%         flat_pwr = T_plot.cond_id == string(1);
%         low_pwr = T_plot.cond_id == string(2);
%         med_pwr = T_plot.cond_id == string(3);
%         high_pwr = T_plot.cond_id == string(4);
%         for meas_i = 1:length(measure_name_plot)
%             T_plot.med_sub_flat(med_pwr) = T_plot.(measure_name_plot{meas_i})(med_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%             T_plot.high_sub_flat(high_pwr) = T_plot.(measure_name_plot{meas_i})(high_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%             T_plot.low_sub_flat(low_pwr) = T_plot.(measure_name_plot{meas_i})(low_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%             %-
%             % T_plot.med_sub_flat(med_pwr) = T_plot.(measure_name_plot{meas_i})(med_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%             % T_plot.high_sub_flat(high_pwr) = T_plot.(measure_name_plot{meas_i})(high_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%             % T_plot.low_sub_flat(low_pwr) = T_plot.(measure_name_plot{meas_i})(low_pwr)-T_plot.(measure_name_plot{meas_i})(flat_pwr);
%         end
%         T_FOOOF_TABLE(inds,:) = T_plot;
%     end
% end

%% ===================================================================== %%
%## TOPO & DIPOLE PLOTS
for k_i = 1:length(clusters)
    cl_i = double(string(clusters(k_i)));
    %## ANATOMY
    atlas_name = atlas_name_store{k_i};
    %## AXES LIMITS
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    %## ALIGNMENT
    GAP = 0.05;
    LEFT_D = 2;
    BOT_D = 6.75;
    %## topo plot 
    im_resize = 5;
    subplot(3,4,1)
    std_topoplot_CL(STUDY,cl_i,'together');
    colormap(linspecer); %colormap_ersp)
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'color','w')
    fig_i.Children(1).Title.String = sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    fig_i.Children(1).Title.Interpreter = 'none';
    fig_i.Children(1).FontSize = 12; %PLOT_STRUCT.font_size;
    fig_i.Children(1).FontName = PLOT_STRUCT.font_name;
    fig_i.Children(1).FontWeight = 'bold';
    fig_i.Children(1).OuterPosition = [0,0,1,1];
    fig_i.Children(1).Units = 'Inches';
    fig_i.Children(1).Position = [0.5,BOT_D-0.175,0.225*im_resize,0.25*im_resize];  %[left bottom width height]
    %## dipole plots (terrain)
    %-
    dpi = 1000;
    target_sz = 1.25; %inch
    target_dim = 1; %2==width, 1==height
    im1 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cl_i)]);
    
    if target_dim==1
        scale = target_sz/(size(im1,1)/dpi);
        im1 = imresize(im1,scale);
    elseif target_dim==2
        scale = target_sz/(size(im1,2)/dpi);
        im1 = imresize(im1,scale);
    end
    
    im2 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cl_i)]);
    im2(1:300,:,:) = [];
    im2(end-300:end,:,:)=[];
    if target_dim==1
        scale = target_sz/(size(im2,1)/dpi);
        im2 = imresize(im2,scale);
    elseif target_dim==2
        scale = target_sz/(size(im2,2)/dpi);
        im2 = imresize(im2,scale);
    end

    im3 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cl_i)]);
%         im3(1,:,:) = [];
    if target_dim==1
        scale = target_sz/(size(im3,1)/dpi);
        im3 = imresize(im3,scale);
    elseif target_dim==2
        scale = target_sz/(size(im3,2)/dpi);
        im3 = imresize(im3,scale);
    end
    
    %-
    szw1 = (size(im1,2)/dpi); %/PG_SIZE(1); %width
    szh1 = (size(im1,1)/dpi); %/PG_SIZE(2); %+0.0001; %height
    szw2 = (size(im2,2)/dpi); %/PG_SIZE(1); %width
    szh2 = (size(im2,1)/dpi); %+0.05; %/PG_SIZE(2); %+0.01;
    szw3 = (size(im3,2)/dpi); %/PG_SIZE(1);
    szh3 = (size(im3,1)/dpi); %/PG_SIZE(2);
    szw = max([szw1,szw2,szw3]);
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D,BOT_D,szw,szh1],'PositionConstraint','outerposition');
    imshow(im1,'border','tight');
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+((szw-szw1)-(szw-szw2)),BOT_D,szw,szh2],'PositionConstraint','outerposition');
    imshow(im2,'border','tight');
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+szw2+0.05,BOT_D,szw,szh3],'PositionConstraint','outerposition');
    imshow(im3,'border','tight');
    %##
    % exportgraphics(fig,[save_dir filesep sprintf('TOPO_DIP_cl%i.tiff',cl_i)],'Resolution',1000)
    exportgraphics(fig,[save_dir filesep sprintf('TOPO_DIP_cl%i.tiff',cl_i)],'Resolution',300)
    close(fig);
end
%% ================================================================= %%
%## PSDS
PSD_BOTTOM = 0.7;
im_resize = 0.5;
GROUP_EDGECOLOR = {};
GROUP_LINESTYLE = {};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
for k_i = 1:length(clusters)
    %##
    cl_i = double(string(clusters(k_i)));
    atlas_name = atlas_name_store{k_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    for j = 1:length(groups)
        %## PLOT SPEEED PSDS
        horiz_shift = 0;
        for des_i = 1:length(designs)
            switch des_i
                case 1
                    color_dark = COLORS_MAPS_TERRAIN;
                    color_light = COLORS_MAPS_TERRAIN;
                    GROUP_CMAP_OFFSET = [0,0.1,0.1];
                    xtick_label_g = {'flat','low','med','high'};
                case 2
                    color_dark = COLOR_MAPS_SPEED;
                    color_light = COLOR_MAPS_SPEED+0.15;
                    GROUP_CMAP_OFFSET = [0.15,0,0];
                    xtick_label_g = {'0.25','0.50','0.75','1.0'};
            end
            %## non-fooof psd (speed)
            axes();
            hold on;
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));
                Pa.EdgeColor = "none";
            end
            axs = [];
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                axs = [axs, ax];
            end 
            %- plot the aperiodic line
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                aperiodic_fit = fooof_apfit_store{des_i}{cl_i}{i,j}';
                dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
            end
            ax = gca;
            xlim([4 40]);
            ylim([-30 -5]);
            [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond_org{des_i}{cl_i}{1}(:,2), 'background', 'Frequency(Hz)');
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            title('PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.08+horiz_shift,PSD_BOTTOM-vert_shift,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
            hold off;
        
            %## fooof psd (speed)
            axes();
            hold on;
            axs = [];
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));
                Pa.EdgeColor = "none";
            end
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                axs = [axs, ax];
            end
            %-
            ax = gca;
            [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond{des_i}{cl_i}{j}(:,2), 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
        %     set(ax,'LineWidth',2)
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            %- legend
            legend([axs,dash],'FontSize',9,'FontName',PLOT_STRUCT.font_name);
            [lg1,icons,plots,txt] = legend('boxoff');
            set(lg1,'Position',[0.20+0.3*im_resize+horiz_shift,PSD_BOTTOM+0.025-vert_shift,0.2,0.1]);
            lg1.ItemTokenSize(1) = 18;
            %-
            title('Flattened PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            carry_ov = 0.12+0.3*im_resize;
            set(ax,'Position',[carry_ov+horiz_shift,PSD_BOTTOM-vert_shift,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
        %         icons(2).XData = [0.05 0.1];
            horiz_shift = horiz_shift + carry_ov + 0.25*im_resize + 0.1;
        end
        %## TITLE
        annotation('textbox',[0.5-0.1,PSD_BOTTOM-vert_shift-0.05+0.25*im_resize,0.2,0.2],...
            'String',string(group_chars(j)),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            'FontSize',14,'FontWeight','Bold','Units','normalized');
        % close(fig);
        vert_shift = vert_shift + 0.25*im_resize+0.1;
    end
    hold off;
    % exportgraphics(fig,[save_dir filesep sprintf('Group_Violins_cl%i.tiff',cl_i)],'Resolution',1000)
    exportgraphics(fig,[save_dir filesep sprintf('Group_PSDs_cl%i.tiff',cl_i)],'Resolution',300)
    close(fig);
end
%% ================================================================= %%
% FOOOF_TABLE = T_FOOOF_TABLE;
%## VIOLIN PLOTS
im_resize= 0.9;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
for k_i = 1:length(clusters)
    %
    atlas_name = ''; %atlas_name_store{k_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    cl_i = double(string(clusters(k_i)));
    %## violin plot's theta/alpha/beta (speed)
    %-
    max_val = zeros(length(measure_name_plot),length(designs));
    % STATS_STRUCT = struct('anova',{{}},...
    %                       'anova_grp',{{}},...
    %                       'pvals',{{}},...
    %                       'pvals_pairs',{{}},...
    %                       'pvals_grp',{{}},...
    %                       'pvals_grp_pairs',{{}},...
    %                       'regress_pval',{{}},...
    %                       'regress_line',{{}},...
    %                       'r2_coeff',{{}},...
    %                       'regress_xvals',0);
    %-
    DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{[]},...
                          'regress_xvals',0,...
                          'subject_char',[],... % this option when filled prints removal of nan() info
                          'group_order',categorical({''}),...
                          'do_include_intercept',false); 
    %
    %##
    STATS_STRUCT = DEFAULT_STATS_STRUCT;
    cnt = 1;
    for des_i = 1:length(designs)
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
            T_plot = FOOOF_TABLE(inds,:);
            STATS_STRUCT(cnt) = DEFAULT_STATS_STRUCT;
            for k = 1:length(groups)
                inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(k);
                T_stats_plot = psd_feature_stats(inds,:);
                switch i 
                    case 1
                        aa = T_stats_plot.theta_anova;
                        c2s = T_stats_plot.theta_cond2_pval;
                        c3s = T_stats_plot.theta_cond3_pval;
                        c4s = T_stats_plot.theta_cond4_pval;
                        rs = T_stats_plot.Th_num_pval;
                        rls = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                        r2 = T_stats_plot.Th_num_R2;
                        norm_p = T_stats_plot.theta_lilnorm_p;
                    case 2
                        aa = T_stats_plot.alpha_anova;
                        c2s = T_stats_plot.alpha_cond2_pval;
                        c3s = T_stats_plot.alpha_cond3_pval;
                        c4s = T_stats_plot.alpha_cond4_pval;
                        rs = T_stats_plot.A_num_pval;
                        rls = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                        r2 = T_stats_plot.A_num_R2;
                        norm_p = T_stats_plot.alpha_lilnorm_p;
                    case 3
                        aa = T_stats_plot.beta_anova;
                        c2s = T_stats_plot.beta_cond2_pval;
                        c3s = T_stats_plot.beta_cond3_pval;
                        c4s = T_stats_plot.beta_cond4_pval;
                        rs = T_stats_plot.B_num_pval;
                        rls = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                        r2 = T_stats_plot.B_num_R2;
                        norm_p = T_stats_plot.beta_lilnorm_p;
                end
                if des_i == 1
                    STATS_STRUCT(cnt).anova{k}=aa;
                    STATS_STRUCT(cnt).pvals{k}=[1,c2s,c3s,c4s];
                    STATS_STRUCT(cnt).pvals_pairs{k}={[1,1],[1,2],[1,3],[1,4]};
                end
                if des_i == 2
                    STATS_STRUCT(cnt).anova{k}=aa;
                    STATS_STRUCT(cnt).regress_pval{k}=rs;
                    STATS_STRUCT(cnt).regress_line{k}=rls;
                    STATS_STRUCT(cnt).r2_coeff(k)=r2;
                    STATS_STRUCT(cnt).regress_xvals=(0:5)*0.25;
                end
            end
            stat_add = (max(T_plot.(measure_name))-min(T_plot.(measure_name)))*0.2;
            for k = 1:length(groups)
                if des_i == 1
                    % tm = max(T_plot.(measure_name))+sum([STATS_STRUCT(cnt).pvals{k}(:)]<0.05)*stat_add;
                    tm = max(T_plot.(measure_name))+2*std(T_plot.(measure_name));
                end
                if des_i == 2
                    % tm = max(T_plot.(measure_name))+sum([STATS_STRUCT(cnt).regress_pval{k}]<0.05)*(stat_add*2);
                    tm = max(T_plot.(measure_name))+2*std(T_plot.(measure_name));
                end
            end
            max_val(i,des_i) = tm;
            cnt = cnt + 1;
        end
    end
    max_vals = max(max_val,[],2)+2*std(max_val,[],2);
    %-
    vert_shift = 0;
    cnt = 1;
    for des_i = 1:length(designs)   
        inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
        T_plot = FOOOF_TABLE(inds,:);
        switch des_i
            case 1
                
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
        end
        horiz_shift= 0;
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            tmp_stats = STATS_STRUCT(cnt);
            % figure;
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_PARAMS = struct('color_map',color_dark,...
                'cond_labels',unique(T_plot.cond_char),'group_labels',unique(T_plot.group_char),...
                'cond_offsets',cond_offsets,...
                'group_offsets',[0.125,0.475,0.812],...
                'y_label','10*log_{10}(Flattened PSD)',...
                'title',title_plot{i},'font_size',7,'ylim',[min(T_plot.(measure_name))-0.3,max_vals(i)],...
                'font_name','Arial','x_label',x_label,'do_combine_groups',false);
            % ax = axes();
            % figfig = figure();
            axax = group_violin(T_plot,measure_name,'cond_id','group_id',...
                fig,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_PARAMS,...
                'STATS_STRUCT',tmp_stats);
            set(axax,'OuterPosition',[0,0,1,1]);
            set(axax,'Position',[0.08+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
            hold off;
            %- iterate
            horiz_shift = horiz_shift + 0.3*im_resize+0.05;
            cnt = cnt + 1;
        end
        %## TITLE
        annotation('textbox',[0.5-0.1,VIOLIN_BOTTOM+vert_shift-0.075+0.225*im_resize,0.2,0.2],...
            'String',string(design_chars{des_i}),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            'FontSize',14,'FontWeight','Bold','Units','normalized');
        vert_shift = vert_shift - (0.1+0.225*im_resize);
    end
    hold off;
    % exportgraphics(fig,[save_dir filesep sprintf('Group_Violins_cl%i.tiff',cl_i)],'Resolution',1000)
    exportgraphics(fig,[save_dir filesep sprintf('Group_Violins_cl%i.tiff',cl_i)],'Resolution',300)
    close(fig);
end
%% ===================================================================== %%
%## COMBINED PLOT FOR ALL (NO GROUPING)
CLUSTERS_TO_PLOT = main_cl_inds(2:end); %valid_clusters(1:end-1); %main_cl_inds(2:end);
VIOLIN_BOTTOM = 0.375;
PSD_BOTTOM = 0.575;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
AX_H  = 0.2;
AX_W = 0.25;
for k_i = 1:length(CLUSTERS_TO_PLOT)
    atlas_name = atlas_name_store{k_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    cl_i = double(string(clusters(k_i)));
    %## ALIGNMENT
    GAP = 0.05;
    LEFT_D = 2;
    BOT_D = 6.75;
    %## topo plot 
    im_resize = 5;
    subplot(3,4,1)
    std_topoplot_CL(STUDY,cl_i,'together');
    colormap(linspecer); %colormap_ersp)
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'color','w')
    fig_i.Children(1).Title.String = sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    fig_i.Children(1).Title.Interpreter = 'none';
    fig_i.Children(1).FontSize = 12; %PLOT_STRUCT.font_size;
    fig_i.Children(1).FontName = PLOT_STRUCT.font_name;
    fig_i.Children(1).FontWeight = 'bold';
    fig_i.Children(1).OuterPosition = [0,0,1,1];
    fig_i.Children(1).Units = 'Inches';
    fig_i.Children(1).Position = [0.5,BOT_D-0.175,0.225*im_resize,0.25*im_resize];  %[left bottom width height]
    %## dipole plots (terrain)
    %-
    dpi = 1000;
    target_sz = 1.25; %inch
    target_dim = 1; %2==width, 1==height
    im1 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cl_i)]);
    
    if target_dim==1
        scale = target_sz/(size(im1,1)/dpi);
        im1 = imresize(im1,scale);
    elseif target_dim==2
        scale = target_sz/(size(im1,2)/dpi);
        im1 = imresize(im1,scale);
    end
    
    im2 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cl_i)]);
    im2(1:300,:,:) = [];
    im2(end-300:end,:,:)=[];
    if target_dim==1
        scale = target_sz/(size(im2,1)/dpi);
        im2 = imresize(im2,scale);
    elseif target_dim==2
        scale = target_sz/(size(im2,2)/dpi);
        im2 = imresize(im2,scale);
    end

    im3 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cl_i)]);
%         im3(1,:,:) = [];
    if target_dim==1
        scale = target_sz/(size(im3,1)/dpi);
        im3 = imresize(im3,scale);
    elseif target_dim==2
        scale = target_sz/(size(im3,2)/dpi);
        im3 = imresize(im3,scale);
    end
    
    %-
    szw1 = (size(im1,2)/dpi); %/PG_SIZE(1); %width
    szh1 = (size(im1,1)/dpi); %/PG_SIZE(2); %+0.0001; %height
    szw2 = (size(im2,2)/dpi); %/PG_SIZE(1); %width
    szh2 = (size(im2,1)/dpi); %+0.05; %/PG_SIZE(2); %+0.01;
    szw3 = (size(im3,2)/dpi); %/PG_SIZE(1);
    szh3 = (size(im3,1)/dpi); %/PG_SIZE(2);
    szw = max([szw1,szw2,szw3]);
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D,BOT_D,szw,szh1],'PositionConstraint','outerposition');
    imshow(im1,'border','tight');
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+((szw-szw1)-(szw-szw2)),BOT_D,szw,szh2],'PositionConstraint','outerposition');
    imshow(im2,'border','tight');
    %-
    axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+szw2+0.05,BOT_D,szw,szh3],'PositionConstraint','outerposition');
    imshow(im3,'border','tight');
    %## PSDS
    im_resize = 0.5;
    vert_shift = 0;
    if length(groups)>1
        error('This is for aggregate stats only. Use plots above for multigroup plotting');
    end
    for j = 1:length(groups)
        %## PLOT SPEEED PSDS
        horiz_shift = 0;
        for des_i = 1:length(designs)
            switch des_i
                case 1
                    color_dark = COLORS_MAPS_TERRAIN;
                    color_light = COLORS_MAPS_TERRAIN;
                    GROUP_CMAP_OFFSET = [0,0.1,0.1];
                    xtick_label_g = {'flat','low','med','high'};
                case 2
                    color_dark = COLOR_MAPS_SPEED;
                    color_light = COLOR_MAPS_SPEED+0.15;
                    GROUP_CMAP_OFFSET = [0.15,0,0];
                    xtick_label_g = {'0.25','0.50','0.75','1.0'};
            end
            %## non-fooof psd (speed)
            axes();
            hold on;
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:),color_light(i,:));
                Pa.EdgeColor = "none";
            end
            axs = [];
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                axs = [axs, ax];
            end 
            %- plot the aperiodic line
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                aperiodic_fit = fooof_apfit_store{des_i}{cl_i}{i,j}';
                dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
            end
            ax = gca;
            xlim([4 40]);
            ylim([-30 -5]);
            [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond_org{des_i}{cl_i}{1}(:,2), 'background', 'Frequency(Hz)');
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            title('PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.08+horiz_shift,PSD_BOTTOM-vert_shift,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
            hold off;
        
            %## fooof psd (speed)
            axes();
            hold on;
            axs = [];
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:),color_light(i,:));
                Pa.EdgeColor = "none";
            end
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                axs = [axs, ax];
            end
            %-
            ax = gca;
            [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond{des_i}{cl_i}{j}(:,2), 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
        %     set(ax,'LineWidth',2)
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            %- legend
            legend([axs,dash],'FontSize',9,'FontName',PLOT_STRUCT.font_name);
            [lg1,icons,plots,txt] = legend('boxoff');
            set(lg1,'Position',[0.20+0.3*im_resize+horiz_shift,PSD_BOTTOM+0.025-vert_shift,0.2,0.1]);
            lg1.ItemTokenSize(1) = 18;
            %-
            title('Flattened PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            carry_ov = 0.12+0.3*im_resize;
            set(ax,'Position',[carry_ov+horiz_shift,PSD_BOTTOM-vert_shift,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
        %         icons(2).XData = [0.05 0.1];
            horiz_shift = horiz_shift + carry_ov + 0.25*im_resize + 0.1;
        end
        %## TITLE
        % annotation('textbox',[0.5-0.1,PSD_BOTTOM-vert_shift-0.05+0.25*im_resize,0.2,0.2],...
        %     'String',string(group_chars(j)),'HorizontalAlignment','center',...
        %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
        %     'FontSize',14,'FontWeight','Bold','Units','normalized');
        % close(fig);
        vert_shift = vert_shift + 0.25*im_resize+0.1;
    end
    hold off;
    %## violin plot's theta/alpha/beta (speed)
    im_resize = 0.5;
    max_val = zeros(length(measure_name_plot),length(designs));
    STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{{}},...
                          'regress_xvals',0);
    %##
    cnt = 1;
    for des_i = 1:length(designs)
        
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
            T_plot = FOOOF_TABLE(inds,:);
            
            for k = 1:length(groups)
                inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(k);
                T_stats_plot = psd_feature_stats(inds,:);
                switch i 
                    case 1
                        aa = T_stats_plot.theta_anova;
                        c2s = T_stats_plot.theta_cond2_pval;
                        c3s = T_stats_plot.theta_cond3_pval;
                        c4s = T_stats_plot.theta_cond4_pval;
                        rs = T_stats_plot.Th_num_pval;
                        rls = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                        r2 = T_stats_plot.Th_num_R2;
                        norm_p = T_stats_plot.theta_lilnorm_p;
                    case 2
                        aa = T_stats_plot.alpha_anova;
                        c2s = T_stats_plot.alpha_cond2_pval;
                        c3s = T_stats_plot.alpha_cond3_pval;
                        c4s = T_stats_plot.alpha_cond4_pval;
                        rs = T_stats_plot.A_num_pval;
                        rls = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                        r2 = T_stats_plot.A_num_R2;
                        norm_p = T_stats_plot.alpha_lilnorm_p;
                    case 3
                        aa = T_stats_plot.beta_anova;
                        c2s = T_stats_plot.beta_cond2_pval;
                        c3s = T_stats_plot.beta_cond3_pval;
                        c4s = T_stats_plot.beta_cond4_pval;
                        rs = T_stats_plot.B_num_pval;
                        rls = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                        r2 = T_stats_plot.B_num_R2;
                        norm_p = T_stats_plot.beta_lilnorm_p;
                end
                if des_i == 1
                    STATS_STRUCT(cnt).anova{k}=aa;
                    STATS_STRUCT(cnt).pvals{k}=[1,c2s,c3s,c4s];
                    STATS_STRUCT(cnt).pvals_pairs{k}={[1,1],[1,2],[1,3],[1,4]};
                end
                if des_i == 2
                    STATS_STRUCT(cnt).anova{k}=aa;
                    STATS_STRUCT(cnt).regress_pval{k}=rs;
                    STATS_STRUCT(cnt).regress_line{k}=rls;
                    STATS_STRUCT(cnt).r2_coeff(k)=r2;
                    STATS_STRUCT(cnt).regress_xvals=(0:5)*0.25;
                end
            end
            stat_add = (max(T_plot.(measure_name))-min(T_plot.(measure_name)))*0.2;
            for k = 1:length(groups)
                if des_i == 1
                    tm = max(T_plot.(measure_name))+sum([STATS_STRUCT(cnt).pvals{k}(:)]<0.05)*stat_add;
                end
                if des_i == 2
                    tm = max(T_plot.(measure_name))+sum([STATS_STRUCT(cnt).regress_pval{k}]<0.05)*(stat_add*2);
                end
            end
            max_val(i,des_i) = tm;
            cnt = cnt + 1;
        end
    end
    max_vals = max(max_val,[],2);
    %-
    vert_shift = 0;
    cnt = 1;
    horiz_shift= 0;
    for des_i = 1:length(designs)   
        inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
        T_plot = FOOOF_TABLE(inds,:);
        switch des_i
            case 1
                
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
        end
        % horiz_shift= 0;
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            tmp_stats = STATS_STRUCT(cnt);
            % figure;
            % VIOLIN_PARAMS = {'width',0.08,...
            %     'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            %     'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            %     'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            %     'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            if length(groups) == 1
                g_char = [];
            else
                g_char = unique(T_plot.group_char);
            end
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_PARAMS = struct('color_map',color_dark,...
                'cond_labels',unique(T_plot.cond_char),'group_labels',unique(T_plot.group_char),...
                'cond_offsets',cond_offsets,...
                'group_offsets',[0.125,0.475,0.812],...
                'y_label','10*log_{10}(Flattened PSD)',...
                'title',title_plot{i},'font_size',7,'ylim',[min(T_plot.(measure_name))-0.3,max_vals(i)],...
                'font_name','Arial','x_label',x_label,'do_combine_groups',false);
            axax = group_violin(T_plot,measure_name,'cond_id','group_id',...
                fig,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_PARAMS',PLOT_PARAMS,...
                'STATS_STRUCT',tmp_stats);
            if i > 1
                ylabel(axax,'');
            end
            set(axax,'OuterPosition',[0,0,1,1]);
            set(axax,'Position',[0.07+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
            hold off;
            %- iterate
            horiz_shift = horiz_shift + 0.3*im_resize;
            cnt = cnt + 1;
        end
        %## TITLE
        % annotation('textbox',[0.5-0.1,VIOLIN_BOTTOM+vert_shift-0.075+0.225*im_resize,0.2,0.2],...
        %     'String',string(design_chars{des_i}),'HorizontalAlignment','center',...
        %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
        %     'FontSize',14,'FontWeight','Bold','Units','normalized');
        horiz_shift = horiz_shift - 0.3*im_resize+0.175;
    end
    hold off;
    %## SAVE TIFF
    exportgraphics(fig,[save_dir filesep sprintf('cl%i_topo-dips-psd-violins.tiff',cl_i)],'Resolution',1000)
    % exportgraphics(fig,[save_dir filesep sprintf('cl%i_topo-dips-psd-violins.tiff',k)],'Resolution',300)
    close(fig);
end
%% ================================================================= %%
%{
%## PSDS
PSD_BOTTOM = 0.7;
im_resize = 0.5;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
for k_i = 1:length(clusters)
    %##
    cl_i = double(string(clusters(k_i)));
    atlas_name = atlas_name_store{k_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    for j = 1:length(groups)
        %## PLOT SPEEED PSDS
        horiz_shift = 0;
        for des_i = 1:length(designs)
            switch des_i
                case 1
                    color_dark = COLORS_MAPS_TERRAIN;
                    color_light = COLORS_MAPS_TERRAIN;
                    GROUP_CMAP_OFFSET = [0,0.1,0.1];
                    xtick_label_g = {'flat','low','med','high'};
                case 2
                    color_dark = COLOR_MAPS_SPEED;
                    color_light = COLOR_MAPS_SPEED+0.15;
                    GROUP_CMAP_OFFSET = [0.15,0,0];
                    xtick_label_g = {'0.25','0.50','0.75','1.0'};
            end
            %## non-fooof psd (speed)
            axes();
            hold on;
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                switch j
                    case 1
                        [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:),color_light(i,:));
                        Pa.EdgeColor = "none";
                        
                    case 2
                        [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:)+GROUP_CMAP_OFFSET,color_light(i,:)+GROUP_CMAP_OFFSET);
                        Pa.EdgeColor = color_light(i,:)+GROUP_CMAP_OFFSET;
                    case 3
                        [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:)+GROUP_CMAP_OFFSET,color_light(i,:)+GROUP_CMAP_OFFSET);
                        Pa.EdgeColor = color_light(i,:)+GROUP_CMAP_OFFSET;

                end
            end
            axs = [];
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                data = spec_data_original{des_i}{cl_i}{i,j}';
                switch j
                    case 1
                        ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                    case 2
                        ax = plot(fooof_freq,mean(data),'color',color_dark(i,:)+GROUP_CMAP_OFFSET,'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                    otherwise
                end
                axs = [axs, ax];
            end 
            %- plot the aperiodic line
            for i = 1:size(spec_data_original{des_i}{cl_i},1)
                aperiodic_fit = fooof_apfit_store{des_i}{cl_i}{i,j}';
                switch j
                    case 1
                        dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
                    case 2
                        dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:)+GROUP_CMAP_OFFSET,'linestyle','-.','linewidth',2,'displayname','ap. fit');
                    otherwise
                end
            end
            ax = gca;
            xlim([4 40]);
            ylim([-30 -5]);
            switch j
                case 1
                    [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond_org{des_i}{cl_i}{1}(:,2), 'background', 'Frequency(Hz)');
                case 2
                    [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond_org{des_i}{cl_i}{1}(:,2), 'background', 'Frequency(Hz)');
            end
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            title('PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.08+horiz_shift,PSD_BOTTOM-vert_shift,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
            hold off;
        
            %## fooof psd (speed)
            axes();
            hold on;
            axs = [];
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                switch j
                    case 1
                        [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:),color_light(i,:));
                        Pa.EdgeColor = "none";
                        
                    case 2
                        [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                            color_dark(i,:)+GROUP_CMAP_OFFSET,color_light(i,:)+GROUP_CMAP_OFFSET);
                        Pa.EdgeColor = color_light(i,:)+GROUP_CMAP_OFFSET;
                        % Pa.FaceAlpha = 0.2;
                        Pa.LineStyle = ":";
                        Pa.LineWidth = 1;
                    otherwise
                end
            end
            for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
                data = fooof_diff_store{des_i}{cl_i}{i,j}';
                switch j
                    case 1
                        ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                    case 2
                        ax = plot(fooof_freq,mean(data),'color',color_dark(i,:)+GROUP_CMAP_OFFSET,'linewidth',2,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
                    otherwise
                end
                axs = [axs, ax];
            end
            %-
            ax = gca;
            switch j
                case 1
                    [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond{des_i}{cl_i}{j}(:,2), 'background', 'Frequency(Hz)');
                case 2
                    [axsignif,Pa] = highlight_CL(ax, fooof_freq, pcond{des_i}{cl_i}{j}(:,2), 'background', 'Frequency(Hz)');
            end
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('Frequency(Hz)');
            ylabel('10*log_{10}(Power)');
            xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
        %     set(ax,'LineWidth',2)
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
                'FontWeight','bold')
            %- legend
            legend([axs,dash],'FontSize',9,'FontName',PLOT_STRUCT.font_name);
            [lg1,icons,plots,txt] = legend('boxoff');
            set(lg1,'Position',[0.20+0.3*im_resize+horiz_shift,PSD_BOTTOM+0.025-vert_shift,0.2,0.1]);
            lg1.ItemTokenSize(1) = 18;
            %-
            title('Flattened PSD')
            set(ax,'OuterPosition',[0,0,1,1]);
            carry_ov = 0.12+0.3*im_resize;
            set(ax,'Position',[carry_ov+horiz_shift,PSD_BOTTOM-vert_shift,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
        %         icons(2).XData = [0.05 0.1];
            horiz_shift = horiz_shift + carry_ov + 0.25*im_resize + 0.1;
        end
        %## TITLE
        annotation('textbox',[0.5-0.1,PSD_BOTTOM-vert_shift-0.05+0.25*im_resize,0.2,0.2],...
            'String',string(group_chars(j)),'HorizontalAlignment','center',...
            'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            'FontSize',14,'FontWeight','Bold','Units','normalized');
        % close(fig);
        vert_shift = vert_shift + 0.25*im_resize+0.1;
    end
    hold off;
    % exportgraphics(fig,[save_dir filesep sprintf('Group_Violins_cl%i.tiff',cl_i)],'Resolution',1000)
    exportgraphics(fig,[save_dir filesep sprintf('Group_PSDs_cl%i.tiff',cl_i)],'Resolution',300)
    % close(fig);
end
%}