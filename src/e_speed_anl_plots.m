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
clear java
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
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
save_dir = [spec_data_dir filesep 'psd_calcs' filesep 'no_spca_psd_data'];
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
tmp = load([save_dir filesep 'STATS_TRACK_STRUCT_speedlin.mat']);
STATS_TRACK_STRUCT = tmp.STATS_TRACK_STRUCT;
STATS_TRACK_STRUCT = struct2table(STATS_TRACK_STRUCT);
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
clusters = unique(FOOOF_TABLE.cluster_id);
for j = 1:length(clusters) % For each cluster requested
    cluster_i = double(string(clusters(j)));
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
%% ===================================================================== %%
clusters = unique(FOOOF_TABLE.cluster_id);
% [STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
txt_store = cell(length(clusters),1);
atlas_name_store = cell(le
    %%ngth(clusters),1);
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
save([save_dir filesep 'anatomy_chars.mat'],'txt_store','atlas_name_store');
%% ==================================================================== %%
%## MIM KINEMATICS
% meas_names_imu = {'nanmean_APexc_mean','nanmean_MLexc_mean','nanmean_APexc_COV','nanmean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
% meas_names_ls = {'nanmean_StepDur','nanmean_StepDur_cov','nanmean_GaitCycleDur_cov','nanmean_GaitCycleDur',};
%##
SCATTER_BOTTOM = 0.65;
IM_RESIZE = 0.7;
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
%% ===================================================================== %%
%## SPEED MANUSCRIPT GROUP PLOT
% CLUSTERS_TO_PLOT = double(string(clusters)); %main_cl_inds(2:end); %valid_clusters(1:end-1); %main_cl_inds(2:end);
% VIOLIN_BOTTOM = 0.375;
% PSD_BOTTOM = 0.575;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%- 
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
horiz_shift = 0;
vert_shift = 0;
%- 
AX_FONT_NAME = 'Arial';
AX_W = 0.35;
AX_H = 0.25;
AX_SUBTITLE_YOFFSET = -0.065;
AX_SUBTITLE_XOFFSET = 0.01;
AX_HORIZ_SHIFT = 0.08;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.09;
AX_TXT_SIZE = 10;
%-
% LINE_WIDTH_REFF = 1;
%- 
DIP_IM_DPI = 1000;
% AX_INIT_HORIZ_TOPO = 0.5;
% AX_INIT_VERT_DIP = 7.275; % inches for some reason, maybe a bug with normal units
AX_INIT_HORIZ_TOPO = 0.085;
AX_INIT_VERT_TOPO = 0.75;
AX_INIT_VERT_DIP = 0.8; %7.2; % inches for some reason, maybe a bug with normal units
AX_INIT_HORIZ_DIP = 0.3846; %2.5; % inches for some reason, maybe a bug with normal units
%- spec plots
XLIM_PSD = [4,40];
XFREQ_LINES = [3,8,13,30];
%* psd
% AX_INIT_VERT_PSD = 0.575;
AX_INIT_VERT_PSD = 0.585;
LEG_HORIZ_SHIFT_PSD = 0.15;
LEG_VERT_SHIFT_PSD =  0.06;
LINE_ALPHA_PSD = 0.7;
LINE_WIDTH_PSD = 2;
LINE_WIDTH_APPSD = 1;
%* fooof
% AX_INIT_VERT_PSDFF = 0.35;
AX_INIT_VERT_PSDFF = 0.35;
LEG_HORIZ_SHIFT_PSDFF = 0.15;
LEG_VERT_SHIFT_PSDFF =  0.05;
LINE_ALPHA_PSDFF = 0.7;
LINE_WIDTH_PSDFF = 2;
%* vio
AX_INIT_VERT_VIO = 0.08;
%-
LAB_A_YOFFSET = -0.17;
LAB_A_XOFFSET = -0.125;
LAB_B_YOFFSET = 0.065;
LAB_B_XOFFSET = -0.125;
LAB_C_YOFFSET = 0.075;
LAB_C_XOFFSET = -0.125;
LAB_D_YOFFSET = 0.08;
LAB_D_XOFFSET = -0.125;
FIGURE_POSITION =[0,0,6.5,9];
FONT_NAME = 'Arial';
FONT_WEIGHT_PSD = 'normal';
TOPO_FONTSIZE = 7;
FONT_SIZE_PSD = 10;
FONT_SIZE_PSD_LEG = 10;
FONT_SIZE_VIO = 10;
FONT_SIZE_VIO_REG = 6;
% PSD_GROUP_TITLE_FONTSIZE = 12;
PSD_GROUP_TITLE_FONTSIZE = 10;
LEG_TOKEN_SIZE_PSD = 20;
%-
des_i = 2;
cl_i = 3;
s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
for i = 1:length(STUDY.design(des_i).variable)
    if strcmp(STUDY.design(des_i).variable(i).label,'cond')
        c_chars = STUDY.design(des_i).variable(i).value;
    elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
        g_chars = STUDY.design(des_i).variable(i).value;
    end
end
cluster_titles = {'Precuneus','Right Posterior Parietal',...
    'Left Occipital','Left Supplementary Motor','Left Sensorimotor','Left Posterior Parietal',...
    'Eye','Left Temporal','Mid/Posterior Cingulate','Right Sensorimotor','Right Temporal'};
psd_ylimits = {[-31.5,-10],[-31.5,-15], [-30,-12.5], [-32.5,-15], [-32.5,-15],...
    [-30,-12.5], [-30,-10], [-30,-10], [-30,-10], [-32.5,-15], [-30,-10]};
c_chars = {'0.25m/s','0.5m/s','0.75m/s','1.0m/s'};
% g_chars_topo = {'Young Adult',{'Older Adult','High Mobility'},{'Older Adult','Low Mobility'}};
g_chars_topo = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
g_chars_subp = {'YA','OHMA','OLMA'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};

%## TOPO & DIPOLE PLOTS
for k_i = 1:length(clusters)
    %%
    cl_i = double(string(clusters(k_i)));
    %## ANATOMY
    atlas_name = cluster_titles{k_i}; %atlas_name_store{k_i};
    %## AXES LIMITS
    fig = figure('color','white','renderer','Painters');
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    % hh = sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    % hh.Position
    set(fig,'Units','inches','Position',FIGURE_POSITION)
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:});
    hold on;
    %## ALIGNMENT
    % GAP = 0.05;
    IM_RESIZE = 0.8;
    %## TOPO PLOT
    % AX_INIT_HORIZ_DIP = 2.25; % inches for some reason, maybe a bug with normal units
    % AX_INIT_VERT_DIP = 7.2; % inches for some reason, maybe a bug with normal units
    % subplot(3,4,1)
    axes();
    set(gca,AXES_DEFAULT_PROPS{:});
    std_topoplot_CL(STUDY,cl_i,'together');
    colormap(linspecer); %colormap_ersp)
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'color','w')
    g_counts = cell(length(groups),1);
    for group_i = 1:length(groups)
        g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
        if length(g_chars_topo{group_i}) == 1 || ischar(g_chars_topo{group_i})
            g_counts{group_i} =sprintf('%s N=%i',g_chars_topo{group_i},sum(g_inds));
        else
            g_counts{group_i} =sprintf('%s\n%s N=%i',g_chars_topo{group_i}{1},g_chars_topo{group_i}{2},sum(g_inds));
        end
    end
    fig_i.Children(1).Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    fig_i.Children(1).Title.Interpreter = 'none';
    fig_i.Children(1).FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    fig_i.Children(1).FontName = FONT_NAME;
    fig_i.Children(1).FontWeight = 'bold';
    fig_i.Children(1).OuterPosition = [0,0,1,1];
    fig_i.Children(1).Units = 'Normalized';
    fig_i.Children(1).Position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,0.225*IM_RESIZE,0.25*IM_RESIZE];  %[left bottom width height]

    %## DIPOLE PLOTS
    %-
    target_sz = 1.25; %inch
    target_dim = 1; %2==width, 1==height
    im1 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cl_i)]);
    im2 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cl_i)]);
    im3 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cl_i)]);
    im1d1 = size(im1,1)/DIP_IM_DPI;
    im1d2 = size(im1,2)/DIP_IM_DPI;
    im2d1 = size(im2,1)/DIP_IM_DPI;
    im2d2 = size(im2,2)/DIP_IM_DPI;
    im3d1 = size(im3,1)/DIP_IM_DPI;
    im3d2 = size(im3,2)/DIP_IM_DPI;
    %## SCALE
    %- 1
    if target_dim==1
        scale = target_sz/im1d1;
        im1 = imresize(im1,scale);
    elseif target_dim==2
        scale = target_sz/im1d2;
        im1 = imresize(im1,scale);
    end
    %- 2
    im2(1:300,:,:) = [];
    im2(end-300:end,:,:)=[];
    if target_dim==1
        scale = target_sz/im2d1;
        im2 = imresize(im2,scale);
    elseif target_dim==2
        scale = target_sz/im2d2;
        im2 = imresize(im2,scale);
    end
    %- 3
    if target_dim==1
        scale = target_sz/im3d1;
        im3 = imresize(im3,scale);
    elseif target_dim==2
        scale = target_sz/im3d2;
        im3 = imresize(im3,scale);
    end
    %- ensure proper crops (automatic, not complete)
    % im1d1 = size(im1,1);
    % im1d2 = size(im1,2);
    % im2d1 = size(im2,1);
    % im2d2 = size(im2,2);
    % im3d1 = size(im3,1);
    % im3d2 = size(im3,2);
    % if im1d1 ~= im2d1 || im1d1 ~= im3d1 || im2d1 ~= im3d1
    %     max([im1d1,im2d1,im3d1])
    %     min([im1d1,im2d1,im3d1])
    % end
    %- ensure proper crops, helps fixes alignment issues (manual)
    %{
    %- use to determine padding
    im1d1 = size(im1,1);
    im1d2 = size(im1,2);
    im2d1 = size(im2,1);
    im2d2 = size(im2,2);
    im3d1 = size(im3,1);
    im3d2 = size(im3,2);
    %}
    im1(1:40,:,:) = [];
    im2 = padarray(im2,[31,0],0,"post");
    im2 = padarray(im2,[10,0],0,"both");
    im3(1:40,:,:) = [];
    im1(im1==255) = 0;
    im2(im2==255) = 0;
    im3(im3==225) = 0;
    %-
    % szw1 = (size(im1,2)/DIP_IM_DPI); %/PG_SIZE(1); %width
    % szh1 = (size(im1,1)/DIP_IM_DPI); %/PG_SIZE(2); %+0.0001; %height
    % szw2 = (size(im2,2)/DIP_IM_DPI); %/PG_SIZE(1); %width
    % szh2 = (size(im2,1)/DIP_IM_DPI); %+0.05; %/PG_SIZE(2); %+0.01;
    % szw3 = (size(im3,2)/DIP_IM_DPI); %/PG_SIZE(1);
    % szh3 = (size(im3,1)/DIP_IM_DPI); %/PG_SIZE(2);
    %-
    % AX_INIT_HORIZ_DIP = AX_INIT_HORIZ_DIP/PG_SIZE(1);
    % AX_INIT_VERT_DIP = AX_INIT_VERT_DIP/PG_SIZE(2);
    szw1 = (size(im1,2)/DIP_IM_DPI)/PG_SIZE(1); %width
    szh1 = (size(im1,1)/DIP_IM_DPI)/PG_SIZE(2); %+0.0001; %height
    szw2 = (size(im2,2)/DIP_IM_DPI)/PG_SIZE(1); %width
    szh2 = (size(im2,1)/DIP_IM_DPI)/PG_SIZE(2); %+0.01;
    szw3 = ((size(im3,2)/DIP_IM_DPI)/PG_SIZE(2)); %((size(im3,2)/DIP_IM_DPI)/PG_SIZE(1))/
    szh3 = (size(im3,1)/DIP_IM_DPI)/PG_SIZE(2);
    szw = max([szw1,szw2,szw3]);
    szh = max([szh1,szh2,szh3]);
    % szsz = max([szw,szh]);
    szh1 = szh1*(szh1/szh);
    % szw1 = szw1*(szw1/szw);
    szh2 = szh2*(szh2/szh);
    % szw2 = szw2*(szw/szh);
    szh3 = szh3*(szh3/szh);
    szw3 = szw3*(szw/szh);
    %-
    % axes('Units','normalized','OuterPosition',[0,0,1,1],...
    %     'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,szw1,szh1],'PositionConstraint','outerposition');
    axes('Units','normalized','OuterPosition',[0,0,1,1],...
        'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,szw1,szh1]);
    set(gca,AXES_DEFAULT_PROPS{:});
    imshow(im1,'border','tight');
    %-
    % axes('Units','normalized','OuterPosition',[0,0,1,1],...
    %     'Position',[AX_INIT_HORIZ_DIP+szw1+((szw-szw1)-(szw-szw2)),AX_INIT_VERT_DIP,szw1,szh2],'PositionConstraint','outerposition');
    axes('Units','normalized','OuterPosition',[0,0,1,1],...
        'Position',[AX_INIT_HORIZ_DIP+szw1+0.0025,AX_INIT_VERT_DIP,szw2,szh2]);
    set(gca,AXES_DEFAULT_PROPS{:});
    imshow(im2,'border','tight');
    
    %-
    % axes('Units','normalized','OuterPosition',[0,0,1,1],...
    %     'Position',[AX_INIT_HORIZ_DIP+szw1+szw2+0.05,AX_INIT_VERT_DIP,szw,szh3],'PositionConstraint','outerposition');
    axes('Units','normalized','OuterPosition',[0,0,1,1],...
        'Position',[AX_INIT_HORIZ_DIP+szw1+szw2-0.005,AX_INIT_VERT_DIP,szw3,szh3]);
    set(gca,AXES_DEFAULT_PROPS{:});
    imshow(im3,'border','tight');
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1],...
        'String','A)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    %%
    %## (PSDS) ====================================================== %%
    % im_resize = 0.5;
    GROUP_EDGECOLOR = {};
    GROUP_LINESTYLE = {};
    
    %## PLOT NONFOOOF PSDS (SPEED) ===================================== %%
    AX_W = 0.325;
    AX_H = 0.25;
    IM_RESIZE = 0.7;
    AX_HORIZ_SHIFT = 0.1;
    horiz_shift = 0;
    %## AUTO SET PSD YLIM
    tmp_data = cell(length(groups),1);
    des_i = 2;
    cnt = 1;
    for j = 1:length(groups)
        for i = 1:size(spec_data_original{des_i}{cl_i},1)
            tmp_data{cnt} = spec_data_original{des_i}{cl_i}{i,j}';
            cnt = cnt + 1;
        end
    end
    tmp_data = cat(1,tmp_data{:});
    Y_LIM_GROUP = [round(prctile(tmp_data,6,'all'),1),prctile(tmp_data,94,'all')];
    %- end auto set psd ylim
    for j = 1:length(groups)
        des_i = 2;
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
                xtick_label_g = c_chars; %{'0.25','0.50','0.75','1.0'};
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
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),...
                'linewidth',LINE_WIDTH_PSD,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
            set(ax,'Color',[color_dark(i,:),LINE_ALPHA_PSD]);
            axs = [axs, ax];
        end 
        %- plot the aperiodic line
        dash = [];
        for i = 1:size(spec_data_original{des_i}{cl_i},1)
            aperiodic_fit = fooof_apfit_store{des_i}{cl_i}{i,j}';
            tmp = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),...
                'linestyle','-.','linewidth',LINE_WIDTH_APPSD,'displayname','ap. fit');
            set(tmp,'Color',[color_dark(i,:),LINE_ALPHA_PSD]);
            if i == 1
                dash = tmp;
            end
        end
        %- setting xlim & ylims before adding sig. shading
        xlim(gca,XLIM_PSD);
        % ylim(gca,psd_ylimits{k_i}(:));
        ylim(gca,Y_LIM_GROUP)
        %- sig. shading
        [axsignif,Pa] = highlight_CL(gca, fooof_freq, pcond_org{des_i}{cl_i}{1}(:,2), 'background', 'Frequency(Hz)');
        %- freq. band lines
        plot(XLIM_PSD,[0 0],'--','color','black');
        for xx = 1:length(XFREQ_LINES)
            xline(XFREQ_LINES(xx),'--');
        end
        %- fonts & ax position
        set(gca,'FontName',FONT_NAME,'FontSize',FONT_SIZE_PSD,...
            'FontWeight',FONT_WEIGHT_PSD)
        set(gca,'OuterPosition',[0,0,1,1]);
        %- legend
        if j == length(groups)
            legend([axs,dash],'FontSize',FONT_SIZE_PSD_LEG,'FontName',FONT_NAME);
            [lg1,icons,plots,txt] = legend('boxoff');
            % set(lg1,'Position',[0.20+0.3*im_resize+horiz_shift,PSD_BOTTOM+0.025-vert_shift,0.2,0.1]);
            set(lg1,'Position',[LEG_HORIZ_SHIFT_PSD+AX_INIT_HORIZ+horiz_shift,...
                LEG_VERT_SHIFT_PSD+AX_INIT_VERT_PSD+vert_shift,lg1.Position(3),lg1.Position(4)]);
            lg1.ItemTokenSize(1) = LEG_TOKEN_SIZE_PSD;
        end
        %## LABELS
        % xlabel('Frequency(Hz)','FontWeight','bold');
        xlabel('');
        if j == 1
            ylabel('10*log_{10}(PSD)','FontWeight','bold','FontSize',FONT_SIZE_PSD);
        else
            ylabel('','FontWeight','bold','FontSize',FONT_SIZE_PSD);
        end
        %## TITLE
        % title('PSD','FontSize',10)
        title([sprintf('%s',string(g_chars_subp(j))),' PSD'],...
            'FontSize',PSD_GROUP_TITLE_FONTSIZE,'FontWeight','Bold');
        % annotation('textbox',[AX_INIT_HORIZ+horiz_shift+AX_SUBTITLE_XOFFSET,AX_INIT_VERT_PSD-vert_shift+AX_SUBTITLE_YOFFSET+AX_H*IM_RESIZE,0.2,0.2],...
        %     'String',string(string(g_chars_subp(j))),'HorizontalAlignment','center',...
        %     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        %     'FontSize',PSD_GROUP_TITLE_FONTSIZE,'FontWeight','Bold','Units','normalized');
        set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_PSD+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_B_XOFFSET+(0.1/2),AX_INIT_VERT_PSD+LAB_B_YOFFSET+(0.1/2),.1,.1],...
        'String','B)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    %## PLOT FOOOF PSDS (SPEED) ===================================== %%
    % AX_W = 0.35;
    % AX_H = 0.25;
    IM_RESIZE = 0.55;
    AX_H  = 0.3;
    AX_W = 0.414;
    AX_HORIZ_SHIFT = 0.125;
    horiz_shift = 0;
    vert_shift = 0;
    %## AUTO SET PSD YLIM
    tmp_data = cell(length(groups),1);
    des_i = 2;
    cnt = 1;
    for j = 1:length(groups)
        for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
            tmp_data{cnt} = fooof_diff_store{des_i}{cl_i}{i,j}';
            cnt = cnt + 1;
        end
    end
    tmp_data = cat(1,tmp_data{:});
    Y_LIM_GROUP = [round(prctile(tmp_data,1,'all'),1),prctile(tmp_data,99,'all')+prctile(tmp_data,99,'all')*0.01];
    %- end auto set psd ylim
    for j = 1:length(groups)
        des_i = 2;
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
                xtick_label_g = c_chars; %{'0.25','0.50','0.75','1.0'};
        end
        %- PLOT
        axes();
        hold on;
        axs = [];
        %- standard error shading
        for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
            data = fooof_diff_store{des_i}{cl_i}{i,j}';
            std_error = std(data)/sqrt(size(data,1));
            [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),mean(data)-std_error,mean(data)+std_error,...
                color_dark(i,:),color_light(i,:));
            Pa.EdgeColor = "none";
        end
        for i = 1:size(fooof_diff_store{des_i}{cl_i},1)
            data = fooof_diff_store{des_i}{cl_i}{i,j}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),...
                'linewidth',LINE_WIDTH_PSDFF,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
            set(ax,'Color',[color_dark(i,:),LINE_ALPHA_PSDFF]);
            axs = [axs, ax];
        end
        %- set limits before sig. shading
        ylim(Y_LIM_GROUP)
        xlim(XLIM_PSD);
        ax = gca;
        [axsignif,Pa] = highlight_CL(gca, fooof_freq, pcond{des_i}{cl_i}{j}(:,2), 'background', 'Frequency(Hz)');
        plot(XLIM_PSD,[0 0],'--','color','black');     
        for xx = 1:length(XFREQ_LINES)
            xline(XFREQ_LINES(xx),'--');
        end
        set(gca,'FontName',FONT_NAME,'FontSize',FONT_SIZE_PSD,...
            'FontWeight',FONT_WEIGHT_PSD)
        %- legend
        if j == length(groups)
            legend(axs,'FontSize',FONT_SIZE_PSD_LEG,'FontName',FONT_NAME);
            [lg2,icons,plots,txt] = legend('boxoff');
            set(lg2,'Position',[LEG_HORIZ_SHIFT_PSDFF+AX_INIT_HORIZ+horiz_shift,...
                LEG_VERT_SHIFT_PSDFF+AX_INIT_VERT_PSDFF+vert_shift,lg2.Position(3),lg2.Position(4)]);
            lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE_PSD;
        end
        %## LABELS
        xlabel('Frequency(Hz)','FontWeight','bold','FontSize',FONT_SIZE_PSD);
        if j == 1
            ylabel('10*log_{10}(PSD - AP. Fit)','FontWeight','bold','FontSize',FONT_SIZE_PSD);
        else
            ylabel('','FontWeight','bold','FontSize',FONT_SIZE_PSD);
        end
        %## TITLE
        % title('Flattened PSD','FontSize',10)
        title([sprintf('%s',string(g_chars_subp(j))),' Flattened PSD'],'FontSize',PSD_GROUP_TITLE_FONTSIZE,'FontWeight','Bold');
        set(gca,'OuterPosition',[0,0,1,1]);
        % carry_ov = 0.12+0.3*im_resize;
    %     set(ax,'Position',[carry_ov+horiz_shift,PSD_BOTTOM-vert_shift,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
    % %         icons(2).XData = [0.05 0.1];
    %     horiz_shift = horiz_shift + carry_ov + 0.25*im_resize + 0.1;
        set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_PSDFF+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        
        %## TITLE
        % annotation('textbox',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_PSDFOOOF-vert_shift+AX_SUBTITLE_OFFSET+AX_H*IM_RESIZE,0.2,0.2],...
        %     'String',string(group_chars(j)),'HorizontalAlignment','center',...
        %     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        %     'FontSize',14,'FontWeight','Bold','Units','normalized');
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_C_XOFFSET+(0.1/2),AX_INIT_VERT_PSDFF+LAB_C_YOFFSET+(0.1/2),.1,.1],...
        'String','C)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    hold on;
    %## VIOLIN PLOTS) ================================================== %%
    % im_resize= 0.9;
    IM_RESIZE = 0.9;
    AX_H  = 0.2;
    AX_W = 0.275;
    AX_HORIZ_SHIFT = 0.06;
    %## violin plot's theta/alpha/beta (speed)
    %-
    max_val = zeros(length(measure_name_plot),length(designs));
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
    des_i = 2;
    % inds = STATS_TRACK_STRUCT.design == designs(des_i) & STATS_TRACK_STRUCT.cluster == clusters(k_i);
    % tmp_table = STATS_TRACK_STRUCT(inds,:);
    % prc_ylim = [prctile(tmp_fooof_t.(measure_name_plot{meas_i}),5),prctile(tmp_fooof_t.(measure_name_plot{meas_i}),95)];
    for meas_i = 1:length(measure_name_plot)
        
        % measure_name_plot{meas_i} = measure_name_plot{meas_i};
        inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
        tmp_fooof_t = FOOOF_TABLE(inds,:);
        prc_ylim = [floor(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),1))+floor(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),99))*.3,...
            ceil(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),99))+ceil(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),99))*.1];
        STATS_STRUCT(cnt) = DEFAULT_STATS_STRUCT;
        for group_i = 1:length(groups)
            inds = STATS_TRACK_STRUCT.design == designs(des_i) & STATS_TRACK_STRUCT.cluster == clusters(k_i) &...
                STATS_TRACK_STRUCT.group == groups(group_i) & strcmp(STATS_TRACK_STRUCT.measure_tag,measure_name_plot{meas_i});
            tmp_table = STATS_TRACK_STRUCT(inds,:);
            switch des_i
                case 1
                    aa = tmp_table.anova_terr_p{1};
                    c2s = tmp_table.lme_terr_p{1}{1}(1);
                    c3s = tmp_table.lme_terr_p{1}{1}(2);
                    c4s = tmp_table.lme_terr_p{1}{1}(3);
                    rs = [];
                    rls = [tmp_table.lme_inter_coeff, tmp_table.lme_terr_coeff{1}{1}(1),...
                        tmp_table.lme_terr_coeff{1}{1}(2), tmp_table.lme_terr_coeff{1}{1}(3)]; 
                    r2 = tmp_table.R2{1};
                    norm_p = tmp_table.norm_test_p;
                    STATS_STRUCT(cnt).anova{group_i}=aa;
                    STATS_STRUCT(cnt).pvals{group_i}=[1,c2s,c3s,c4s];
                    STATS_STRUCT(cnt).pvals_pairs{group_i}={[1,1],[1,2],[1,3],[1,4]};
                case 2
                    aa = tmp_table.anova_speed_p{1};
                    c2s = [];
                    c3s = [];
                    c4s = [];
                    rs = tmp_table.lme_speed_p{1};
                    rls = [tmp_table.lme_inter_coeff{1}, tmp_table.lme_speed_coeff{1}]; 
                    r2 = tmp_table.R2{1};
                    norm_p = tmp_table.norm_test_p;
                    STATS_STRUCT(cnt).anova{group_i}=aa;
                    STATS_STRUCT(cnt).regress_pval{group_i}=rs;
                    STATS_STRUCT(cnt).regress_line{group_i}=rls;
                    STATS_STRUCT(cnt).r2_coeff(group_i)=r2;
                    STATS_STRUCT(cnt).regress_xvals=(0:5)*0.25;
            end
        end
        stat_add = (max(tmp_fooof_t.(measure_name_plot{meas_i}))-min(tmp_fooof_t.(measure_name_plot{meas_i})))*0.2;
        for group_i = 1:length(groups)
            switch des_i
                case 1
                    tm = ceil(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),97));
                    % tm = max(T_plot.(measure_name_plot{meas_i}))+sum([STATS_STRUCT(cnt).pvals{k}(:)]<0.05)*stat_add;
                    % tm = max(tmp_fooof_t.(measure_name_plot{meas_i}))+1*std(tmp_fooof_t.(measure_name_plot{meas_i}));
                case 2
                    tm = ceil(prctile(tmp_fooof_t.(measure_name_plot{meas_i}),97));
                    % tm = max(T_plot.(measure_name_plot{meas_i}))+sum([STATS_STRUCT(cnt).regress_pval{k}]<0.05)*(stat_add*2);
                    % tm = max(tmp_fooof_t.(measure_name_plot{meas_i}))+1*std(tmp_fooof_t.(measure_name_plot{meas_i}));
            end
        end
        max_val(meas_i,des_i) = tm;
        cnt = cnt + 1;
    end
    max_vals = max(max_val,[],2)+1*std(max_val,[],2);
    %-
    vert_shift = 0;
    cnt = 1;   
    inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
    tmp_fooof_t = FOOOF_TABLE(inds,:);
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
    for meas_i = 1:length(measure_name_plot)
        ax = axes();
        tmp_stats = STATS_STRUCT(cnt);
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_STRUCT = struct('color_map',color_dark,...
            'cond_labels',unique(tmp_fooof_t.cond_char),'group_labels',categorical(g_chars_subp),...
            'cond_offsets',cond_offsets,...
            'group_offsets',[0.125,0.475,0.812],...
            'y_label','10*log_{10}(Flattened PSD)',...
            'title',title_plot{meas_i},'font_size',FONT_SIZE_VIO,'ylim',[min(tmp_fooof_t.(measure_name_plot{meas_i}))-0.5,max_vals(meas_i)],...
            'font_name','Arial','x_label',x_label,'do_combine_groups',false,...
            'regresslab_txt_size',FONT_SIZE_VIO_REG);
        % PLOT_PARAMS = struct('color_map',color_dark,...
        %     'cond_labels',unique(tmp_fooof_t.cond_char),'group_labels',unique(tmp_fooof_t.group_char),...
        %     'cond_offsets',cond_offsets,...
        %     'group_offsets',[0.125,0.475,0.812],...
        %     'y_label','10*log_{10}(Flattened PSD)',...
        %     'title',title_plot{meas_i},'font_size',PLOT_STRUCT.font_size,'ylim',prc_ylim,...
        %     'font_name','Arial','x_label',x_label,'do_combine_groups',false,...
        %     'regresslab_txt_size',9);
        % ax = axes();
        % figfig = figure();
        ax = group_violin(tmp_fooof_t,measure_name_plot{meas_i},'cond_id','group_id',...
            ax,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',tmp_stats);
        if meas_i ~= 1
            ylabel('');
        end
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        %## TITLE
        % annotation('textbox',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO-vert_shift+AX_SUBTITLE_OFFSET+AX_H*IM_RESIZE,0.2,0.2],...
        % 'String',string(group_chars(j)),'HorizontalAlignment','center',...
        % 'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        % 'FontSize',14,'FontWeight','Bold','Units','normalized');
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        % set(axax,'Position',[0.08+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        % hold off;
        % %- iterate
        % horiz_shift = horiz_shift + 0.3*im_resize+0.05;
        %- set ylabel & title
        AXES_FONT_SIZE_VIO = 10;
        GROUP_LAB_FONTSIZE = AXES_FONT_SIZE_VIO;
        GROUP_LAB_FONTWEIGHT = 'bold ';
        XLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
        YLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
        XTICK_FONTSIZE = AXES_FONT_SIZE_VIO;
        XLAB_FONTWEIGHT = 'bold';
        YLAB_FONTWEIGHT = 'bold';
        TITLE_FONTSIZE = AXES_FONT_SIZE_VIO;
        TITLE_FONTWEIGHT = 'bold';
        XLABEL_OFFSET = -.05;
        GROUP_LAB_YOFFSET = -0.275;
        % ax = gca;
        ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
        ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
        ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
        ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
        ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
        ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
        % set(ax,'FontName','Arial','FontSize',XTICK_FONTSIZE,'FontWeight','normal')
        % yticks(ax,'FontSize',XTICK_FONTSIZE)
        yt = yticks(ax);
        xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
        pos1=get(xlh,'Position');
        pos1(1,2)=pos1(1,2)+XLABEL_OFFSET;
        set(xlh,'Position',pos1);
        if meas_i == 1
            ylabel(ax,'10*log_{10}(PSD - AP. Fit)','FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
        else
            ylabel(ax,'','FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
        end
        title(ax,PLOT_STRUCT.title,...
            'FontSize',TITLE_FONTSIZE,'FontWeight',TITLE_FONTWEIGHT);
        % fontsize(ax,10,"points");
        % ylim(gca,PLOT_STRUCT.ylim);
        cnt = cnt + 1;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_D_XOFFSET+(0.1/2),AX_INIT_VERT_VIO+LAB_D_YOFFSET+(0.1/2),.1,.1],...
        'String','D)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    %## TITLE
    % annotation('line',[AX_INIT_HORIZ+0.0775,1-AX_INIT_HORIZ-0.0875],repmat(AX_INIT_VERT_VIO+AX_H*IM_RESIZE+0.035,[1,2]),...
    %     'Color','black','Units','normalized','LineStyle','-','LineWidth',2);
    % annotation('rectangle',[AX_INIT_HORIZ-0.075,AX_INIT_VERT_VIO-0.075,1-(AX_INIT_HORIZ-0.075),AX_H*IM_RESIZE+0.1],...
    %     'Color','black','Units','normalized','LineStyle','-','LineWidth',3);
    % vert_shift = vert_shift - (0.1+0.225*im_resize);
    hold off;
    exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i.tiff',cl_i)],'Resolution',1000)
    exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i.png',cl_i)],'Resolution',300)
    % exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i.tiff',cl_i)],'Resolution',300)
    % close(fig);
end
