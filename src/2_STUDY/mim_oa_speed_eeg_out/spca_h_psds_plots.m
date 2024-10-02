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
SUB_GROUP_FNAME = ['group_spec' filesep 'split_band_test'];
% SUB_GROUP_FNAME = 'all_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## SET STUDY PATHS
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
% if ~isempty(SUB_GROUP_FNAME)
%     spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
% else
%     spec_data_dir = [cluster_dir filesep 'spec_data'];
% end
% % save_dir = [spec_data_dir filesep 'psd_calcs'];
% save_dir = [spec_data_dir filesep 'psd_calcs' filesep 'split_band_test'];
% if ~exist(save_dir,'dir')
%     mkdir(save_dir);
% end
if ~isempty(SUB_GROUP_FNAME)
    save_dir = [cluster_dir filesep 'psd_calcs' filesep SUB_GROUP_FNAME];
else
    save_dir = [cluster_dir filesep 'psd_calcs'];
end
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
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
% cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
%% RE-POP PARAMS
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%% (LOAD EXISTING TALBES && FORMAT STUDY)
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
% tmp = load([save_dir filesep 'STATS_TRACK_STRUCT_speedlin.mat']);
tmp = load([save_dir filesep 'STATS_TRACK_STRUCT_speedlin_alt.mat']);
% STATS_TRACK_STRUCT = tmp.STATS_TRACK_STRUCT;
% STATS_TRACK_STRUCT = struct2table(STATS_TRACK_STRUCT);
vio_stats_struct = tmp.stats_struct;
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
%% (TOPO PLOTS) ======================================================== %%
CALC_STRUCT = struct('cluster_inds',(2:length(STUDY.cluster)),...
    'save_inf',true,...
    'recalculate',true);
[STUDY,dipfit_structs,topo_cells] = eeglab_get_topodip(STUDY,...
    'CALC_STRUCT',CALC_STRUCT,...
    'ALLEEG',ALLEEG);
%% (ANATOMY) =========================================================== %%
addpath([PATHS.submods_dir filesep 'AAL3']);
ANATOMY_STRUCT = struct('atlas_fpath',{{[PATHS.submods_dir filesep 'AAL3' filesep 'AAL3v1.nii']}},...
    'group_chars',{unique({STUDY.datasetinfo.group})},...
    'cluster_inds',2:length(STUDY.cluster),...
    'anatomy_calcs',{{'all aggregate','all centroid'}},... % ('all calcs','group centroid','all centroid','group aggregate','all aggregate')
    'save_dir',cluster_dir,...
    'save_inf',true);
[STUDY,anat_struct,~,~,txt_out] = eeglab_get_anatomy(STUDY,...
    'ANATOMY_STRUCT',ANATOMY_STRUCT);
%-
% atlas_char = 'AAL3v1.nii';
% inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'aggregate label for all');
% at_out={anat_struct(inds).anatomy_label};
% cl_tmp = STUDY.cluster;
% [cl_tmp(2:length(at_out)).agg_anat_all] = at_out{:};
%-
% inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'centroid label for all');
% at_out={anat_struct(inds).anatomy_label};
% at_ind = [anat_struct(inds).cluster];
% cl_tmp = STUDY.cluster;
% [cl_tmp(at_ind).ctr_anat_all] = at_out{:};
%% (LOAD) ============================================================== %%
atlas_char = 'AAL3v1.nii';
inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'centroid label for all');
at_out = {anat_struct(inds).anatomy_label};
atlas_name_store = at_out;
% cl_struct = STUDY.cluster;
% par_save(cl_struct,cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
%% ===================================================================== %%
%## SPEED MANUSCRIPT GROUP PLOT
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
meas_names = {'theta_avg_power','alpha_avg_power','beta_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% measure_name_plot = {'theta_avg_power','alpha2_avg_power','beta2_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \alpha','Mean \beta'};
% measure_name_plot = {'med_sub_flat','low_sub_flat','high_sub_flat'};
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
FIGURE_POSITION =[1,1,6.5,9];
PG_SIZE = [6.5,9];
FONT_NAME = 'Arial';
FONT_WEIGHT_PSD = 'normal';

YLABEL_FONT_SIZE_PSD = 10;
YLABEL_FONT_SIZE_PSDFF = 10;
XLABEL_FONT_SIZE_PSDFF = 10;
TOPO_FONTSIZE = 7;
FONT_SIZE_PSD = 9;
FONT_SIZE_PSD_LEG = 9;
% PSD_GROUP_TITLE_FONTSIZE = 12;
PSD_GROUP_TITLE_FONTSIZE = 10;
LEG_TOKEN_SIZE_PSD = 20;
%##    
VIO_PLOT_STRUCT = struct('color_map',[],...
            'cond_labels',{{'0.25','0.50','0.75','1.0'}},...
            'cond_offsets',[-0.35,-0.1,0.15,0.40],...
            'group_labels',{{'YA','OHMA','OLMA'}},...
            'group_offsets',[0.125,0.475,0.812],...
            'group_lab_yoffset',-0.23,...
            'group_lab_fontweight','normal',...
            'group_lab_fontsize',10,...
            'y_label','unit',...
            'y_label_fontsize',10,...
            'y_label_fontweight','bold',...
            'ylim',[],...
            'x_label','Speed (m/s)',...
            'x_label_fontsize',10,...
            'x_label_fontweight','bold',...
            'x_label_yoffset',-.14,...
            'xlim',[],...
            'title',{{''}},...
            'title_fontsize',10,...
            'title_fontweight','bold',...
            'font_size',9,...
            'font_name','Arial',...
            'do_combine_groups',false,...
            'regresslab_txt_size',6,...
            'ax_position',[0,0,1,1],...
            'ax_line_width',1);
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
% cluster_new_labs = {'Precuneus_R','Postcentral_R','Calcarine_L/Lingaul_L','Frontal_Sup_2_L',...
%     'Precentral_L','Parietal_Inf_L','Putamen_R/Caudate_R/Frontal_Mid_2_R','Insula_L',...
%     'Supp_Motor_Area_R/error','Frontal_Sup_2_R','Rolandic_Oper_R'};
% cluster_titles = {'Precuneus','Right Posterior Parietal',...
%     'Left Occipital','Left Supplementary Motor','Left Sensorimotor','Left Posterior Parietal',...
%     'Eye','Left Temporal','Mid/Posterior Cingulate','Right Sensorimotor','Right Temporal'};
% cluster_titles = {'Precuneus','Right Sensorimotor',...
%     'Left Occipital','Left Supplementary Motor','Left Sensorimotor','Left Posterior Parietal',...
%     'Eye','Left Temporal','Mid Cingulate','Right Supplementary Motor','Right Temporal'};
% cluster_titles = atlas_name_store;
%- (09/04/2024) ICLabel Chosen Brain Areas
% cluster_titles = {'Precuneus','Right Supplementary Motor',...
%     'Left Sensorimotor','Left Occipital','Right Temporal','Right Sensorimotor',...
%     'Left Temporal','Mid Cingulate','Left Posterior Parietal','Right Posterior Parietal','Left Supplementary Motor'};
%- (09/8/2024) ICLabel & kmeans bug fix
cluster_titles = {'Right Occipital','Left Occipital','Mid Cingulate',...
    'Right Sensorimotor','Left Supplementary','Precuneus','Left Temporal','Left Sensorimotor',...
    'Right Posterior Parietal','Left Posterior Parietal','Right Temporal'};
%## 
% psd_ylimits = {[-31.5,-10],[-32.5,-15],...
%     [-30,-12.5], [-32.5,-15], [-32.5,-15],[-30,-12.5],...
%     [-30,-10],[-30,-10],[-30,-10],[-32.5,-15],[-30,-10]};
% psd_ylimits = {[-31.5,-10],[-32.5,-15],...
%     [-37.5,-15], [-32.5,-15], [-32.5,-15],[-30,-12.5],...
%     [-30,-10],[-30,-10],[-30,-10],[-32.5,-15],[-30,-10]};
% violin_ylimits_theta = {[-1,6],[-1.5,4],...
%     [-1.5,5], [], [],[],...
%     [],[],[],[],[]};
% violin_ylimits_alpha = {[-1,12.5],[-1,12.5],...
%     [-1,12.5], [], [],[],...
%     [],[],[],[],[]};
% violin_ylimits_beta = {[-1,7],[-1,9],...
%     [-1.5,7], [], [],[],...
%     [],[],[],[],[]};
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
% g_chars_topo = {'Young Adult',{'Older Adult','High Mobility'},{'Older Adult','Low Mobility'}};
g_chars_topo = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
g_chars_subp = {'YA','OHMA','OLMA'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
dip_dir = [cluster_dir filesep 'topo_dip_inf'];
%%
des_i = 2;
%## TOPO & DIPOLE PLOTS
for k_i = 1:length(clusters)
    %%
    cl_i = double(string(clusters(k_i)));
    %## ANATOMY
    atlas_name = cluster_titles{k_i};
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
    for g_i = 1:length(groups)
        g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
        if length(g_chars_topo{g_i}) == 1 || ischar(g_chars_topo{g_i})
            g_counts{g_i} =sprintf('%s N=%i',g_chars_topo{g_i},sum(g_inds));
        else
            g_counts{g_i} =sprintf('%s\n%s N=%i',g_chars_topo{g_i}{1},g_chars_topo{g_i}{2},sum(g_inds));
        end
    end
    % fig_i.Children(1).Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    % fig_i.Children(1).Title.Interpreter = 'none';
    % fig_i.Children(1).FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    % fig_i.Children(1).FontName = FONT_NAME;
    % fig_i.Children(1).FontWeight = 'bold';
    % fig_i.Children(1).OuterPosition = [0,0,1,1];
    % fig_i.Children(1).Units = 'Normalized';
    % fig_i.Children(1).Position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,0.225*IM_RESIZE,0.25*IM_RESIZE];  %[left bottom width height]
    obj_i = findobj(fig_i,'type','Axes');
    obj_i(1).Title.String = g_counts;
    obj_i(1).Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    obj_i(1).Title.Interpreter = 'none';
    obj_i(1).FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    obj_i(1).FontName = FONT_NAME;
    obj_i(1).FontWeight = 'bold';
    obj_i(1).OuterPosition = [0,0,1,1];
    obj_i(1).Units = 'Normalized';
    obj_i(1).Position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,0.225*IM_RESIZE,0.25*IM_RESIZE];  %[left bottom width height]

    %## DIPOLE PLOTS
    %-
    target_sz = 1.25; %inch
    target_dim = 1; %2==width, 1==height
    
    im1 = imread([dip_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cl_i)]);
    im2 = imread([dip_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cl_i)]);
    im3 = imread([dip_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cl_i)]);
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
    % im1(im1==255) = 0;
    % im2(im2==255) = 0;
    % im3(im3==225) = 0;
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
                'linestyle','-.','linewidth',LINE_WIDTH_APPSD,'displayname','AP. Fit');
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
            ylabel('10*log_{10}(PSD)','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSD);
        else
            ylabel('','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSD);
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
        % des_i = 2;
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
        %##
        xlabel('Frequency(Hz)','FontWeight','bold','FontSize',XLABEL_FONT_SIZE_PSDFF);
        if j == 1
            ylabel('10*log_{10}(PSD) - AP. Fit','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSDFF);
        else
            ylabel('','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSDFF);
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
    IM_RESIZE = 0.9;
    AX_H  = 0.2;
    AX_W = 0.275;
    AX_HORIZ_SHIFT = 0.06;
    % figure;
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
            speed_xvals = []
    end
    %## violin plot's theta/alpha/beta (speed)
    DEFAULT_STATS_STRUCT = struct('anova',{{}},...
        'anova_grp',{{}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{{}},...
        'pvals_grp_pairs',{{}},...
        'regress_pval',{{}},...
        'regress_line',{{}},...
        'line_type',{'best_fit'},... % ('best_fit' | 'means')
        'regress_xvals',speed_xvals,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'display_stats_char',true,...
        'stats_char',{{}},...
        'bracket_conn_yshift',[1,1,1],...
        'bracket_rawshifty_upper',0,...
        'bracket_rawshifty_lower',0,...
        'grp_sig_offset_x',[0,0,0],... %zeros(length(unique(tmp_table.(GROUP_TABLE_VAR)))),...
        'grp_sig_offset_y',[0,0,0]);
    %
    %##
    STATS_STRUCT = DEFAULT_STATS_STRUCT;
    cnt = 1;
    vert_shift = 0;    
    horiz_shift= 0;
    prc_ylim = zeros(length(meas_names),2);
    for meas_i = 1:length(meas_names)
        % measure_name_plot{meas_i} = measure_name_plot{meas_i};
        inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
        tmp_fooof_t = FOOOF_TABLE(inds,:);
        prc_ylim(meas_i,:) = [floor(prctile(tmp_fooof_t.(meas_names{meas_i}),1))-floor(std(tmp_fooof_t.(meas_names{meas_i}))),...
            ceil(prctile(tmp_fooof_t.(meas_names{meas_i}),99))+ceil(std(tmp_fooof_t.(meas_names{meas_i})))*1.5];
        disp(prc_ylim)
        STATS_STRUCT(cnt) = DEFAULT_STATS_STRUCT;
        for g_i = 1:length(groups)
            inds = vio_stats_struct.design_tag == designs(des_i) & vio_stats_struct.cluster_tag == clusters(k_i) &...
                vio_stats_struct.group_tag == groups(g_i) & vio_stats_struct.measure_tag==meas_names{meas_i};
            tmp_table = vio_stats_struct(inds,:);
            switch des_i
                case 1
                    anova_p = tmp_table.anova_preds_p{2};
                    terr_p = tmp_table.mod_preds_p(2:4);
                    speed_r2 = tmp_table.mod_r2;
                    coeffs = tmp_table.mod_preds_coeff(1:4);

                    aa = tmp_table.anova_terr_p{1};
                    c2s = tmp_table.lme_terr_p{1}{1}(1);
                    c3s = tmp_table.lme_terr_p{1}{1}(2);
                    c4s = tmp_table.lme_terr_p{1}{1}(3);
                    rs = [];
                    rls = [tmp_table.lme_inter_coeff, tmp_table.lme_terr_coeff{1}{1}(1),...
                        tmp_table.lme_terr_coeff{1}{1}(2), tmp_table.lme_terr_coeff{1}{1}(3)]; 
                    r2 = tmp_table.R2{1};
                    norm_p = tmp_table.norm_test_p;
                    STATS_STRUCT(cnt).anova{g_i}=aa;
                    STATS_STRUCT(cnt).pvals{g_i}=[1,c2s,c3s,c4s];
                    STATS_STRUCT(cnt).pvals_pairs{g_i}={[1,1],[1,2],[1,3],[1,4]};
                case 2
                    anova_p = tmp_table.anova_preds_p{1}(2);
                    reg_p = tmp_table.mod_preds_p{1}(2);
                    reg_lin = [tmp_table.mod_preds_coeff{:}];
                    r2 = tmp_table.mod_r2{1};
                    aa = anova_p; %tmp_table.anova_speed_p{1};
                    rs = reg_p; %tmp_table.lme_speed_p{1};
                    rls = reg_lin; %[tmp_table.lme_inter_coeff{1}, tmp_table.lme_speed_coeff{1}]; 
                    % r2 = r2; %tmp_table.R2{1};
                    if aa > 0.01 && aa < 0.05
                        % str = sprintf('* %0.1g\{times}<speed>+%0.1g\nR^2=%0.2g',rls(2),rls(1),r2);
                        % str = [sprintf('* %0.1g',rls(2)),'\times<speed>+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                        str = [sprintf('* %0.1g',rls(2)),'x+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                    elseif aa <= 0.01 && aa > 0.001
                        % str = sprintf('* %0.1g\{times}<speed>+%0.1g\nR^2=%0.2g',rls(2),rls(1),r2);
                        % str = [sprintf('* %0.1g',rls(2)),'\times<speed>+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                        str = [sprintf('** %0.1g',rls(2)),'x+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                    elseif aa <= 0.001
                        % str = [sprintf('* %0.1g',rls(2)),'\times<speed>+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                        str = [sprintf('*** %0.1g',rls(2)),'x+',sprintf('%0.1g\nR^2=%0.2g',rls(1),r2)];
                    else
                        str = '';
                    end
                    STATS_STRUCT(cnt).line_type = 'best_fit';
                    STATS_STRUCT(cnt).stats_char{g_i} = str;
                    STATS_STRUCT(cnt).display_stats_char = true;
                    norm_p = tmp_table.norm_test_p;
                    STATS_STRUCT(cnt).anova{g_i}=aa;
                    STATS_STRUCT(cnt).regress_pval{g_i}=rs;
                    STATS_STRUCT(cnt).regress_line{g_i}=rls;
                    % STATS_STRUCT(cnt).r2_coeff(group_i)=r2;
                    STATS_STRUCT(cnt).regress_xvals=(0:5)*0.25;
            end
        end
        %## PLOT
        ax = axes();
        tmp_stats = STATS_STRUCT(cnt);
        VIO_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        %-
        VIO_PLOT_STRUCT.color_map = color_dark;
        VIO_PLOT_STRUCT.cond_labels = xtick_label_g;
        VIO_PLOT_STRUCT.title = title_plot(meas_i);
        if meas_i == 1
            VIO_PLOT_STRUCT.y_label ='10*log_{10}(PSD) - AP. Fit';
        else
            VIO_PLOT_STRUCT.y_label ='';
        end
        VIO_PLOT_STRUCT.ylim = prc_ylim(meas_i,:);
        VIO_PLOT_STRUCT.ax_position = [AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        %-
        ax = group_violin(tmp_fooof_t,meas_names{meas_i},'cond_id','group_id',...
            ax,...
            'VIOLIN_PARAMS',VIO_PARAMS,...
            'PLOT_STRUCT',VIO_PLOT_STRUCT,...
            'STATS_STRUCT',tmp_stats);
        if meas_i ~= 1
            ylabel('');
        end
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        cnt = cnt + 1;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_D_XOFFSET+(0.1/2),AX_INIT_VERT_VIO+LAB_D_YOFFSET+(0.1/2),.1,.1],...
        'String','D)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    hold off;
    drawnow;
    exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i_des%i.tiff',cl_i,des_i)],'Resolution',1000)
    % exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i_des%i.png',cl_i,des_i)],'Resolution',300)
    savefig(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i_des%i.fig',cl_i,des_i)])
    % exportgraphics(fig,[save_dir filesep sprintf('Group_speed_violin_psds_cl%i.tiff',cl_i)],'Resolution',300)
    % close(fig);
end

%% MULTI-CLUSTER PLOT OF ALL SUBJECTS ================================== %%
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
conditions = unique(FOOOF_TABLE.cond_char);
groups = unique(FOOOF_TABLE.group_char);
%- (09/8/2024) ICLabel & kmeans bug fix
cluster_titles = {'Right Occipital','Left Occipital','Mid Cingulate',...
    'Right Sensorimotor','Left Supplementary','Precuneus','Left Temporal','Left Sensorimotor',...
    'Right Posterior Parietal','Left Posterior Parietal','Right Temporal'};
%-
des_i = 2;
cl_i = 8;
clust_i = double(string(clusters(cl_i)));
%-
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_subp = {'YA','OHMA','OLMA'};
AX_W = 0.4;
AX_H = 0.25;
IM_RESIZE = 0.5;
HZ_DIM = 3;
VERTICAL_SHIFT =  0.45;
HORIZONTAL_SHIFT = 0.5;
AX_INIT_HORIZ = 0.175;
AX_INIT_VERTICAL = 0.75;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%## PLOT ============================================================ %%
fig = figure('color','white','renderer','Painters');
sgtitle(sprintf('%s',cluster_titles{cl_i}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
set(gca,AXES_DEFAULT_PROPS{:})
hold on;
%-
horiz_shift = AX_INIT_HORIZ;
hz = 0;
vert_shift = AX_INIT_VERTICAL;
%##
switch des_i
    case 1
        color_dark = COLOR_MAPS_TERRAIN;
        color_light = COLOR_MAPS_TERRAIN;
        GROUP_CMAP_OFFSET = [0,0.1,0.1];
        xtick_label_g = {'flat','low','med','high'};
    case 2
        color_dark = COLOR_MAPS_SPEED;
        color_light = COLOR_MAPS_SPEED+0.15;
        GROUP_CMAP_OFFSET = [0.15,0,0];
        xtick_label_g = {'0.25','0.50','0.75','1.0'};
end

for c_i = 1:4
    cond_i = c_i; %double(string(conditions(c_i)));
    cond_char = string(conditions(c_i));
    for g_i = 1:length(groups)
        group_i = g_i; %double(string(groups(g_i)));
        %-
        hz = hz + 1;
        ax = axes();
        hold on;
        plot([0 40],[0 0],'--','color','black');
        xlim([4 40]);
        ylim([-2 10]);
        %-
        fooof_psd = fooof_diff_store{des_i}{clust_i}{cond_i,group_i}';
        fooof_psd_mean = mean(fooof_psd);
        subjs = plot(fooof_freq,fooof_psd,'color',[0,0,0,0.15],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
        mean_plot = plot(fooof_freq,fooof_psd_mean,'color',color_dark(cond_i,:),'linestyle','-','linewidth',4,'displayname','orig. subj psd');
        %-
        xlabel('Frequency(Hz)');
        if hz ~= 1
            ylabel('');
        else
            ylabel('10*log_{10}(PSD) - AP. Fit');
        end
        set(ax,'FontName','Arial',...
            'FontSize',12,...
            'FontWeight','bold')
        xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
        set(ax,'FontName','Arial','FontSize',10,...
            'FontWeight','bold')
        title(sprintf('%s',g_chars_subp{group_i}),'FontWeight','normal')
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[horiz_shift,vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        %## TITLE CONDITION
        annotation('textbox',[0.5-(0.1/2),vert_shift-(0.1/2)+AX_H*IM_RESIZE,.1,.1],...
            'String',sprintf('%s',c_chars{cond_i}),'HorizontalAlignment','left',...
            'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
            'FontSize',12,'FontWeight','Bold','Units','normalized');
        if hz < HZ_DIM
            horiz_shift = horiz_shift + HORIZONTAL_SHIFT*IM_RESIZE;
        else
            vert_shift = vert_shift - VERTICAL_SHIFT*IM_RESIZE;
            horiz_shift = AX_INIT_HORIZ;
            hz = 0;
        end
    end
end
hold off;
exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_1group_comparison_cluster_fooofs_subjs_means.tiff',cl_i,des_i,cond_i,group_i)],'Resolution',1000)
% close(fig);
