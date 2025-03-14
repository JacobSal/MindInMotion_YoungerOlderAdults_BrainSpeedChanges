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
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SRC_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
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
    SRC_DIR = fileparts(SCRIPT_DIR);    
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%
%- colors
cmaps_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmaps_terrain = [cmaps_terrain(3,:);custom_yellow;cmaps_terrain(4,:);cmaps_terrain(2,:)];
cmaps_speed = linspecer(4*3);
cmaps_speed = [cmaps_speed(1,:);cmaps_speed(2,:);cmaps_speed(3,:);cmaps_speed(4,:)];
%- statistics & conditions
SPEED_VALS = {'0.25','0.5','0.75','1.0';
              '0p25','0p5','0p75','1p0'};
TERRAIN_VALS = {'flat','low','med','high'};
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
    'freqrange',[1,200]);
%## FOOOF
% settings = struct('peak_width_limits',[1,8],...
%     'min_peak_height',0.05,...
%     'max_n_peaks',3);
settings = struct('peak_width_limits',[1,8],...
    'min_peak_height',0.05,...
    'max_n_peaks',5);
f_range = [3, 40];
theta_band_lims = [4, 8];
alpha_band_lims = [8 12];
beta_band_lims  = [12 30];
alpha1_band_lims = [8,10.5];
alpha2_band_lims = [10.5,13];
beta1_band_lims = [13,20];
beta2_band_lims = [20,30];
%% (PATHS) ============================================================= %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load study file
study_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME)];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [study_fpath filesep '__iclabel_cluster_kmeansalt_rb5'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- save directory 
ANALYSIS_DNAME = 'spca_fooof_psd_anl';
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
% end
%-
% cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
% STUDY.cluster = cl_struct;
% [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
% CLUSTER_PICS = main_cl_inds;
%%
%(01/16/2025) JS, this will remove all topo information from the cluster
%structure. Moving cluster loading after the std_makedesign.m call.
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}}};
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%## LOAD CLUSTER DATA
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%% (LOAD EXISTING TALBES && FORMAT STUDY)
tmp = load([save_dir filesep 'psd_feature_table.mat']);
fooof_table = tmp.FOOOF_TABLE;
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
%% (LOAD) ============================================================== %%
% atlas_char = 'AAL3v1.nii';
% inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'centroid label for all');
% at_out = {anat_struct(inds).anatomy_label};
% atlas_name_store = at_out;
% cl_struct = STUDY.cluster;
% par_save(cl_struct,cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
%% ===================================================================== %%
%## SPEED MANUSCRIPT GROUP PLOT
designs = unique(fooof_table.design_id);
clusters = unique(fooof_table.cluster_id);
groups = unique(fooof_table.group_id);
meas_names = {'theta_avg_power','alpha_avg_power','beta_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% measure_name_plot = {'theta_avg_power','alpha2_avg_power','beta2_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \alpha','Mean \beta'};
% measure_name_plot = {'med_sub_flat','low_sub_flat','high_sub_flat'};
%-
TITLE_TXT_SIZE = 14;
y_shift = 0;
%- 
AX_FONT_NAME = 'Arial';
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
AX_INIT_VERT_TOPO = 0.765;
AX_INIT_VERT_DIP = 0.83; %0.79; %0.8 --- %7.2; % inches for some reason, maybe a bug with normal units
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
AX_INIT_VERT_PSDFF = 0.365; %0.37
LEG_HORIZ_SHIFT_PSDFF = 0.15;
LEG_VERT_SHIFT_PSDFF =  0.05;
LINE_ALPHA_PSDFF = 0.7;
LINE_WIDTH_PSDFF = 2;
%* vio
AX_INIT_VERT_VIO = 0.09; %0.08
%-
LAB_A_YOFFSET = -0.16;
LAB_A_XOFFSET = -0.125;
LAB_B_YOFFSET = 0.065;
LAB_B_XOFFSET = -0.125;
LAB_C_YOFFSET = 0.06; %0.075
LAB_C_XOFFSET = -0.125;
LAB_D_YOFFSET = 0.09;
LAB_D_XOFFSET = -0.125;
FIGURE_POSITION =[1,1,6.5,9];
PG_SIZE = [6.5,9];
FONT_NAME = 'Arial';
FONT_WEIGHT_PSD = 'normal';

YLABEL_FONT_SIZE_PSD = 10;
YLABEL_FONT_SIZE_PSDFF = 10;
XLABEL_FONT_SIZE_PSDFF = 10;
TOPO_FONTSIZE = 8;
FONT_SIZE_PSD = 9;
FONT_SIZE_PSD_LEG = 9;
% PSD_GROUP_TITLE_FONTSIZE = 12;
PSD_GROUP_TITLE_FONTSIZE = 10;
LEG_TOKEN_SIZE_PSD = 20;
%## VIOLIN PLOT STRUCT
VIO_PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{'0.25','0.50','0.75','1.0'}},...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_labels',{{'YA','OHFA','OLFA'}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.23,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label','',...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label','Speed (m/s)',...
    'x_label_fontsize',10,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.155,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',10,...
    'title_fontweight','bold',...
    'font_size',9,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'regresslab_txt_size',7,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
%-
des_i = 2;
cl_i = 3;
s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
desdes = cat(1,STUDY.design.variable);
% c_chars = desdes(strcmp({desdes.label},'cond'));
% c_chars = {c_chars.value};
g_chars = desdes(strcmp({desdes.label},'group'));
g_chars = {g_chars.value};
g_chars = g_chars{1};

%- rb10
% cluster_titles = {'','','Right Sensorimotor','Mid Cingulate','Left Temporal','Right Occipital',...
%     'Right Premotor','Left Occipital','Left Sensorimotor','Right Posterior Parietal',...
%     'Left Posterior Parietal',...
%     'Left Supplementary Motor','Right Temporal'};
%- rb5
cluster_titles = {'Anterior Cingulate','Right Sensorimotor','Left Occipital','Left Supplementary Motor',...
    'Mid Cingulate','Cuneus','Left Sensorimotor','Left Posterior Parietal',...
    'Right Temporal',...
    'Right Occipital','Right Posterior Parietal'};
output_titles = {'lsm','rppa','midc','rcun','rsm','lsupm','rocp','locp','ltemp','lppa','rtemp'};
fig_n = 1:length(cluster_titles); %[0,0,7,10,11,10,8,12,0,0,0,9,0];
%##
violin_ylimits = {{[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]}};
% violin_ylimits = {{[-1.5,7],[-1.5,7],[-1.5,5],... % theta
%         [-2,8],[-1.5,5.5],[-1.5,6],[],...
%         [],...
%         [],[-2,8],[-1,5]};...
%             {[-2.5,16],[-2.5,16],[-2,10],... % theta
%         [-2.5,16],[-2.5,16],[-2,7],[],...
%         [],...
%         [],[-2.5,16],[-2.5,11]};...
%             {[-2,10.5],[-2,10.5],[-2,8],... % theta
%         [-2,8],[-2,10.5],[-2,7],[],...
%         [],...
%         [],[-2,8],[-1,6]}};
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
% g_chars_topo = {'Young Adult',{'Older Adult','High Mobility'},{'Older Adult','Low Mobility'}};
g_chars_topo = {'Young Adults','Older High Functioning Adults','Older Low Functioning Adults'};
g_chars_subp = {'YA','OHFA','OLFA'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%-
% dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'cl3-cl13'];
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'figure_gen'];
%-
% cluster_inds_plot = [3,4,6,7,9,10,11,12]
cluster_inds_plot = [3,4,6,7,8,9,10,13];

%%
des_i = 2;
%## TOPO & DIPOLE PLOTS
for ii = 1:length(cluster_inds_plot)
    %%
    k_i = find(cluster_inds_plot(ii) == double(string(clusters)));
    % k_i = cluster_inds(ii);
    cl_i = double(string(clusters(k_i)));
    %## ANATOMY
    atlas_name = cluster_titles{k_i};
    %## AXES LIMITS
    fig = figure('color','white','Renderer','Painters');
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
    p_sz = get(fig,'Position');
    % set(fig,'Units','normalized');
    set(gca,AXES_DEFAULT_PROPS{:});
    % fig_ax = gca;
    hold on;
    %## ALIGNMENT
    IM_RESIZE = 0.225;
    %## TOPO PLOT
    [~,h] = std_topoplot_CL(STUDY,cl_i,'together');
    set(h,'color','w')
    set(h,'Units','normalized');    
    ax = get(h,'CurrentAxes');
    colormap(h,linspecer);   
    pos = get(ax,'Position');
    opos = get(ax,'OuterPosition');
    ac = get(ax,'Children');
    ax1 = axes('Parent',fig,...
        'DataAspectRatio',get(ax,'DataAspectRatio'));
    copyobj(ac,ax1);
    set(ax1,AXES_DEFAULT_PROPS{:});
    colormap(ax1,linspecer)
    g_counts = cell(length(groups),1);
    for g_i = 1:length(groups)
        g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
        if length(g_chars_topo{g_i}) == 1 || ischar(g_chars_topo{g_i})
            g_counts{g_i} =sprintf('%s N=%i',g_chars_topo{g_i},sum(g_inds));
        else
            g_counts{g_i} =sprintf('%s\n%s N=%i',g_chars_topo{g_i}{1},g_chars_topo{g_i}{2},sum(g_inds));
        end
    end
    ax1.Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    ax1.Title.Interpreter = 'none';
    ax1.FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    ax1.FontName = FONT_NAME;
    ax1.FontWeight = 'bold';
    ax1.OuterPosition = [0,0,1,1];
    ax1.Units = 'Normalized';
    pp = get(ax1,'Position');    
    ax1.Position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,pp(3)*IM_RESIZE,pp(4)*IM_RESIZE];  %[left bottom width height]
    close(h)
    %## FIGURE IMPLEMENTATION ========================================== %%
    % AX_INIT_HORIZ_DIP = 0.075;
    % AX_INIT_VERT_DIP = 0.6;
    %-
    IM_RESIZE = 1.1;
    % AX1_ADJ = 1.05;
    % AX2_ADJ = 1.55;
    % AX3_ADJ = 1.2725;
    % %-
    % AX1_X_O = 0; %0.19;
    % AX1_Y_O = -0.2; %-0.17;
    % AX2_X_O = 0.11; %0.19;
    % AX2_Y_O = -0.175; %-0.17;
    % AX3_X_O = 0.15; %0.24;
    % AX3_Y_O = -0.11;
    %-
    tmp = openfig([dip_dir filesep sprintf('%i_dipplot_alldipspc_top.fig',cl_i)]);
    cnt = length(tmp.Children(end).Children);
    try
        while cnt > 3
            % tmp.Children(end).Children(cnt).LineWidth = 0.4;
            tmp.Children(end).Children(cnt).SizeData = 7; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
            tmp.Children(end).Children(cnt).LineWidth = 0.1; %'none';
            tmp.Children(end).Children(cnt).CData = [0,0,0]; %tmp.Children(end).Children(cnt-1).MarkerFaceColor;
            % fig.Children(end).Children(cnt-1).MarkerFaceColor = DIPPLOT_STRUCT.color{1}; %fig.Children(end).Children(cnt-1+1).Color;
            tmp.Children(end).Children(cnt).MarkerFaceAlpha = 0.75;
            cnt = cnt - 1;
        end
    catch
        fprintf('not scatter 3d');
    end
    ax = get(tmp,'CurrentAxes');
    view(ax,[90,0])
    daspect = get(ax,'DataAspectRatio');
    ac = get(ax,'Children');
    XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
    YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
    ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];
    %## AXES 1
    ax1 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
    copyobj(ac,ax1);
    set(ax1,'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
    ax1.XRuler.FirstCrossoverValue = 0;
    ax1.YRuler.FirstCrossoverValue = 0;
    ax1.ZRuler.FirstCrossoverValue = 0;
    view([0,90])
    camzoom(ax1,1.1^2)
    set(ax1,'YLim',YLIM)
    set(ax1,'XLim',XLIM)
    drawnow;
    pp = get(ax1,'Position');
    %-
    AX1_ADJ = 1.02;
    AX1_Y_O = -0.145; %-0.17;
    AX1_X_O = -0.04; %-0.17;
    ax1_x = pp(3)*IM_RESIZE*AX1_ADJ/p_sz(3);
    ax1_y = pp(4)*IM_RESIZE*AX1_ADJ/p_sz(3);
    set(ax1,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x*AX1_X_O,AX_INIT_VERT_DIP+ax1_y*AX1_Y_O,ax1_x,ax1_y]);
    % ax1_x = ax1_x;
    % ax1_h = ax1_h*AX1_W_O;
    %## AXES 2
    ax2 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
    copyobj(ac,ax2);
    %-
    set(ax2,'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
    ax2.XRuler.FirstCrossoverValue = 0;
    ax2.YRuler.FirstCrossoverValue = 0;
    ax2.ZRuler.FirstCrossoverValue = 0;
    set(ax2,'View',[90,0]);
    camzoom(ax2,1.1^2)
    set(ax2,'YLim',YLIM)
    set(ax2,'ZLim',ZLIM)
    drawnow;
    pp2 = get(ax2,'Position');
    AX2_ADJ = 1.55;
    AX2_Y_O = -0.175;
    AX2_X_O = 0.15;
    ax2_x = pp2(3)*IM_RESIZE*AX2_ADJ/p_sz(3);
    ax2_y = pp2(4)*IM_RESIZE*AX2_ADJ/p_sz(4);    
    set(ax2,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+(ax2_x*AX2_X_O),AX_INIT_VERT_DIP+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);    
    %## AXES 3
    ax3 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
    copyobj(ac,ax3);
    %-
    set(ax3,'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
    ax3.XRuler.FirstCrossoverValue = 0;
    ax3.YRuler.FirstCrossoverValue = 0;
    ax3.ZRuler.FirstCrossoverValue = 0;
    set(ax3,'view',[0,0])
    set(ax3,'ZLim',ZLIM)
    set(ax3,'XLim',XLIM)
    camzoom(ax3,1.1^2)
    pp3 = get(ax3,'Position');
    AX3_ADJ = 1.2725;
    AX3_Y_O = -0.11;
    AX3_X_O = 0.22;
    ax3_x = pp3(3)*IM_RESIZE*AX3_ADJ/p_sz(3);
    ax3_y = pp3(4)*IM_RESIZE*AX3_ADJ/p_sz(4);    
    set(ax3,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),AX_INIT_VERT_DIP+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
    close(tmp)
    %## LETTER
    annotation(fig,'textbox',[AX_INIT_HORIZ+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1],...
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
    x_shift = 0;
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
                color_dark = cmaps_speed;
                color_light = cmaps_speed+0.15;
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
            set(lg1,'Position',[LEG_HORIZ_SHIFT_PSD+AX_INIT_HORIZ+x_shift,...
                LEG_VERT_SHIFT_PSD+AX_INIT_VERT_PSD+y_shift,lg1.Position(3),lg1.Position(4)]);
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
        set(gca,'Position',[AX_INIT_HORIZ+x_shift,AX_INIT_VERT_PSD+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
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
    x_shift = 0;
    y_shift = 0;
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
                color_dark = cmaps_speed;
                color_light = cmaps_speed+0.15;
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
            set(lg2,'Position',[LEG_HORIZ_SHIFT_PSDFF+AX_INIT_HORIZ+x_shift,...
                LEG_VERT_SHIFT_PSDFF+AX_INIT_VERT_PSDFF+y_shift,lg2.Position(3),lg2.Position(4)]);
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
        set(gca,'Position',[AX_INIT_HORIZ+x_shift,AX_INIT_VERT_PSDFF+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        
        %## TITLE
        % annotation('textbox',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_PSDFOOOF-vert_shift+AX_SUBTITLE_OFFSET+AX_H*IM_RESIZE,0.2,0.2],...
        %     'String',string(group_chars(j)),'HorizontalAlignment','center',...
        %     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        %     'FontSize',14,'FontWeight','Bold','Units','normalized');
        x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_C_XOFFSET+(0.1/2),AX_INIT_VERT_PSDFF+LAB_C_YOFFSET+(0.1/2),.1,.1],...
        'String','C)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    hold on;
    %## VIOLIN PLOTS) ================================================== %%
    IM_RESIZE = 0.9;
    AX_H  = 0.22;
    AX_W = 0.28;
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
            color_dark = cmaps_speed; %color.speed;
            color_light = cmaps_speed+0.15; %color.speed_shade;
            xtick_label_g = {'0.25','0.50','0.75','1.0'};
            x_label = 'speed (m/s)';
            cond_offsets = [-0.35,-0.1,0.15,0.40];
            speed_xvals = (0:5)*0.25;
    end
    %## VIOLIN PLOTS
    VIO_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,...
        'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,...
        'Bandwidth',0.15,...
        'QuartileStyle','shadow',...
        'HalfViolin','full',...
        'DataStyle','scatter',...
        'MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],...
        'ViolinAlpha',{0.2 0.3}};
    DEF_STATS_STRUCT = struct('anova',{{}},...
        'anova_grp',{{}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{{}},...
        'pvals_grp_pairs',{{}},...
        'regress_pval',{{}},...
        'regress_line',{{}},...
        'line_type',{'best_fit'},... % ('best_fit' | 'means')
        'regress_xvals',0,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'display_stats_char',false,... 
        'stats_char',{{}});
    DEF_BRACKET_STRUCT = struct('sig_sign','+',...
        'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
        'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
        'bracket_conn',[],...
        'conn_offset_y_upper',[],...
        'bracket_offset_y_upper',0,...
        'bracket_offset_y_lower',0,...
        'sig_levels',[0.05,0.01,0.001],...
        'sig_offset_x',0,...
        'sig_offset_y',[]);
    DEF_SIGLINE_STRUCT = struct('sig_sign','*',...
        'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
        'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
        'conn_y',[],...
        'conn_offset_y',[],...
        'sig_levels',[0.05,0.01,0.001],...
        'sig_offset_x',0,...
        'sig_offset_y',0); 
    %-
    STATS_STRUCT = DEF_STATS_STRUCT;
    cnt = 1;
    y_shift = 0;    
    x_shift= 0;
    prc_ylim = zeros(length(meas_names),2);
    
    for meas_i = 1:length(meas_names)
        % measure_name_plot{meas_i} = measure_name_plot{meas_i};
        inds = fooof_table.design_id == num2str(des_i) & fooof_table.cluster_id == num2str(cl_i);
        tmp_fooof_t = fooof_table(inds,:);
        if isempty(violin_ylimits{meas_i}{k_i})
            prc_ylim(meas_i,:) = [floor(prctile(tmp_fooof_t.(meas_names{meas_i}),1))-floor(std(tmp_fooof_t.(meas_names{meas_i}))),...
                ceil(prctile(tmp_fooof_t.(meas_names{meas_i}),99))+ceil(std(tmp_fooof_t.(meas_names{meas_i})))*1.5];
            disp(prc_ylim)
        else
            prc_ylim(meas_i,:) = violin_ylimits{meas_i}{k_i};
        end
        STATS_STRUCT(cnt) = DEF_STATS_STRUCT;
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
        VIO_PLOT_STRUCT.ax_position = [AX_INIT_HORIZ+x_shift,AX_INIT_VERT_VIO+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        %-
        ax = group_violin(tmp_fooof_t,meas_names{meas_i},'cond_id','group_id',...
            ax,...
            'VIOLIN_PARAMS',VIO_PARAMS,...
            'PLOT_STRUCT',VIO_PLOT_STRUCT,...
            'STATS_STRUCT',tmp_stats,...
            'BRACKET_STRUCT',DEF_BRACKET_STRUCT,...
            'SIGLINE_STRUCT',DEF_SIGLINE_STRUCT);
        if meas_i ~= 1
            ylabel('');
        end
        x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        cnt = cnt + 1;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_D_XOFFSET+(0.1/2),AX_INIT_VERT_VIO+LAB_D_YOFFSET+(0.1/2),.1,.1],...
        'String','D)','HorizontalAlignment','left',...
        'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    hold off;
    drawnow;
    exportgraphics(fig,[save_dir filesep sprintf('figure_%i_%s.tif',fig_n(k_i),output_titles{k_i})],'Resolution',600)
    % exportgraphics(fig,[save_dir filesep sprintf('figure_%i_%s.pdf',fig_n(k_i),output_titles{k_i})],'ContentType','vector')
    close(fig);
end

%% MULTI-CLUSTER PLOT OF ALL SUBJECTS ================================== %%
designs = unique(fooof_table.design_id);
clusters = unique(fooof_table.cluster_id);
conditions = unique(fooof_table.cond_char);
groups = unique(fooof_table.group_char);
%- (09/8/2024) ICLabel & kmeans bug fix
% cluster_titles = {'Right Occipital','Left Occipital','Mid Cingulate',...
%     'Right Sensorimotor','Left Supplementary','Precuneus','Left Temporal','Left Sensorimotor',...
%     'Right Posterior Parietal','Left Posterior Parietal','Right Temporal'};
%-
cluster_titles = {'Left Sensorimotor','Right Posterior Parietal','Mid Cingulate',...
    'Right Cuneus','Right Sensorimotor','Left Supplemeantary Motor','Right Occipital',...
    'Left Occipital',...
    'Left Temporal','Left Posterior Parietal','Right Temporal'};
%-
des_i = 2;
cl_i = 1;
clust_i = double(string(clusters(cl_i)));
%-
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_subp = {'YA','OHFA','OLFA'};
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
x_shift = AX_INIT_HORIZ;
hz = 0;
y_shift = AX_INIT_VERTICAL;
%##
switch des_i
    case 1
        color_dark = cmaps_terrain;
        color_light = cmaps_terrain;
        GROUP_CMAP_OFFSET = [0,0.1,0.1];
        xtick_label_g = {'flat','low','med','high'};
    case 2
        color_dark = cmaps_speed;
        color_light = cmaps_speed+0.15;
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
        ylim([-2 13]);
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
        set(ax,'Position',[x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        %## TITLE CONDITION
        annotation('textbox',[0.5-(0.1/2),y_shift-(0.1/2)+AX_H*IM_RESIZE,.1,.1],...
            'String',sprintf('%s',c_chars{cond_i}),'HorizontalAlignment','left',...
            'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
            'FontSize',12,'FontWeight','Bold','Units','normalized');
        if hz < HZ_DIM
            x_shift = x_shift + HORIZONTAL_SHIFT*IM_RESIZE;
        else
            y_shift = y_shift - VERTICAL_SHIFT*IM_RESIZE;
            x_shift = AX_INIT_HORIZ;
            hz = 0;
        end
    end
end
hold off;
% exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_1group_comparison_cluster_fooofs_subjs_means.tiff',cl_i,des_i,cond_i,group_i)],'ContentType','image','Resolution',600)
exportgraphics(fig,[save_dir filesep sprintf('figure_subject_fooof_example.pdf')],'ContentType','vector')
exportgraphics(fig,[save_dir filesep sprintf('figure_subject_fooof_example.tiff')],'ContentType','image','Resolution',1000)
% export_fig(fig,[save_dir filesep sprintf('figure_subject_fooof_example.pdf')],'ContentType','vector','Resolution',1000)

% close(fig);
%% (SUBJECT DIPOLES PLOTS) ============================================= %%
%-
AX_INIT_HORIZ_DIP = 0.075;
AX_INIT_VERT_DIP = 0.6;
AX_INIT_VERT_DIP2 = 0.775;
%## FIGURE IMPLEMENTATION
%-
% AX1_RESIZE = 0.25;
% AX2_RESIZE = 0.39;
% AX3_RESIZE = 0.302;
IM_RESIZE = 1.5;
AX2_ADJ = 1.52;
AX3_ADJ = 1.255;
%-
AX2_X_O = 0.19;
AX2_Y_O = -0.17;
AX3_X_O = 0.24;
AX3_Y_O = -0.1;
%-
FONT_NAME = 'Arial';
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
% dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'figure_gen'];
%-
% cluster_titles = {'','','Left Sensorimotor','Right Posterior Parietal','Mid Cingulate',...
%     'Precuneus','Right Sensorimotor','Left Supplemeantary Motor','Right Occipital',...
%     'Left Occipital',...
%     'Left Temporal','Left Posterior Parietal','Right Temporal'};
% cluster_titles = {'','','Right Posterior Parietal','Mid Cingulate','Left Temporal','Right Occipital',...
%     'Right Supplementary Motor','Left Occipital','Left Sensorimotor','Cuneus',...
%     'Left Posterior Parietal',...
%     'Left Supplementary Motor','Right Temporal'};
% cluster_titles = {'','','Right Sensorimotor','Mid Cingulate','Left Temporal','Right Occipital',...
%     'Right Premotor','Left Occipital','Left Sensorimotor','Right Posterior Parietal',...
%     'Left Posterior Parietal',...
%     'Left Supplementary Motor','Right Temporal'};
cluster_colors = {[1 1 1],...        % White
            [1 1 0]...             % Yellow
            [221,52,151]/255,...    % Pink
            [1 0 0],...             % Red
            [250 140 0]/255,...     % oragne
            [210 173 255]/255,...   % purple
            [0.5 0.5 0],...         % Olive
            [0.5 0 0.5],...         % Purple
            [0.5 0 0],...           % Maroon
            [0 1 1],...             % Aqua
            [0 1 0],...             % Lime
            [0 0.5 0.5],...         % Teal
            [0 0.5 0],...           % Green
            [0 0 1],...             % Blue
            [0 0 0.5],...           % Navy
            [0.8 0.8 0.8]};          % Gray
cluster_colors = cluster_colors([4 11 14 2 13 10 5 6 15 16 1 7 9 3]);
% cluster_inds_plot = [3,4,6,7,9,10,11,12];
cluster_titles = {'','','Anterior Cingulate','Right Sensorimotor','Left Occipital','Left Supplementary Motor',...
    'Mid Cingulate','Cuneus','Left Sensorimotor','Left Posterior Parietal',...
    'Right Temporal',...
    'Right Occipital','Right Posterior Parietal'};
dip_dir = [cluster_k_dir filesep 'topo_dip_inf' filesep 'figure_gen'];
%-
% cluster_inds_plot = [3,4,6,7,9,10,11,12]
cluster_inds_plot = [3,4,6,7,8,9,10,13];

%%
%## PLOT ============================================================ %%
fig = figure('color','white','Renderer','painters');
% sgtitle(sprintf('%s',cluster_titles{cl_i}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
p_sz = get(fig,'Position');
set(gca,AXES_DEFAULT_PROPS{:})
hold on;
%## NEXT ROW =========================================================== %%
tmp = openfig([dip_dir filesep 'dipplot_alldipspc_top.fig']);
cnt = length(tmp.Children(end).Children);
while cnt > 3
    % tmp.Children(end).Children(cnt-1).LineWidth = 0.02; %'none';
    % tmp.Children(end).Children(cnt-1).CData = [0,0,0]; %tmp.Children(end).Children(cnt-1).MarkerFaceColor;
    % % fig.Children(end).Children(cnt-1).MarkerFaceColor = DIPPLOT_STRUCT.color{1}; %fig.Children(end).Children(cnt-1+1).Color;
    % tmp.Children(end).Children(cnt-1).MarkerFaceAlpha = 0.85;
    % tmp.Children(end).Children(cnt-1).SizeData = 6; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
    % cnt = cnt - 2;
    tmp.Children(end).Children(cnt).LineWidth = 0.02; %'none';
    tmp.Children(end).Children(cnt).CData = [0,0,0]; %tmp.Children(end).Children(cnt-1).MarkerFaceColor;
    % fig.Children(end).Children(cnt-1).MarkerFaceColor = DIPPLOT_STRUCT.color{1}; %fig.Children(end).Children(cnt-1+1).Color;
    tmp.Children(end).Children(cnt).MarkerFaceAlpha = 0.85;
    tmp.Children(end).Children(cnt).SizeData = 6; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
    cnt = cnt - 1;
end
ax = get(tmp,'CurrentAxes');
view(ax,[90,0])
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');
XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];
dx = abs(ac(1).XData(1,2) - ac(1).XData(1,1));
dy = abs(ac(1).YData(1,1) - ac(1).YData(2,1));
dz = abs(ac(2).ZData(1,1) - ac(2).ZData(2,1));
%## AXES 1
ax1 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax1);
set(ax1,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax1,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax1.XRuler.FirstCrossoverValue = 0;
ax1.YRuler.FirstCrossoverValue = 0;
ax1.ZRuler.FirstCrossoverValue = 0;
view([0,90])
camzoom(ax1,1.1^2)
set(ax1,'YLim',YLIM)
set(ax1,'XLim',XLIM)
drawnow;
pp = get(ax1,'Position');
ax1_x = pp(3)*IM_RESIZE/p_sz(3);
ax1_y = pp(4)*IM_RESIZE/p_sz(4);
% ax1.SortMethod='ChildOrder';
set(ax1,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,ax1_x,ax1_y]);
%## AXES 2
ax2 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax2);
%-
set(ax2,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax2,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax2.XRuler.FirstCrossoverValue = 0;
ax2.YRuler.FirstCrossoverValue = 0;
ax2.ZRuler.FirstCrossoverValue = 0;
set(ax2,'View',[90,0]);
camzoom(ax2,1.1^2)
set(ax2,'YLim',YLIM)
set(ax2,'ZLim',ZLIM)
% set(ax2,'XLim',XLIM)
drawnow;
pp2 = get(ax2,'Position');
% ax2.SortMethod='ChildOrder';
ax2_x = pp2(3)*IM_RESIZE*AX2_ADJ/p_sz(3);
ax2_y = pp2(4)*IM_RESIZE*AX2_ADJ/p_sz(4);
% AX2_W_O = 0.19;
% AX2_H_O = -0.1658;
set(ax2,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+(ax2_x*AX2_X_O),AX_INIT_VERT_DIP+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);

%## AXES 3
ax3 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax3);
%-
set(ax3,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax3,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax3.XRuler.FirstCrossoverValue = 0;
ax3.YRuler.FirstCrossoverValue = 0;
ax3.ZRuler.FirstCrossoverValue = 0;
set(ax3,'view',[0,0])
set(ax3,'ZLim',ZLIM)
set(ax3,'XLim',XLIM)
camzoom(ax3,1.1^2)
pp3 = get(ax3,'Position');

ax3_x = pp3(3)*IM_RESIZE*AX3_ADJ/p_sz(3);
ax3_y = pp3(4)*IM_RESIZE*AX3_ADJ/p_sz(4);
% AX3_W_O = 0.24;
% AX3_H_O = -0.1;
% ax3.SortMethod='ChildOrder';
set(ax3,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),AX_INIT_VERT_DIP+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
close(tmp)
%## FIRST ROW ========================================================== %%
%##
tmp = openfig([dip_dir filesep 'dipplot_avgdipspc_top.fig']);
cnt = length(tmp.Children(end).Children);
surf = tmp.Children(end).Children(1);
while cnt > 3
    hold on;
    tmp.Children(end).Children(cnt).Marker = 'o'; 
    tmp.Children(end).Children(cnt).LineWidth = 0.05; %'none';
    tmp.Children(end).Children(cnt).MarkerFaceAlpha = 0.85;
    tmp.Children(end).Children(cnt).SizeData = 100; %fig.Children(end).Children(cnt+3+1).MarkerSize*2.6;
    cnt = cnt - 1;
end
%## DEBUG
% fig = figure('color','white','Renderer','painters');
% % sgtitle(sprintf('%s',cluster_titles{cl_i}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
% set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
% set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
% set(gca,AXES_DEFAULT_PROPS{:})
%## END DEBUG
% set(tmp,'color','white')
ax = get(tmp,'CurrentAxes');
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');
XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];
dx = abs(ac(1).XData(1,2) - ac(1).XData(1,1));
dy = abs(ac(1).YData(1,1) - ac(1).YData(2,1));
dz = abs(ac(2).ZData(1,1) - ac(2).ZData(2,1));
%-
% IM_RESIZE = 0.25;
%## AXES 1 
ax1 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax1);
%-
set(ax1,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax1,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax1.XRuler.FirstCrossoverValue = 0;
ax1.YRuler.FirstCrossoverValue = 0;
ax1.ZRuler.FirstCrossoverValue = 0;
view([0,90])
camzoom(ax1,1.1^2)
set(ax1,'YLim',YLIM)
set(ax1,'XLim',XLIM)
drawnow;
pp = get(ax1,'Position');
ax1_x = pp(3)*IM_RESIZE/p_sz(3);
ax1_y = pp(4)*IM_RESIZE/p_sz(4);
% ax1.SortMethod='ChildOrder';
set(ax1,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP2,ax1_x,ax1_y]);
%## AXES 2
ax2 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax2);
%-
set(ax2,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax2,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax2.XRuler.FirstCrossoverValue = 0;
ax2.YRuler.FirstCrossoverValue = 0;
ax2.ZRuler.FirstCrossoverValue = 0;
set(ax2,'View',[90,0]);
camzoom(ax2,1.1^2)
set(ax2,'YLim',YLIM)
set(ax2,'ZLim',ZLIM)
% set(ax2,'XLim',XLIM)
drawnow;
pp2 = get(ax2,'Position');
% ax2.SortMethod='ChildOrder';
ax2_x = pp2(3)*IM_RESIZE*AX2_ADJ/p_sz(3); 
ax2_y = pp2(4)*IM_RESIZE*AX2_ADJ/p_sz(4); 
% AX2_W_O = 0.19;
% AX2_H_O = -0.17;
set(ax2,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+(ax2_x*AX2_X_O),AX_INIT_VERT_DIP2+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);
%## AXES 3
ax3 = axes('Parent',fig,'DataAspectRatio',daspect,'Units','normalized');
copyobj(ac,ax3);
%-
set(ax3,'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
% set(ax3,'box','on','xtick',[],'ytick',[],...
%     'ztick',[],'xcolor',[0,0,0],'ycolor',[0,0,0],'zcolor',[0,0,0]);
ax3.XRuler.FirstCrossoverValue = 0;
ax3.YRuler.FirstCrossoverValue = 0;
ax3.ZRuler.FirstCrossoverValue = 0;
set(ax3,'view',[0,0])
set(ax3,'ZLim',ZLIM)
set(ax3,'XLim',XLIM)
camzoom(ax3,1.1^2)
pp3 = get(ax3,'Position');
% AX3_W_ADJ = 1.255;
ax3_x = pp3(3)*IM_RESIZE*AX3_ADJ/p_sz(3);
ax3_y = pp3(4)*IM_RESIZE*AX3_ADJ/p_sz(4);
% AX3_W_O = 0.24;
% AX3_H_O = -0.1;
% ax3.SortMethod='ChildOrder';
set(ax3,'OuterPosition',[0,0,1,1],'Position',[AX_INIT_HORIZ_DIP+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),AX_INIT_VERT_DIP2+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
close(tmp)
%## ANNOTATIONS ======================================================== %%
%-
ANN_H = 0.875;
annotation('textbox',[0.15,ANN_H,.1,.1],...
    'String','Transverse','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
%-
annotation('textbox',[0.43,ANN_H,.1,.1],...
    'String','Sagittal','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
%-
annotation('textbox',[0.745,ANN_H,.1,.1],...
    'String','Coronal','HorizontalAlignment','center',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',12,'FontWeight','Bold','Units','normalized');
%- UPDATE THIS CODE
%## DEBUG
% fig = figure('color','white','Renderer','painters');
% % sgtitle(sprintf('%s',cluster_titles{cl_i}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
% set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
% set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
% p_sz = get(fig,'Position');
% set(gca,AXES_DEFAULT_PROPS{:})
hold on;
%## END DEBUG
bc_shift_x = -0.035;
bc_shift_y = -0.07;
p_sz = get(fig,'Position');
circle_dim = 0.2;
circle_dim_x = circle_dim/(p_sz(3));
circle_dim_y = circle_dim/(p_sz(4));
txt_dim_y = 1/p_sz(3);
txt_dim_x = 2.25/p_sz(4);
txt_shift_x = 0.35/p_sz(4);
txt_shift_y = -0.43/p_sz(3);
HZ_DIM = 3;
AX_INIT_HORIZ = AX_INIT_HORIZ_DIP;
HORIZ_SHIFT_VAL = 13.5/p_sz(4);
VERT_SHIFT_VAL = -2/p_sz(3);
hz = 0;
x_shift = AX_INIT_HORIZ;
y_shift = AX_INIT_VERT_DIP;
aa_bc = [];
aa_txt = [];
for bi = 1:length(cluster_inds_plot)
    hz = hz + 1;
    ANN_H = 0.875;
    tmp = annotation(fig,'ellipse',[bc_shift_x+x_shift,bc_shift_y+y_shift,circle_dim_x,circle_dim_y],...
        'Units','normalized',...
        'FaceColor',cluster_colors{cluster_inds_plot(bi)});
    aa_bc = [aa_bc, tmp];
    tmp = annotation('textbox',[bc_shift_x+x_shift+txt_shift_x,bc_shift_y+y_shift+txt_shift_y,txt_dim_x,txt_dim_y],...
        'String',sprintf('%s',cluster_titles{cluster_inds_plot(bi)}),'HorizontalAlignment','center',...
        'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',12,'FontWeight','Bold','Units','normalized');
    aa_txt = [aa_txt, tmp];
    if hz < HZ_DIM
        x_shift = x_shift + HORIZ_SHIFT_VAL*txt_dim_x;
    else
        y_shift = y_shift + VERT_SHIFT_VAL*txt_dim_y;
        x_shift = AX_INIT_HORIZ;
        hz = 0;
    end    
end
%- edits
aa_bc(end).Position(1) = aa_bc(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
aa_txt(end).Position(1) = aa_txt(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
hold off;
% export_fig([save_dir filesep 'figure_group_dipoleplots_export_fig.pdf'],'-pdf','-nocrop','-painters')
% export_fig([save_dir filesep 'figure_group_dipoleplots_export_fig.pdf'],'pdf')
% exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_1group_comparison_cluster_fooofs_subjs_means.tiff',cl_i,des_i,cond_i,group_i)],'ContentType','image','Resolution',1000)
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.pdf'],'ContentType','vector')
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.tiff'],'ContentType','image','Resolution',600)
% print(fig,'-depsc','-tiff','-painters',[save_dir filesep 'figure_group_dipoleplots_Test_vect.eps'])
% print('-vector','-depsc',[save_dir filesep 'figure_group_dipoleplots_Test_vect.eps']) 
% print(fig,'-vector','-dpdf',[save_dir filesep 'figure_group_dipoleplots_Test_vect.pdf']) 
% print('-vector', '-dsvg', [save_dir filesep 'figure_group_dipoleplots_Test_vect.svg']) 

% close(fig);