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
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(STUDY_DIR);
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
    SRC_DIR = fileparts(STUDY_DIR);
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR);
addpath(SRC_DIR);
addpath(STUDY_DIR);
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
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
% studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load study file
% study_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME)];
%- override
study_fpath = 'R:\Ferris-Lab\jsalminen\Experiments_Data\MIND_IN_MOTION_PRJ\mim_proc_figs_saves\04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [study_fpath filesep 'iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- save directory 
ANALYSIS_DNAME = ['psd_calcs' filesep 'group_spec' filesep 'split_band_test'];
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
%%
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
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
%-- r-stats
r_stats_dir = 'M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\MindInMotion_YoungerOlderAdults_BrainSpeedChanges\src\r_scripts\eeg_speed_lmes';
RSTATS_IMPORT = readtable([r_stats_dir filesep sprintf('02202025_lme_eeg_kin_speed_manu_tests_stats.xlsx')], ...
    "FileType","spreadsheet","UseExcel",true);
%--
tmp = load([save_dir filesep 'fooof_pcond.mat']);
pcond = tmp.pcond;
tmp = load([save_dir filesep 'org_pcond.mat']);
pcond_org = tmp.pcond_org;
tmp = load([save_dir filesep 'fooof_diff_store.mat']);
fooof_diff_store = tmp.fooof_diff_store;
tmp = load([save_dir filesep 'fooof_apfit_store.mat']);
fooof_apfit_store = tmp.fooof_apfit_store;
tmp = load([save_dir filesep 'spec_data_original.mat']);
spec_data_original = tmp.spec_data_original;
tmp = load([save_dir filesep 'fooof_results.mat']);
fooof_results = tmp.fooof_results;
fooof_freq = fooof_results{1}{1,1}{1}.freqs;

%% ===================================================================== %%
%-
des_i = 2;
cl_n = 3;
s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_n).sets).subject};
desdes = cat(1,STUDY.design.variable);
% c_chars = desdes(strcmp({desdes.label},'cond'));
% c_chars = {c_chars.value};
g_chars = desdes(strcmp({desdes.label},'group'));
g_chars = {g_chars.value};
g_chars = g_chars{1};
%--
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
%-- rb3
cluster_titles = {'Left Sensorimotor','Right Posterior Parietal','Mid Cingulate', ...
    'Right Cuneus','Right Sensorimotor','Left Supplementary Motor','Right Occipital', ...
    'Left Occipital','Left Temporal','Left Posterior Parietal','Right Temporal'};
%##
%--
dip_dir = [cluster_k_dir filesep 'topo_dip_inf'];
%--
cluster_inds_plot = [3,4,5,6,7,8,12];
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
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[], ...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
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
hold on;
%- UPDATE THIS CODE
%## DEBUG
% fig = figure('color','white','Renderer','painters');
% % sgtitle(sprintf('%s',cluster_titles{cl_i}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
% set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
% set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
% p_sz = get(fig,'Position');
% set(gca,AXES_DEFAULT_PROPS{:})
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
for cl_i = 1:length(cluster_inds_plot)
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};
    %--
    hz = hz + 1;
    ANN_H = 0.875;
    tmp = annotation(fig,'ellipse',[bc_shift_x+x_shift,bc_shift_y+y_shift,circle_dim_x,circle_dim_y],...
        'Units','normalized',...
        'FaceColor',cluster_colors{cluster_inds_plot(cl_i)});
    aa_bc = [aa_bc, tmp];
    tmp = annotation('textbox',[bc_shift_x+x_shift+txt_shift_x,bc_shift_y+y_shift+txt_shift_y,txt_dim_x,txt_dim_y],...
        'String',sprintf('%s',atlas_name), ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',12, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    aa_txt = [aa_txt, tmp];
    if hz < HZ_DIM
        x_shift = x_shift + HORIZ_SHIFT_VAL*txt_dim_x;
    else
        y_shift = y_shift + VERT_SHIFT_VAL*txt_dim_y;
        x_shift = AX_INIT_HORIZ;
        hz = 0;
    end    
end
%-- edits (03/03/2025)
aa_bc(end).Position(1) = aa_bc(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
aa_txt(end).Position(1) = aa_txt(end).Position(1) + HORIZ_SHIFT_VAL*txt_dim_x;
hold off;
%## SAVE
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.pdf'],'ContentType','vector')
exportgraphics(fig,[save_dir filesep 'figure_group_dipoleplots.tiff'],'ContentType','image','Resolution',600)
