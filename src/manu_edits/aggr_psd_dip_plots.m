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
% r_stats_fn = sprintf('02202025_lme_eeg_kin_speed_manu_tests_stats_oldstudy.xlsx');
r_stats_fn = sprintf('02202025_lme_eeg_kin_speed_manu_tests_stats.xlsx');
RSTATS_IMPORT = readtable([r_stats_dir filesep r_stats_fn], ...
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
%## SPEED MANUSCRIPT GROUP PLOT
designs = unique(fooof_table.design_id);
clusters = unique(fooof_table.cluster_id);
groups = unique(fooof_table.group_id);
EEG_MEASURES = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
PLOT_TITLES = {'Mean \theta','Mean \alpha','Mean \beta'};
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
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
% TOPO_FONTSIZE = 8;
FONT_SIZE_PSD = 9;
FONT_SIZE_PSD_LEG = 9;
% PSD_GROUP_TITLE_FONTSIZE = 12;
PSD_GROUP_TITLE_FONTSIZE = 10;
LEG_TOKEN_SIZE_PSD = 20;
%## VIOLIN PLOT STRUCT
VIO_PLOT_STRUCT = struct('color_map',[],...
    'cond_labels',{{}},...
    'cond_offsets',[-0.35,-0.10,0.15,0.4],...
    'do_group_labels',true, ...
    'group_labels',{{}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.285,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label',{''},...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label',{''},...
    'x_label_fontsize',9,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.1,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',10,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',75);
%## VIOLIN PLOTS
VIOLIN_STRUCT = struct('Width',0.15,...
    'ShowWhiskers',false,...
    'ShowNotches',false,...
    'ShowBox',true,...
    'ShowMedian',true,...
    'Bandwidth',0.1,...
    'QuartileStyle','shadow',...
    'HalfViolin','full',...
    'DataStyle','scatter',...
    'MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],...
    'ViolinAlpha',{{0.2 0.3}},...
    'do_plot_outlier_marks',true,...
    'use_raw_bandwidth',false);
STATS_STRUCT = struct('anova',{{}},...
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
BRACKET_STRUCT = struct('sig_sign','+',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
SIGLINE_STRUCT = struct('sig_sign','*',...
    'line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_offset_x',0,...
    'sig_offset_y',0); 
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

%- rb5
cluster_titles = {'Left Sensorimotor','Right Posterior Parietal','Mid Cingulate', ...
    'Right Cuneus','Right Sensorimotor','Left Supplementary Motor','Right Occipital', ...
    'Left Occipital','Left Temporal','Left Posterior Parietal','Right Temporal'};
out = cellfun(@(x) regexp(x,'(.).*\s(...)','tokens'),cluster_titles);
output_titles = cellfun(@(x) strjoin(x,''),out,'UniformOutput',false);
fig_n = 1:length(cluster_titles);

%##
violin_ylimits = {{[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]};...
            {[],[],[],[],[],[],[],[],[],[],[],[]}};
%--
speed_xvals = (0:5)*0.25;
c_chars = {'0.25 m/s','0.50 m/s','0.75 m/s','1.0 m/s'};
g_chars_topo = {'Young Adults','Older High Functioning Adults','Older Low Functioning Adults'};
g_chars_subp = {'YA','OHFA','OLFA'};
dip_dir = [cluster_k_dir filesep 'topo_dip_inf'];
%--
% cluster_inds_plot = [3,4,5,6,7,8,12];
cluster_inds_plot = (3:13);
%% ===================================================================== %%
%## LEFT SENSORIMOTOR
% cl_i = 1;
des_i = 2;
for cl_i = 1:length(cluster_inds_plot)
    %%
    %## TOPO & DIPOLE PLOTS
    %-- initiate params
    cl_ii = find(cluster_inds_plot(cl_i) == double(string(clusters)));
    cl_n = double(string(clusters(cl_ii)));
    atlas_name = cluster_titles{cl_ii};
    
    %## INITIATE FIGURE
    fig = figure('color','white', ...
        'Renderer','Painters');
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    set(fig,'Units','inches', ...
        'Position',FIGURE_POSITION, ...
        'PaperUnits','inches', ...
        'PaperSize',[1 1], ...
        'PaperPosition',[0 0 1 1])
    p_sz = get(fig,'Position');
    set(gca,AXES_DEFAULT_PROPS{:});
    hold on;
    
    %## TOPOGRAPHY PLOT
    IM_RESIZE = 0.225;
    ax_position = [AX_INIT_HORIZ_TOPO,AX_INIT_VERT_TOPO,0,0];
    plot_topography(fig,STUDY,cl_n, ...
        groups,g_chars,g_chars_topo, ...
        ax_position,IM_RESIZE)
    
    %## DIPOLE PLOT
    IM_RESIZE = 1.1;
    dip_fig_path = [dip_dir filesep sprintf('%i_dipplot_alldipspc_top.fig',cl_n)];
    ax_position = [AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,0,0];
    label_position = [AX_INIT_HORIZ+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1];
    plot_dipole_slices(fig,dip_fig_path,IM_RESIZE, ...
        p_sz,ax_position,label_position)
    
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
        for i = 1:size(spec_data_original{des_i}{cl_n},1)
            tmp_data{cnt} = spec_data_original{des_i}{cl_n}{i,j}';
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
        for i = 1:size(spec_data_original{des_i}{cl_n},1)
            data = spec_data_original{des_i}{cl_n}{i,j}';
            [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));
            Pa.EdgeColor = "none";
        end
        axs = [];
        for i = 1:size(spec_data_original{des_i}{cl_n},1)
            data = spec_data_original{des_i}{cl_n}{i,j}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),...
                'linewidth',LINE_WIDTH_PSD,'LineStyle','-','displayname',sprintf('%s',xtick_label_g{i}));
            set(ax,'Color',[color_dark(i,:),LINE_ALPHA_PSD]);
            axs = [axs, ax];
        end 
        %- plot the aperiodic line
        dash = [];
        for i = 1:size(spec_data_original{des_i}{cl_n},1)
            aperiodic_fit = fooof_apfit_store{des_i}{cl_n}{i,j}';
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
        [axsignif,Pa] = highlight_CL(gca, fooof_freq, pcond_org{des_i}{cl_n}{1}(:,2), 'background', 'Frequency(Hz)');
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
        xlabel('');
        if j == 1
            ylabel('10*log_{10}(PSD)','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSD);
        else
            ylabel('','FontWeight','bold','FontSize',YLABEL_FONT_SIZE_PSD);
        end
        %## TITLE
        title([sprintf('%s',string(g_chars_subp(j))),' PSD'],...
            'FontSize',PSD_GROUP_TITLE_FONTSIZE,'FontWeight','Bold');
        set(gca,'Position',[AX_INIT_HORIZ+x_shift,AX_INIT_VERT_PSD+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_B_XOFFSET+(0.1/2),AX_INIT_VERT_PSD+LAB_B_YOFFSET+(0.1/2),.1,.1],...
        'String','B)', ...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    
    %## PLOT FOOOF PSDS (SPEED) ===================================== %%
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
        for i = 1:size(fooof_diff_store{des_i}{cl_n},1)
            tmp_data{cnt} = fooof_diff_store{des_i}{cl_n}{i,j}';
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
                xtick_label_g = c_chars;
        end
        %- PLOT
        axes();
        hold on;
        axs = [];
        %- standard error shading
        for i = 1:size(fooof_diff_store{des_i}{cl_n},1)
            data = fooof_diff_store{des_i}{cl_n}{i,j}';
            std_error = std(data)/sqrt(size(data,1));
            [Pa,Li] = JackKnife_sung(fooof_freq,mean(data),mean(data)-std_error,mean(data)+std_error,...
                color_dark(i,:),color_light(i,:));
            Pa.EdgeColor = "none";
        end
        for i = 1:size(fooof_diff_store{des_i}{cl_n},1)
            data = fooof_diff_store{des_i}{cl_n}{i,j}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),...
                'linewidth',LINE_WIDTH_PSDFF, ...
                'LineStyle','-', ...
                'displayname',sprintf('%s',xtick_label_g{i}));
            set(ax,'Color',[color_dark(i,:),LINE_ALPHA_PSDFF]);
            axs = [axs, ax];
        end
        %- set limits before sig. shading
        ylim(Y_LIM_GROUP)
        xlim(XLIM_PSD);
        ax = gca;
        [axsignif,Pa] = highlight_CL(gca, fooof_freq, pcond{des_i}{cl_n}{j}(:,2), ...
            'background','Frequency(Hz)');
        plot(XLIM_PSD,[0 0],'--','color','black');     
        for xx = 1:length(XFREQ_LINES)
            xline(XFREQ_LINES(xx),'--');
        end
        set(gca,'FontName',FONT_NAME, ...
            'FontSize',FONT_SIZE_PSD,...
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
        xlabel('Frequency(Hz)', ...
            'FontWeight','bold', ...
            'FontSize',XLABEL_FONT_SIZE_PSDFF);
        if j == 1
            ylabel('10*log_{10}(PSD) - AP. Fit', ...
                'FontWeight','bold', ...
                'FontSize',YLABEL_FONT_SIZE_PSDFF);
        else
            ylabel('','FontWeight','bold', ...
                'FontSize',YLABEL_FONT_SIZE_PSDFF);
        end
        %## TITLE
        % title('Flattened PSD','FontSize',10)
        title([sprintf('%s',string(g_chars_subp(j))),' Flattened PSD'], ...
            'FontSize',PSD_GROUP_TITLE_FONTSIZE, ...
            'FontWeight','Bold');
        set(gca,'OuterPosition',[0,0,1,1]);
        set(gca,'Position',[AX_INIT_HORIZ+x_shift,AX_INIT_VERT_PSDFF+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);
        x_shift = x_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    end
    %## LETTER
    annotation('textbox',[AX_INIT_HORIZ+LAB_C_XOFFSET+(0.1/2),AX_INIT_VERT_PSDFF+LAB_C_YOFFSET+(0.1/2),.1,.1],...
        'String','C)', ...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top', ...
        'LineStyle','none', ...
        'FontName',FONT_NAME,...
        'FontSize',14, ...
        'FontWeight','Bold', ...
        'Units','normalized');
    hold on;

    %## VIOLIN PLOTS) ================================================== %%
    IM_RESIZE = 0.9;
    AX_H  = 0.22;
    AX_W = 0.28;
    AX_HORIZ_SHIFT = 0.06;
    g_chars_vio = {'Y','OHF','OLF'};
    % VIO_PLOT_STRUCT.regresslab_txt_size = 9;
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
    cnt = 1;
    y_shift = 0;    
    x_shift= 0;
    prc_ylim = zeros(length(EEG_MEASURES),2);
    inds = fooof_table.design_id == num2str(des_i) & fooof_table.cluster_id == num2str(cl_n);
    tmp_fooof_t = fooof_table(inds,:); 
    for m_i = 1:length(EEG_MEASURES)
        if isempty(violin_ylimits{m_i}{cl_ii})
            prc_ylim(m_i,:) = [floor(prctile(tmp_fooof_t.(EEG_MEASURES{m_i}),1))-floor(std(tmp_fooof_t.(EEG_MEASURES{m_i}))),...
                ceil(prctile(tmp_fooof_t.(EEG_MEASURES{m_i}),99))+ceil(std(tmp_fooof_t.(EEG_MEASURES{m_i})))*1.5];
            disp(prc_ylim)
        else
            prc_ylim(m_i,:) = violin_ylimits{m_i}{cl_ii};
        end
    end
    for meas_i = 1:length(EEG_MEASURES)
        %## EXTRACT STATS INFO
        COEFF_CHARS_INT = {'(Intercept)','speed_cond_num','group_charH2000''s','group_charH3000''s', ...
            'speed_cond_num:group_charH2000''s','speed_cond_num:group_charH3000''s'};
        ANV_CHARS_INT = {'(Intercept)','speed_cond_num','group_char','speed_cond_num:group_char'};
        ANV_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_char'};
        COEFF_CHARS_GROUP = {'(Intercept)','speed_cond_num','group_charH2000''s','group_charH3000''s'};
        %--
        tmp_stats = RSTATS_IMPORT.cluster_num==double(string(cl_n)) &...
            strcmp(RSTATS_IMPORT.model_char,'speed_group_intact_all') &...
            strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
            strcmp(RSTATS_IMPORT.group_char,'all');
        tmp_stats = RSTATS_IMPORT(tmp_stats,:);
        %--
        tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
        tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
        %-- anova p-values
        anvs = zeros(length(ANV_CHARS_INT),1);
        for cc = 1:length(ANV_CHARS_INT)
            ind = strcmp(ANV_CHARS_INT{cc},tmp_ac);
            if ~isempty(ind)
                anvs(cc) = tmp_anv(ind);              
            else
                fprintf("Coefficient %s not found.\n",ANV_CHARS_INT{cc})
            end            
        end
        
        %## CHECK FOR SIGNIFICANT INTERACTION
        if anvs(4) > 0.05
            %## USE NON-INTERACTION MODEL
            tmp_stats = RSTATS_IMPORT.cluster_num==double(string(cl_n)) &...
                strcmp(RSTATS_IMPORT.model_char,'speed_group_all') &...
                strcmp(RSTATS_IMPORT.freq_band_char,EEG_MEASURES{meas_i}) &...
                strcmp(RSTATS_IMPORT.group_char,'all');
            tmp_stats = RSTATS_IMPORT(tmp_stats,:);
            %--
            tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
            tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
            %-- anova p-values
            anvs = zeros(length(ANV_CHARS_GROUP),1);
            for cc = 1:length(ANV_CHARS_GROUP)
                ind = strcmp(ANV_CHARS_GROUP{cc},tmp_ac);
                if ~isempty(ind)
                    anvs(cc) = tmp_anv(ind);              
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_GROUP{cc})
                end            
            end
            tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
            tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
            tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
            tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
            tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
            % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
            tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
            tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
            %-- model coefficients       
            coeffs = zeros(length(COEFF_CHARS_GROUP),1);        
            for cc = 1:length(COEFF_CHARS_GROUP)
                ind = strcmp(COEFF_CHARS_GROUP{cc},tmp_cc);
                if ~isempty(ind)
                    coeffs(cc) = tmp_coeffs(ind);              
                else
                    fprintf("Coefficient %s not found.\n",COEFF_CHARS_GROUP{cc})
                end            
            end        
            %-- cohens f^2 values
            fsq_chars = strcmp(ANV_CHARS_GROUP,'(Intercept)');
            fsq_chars = ANV_CHARS_GROUP(~fsq_chars);
            fsqs = zeros(length(fsq_chars),1);
            for cc = 1:length(fsq_chars)
                ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
                if ~isempty(ind)
                    fsqs(cc) = tmp_fsq(ind);
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_GROUP{cc})
                end
            end
            %-- confidence intervals
            ci_chars = strcmp(g_chars,'(Intercept)');
            ci_chars = g_chars(~ci_chars);
            cis = zeros(length(ci_chars),1,2);
            for cc = 1:length(ci_chars)
                ind = strcmp(ci_chars{cc},tmp_ci_chars);
                if ~isempty(ind)
                    cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
                else
                    fprintf("Coefficient %s not found.\n",g_chars{cc})
                end
            end
            %--
            % ran_effs_char = strsplit(tmp_stats.ran_effs_char{1},',');
            % ran_effs_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
            %--
            if anvs(3) < 0.05 && anvs(3) > 0.01
                strg = '^{+}';
            elseif anvs(3) <= 0.01 && anvs(3) > 0.001
                strg = '^{++}';
            elseif anvs(3) <= 0.001
                strg = '^{+++}';
            else
                strg = '^{ns}';
            end
            %--
            if anvs(2) < 0.05 && anvs(2) > 0.01
                strs = '^{*}';
            elseif anvs(2) <= 0.01 && anvs(2) > 0.001
                strs = '^{**}';
            elseif anvs(2) <= 0.001
                strs = '^{***}';
            else
                strs = '^{ns}';
            end

            %## ASSIGN STATS
            str = {[sprintf('%sf_{s}^{2}=%1.2f    %sf_{g}^{2}=%1.2f\nR^2=%1.2f', ...
                strs,fsqs(1), ...
                strg,fsqs(2), ...
                tmp_stats.r2_c_int)],'',''};
            %--
            if anvs(3) > 0.05 && anvs(2) > 0.05
                chkd = false;
            else
                chkd = true;
            end
            CONFINT_STRUCT = struct('do_display',chkd, ...
                    'y_bnds',cis, ...
                    'x_vals',repmat(0,[3,1]), ...
                    'errbar_struct',struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
                        'err_bar_width',0.25));
            tssc = struct('var_type','continuous', ...
                'anova',anvs(2), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[[coeffs(1),coeffs(2)]; ...
                    [coeffs(1)+coeffs(3),coeffs(2)]; ...
                    [coeffs(1)+coeffs(4),coeffs(2)]], ... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',(0:5)*0.25,... % continuous predictors
                'order',{{}});
            tssg = struct('var_type','categorical', ...
                'anova',anvs(3), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[],... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',0,... % continuous predictors
                'order',{{}});
            STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
                'cond_stats',tssc, ...
                'group_stats',tssg, ...
                'stats_char',struct('str',{str}, ...
                    'offsets',[-0.1,-0.05], ...
                    'font_size',9, ...
                    'do_display',true));
        else
            %## USE INTERACTION MODEL
            tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
            tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
            tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
            tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
            tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
            % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
            tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
            tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
            %-- model coefficients       
            coeffs = zeros(length(COEFF_CHARS_INT),1);        
            for cc = 1:length(COEFF_CHARS_INT)
                ind = strcmp(COEFF_CHARS_INT{cc},tmp_cc);
                if ~isempty(ind)
                    coeffs(cc) = tmp_coeffs(ind);              
                else
                    fprintf("Coefficient %s not found.\n",COEFF_CHARS_INT{cc})
                end            
            end        
            %-- cohens f^2 values
            fsq_chars = strcmp(ANV_CHARS_INT,'(Intercept)');
            fsq_chars = ANV_CHARS_INT(~fsq_chars);
            fsqs = zeros(length(fsq_chars),1);
            for cc = 1:length(fsq_chars)
                ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
                if ~isempty(ind)
                    fsqs(cc) = tmp_fsq(ind);
                else
                    fprintf("Coefficient %s not found.\n",ANV_CHARS_INT{cc})
                end
            end
            %-- confidence intervals
            ci_chars = strcmp(g_chars,'(Intercept)');
            ci_chars = g_chars(~ci_chars);
            cis = zeros(length(ci_chars),1,2);
            for cc = 1:length(ci_chars)
                ind = strcmp(ci_chars{cc},tmp_ci_chars);
                if ~isempty(ind)
                    cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
                else
                    fprintf("Coefficient %s not found.\n",g_chars{cc})
                end
            end
            %--
            if anvs(4) < 0.05 && anvs(4) > 0.01
                stri = '*';
            elseif anvs(4) <= 0.01 && anvs(4) > 0.001
                stri = '**';
            elseif anvs(4) <= 0.001
                stri = '***';
            else
                stri = '^{ns}';
            end
            %## ASSIGN STATS
            str = {sprintf('%sf_{s:g}^{2}=%1.2f\nR^2=%1.2f', ...
                stri,fsqs(3),tmp_stats.r2_c_int),'',''};
            txt_sz = 9;
            offs = [-0.1,-0.05];
            %--
            % str = {sprintf('%sm_{ya}=%1.2f  m_{ohf}=%1.2f  m_{olf}=%1.2f\nR^2=%1.2f', ...
            %     stri,coeffs(2), ...
            %     coeffs(2)+coeffs(5), ...
            %     coeffs(2)+coeffs(6), ...
            %     tmp_stats.r2_c_int),'',''};
            % txt_sz = 9;
            % offs = [-0.11,-0.05];
            %--
            % str = {sprintf('%sy=(%1.1f)x+(%1.1f)\nR^2=%1.2f',stri,coeffs(2),coeffs(1),tmp_stats.r2_c_int), ...
            %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(5),coeffs(1)+coeffs(3)), ...
            %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(6),coeffs(1)+coeffs(4))};
            % txt_sz = 7;
            % offs = [-0.11,-0.05];
            %--
            if anvs(4) > 0.05 && anvs(3) > 0.05
                chkd = false;
            else
                chkd = true;
            end
            CONFINT_STRUCT = struct('do_display',chkd, ...
                    'y_bnds',cis, ...
                    'x_vals',repmat(0,[3,1]), ...
                    'errbar_struct',struct('line_specs',{{'LineStyle','-','LineWidth',2,'Color','k'}}, ...
                        'err_bar_width',0.25));
            tssc = struct('var_type','continuous', ...
                'anova',anvs(4), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[[coeffs(1),coeffs(2)]; ...
                    [coeffs(1)+coeffs(3),coeffs(2)+coeffs(5)]; ...
                    [coeffs(1)+coeffs(4),coeffs(2)+coeffs(6)]], ...
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',(0:5)*0.25,... % continuous predictors
                'order',{{}});
            tssg = struct('var_type','categorical', ...
                'anova',anvs(4), ...
                'multc_pvals',[],...
                'multc_pvals_pairs',[],...
                'regress_line',[],... 
                'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
                'regress_xvals',0,... % continuous predictors
                'order',{{}});
            STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
                'cond_stats',tssc, ...
                'group_stats',tssg, ...
                'stats_char',struct('str',{str}, ...
                    'offsets',[-0.1,-0.05], ...
                    'font_size',txt_sz, ...
                    'do_display',true));
        end

        %## PLOT
        ax = axes();
        %-- set parameters
        VIO_PLOT_STRUCT.group_labels = g_chars_subp;
        VIO_PLOT_STRUCT.color_map = color_dark;
        VIO_PLOT_STRUCT.cond_labels = xtick_label_g;
        VIO_PLOT_STRUCT.title = PLOT_TITLES(meas_i);
        if meas_i == 1
            VIO_PLOT_STRUCT.y_label ='10*log_{10}(PSD) - AP. Fit';
        else
            VIO_PLOT_STRUCT.y_label ='';
        end
        VIO_PLOT_STRUCT.ylim = prc_ylim(meas_i,:);
        VIO_PLOT_STRUCT.ax_position = [AX_INIT_HORIZ+x_shift,AX_INIT_VERT_VIO+y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
        %-- group violin plot
        ax = group_violin(tmp_fooof_t,EEG_MEASURES{meas_i},'cond_id','group_id',...
            ax,...
            'VIOLIN_STRUCT',VIOLIN_STRUCT,...
            'PLOT_STRUCT',VIO_PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT,...
            'BRACKET_STRUCT',BRACKET_STRUCT,...
            'SIGLINE_STRUCT',SIGLINE_STRUCT, ...
            'CONFINT_STRUCT',CONFINT_STRUCT);
        % ax = group_violin(tmp_fooof_t,EEG_MEASURES{meas_i},'cond_id','group_id',...
        %     ax,...
        %     'VIOLIN_PARAMS',VIOLIN_STRUCT,...
        %     'PLOT_STRUCT',VIO_PLOT_STRUCT,...
        %     'STATS_STRUCT',tmp_stats,...
        %     'BRACKET_STRUCT',BRACKET_STRUCT,...
        %     'SIGLINE_STRUCT',SIGLINE_STRUCT);
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
    exportgraphics(fig,[save_dir filesep sprintf('figure_update_%i_%s.tif',fig_n(cl_ii),output_titles{cl_ii})],'Resolution',600)
    exportgraphics(fig,[save_dir filesep sprintf('figure_update_%i_%s.pdf',fig_n(cl_ii),output_titles{cl_ii})],'ContentType','vector')
    % close(fig);
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
cl_n = 1;
clust_i = double(string(clusters(cl_n)));
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
sgtitle(sprintf('%s',cluster_titles{cl_n}),'FontSize',14,'FontName','Arial','FontWeight','Bold')
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
