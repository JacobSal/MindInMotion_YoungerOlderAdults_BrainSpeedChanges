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
ADD_CLEANING_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
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
%## Hard Defines
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
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- study info
SUB_GROUP_FNAME = 'group_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## LOAD STUDY
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
save_dir = [spec_data_dir filesep 'psd_calcs' filesep 'split_band_test'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%%
%## RE-POP PARAMS
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
FOOOF_TABLE.('log_speed_div_stat') = log10(FOOOF_TABLE.speed_div_stat);
%-
% tmp = load([save_dir filesep 'STATS_TRACK_STRUCT_speedlin.mat']);
% tmp = load([save_dir filesep 'STATS_TRACK_STRUCT_speedquad.mat']);
% psd_feature_stats = tmp.psd_feature_stats;
% psd_feature_stats = tmp.STATS_TRACK_STRUCT;
% psd_feature_stats = struct2table(psd_feature_stats);
%-
tmp = load([save_dir filesep 'fooof_results.mat']);
fooof_results = tmp.fooof_results;
fooof_freq = fooof_results{1}{1,1}{1}.freqs;
%## STATS
iter = 200; % in eeglab, the fdr stats will automatically *20
try
    STUDY.etc = rmfield(STUDY.etc,'statistics');
end
STUDY = pop_statparams(STUDY, 'groupstats','off','condstats', 'on',...
            'method','perm',...
            'singletrials','off','mode','fieldtrip','fieldtripalpha',NaN,...
            'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',iter*20);
stats = STUDY.etc.statistics;
stats.paired{1} = 'on'; % Condition stats
stats.paired{2} = 'off'; % Group stats
%% ===================================================================== %%
%## PARAMS
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
% measure_name_plot = {'theta_avg_power','alpha_avg_power','beta_avg_power','tbr_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% title_plot = {'Mean \theta','Mean \alpha','Mean \beta','\beta//\theta'};
MEASURES_ANALYZE = {'theta_avg_power','beta_avg_power','tbr_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \beta','\beta//\theta'};
%% ===================================================================== %%
if ~ispc
    tmp = strsplit(path,':');
else
    tmp = strsplit(path,';');
end
%*
b1 = regexp(tmp,'fieldtrip','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_fieldtrip = b2{1}(1:b1{1});
    fprintf('fieldtrip path: %s\n',path_fieldtrip);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('fieldtrip path not found.\n');
    end
end
%*
b1 = regexp(tmp,'AAL3','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_aal3 = b2{1}(1:b1{1});
    fprintf('ALL3 path: %s\n',path_aal3);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('AAL3 path not found.\n');
    end
end
%- atlas paths
% ATLAS_FPATHS = {[path_fieldtrip filesep 'template',...
%         filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
%     [path_fieldtrip filesep 'template' filesep 'atlas'...
%         filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
%     [path_fieldtrip filesep 'template' filesep 'atlas'...
%         filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
% ATLAS_FPATHS = {[path_aal3 filesep 'AAL3v1.nii']}; 
ATLAS_FPATHS = {[path_aal3 filesep 'AAL3v1_1mm.nii']}; 
%%
clusters = unique(FOOOF_TABLE.cluster_id);
% [STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
txt_store = cell(length(clusters),1);
atlas_name_store = cell(length(clusters),1);
f = fopen([save_dir filesep 'anatomy_output.txt'],'w');
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
    sprintf('CENTROID: Dip [%0.1f,%0.1f,%0.1f]\n\n',STUDY.cluster(k).centroid.dipole.posxyz)];
    atlas_name_store{k_i} = sprintf('CL%i: %s\n',k,atlas_name);
end
txt_store = txt_store(~cellfun(@isempty,txt_store));
cellfun(@(x) fprintf(f,x),txt_store);
fclose(f);
par_save(atlas_name_store,save_dir,'atlas_names.mat');
%% (LOAD) ============================================================== %%
tmp = load([save_dir filesep 'anatomy_chars.mat']);
atlas_name_store = tmp.atlas_name_store;
%% ===================================================================== %%
% FOOOF_TABLE = FOOOF_TABLE;
%## MEASURES TO COMPUTE (LIMIT TO 3 WHEN POSSIBLE)
MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
MEASURE_SYMBOLS = {'\theta','\alpha','\beta'};
%## (STATS STRUCT) ====================================================== %%
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'measure_tag',categorical({''}),...
    'measure_char',{''},...
    'design_tag',categorical({''}),...
    'design_num',[],...
    'cluster_num',[],...
    'mod_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{''},...
    'anova_preds_p',{[]},...
    'anova_preds_stat',{[]},...
    'anova_preds_df',{[]},...
    'mod_preds_terms',{{''}},...
    'mod_preds_p',[],...
    'mod_preds_stat',[],...
    'mod_preds_coeff',[],...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',[],...
    'multi_comp_p',[],...
    'multi_comp_coeff',[],...
    'multi_comp_lci_uci',[],...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_struct = DEF_STATS_TRACK_STRUCT;
cnts = 1;
%-
TERRAIN_DES_ID = 1;
SPEED_DES_ID = 2;
% MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
% LOG_MEASURE_NAMES = {'log_alpha_avg_power','log_theta_avg_power','log_beta_avg_power'};
%% (STAT) PREDICTORS: SPEED_DIFF_STAT. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
% SPEED_MEAS_ANL = 'speed_diff_stat';
SPEED_MEAS_ANL = 'speed_div_stat';
% SPEED_MEAS_ANL = 'speed_diffdiv_stat';
%-
IM_RESIZE = 1;
TITLE_TXT_SIZE = 14;
%-
AX_FONT_NAME = 'Arial';
AX_W = 0.35;
AX_H = 0.25;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT = 0.65;
g_chars_subp = {'YA','OHMA','OLMA'};

tmp_savedir = [save_dir filesep sprintf('lme_P%s_stat-Reeg_linear',SPEED_MEAS_ANL)];
mkdir(tmp_savedir);
%-
%-
groups = unique(FOOOF_TABLE.group_id);
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for cl_i = 1:length(clusters)
    atlas_name = atlas_name_store{cl_i};
    fig = figure('color','white');
    sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    horiz_shift = 0;
    plot_leg = [];
    % figure;
    %##
    cond_plot_store = [];
    group_plot_store = [];
    %- data
    inds = FOOOF_TABLE.design_id == designs(des_i) & FOOOF_TABLE.cluster_id == clusters(cl_i);
    tmp_data = FOOOF_TABLE(inds,:);
    tmp_data.(SPEED_MEAS_ANL) = log10(tmp_data.(SPEED_MEAS_ANL));
    y_lim_calc = [prctile(tmp_data.(SPEED_MEAS_ANL),3)-std(tmp_data.(SPEED_MEAS_ANL)),...
        prctile(tmp_data.(SPEED_MEAS_ANL),97)+std(tmp_data.(SPEED_MEAS_ANL))];
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
    ax = axes();
    VIOLIN_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
    PLOT_STRUCT = struct('color_map',color_dark,...
        'cond_labels',categorical({'.'}),'group_labels',categorical(g_chars_subp),...
        'cond_offsets',[0],...
        'group_offsets',[0.125,0.475,0.812],...
        'y_label',SPEED_MEAS_ANL,...
        'title','Terrain Treadmill Speed by Group','font_size',10,'ylim',y_lim_calc,...
        'font_name','Arial','x_label',{''},'do_combine_groups',false,...
        'regresslab_txt_size',8);
    axax = group_violin(tmp_data,SPEED_MEAS_ANL,'design_id','group_id',...
        ax,...
        'VIOLIN_PARAMS',VIOLIN_PARAMS,...
        'PLOT_STRUCT',PLOT_STRUCT);
    set(axax,'OuterPosition',[0,0,1,1]);
    set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]

    hold off;
    %##
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_eeg-%s_lin.tiff',string(clusters(cl_i)),SPEED_MEAS_ANL)],'Resolution',300)
    % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_g%i_eeg-speeddiffstat_linetestp.tiff',string(clusters(cl_i)),groups(group_i))],'Resolution',300)
    % close(fig)
end
%% (STAT) PREDICTORS: SPEED_DIFF_STAT. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
g_i = 1;
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% SPEED_MEAS_ANL = 'speed_diff_stat';
% SPEED_MEAS_ANL = 'speed_div_stat';
SPEED_MEAS_ANL = 'log_speed_div_stat';
% SPEED_MEAS_ANL = 'speed_diffdiv_stat';
%-
COLORS = linspecer(2);
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_FONT_NAME = 'Arial';
AX_W = 0.35;
AX_H = 0.25;
AX_HORIZ_SHIFT = 0.1;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT = 0.65;
AX_TXT_SIZE = 10;
%-
LINE_WIDTH_REFF = 2;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = true;
REG_TXT_SIZE = 8;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.125;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep sprintf('lme_P%s_stat-Reeg_linear',SPEED_MEAS_ANL)];
mkdir(tmp_savedir);
%-
%-
groups = unique(FOOOF_TABLE.group_id);
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for cl_i = 1:length(clusters)
    for g_i = 1:length(groups)
        %##
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        horiz_shift = 0;
        plot_leg = [];
        for meas_i = 1:length(MEASURES_ANALYZE)
            tmptmp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == clusters(cl_i) &...
                FOOOF_TABLE.design_id == designs(des_i),:);
            y_lim_calc = [prctile(tmptmp_table.(MEASURES_ANALYZE{meas_i}),3)-std(tmptmp_table.(MEASURES_ANALYZE{meas_i})),...
                        prctile(tmptmp_table.(MEASURES_ANALYZE{meas_i}),97)+std(tmptmp_table.(MEASURES_ANALYZE{meas_i}))];
            % figure;
            %##
            cond_plot_store = [];
            group_plot_store = [];
            %- data
            tmptmp_table = FOOOF_TABLE(FOOOF_TABLE.group_id == groups(g_i) & FOOOF_TABLE.cluster_id == clusters(cl_i) &...
                FOOOF_TABLE.design_id == designs(des_i),:);
            % tmptmp_table = table(categorical(string(tmptmp_table.subj_char)),double(tmptmp_table.(MEASURES_ANALYZE{meas_i})),...
            %     double(string(tmptmp_table.cond_char)),...
            %     categorical(string(tmptmp_table.group_id)),'VariableNames',{'subj_char',MEASURES_ANALYZE{meas_i},'cond_id','group_id'});
            tmptmp_table = table(categorical(string(tmptmp_table.subj_char)),double(tmptmp_table.(MEASURES_ANALYZE{meas_i})),...
                double(string(tmptmp_table.(SPEED_MEAS_ANL))),double(tmptmp_table.(SPEED_MEAS_ANL)),...
                categorical(string(tmptmp_table.group_id)),'VariableNames',{'subj_char',MEASURES_ANALYZE{meas_i},'cond_id',SPEED_MEAS_ANL,'group_id'});
    
            %## LINEAR MODEL
            % mod_lme = sprintf('%s ~ 1 + cond_char + (1|subj_char)',MEASURES_ANALYZE{meas_i});
            % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
            % stats_out = fitlme(tmptmp_table,mod_lme,'FitMethod','REML');
            % anova_out = anova(stats_out);
            mod_lme = sprintf('%s ~ 1 + cond_id + (1|subj_char)',MEASURES_ANALYZE{meas_i});
            % mod_lme = sprintf('%s ~ 1 + cond_id^2 + cond_id + (1|subj_char)',MEASURES_ANALYZE{meas_i});
            % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
            stats_out = fitlme(tmptmp_table,mod_lme,'FitMethod','REML');
            % anova_out = anova(stats_out);
            [p,t,anova_out,terms] = anovan(tmptmp_table.(MEASURES_ANALYZE{meas_i}),{tmptmp_table.cond_id},...
            'sstype',3,'varnames',{'cond_id'},'model','linear','Display','off');
            %## GATHER STATS
            %- test normality
            [norm_h,norm_p] = lillietest(stats_out.residuals);
            %- get effects
            [~,bnames,~] = stats_out.fixedEffects();
            [~,brnames,bretable] = stats_out.randomEffects();
            %- intercept only model
            % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            R2 = stats_out.Rsquared.Adjusted;
            %- intercept only model
            altmod_out = sprintf('%s ~ 1 + (1|subj_char)',MEASURES_ANALYZE{meas_i});
            altstats_out = fitlme(tmptmp_table,altmod_out);
            %- alternative f2?
            R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
            R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
            alt_f2 = (R22-R21)/(1-R22);
            %- populate struct
            stats_struct(cnts).stat_test_mod = mod_lme;
            stats_struct(cnts).measure_tag = categorical({MEASURES_ANALYZE{meas_i}});
            stats_struct(cnts).design_tag = categorical(des_i);
            stats_struct(cnts).mod_tag = categorical(2);
            % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
            stats_struct(cnts).mod_resp_terms = MEASURES_ANALYZE{meas_i};
            stats_struct(cnts).anova_preds_terms = t(:,1)';
            tmp = t(:,7)';
            tmp = tmp(~cellfun(@isempty,tmp));
            stats_struct(cnts).anova_preds_p = tmp;
            tmp = t(:,6)';
            tmp = tmp(~cellfun(@isempty,tmp));
            stats_struct(cnts).anova_preds_stat = tmp;
            tmp = t(:,3)';
            tmp = tmp(~cellfun(@isempty,tmp));
            stats_struct(cnts).anova_preds_df =tmp;
            stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
            stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
            stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
            stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
            stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
            stats_struct(cnts).norm_test_p = norm_p;
            stats_struct(cnts).norm_test_h = norm_h;
            stats_struct(cnts).effect_size = alt_f2;
            stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
            cnts = cnts + 1;
            stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
            %## PLOT =============================================================== %%
            tmp = t(:,7)';
            anova_p_var = tmp{2};
            lme_speed_coeff = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_id')));
            lme_inter_coeff = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            rnd_eff = bretable;
            %## SCATTER
            %- plot
            % figure;
            axes();
            hold on;
            data = tmptmp_table;
            [~,inds] = sort(data.(MEASURES_ANALYZE{meas_i}));
            data = data(inds,:);
            % ss = scatter(data,'speed_diff_stat',MEASURES_ANALYZE{meas_i},'DisplayName',sprintf('%s',cond_chars(cond_i)));
            ss = scatter(data,SPEED_MEAS_ANL,MEASURES_ANALYZE{meas_i},'DisplayName','rel. speed');
            ss.CData = COLORS(1,:);
            ss.SizeData = 15;
            ss.MarkerEdgeAlpha = 0.6;
            %## LINEAR MODEL FIT
            hold on;
            data = tmptmp_table;
            [vals,inds] = sort(data.(SPEED_MEAS_ANL));
            data = data(inds,:);
            if anova_p_var < ALPHA 
                %## PLOT RANDOM EFFECTS
                subjects = unique(rnd_eff.Level);
                rnd_colors = linspecer(length(subjects));
                % cnt = 1;
                inter = '(Intercept)';
                x1 = 'cond_id';
                for subj_i = 1:length(subjects)
                    inds = strcmp(rnd_eff.Level, subjects{subj_i});
                    est = rnd_eff.Estimate(inds);
                    levels = rnd_eff.Name(inds);
                    interv = est(strcmp(levels,inter));
                    x1v = est(strcmp(levels,x1));
                    % x = unique(data.(SPEED_MEAS_ANL));
                    subj_inds = data.subj_char == categorical(subjects(subj_i));
                    x = data(subj_inds,:);
                    x = x.(SPEED_MEAS_ANL);
                    %## LINEAR
                    y = lme_speed_coeff.*x + lme_inter_coeff + interv;
                    pps = plot(x,y,...
                        'DisplayName',sprintf('subj. %s',subjects{subj_i}),...
                        'LineWidth',LINE_WIDTH_REFF);
                    pps.Color = rnd_colors(subj_i,:)*0.8;
                    pps.LineStyle = ':';
                end
                %## PLOT MAIN EFFECTS
                x = data.(SPEED_MEAS_ANL);
                %## LINEAR
                y = lme_speed_coeff.*x + lme_inter_coeff;                
                pp = plot(x,y,...
                    'DisplayName',sprintf('pvalue_{%s}=(%0.2f)','sp',anova_p_var),...
                    'LineWidth',LINE_WIDTH_MEFF);
                pp.Color = COLORS(2,:);
                %##
                if DO_PLOT_R2
                    eq = '';
                    if anova_p_var > 0.01 && anova_p_var < ALPHA
                        annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                            'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                            'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                            'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                            'BackgroundColor','none');
                    elseif anova_p_var <= 0.01 && anova_p_var > 0.001
                        annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                            'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                            'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                            'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                            'BackgroundColor','none');
                    elseif anova_p_var <= 0.001
                        annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                            'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                            'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                            'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                            'BackgroundColor','none');
                    end
                end
                plot_leg = [plot_leg, pp];
            end
            %## AX EDITS
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('')
            end
            % xlabel('\DeltaSpeed_{OG} (%)');
            % xlabel('log_{10}(Speed_{ms}/OG_{ms})')
            % xlabel('$\log_{10}$$\frac{Trial_{ms}}{Subject_{ms}}$','Interpreter','latex');
            xlabel('$\log_{10} ( \frac{Trial_{ms}}{Subject_{ms}} )$','Interpreter','latex','FontSize',14,'FontWeight','bold');
            set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            % xlab = xlabel(['$(log_{10})$','$(\frac{x}{x})$','(%)'],'FontSize',8,'Interpreter','latex');
            title(MEASURE_NAME_LABS{meas_i});
            ylim(y_lim_calc)
            % xlim([-1,max(data.(SPEED_MEAS_ANL))])
            %## legend
            if meas_i == length(MEASURES_ANALYZE)
                %- lg3
                if ~isempty(plot_leg)
                    legend(gca,plot_leg);
                    [lg3,~,~,~] = legend('boxoff');
                    set(lg3,'Orientation','horizontal')
                    set(lg3,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE);
                    set(lg3,'NumColumns',4);
                    % set(lg1,'Position',[0.1,AX_INIT_VERT+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                        AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.05,lg3.Position(3),lg3.Position(4)]);
                    lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
            end
            set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            if mod(meas_i,AX_MAX) == 0
                horiz_shift = 0;
                vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
            end
        end
        hold off;
        %##
        % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_g%i_eeg-speeddiffstat_quadtestp.tiff',string(clusters(cl_i)),groups(group_i))],'Resolution',300)
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_g%i_%s_linetestp.tiff',string(clusters(cl_i)),groups(g_i),SPEED_MEAS_ANL)],'Resolution',300)
        close(fig)
    end
end

%% (LINE PLOT) PREDICTORS: SPEED_DIFF_STAT. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%- model specific
MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
MEASURE_SYMBOLS = {'\theta','\alpha','\beta'};
% SPEED_MEAS_ANL = 'speed_diffdiv_stat';
%- collect
DEF_FIT_STRUCT = struct('subject',categorical({''}),...
                    'cluster',categorical({''}),...
                    'design',categorical({''}),...
                    'group',categorical({''}),...
                    'group_char',categorical({''}),...
                    'measure',categorical({''}),...
                    'measure_char',{{''}},...
                    'b1',0,...
                    'b2',0,...
                    'b0',0,...
                    'R2',0,...
                    'llh',0,...
                    'x_fit',zeros(1,2),...
                    'y_fit',zeros(1,2),...
                    'x_raw',zeros(2,1),...
                    'y_raw',zeros(2,1));
FIT_STRUCT = DEF_FIT_STRUCT;
cnt = 1;
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
% group_i = 1;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
group_chars = unique(FOOOF_TABLE.group_char);
%-
COLORS = linspecer(2);
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
% ALPHA = 0.05;
% DO_PLOT_RANDOM_VARS = false;
%-
AX_FONT_NAME = 'Arial';
AX_W = 0.35;
AX_H = 0.25;
AX_HORIZ_SHIFT = 0.08;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.09;
AX_INIT_VERT = 0.75;
AX_TXT_SIZE = 10;
%-
LINE_WIDTH_REFF = 1;
% LINE_WIDTH_MEFF = 3;
%-
% DO_PLOT_R2 = false;
% REG_TXT_SIZE = 8;
% REG_HORIZ_SHIFT = 0.05;
% REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.125;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep sprintf('lmquad_P%s_nostat-Reeg_linear',SPEED_MEAS_ANL)];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for g_i = 1:length(groups)
        
        fig = figure('color','white');
        
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        horiz_shift = 0;
        plot_leg = [];
        for meas_i = 1:length(MEASURES_ANALYZE)
            %##
            cond_plot_store = [];
            group_plot_store = [];
            %- get flat eeg power
            inds = (FOOOF_TABLE.cluster_id == clusters(cl_i) &...
                FOOOF_TABLE.group_id == groups(g_i) & FOOOF_TABLE.cond_char == 'flat');
            flat_vals = FOOOF_TABLE(inds,:);
            %- data
            if ~isempty(g_i)
                inds = FOOOF_TABLE.design_id == designs(des_i) & FOOOF_TABLE.cluster_id == clusters(cl_i)...
                    & FOOOF_TABLE.group_id == groups(g_i);
            else
                inds = FOOOF_TABLE.design_id == designs(des_i) & FOOOF_TABLE.cluster_id == clusters(cl_i);
            end
            tmp_data = FOOOF_TABLE(inds,:);
            %## TITLE
            subjects = unique(tmp_data.subj_id);
            atlas_name = sprintf('%s, N=%i) %s',group_chars(g_i),length(subjects),atlas_name_store{cl_i});
            sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
            %## NORMALIZE MEASURE
            for subj_i = 1:length(subjects)
                %- grab subject and its data
                indsd = tmp_data.subj_id == subjects(subj_i);
                datad = tmp_data(indsd,:);
                %- grab subjects data on flat terrain
                inds_flat = flat_vals.subj_id == subjects(subj_i);
                data_flat = flat_vals(inds_flat,:);
                %- get eeg data and normalize to flat's eeg data.
                indv = strcmp(tmp_data.Properties.VariableNames,MEASURES_ANALYZE{meas_i});
                % tmp_data{indsd,indv} = datad.(MEASURES_ANALYZE{meas_i}) - data_flat.(MEASURES_ANALYZE{meas_i});
                tmp_data{indsd,indv} = datad.(MEASURES_ANALYZE{meas_i}) - data_flat.(MEASURES_ANALYZE{meas_i});
                %- check if overground speed matches flat speed & if eeg is
                %same
                chk = any(data_flat.speed_ms == double(string(datad.cond_char)));
                if ~any(chk)
                    %- grab eeg data and subject's speed normalized to
                    %overground
                    indvs = strcmp(tmp_data.Properties.VariableNames,SPEED_MEAS_ANL);
                    tmp_data(end+1,:) = datad(1,:);
                    tmp_data{end,indv} = data_flat.(MEASURES_ANALYZE{meas_i}) - data_flat.(MEASURES_ANALYZE{meas_i});
                    tmp_data{end,indvs} = 0;
                end                
            end
            %## ADD FLAT
            % for subj_i = 1:length(subjects)
            %     %- grab subject and its data
            %     indsd = tmp_data.subj_id == subjects(subj_i);
            %     datad = tmp_data(indsd,:);
            %     %- grab subjects data on flat terrain
            %     inds_flat = flat_vals.subj_id == subjects(subj_i);
            %     data_flat = flat_vals(inds_flat,:);
            %     %- get eeg data.
            %     indv = strcmp(tmp_data.Properties.VariableNames,MEASURES_ANALYZE{meas_i});
            %     %- check if overground speed matches flat speed & if eeg is
            %     %same
            %     chk = any(data_flat.speed_ms == double(string(datad.cond_char)));
            %     if ~any(chk)
            %         %- grab FLAT eeg data and subject's
            %         indvs = strcmp(tmp_data.Properties.VariableNames,SPEED_MEAS_ANL);
            %         tmp_data(end+1,:) = datad(1,:);
            %         tmp_data{end,indv} = data_flat.(MEASURES_ANALYZE{meas_i});
            %         tmp_data{end,indvs} = 0;
            %     end             
            % end
            %-
            y_lim_calc = [min(tmp_data.(MEASURES_ANALYZE{meas_i}))-std(tmp_data.(MEASURES_ANALYZE{meas_i})),...
                max(tmp_data.(MEASURES_ANALYZE{meas_i}))+std(tmp_data.(MEASURES_ANALYZE{meas_i}))];
            %## SCATTER
            %- plot
            axes();
            hold on;
            data = tmp_data;
            ss = scatter(data,SPEED_MEAS_ANL,MEASURES_ANALYZE{meas_i},'DisplayName','rel. speed');
            ss.CData = COLORS(1,:);
            ss.SizeData = 15;
            ss.MarkerEdgeAlpha = 0.6;
            
            %## LINEAR MODEL FIT
            hold on;
            % data = tmp_data;
            % [vals,inds] = sort(data.(SPEED_MEAS_ANL));
            % data = data(inds,:);
            % rnd_eff = tmp_stats.lme_rnd_effects{1}{1};
            subjects = unique(tmp_data.subj_id);
            rnd_colors = linspecer(size(subjects,1));
            for subj_i = 1:length(subjects)
                inds = tmp_data.subj_id == subjects(subj_i);
                data = tmp_data(inds,:);
                %- raw data
                % x = data.(SPEED_MEAS_ANL);
                % y = data.(MEASURES_ANALYZE{meas_i});
                % pps = plot(x,y,...
                %     'DisplayName',sprintf('subj. %s',subjects(subj_i)),...
                %     'LineWidth',LINE_WIDTH_REFF);
                % % pp.Color = rnd_colors(subj_i)*0.8;
                % pps.Color = [COLORS(2,:),0.5];
                % pps.LineStyle = '-';
                %- fit quad
                mod_lm = sprintf('%s ~ 1 + %s + %s^2',MEASURES_ANALYZE{meas_i},SPEED_MEAS_ANL,SPEED_MEAS_ANL);
                lm_out = fitlm(data,mod_lm);
                b0 = lm_out.Coefficients.Estimate(strcmp(lm_out.Coefficients.Properties.RowNames,'(Intercept)'));
                b1 = lm_out.Coefficients.Estimate(strcmp(lm_out.Coefficients.Properties.RowNames,SPEED_MEAS_ANL));
                b2 = lm_out.Coefficients.Estimate(strcmp(lm_out.Coefficients.Properties.RowNames,sprintf('%s^2',SPEED_MEAS_ANL)));
                R2 = lm_out.Rsquared.Adjusted;
                llh = lm_out.LogLikelihood;
                %- plot model fit with raw x values
                % x = sort(data.(SPEED_MEAS_ANL));
                % y = b1*x + b2*(x.^2) + b0;
                %- plot model fit with smooth x values
                x = data.(SPEED_MEAS_ANL);
                x = linspace(min(x),max(x),50);
                y = b1*x + b2*(x.^2) + b0;
                pps = plot(x,y,...
                    'DisplayName',sprintf('subj. %s',subjects(subj_i)),...
                    'LineWidth',LINE_WIDTH_REFF);
                pps.Color = [rnd_colors(subj_i,:),0.5];
                % pps.Color = [COLORS(2,:),0.5];
                pps.LineStyle = '-';
                %- store
                FIT_STRUCT(cnt).subject = subjects(subj_i);
                FIT_STRUCT(cnt).cluster = clusters(cl_i);
                FIT_STRUCT(cnt).design = designs(des_i);
                FIT_STRUCT(cnt).group = groups(g_i);
                FIT_STRUCT(cnt).group_char = group_chars(g_i);
                FIT_STRUCT(cnt).measure_char = MEASURES_ANALYZE{meas_i};
                FIT_STRUCT(cnt).measure = categorical(meas_i);
                FIT_STRUCT(cnt).b1 = b1;
                FIT_STRUCT(cnt).b2 = b2;
                FIT_STRUCT(cnt).b0 = b0;
                FIT_STRUCT(cnt).R2 = R2;
                FIT_STRUCT(cnt).llh = llh;
                FIT_STRUCT(cnt).x_raw = data.(SPEED_MEAS_ANL);
                FIT_STRUCT(cnt).y_raw = data.(MEASURES_ANALYZE{meas_i});
                FIT_STRUCT(cnt).x_fit = x;
                FIT_STRUCT(cnt).y_fit = y;
                FIT_STRUCT(cnt+1) = DEF_FIT_STRUCT;
                cnt = cnt + 1;
            end
            %## AX EDITS
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('')
            end
            xlabel('\DeltaSpeed_{OG} (%)');
            title(MEASURE_NAME_LABS{meas_i});
            set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            ylim(y_lim_calc)
            % xlim([-1,max(data.(SPEED_MEAS_ANL))])
            %## legend
            if meas_i == length(MEASURES_ANALYZE)
                %- lg3
                if ~isempty(plot_leg)
                    legend(gca,plot_leg);
                    [lg3,~,~,~] = legend('boxoff');
                    set(lg3,'Orientation','horizontal')
                    set(lg3,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE);
                    set(lg3,'NumColumns',4);
                    % set(lg1,'Position',[0.1,AX_INIT_VERT+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                        AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.05,lg3.Position(3),lg3.Position(4)]);
                    lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
            end
            set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            if mod(meas_i,AX_MAX) == 0
                horiz_shift = 0;
                vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_eeg-speeddiffstat_smooth_normal.tiff',string(clusters(cl_i)),string(group_chars(g_i)))],'Resolution',450)
        
        % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_eeg-speeddiffstat_smooth_nonormal.tiff',string(clusters(cl_i)),string(group_chars(group_i)))],'Resolution',450)
        close(fig)
    end
end
FIT_STRUCT = struct2table(FIT_STRUCT);
FIT_STRUCT(isundefined(FIT_STRUCT.group),:) = [];
save([tmp_savedir filesep 'FIT_STRUCT_nonormal.m','FIT_STRUCT'])

%% (MAIN PLOTTING LOOP)
%## PARAMS
%- model specific
MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
MEASURE_SYMBOLS = {'\theta','\alpha','\beta'};
% SPEED_MEAS_ANL = 'speed_diffdiv_stat';
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 1;
% group_i = 1;
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
group_chars = unique(FOOOF_TABLE.group_char);
%-
COLORS = linspecer(2);
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
ALPHA = 0.05;
%-
AX_FONT_NAME = 'Arial';
AX_W = 0.35;
AX_H = 0.25;
AX_HORIZ_SHIFT = 0.08;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.09;
AX_INIT_VERT = 0.65;
AX_TXT_SIZE = 10;
%-
LINE_WIDTH_REFF = 1;
REG_TXT_SIZE = 25;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.125;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep sprintf('lmquad_P%s_meanstat_manova-Reeg_linear',SPEED_MEAS_ANL)];
mkdir(tmp_savedir);
%-
groups = unique(FIT_STRUCT.group);
designs = unique(FIT_STRUCT.design);
clusters = unique(FIT_STRUCT.cluster);
measures = unique(FIT_STRUCT.measure);
for cl_i = 1:length(clusters)
    %%
    fig = figure('color','white');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    horiz_shift = 0;
    plot_leg = [];
    for meas_i = 1:length(MEASURES_ANALYZE)
        %## IMPLEMENTATION
        % GROUPINGS = {1,[2,3]};
        GROUPINGS = {1,2,3};
        groups_plot = 1:length(GROUPINGS);
        inds = FIT_STRUCT.cluster == clusters(cl_i) &...
            FIT_STRUCT.design == designs(des_i) &...
            FIT_STRUCT.measure == measures(meas_i);
        tmp_table = FIT_STRUCT(inds,:);
        %## TITLE
        subjects = unique(tmp_table.subject);
        atlas_name = sprintf('N=%i) %s',length(subjects),atlas_name_store{cl_i});
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        %## MANOVA
        mod_in = 'b0,b1,b2 ~ 1+group';
        % mod_in = 'b1,b2 ~ 1+group';
        man_out = manova(tmp_table,mod_in);
        stats_out = man_out.stats;
        chk = stats_out.pValue < ALPHA;
        if any(chk)
            anova_p_var = stats_out.pValue(chk);
        else
            anova_p_var = 1;
        end
        % groupmeans(man_out,"b1")
        %## CHANGE GROUPINGS?
        if length(GROUPINGS) ~= length(groups)
            indsn = zeros(size(tmp_table,1),1);
            for g_i = 1:length(GROUPINGS)
                g1 = GROUPINGS{g_i};
                for g_j = 1:length(g1)
                    indsg = tmp_table.group == groups(g1(g_j));
                    % indsn(indsg) = g_i;
                    tmp_table.group(indsg,:) = categorical(g_i);
                end
            end
        end
        %## PLOT
        xname = 'x_raw';
        yname = 'y_raw';
        colors = linspecer(length(groups_plot));
        colors_mean = colors*0.8;
        axes();
        hold on;
        % plot_leg = [];
        for g_i = 1:length(groups_plot)
            indsg = tmp_table.group == categorical(groups_plot(g_i));
            tmp = tmp_table(indsg,:);
            xls = cellfun(@length,tmp.(xname));
            xn = zeros(length(tmp.(xname)),max(xls));
            x = cat(1,tmp.(xname){:});
            xu = unique(x);
            % yls = cellfun(@length,tmp.(yname));
            % yls = max(yls);
            %-
            % xn = cat(1,tmp.(xname){:});
            yn = nan(length(tmp.(yname)),length(xu));
            for subj_i = 1:length(tmp.(yname))
                yn(subj_i,:) = tmp.b0(subj_i) + tmp.b1(subj_i)*xu + tmp.b2(subj_i)*xu.^2;
                % tmpz = zeros(size(xu));
                % [~,xinds] = setdiff(xu,tmp.(xname){subj_i});
                % tmpz(xinds) = 1;
                % yn(subj_i,~logical(tmpz)) = tmp.(xname){subj_i};
            end
            % figure;
            xm = xu';
            ym = mean(yn,1);
            ys = std(yn,[],1);
            %-
            xf = [xm, fliplr(xm)];
            yf = [(ym + ys), fliplr((ym - ys))];
            ff = fill(xf, yf,'g');
            hold on;
            ff.FaceColor = colors(g_i,:);
            ff.FaceAlpha = 0.3;
            ff.EdgeColor = 'none';
            pp = plot(xm,ym,'Color',colors_mean(g_i,:),'DisplayName',sprintf('Group: %s',string(group_chars(g_i))),...
                'LineWidth',3);
            if meas_i == 1
                plot_leg = [plot_leg, pp];
            end
        end
        xline(gca,0,'k--')
        yline(gca,0,'k--')
        if anova_p_var < ALPHA
            eq = '';
            if anova_p_var > 0.01 && anova_p_var < ALPHA
                annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                    'String','*','HorizontalAlignment','center',...
                    'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                    'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                    'BackgroundColor','none');
            elseif anova_p_var <= 0.01 && anova_p_var > 0.001
                annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                    'String','**','HorizontalAlignment','center',...
                    'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                    'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                    'BackgroundColor','none');
            elseif anova_p_var <= 0.001
                annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                    'String','***','HorizontalAlignment','center',...
                    'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                    'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                    'BackgroundColor','none');
            end
        end
        %## AX EDITS
        if meas_i == 1
            ylabel('10*log_{10}(Flattened PSD)');
        else
            ylabel('')
        end
        xlabel('\DeltaSpeed_{OG} (%)');
        % xlim([-1,max(data.(SPEED_MEAS_ANL))])
        ylim([-5,5]);
        title(MEASURE_NAME_LABS{meas_i});
        set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
        %## LEGEND
        if meas_i == length(MEASURES_ANALYZE)
            %- lg3
            if ~isempty(plot_leg)
                legend(gca,plot_leg);
                [lg3,~,~,~] = legend('boxoff');
                set(lg3,'Orientation','horizontal')
                set(lg3,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE);
                set(lg3,'NumColumns',4);
                % set(lg1,'Position',[0.1,AX_INIT_VERT+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                    AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.05,lg3.Position(3),lg3.Position(4)]);
                lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
            end
        end
        set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        if mod(meas_i,AX_MAX) == 0
            horiz_shift = 0;
            vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
        end
    end
    hold off;
    %##
    % exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_eeg-speeddiffstat_smooth_nonormal.tiff',string(clusters(cl_i)))],'Resolution',450)
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_eeg-speeddiffstat_smooth_normal.tiff',string(clusters(cl_i)),SPEED_MEAS_ANL)],'Resolution',450)

    close(fig)
end
%% (CURVE STATS) Individual Fits
%{
%## CONCEPT
%{
figure;
hold on;
b1 = 2;
b2 = 2;
b0 = 0;
x = linspace(-1,1,50);
y = b1*x + b2*(x.^2) + b0;
plot(x,y,'DisplayName','+b2,+b1');
b1 = -2;
b2 = 2;
b0 = 0;
x = linspace(-1,1,50);
y = b1*x + b2*(x.^2) + b0;
plot(x,y,'DisplayName','+b2,-b1');
b1 = 2;
b2 = -2;
b0 = 0;
x = linspace(-1,1,50);
y = b1*x + b2*(x.^2) + b0;
plot(x,y,'DisplayName','-b2,+b1');
hold off;
legend();
%}

%## IMPLEMENTATION
groups = unique(FIT_STRUCT.group);
designs = unique(FIT_STRUCT.design);
clusters = unique(FIT_STRUCT.cluster);
measures = unique(FIT_STRUCT.measure);
cl_i = 1;
des_i = 1;
meas_i = 1;
% GROUPINGS = {1,[2,3]};
GROUPINGS = {1,2,3};
groups_plot = 1:length(GROUPINGS);
inds = FIT_STRUCT.cluster == clusters(cl_i) &...
    FIT_STRUCT.design == designs(des_i) &...
    FIT_STRUCT.measure == measure(meas_i);
tmp_table = FIT_STRUCT(inds,:);
%## MANOVA
%-
mod_in = 'b0,b1,b2 ~ 1+group';
man_out = manova(tmp_table,mod_in);
stats_out = man_out.stats;
anova_p_var = stats_out.pValue;
% groupmeans(man_out,"b1")
%## CHANGE GROUPINGS?
if length(GROUPINGS) ~= length(groups)
    indsn = zeros(size(tmp_table,1),1);
    for g_i = 1:length(GROUPINGS)
        g1 = GROUPINGS{g_i};
        for g_j = 1:length(g1)
            indsg = tmp_table.group == groups(g1(g_j));
            % indsn(indsg) = g_i;
            tmp_table.group(indsg,:) = categorical(g_i);
        end
    end
end
%## PLOT
xname = 'x_raw';
yname = 'y_raw';
colors = linspecer(length(groups_plot));
colors_mean = colors*0.8;
figure;
hold on;
plot_leg = [];
for g_i = 1:length(groups_plot)
    indsg = tmp_table.group == categorical(groups_plot(g_i));
    tmp = tmp_table(indsg,:);
    xls = cellfun(@length,tmp.(xname));
    xn = zeros(length(tmp.(xname)),max(xls));
    x = cat(1,tmp.(xname){:});
    xu = unique(x);
    % yls = cellfun(@length,tmp.(yname));
    % yls = max(yls);
    %-
    % xn = cat(1,tmp.(xname){:});
    yn = nan(length(tmp.(yname)),length(xu));
    for subj_i = 1:length(tmp.(yname))
        yn(subj_i,:) = tmp.b0(subj_i) + tmp.b1(subj_i)*xu + tmp.b2(subj_i)*xu.^2;
        % tmpz = zeros(size(xu));
        % [~,xinds] = setdiff(xu,tmp.(xname){subj_i});
        % tmpz(xinds) = 1;
        % yn(subj_i,~logical(tmpz)) = tmp.(xname){subj_i};
    end
    % figure;
    xm = xu';
    ym = mean(yn,1);
    ys = std(yn,[],1);
    %-
    xf = [xm, fliplr(xm)];
    yf = [(ym + ys), fliplr((ym - ys))];
    ff = fill(xf, yf,'g');
    hold on;
    ff.FaceColor = colors(g_i,:);
    ff.FaceAlpha = 0.3;
    ff.EdgeColor = 'none';
    pp = plot(xm,ym,'Color',colors_mean(g_i,:),'DisplayName',sprintf('Group: %i',g_i),...
        'LineWidth',3);
    plot_leg = [plot_leg, pp];
end
xlim([-1,1]);
ylim([-5,5]);
legend(plot_leg);
hold off;
%}