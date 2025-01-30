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
        SRC_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        % SRC_DIR = getenv('STUDY_DIR');
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
    SRC_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
end
%## Add Study, Src, & Script Paths
addpath(SCRIPT_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%-
cmap_terrain = linspecer(4);
custom_yellow = [254,223,0]/255;
cmap_terrain = [cmap_terrain(3,:);custom_yellow;cmap_terrain(4,:);cmap_terrain(2,:)];
cmap_speed = linspecer(4*3);
cmap_speed = [cmap_speed(1,:);cmap_speed(2,:);cmap_speed(3,:);cmap_speed(4,:)];
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
%% (PATHS) ============================================================= %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
% STUDY_FNAME = 'kin_eeg_epoch_study';
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
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',save_dir);
% end
% %## LOAD STUDY
if ~ispc
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s_UNIX.study',STUDY_FNAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep sprintf('%s.study',STUDY_FNAME)]);
    STUDY = tmp.STUDY;
end
% %-
% cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
% STUDY.cluster = cl_struct;
% [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
% CLUSTER_PICS = main_cl_inds;
% %-
% save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
% if ~exist(save_dir,'dir')
%     mkdir(save_dir);
% end
%% MEASURES TO ANALYZE ================================================= %%
%% (LOAD STATISTICS & DATA EXCEL SHEET FROM R) ========================= %%
%-
tmp_dir = [cluster_k_dir];
% CLIN_TABLE = readtable([tmp_dir filesep 'fooof_clins_table_nans.xlsx'], ...
%     "FileType","spreadsheet","UseExcel",true);
sppb_table = readtable('M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\_data\MIM_dataset\_studies\subject_mgmt\sppb_analysis.xlsx', ...
    "FileType","spreadsheet","UseExcel",true);
% for i = 1:size(CLIN_TABLE,1)
%     if CLIN_TABLE.rej_cat_other(i)==1
%         CLIN_TABLE.subj_cat(i) = 1;
%     elseif CLIN_TABLE.rej_cat_other(i)==2 && CLIN_TABLE.cat_char_2(i)==3
%         CLIN_TABLE.subj_cat(i) = 3;
%     elseif CLIN_TABLE.rej_cat_other(i)==0
%         CLIN_TABLE.subj_cat(i) = 0;
%     else
%         CLIN_TABLE.subj_cat(i) = 2;
%     end
% end
% for i = 1:size(CLIN_TABLE,1)
%     if CLIN_TABLE.rej_cat_other(i)==1
%         CLIN_TABLE.subj_cat(i) = 1;
%     elseif CLIN_TABLE.rej_cat_other(i)==2 && CLIN_TABLE.rej_cat_cond(i)==1 && CLIN_TABLE.rej_cat_trial(i)==4
%         CLIN_TABLE.subj_cat(i) = 4;
%     elseif CLIN_TABLE.rej_cat_other(i)==2 && CLIN_TABLE.rej_cat_cond(i)==3 && CLIN_TABLE.rej_cat_trial(i)==4
%         CLIN_TABLE.subj_cat(i) = 3;
%     elseif CLIN_TABLE.rej_cat_other(i)==0 
%         CLIN_TABLE.subj_cat(i) = 0;
%     else
%         CLIN_TABLE.subj_cat(i) = 2;
%     end
% end
% for i = 1:size(CLIN_TABLE,1)
%     if ~CLIN_TABLE.is_processable(i)
%         CLIN_TABLE.subj_cat(i) = 0;
%     else
%         CLIN_TABLE.subj_cat(i) = 1*(CLIN_TABLE.is_considered(i));
%         CLIN_TABLE.subj_cat(i) = CLIN_TABLE.subj_cat(i)+3*CLIN_TABLE.is_in_study(i);
%         CLIN_TABLE.subj_cat(i) = CLIN_TABLE.subj_cat(i)+5*CLIN_TABLE.all_speed_conds(i);
%         CLIN_TABLE.subj_cat(i) = CLIN_TABLE.subj_cat(i)+7*CLIN_TABLE.all_speed_trials(i); 
%     end
% end
% CLIN_TABLE.group_code = categorical(CLIN_TABLE.group_code);
% CLIN_TABLE.subj_cat = categorical(CLIN_TABLE.subj_cat);
% %-
% fprintf('Subjs total: %i\n',size(CLIN_TABLE,1));
% fprintf('Subjs possible (non-fu, non-case, has mri): %i\n',sum(CLIN_TABLE.is_processable));
% fprintf('Subjs considered (non-fu, non-case, has mri, has rest): %i\n\n',sum(CLIN_TABLE.is_considered));
% %-
% chk0 = CLIN_TABLE.subj_cat==categorical(15);
% chk1 = (CLIN_TABLE.subj_cat==categorical(4));
% chk2 = (CLIN_TABLE.subj_cat==categorical(9));
% chk3 = (CLIN_TABLE.subj_cat==categorical(16));
% fprintf('Subjs in study total: %i\n',sum(CLIN_TABLE.is_in_study));
% fprintf('Subjs not considered, in study, all trials, all conds (error): %i\n',sum(chk0));
% fprintf('Subjs in study not all conds (error): %i\n',sum(chk1));
% fprintf('Subjs in study, all conds, not all trials: %i\n',sum(chk2));
% fprintf('Subjs in study, all conds, all trials: %i\n\n',sum(chk3));
% inds = (chk0 | chk1 | chk2 | chk3);
% CLIN_TABLE.subj_cat(chk0|chk1) = categorical(0);
% CLIN_TABLE.subj_cat(chk2) = categorical(2);
% CLIN_TABLE.subj_cat(chk3) = categorical(3);
% %-
% chk4 = CLIN_TABLE.subj_cat==categorical(0);
% chk5 = CLIN_TABLE.subj_cat==categorical(12);
% chk6 = CLIN_TABLE.subj_cat==categorical(5);
% chk7 = CLIN_TABLE.subj_cat==categorical(1);
% chk8 = CLIN_TABLE.subj_cat==categorical(13);
% chk9 = CLIN_TABLE.subj_cat==categorical(6);
% % chk10 = CLIN_TABLE.subj_cat==categorical(8)
% fprintf('Subjs not processable (no mri, other): %i\n',sum(chk4));
% fprintf('Subjs not considered, but all trials & conds: %i\n',sum(chk5));
% fprintf('Subjs not considered, not in study, all conds: %i\n',sum(chk6));
% fprintf('Subjs considered, not in study, not all conds, not all trials: %i\n',sum(chk7));
% fprintf('Subjs considered, not in study, but all trials & conds: %i\n',sum(chk8));
% fprintf('Subjs considered, but all conds: %i\n',sum(chk9));
% inds = (chk4 | chk5 | chk6 | chk7 | chk8 | chk9);
% CLIN_TABLE.subj_cat(inds) = categorical(1);
% inds = chk7;
% CLIN_TABLE.subj_cat(inds) = categorical(4);
% %-
% CLIN_TABLE.subj_cat = string(CLIN_TABLE.subj_cat);
% CLIN_TABLE.subj_cat = categorical(CLIN_TABLE.subj_cat);
% %-
% subj_cat_chars = {'not studied','missing trials','all trials','not studied missing conds'};
% % subj_cat_chars = {'np','p+ns+ncond+ntri','p+nc+ns+acond+ntri','p+c+ns+ncond+ntri','p+c+s+acond+ntri','p+nc+acond+atri','p+c+ns+acond+atri','p+c+s+acond+atri'};
% inds = isnan([CLIN_TABLE.sppb_12]);
% % CLIN_TABLE.sppb_12(inds) = 0;
% CLIN_TABLE = CLIN_TABLE(~inds,:);
%-
% inds = ~(CLIN_TABLE.subj_cat == categorical(0))
% CLIN_TABLE = CLIN_TABLE(inds,:);

%## MANUAL ENTRY
reject_table = readtable('M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\_data\MIM_dataset\_studies\subject_mgmt\subject_errors.xlsx', ...
    "Sheet","subj_errors_01132024","FileType","spreadsheet","UseExcel",true);
subj_cat_chars = {'all conds','missing conds','other'};
not_all_cond_subjs = {'H3046','H3047','H3053','H3063','H3073', ...
    'H3092','NH3023','NH3025','NH3028','NH3036','NH3051','NH3056','NH3071',...
    'NH3082','NH3123'};
all_conds_subjs = {STUDY.datasetinfo.subject};
% sppb_table.condition_status = categorical(repmat({''},size(sppb_table,1),1));
sppb_table.condition_status = cell(size(sppb_table,1),1);
for i = 1:size(sppb_table,1)
    if any(strcmp(sppb_table.subject_code(i),all_conds_subjs))
        sppb_table.condition_status{i} = 'all conditions';
    elseif any(strcmp(sppb_table.subject_code(i),not_all_cond_subjs))
        sppb_table.condition_status{i} = 'missing conditions';
    elseif any(strcmp(sppb_table.subject_code(i),reject_table.subject_name))
        sppb_table.condition_status{i} = 'other reject';
    else
        sppb_table.condition_status{i} = 'not considered';
    end
end
% inds = sppb_table.condition_status == categorical({'all conditions'}) | ...
%     sppb_table.condition_status == categorical({'missing conditions'}) | ...
%     sppb_table.condition_status == categorical({'other reject'});
inds = (sppb_table.condition_status == categorical({'all conditions'}) | ...
    sppb_table.condition_status == categorical({'missing conditions'}) | ...
    sppb_table.condition_status == categorical({'other reject'})) & ...
    any(sppb_table.group_code == [2,3],2);
% inds = (sppb_table.condition_status == categorical({'all conditions'}) | ...
%     sppb_table.condition_status == categorical({'missing conditions'})) & ...
%     any(sppb_table.group_code == [2,3],2);
sppb_table = sppb_table(inds,:);
sppb_table.condition_status = categorical(sppb_table.condition_status);

%% ===================================================================== %%
%## PARAMETERS
cmaps_speed = linspecer(3);
%-
SAVE_RES = 300;
TITLE_TXT_SIZE = 14;
IM_RESIZE = 1.3;
AX_W = 0.4;
AX_H = 0.25;
AX_FONT_NAME = 'Arial';
AX_X_SHIFT = 1.2;
AX_Y_SHIFT = 1.55;
AX_INIT_X = 0.1;
AX_INIT_Y = 0.6;
%## 
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
%## VIOLIN PLOTS
PLOT_STRUCT = struct('color_map',cmaps_speed,...
    'cond_labels',{subj_cat_chars},... %cellstr(unique(CLIN_TABLE.subj_cat))
    'cond_offsets',[-0.05,0.05,0.15],...
    'do_group_labels',false, ...
    'group_labels',{{}},...
    'group_offsets',0,...
    'group_lab_yoffset',-0.1,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',10,...
    'y_label','SPPB Score Total',...
    'y_label_fontsize',10,...
    'y_label_fontweight','bold',...
    'ylim',[],...
    'x_label','Rejection Status',...
    'x_label_fontsize',10,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.05,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',10,...
    'title_fontweight','bold',...
    'font_size',9,...
    'font_name','Arial',...
    'do_combine_groups',true,...
    'regresslab_txt_size',10,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1,...
    'xtick_angle',0);
VIOLIN_STRUCT = struct('Width',0.03,...
    'ShowWhiskers',false,...
    'ShowNotches',false,...
    'ShowBox',true,...
    'ShowMedian',true,...
    'Bandwidth',0.075,...
    'QuartileStyle','shadow',...
    'HalfViolin','full',...
    'DataStyle','scatter',...
    'MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],...
    'ViolinAlpha',{{0.2 0.3}},...
    'do_plot_outlier_marks',true,...
    'use_raw_bandwidth',false);
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
    'text_specs',{{'FontSize',12,'FontName','Arial','FontWeight','bold'}},...
    'conn_y',[],...
    'conn_offset_y',[],...
    'sig_levels',[0.05,0.01,0.001],...
    'sig_offset_x',0,...
    'sig_offset_y',0);
%% ===================================================================== %%
tmp_savedir = save_dir;
mkdir(tmp_savedir);
%## SET SHIFTS
x_cnt = 1;
y_shift = AX_INIT_Y;
x_shift = AX_INIT_X;
%## INITIATE FIGURE
fig = figure('color','white');
sgtitle('SPPB Scores','FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
set(fig, ...
    'Units','inches', ...
    'Position',[0.5,0.5,6.5,9])
set(fig, ...
    'PaperUnits','inches', ...
    'PaperSize',[1 1], ...
    'PaperPosition',[0 0 1 1])
hold on;
set(gca,AXES_DEFAULT_PROPS{:})

%## STATS
modo = fitlm(sppb_table,'sppb_12 ~ 1+condition_status');
anvo = anova(modo);
[p,t,anova_out,terms] = anovan([sppb_table.sppb_12],{sppb_table.condition_status},...
    'sstype',3,'varnames',{'condition_status'},'model','linear','Display','off');
[comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',1,...
    'display','off','Alpha',0.05,'CriticalValueType','bonferroni'); % compaarisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
disp(modo)
% [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
adj_p = comparisons(:,6);
%-
if t{2,7} < 0.001
    str = '***';
elseif t{2,7} > 0.001 && t{2,7} <= 0.01
    str = '**';
elseif t{2,7} < 0.05
    str = '*';
else
    str = '';
end
tmp_stats_struct = struct('anova',{{t{2,7}}},...
    'anova_grp',{{}},...
    'pvals',{{adj_p}},...
    'pvals_pairs',{{num2cell(comparisons(:,1:2),2)}},...
    'pvals_grp',{{}},...
    'pvals_grp_pairs',{{}},...
    'regress_pval',{{}},...
    'regress_line',{{}},...
    'line_type',{'best_fit'},... % ('best_fit' | 'means')
    'regress_xvals',[],...
    'subject_char',[],... % this option when filled prints removal of nan() info
    'group_order',categorical({''}),...
    'display_stats_char',true,... 
    'stats_char',{{sprintf('%sF = %0.2f\n',str,t{2,6})}},...
    'stats_char_offsets',[0.05,0]);
%## PLOT
ax = axes();
%-
PLOT_STRUCT.y_label ='SPPB Total Score';
PLOT_STRUCT.ylim = [0,13];
PLOT_STRUCT.ax_position = [x_shift,y_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];
%-
% ax = group_violin(CLIN_TABLE,'sppb_12','design_id','group_char',...
%     ax,...
%     'VIOLIN_PARAMS',VIO_PARAMS,...
%     'PLOT_STRUCT',VIO_PLOT_STRUCT,...
%     'STATS_STRUCT',tmp_stats_struct,...
%     'BRACKET_STRUCT',DEF_BRACKET_STRUCT,...
%     'SIGLINE_STRUCT',DEF_SIGLINE_STRUCT);


ax = group_violin(sppb_table,'sppb_12','condition_status','dummy_var',...
    ax,...
    'VIOLIN_STRUCT',VIOLIN_STRUCT,...
    'PLOT_STRUCT',PLOT_STRUCT,...
    'STATS_STRUCT',tmp_stats_struct,...
    'BRACKET_STRUCT',DEF_BRACKET_STRUCT,...
    'SIGLINE_STRUCT',DEF_SIGLINE_STRUCT);
%## AX SHIFT
if x_cnt < X_DIM
    x_shift = x_shift + AX_X_SHIFT*IM_RESIZE*AX_W;
else
    y_shift = y_shift - AX_Y_SHIFT*IM_RESIZE*AX_H;
    x_shift = AX_INIT_X;
    x_cnt = 0;
end
x_cnt = x_cnt + 1;
exportgraphics(fig,[tmp_savedir filesep sprintf('sppb_plot_allsubj.tiff')], ...
    'Resolution',SAVE_RES)
% close(fig)