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
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
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
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',spec_data_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',spec_data_dir);
% end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
save_dir = [spec_data_dir filesep 'psd_calcs'];
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
%{
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
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for i = 1:length(designs)
    des_i = string(designs(i));
    for j = 1:length(clusters) % For each cluster requested
        cl_i = double(string(clusters(j)));
        if isempty(tmp_study.cluster(cl_i).topo) || RE_CALC
            inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == string(cl_i));
            sets_i = unique([FOOOF_TABLE.subj_cl_ind(inds)]);
            tmp_study.cluster(cl_i).sets = tmp_study.cluster(cl_i).sets(sets_i);
            tmp_study.cluster(cl_i).comps = tmp_study.cluster(cl_i).comps(sets_i);
            tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
            STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
            STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
            STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
            STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
            STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
        end
    end
end
%}
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
measure_name_plot = {'theta_avg_power','alpha_avg_power','beta_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \alpha','Mean \beta'};
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
kin_savedir = [save_dir filesep 'kinematics'];
mkdir(kin_savedir);
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
design_chars = {'terrain','speed'};
group_chars = unique(FOOOF_TABLE.group_char);
conditions = unique(FOOOF_TABLE.cond_char);
xlsx_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'raw_data_vis'];
imu_table = [xlsx_fpath filesep 'imu_table_meantrial.xlsx'];
ls_table = [xlsx_fpath filesep 'ls_table_meantrial.xlsx'];
% imu_table = [xlsx_fpath filesep 'imu_table_out.xlsx'];
% ls_table = [xlsx_fpath filesep 'ls_table_out.xlsx'];
imu_table = readtable(imu_table);
ls_table = readtable(ls_table);
PLOT_STRUCT = struct('color_map',linspecer(4),...
                'cond_labels',unique(FOOOF_TABLE.cond_char),'group_labels',unique(FOOOF_TABLE.group_char),...
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label','10*log_{10}(Flattened PSD)',...
                'title','','font_size',9,'y_lim',[-1,15],...
                'font_name','Arial','x_label','');
%% ================================================================== %%

%## ATTACH KINEMATICS TO SPECTRUM TABLE
NON_KINS = {'speed_ms','subj_id','subj_char','comp_id','design_id',...
    'cond_id','cond_char','group_id','cluster_id','aperiodic_exp',...
    'aperiodic_offset','central_freq','power','r_squared',...
    'theta_avg_power','alpha_avg_power','beta_avg_power','theta',...
    'alpha','beta','theta_center','alpha_center','beta_center',...
    'group_char','SubjectName','TrialName','SubjectCategory','terrain_speed'};
%## BIG BIG MODEL
varnames = imu_table.Properties.VariableNames;
total_miss = {};
%-
for var_i = 1:length(varnames)
    for i = 1:size(FOOOF_TABLE,1)
        ind_sub = imu_table.SubjectName == string(FOOOF_TABLE.subj_char(i));
        tmp = strsplit(string(FOOOF_TABLE.cond_char(i)),'.');
        tmp = strjoin(tmp,'p');
        tt = imu_table(ind_sub,:);
        ind_cond = tt.TrialName == tmp;
        if ~any(strcmp(varnames{var_i},NON_KINS))
            vals = tt.(varnames{var_i})(ind_cond);
            if length(vals) > 1
                % fprintf('Multiple entries for this item on subject %s.\n',FOOOF_TABLE.subj_char(i));
                ii = isnan(vals);
                if all(ii)
                    FOOOF_TABLE.(varnames{var_i})(i) = nan();
                else
                    vals = vals(~ii);
                    FOOOF_TABLE.(varnames{var_i})(i) = vals(1);
                end
            elseif isempty(vals)
                fprintf('Subject %s does not have an entry for this measure.\n',FOOOF_TABLE.subj_char(i));
                if isempty(total_miss) || ~any(total_miss == FOOOF_TABLE.subj_char(i))
                    total_miss = [total_miss FOOOF_TABLE.subj_char(i)];
                end
            else
                ii = isnan(vals);
                if all(ii)
                    FOOOF_TABLE.(varnames{var_i})(i) = nan();
                else
                    vals = vals(~ii);
                    FOOOF_TABLE.(varnames{var_i})(i) = vals(1);
                end
            end
        end
    end
end
disp(total_miss);
%## BIG BIG MODEL
varnames = ls_table.Properties.VariableNames;
total_miss = {};
%-
for var_i = 1:length(varnames)
    for i = 1:size(FOOOF_TABLE,1)
        ind_sub = ls_table.SubjectName == string(FOOOF_TABLE.subj_char(i));
        tmp = strsplit(string(FOOOF_TABLE.cond_char(i)),'.');
        tmp = strjoin(tmp,'p');
        tt = ls_table(ind_sub,:);
        ind_cond = tt.TrialName == tmp;
        if ~any(strcmp(varnames{var_i},NON_KINS))
            vals = tt.(varnames{var_i})(ind_cond);
            if length(vals) > 1
                % fprintf('Multiple entries for this item on subject %s.\n',FOOOF_TABLE.subj_char(i));
                ii = isnan(vals);
                if all(ii)
                    FOOOF_TABLE.(varnames{var_i})(i) = nan();
                else
                    vals = vals(~ii);
                    FOOOF_TABLE.(varnames{var_i})(i) = vals(1);
                end
            elseif isempty(vals)
                fprintf('Subject %s does not have an entry for this measure.\n',FOOOF_TABLE.subj_char(i));
                if isempty(total_miss) || ~any(total_miss == FOOOF_TABLE.subj_char(i))
                    total_miss = [total_miss FOOOF_TABLE.subj_char(i)];
                end
            else
                ii = isnan(vals);
                if all(ii)
                    FOOOF_TABLE.(varnames{var_i})(i) = nan();
                else
                    vals = vals(~ii);
                    FOOOF_TABLE.(varnames{var_i})(i) = vals(1);
                end
            end
        end
    end
end
disp(total_miss);
writetable(FOOOF_TABLE,[save_dir filesep 'fooof_kinematics_table.xlsx'])
par_save(FOOOF_TABLE,save_dir,'fooof_kinematics_table.mat');

%% ===================================================================== %%
% FOOOF_TABLE = readtable([save_dir filesep 'fooof_kinematics_table.xlsx']);
FOOOF_TABLE = par_load(save_dir,'fooof_kinematics_table.mat');
%-
NON_KINS1 = {'subj_cl_ind','central_freq_1','central_freq_2','central_freq_3','power_1','power_2','power_3','theta_1',...
    'theta_2','theta_3','theta_4','alpha_1','alpha_2','alpha_3','alpha_4','alpha_5','alpha_6',...
    'beta_1','beta_2','beta_3','beta_4','beta_5','beta_6','theta_center_1','theta_center_2',...
    'alpha_center_1','alpha_center_2','beta_center_1','beta_center_2'};
NON_KINS2 = {'speed_ms','subj_char','comp_id',...
    'cond_char','aperiodic_exp',...
    'aperiodic_offset','central_freq','power','r_squared','theta',...
    'alpha','beta','theta_center','alpha_center','beta_center',...
    'group_char','SubjectName','TrialName','SubjectCategory','terrain_speed'};
MAIN_EFFS = {'cond_id','design_id','group_id','cluster_id','subj_id'};
MAIN_RESP = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
VARS_CHK = [MAIN_EFFS,MAIN_RESP];
%- (04/01/2024) JS, the purpose of this section is to reduce the datamatrix
% (aka kinematics & clinical outcomes) to a rank non-deficient state 
% (rank = 0). Not sure if its appropriate to also include the effects with
% the responses or to only perform it on the effects.
%##
inds = [];
for i = 1:length(FOOOF_TABLE.Properties.VariableNames)
    if ~any(strcmp(FOOOF_TABLE.Properties.VariableNames{i},[NON_KINS1,NON_KINS2,VARS_CHK]))
        inds = [inds, i];
    end
end
vv_main_anl = {'mean_APexc_mean','mean_APexc_COV','mean_MLexc_COV','mean_StepDur'};
varnames = FOOOF_TABLE.Properties.VariableNames(inds);
inds_l = contains(varnames,'_L');
inds_r = contains(varnames,'_R');
varnames = varnames(~(inds_l|inds_r));
varnames_labs = cell(length(varnames),1);
varnames_hold = varnames;
for i = 1:length(varnames)
    tmp = strsplit(varnames{i},'_');
    chk = contains(tmp,'mean');
    i1 = find(chk);
    i2 = find(~chk);
    % tmp = strjoin([tmp(i1(1)),tmp(i2)],' ');
    tmp = strjoin([tmp(i2)],' ');
    varnames_labs{i} = tmp;
end
TMP_FOOOF_T = FOOOF_TABLE;
%## HARDCODED VARIABLES OF INTEREST
varnames = {'mean_APexc_COV','mean_APexc_mean','mean_MLexc_COV',...
            'mean_MLexc_mean','mean_StepDur','mean_UDexc_COV',...
            'mean_UDexc_mean','mean_StanceDur','mean_GaitCycleDur','mean_PeakUpDownVel_mean'};
varnames_labs = {'APexc COV','APexc','MLexc COV','MLexc','Step Dur','UDexc COV','UDexc','Stance Dur','GaitCycle Dur','Peak UpDown Vel'};
%% PREDICTORS: SPEED. RESPONSE: BRAIN ACTIVITY, TEST FOR LINEARITY
% STATS_OUT = [];
% im_resize= 1.2;
VIOLIN_BOTTOM = 0.65;
% AX_H  = 0.2;
% AX_W = 0.25;
DO_PLOT_GROUPS = true;
% REGRESS_TXT_SIZE = 8;
% REGRESS_TXT_XMULTI = 0.9;
% REGRESS_TXT_YMULTI = 0.9;
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
im_resize = 0.8;
AX_W = 0.35;
AX_H = 0.25;
% GROUP_SHORTS = {'YA','HO','FO'};
% GROUP_MARKERS = {'o','x','*'};
tmp_savedir = [save_dir filesep 'Pspeed-Reeg-linearity'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %
    atlas_name = atlas_name_store{cl_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    for des_i = 2 %## JUST SPEED
        %##
        horiz_shift = 0;
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
        for meas_i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{meas_i};
            inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
            T_vals_plot = TMP_FOOOF_T(inds,:);
            T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
            subjects = unique(T_vals_plot.subj_char);
            conds = unique(T_vals_plot.cond_id);
            %##
            STATS_STRUCT = struct();
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_STRUCT = struct('color_map',color_dark,...
                'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
                'cond_offsets',[-0.35,-0.1,0.15,0.4],'group_offsets',[0.125,0.475,0.812],...
                'y_label',measure_name,...
                'title',measure_name,'font_size',8,'ylim',[min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+3*std(T_vals_plot.(measure_name))],...
                'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);

            hold on;
            axax = group_violin(T_vals_plot,measure_name,'cond_char','group_char',...
                fig,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_STRUCT,...
                'STATS_STRUCT',STATS_STRUCT);
            %-
            hold on;
            cnt_g = 1;
            xtxt_shft = 0;
            ytxt_shft = 0; 
            xticks = get(axax,'XTick');
            for group_i = 1:length(groups)
                if DO_PLOT_GROUPS
                    inds = T_vals_plot.group_char == string(group_chars(group_i));
                    data = T_vals_plot(inds,:);
                else
                    data = T_vals_plot;
                end
                x_plot = xticks(cnt_g:cnt_g+length(conds)-1);
                y = data.(measure_name);
                x = data.cond_char;
                p2 = polyfit(x,y,2);
                v2 = polyval(p2,x);              % Evaluate
                p1 = polyfit(x,y,1);
                v1 = polyval(p1,x);              % Evaluate
                TSE = sum((v2 - y).^2)-sum((v1 - y).^2);         % Total Squared Error
                %-
                y_lim_calc = [min(y)-std(y),max(y)+3*std(y)];
                x = unique(data.cond_char);
                y = [];
                for i = 1:length(x)
                    y(i) = p2(1)*x(i)^2 + p2(2)*x(i) + p2(3);
                end
                plot(x_plot,y,...
                    'DisplayName',sprintf('polyfit^2'),...
                    'LineWidth',2,'LineStyle','-.','Color','k');
                hold on;
                x = unique(data.cond_char);
                y = [];
                for i = 1:length(x)
                    y(i) = p1(1)*x(i) + p1(2);
                end
                plot(x_plot,y,...
                    'DisplayName',sprintf('polyfit^2'),...
                    'LineWidth',2,'LineStyle','--','Color','k');
                hold on;
                eq = sprintf('y=(%0.1g)*x^2+(%0.1g)*x+(%0.1g)',p2(1),p2(2),p2(3));
                text(0.1+xtxt_shft,VIOLIN_BOTTOM-vert_shift+0.025+ytxt_shft,sprintf('TSE_{quad-lin}= %0.2g',TSE),...
                    'FontSize',6,...
                    'FontName',PLOT_STRUCT.font_name,...
                    'FontWeight','bold','Units','normalized');
                cnt_g = cnt_g + length(conds);
                xtxt_shft = xtxt_shft + 0.25;
                ytxt_shft = ytxt_shft + 0.1;
                if ~DO_PLOT_GROUPS
                    break;
                end
            end
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('');
            end
            xlabel('Speed (m/s)');
            title(MEASURE_NAME_LABS{meas_i});
            set(gca,'FontWeight','bold');
            ylim(y_lim_calc)
            set(gca,'OuterPosition',[0,0,1,1]);
            set(gca,'Position',[0.075+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
            % hold off;
            % set(gca,'Position',[horiz_shift+0.05,SCATTER_BOTTOM-vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_H*im_resize + 0.115;
            hold on;
        end
        if DO_PLOT_GROUPS
            exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_speed-eeg-group_linearity.tiff',string(clusters(cl_i)))],'Resolution',300)
        else
            exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_speed-eeg_linearity.tiff',string(clusters(cl_i)))],'Resolution',300)
        end
        close(fig)
   end
end
%% PREDICTORS: SPEED CONDITION, RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 1.2;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
DO_PLOT_GROUPS = false;
tmp_savedir = [save_dir filesep 'Pspeed-Rkin'];
mkdir(tmp_savedir);
for var_i = 1:length(varnames)
    %
    vert_shift = 0;
    for des_i = 2 %## JUST SPEED
        %##
        horiz_shift = 0;
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
        inds = TMP_FOOOF_T.design_id == designs(des_i);
        T_vals_plot = TMP_FOOOF_T(inds,:);
        subjects = unique(T_vals_plot.subj_char);
        conds = unique(T_vals_plot.cond_id);
        % groups = unique(T_vals_plot.group_id);
        t_tmp = [];
        for i = 1:length(subjects)
            ii = find(T_vals_plot.subj_char == subjects(i));
            tt = T_vals_plot(ii,:);
            for j = 1:length(conds)
                jj = find(tt.cond_id == conds(j));
                t_tmp = [t_tmp; tt(jj(1),:)];
            end
        end
        T_vals_plot = table(categorical(string(t_tmp.cond_char)),t_tmp.(varnames{var_i}),categorical(string(t_tmp.group_char)),...
           'VariableNames',{'cond_char',varnames{var_i},'group_char'});
        % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        try
            mod = sprintf('%s ~ 1 + %s',varnames{var_i},'cond_char');
            % stats_out = fitlme(T_vals_plot,mod);
            stats_out = fitlm(T_vals_plot,mod);
            % anova_out = anova(stats_out);
            out = anova(T_vals_plot,mod,'SumOfSquaresType',"three",'CategoricalFactors',{'cond_char'},...
                'ModelSpecification','linear');
            anova_out = out.stats();
            % anova_out = anovan(double(T_vals_plot.(varnames{var_i})),{T_vals_plot.cond_char},...
            %             'model','linear',...
            %             'model',1,...
            %             'sstype',3,...
            %             'varnames',strvcat('speed'));

            %## PRINT TABLES
            t = sprintf_table(anova_out);
            t.saveToFile([tmp_savedir filesep sprintf('%s_kinematics-speed_ANOVA.txt',varnames{var_i})]);
            t = sprintf_table(stats_out.Coefficients);
            t.saveToFile([tmp_savedir filesep sprintf('%s_kinematics-speed_LM.txt',varnames{var_i})]);
            %-
            R2 = stats_out.Rsquared.Adjusted;
            anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
            pval_0p5 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
            pval_0p75 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
            pval_1p0 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
        catch e
            fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
            R2 = 0;
            pval = 1;
            slope = 0;
            inter = 0;
        end
        %##
        STATS_STRUCT = struct('anova',{{}},...
                      'anova_grp',{{}},...
                      'pvals',{{}},...
                      'pvals_pairs',{{}},...
                      'pvals_grp',{{}},...
                      'pvals_grp_pairs',{{}},...
                      'regress_pval',{{}},...
                      'regress_line',{{}},...
                      'r2_coeff',{[]},...
                      'regress_xvals',0);
        if DO_PLOT_GROUPS
            for gg = 1:length(groups)
                STATS_STRUCT.anova{gg}=anova_p_var;
                STATS_STRUCT.pvals_pairs{gg}={[1,1],[1,2],[1,3],[1,4]};
                STATS_STRUCT.pvals{gg}=[pval_inter,pval_0p5,pval_0p75,pval_1p0];
            end
        else
            STATS_STRUCT.anova{1}=anova_p_var;
            STATS_STRUCT.pvals_pairs{1}={[1,1],[1,2],[1,3],[1,4]};
            STATS_STRUCT.pvals{1}=[pval_inter,pval_0p5,pval_0p75,pval_1p0];
        end
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_STRUCT = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.35,-0.1,0.15,0.4],'y_label',varnames_labs{var_i},...
            'title',varnames_labs{var_i},'font_size',10,'group_offsets',[0.125,0.475,0.812],...
            'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+3*std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);
        fig = figure('color','white','renderer','Painters');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_char',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[tmp_savedir filesep sprintf('%s_kinematics-speed-factor_grouped.tiff',varnames{var_i})],'Resolution',300)
        % exportgraphics(fig,[tmp_savedir filesep sprintf('%s_kinematics-speed-contin_grouped.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% ===================================================================== %%
%{
pval_out = zeros(length(STATS_OUT),1);
lens = zeros(length(STATS_OUT),1);
for i = 1:length(STATS_OUT)
    pval_out(i) = STATS_OUT(i).anova{1};
    lens(i) = length(varnames{i});
end
[val,ind] = sort(pval_out);
varnames = varnames(ind);
format = sprintf('\n\n   %%s%%%is\t\t%%s\n',max(lens)-length('Name'));
fid = fopen([kin_savedir filesep 'kinematics_pvals.txt']);
fprintf(fid,format,'Name','','pValue');
for i = 1:length(varnames)
    format = sprintf('%i) %%s%%%is\t\t%%0.2g\n',i,max(lens)-length(varnames{i}));
    fprintf(fid,format,varnames{i},'',val(i));
end
fclose(fid);
varnames = varnames(1:3);
varnames = unique([varnames, vv_main_anl]);
inds = cellfun(@(x) find(strcmp(x,varnames_hold)),varnames);
varnames_labs = varnames_labs(inds);
%}
%% PREDICTORS: SPEED factor, GROUP, & INTERACTION; RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 1.2;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
DO_PLOT_GROUPS = true;
tmp_savedir = [save_dir filesep 'Pspeed-Pgroup-Pinter-Rkin'];
mkdir(tmp_savedir);
for var_i = 1:length(varnames)
    vert_shift = 0;
    for des_i = 2 %## JUST SPEED
        %##
        horiz_shift = 0;
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
        inds = TMP_FOOOF_T.design_id == designs(des_i);
        T_vals_plot = TMP_FOOOF_T(inds,:);
        subjects = unique(T_vals_plot.subj_char);
        conds = unique(T_vals_plot.cond_char);
        % groups = unique(T_vals_plot.group_id);
        t_tmp = [];
        for i = 1:length(subjects)
            ii = find(T_vals_plot.subj_char == subjects(i));
            tt = T_vals_plot(ii,:);
            for j = 1:length(conds)
                jj = find(tt.cond_char == conds(j));
                t_tmp = [t_tmp; tt(jj(1),:)];
            end
        end
        T_vals_plot = table(categorical(string(t_tmp.cond_char)),t_tmp.(varnames{var_i}),categorical(string(t_tmp.group_char)),...
           'VariableNames',{'cond_char',varnames{var_i},'group_char'});
        % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        try
            mod = sprintf('%s ~ 1 + cond_char + group_char + cond_char*group_char',varnames{var_i});
            % stats_out = fitlme(T_vals_plot,mod);
            stats_out = fitlm(T_vals_plot,mod);
            anova_out = anova(T_vals_plot,mod,'SumOfSquaresType',"three",'CategoricalFactors',{'group_char','cond_char'},...
                'ModelSpecification','interactions');
            anova_out = anova_out.stats();
            % anova_out = anovan(double(T_vals_plot.(varnames{var_i})),{T_vals_plot.group_char,T_vals_plot.cond_char},...
            %             'model','interaction',...
            %             'model',2,...
            %             'sstype',3,...
            %             'varnames',strvcat('group','speed'));
            R2 = stats_out.Rsquared.Adjusted;
            %## PRINT TABLES
            t = sprintf_table(anova_out);
            t.saveToFile([tmp_savedir filesep sprintf('%s_ANOVA.txt',varnames{var_i})]);
            t = sprintf_table(stats_out.Coefficients);
            t.saveToFile([tmp_savedir filesep sprintf('%s_LM.txt',varnames{var_i})]);
            %-
            anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
            anova_p_group = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'group_char'));
            anova_p_inter = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char:group_char'));
            %-
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
            pval_0p5 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
            pval_0p75 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
            pval_1p0 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
            %-
            pval_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H2000''s'));
            pval_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s'));
            %-
            pval_var_0p5_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5:group_char_H2000''s'));
            pval_var_0p75_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75:group_char_H2000''s'));
            pval_var_1p0_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0:group_char_H2000''s'));
            %-
            pval_var_0p5_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5:group_char_H3000''s'));
            pval_var_0p75_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75:group_char_H3000''s'));
            pval_var_1p0_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0:group_char_H3000''s'));

            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id')));
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id')));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
        catch e
            fprintf('Error. Variable %s\n\n%s\n',string(varnames{var_i}),getReport(e))
            R2 = 0;
            pval = 1;
            slope = 0;
            inter = 0;
        end
        DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{[]},...
                          'regress_xvals',0);
        STATS_STRUCT.anova_grp = {anova_p_group,anova_p_group,anova_p_group};
        STATS_STRUCT.pvals_grp_pairs = {[1,1],[1,2],[1,3]};
        STATS_STRUCT.pvals_grp = {1,pval_g2,pval_g3};
        STATS_STRUCT.anova={anova_p_var,anova_p_inter,anova_p_inter};
        STATS_STRUCT.pvals = {[1,pval_0p5,pval_0p75,pval_1p0],...
            [1,pval_var_0p5_g2,pval_var_0p75_g2,pval_var_1p0_g2],...
            [1,pval_var_0p5_g3,pval_var_0p75_g3,pval_var_1p0_g3]};
        STATS_STRUCT.pvals_pairs = {{[1,1],[1,2],[1,3],[1,4]},...
            {[1,1],[1,2],[1,3],[1,4]},...
            {[1,1],[1,2],[1,3],[1,4]}};
        % STATS_STRUCT.regress_pval={pval_var};
        % STATS_STRUCT.regress_line={[inter_mn,slope_var]};
        % STATS_STRUCT.r2_coeff=R2;
        % STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_STRUCT = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},'group_offsets',[0.125,0.475,0.812],...
            'title',varnames_labs{var_i},'font_size',10,'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+3*std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);
        % ax = axes();
        fig = figure('color','white','renderer','Painters');
        % sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_char',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[tmp_savedir filesep sprintf('%s_kinematics-speed-group-interact.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% PREDICTORS: SPEED Factor, GROUP; RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 1.2;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
DO_PLOT_GROUPS = true;
tmp_savedir = [save_dir filesep 'Pspeed-Pgroup-Rkin'];
mkdir(tmp_savedir);
for var_i = 1:length(varnames)
    vert_shift = 0;
    for des_i = 2 %## JUST SPEED
        %##
        horiz_shift = 0;
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
        inds = TMP_FOOOF_T.design_id == designs(des_i);
        T_vals_plot = TMP_FOOOF_T(inds,:);
        subjects = unique(T_vals_plot.subj_char);
        conds = unique(T_vals_plot.cond_char);
        % groups = unique(T_vals_plot.group_id);
        t_tmp = [];
        for i = 1:length(subjects)
            ii = find(T_vals_plot.subj_char == subjects(i));
            tt = T_vals_plot(ii,:);
            for j = 1:length(conds)
                jj = find(tt.cond_char == conds(j));
                t_tmp = [t_tmp; tt(jj(1),:)];
            end
        end
        T_vals_plot = table(categorical(string(t_tmp.cond_char)),t_tmp.(varnames{var_i}),categorical(string(t_tmp.group_char)),...
           'VariableNames',{'cond_char',varnames{var_i},'group_char'});
        % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        try
            mod = sprintf('%s ~ 1 + cond_char + group_char',varnames{var_i});
            % stats_out = fitlme(T_vals_plot,mod);
            stats_out = fitlm(T_vals_plot,mod);
            % anova_out = anova(stats_out);
            out = anova(T_vals_plot,mod,'SumOfSquaresType',"three",'CategoricalFactors',{'group_char','cond_char'},...
                'ModelSpecification','linear');
            anova_out = out.stats();
            % anova_out = anovan(double(T_vals_plot.(varnames{var_i})),{T_vals_plot.group_char,T_vals_plot.cond_char},...
            %             'model','interaction',...
            %             'model',1,...
            %             'sstype',3,...
            %             'varnames',strvcat('group','speed'));
            
            %## PRINT TABLES
            t = sprintf_table(anova_out);
            t.saveToFile([tmp_savedir filesep sprintf('%s_ANOVA.txt',varnames{var_i})]);
            t = sprintf_table(stats_out.Coefficients);
            t.saveToFile([tmp_savedir filesep sprintf('%s_LM.txt',varnames{var_i})]);

            R2 = stats_out.Rsquared.Adjusted;
            anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
            anova_p_group = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'group_char'));
            %- 
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
            pval_0p5 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
            pval_0p75 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
            pval_1p0 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
            %-
            pval_grp2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H2000''s'));
            pval_grp3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s'));
            %-
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
        catch e
            fprintf('Error. Variable %s\n\n%s\n',string(varnames{var_i}),getReport(e))
            R2 = 0;
            pval = 1;
            slope = 0;
            inter = 0;
        end
        STATS_STRUCT = struct('anova',{{}},...
                              'anova_grp',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'pvals_grp',{{}},...
                              'pvals_grp_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{[]},...
                              'regress_xvals',0);
        STATS_STRUCT.anova_grp = {anova_p_group,anova_p_group,anova_p_group};
        STATS_STRUCT.pvals_grp_pairs = {[1,1],[1,2],[1,3]};
        STATS_STRUCT.pvals_grp = {1,pval_grp2,pval_grp3};
        STATS_STRUCT.anova={anova_p_var,anova_p_var,anova_p_var};
        STATS_STRUCT.pvals = {[1,pval_0p5,pval_0p75,pval_1p0],...
            [1,1,1,1],...
            [1,1,1,1]};
        STATS_STRUCT.pvals_pairs = {{[1,1],[1,2],[1,3],[1,4]},...
            {[1,1],[1,2],[1,3],[1,4]},...
            {[1,1],[1,2],[1,3],[1,4]}};
        % STATS_STRUCT.regress_pval={pval_var};
        % STATS_STRUCT.regress_line={[inter_mn,slope_var]};
        % STATS_STRUCT.r2_coeff=R2;
        % STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_STRUCT = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},'group_offsets',[0.125,0.475,0.812],...
            'title',varnames_labs{var_i},'font_size',10,'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+3*std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);
        % ax = axes();
        fig = figure('color','white','renderer','Painters');
        % sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_char',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_STRUCT',PLOT_STRUCT,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[tmp_savedir filesep sprintf('%s_kinematics-speed-group-interact.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% PREDICTORS: SPEED Factor, Kinematics, INTERACTION. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
% REGRESS_TXT_XMULTI = 0.9;
% REGRESS_TXT_YMULTI = 1.0;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25*3;
SCATTER_BOTTOM = 0.175;
GROUP_SHORTS = {'YA','HO','FO'};
GROUP_MARKS = {'o','x','^'};
GROUP_LINESTYLES = {'-','-.',':'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
PLOT_STRUCT = struct('color_map',[],...
        'cond_labels',unique(TMP_FOOOF_T.cond_char),'group_labels',unique(TMP_FOOOF_T.group_char),...
        'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',[],...
        'title',[],'font_size',10,'ylim',[],...
        'font_name','Arial','x_label','speed','do_combine_groups',true);
tmp_savedir = [save_dir filesep 'Pspeedf-Pkin-Pinter-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %%
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white','renderer','Painters');
        sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,3])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        for des_i = 2
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
            horiz_shift = 0;
            stats_store = [];
            for meas_i = 1:length(measure_name_plot)
                measure_name = measure_name_plot{meas_i};
                %##
                cond_plot_store = [];
                % inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(group_i);
                % T_stats_plot = psd_feature_stats(inds,:);
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
                T_vals_plot = TMP_FOOOF_T(inds,:);
                % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                T_vals_plot = table(double(T_vals_plot.(measure_name)),categorical(string(T_vals_plot.cond_char)),categorical(string(T_vals_plot.group_char)),double(T_vals_plot.(varnames{var_i})),...
                    'VariableNames',{measure_name,'cond_char','group_char',varnames{var_i}});
                try
                    mod = sprintf('%s ~ 1 + cond_char + %s + cond_char:%s',measure_name,varnames{var_i},varnames{var_i});
                    % stats_out = fitlme(T_vals_plot,mod);
                    stats_out = fitlm(T_vals_plot,mod);
                    % anova_out = anova(stats_out);
                    out = anova(T_vals_plot,mod,'SumOfSquaresType',"three",'CategoricalFactors',{'cond_char'},...
                            'ModelSpecification','interactions');
                    anova_out = out.stats();
                    % anova_out = anovan(double(T_vals_plot.(measure_name)),{T_vals_plot.group_char,T_vals_plot.cond_char},...
                    %         'model','interaction',...
                    %         'model',2,...
                    %         'sstype',3,...
                    %         'varnames',strvcat('group','speed'));
                    R2 = stats_out.Rsquared.Adjusted;
                    %## PRINT TABLES
                    t = sprintf_table(anova_out);
                    t.saveToFile([tmp_savedir filesep sprintf('cl%s_%s_%s_ANOVA.txt',string(clusters(cl_i)),varnames{var_i},measure_name)]);
                    t = sprintf_table(stats_out.Coefficients);
                    t.saveToFile([tmp_savedir filesep sprintf('cl%s_%s_%s_LM.txt',string(clusters(cl_i)),varnames{var_i},measure_name)]);
                    %-
                    anova_p_cond = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    anova_p_inter = anova_out.pValue(strcmp(anova_out.Properties.RowNames,sprintf('cond_char:%s',varnames{var_i})));
                    %-
                    pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                    pval_0p5 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
                    pval_0p75 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
                    pval_1p0 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    pval_0p5_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_0.5:%s',varnames{var_i})));
                    pval_0p75_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_0.75:%s',varnames{var_i})));
                    pval_1p0_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_1.0:%s',varnames{var_i})));
                    %-
                    slope_inter = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                    slope_0p5 = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
                    slope_0p75 = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
                    slope_1p0 = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
                    slope_var = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_0p5_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_0.5:%s',varnames{var_i}))));
                    slope_0p75_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_0.75:%s',varnames{var_i}))));
                    slope_1p0_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('cond_char_1.0:%s',varnames{var_i}))));
                catch e
                    fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
                    R2 = 0;
                    pval = 1;
                    slope = 0;
                    inter = 0;
                end
                %## SCATTER
                axes();
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    inds = T_vals_plot.cond_char==loc_cond_chars(cond_i);
                    data = T_vals_plot(inds,:);
                    [vals,inds] = sort(data.(varnames{var_i}));
                    data = data(inds,:);
                    ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',GROUP_SHORTS{group_i}));
                    ss.CData = color_dark(cond_i,:);
                    ss.SizeData = 15;
                    ss.Marker = GROUP_MARKS{group_i};
                    if meas_i == 1
                        cond_plot_store = [cond_plot_store, ss];
                    end
                end
                %## LINEAR MODEL FIT
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    switch cond_i
                        case 1
                            c2 = 0;
                            c3 = 0;
                            c4 = 0;
                        case 2
                            c2 = 1;
                            c3 = 0;
                            c4 = 0;
                        case 3
                            c2 = 0;
                            c3 = 1;
                            c4 = 0;
                        case 4
                            c2 = 0;
                            c3 = 0;
                            c4 = 1;
                    end
                    inds = T_vals_plot.cond_char==loc_cond_chars(cond_i);
                    data = T_vals_plot(inds,:);
                    [vals,inds] = sort(data.(varnames{var_i}));
                    data = data(inds,:);
                    if anova_p_var < 0.1 || anova_p_cond < 0.1 || anova_p_inter
                        x = unique(data.(varnames{var_i}));
                        y = [];
                        for i = 1:length(x)
                            y(i) = x(i)*slope_var + slope_0p5*c2 + slope_0p75*c3 + slope_1p0*c4 + slope_0p5_var*x(i)*c2 + slope_0p75_var*x(i)*c3 + slope_1p0_var*x(i)*c4 + slope_inter ;
                        end
                        pp = plot(x,y,...
                            'DisplayName',sprintf('p_{%s}=(%0.2f,%0.2f,%0.2f)',MEASURE_NAME_LABS{meas_i},anova_p_var,anova_p_cond,anova_p_inter),...
                            'LineWidth',2);
                        pp.LineStyle = GROUP_LINESTYLES{group_i};
                        pp.Color = color_light(cond_i,:);
                        if cond_i == 1
                            % eq = sprintf('y=(%0.1g)*x+(%0.1g)*%0.2f+(%0.1g)*x+(%0.1g)',slope_var,slope_cnd,loc_cond_chars(cond_i),slope_grp,inter_mn);
                            % eq = sprintf('y=(%0.1g)*x+(%0.1g)*c_i\n%6s+(%0.1g)',slope_var,slope_cnd,'',inter_mn);
                            eq = '';
                            x_txt = 0.15;
                            y_txt = 0.7;
                            if anova_p_var > 0.01 & anova_p_var < 0.05
                                text(x_txt,y_txt,sprintf('* %s\nR^2=%0.2g',eq,R2),...
                                    'FontSize',REGRESS_TXT_SIZE,...
                                    'FontName',PLOT_STRUCT.font_name,...
                                    'FontWeight','bold','Units','normalized');
                            elseif anova_p_var <= 0.01 & anova_p_var > 0.001 
                                text(x_txt,y_txt,sprintf('** %s\nR^2=%0.2g',eq,R2),...
                                    'FontSize',REGRESS_TXT_SIZE,...
                                    'FontName',PLOT_STRUCT.font_name,...
                                    'FontWeight','bold','Units','normalized');
                            else
                                text(x_txt,y_txt,sprintf('*** %s\nR^2=%0.2g',eq,R2),...
                                    'FontSize',REGRESS_TXT_SIZE,...
                                    'FontName',PLOT_STRUCT.font_name,...
                                    'FontWeight','bold','Units','normalized');
                            end
                            stats_store = [stats_store, pp];
                        end
                    end
                end
                ylabel('10*log_{10}(Flattened PSD)');
                xlabel(varnames_labs{var_i});
                title(MEASURE_NAME_LABS{meas_i});
                set(gca,'FontWeight','bold');
                ylim(y_lim_calc)
                %## legend
                if meas_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',double(string(loc_cond_chars(cnt))));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName','Arial','FontSize',9)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Position',[0.1,0.84,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = 18;
                elseif meas_i == 2
                    % %- lg1
                    % legend(gca,group_plot_store);
                    % [lg1,icons,plots,txt] = legend('boxoff');
                    % set(lg1,'Orientation','horizontal')
                    % set(lg1,'FontName','Arial','FontSize',9)
                    % % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift+0.025,lg1.Position(3),lg1.Position(4)]);
                    % lg1.ItemTokenSize(1) = 18;
                elseif meas_i == 3
                    %- lg3
                    if ~isempty(stats_store)
                        legend(gca,stats_store);
                        [lg3,~,~,~] = legend('boxoff');
                        set(lg3,'Orientation','horizontal')
                        set(lg3,'FontName','Arial','FontSize',9)
                        % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                        set(lg3,'Position',[0.1,0.775,lg3.Position(3),lg3.Position(4)]);
                        lg3.ItemTokenSize(1) = 18;
                    end
                end
                set(gca,'Position',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
                horiz_shift = horiz_shift + AX_W*im_resize + 0.05;
            end
            %## TITLE
            % annotation('textbox',[0.5-0.1,SCATTER_BOTTOM-vert_shift-0.05+AX_H*im_resize,0.2,0.2],...
            %     'String',string(design_chars(des_i)),'HorizontalAlignment','center',...
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            % vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_eeg-kin-speed-inter.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        % close(fig)
    end
end

%% PREDICTORS: SPEED Factor, GROUP, INTERACTION. RESPONSE: BRAIN ACTIVITY, STATS TEST
im_resize = 0.8;
AX_W = 0.35;
AX_H = 0.25*3;
DO_PLOT_GROUPS = true;
VIOLIN_BOTTOM = 0.2;
% GROUP_SHORTS = {'HO','FO'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
tmp_savedir = [save_dir filesep 'Pspeedf-PGroup-Pinter-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    atlas_name = atlas_name_store{cl_i};
    fig = figure('color','white','renderer','Painters');
    sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
    set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    hold on;
    set(gca,AXES_DEFAULT_PROPS{:})
    vert_shift = 0;
    for des_i = 2
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
        horiz_shift = 0;
        stats_store = [];
        for meas_i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{meas_i};
            %##
            cond_plot_store = [];
            group_plot_store = [];
            % inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(group_i);
            % T_stats_plot = psd_feature_stats(inds,:);
            inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
            T_vals_plot = TMP_FOOOF_T(inds,:);
            % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
            loc_cond_chars = unique(T_vals_plot.cond_char);
            y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
            T_vals_plot = table(double(T_vals_plot.(measure_name)),categorical(string(T_vals_plot.cond_char)),categorical(string(T_vals_plot.group_char)),...
                'VariableNames',{measure_name,'cond_char','group_char'});
            try
                mod = sprintf('%s ~ 1 + cond_char + group_char + cond_char*group_char',measure_name);
                % stats_out = fitlme(T_vals_plot,mod);
                stats_out = fitlm(T_vals_plot,mod);
                % anova_out = anova(stats_out);
                out = anova(T_vals_plot,mod,'SumOfSquaresType',"three");
                anova_out = out.stats();
                % anova_out = anovan(double(T_vals_plot.(measure_name)),{T_vals_plot.group_char,T_vals_plot.cond_char},...
                %         'model','interaction',...
                %         'model',2,...
                %         'sstype',3,...
                %         'varnames',strvcat('group','speed'));
                R2 = stats_out.Rsquared.Adjusted;
                %## PRINT TABLES
                t = sprintf_table(anova_out);
                t.saveToFile([tmp_savedir filesep sprintf('cl%s_%s_ANOVA.txt',string(clusters(cl_i)),measure_name)]);
                t = sprintf_table(stats_out.Coefficients);
                t.saveToFile([tmp_savedir filesep sprintf('cl%s_%s_LM.txt',string(clusters(cl_i)),measure_name)]);
                %-
                anova_p_cond = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
                anova_p_group = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'group_char'));
                anova_p_inter = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char:group_char'));
                
                pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                pval_0p5 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5'));
                pval_0p75 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75'));
                pval_1p0 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0'));
                pval_grp2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H2000''s'));
                pval_grp3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s'));
                %-
                pval_var_0p5_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5:group_char_H2000''s'));
                pval_var_0p75_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75:group_char_H2000''s'));
                pval_var_1p0_g2 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0:group_char_H2000''s'));
                %-
                pval_var_0p5_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.5:group_char_H3000''s'));
                pval_var_0p75_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_0.75:group_char_H3000''s'));
                pval_var_1p0_g3 = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char_1.0:group_char_H3000''s'));
            catch e
                fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
                R2 = 0;
                pval = 1;
                slope = 0;
                inter = 0;
            end
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
            ap = 1;
            if anova_p_group < 0.05
                ap = anova_p_group;
            end
            if anova_p_var < 0.05
                ap = anova_p_var;
            end
            if anova_p_inter < 0.05
                ap = anova_p_inter;
            end
            STATS_STRUCT.anova_grp = {anova_p_group,anova_p_group,anova_p_group};
            STATS_STRUCT.pvals_grp_pairs = {[1,1],[1,2],[1,3]};
            STATS_STRUCT.pvals_grp = {1,pval_grp2,pval_grp3};
            STATS_STRUCT.anova={anova_p_var,anova_p_var,anova_p_var};
            STATS_STRUCT.pvals = {[1,pval_0p5,pval_0p75,pval_1p0],...
                [1,pval_var_0p5_g2,pval_var_0p75_g2,pval_var_1p0_g2],...
                [1,pval_var_0p5_g3,pval_var_0p75_g3,pval_var_1p0_g3]};
            STATS_STRUCT.pvals_pairs = {{[1,1],[1,2],[1,3],[1,4]},...
                {[1,1],[1,2],[1,3],[1,4]},...
                {[1,1],[1,2],[1,3],[1,4]}};
            % STATS_STRUCT.anova={ap,ap,ap};
            % STATS_STRUCT.pvals = {[1,pval_var_0p5,pval_var_0p75,pval_var_1p0],...
            %     [1,pval_var_0p5_g2,pval_var_0p75_g2,pval_var_1p0_g2],...
            %     [1,pval_var_0p5_g3,pval_var_0p75_g3,pval_var_1p0_g3]};
            % STATS_STRUCT.pvals_pairs = {{[1,1],[1,2],[1,3],[1,4]},...
            %     {[1,1],[1,2],[1,3],[1,4]},...
            %     {[1,1],[1,2],[1,3],[1,4]}};
            %-
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_STRUCT = struct('color_map',color_dark,...
                'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label','10*log_{10}(Flattened PSD)','group_offsets',[0.125,0.475,0.812],...
                'title',MEASURE_NAME_LABS{meas_i},'font_size',10,'ylim',[min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+3*std(T_vals_plot.(measure_name))],...
                'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);
            axax = group_violin(T_vals_plot,measure_name,'cond_char','group_char',...
                fig,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_STRUCT,...
                'STATS_STRUCT',STATS_STRUCT);
            if meas_i > 1
                ylabel(axax,'');
            end
            set(axax,'OuterPosition',[0,0,1,1]);
            set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_H*im_resize + 0.1;
        end
        %## TITLE
        vert_shift = vert_shift + AX_H*im_resize+0.1;
    end
     hold off;
    %##
    exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_eeg-speed-group-interact.tiff',string(clusters(cl_i)))],'Resolution',300)
    close(fig)
end

%% PREDICTORS: KINEMATICS, GROUP, & INTERACTION. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
REGRESS_TXT_XMULTI = 0.9;
REGRESS_TXT_YMULTI = 1.0;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25;
GROUP_SHORTS = {'YA','HO','FO'};
GROUP_MARKS = {'o','x','^'};
GROUP_LINESTYLES = {'-','-.',':'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
PLOT_STRUCT = struct('color_map',[],...
        'cond_labels',unique(TMP_FOOOF_T.cond_char),'group_labels',unique(TMP_FOOOF_T.group_char),...
        'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',[],...
        'title',[],'font_size',10,'ylim',[],...
        'font_name','Arial','x_label','speed','do_combine_groups',true);
tmp_savedir = [save_dir filesep 'Pkin-PGroup-Pinter-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %%
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white','renderer','Painters');
        sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        for des_i = 2
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
            horiz_shift = 0;
            stats_store = [];
            cond_plot_store = [];
            for meas_i = 1:length(measure_name_plot)
                measure_name = measure_name_plot{meas_i};
                %##
                group_plot_store = [];
                % inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(group_i);
                % T_stats_plot = psd_feature_stats(inds,:);
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
                T_vals_plot = TMP_FOOOF_T(inds,:);
                %- continuous or categorical
                % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                %-
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                x_txt = min(T_vals_plot.(varnames{var_i}))*REGRESS_TXT_XMULTI+std(T_vals_plot.(varnames{var_i}));
                y_txt = max(T_vals_plot.(measure_name))*REGRESS_TXT_YMULTI+std(T_vals_plot.(measure_name));
                try
                    mod = sprintf('%s ~ 1 + %s + group_char + %s*group_char',measure_name,varnames{var_i},varnames{var_i});
                    stats_out = fitlm(T_vals_plot,mod);
                    % anova_out = anova(stats_out);
                    anova_out = anova(T_vals_plot,mod,'SumOfSquaresType',"three",'CategoricalFactors',{'group_char'},...
                        'ModelSpecification','linear');
                    anova_out = anova_out.stats();
                    % anova_out = anovan(double(T_vals_plot.(measure_name)),{T_vals_plot.group_char,double(T_vals_plot.(varnames{var_i}))},...
                    %     'continuous',[2],...
                    %     'model','interaction',...
                    %     'model',2,...
                    %     'sstype',3,...
                    %     'varnames',strvcat('group',varnames{var_i}));
                
                    R2 = stats_out.Rsquared.Adjusted;
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    anova_p_grp = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'group_char'));
                    anova_p_inter = anova_out.pValue(strcmp(anova_out.Properties.RowNames,sprintf('group_char:%s',varnames{var_i})));
                    
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    pval_grp2 = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H2000''s')));
                    pval_grp3 = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
                    pval_var_grp2 = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('group_char_H2000''s:%s',varnames{var_i}))));
                    pval_var_grp3 = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('group_char_H3000''s:%s',varnames{var_i}))));
                    
                    slope_var = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_grp2 = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H2000''s')));
                    slope_grp3 = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
                    slope_var_grp2 = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('group_char_H2000''s:%s',varnames{var_i}))));
                    slope_var_grp3 = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,sprintf('group_char_H3000''s:%s',varnames{var_i}))));
                    inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                    %## FPRINTF
                    fid = fopen([tmp_savedir filesep sprintf('cl%s_%s_eeg-kin-group-interact_ANOVA.txt',string(clusters(cl_i)),varnames{var_i})],'wt');
                    %- converst Summary anova to char array
                    txt = evalc('anova_out');
                    fprintf(fid,'ANOVA EVAL\n');
                    fprintf(fid,'%s',txt);
                    fprintf(fid,'\n');
                    %- Convert summary to char array
                    fprintf(fid,'FITGLME EVAL\n');
                    txt = evalc('stats_out');
                    fprintf(fid,'%s',txt);
                    fclose(fid);
                catch e
                    fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
                    R2 = 0;
                    pval = 1;
                    slope = 0;
                    inter = 0;
                end
                %## SCATTER
                axes();
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',GROUP_SHORTS{group_i}));
                        ss.CData = color_dark(cond_i,:);
                        ss.SizeData = 15;
                        ss.Marker = GROUP_MARKS{group_i};
                        if group_i == 1 && meas_i == 1
                            cond_plot_store = [cond_plot_store, ss];
                        end
                        if cond_i == 1
                            group_plot_store = [group_plot_store, ss];
                        end
                    end
                end
                %## LINEAR MODEL FIT
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        switch group_i
                            case 1
                                g2 = 0;
                                g3 = 0;
                            case 2
                                g2 = 1;
                                g3 = 0;
                            case 3
                                g2 = 0;
                                g3 = 1;
                        end
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        if anova_p_var < 0.1 || anova_p_grp < 0.1
                            x = unique(data.(varnames{var_i}));
                            y = [];
                            for i = 1:length(x)
                                % y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*x(i) + inter_mn;
                                y(i) = x(i)*slope_var + slope_grp2*g2 + slope_grp3*g3 + slope_var_grp2*x(i)*g2 + slope_var_grp3*x(i)*g3 + inter_mn ;
                            end
                            pp = plot(x,y,...
                                'DisplayName',sprintf('p-vals=(%0.2f,%0.2f,%0.2f)',anova_p_var,anova_p_grp,anova_p_inter),...
                                'LineWidth',2);
                            pp.LineStyle = GROUP_LINESTYLES{group_i};
                            pp.Color = color_light(cond_i,:);
                            if cond_i == 1
                                % eq = sprintf('y=(%0.1g)*x+(%0.1g)*%0.2f+(%0.1g)*x+(%0.1g)',slope_var,slope_cnd,loc_cond_chars(cond_i),slope_grp,inter_mn);
                                % eq = sprintf('y=(%0.1g)*x+(%0.1g)*c_i\n%6s+(%0.1g)',slope_var,slope_cnd,'',inter_mn);
                                eq = '';
                                if anova_p_var > 0.01 & anova_p_var < 0.05
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif anova_p_var <= 0.01 & anova_p_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                end
                                if group_i == 1
                                    stats_store = [stats_store, pp];
                                end
                            end
                        end
                    end
                end
                ylabel('10*log_{10}(Flattened PSD)');
                xlabel(varnames_labs{var_i});
                title(MEASURE_NAME_LABS{meas_i});
                set(gca,'FontWeight','bold');
                ylim(y_lim_calc)
                %## legend
                if meas_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',double(string(loc_cond_chars(cnt))));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName','Arial','FontSize',9)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Position',[0.1-0.027,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = 18;
                elseif meas_i == 2
                    %- lg1
                    legend(gca,group_plot_store);
                    [lg1,icons,plots,txt] = legend('boxoff');
                    set(lg1,'Orientation','horizontal')
                    set(lg1,'FontName','Arial','FontSize',9)
                    % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift+0.025,lg1.Position(3),lg1.Position(4)]);
                    lg1.ItemTokenSize(1) = 18;
                elseif meas_i == 3
                    %- lg3
                    if ~isempty(stats_store)
                        legend(gca,stats_store);
                        [lg3,~,~,~] = legend('boxoff');
                        set(lg3,'Orientation','horizontal')
                        set(lg3,'FontName','Arial','FontSize',9)
                        % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                        set(lg3,'Position',[0.1+0.01,SCATTER_BOTTOM+AX_W*im_resize-vert_shift-0.03,lg3.Position(3),lg3.Position(4)]);
                        lg3.ItemTokenSize(1) = 18;
                    end
                end
                set(gca,'Position',[horiz_shift+0.05,SCATTER_BOTTOM-vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
                horiz_shift = horiz_shift + AX_H*im_resize + 0.125;
            end
            %## TITLE
            % annotation('textbox',[0.5-0.1,SCATTER_BOTTOM-vert_shift-0.05+AX_H*im_resize,0.2,0.2],...
            %     'String',string(design_chars(des_i)),'HorizontalAlignment','center',...
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_%s_eeg-kin-group-inter.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end

%% PREDICTORS: KINEMATICS, GROUP. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
REGRESS_TXT_XMULTI = 0.9;
REGRESS_TXT_YMULTI = 1.0;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25;
GROUP_SHORTS = {'HO','FO'};
PLOT_STRUCT = struct('color_map',color_dark,...
        'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
        'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames{var_i},...
        'title',varnames{var_i},'font_size',10,'ylim',[],...
        'font_name','Arial','x_label','speed','do_combine_groups',true);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %##
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white','renderer','Painters');
        sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        for des_i = 2
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
            horiz_shift = 0;
            stats_store = [];
            cond_plot_store = [];
            for meas_i = 1:length(measure_name_plot)
                measure_name = measure_name_plot{meas_i};
                %##
                group_plot_store = [];
                % inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(group_i);
                % T_stats_plot = psd_feature_stats(inds,:);
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
                T_vals_plot = TMP_FOOOF_T(inds,:);
                T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                x_txt = min(T_vals_plot.(varnames{var_i}))*REGRESS_TXT_XMULTI+std(T_vals_plot.(varnames{var_i}));
                y_txt = max(T_vals_plot.(measure_name))*REGRESS_TXT_YMULTI+std(T_vals_plot.(measure_name));
                try
                    mod = sprintf('%s ~ 1 + %s + group_char',measure_name,varnames{var_i});
                    stats_out = fitlm(T_vals_plot,mod);
                    anova_out = anova(stats_out);
                    R2 = stats_out.Rsquared.Adjusted;
                    % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,varnames{var_i}));
                    % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,varnames{var_i})));
                    % inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    pval_grp = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s'));
                    tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i})));
                    slope_grp = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
                    % slope_cnd = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
                    inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                catch e
                    fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
                    R2 = 0;
                    pval = 1;
                    slope = 0;
                    inter = 0;
                end
                %##
                axes();
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',GROUP_SHORTS{group_i}));
                        if group_i == 1
                            ss.Marker = 'o';
                            ss.CData = color_light(cond_i,:);
                            ss.SizeData = 15;
                            if meas_i == 1
                                cond_plot_store = [cond_plot_store,ss];
                            end
                        else
                            ss.Marker = 'x';
                            ss.CData = color_light(cond_i,:)+GROUP_CMAP_OFFSET;
                            ss.SizeData = 15;
                            % ss.MarkerSize = 5;
                        end
                        if cond_i == 1
                            group_plot_store = [group_plot_store, ss];
                        end
                    end
                end
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        if pval_var < 0.1 || pval_grp < 0.1
                            x = unique(data.(varnames{var_i}));
                            y = [];
                            for i = 1:length(x)
                                % y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*x(i) + inter_mn;
                                y(i) = x(i)*slope_var + slope_grp*(group_i-1) + inter_mn;
                            end
                            pp = plot(x,y,...
                                'DisplayName',sprintf('p-vals=(%0.2f,%0.2f)',GROUP_SHORTS{group_i},pval_var,pval_grp),...
                                'LineWidth',2);
                            if group_i == 1
                                pp.LineStyle = '-';
                                pp.Color = color_dark(cond_i,:);
                            else
                                pp.LineStyle = '-.';
                                pp.Color = color_dark(cond_i,:)+GROUP_CMAP_OFFSET;
                            end
                            if cond_i == 1
                                % eq = sprintf('y=(%0.1g)*x+(%0.1g)*%0.2f+(%0.1g)*x+(%0.1g)',slope_var,slope_cnd,loc_cond_chars(cond_i),slope_grp,inter_mn);
                                eq = sprintf('y=(%0.1g)*x+(%0.1g)*g_i\n%6s+(%0.1g)',slope_var,slope_grp,'',inter_mn);
                                if pval_var > 0.01 & pval_var < 0.05
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                end
                                if group_i == 1
                                    stats_store = [stats_store, pp];
                                end
                            end
                        end
                    end
                end
                ylabel('10*log_{10}(Flattened PSD)');
                xlabel(varnames_labs{var_i});
                title(MEASURE_NAME_LABS{meas_i});
                set(gca,'FontWeight','bold');
                ylim(y_lim_calc)
                %## legend
                if meas_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',loc_cond_chars(cnt));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName','Arial','FontSize',9)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Position',[0.1-0.027,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = 18;
                elseif meas_i == 2
                    %- lg1
                    legend(gca,group_plot_store);
                    [lg1,icons,plots,txt] = legend('boxoff');
                    set(lg1,'Orientation','horizontal')
                    set(lg1,'FontName','Arial','FontSize',9)
                    % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift+0.025,lg1.Position(3),lg1.Position(4)]);
                    lg1.ItemTokenSize(1) = 18;
                elseif meas_i == 3
                    %- lg3
                    if ~isempty(stats_store)
                        legend(gca,stats_store);
                        [lg3,~,~,~] = legend('boxoff');
                        set(lg3,'Orientation','horizontal')
                        set(lg3,'FontName','Arial','FontSize',9)
                        % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                        set(lg3,'Position',[0.1+0.01,SCATTER_BOTTOM+AX_W*im_resize-vert_shift-0.03,lg3.Position(3),lg3.Position(4)]);
                        lg3.ItemTokenSize(1) = 18;
                    end
                end
                set(gca,'Position',[horiz_shift+0.05,SCATTER_BOTTOM-vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
                horiz_shift = horiz_shift + AX_H*im_resize + 0.125;
                % if meas_i == 1
                %     ylabel(measure_name,'Interpreter','none');
                % else
                %     ylabel('');
                % end
            end
            %## TITLE
            % annotation('textbox',[0.5-0.1,SCATTER_BOTTOM-vert_shift-0.05+AX_H*im_resize,0.2,0.2],...
            %     'String',string(design_chars(des_i)),'HorizontalAlignment','center',...
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[kin_savedir filesep sprintf('cl%s_%s_group_kin-eeg.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end

%% PREDICTORS: KINEMATICS, GROUP, COND. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
REGRESS_TXT_XMULTI = 0.9;
REGRESS_TXT_YMULTI = 1.0;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25;
GROUP_SHORTS = {'HO','FO'};
PLOT_STRUCT = struct('color_map',color_dark,...
        'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
        'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames{var_i},...
        'title',varnames{var_i},'font_size',10,'ylim',[],...
        'font_name','Arial','x_label','speed','do_combine_groups',true);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %##
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white','renderer','Painters');
        sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        for des_i = 2
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
            horiz_shift = 0;
            stats_store = [];
            cond_plot_store = [];
            for meas_i = 1:length(measure_name_plot)
                measure_name = measure_name_plot{meas_i};
                %##
                group_plot_store = [];
                % inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(group_i);
                % T_stats_plot = psd_feature_stats(inds,:);
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
                T_vals_plot = TMP_FOOOF_T(inds,:);
                T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                x_txt = min(T_vals_plot.(varnames{var_i}))*REGRESS_TXT_XMULTI+std(T_vals_plot.(varnames{var_i}));
                y_txt = max(T_vals_plot.(measure_name))*REGRESS_TXT_YMULTI+std(T_vals_plot.(measure_name));
                try
                    mod = sprintf('%s ~ 1 + %s + group_char + cond_char',measure_name,varnames{var_i});
                    stats_out = fitlm(T_vals_plot,mod);
                    anova_out = anova(stats_out);
                    R2 = stats_out.Rsquared.Adjusted;
                    % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,varnames{var_i}));
                    % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,varnames{var_i})));
                    % inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    pval_grp = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
                    pval_cnd = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
                    tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i})));
                    slope_grp = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
                    slope_cnd = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
                    inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
                catch e
                    fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
                    R2 = 0;
                    pval = 1;
                    slope = 0;
                    inter = 0;
                end
                %##
                axes();
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',GROUP_SHORTS{group_i}));
                        if group_i == 1
                            ss.Marker = 'o';
                            ss.CData = color_light(cond_i,:);
                            ss.SizeData = 15;
                            if meas_i == 1 
                                cond_plot_store = [cond_plot_store, ss];
                            end
                        else
                            ss.Marker = 'x';
                            ss.CData = color_light(cond_i,:)+GROUP_CMAP_OFFSET;
                            ss.SizeData = 15;
                        end
                        if cond_i == 1
                            group_plot_store = [group_plot_store, ss];
                        end
                    end
                end
                hold on;
                for cond_i = 1:length(loc_cond_chars)
                    for group_i = 1:length(groups)
                        inds = T_vals_plot.cond_char==loc_cond_chars(cond_i) & T_vals_plot.group_id==string(group_i);
                        data = T_vals_plot(inds,:);
                        [vals,inds] = sort(data.(varnames{var_i}));
                        data = data(inds,:);
                        if pval_var < 0.1 || pval_grp < 0.1
                            x = unique(data.(varnames{var_i}));
                            y = [];
                            for i = 1:length(x)
                                % y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*x(i) + inter_mn;
                                y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*(group_i-1) + inter_mn;
                            end
                            pp = plot(x,y,...
                                'DisplayName',sprintf('p-vals=\n%6s(%0.2f,%0.2f,%0.2f)',GROUP_SHORTS{group_i},'',pval_var,pval_grp,pval_cnd),...
                                'LineWidth',2);
                            if group_i == 1
                                pp.LineStyle = '-';
                                pp.Color = color_dark(cond_i,:)+GROUP_CMAP_OFFSET;
                            else
                                pp.LineStyle = '-.';
                                pp.Color = color_dark(cond_i,:)+GROUP_CMAP_OFFSET;
                            end
                            if cond_i == 1
                                % eq = sprintf('y=(%0.1g)*x+(%0.1g)*%0.2f+(%0.1g)*x+(%0.1g)',slope_var,slope_cnd,loc_cond_chars(cond_i),slope_grp,inter_mn);
                                eq = sprintf('y=(%0.1g)*x+(%0.1g)*c_i\n%6s+(%0.1g)*g_i+(%0.1g)',slope_var,slope_cnd,'',slope_grp,inter_mn);
                                if pval_var > 0.01 & pval_var < 0.05
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                end
                                if group_i == 1
                                    stats_store = [stats_store, pp];
                                end
                            end
                        end
                    end
                end
                ylabel('10*log_{10}(Flattened PSD)');
                xlabel(varnames_labs{var_i});
                title(MEASURE_NAME_LABS{meas_i});
                set(gca,'FontWeight','bold');
                ylim(y_lim_calc)
                %## legend
                if meas_i == 1
                    %- lg2
                    legend(gca,cond_plot_store);
                    [lg2,icons,plots,txt]  = legend('boxoff');
                    tmp = get(lg2,'String');
                    cnt = 1;
                    for i = 1:length(cond_plot_store)
                        tmp{i} = sprintf('%0.2g',loc_cond_chars(cnt));
                        cnt = cnt + 1;
                    end
                    set(lg2,'String',tmp,'FontName','Arial','FontSize',9)
                    set(lg2,'Orientation','horizontal')
                    set(lg2,'Position',[0.1-0.027,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg2.Position(3),lg2.Position(4)]);
                    lg2.ItemTokenSize(1) = 18;
                elseif meas_i == 2
                    %- lg1
                    legend(gca,group_plot_store);
                    [lg1,icons,plots,txt] = legend('boxoff');
                    set(lg1,'Orientation','horizontal')
                    set(lg1,'FontName','Arial','FontSize',9)
                    % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift+0.025,lg1.Position(3),lg1.Position(4)]);
                    lg1.ItemTokenSize(1) = 18;
                elseif meas_i == 3
                    %- lg3
                    if ~isempty(stats_store)
                        legend(gca,stats_store);
                        [lg3,~,~,~] = legend('boxoff');
                        set(lg3,'Orientation','horizontal')
                        set(lg3,'FontName','Arial','FontSize',9)
                        % set(lg1,'Position',[0.1,SCATTER_BOTTOM+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                        set(lg3,'Position',[0.1+0.01,SCATTER_BOTTOM+AX_W*im_resize-vert_shift-0.05,lg3.Position(3),lg3.Position(4)]);
                        lg3.ItemTokenSize(1) = 18;
                    end
                end
                set(gca,'Position',[horiz_shift+0.05,SCATTER_BOTTOM-vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
                horiz_shift = horiz_shift + AX_H*im_resize + 0.125;
            end
            %## TITLE
            % annotation('textbox',[0.5-0.1,SCATTER_BOTTOM-vert_shift-0.05+AX_H*im_resize,0.2,0.2],...
            %     'String',string(design_chars(des_i)),'HorizontalAlignment','center',...
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_STRUCT.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[kin_savedir filesep sprintf('cl%s_%s_group_cond_kin-eeg.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end
