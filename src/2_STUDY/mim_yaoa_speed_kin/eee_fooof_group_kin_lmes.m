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
% measure_name_plot = {'theta_avg_power','alpha_avg_power','beta_avg_power','tbr_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% title_plot = {'Mean \theta','Mean \alpha','Mean \beta','\beta//\theta'};
MEASURES_ANALYZE = {'theta_avg_power','beta_avg_power','tbr_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'Mean \theta','Mean \beta','\beta//\theta'};
%% ===================================================================== %%
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
%% ================================================================== %%
%## CREATE EEG-KINEMATICS TABLE
%{
%## LOAD EEG & KINEMATIC DATA
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
xlsx_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'raw_data_vis'];
imu_table = [xlsx_fpath filesep 'imu_table_meantrial.xlsx'];
ls_table = [xlsx_fpath filesep 'ls_table_meantrial.xlsx'];
% imu_table = [xlsx_fpath filesep 'imu_table_out.xlsx'];
% ls_table = [xlsx_fpath filesep 'ls_table_out.xlsx'];
imu_table = readtable(imu_table);
ls_table = readtable(ls_table);
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
writetable(FOOOF_TABLE,[save_dir filesep 'fooof_kinematics_table_nans.xlsx'])
par_save(FOOOF_TABLE,save_dir,'fooof_kinematics_table_nans.mat');
%}
%% ===================================================================== %%
% FOOOF_TABLE = readtable([save_dir filesep 'fooof_kinematics_table.xlsx']);
% FOOOF_TABLE = par_load(save_dir,'fooof_kinematics_table.mat');
% FOOOF_TABLE = readtable([save_dir filesep 'fooof_kinematics_table_nans.xlsx']);
FOOOF_TABLE = par_load(save_dir,'fooof_kinematics_table_nans.mat');
%- calculate theta-beta ratios
FOOOF_TABLE.btr_avg_power = FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power;
FOOOF_TABLE.tbr_avg_power = FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power;
FOOOF_TABLE.logbtr_avg_power = real(log10(FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power));
FOOOF_TABLE.logtbr_avg_power = real(log10(FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power));
% FOOOF_TABLE.logbtr_norm_avg_power = log10(FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power+max(abs(FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power)));
% FOOOF_TABLE.logtbr_norm_avg_power = log10(FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power+max(abs(FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power)));
% FOOOF_TABLE.logbtr_norm_avg_power = log10((FOOOF_TABLE.beta_avg_power+abs(min(FOOOF_TABLE.beta_avg_power)))./(FOOOF_TABLE.theta_avg_power+abs(min(FOOOF_TABLE.theta_avg_power))));
% FOOOF_TABLE.logtbr_norm_avg_power = log10((FOOOF_TABLE.theta_avg_power+abs(min(FOOOF_TABLE.theta_avg_power)))./(FOOOF_TABLE.beta_avg_power+abs(min(FOOOF_TABLE.beta_avg_power))));
FOOOF_TABLE.logbtr_norm_avg_power = real(log10(FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power).*conj(log10(FOOOF_TABLE.beta_avg_power./FOOOF_TABLE.theta_avg_power)));
FOOOF_TABLE.logtbr_norm_avg_power = real(log10(FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power).*conj(log10(FOOOF_TABLE.theta_avg_power./FOOOF_TABLE.beta_avg_power)));
TMP_FOOOF_T = FOOOF_TABLE;

%% STATISTICS TRACKING STRUCTURE
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'resp_terms',{{''}},...
    'pred_terms',{{''}},...
    'rnd_terms',{{''}},...
    'anova_grp_p',[],...
    'anova_kin_p',[],...
    'anova_eeg_p',[],...
    'anova_speed_p',[],...
    'anova_intact_p',[],...
    'anova_inter_p',[],...
    'lme_grp_p',[],...
    'lme_kin_p',[],...
    'lme_eeg_p',[],...
    'lme_speed_p',[],...
    'lme_intact_p',[],...
    'lme_inter_p',[],...
    'lme_grp_coeff',[],...
    'lme_kin_coeff',[],...
    'lme_eeg_coeff',[],...
    'lme_speed_coeff',[],...
    'lme_intact_coeff',[],...
    'lme_inter_coeff',[],...
    'lme_rnd_effects',{{}},...
    'R2',[],...
    'norm_test_h',[],...
    'norm_test_p',[]);
STATS_TRACK_STRUCT = DEF_STATS_TRACK_STRUCT;
cntsts = 1;
%## MEASURES TO COMPUTE (LIMIT TO 3 WHEN POSSIBLE)
%- measures: theta_avg_power, alpha_avg_power, theta_avg_power,
%   tbr_avg_power, btr_avg_power, logbtr_avg_power, logtbr_avg_power
% (NOTE): logs are real(log(measure)), this may be improved in some manner
%latter. 
MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
MEASURE_SYMBOLS = {'\theta','\alpha','\beta'};
%-
% MEASURES_ANALYZE = {'theta_avg_power','beta_avg_power','logtbr_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% MEASURE_NAME_LABS = {'Mean \theta','Mean \beta','real(log_{10}(\theta//\beta))'};
%-
% MEASURES_ANALYZE = {'theta_avg_power','tbr_avg_power','logtbr_norm_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% MEASURE_NAME_LABS = {'Mean \theta','(\theta/\beta)','conj(log_{10}(\theta/\beta))'};
%-
% MEASURES_ANALYZE = {'theta_avg_power','alpha_avg_power','beta_avg_power','tbr_avg_power',...
%     'logtbr_avg_power','logtbr_norm_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
% MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta','\theta/\beta',...
%     'Re(log(\theta/\beta))','conj(log_{10}(\theta/\beta))'};
%## HARDCODED VARIABLES OF INTEREST
% varnames = {'mean_APexc_COV','mean_APexc_mean','mean_MLexc_COV',...
%             'mean_MLexc_mean','mean_StepDur','mean_UDexc_COV',...
%             'mean_UDexc_mean','mean_StanceDur','mean_GaitCycleDur','mean_PeakUpDownVel_mean'};
% varnames_labs = {'APexc COV','APexc','MLexc COV','MLexc','Step Dur','UDexc COV','UDexc','Stance Dur','GaitCycle Dur','Peak UpDown Vel'};
%-
varnames = {'mean_APexc_COV','mean_MLexc_COV',...
            'mean_StepDur_cov','mean_StepDur','mean_UDexc_COV',...
            'mean_StanceDur'};
varnames_labs = {'AP Exc. COV','ML Exc. COV','Step Duration COV','Step Duration (s)','UD Exc. COV','Stance Duration (s)'};
% varnames = {'mean_StepDur','mean_GaitCycleDur','mean_SwingDur',...
%             'mean_StanceDur','mean_SingleSupport','mean_TotalDS'};
% varnames_labs = {'mean_StepDur','mean_GaitCycleDur','mean_SwingDur',...
%             'mean_StanceDur','mean_SingleSupport','mean_TotalDS'};
%## GROUP VARIABLE?
group_i = 1;
%% (SIMPLE MEDIATION MODEL) PREDICTORS: KINEMATIC. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
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
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = false;
REG_TXT_SIZE = 8;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.125;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep 'lme_Pkin-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
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
        stats_store = [];
        for meas_i = 1:length(MEASURES_ANALYZE)
            measure_name = MEASURES_ANALYZE{meas_i};
            %##
            cond_plot_store = [];
            group_plot_store = [];
            if ~isempty(group_i)
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i) & TMP_FOOOF_T.cluster_id == groups(group_i);
            else
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
            end
            T_vals_plot = TMP_FOOOF_T(inds,:);
            %-
            y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
            T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
            %## REMOVE NANS
            cond_chars = unique(T_vals_plot.cond_char);
            subjects = unique(T_vals_plot.subj_char);
            t_tmp = [];
            for i = 1:length(subjects)
                ii = find(T_vals_plot.subj_char == subjects(i));
                tt = T_vals_plot(ii,:);
                for j = 1:length(cond_chars)
                    jj = find(tt.cond_char == cond_chars(j));
                    t_tmp = [t_tmp; tt(jj(1),:)];
                end
            end
            inds = isnan(t_tmp.(varnames{var_i}));
            if any(inds)
                t_tmp = t_tmp(~inds,:);
            end
            T_vals_plot = t_tmp;
            %- end remove
            %## LINEAR MODEL
            mod_lme = sprintf('%s ~ 1 + %s + (1|subj_char)',measure_name,varnames{var_i});
            stats_out = fitlme(T_vals_plot,mod_lme); %,'DummyVarCoding','effects');
            anova_out = anova(stats_out);
            %- intercept only model
            % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            R2 = stats_out.Rsquared.Adjusted;
            %## GRAB COEFFICIENTS & PVALUES
            anova_p_var = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            %-
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            %-
            slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i}))));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            %## GATHER STATS
            [norm_h,norm_p] = lillietest(stats_out.residuals);
            [b,bnames,fetable] = stats_out.fixedEffects();
            [br,brnames,bretable] = stats_out.randomEffects();
            STATS_TRACK_STRUCT(cntsts).stat_test_mod = {mod_lme};
            STATS_TRACK_STRUCT(cntsts).resp_terms = {measure_name};
            STATS_TRACK_STRUCT(cntsts).pred_terms = bnames.Name;
            STATS_TRACK_STRUCT(cntsts).rnd_terms = {''};
            %-
            STATS_TRACK_STRUCT(cntsts).anova_kin_p = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).anova_inter_p = anova_out.pValue(strcmp(anova_out.Term,'(Intercept)'));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_kin_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            STATS_TRACK_STRUCT(cntsts).lme_rnd_effects = {bretable};
            STATS_TRACK_STRUCT(cntsts).norm_test_h = norm_h;
            STATS_TRACK_STRUCT(cntsts).norm_test_p = norm_p;
            STATS_TRACK_STRUCT(cntsts).R2 = R2;
            cntsts = cntsts + 1;
            STATS_TRACK_STRUCT(cntsts) = DEF_STATS_TRACK_STRUCT;
            %- end gather stats
            %## SCATTER
            %- plot
            axes();
            hold on;
            data = T_vals_plot;
            [~,inds] = sort(data.(varnames{var_i}));
            data = data(inds,:);
            ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',cond_chars(cond_i)));
            ss.CData = COLORS(1,:);
            ss.SizeData = 15;
            ss.MarkerEdgeAlpha = 0.6;
            %## LINEAR MODEL FIT
            hold on;
            data = T_vals_plot;
            [vals,inds] = sort(data.(varnames{var_i}));
            data = data(inds,:);
            if anova_p_var < ALPHA
                %## PLOT RANDOM EFFECTS
                if DO_PLOT_RANDOM_VARS
                    rnd_colors = linspecer(size(bretable,1));
                    for subj_i = 1:size(bretable,1)
                        x = unique(data.(varnames{var_i}));
                        y = [];
                        for i = 1:length(x)
                            y(i) = x(i)*slope_var + bretable.Estimate(subj_i) + inter_mn ;
                        end
                        pps = plot(x,y,...
                            'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                            'LineWidth',LINE_WIDTH_REFF);
                        % pp.Color = rnd_colors(subj_i)*0.8;
                        pps.Color = [COLORS(2,:),0.35];
                        pps.LineStyle = ':';
                    end
                end
                %## PLOT MAIN EFFECTS
                x = unique(data.(varnames{var_i}));
                y = [];
                for i = 1:length(x)
                    y(i) = x(i)*slope_var + inter_mn ;
                end
                pp = plot(x,y,...
                    'DisplayName',sprintf('pvalue_{%s}=(%0.2f)','kin',anova_p_var),...
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
                stats_store = [stats_store, pp];
            end
            %## AX EDITS
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('')
            end
            xlabel(varnames_labs{var_i});
            title(MEASURE_NAME_LABS{meas_i});
            set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            ylim(y_lim_calc)
            %## legend
            if meas_i == length(MEASURES_ANALYZE)
                %- lg3
                if ~isempty(stats_store)
                    legend(gca,stats_store);
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
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_eeg-kin.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end
%% (SIMPLE MEDIATION MODEL) PREDICTORS: KINEMATICS, RESPONSE: SPEED, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
%-
COLORS = linspecer(2);
IM_RESIZE = 0.7;
TITLE_TXT_SIZE =14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_W = 0.35;
AX_H = 0.25;
AX_TXT_SIZE = 10;
AX_FONT_NAME = 'Arial';
AX_HORIZ_SHIFT = 0.1;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT = 0.65;
%-
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = false;
REG_TXT_SIZE = 8;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.125;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep 'lme_Pkin-Rspeed'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %##
        atlas_name = atlas_name_store{cl_i};
        fig = figure('color','white','renderer','Painters');
        sgtitle(atlas_name,'FontName',AX_FONT_NAME,'FontSize',TITLE_TXT_SIZE,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        vert_shift = 0;
        horiz_shift = 0;
        stats_store = [];
        measure_name = MEASURES_ANALYZE{meas_i};
        %##
        cond_plot_store = [];
        group_plot_store = [];
        inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
        T_vals_plot = TMP_FOOOF_T(inds,:);
        %-
        y_lim_calc = [0,1.25];
        %-
        T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        %## REMOVE NANS
        cond_chars = unique(T_vals_plot.cond_char);
        subjects = unique(T_vals_plot.subj_char);
        t_tmp = [];
        for i = 1:length(subjects)
            ii = find(T_vals_plot.subj_char == subjects(i));
            tt = T_vals_plot(ii,:);
            for j = 1:length(cond_chars)
                jj = find(tt.cond_char == cond_chars(j));
                t_tmp = [t_tmp; tt(jj(1),:)];
            end
        end
        inds = isnan(t_tmp.(varnames{var_i}));
        if any(inds)
            t_tmp = t_tmp(~inds,:);
        end
        T_vals_plot = t_tmp;
        %- end remove
        %## LINEAR MODEL
        mod_lme = sprintf('%s ~ 1 + %s + (1|subj_char)','cond_char',varnames{var_i});
        stats_out = fitlme(T_vals_plot,mod_lme); %,'DummyVarCoding','effects');
        anova_out = anova(stats_out);
        %- intercept only model
        % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
        % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
        % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
        R2 = stats_out.Rsquared.Adjusted;
        %## GRAB COEFFICIENTS & PVALUES
        anova_p_var = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
        %-
        pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
        pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
        %-
        slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i}))));
        inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
        %## GATHER STATS
        [norm_h,norm_p] = lillietest(stats_out.residuals);
        [b,bnames,fetable] = stats_out.fixedEffects();
        [br,brnames,bretable] = stats_out.randomEffects();
        STATS_TRACK_STRUCT(cntsts).stat_test_mod = {mod_lme};
        STATS_TRACK_STRUCT(cntsts).resp_terms = {measure_name};
        STATS_TRACK_STRUCT(cntsts).pred_terms = bnames.Name;
        STATS_TRACK_STRUCT(cntsts).rnd_terms = {''};
        %-
        STATS_TRACK_STRUCT(cntsts).anova_kin_p = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
        STATS_TRACK_STRUCT(cntsts).anova_inter_p = anova_out.pValue(strcmp(anova_out.Term,'(Intercept)'));
        %-
        STATS_TRACK_STRUCT(cntsts).lme_kin_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
        STATS_TRACK_STRUCT(cntsts).lme_inter_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
        %-
        STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
        STATS_TRACK_STRUCT(cntsts).lme_inter_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
        STATS_TRACK_STRUCT(cntsts).lme_rnd_effects = {bretable};
        STATS_TRACK_STRUCT(cntsts).norm_test_h = norm_h;
        STATS_TRACK_STRUCT(cntsts).norm_test_p = norm_p;
        STATS_TRACK_STRUCT(cntsts).R2 = R2;
        cntsts = cntsts + 1;
        STATS_TRACK_STRUCT(cntsts) = DEF_STATS_TRACK_STRUCT;
        %- end gather stats
        %## SCATTER
        %- plot
        axes();
        hold on;
        data = T_vals_plot;
        [~,inds] = sort(data.(varnames{var_i}));
        data = data(inds,:);
        ss = scatter(data,varnames{var_i},'cond_char','DisplayName',sprintf('%s',cond_chars(cond_i)));
        ss.CData = COLORS(1,:);
        ss.SizeData = 15;
        ss.MarkerEdgeAlpha = 0.6;
        %## LINEAR MODEL FIT
        hold on;
        data = T_vals_plot;
        [vals,inds] = sort(data.(varnames{var_i}));
        data = data(inds,:);
        if anova_p_var < ALPHA
            %## PLOT RANDOM EFFECTS
            if DO_PLOT_RANDOM_VARS
                rnd_colors = linspecer(size(bretable,1));
                for subj_i = 1:size(bretable,1)
                    x = unique(data.(varnames{var_i}));
                    y = [];
                    for i = 1:length(x)
                        y(i) = x(i)*slope_var + bretable.Estimate(subj_i) + inter_mn ;
                    end
                    pps = plot(x,y,...
                        'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                        'LineWidth',LINE_WIDTH_REFF);
                    % pp.Color = rnd_colors(subj_i)*0.8;
                    pps.Color = [COLORS(2,:),0.35];
                    pps.LineStyle = ':';
                end
            end
            %## PLOT MAIN EFFECTS
            x = unique(data.(varnames{var_i}));
            y = [];
            for i = 1:length(x)
                y(i) = x(i)*slope_var + inter_mn ;
            end
            pp = plot(x,y,...
                'DisplayName',sprintf('pvalue_{%s}=(%0.2f)','kin',anova_p_var),...
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
            stats_store = [stats_store, pp];
        end
        %## AX EDITS
        if meas_i == 1
            ylabel('Speed (m/s)');
        else
            ylabel('')
        end
        xlabel(varnames_labs{var_i});
        title(MEASURE_NAME_LABS{meas_i});
        set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
        ylim(y_lim_calc)
        %## LEGEND
        if meas_i == length(MEASURES_ANALYZE)
            %- lg3
            if ~isempty(stats_store)
                legend(gca,stats_store);
                [lg3,~,~,~] = legend('boxoff');
                set(lg3,'Orientation','horizontal')
                set(lg3,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE);
                set(lg3,'NumColumns',4);
                set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                    AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.05,lg3.Position(3),lg3.Position(4)]);
                lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
            end
        end
        %## AX POSITION
        set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
        horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
        if mod(meas_i,AX_MAX) == 0
            horiz_shift = 0;
            vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
        end
        hold off;
        %## SAVE
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_speed.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end
%% (SIMPLE MEDIATION MODEL) PREDICTORS: KINEMATIC, SPEED. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
%-
% COLORS = linspecer(4);
COLORS = COLOR_MAPS_SPEED;
IM_RESIZE = 0.7;
TITLE_TXT_SIZE =14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_W = 0.35;
AX_H = 0.25;
AX_TXT_SIZE = 10;
AX_FONT_NAME = 'Arial';
AX_HORIZ_SHIFT = 0.1;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT = 0.65;
%-
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = false;
REG_TXT_SIZE = 8;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.135;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep 'lme_Pkin-Pspeed-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
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
        stats_store = [];
        for meas_i = 1:length(MEASURES_ANALYZE)
            measure_name = MEASURES_ANALYZE{meas_i};
            %##
            cond_plot_store = [];
            group_plot_store = [];
            inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
            T_vals_plot = TMP_FOOOF_T(inds,:);
            %-
            y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
            T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
            %## REMOVE NANS
            cond_chars = unique(T_vals_plot.cond_char);
            subjects = unique(T_vals_plot.subj_char);
            t_tmp = [];
            for i = 1:length(subjects)
                ii = find(T_vals_plot.subj_char == subjects(i));
                tt = T_vals_plot(ii,:);
                for j = 1:length(cond_chars)
                    jj = find(tt.cond_char == cond_chars(j));
                    t_tmp = [t_tmp; tt(jj(1),:)];
                end
            end
            inds = isnan(t_tmp.(varnames{var_i}));
            if any(inds)
                t_tmp = t_tmp(~inds,:);
            end
            T_vals_plot = t_tmp;
            %- end remove
            %## LINEAR MODEL
            mod_lme = sprintf('%s ~ 1 + %s + cond_char + (1|subj_char)',measure_name,varnames{var_i});
            stats_out = fitlme(T_vals_plot,mod_lme); %,'DummyVarCoding','effects');
            anova_out = anova(stats_out);
            %- intercept only model
            % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            R2 = stats_out.Rsquared.Adjusted;
            %## GRAB COEFFICIENTS & PVALUES
            anova_p_var = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            anova_p_cond = anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
            %-
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            pval_cond = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            %-
            slope_cond = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i}))));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            %## GATHER STATS
            [norm_h,norm_p] = lillietest(stats_out.residuals);
            [b,bnames,fetable] = stats_out.fixedEffects();
            [br,brnames,bretable] = stats_out.randomEffects();
            STATS_TRACK_STRUCT(cntsts).stat_test_mod = {mod_lme};
            STATS_TRACK_STRUCT(cntsts).resp_terms = {measure_name};
            STATS_TRACK_STRUCT(cntsts).pred_terms = bnames.Name;
            STATS_TRACK_STRUCT(cntsts).rnd_terms = {''};
            %-
            STATS_TRACK_STRUCT(cntsts).anova_speed_p = anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
            STATS_TRACK_STRUCT(cntsts).anova_kin_p = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).anova_inter_p = anova_out.pValue(strcmp(anova_out.Term,'(Intercept)'));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_speed_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            STATS_TRACK_STRUCT(cntsts).lme_kin_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            STATS_TRACK_STRUCT(cntsts).lme_rnd_effects = {bretable};
            STATS_TRACK_STRUCT(cntsts).norm_test_h = norm_h;
            STATS_TRACK_STRUCT(cntsts).norm_test_p = norm_p;
            STATS_TRACK_STRUCT(cntsts).R2 = R2;
            cntsts = cntsts + 1;
            STATS_TRACK_STRUCT(cntsts) = DEF_STATS_TRACK_STRUCT;
            %- end gather stats
            %## SCATTER
            %- plot
            axes();
            hold on;
            for cond_i = 1:length(cond_chars)
                inds = T_vals_plot.cond_char==double(string(cond_chars(cond_i)));
                data = T_vals_plot(inds,:);
                [vals,inds] = sort(data.(varnames{var_i}));
                data = data(inds,:);
                ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%0.2f',cond_chars(cond_i)));
                ss.CData = COLORS(cond_i,:);
                ss.SizeData = 15;
                ss.MarkerEdgeAlpha = 0.6;
                if meas_i == 1
                    cond_plot_store = [cond_plot_store, ss];
                end
            end
            %## LINEAR MODEL FIT
            hold on;
            for cond_i = 1:length(cond_chars)
                inds = T_vals_plot.cond_char==double(string(cond_chars(cond_i)));
                data = T_vals_plot(inds,:);
                [vals,inds] = sort(data.(varnames{var_i}));
                data = data(inds,:);
                if anova_p_var < ALPHA || anova_p_cond < ALPHA
                    %## PLOT RANDOM EFFECTS
                    if DO_PLOT_RANDOM_VARS
                        rnd_colors = linspecer(size(bretable,1));
                        for subj_i = 1:size(bretable,1)
                            x = unique(data.(varnames{var_i}));
                            xs = double(string(cond_chars(cond_i)));
                            y = [];
                            for i = 1:length(x)
                                y(i) = x(i)*slope_var + xs*slope_cond + inter_mn + bretable.Estimate(subj_i);
                            end
                            pps = plot(x,y,...
                                'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                                'LineWidth',LINE_WIDTH_REFF);
                            % pp.Color = rnd_colors(subj_i)*0.8;
                            pps.Color = [COLORS(cond_i,:)*0.5,0.35];
                            pps.LineStyle = ':';
                        end
                    end
                    %## PLOT MAIN EFFECTS
                    x = unique(data.(varnames{var_i}));
                    xs = double(string(cond_chars(cond_i)));
                    y = [];
                    for i = 1:length(x)
                        y(i) = x(i)*slope_var + xs*slope_cond + inter_mn ;
                    end
                    pp = plot(x,y,...
                        'DisplayName',sprintf('pval_{kin,%s}=%0.2f\npval_{speed,%s}=%0.2f',...
                        MEASURE_SYMBOLS{meas_i},anova_p_var,MEASURE_SYMBOLS{meas_i},anova_p_cond),...
                        'LineWidth',LINE_WIDTH_MEFF);
                    pp.Color = COLORS(cond_i,:)*0.8;
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
                    if cond_i == 1
                        stats_store = [stats_store, pp];
                    end
                end
            end
            %## AX EDITS
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('')
            end
            xlabel(varnames_labs{var_i});
            title(MEASURE_NAME_LABS{meas_i});
            set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            ylim(y_lim_calc)
            %## LEGEND
            if meas_i == 1
                %- lg2
                legend(gca,cond_plot_store);
                [lg2,icons,plots,txt]  = legend('boxoff');
                tmp = get(lg2,'String');
                cnt = 1;
                for i = 1:length(cond_plot_store)
                    tmp{i} = sprintf('%0.2g',double(string(cond_chars(cnt))));
                    cnt = cnt + 1;
                end
                set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                set(lg2,'Orientation','horizontal')
                set(lg2,'Units','normalized')
                set(lg2,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                    AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
                lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
            elseif meas_i == length(MEASURES_ANALYZE)
                %- lg3
                if ~isempty(stats_store)
                    legend(gca,stats_store);
                    [lg3,~,~,~] = legend('boxoff');
                    set(lg3,'Orientation','horizontal')
                    set(lg3,'FontName',AX_FONT_NAME,'FontSize',9);
                    set(lg3,'NumColumns',4);
                    set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                        AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.06,lg3.Position(3),lg3.Position(4)]);
                    lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
            end
            %## AX POSITION
            set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            if mod(meas_i,AX_MAX) == 0
                horiz_shift = 0;
                vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_kin_speed.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end
%% (SIMPLE MEDIATION MODEL) PREDICTORS: KINEMATIC, SPEED, INTERACT. RESPONSE: BRAIN ACTIVITY, STATS TEST
%## PARAMS
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
des_i = 2;
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
%-
COLORS = COLOR_MAPS_SPEED;
IM_RESIZE = 0.7;
TITLE_TXT_SIZE = 14;
ALPHA = 0.05;
DO_PLOT_RANDOM_VARS = false;
%-
AX_W = 0.35;
AX_H = 0.25;
AX_TXT_SIZE = 10;
AX_FONT_NAME = 'Arial';
AX_HORIZ_SHIFT = 0.1;
AX_VERT_SHIFT = 0.1250;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT = 0.65;
%-
LINE_WIDTH_REFF = 1;
LINE_WIDTH_MEFF = 3;
%-
DO_PLOT_R2 = false;
REG_TXT_SIZE = 8;
REG_HORIZ_SHIFT = 0.05;
REG_VERT_SHIFT = 0.05;
%-
LEG_HORIZ_SHIFT = 0;
LEG_VERT_SHIFT =  0.15;
LEG_TXT_SIZE = 9;
LEG_TOKEN_SIZE = 20;
AX_MAX = 3;
tmp_savedir = [save_dir filesep 'lme_Pkin-Pspeed-PIntact-Reeg'];
mkdir(tmp_savedir);
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
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
        stats_store = [];
        for meas_i = 1:length(MEASURES_ANALYZE)
            measure_name = MEASURES_ANALYZE{meas_i};
            %##
            cond_plot_store = [];
            group_plot_store = [];
            if ~isempty(group_i)
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i) & TMP_FOOOF_T.cluster_id == groups(group_i);
            else
                inds = TMP_FOOOF_T.design_id == designs(des_i) & TMP_FOOOF_T.cluster_id == clusters(cl_i);
            end
            T_vals_plot = TMP_FOOOF_T(inds,:);
            %-
            y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
            T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
            %## REMOVE NANS
            cond_chars = unique(T_vals_plot.cond_char);
            subjects = unique(T_vals_plot.subj_char);
            t_tmp = [];
            for i = 1:length(subjects)
                ii = find(T_vals_plot.subj_char == subjects(i));
                tt = T_vals_plot(ii,:);
                for j = 1:length(cond_chars)
                    jj = find(tt.cond_char == cond_chars(j));
                    t_tmp = [t_tmp; tt(jj(1),:)];
                end
            end
            inds = isnan(t_tmp.(varnames{var_i}));
            if any(inds)
                t_tmp = t_tmp(~inds,:);
            end
            T_vals_plot = t_tmp;
            %- end remove
            %## LINEAR MODEL
            mod_lme = sprintf('%s ~ 1 + %s + cond_char + %s:cond_char + (1|subj_char)',measure_name,varnames{var_i},varnames{var_i});
            stats_out = fitlme(T_vals_plot,mod_lme); %,'DummyVarCoding','effects');
            anova_out = anova(stats_out);
            %- intercept only model
            % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            R2 = stats_out.Rsquared.Adjusted;
            %## GRAB COEFFICIENTS & PVALUES
            anova_p_var = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            anova_p_cond = anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
            anova_p_intact = anova_out.pValue(strcmp(anova_out.Term,sprintf('cond_char:%s',varnames{var_i})));
            %-
            pval_inter = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            pval_cond = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            pval_intact = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('cond_char:%s',varnames{var_i})));
            %-
            slope_cond = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i}))));
            slope_intact = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('cond_char:%s',varnames{var_i}))));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            %## GATHER STATS
            [norm_h,norm_p] = lillietest(stats_out.residuals);
            [b,bnames,fetable] = stats_out.fixedEffects();
            [br,brnames,bretable] = stats_out.randomEffects();
            STATS_TRACK_STRUCT(cntsts).stat_test_mod = {mod_lme};
            STATS_TRACK_STRUCT(cntsts).resp_terms = {measure_name};
            STATS_TRACK_STRUCT(cntsts).pred_terms = bnames.Name;
            STATS_TRACK_STRUCT(cntsts).rnd_terms = {''};
            %-
            STATS_TRACK_STRUCT(cntsts).anova_speed_p = anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
            STATS_TRACK_STRUCT(cntsts).anova_kin_p = anova_out.pValue(strcmp(anova_out.Term,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).anova_inter_p = anova_out.pValue(strcmp(anova_out.Term,'(Intercept)'));
            STATS_TRACK_STRUCT(cntsts).anova_intact_p = anova_out.pValue(strcmp(anova_out.Term,sprintf('cond_char:%s',varnames{var_i})));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_speed_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            STATS_TRACK_STRUCT(cntsts).lme_kin_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            STATS_TRACK_STRUCT(cntsts).lme_intact_p = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,sprintf('cond_char:%s',varnames{var_i})));
            %-
            STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            STATS_TRACK_STRUCT(cntsts).lme_kin_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_inter_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'));
            STATS_TRACK_STRUCT(cntsts).lme_intact_coeff = stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,sprintf('cond_char:%s',varnames{var_i})));
            STATS_TRACK_STRUCT(cntsts).lme_rnd_effects = {bretable};
            STATS_TRACK_STRUCT(cntsts).norm_test_h = norm_h;
            STATS_TRACK_STRUCT(cntsts).norm_test_p = norm_p;
            STATS_TRACK_STRUCT(cntsts).R2 = R2;
            cntsts = cntsts + 1;
            STATS_TRACK_STRUCT(cntsts) = DEF_STATS_TRACK_STRUCT;
            %- end gather stats
            %## SCATTER
            %- plot
            axes();
            hold on;
            for cond_i = 1:length(cond_chars)
                inds = T_vals_plot.cond_char==double(string(cond_chars(cond_i)));
                data = T_vals_plot(inds,:);
                [vals,inds] = sort(data.(varnames{var_i}));
                data = data(inds,:);
                ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%0.2f',cond_chars(cond_i)));
                ss.CData = COLORS(cond_i,:);
                ss.SizeData = 15;
                % ss.MarkerFaceAlpha = 'flat';
                % ss.AlphaData = repmat(0.3,[size(data,1),1]);
                ss.MarkerEdgeAlpha = 0.6;
                if meas_i == 1
                    cond_plot_store = [cond_plot_store, ss];
                end
            end
            %## LINEAR MODEL FIT
            hold on;
            for cond_i = 1:length(cond_chars)
                inds = T_vals_plot.cond_char==double(string(cond_chars(cond_i)));
                data = T_vals_plot(inds,:);
                [vals,inds] = sort(data.(varnames{var_i}));
                data = data(inds,:);
                if anova_p_intact < ALPHA
                    %## PLOT RANDOM EFFECTS
                    if DO_PLOT_RANDOM_VARS
                        rnd_colors = linspecer(size(bretable,1));
                        for subj_i = 1:size(bretable,1)
                            x = unique(data.(varnames{var_i}));
                            xs = double(string(cond_chars(cond_i)));
                            y = [];
                            for i = 1:length(x)
                                y(i) = x(i)*slope_var + xs*slope_cond + xs*x(i)*slope_intact + bretable.Estimate(subj_i) + inter_mn ;
                            end
                            pps = plot(x,y,...
                                'DisplayName',sprintf('subj. %s',bretable.Level{subj_i}),...
                                'LineWidth',LINE_WIDTH_REFF);
                            % pp.Color = rnd_colors(subj_i)*0.8;
                            pps.Color = [COLORS(cond_i,:)*0.5,0.35];
                            pps.LineStyle = ':';
                        end
                    end
                    %## PLOT MAIN EFFECTS
                    x = unique(data.(varnames{var_i}));
                    xs = double(string(cond_chars(cond_i)));
                    y = [];
                    for i = 1:length(x)
                        y(i) = x(i)*slope_var + xs*slope_cond + xs*x(i)*slope_intact + inter_mn ;
                    end
                    pp = plot(x,y,...
                        'DisplayName',sprintf('pval_{kin,%s}=%0.2f\npval_{speed,%s}=%0.2f\npval_{intact,%s}=%0.2f',...
                        MEASURE_SYMBOLS{meas_i},anova_p_var,MEASURE_SYMBOLS{meas_i},anova_p_cond,MEASURE_SYMBOLS{meas_i},anova_p_intact),...
                        'LineWidth',LINE_WIDTH_MEFF);
                    pp.Color = COLORS(cond_i,:)*0.8;
                    
                    %## GIVE CORRELATION COEFFICIENT
                    if DO_PLOT_R2
                        eq = '';
                        if anova_p_intact > 0.01 && anova_p_intact < ALPHA
                            annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                                'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                'BackgroundColor','none');
                        elseif anova_p_intact <= 0.01 && anova_p_intact > 0.001
                            annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                                'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                'BackgroundColor','none');
                        elseif anova_p_intact <= 0.001
                            annotation('textbox',[horiz_shift+REG_HORIZ_SHIFT*IM_RESIZE,AX_INIT_VERT+vert_shift+REG_VERT_SHIFT*IM_RESIZE,0.2,0.2],...
                                'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                'VerticalAlignment','middle','LineStyle','none','FontName',AX_FONT_NAME,...
                                'FontSize',REG_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                'BackgroundColor','none');
                        end
                    end
                    if cond_i == 1
                        stats_store = [stats_store, pp];
                    end
                end
            end
            %## AX EDITS
            if meas_i == 1
                ylabel('10*log_{10}(Flattened PSD)');
            else
                ylabel('')
            end
            xlabel(varnames_labs{var_i});
            title(MEASURE_NAME_LABS{meas_i});
            set(gca,'FontWeight','bold','FontSize',AX_TXT_SIZE);
            ylim(y_lim_calc)
            %## LEGEND
            if meas_i == 1
                %- lg2
                legend(gca,cond_plot_store);
                [lg2,icons,plots,txt]  = legend('boxoff');
                tmp = get(lg2,'String');
                cnt = 1;
                for i = 1:length(cond_plot_store)
                    tmp{i} = sprintf('%0.2g',double(string(cond_chars(cnt))));
                    cnt = cnt + 1;
                end
                set(lg2,'String',tmp,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE)
                set(lg2,'Orientation','horizontal')
                set(lg2,'Units','normalized')
                set(lg2,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                    AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE,lg2.Position(3),lg2.Position(4)]);
                lg2.ItemTokenSize(1) = LEG_TOKEN_SIZE;
            elseif meas_i == length(MEASURES_ANALYZE)
                %- lg3
                if ~isempty(stats_store)
                    legend(gca,stats_store);
                    [lg3,~,~,~] = legend('boxoff');
                    set(lg3,'Orientation','horizontal')
                    set(lg3,'FontName',AX_FONT_NAME,'FontSize',LEG_TXT_SIZE);
                    set(lg3,'NumColumns',4);
                    % set(lg1,'Position',[0.1,AX_INIT_VERT+AX_W*im_resize-vert_shift,lg1.Position(3),lg1.Position(4)]);
                    set(lg3,'Position',[AX_INIT_HORIZ+LEG_HORIZ_SHIFT*IM_RESIZE,...
                        AX_INIT_VERT+AX_H*IM_RESIZE+LEG_VERT_SHIFT*IM_RESIZE-0.075,lg3.Position(3),lg3.Position(4)]);
                    lg3.ItemTokenSize(1) = LEG_TOKEN_SIZE;
                end
            end
            %## AX SHIFT
            set(gca,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            if mod(meas_i,AX_MAX) == 0
                horiz_shift = 0;
                vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
            end
        end
        hold off;
        %##
        exportgraphics(fig,[tmp_savedir filesep sprintf('cl%s_kin%s_plot.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end
%% WRITE STATISTICS TO MAT & XLSX
par_save(STATS_TRACK_STRUCT,save_dir,'stats_track_struct.mat');
%-
tt = struct2table(STATS_TRACK_STRUCT);
tt = removevars(tt,{'lme_rnd_effects'});
tmp = cellfun(@(x) strjoin(x,','),tt.pred_terms,'UniformOutput',false);
tt.pred_terms = tmp;
writetable(tt,[save_dir filesep 'stats_track_struct.xlsx']);