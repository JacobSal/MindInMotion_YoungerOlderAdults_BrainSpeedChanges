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
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
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
%## Add Study & Script Paths
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
study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
%- study info
SUB_GROUP_FNAME = 'group_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 12;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## LOAD STUDY
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
    plot_store_dir = [cluster_dir filesep 'plots_out' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    plot_store_dir = [cluster_dir filesep 'plots_out'];
end
if ~exist(spec_data_dir,'dir')
    error('spec_data dir does not exist');
end
if ~exist(plot_store_dir,'dir')
    mkdir(plot_store_dir);
end
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',spec_data_dir);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',spec_data_dir);
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {STUDY.datasetinfo.subject};
%-
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
CLUSTER_PICKS = main_cl_inds(2:end);
DESIGN_I = 1:length(STUDY.design);
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%%
%## CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
%     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
% allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
%## (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
ERSP_TIMERANGE=[averaged_warpto_events(1), averaged_warpto_events(end)];
STUDY.etc.averaged_warpto_events = averaged_warpto_events;
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
disp(averaged_warpto_events);
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
[STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
txt_store = cell(length(clusters),1);
atlas_name_store = cell(length(clusters),1);
for k_i = 1:length(clusters)
    k = double(string(clusters(k_i)));
    %## ANATOMY
    dip1 = STUDY.cluster(k).all_diplocs;
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
PLOT_PARAMS = struct('color_map',linspecer(4),...
                'cond_labels',unique(FOOOF_TABLE.cond_char),'group_labels',unique(FOOOF_TABLE.group_char),...
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label','10*log_{10}(Flattened PSD)',...
                'title','','font_size',9,'y_lim',[-1,15],...
                'font_name','Arial','x_label','');
%% ================================================================== %%
%{
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
                    FOOOF_TABLE.(varnames{var_i})(i) = 0;
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
                    FOOOF_TABLE.(varnames{var_i})(i) = 0;
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
                    FOOOF_TABLE.(varnames{var_i})(i) = 0;
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
                    FOOOF_TABLE.(varnames{var_i})(i) = 0;
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
%}
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
varnames = {'mean_APexc_COV','mean_APexc_mean','mean_MLexc_COV','mean_MLexc_mean','mean_StepDur','mean_UDexc_COV','mean_UDexc_mean'};
varnames_labs = {'APexc COV','APexc','MLexc COV','MLexc','StepDur','UDexc COV','UDexc'};
%% PREDICTORS: SPEED CONDITION, RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 1.1;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
DO_PLOT_GROUPS = true;
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
       T_vals_plot = t_tmp;
       % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
       try
           mod = sprintf('%s ~ 1 + %s',varnames{var_i},'cond_char');
           % stats_out = fitlme(T_vals_plot,mod);
           stats_out = fitlm(T_vals_plot,mod);
           anova_out = anova(stats_out);
           R2 = stats_out.Rsquared.Adjusted;
           %## FPRINTF
            fid = fopen([kin_savedir filesep sprintf('%s_kinematics-speed_ANOVA.txt',varnames{var_i})],'wt');
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
           % anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
           % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
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
        STATS_STRUCT = struct('anova',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{{}},...
                              'regress_xvals',0);
        if DO_PLOT_GROUPS
            for gg = 1:length(groups)
                STATS_STRUCT.anova{gg}=anova_p_var;
                STATS_STRUCT.regress_pval{gg}=pval_var;
                STATS_STRUCT.regress_line{gg}=[inter_mn,slope_var];
            end
            STATS_STRUCT.r2_coeff=repmat(R2,length(groups),1);
            STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        else
            STATS_STRUCT.anova={anova_p_var};
            STATS_STRUCT.regress_pval={pval_var};
            STATS_STRUCT.regress_line={[inter_mn,slope_var]};
            STATS_STRUCT.r2_coeff=R2;
            STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        end
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        STATS_STRUCT = struct('anova',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{{}},...
                              'regress_xvals',0);
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_PARAMS = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},...
            'title',varnames_labs{var_i},'font_size',10,'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',~DO_PLOT_GROUPS);
        % ax = axes();
        fig = figure('color','white','renderer','Painters');
        % sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_id',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_PARAMS',PLOT_PARAMS,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[kin_savedir filesep sprintf('%s_kinematics-speed_grouped.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% PREDICTORS: SPEED CONDITION, GROUP, & INTERACTION; RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 0.9;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
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
        T_vals_plot = t_tmp;
        % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        try
            mod = sprintf('%s ~ 1 + %s + group_id + %s*group_id',varnames{var_i},'cond_char','cond_char');
            % stats_out = fitlme(T_vals_plot,mod);
            stats_out = fitlm(T_vals_plot,mod);
            anova_out = anova(stats_out);
            R2 = stats_out.Rsquared.Adjusted;
            %## FPRINTF
            fid = fopen([kin_savedir filesep sprintf('%s_kinematics-speed-group-interact_ANOVA.txt',varnames{var_i})],'wt');
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
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Name,'cond_char'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            % inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
            figure; plot(stats_out); plot(stats_out.Formula)
            pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id')));
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_id')));

            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
        catch e
            fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
            R2 = 0;
            pval = 1;
            slope = 0;
            inter = 0;
        end
        STATS_STRUCT = struct('anova',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{{}},...
                              'regress_xvals',0);
        if DO_PLOT_GROUPS
            % STATS_STRUCT.anova={repmat(anova_p_var,length(groups),1)};
            % STATS_STRUCT.regress_pval={repmat(pval_var,length(groups),1)};
            % STATS_STRUCT.regress_line={repmat([inter_mn,slope_var],length(groups),1)};
            % STATS_STRUCT.r2_coeff={repmat(R2,length(groups),1)};
            % STATS_STRUCT.regress_xvals=[0,unique(T_vals_plot.cond_char)',1.25];
            for gg = 1:length(groups)
                STATS_STRUCT.anova{gg}=anova_p_var;
                STATS_STRUCT.regress_pval{gg}=pval_var;
                STATS_STRUCT.regress_line{gg}=[inter_mn,slope_var];
            end
            STATS_STRUCT.r2_coeff=repmat(R2,length(groups),1);
            STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        else
            STATS_STRUCT.anova={anova_p_var};
            STATS_STRUCT.regress_pval={pval_var};
            STATS_STRUCT.regress_line={[inter_mn,slope_var]};
            STATS_STRUCT.r2_coeff=R2;
            STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        end
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        STATS_STRUCT = struct('anova',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{{}},...
                              'regress_xvals',0);
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_PARAMS = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},...
            'title',varnames_labs{var_i},'font_size',10,'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',true);
        % ax = axes();
        fig = figure('color','white','renderer','Painters');
        % sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_id',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_PARAMS',PLOT_PARAMS,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[kin_savedir filesep sprintf('%s_kinematics-speed-group-interact.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% PREDICTORS: SPEED CONDITION, GROUP; RESPONSE: KINEMATICS, STATS TEST
STATS_OUT = [];
im_resize= 0.9;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
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
        T_vals_plot = t_tmp;
        % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
        try
            mod = sprintf('%s ~ 1 + %s + group_id',varnames{var_i},'cond_char');
            % stats_out = fitlme(T_vals_plot,mod);
            stats_out = fitlm(T_vals_plot,mod);
            anova_out = anova(stats_out);
            R2 = stats_out.Rsquared.Adjusted;
            %## FPRINTF
            fid = fopen([kin_savedir filesep sprintf('%s_kinematics-speed-group_ANOVA.txt',varnames{var_i})],'wt');
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
            % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char'));
            % tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Name,'cond_char'));
            % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')));
            % inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
            anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,'cond_char'));
            pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char'));
            slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
            inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'(Intercept)')));
        catch e
            fprintf('Error. Cluster %s\n\n%s\n',string(clusters(cl_i)),getReport(e))
            R2 = 0;
            pval = 1;
            slope = 0;
            inter = 0;
        end
        STATS_STRUCT = struct('anova',{{}},...
                              'pvals',{{}},...
                              'pvals_pairs',{{}},...
                              'regress_pval',{{}},...
                              'regress_line',{{}},...
                              'r2_coeff',{{}},...
                              'regress_xvals',0);
        STATS_STRUCT.anova={anova_p_var};
        STATS_STRUCT.regress_pval={pval_var};
        STATS_STRUCT.regress_line={[inter_mn,slope_var]};
        STATS_STRUCT.r2_coeff=R2;
        STATS_STRUCT.regress_xvals=[0,unique(double(string(T_vals_plot.cond_char)))',1.25];
        STATS_OUT = [STATS_OUT; STATS_STRUCT];
        % figure;
        VIOLIN_PARAMS = {'width',0.1,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
        PLOT_PARAMS = struct('color_map',color_dark,...
            'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
            'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},...
            'title',varnames_labs{var_i},'font_size',10,'ylim',[min(T_vals_plot.(varnames{var_i}))-std(T_vals_plot.(varnames{var_i})),max(T_vals_plot.(varnames{var_i}))+std(T_vals_plot.(varnames{var_i}))],...
            'font_name','Arial','x_label','speed','do_combine_groups',true);
        % ax = axes();
        fig = figure('color','white','renderer','Painters');
        % sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        set(gca,AXES_DEFAULT_PROPS{:})
        axax = group_violin(T_vals_plot,varnames{var_i},'cond_char','group_id',...
            fig,...
            'VIOLIN_PARAMS',VIOLIN_PARAMS,...
            'PLOT_PARAMS',PLOT_PARAMS,...
            'STATS_STRUCT',STATS_STRUCT);
        set(axax,'OuterPosition',[0,0,1,1]);
        set(axax,'Position',[0.1+horiz_shift,VIOLIN_BOTTOM+vert_shift,AX_W*im_resize,AX_H*im_resize]);  %[left bottom width height]
        hold off;
        exportgraphics(fig,[kin_savedir filesep sprintf('%s_kinematics-speed-group-interact.tiff',varnames{var_i})],'Resolution',300)
        close(fig)
        %- iterate
   end
end
%% ===================================================================== %%
%{
pval_out = zeros(length(STATS_OUT),1);
lens = zeros(length(STATS_OUT),1);
for i = 1:length(STATS_OUT)
    pval_out(i) = STATS_OUT(i).regress_pval{1};
    lens(i) = length(varnames{i});
end
[val,ind] = sort(pval_out);
varnames = varnames(ind);
format = sprintf('\n\n   %%s%%%is\t\t%%s\n',max(lens)-length('Name'));
fprintf(format,'Name','','pValue');
for i = 1:length(varnames)
    format = sprintf('%i) %%s%%%is\t\t%%0.2g\n',i,max(lens)-length(varnames{i}));
    fprintf(format,varnames{i},'',val(i));
end

varnames = varnames(1:3);
varnames = unique([varnames, vv_main_anl]);
inds = cellfun(@(x) find(strcmp(x,varnames_hold)),varnames);
varnames_labs = varnames_labs(inds);
%}
%## HARDCODED VARIABLES OF INTEREST
varnames = {'mean_APexc_COV','mean_APexc_mean','mean_MLexc_COV','mean_MLexc_mean','mean_StepDur','mean_UDexc_COV','mean_UDexc_mean'};
varnames_labs = {'APexc COV','APexc','MLexc COV','MLexc','StepDur','UDexc COV','UDexc'};
%% PREDICTORS: KINEMATICS. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
REGRESS_TXT_XMULTI = 0.9;
REGRESS_TXT_YMULTI = 0.9;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25;
GROUP_SHORTS = {'HO','FO'};
MEASURE_NAME_LABS = {'Mean \theta','Mean \alpha','Mean \beta'};
for cl_i = 1:length(clusters)
    %##
    for var_i = 1:length(varnames)
        %##
        %%
        PLOT_PARAMS = struct('color_map',color_dark,...
                'cond_labels',unique(T_vals_plot.cond_char),'group_labels',unique(T_vals_plot.group_char),...
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label',varnames_labs{var_i},...
                'title',varnames_labs{var_i},'font_size',10,'ylim',[],...
                'font_name','Arial','x_label','speed','do_combine_groups',true);
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
                T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                x_txt = min(T_vals_plot.(varnames{var_i}))*REGRESS_TXT_XMULTI+std(T_vals_plot.(varnames{var_i}));
                y_txt = max(T_vals_plot.(measure_name))*REGRESS_TXT_YMULTI+std(T_vals_plot.(measure_name));
                try
                    mod = sprintf('%s ~ 1 + %s',measure_name,varnames{var_i});
                    stats_out = fitlm(T_vals_plot,mod);
                    anova_out = anova(stats_out);
                    R2 = stats_out.Rsquared.Adjusted;
                    % pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,varnames{var_i}));
                    % slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,varnames{var_i})));
                    % inter_mn = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)')));
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i})));
                    % slope_grp = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,'group_char_H3000''s')));
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
                        % ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',group_chars(group_i)));
                        ss = scatter(data,varnames{var_i},measure_name,'DisplayName',sprintf('%s',GROUP_SHORTS{group_i}));
                        if group_i == 1
                            ss.Marker = 'o';
                            ss.CData = color_light(cond_i,:);
                            ss.SizeData = 15;
                            cond_plot_store = [cond_plot_store, ss];
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
                        if pval_var < 0.1
                            x = unique(data.(varnames{var_i}));
                            y = [];
                            for i = 1:length(x)
                                % y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*x(i) + inter_mn;
                                y(i) = x(i)*slope_var + inter_mn;
                            end
                            pp = plot(x,y,...
                                'DisplayName',sprintf('p-val_{%s}=(%0.1g)',GROUP_SHORTS{group_i},pval_var),...
                                'LineWidth',2);
                            % pp = plot(x,y,...
                            %     'DisplayName',sprintf('%s',loc_cond_chars(cond_i)),...
                            %     'LineWidth',2);
                            if group_i == 1
                                pp.LineStyle = '-';
                                pp.Color = color_dark(cond_i,:);
                                % cond_plot_store = [cond_plot_store, pp];
                            else
                                pp.LineStyle = '-.';
                                pp.Color = color_dark(cond_i,:)+GROUP_CMAP_OFFSET;
                            end
                            if cond_i == 1
                                % eq = sprintf('y=(%0.1g)*x+(%0.1g)*%0.2f+(%0.1g)*x+(%0.1g)',slope_var,slope_cnd,loc_cond_chars(cond_i),slope_grp,inter_mn);
                                eq = sprintf('y=(%0.1g)*x+(%0.1g)',slope_var,inter_mn);
                                if pval_var > 0.01 & pval_var < 0.05
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
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
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[kin_savedir filesep sprintf('cl%s_%s_kin-eeg.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        % close(fig)
    end
end
%% PREDICTORS: KINEMATICS, COND, & INTERACTION. RESPONSE: BRAIN ACTIVITY, STATS TEST
REGRESS_TXT_SIZE = 8;
REGRESS_TXT_XMULTI = 0.9;
REGRESS_TXT_YMULTI = 1.0;
im_resize = 0.7;
AX_W = 0.35;
AX_H = 0.25;
GROUP_SHORTS = {'HO','FO'};
PLOT_PARAMS = struct('color_map',color_dark,...
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
                %- continuous or categorical
                % T_vals_plot.cond_char = double(string(T_vals_plot.cond_char));
                %-
                loc_cond_chars = unique(T_vals_plot.cond_char);
                y_lim_calc = [min(T_vals_plot.(measure_name))-std(T_vals_plot.(measure_name)),max(T_vals_plot.(measure_name))+std(T_vals_plot.(measure_name))];
                x_txt = min(T_vals_plot.(varnames{var_i}))*REGRESS_TXT_XMULTI+std(T_vals_plot.(varnames{var_i}));
                y_txt = max(T_vals_plot.(measure_name))*REGRESS_TXT_YMULTI+std(T_vals_plot.(measure_name));
                try
                    mod = sprintf('%s ~ 1 + %s + cond_char + %s*cond_char',measure_name,varnames{var_i},varnames{var_i});
                    stats_out = fitlm(T_vals_plot,mod);
                    anova_out = anova(stats_out);
                    R2 = stats_out.Rsquared.Adjusted;
                    anova_p_var = anova_out.pValue(strcmp(anova_out.Properties.RowNames,varnames{var_i}));
                    pval_var = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    pval_cnd = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Properties.RowNames,'cond_char')));
                    tstat_var = stats_out.Coefficients.tStat(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i}));
                    slope_var = double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Properties.RowNames,varnames{var_i})));
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
                        if pval_var < 0.1 || pval_cnd < 0.1
                            x = unique(data.(varnames{var_i}));
                            y = [];
                            for i = 1:length(x)
                                % y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + slope_grp*x(i) + inter_mn;
                                y(i) = x(i)*slope_var + loc_cond_chars(cond_i)*slope_cnd + inter_mn;
                            end
                            pp = plot(x,y,...
                                'DisplayName',sprintf('p-vals=(%0.2f,%0.2f)',GROUP_SHORTS{group_i},pval_var,pval_cnd),...
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
                                eq = sprintf('y=(%0.1g)*x+(%0.1g)*c_i\n%6s+(%0.1g)',slope_var,slope_cnd,'',inter_mn);
                                if pval_var > 0.01 & pval_var < 0.05
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('* %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
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
            end
            %## TITLE
            % annotation('textbox',[0.5-0.1,SCATTER_BOTTOM-vert_shift-0.05+AX_H*im_resize,0.2,0.2],...
            %     'String',string(design_chars(des_i)),'HorizontalAlignment','center',...
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[kin_savedir filesep sprintf('cl%s_%s_cond_kin-eeg.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
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
PLOT_PARAMS = struct('color_map',color_dark,...
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
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
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
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
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
PLOT_PARAMS = struct('color_map',color_dark,...
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
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                elseif pval_var <= 0.01 & pval_var > 0.001
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
                                        'FontSize',REGRESS_TXT_SIZE,'FontWeight','Bold','Units','normalized',...
                                        'BackgroundColor','none');
                                else
                                    annotation('textbox',[horiz_shift+0.1,SCATTER_BOTTOM-vert_shift+0.05,0.2,0.2],...
                                        'String',sprintf('*** %s\nR^2=%0.2g',eq,R2),'HorizontalAlignment','center',...
                                        'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
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
            %     'VerticalAlignment','middle','LineStyle','none','FontName',PLOT_PARAMS.font_name,...
            %     'FontSize',14,'FontWeight','Bold','Units','normalized');
            vert_shift = vert_shift + AX_H*im_resize+0.1;
        end
        hold off;
        %##
        exportgraphics(fig,[kin_savedir filesep sprintf('cl%s_%s_group_cond_kin-eeg.tiff',string(clusters(cl_i)),varnames{var_i})],'Resolution',300)
        close(fig)
    end
end

%% ===================================================================== %%
%## GROUP PLOTS
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
                    otherwise
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
    close(fig);
end
%% ================================================================= %%
im_resize= 0.9;
VIOLIN_BOTTOM = 0.7;
AX_H  = 0.2;
AX_W = 0.25;
%## VIOLIN PLOTS
for k_i = 1:length(clusters)
    atlas_name = atlas_name_store{k_i};
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
    STATS_STRUCT = struct('anova',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',0,...
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
    for des_i = 1:length(designs)   
        inds = FOOOF_TABLE.design_id == num2str(des_i) & FOOOF_TABLE.cluster_id == num2str(cl_i);
        T_plot = FOOOF_TABLE(inds,:);
        switch des_i
            case 1
                
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
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
                'cond_offsets',[-0.3,-0.1,0.1,0.3],'y_label','10*log_{10}(Flattened PSD)',...
                'title',title_plot{i},'font_size',9,'ylim',[min(T_plot.(measure_name))-0.3,max_vals(i)],...
                'font_name','Arial','x_label',x_label,'do_combine_groups',false);
            % ax = axes();
            % figfig = figure();
            axax = group_violin(T_plot,measure_name,'cond_id','group_id',...
                fig,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_PARAMS',PLOT_PARAMS,...
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