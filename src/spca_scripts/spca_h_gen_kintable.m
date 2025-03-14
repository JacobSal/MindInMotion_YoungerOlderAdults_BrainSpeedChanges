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
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
%- study info
SUB_GROUP_FNAME = ['group_spec' filesep 'split_band_test'];
% SUB_GROUP_FNAME = 'all_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt'];
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt_rb5'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## SET STUDY PATHS
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    save_dir = [cluster_dir filesep 'psd_calcs' filesep SUB_GROUP_FNAME];
else
    save_dir = [cluster_dir filesep 'psd_calcs'];
end
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ================================================================== %%
%## CREATE EEG-KINEMATICS TABLE
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