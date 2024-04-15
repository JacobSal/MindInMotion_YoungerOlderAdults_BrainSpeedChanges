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
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = true;
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
DO_SUBJ_PLOTS = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
%- ERSP STAT PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
%- SPECTRUM PARAMS
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
%- ERSP PARAMS
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
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
%- pull out pertinent study info
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
% subject_chars = {STUDY.datasetinfo.subject};
% fPaths = {STUDY.datasetinfo.filepath};
% fNames = {STUDY.datasetinfo.filename};
CLUSTER_PICKS = main_cl_inds(2:end);
DESIGN_I = 1:length(STUDY.design);
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## RE-POP PARAMS
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
    'groupstats',ERSP_STAT_PARAMS.groupstats,...
    'method',ERSP_STAT_PARAMS.method,...
    'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
    'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
    'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
    'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_TIMERANGE);
%% STUDY STATS
for subj_i=1:length(ALLEEG)
    EEG = ALLEEG(subj_i);
    CAT = EEG.etc.COND_CAT(cond_i);
    fprintf('%s) Number of trials: %i',EEG.subject,CAT.trials);
end
%% FOOOF SETUP & PYTHON
%{
% Check which python is being used
pyversion

% The print out from above should tell you which Python you are calling
%  It should show that you are using Python version 3.X
%  If you are using anaconda, it should show your Python is in the anaconda folder
%  If either of these things are not right, reset which Python you are using, as below

% Set python version to use
%  Note: you must do this first thing after opening Matlab (relaunch if you need to)
%  You should only ever have to run this at most, once.
%  You might need to change the path to where your python or anaconda install is
%    For example, your anaconda folder might be `anaconda3` instead of `anaconda`
%    or your anaconda path might be somewhere else, for example, '/opt/anaconda3/bin/python'
%## MANUALLY SET TO A VERSION LOWER THAT 3.10
pyversion('C:\Users\jsalminen\AppData\Local\Programs\Python\Python37\python.EXE');
% pyversion('C:\Users\jsalminen\AppData\Local\Microsoft\WindowsApps\python37.EXE');
% pyversion('/cygdrive/c/Users/jsalminen/AppData/Local/Microsoft/WindowsApps/python3.8')
pe = pyenv;
%}
%% Setup Fooof
SAVE_DATA = true;
settings = struct();  % Use defaults
settings.peak_width_limits = [1 8];%default [1 6] % Amanda used [1 8] in her paper
settings.min_peak_height = 0.05;
% settings.peak_threshold = 2;
settings.max_n_peaks = 3; % originally set to be 2 - 2023-06-07
% the settings are consitent with fooof on github
f_range = [3, 40];
theta_band = [4, 8];
alpha_band = [8 12];
beta_band  = [12 30];
cond_terrains = {'flat','low','med','high'};
cond_speeds = {'0p25','0p5','0p75','1p0'};
COND_CHARS = {cond_terrains,cond_speeds};
% spec_data_cl2_0p250p50p751p0_subbaselined_commonbase
% spec_data_cl2_flatlowmedhigh_subtractmean
keywords = {'Terrain','Speed'};
% save_path = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\icrej_6\14';
% data_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\icrej_6\14\spec_data';
% save_path = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\icrej_5\14';
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% mkdir(fullfile(save_path, 'Fooof_Plots'))
%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
speed_alleeg = cell(length(ALLEEG),2);
for i = 1:length(ALLEEG)
    ss = ALLEEG(i).subject;
    ind = speed_table.Var1==ss;
    chk1 = strcmp(ALLEEG(i).group,SUB_GROUP_FNAME_REGEX) || isempty(SUB_GROUP_FNAME_REGEX);
    if any(ind) && chk1
        speed_alleeg{i,1} = speed_table.Var1(ind);
        speed_alleeg{i,2} = double(speed_table.Var2(ind));
    end
end
speed_alleeg = speed_alleeg(~cellfun(@isempty,speed_alleeg(:,1)),:);
%% (TABLE) GENERATE FOOOF VALUES ======================================= %%
fooof_group_results_org = cell(1,length(DESIGN_I));
fooof_results = cell(length(DESIGN_I),1);
fooof_diff_store = cell(length(DESIGN_I),1);
fooof_apfit_store = cell(length(DESIGN_I),1);
spec_data_original = cell(length(DESIGN_I),1);
for g = DESIGN_I 
    for k = CLUSTER_PICKS
        if isempty(SUB_GROUP_FNAME_REGEX)
            file_mat = [spec_data_dir filesep sprintf('spec_data_cl%i_%s_subtractmean.mat',k,[COND_CHARS{g}{:}])]; %['readSPEC_',num2str(g),'_',keywords{g},'_subSubjectMean.mat'];
        else
            file_mat = [spec_data_dir filesep sprintf('spec_data_cl%i_%s_%s_subtractmean.mat',k,[COND_CHARS{g}{:}]),SUB_GROUP_FNAME_REGEX];
        end
%         file_mat = [data_fpath filesep sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',k,[COND_CHARS{g}{:}])]; %['readSPEC_',num2str(g),'_',keywords{g},'_subSubjectMean.mat'];
        tmp = par_load(file_mat,[]);
        % Note: spec_subj_mean_stats separate young and older adults
        specdata = {tmp.specdata}';
        specfreqs = {tmp.specfreqs}';
        specfreqs = specfreqs{1};
        %% Run fooof
        % Input should be in linear spacing   
        i_ind = 0;
        %- get subjects in cluster
        s_chars = {STUDY.datasetinfo(STUDY.cluster(k).sets).subject};
        if ~isempty(SUB_GROUP_FNAME_REGEX)
            g_inds = cellfun(@(x) strcmp(x,SUB_GROUP_FNAME_REGEX),{STUDY.datasetinfo(STUDY.cluster(k).sets).group});
        else
            g_inds = 1:length(s_chars);
        end
        cl_chars = s_chars(g_inds);
        cl_inds = find(g_inds);
        cl_comps = cluster_update(k).comps(g_inds);
        cl_speeds = zeros(length(cl_chars),1);
        for i = 1:length(cl_speeds)
            ind = cellfun(@(x) x == categorical(cl_chars(i)),speed_alleeg(:,1));
            cl_speeds(i) = speed_alleeg{ind,2};
        end
        for group = 1:size(specdata,2) % in case there is young and old adult group
            for cond = 1:size(specdata,1) % different level of terrains
                specdata_nolog = 10.^(specdata{cond,group}/10);
                % Run FOOOF
                return_model = true;
                for i = 1:size(specdata{cond,group},2)
                    fooof_results{g}{cond,group}{i} = fooof(specfreqs, specdata_nolog(:,i), f_range, settings, return_model);
%                     fooof_group_results_org{g}{k}(i_ind + i).subjects = s_chars(g_inds);
                    fooof_group_results_org{g}{k}(i_ind + i).speed_ms = cl_speeds(i);
                    fooof_group_results_org{g}{k}(i_ind + i).subID = cl_inds(i); %cluster_update(k).sets(i);
                    fooof_group_results_org{g}{k}(i_ind + i).sub_char = categorical(cl_chars(i)); %cluster_update(k).sets(i);
                    fooof_group_results_org{g}{k}(i_ind + i).compID = cl_comps(i); %cluster_update(k).comps(i);
                    %-
                    fooof_group_results_org{g}{k}(i_ind + i).study = g;%1 = terrain, 2 = speed
                    fooof_group_results_org{g}{k}(i_ind + i).cond = cond;
                    fooof_group_results_org{g}{k}(i_ind + i).group = group;
                    fooof_group_results_org{g}{k}(i_ind + i).cluster = k;
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_exp = fooof_results{g}{cond,group}{i}.aperiodic_params(2);
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_offset = fooof_results{g}{cond,group}{i}.aperiodic_params(1);
                    fooof_group_results_org{g}{k}(i_ind + i).central_freq = fooof_results{g}{cond,group}{i}.peak_params(:,1);
                    fooof_group_results_org{g}{k}(i_ind + i).power = fooof_results{g}{cond,group}{i}.peak_params(:,2);
                    fooof_group_results_org{g}{k}(i_ind + i).r_squared = fooof_results{g}{cond,group}{i}.r_squared;
                    %------------ Compute average power after flatten curve
%                     fooof_diff = fooof_results{g}{cond,group}{i}.power_spectrum - fooof_results{g}{cond,group}{i}.ap_fit;
                    % Super important, the output is already logged, the
                    % only difference is the magnitude by 10
                    fooof_diff = 10*(fooof_results{g}{cond,group}{i}.power_spectrum) - 10*(fooof_results{g}{cond,group}{i}.ap_fit);
                    fooof_freq = fooof_results{g}{cond,group}{i}.freqs;
                    fooof_group_results_org{g}{k}(i_ind + i).theta_avg_power = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_avg_power = mean(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).beta_avg_power = mean(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));

                    % data structure needs to be freq x subject
                    fooof_diff_store{g}{k}{cond}(:,i) = fooof_diff';
                    fooof_apfit_store{g}{k}{cond}(:,i) = 10*(fooof_results{g}{cond,group}{i}.ap_fit);
                    
                    % - store original spec data
                    spec_data_original{g}{k}{cond} = specdata{cond,group}(specfreqs >= f_range(1) & specfreqs <= f_range(2),:);
                end
                i_ind = i_ind + size(specdata{cond,group},2);

            end
        end
    end
end
save([save_dir filesep 'fooof_results.mat'],'fooof_results');
%% Determine peak occurs at different band, don't think I am using this
% Looks like 
% Define frequency bands of interest
for g = DESIGN_I
    for k = 3:length(fooof_group_results_org{g})
        for i = 1:length(fooof_group_results_org{g}{k})
            if ~isempty(fooof_group_results_org{g}{k}(i).central_freq)
                fooof_group_results_org{g}{k}(i).theta = [];
                fooof_group_results_org{g}{k}(i).alpha = [];
                fooof_group_results_org{g}{k}(i).beta = [];
                for j = 1:length(fooof_group_results_org{g}{k}(i).central_freq)
                    cf = fooof_group_results_org{g}{k}(i).central_freq(j);
                    if cf > theta_band(1) & cf <= theta_band(2)
                        fooof_group_results_org{g}{k}(i).theta = [fooof_group_results_org{g}{k}(i).theta; cf fooof_group_results_org{g}{k}(i).power(j)];
                    elseif cf > alpha_band(1) & cf <= alpha_band(2)
                        fooof_group_results_org{g}{k}(i).alpha = [fooof_group_results_org{g}{k}(i).alpha; cf fooof_group_results_org{g}{k}(i).power(j)];
                    elseif cf > beta_band(1) & cf <= beta_band(2)
                        fooof_group_results_org{g}{k}(i).beta = [fooof_group_results_org{g}{k}(i).beta; cf fooof_group_results_org{g}{k}(i).power(j)];
                    end
                end
                if length(fooof_group_results_org{g}{k}(i).theta) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).theta(:,1)-6));
                    temp_power = fooof_group_results_org{g}{k}(i).theta(:,2);
                    fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).theta(indx,1) temp_power(indx)];
                end
                if length(fooof_group_results_org{g}{k}(i).alpha) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).alpha(:,1)-10));
                    temp_power = fooof_group_results_org{g}{k}(i).alpha(:,2);
                    fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).alpha(indx,1) temp_power(indx)];
                end
                if length(fooof_group_results_org{g}{k}(i).beta) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).beta(:,1)-20));
                    temp_power = fooof_group_results_org{g}{k}(i).beta(:,2);
                    fooof_group_results_org{g}{k}(i).beta_center = [fooof_group_results_org{g}{k}(i).beta(indx,1) temp_power(indx)];
                end
            end
        end
    end
end
%## SAVE DATA
save([save_dir filesep 'fooof_results_summary.mat'],'fooof_group_results_org');
save([save_dir filesep 'fooof_diff_store.mat'],'fooof_diff_store');
save([save_dir filesep 'fooof_apfit_store.mat'],'fooof_apfit_store');
save([save_dir filesep 'spec_data_original.mat'],'spec_data_original');
%%
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
for k = CLUSTER_PICKS % For each cluster requested
    if isempty(tmp_study.cluster(k).topo) || RE_CALC
        sets_i = unique([fooof_group_results_org{g}{k}(:).subID]);
        tmp_study.cluster(k).sets = tmp_study.cluster(k).sets(sets_i);
        tmp_study.cluster(k).comps = tmp_study.cluster(k).comps(sets_i);
        tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,k);% Using this custom modified code to allow taking average within participant for each cluster
        STUDY.cluster(k).topox = tmp_study.cluster(k).topox;
        STUDY.cluster(k).topoy = tmp_study.cluster(k).topoy;
        STUDY.cluster(k).topoall = tmp_study.cluster(k).topoall;
        STUDY.cluster(k).topo = tmp_study.cluster(k).topo;
        STUDY.cluster(k).topopol = tmp_study.cluster(k).topopol;
    end
end
    
%% Create table from group results, take mean across participants ICs    
psd_feature_table = table;
for g = DESIGN_I
    for k = 3:length(fooof_group_results_org{g})
        temp_table_C1 = fooof_group_results_org{g}{k};
        for i = 1:length(temp_table_C1)
            if isempty(temp_table_C1(i).alpha)
                temp_table_C1(i).alpha = [NaN, NaN];
            end
            if isempty(temp_table_C1(i).beta)
                temp_table_C1(i).beta = [NaN, NaN];
            end
            if isempty(temp_table_C1(i).alpha_center)
                temp_table_C1(i).alpha_center = [NaN, NaN];
            end

            if isempty(temp_table_C1(i).beta_center)
                temp_table_C1(i).beta_center = [NaN, NaN];
            end

            [~,idx_a] = max(temp_table_C1(i).alpha(:,2));
            [~,idx_b] = max(temp_table_C1(i).beta(:,2));

            psd_feature_table = vertcat(psd_feature_table,table(temp_table_C1(i).subID,temp_table_C1(i).compID,temp_table_C1(i).study,...
                temp_table_C1(i).cond,temp_table_C1(i).group,temp_table_C1(i).cluster,temp_table_C1(i).sub_char,temp_table_C1(i).speed_ms,temp_table_C1(i).aperiodic_exp,temp_table_C1(i).aperiodic_offset,...
                temp_table_C1(i).r_squared,temp_table_C1(i).alpha(idx_a,1),temp_table_C1(i).alpha(idx_a,2),temp_table_C1(i).beta(idx_b,1),temp_table_C1(i).beta(idx_b,2),...
                temp_table_C1(i).alpha_center(1,1),temp_table_C1(i).alpha_center(1,2),temp_table_C1(i).beta_center(1,1),temp_table_C1(i).beta_center(1,2),...
                temp_table_C1(i).theta_avg_power(1,1),temp_table_C1(i).alpha_avg_power(1,1),temp_table_C1(i).beta_avg_power(1,1),...
                'VariableNames',...
                {'subID','compID','study','cond','group','cluster','subj_char','speed_ms','aperiodic_exp','aperiodic_offset','r_squared','alpha_cf','alpha_p',...
                'beta_cf','beta_p','alpha_center','alpha_centerP','beta_center','beta_centerP',...
                'theta_avg_power','alpha_avg_power','beta_avg_power'}));
        end
    end
end
psd_feature_table.subID = categorical(psd_feature_table.subID);
psd_feature_table.cond = categorical(psd_feature_table.cond);
psd_feature_table.cluster = categorical(psd_feature_table.cluster);
psd_feature_table.study = categorical(psd_feature_table.study);

% C1_table.study = categorical(C1_table.study);

grp_C1_table = grpstats(psd_feature_table,["subID","study","cond","cluster"],'nanmedian','DataVars',["r_squared","alpha_cf","alpha_p",...
            "beta_cf","beta_p","alpha_center","alpha_centerP","beta_center","beta_centerP","theta_avg_power","alpha_avg_power","beta_avg_power"]);
grp_C1_table.Properties.VariableNames = {'subID','study','cond','cluster','GroupCount','med_r_squared','med_alpha_cf','med_alpha_p',...
            'med_beta_cf','med_beta_p','med_alpha_center','med_alpha_centerP','med_beta_center','med_beta_centerP',...
            'med_theta_avg_power','med_alpha_avg_power','med_beta_avg_power'};
writetable(psd_feature_table,[save_dir filesep 'fooof_spec_table.xlsx']);
save([save_dir filesep 'psd_feature_table.mat'],'psd_feature_table');
%% (STATISTICS CALCULATIONS) =========================================== %% 
%## Preliminary stats on average power (substracted background)
% Create stats table
psd_feature_stats = table;
for g = DESIGN_I
    for k = 3:length(fooof_diff_store{g})
        %- extract table for cluster and design
        t1 = psd_feature_table(psd_feature_table.cluster == num2str(k) & psd_feature_table.study == num2str(g) ,:);
        t1.cond = categorical(t1.cond);
        t1.log_alpha_avg_power = log(t1.alpha_avg_power+5);
        t1.log_theta_avg_power = log(t1.theta_avg_power+5);
        % t2 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,'2')& grp_C1_table.cond ~= 1,:);
        % t2.cond = categorical(t2.cond);
        %- theta
        lme_theta_avg_power = fitlme(t1,'theta_avg_power ~ cond + (1|subID)');
%         lme_theta_avg_power = fitlme(t1,'theta_avg_power ~ cond + (1|speed_ms)');
        theta_anova = anova(lme_theta_avg_power); [theta_h,theta_p] = lillietest(lme_theta_avg_power.residuals);
        %- alpha
        lme_alpha_avg_power= fitlme(t1,'alpha_avg_power ~ cond + (1|subID)');
%         lme_alpha_avg_power = fitlme(t1,'alpha_avg_power ~ cond + (1|speed_ms)');
        alpha_anova = anova(lme_alpha_avg_power); [alpha_h,alpha_p] = lillietest(lme_alpha_avg_power.residuals);
        %- beta
        lme_beta_avg_power = fitlme(t1,'beta_avg_power ~ cond + (1|subID)');
%         lme_beta_avg_power = fitlme(t1,'beta_avg_power ~ cond + (1|speed_ms)');
        beta_anova = anova(lme_beta_avg_power); [beta_h,beta_p] = lillietest(lme_beta_avg_power.residuals);
%         lme_log_theta_avg_power = fitlme(t1,'log_theta_avg_power ~ cond + (1|subID)');
%         log_Th = anova(lme_log_theta_avg_power);[h,p] = lillietest(lme_log_theta_avg_power.residuals)
%         lme_log_alpha_avg_power= fitlme(t1,'log_alpha_avg_power ~ cond + (1|subID)');
%         log_A = anova(lme_log_alpha_avg_power);[h,p] = lillietest(lme_log_alpha_avg_power.residuals)
        %- add stats for aperiod fit/offsets
        lme_ap_exp = fitlme(t1,'aperiodic_exp ~ cond + (1|speed_ms)');
        ap_exp = anova(lme_ap_exp);
        lme_ap_offset = fitlme(t1,'aperiodic_offset ~ cond + (1|speed_ms)');
        ap_offset = anova(lme_ap_offset);
        %-
        temp_stats = table;
        temp_stats.theta_lilnorm_h = theta_h;
        temp_stats.theta_lilnorm_p = theta_p;
        temp_stats.alpha_lilnorm_h = alpha_h;
        temp_stats.alpha_lilnorm_p = alpha_p;
        temp_stats.beta_lilnorm_h = beta_h;
        temp_stats.beta_lilnorm_p = beta_p;
        %- statistics & labels
        temp_stats.study = g;
        temp_stats.cluster = k;
        temp_stats.theta_anova = theta_anova.pValue(2);
        temp_stats.theta_F = theta_anova.FStat(2);
        temp_stats.theta_F_DF2 = theta_anova.DF2(2);
        temp_stats.alpha_anova = alpha_anova.pValue(2);
        temp_stats.alpha_F = alpha_anova.FStat(2);
        temp_stats.alpha_F_DF2 = alpha_anova.DF2(2);
        temp_stats.beta_anova = beta_anova.pValue(2);
        temp_stats.beta_F = beta_anova.FStat(2);
        temp_stats.beta_F_DF2 = beta_anova.DF2(2);
        %- pvalue
        temp_stats.theta_cond2_pval = lme_theta_avg_power.Coefficients.pValue(2);
        temp_stats.theta_cond3_pval = lme_theta_avg_power.Coefficients.pValue(3);
        temp_stats.theta_cond4_pval = lme_theta_avg_power.Coefficients.pValue(4);

        temp_stats.alpha_cond2_pval = lme_alpha_avg_power.Coefficients.pValue(2);
        temp_stats.alpha_cond3_pval = lme_alpha_avg_power.Coefficients.pValue(3);
        temp_stats.alpha_cond4_pval = lme_alpha_avg_power.Coefficients.pValue(4);

        temp_stats.beta_cond2_pval = lme_beta_avg_power.Coefficients.pValue(2);
        temp_stats.beta_cond3_pval = lme_beta_avg_power.Coefficients.pValue(3);
        temp_stats.beta_cond4_pval = lme_beta_avg_power.Coefficients.pValue(4);
        %- aperiodic fit
        temp_stats.ap_exp_anova = ap_exp.pValue(2);
        temp_stats.ap_exp2 = lme_ap_exp.Coefficients.pValue(2);
        temp_stats.ap_exp3 = lme_ap_exp.Coefficients.pValue(3);
        temp_stats.ap_exp4 = lme_ap_exp.Coefficients.pValue(4);
        %- aperiodic offset
        temp_stats.ap_offset_anova = ap_offset.pValue(2);
        temp_stats.ap_offset2 = lme_ap_offset.Coefficients.pValue(2);
        temp_stats.ap_offset3 = lme_ap_offset.Coefficients.pValue(3);
        temp_stats.ap_offset4 = lme_ap_offset.Coefficients.pValue(4);
        
       
        if g == 2
            % use continuous variable
            t1.cond = double(t1.cond)*0.25;
            %- theta
%             lme_theta_avg_power_num = fitlme(t1,'theta_avg_power ~ cond + (1|speed_ms)');
            lme_theta_avg_power_num = fitlme(t1,'theta_avg_power ~ cond + (1|subID)');
            Th_num = anova(lme_theta_avg_power_num);[theta_h,theta_p] = lillietest(lme_theta_avg_power_num.residuals);
            %- alpha
%             lme_alpha_avg_power_num = fitlme(t1,'alpha_avg_power ~ cond + (1|speed_ms)');
            lme_alpha_avg_power_num = fitlme(t1,'alpha_avg_power ~ cond  + (1|subID)');
            A_num  = anova(lme_alpha_avg_power_num);[alpha_h,alpha_p] = lillietest(lme_alpha_avg_power_num.residuals);
            %- beta
%             lme_beta_avg_power_num = fitlme(t1,'beta_avg_power ~ cond + (1|speed_ms)');
            lme_beta_avg_power_num = fitlme(t1,'beta_avg_power ~ cond + (1|subID)');
            B_num  = anova(lme_beta_avg_power_num);[beta_h,beta_p] = lillietest(lme_beta_avg_power_num.residuals);
            %- lillie test
            temp_stats.theta_lilnorm_h = theta_h;
            temp_stats.theta_lilnorm_p = theta_p;
            temp_stats.alpha_lilnorm_h = alpha_h;
            temp_stats.alpha_lilnorm_p = alpha_p;
            temp_stats.beta_lilnorm_h = beta_h;
            temp_stats.beta_lilnorm_p = beta_p;
            %- 
            temp_stats.Th_num_pval = Th_num.pValue(2);
            temp_stats.A_num_pval = A_num.pValue(2);
            temp_stats.B_num_pval = B_num.pValue(2);
            temp_stats.theta_F_num = Th_num.FStat(2);
            temp_stats.theta_F_DF2_num = Th_num.DF2(2);
            temp_stats.alpha_F_num = A_num.FStat(2);
            temp_stats.alpha_F_DF2_num = A_num.DF2(2);
            temp_stats.beta_F_num = B_num.FStat(2);
            temp_stats.beta_F_DF2_num = B_num.DF2(2);
            
            % store the linear fit for all these data
            temp_stats.Th_slope = lme_theta_avg_power_num.Coefficients.Estimate(2);
            temp_stats.Th_intercept = lme_theta_avg_power_num.Coefficients.Estimate(1);
            temp_stats.A_slope = lme_alpha_avg_power_num.Coefficients.Estimate(2);
            temp_stats.A_intercept = lme_alpha_avg_power_num.Coefficients.Estimate(1);
            temp_stats.B_slope = lme_beta_avg_power_num.Coefficients.Estimate(2);
            temp_stats.B_intercept = lme_beta_avg_power_num.Coefficients.Estimate(1);
            temp_stats.Th_num_R2 = lme_theta_avg_power_num.Rsquared.Adjusted;
            temp_stats.A_num_R2 = lme_alpha_avg_power_num.Rsquared.Adjusted;
            temp_stats.B_num_R2 = lme_beta_avg_power_num.Rsquared.Adjusted;
        else
            %- lillie test
            temp_stats.theta_lilnorm_h = NaN;
            temp_stats.theta_lilnorm_p = NaN;
            temp_stats.alpha_lilnorm_h = NaN;
            temp_stats.alpha_lilnorm_p = NaN;
            temp_stats.beta_lilnorm_h = NaN;
            temp_stats.beta_lilnorm_p = NaN;
            
            temp_stats.Th_num_pval = NaN;
            temp_stats.A_num_pval = NaN;
            temp_stats.B_num_pval = NaN;
            temp_stats.theta_F_num = NaN;
            temp_stats.theta_F_DF2_num = NaN;
            temp_stats.alpha_F_num = NaN;
            temp_stats.alpha_F_DF2_num = NaN;
            temp_stats.beta_F_num = NaN;
            temp_stats.beta_F_DF2_num = NaN;
            
            temp_stats.Th_slope = NaN;
            temp_stats.Th_intercept = NaN;
            temp_stats.A_slope = NaN;
            temp_stats.A_intercept = NaN;
            temp_stats.B_slope = NaN;
            temp_stats.B_intercept = NaN;
            
            temp_stats.Th_num_R2 = NaN;
            temp_stats.A_num_R2 = NaN;
            temp_stats.B_num_R2 = NaN;
        end
        
        psd_feature_stats = vertcat(psd_feature_stats,temp_stats);
    end
end
psd_feature_stats.study = categorical(psd_feature_stats.study);
psd_feature_stats.cluster = categorical(psd_feature_stats.cluster);
writetable(psd_feature_stats,[save_dir filesep 'psd_band_power_stats.xlsx']);
save([save_dir filesep 'psd_band_power_stats.mat'],'psd_feature_stats');
%% Perform time series stats on the flattened curve
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
% fooof_diff_store needs to be freq x subject, and condition by row
design_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
cnt = 1;
for g = DESIGN_I
    for k = 3:length(fooof_diff_store{g})
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(fooof_diff_store{g}{k}', stats);
        %- 
        clust_t{cnt} = k;
        design_t{cnt} = g;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter;
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond{g}{k} = temp_pcond;
        pgroup{g}{k} = temp_pcond;
        pinter{g}{k} = temp_pinter;
        statcond{g}{k} = temp_statcond;
        statgroup{g}{k} = temp_statgroup;
        statinter{g}{k} = temp_statinter;
        for k0 = 1:length(pcond{g}{k})
            pcond{g}{k}{k0}(:,2) = pcond{g}{k}{k0}(:,1)<0.05;    
            pcond_test_t{cnt} = pcond{g}{k}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup{g}{k})
            if ~isempty(pgroup{g}{k}{k0})
                pgroup{g}{k}{k0}(:,2) = pgroup{g}{k}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup{g}{k}{k0}(:,1)<0.05;    
            end
        end
        for k0 = 1:length(pinter{g}{k})
            if ~isempty(pinter{g}{k}{k0})
                pinter{g}{k}{k0}(:,2) = pinter{g}{k}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter{g}{k}{k0}(:,1)<0.05;    
            end
        end
        cnt = cnt + 1;
    end
end
% table_out = table(design_t,clust_t,pcond_t,pcond_test_t,pgroup_t,pgroup_test_t,pinter_t,pinter_test_t,statcond_t,statgroup_t,statinter_t);
% save([save_dir filesep 'fooof_psd_stats.mat'],'table_out');
save([save_dir filesep 'fooof_pcond.mat'],'pcond');
%% ===================================================================== %%
% fooof_diff_store needs to be freq x subject, and condition by row
freq_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
% subj_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
% cond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
design_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
cnt = 1;
for g = DESIGN_I
    for k = 3:length(spec_data_original{g})
%         k
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(spec_data_original{g}{k}', stats);
        %- 
        clust_t{cnt} = k;
        design_t{cnt} = g;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter; %{temp_pinter};
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond_org{g}{k} = temp_pcond;
        pgroup_org{g}{k} = temp_pcond;
        pinter_org{g}{k} = temp_pinter;
        statcond_org{g}{k} = temp_statcond;
        statgroup_org{g}{k} = temp_statgroup;
        statinter_org{g}{k} = temp_statinter;
        for k0 = 1:length(pcond_org{g}{k})
            pcond_org{g}{k}{k0}(:,2) = pcond_org{g}{k}{k0}(:,1)<0.05; 
            pcond_test_t{cnt} = pcond_org{g}{k}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup_org{g}{k})
            if ~isempty(pgroup_org{g}{k}{k0})
                pgroup_org{g}{k}{k0}(:,2) = pgroup_org{g}{k}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup_org{g}{k}{k0}(:,1)<0.05;   
            end
        end
        for k0 = 1:length(pinter_org{g}{k})
            if ~isempty(pinter_org{g}{k}{k0})
                pinter_org{g}{k}{k0}(:,2) = pinter_org{g}{k}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter_org{g}{k}{k0}(:,1)<0.05; %{pinter_org{g}{k}{k0}(:,1)<0.05};  
            end
        end            
        cnt=cnt+1;
    end
end
% table_out = table(design_t,clust_t,pcond_t,pcond_test_t,pgroup_t,pgroup_test_t,pinter_t,pinter_test_t,statcond_t,statgroup_t,statinter_t);
% writetable(table_out,[save_path filesep 'original_psd_stats.xlsx']);
% save([save_dir filesep 'original_psd_stats.mat'],'table_out');
save([save_dir filesep 'org_pcond.mat'],'pcond_org');
% %% LOAD
% tmp = load([save_dir filesep 'org_pcond.mat']);
% pcond_org = tmp.pcond_org;
%% (LOAD EXISTING TALBES && FORMAT STUDY)
tmp = load([save_dir filesep 'psd_feature_table.mat']);
psd_feature_table = tmp.psd_feature_table;
tmp = load([save_dir filesep 'psd_band_power_stats.mat']);
psd_feature_stats = tmp.psd_feature_stats;
tmp = load([save_dir filesep 'fooof_pcond.mat']);
pcond = tmp.pcond;
tmp = load([save_dir filesep 'org_pcond.mat']);
pcond_org = tmp.pcond_org;
tmp = load([save_dir filesep 'fooof_results_summary.mat']);
fooof_group_results_org = tmp.fooof_group_results_org;
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
for g = DESIGN_I
    for k = CLUSTER_PICKS % For each cluster requested
        if isempty(tmp_study.cluster(k).topo) || RE_CALC
            sets_i = unique([fooof_group_results_org{g}{k}(:).subID]);
            tmp_study.cluster(k).sets = tmp_study.cluster(k).sets(sets_i);
            tmp_study.cluster(k).comps = tmp_study.cluster(k).comps(sets_i);
            tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,k);% Using this custom modified code to allow taking average within participant for each cluster
            STUDY.cluster(k).topox = tmp_study.cluster(k).topox;
            STUDY.cluster(k).topoy = tmp_study.cluster(k).topoy;
            STUDY.cluster(k).topoall = tmp_study.cluster(k).topoall;
            STUDY.cluster(k).topo = tmp_study.cluster(k).topo;
            STUDY.cluster(k).topopol = tmp_study.cluster(k).topopol;
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

%% Sanity check - time series plots from aperiodic subtraction
for g = DESIGN_I
    data_min = [min(fooof_diff_store{g}{k}{1});min(fooof_diff_store{g}{k}{2});min(fooof_diff_store{g}{k}{3});...
    min(fooof_diff_store{g}{k}{4});];
    data_max = [max(fooof_diff_store{g}{k}{1});min(fooof_diff_store{g}{k}{2});max(fooof_diff_store{g}{k}{3});...
    max(fooof_diff_store{g}{k}{4});];
    for k = CLUSTER_PICKS
        figure('color','white');
        hold on;
        for i = 1:4
            subplot(4,1,i)
            data = fooof_diff_store{g}{k}{i};
            plot(fooof_freq,data,'color',color.terrain(i,:));
            ylabel('log10(Power)')
            ylim([min(data_min,[],'all') max(data_max,[],'all')]);
        end
        xlabel('Frequency(Hz)');
        title(['Cluster ',num2str(k)]);
        hold off;
    end
end
%% Sanity check: Plot distribution of aperiodic params (exp), central frequency, and goodness of fit
for g = DESIGN_I
    for k = CLUSTER_PICKS
        temp_table = psd_feature_table(psd_feature_table.cluster == num2str(k) & psd_feature_table.study == num2str(g),: );
        temp_table.cond = categorical(temp_table.cond);
        figure();set(gcf,'color','white');
        subplot(2,2,1)
        boxchart(temp_table.cond,temp_table.aperiodic_exp);
        ylabel('Aperodic exponent');

    %     subplot(2,2,2)
    %     boxchart(temp_table.cond,temp_table.center_frequency);
    %     xlabel('Central frequency');ylabel('# Peaks')

        subplot(2,2,3)
        boxchart(temp_table.cond,temp_table.r_squared);
        ylabel('R squared');
    end
end
%% (Paper Figures) ===================================================== %%
% NOTE: run the above code cells before running this
PSD_YLIM_ORG = [-30,-10];
PSD_YLIM_FOOOF = [-0.5,6];
VIOLIN_YLIM_APERIODIC_EXP = [0,5];
VIOLIN_YLIM_APERIODIC_OFFSET = [-20,5];
% PSD_YLIM_FOOOF = [-0.5,6];
% CINGULATE_I = 3;
% CINGULATE_I = {};
% CINGULATE_I = {3};
% CINGULATE_I = {7};
% CINGULATE_I = {11};
% CINGULATE_I = {14};
% CINGULATE_I = {9};
% CINGULATE_I = {10};
CINGULATE_I = {11,13};
theta_3=[-1,6];
alpha_3=[-1.5,6];
beta_3=[-1,9];
psd_ylim_fooof_3 = [];
psd_ylim_orig_3 = [];
%-
% SENSORIMOTOR_I = {7,8};
% SENSORIMOTOR_I = {6,12};
% SENSORIMOTOR_I = {6,10};
% SENSORIMOTOR_I = {9,11};
% SENSORIMOTOR_I = {5,14};
% SENSORIMOTOR_I = {13,8};
% SENSORIMOTOR_I = {6,13};
% SENSORIMOTOR_I = {5,6};
SENSORIMOTOR_I = {4,14};
theta_1=[-1.5,6.8];
alpha_1=[-1.5,16];
beta_1=[-1,10];
psd_ylim_fooof_1 = [];
psd_ylim_orig_1 = [];
%-
% POSTERIORP_I = {12,5};
% POSTERIORP_I = {11,14};
% POSTERIORP_I = {4,8};
% POSTERIORP_I = {8,14};
% POSTERIORP_I = {3,12};
% POSTERIORP_I = {12,7};
% POSTERIORP_I = {11,4};
% POSTERIORP_I = {14,12};
POSTERIORP_I = {10,9};
theta_2=[-1,6];
alpha_2=[-1.5,15.5];
beta_2=[-1,9];
psd_ylim_fooof_2 = [];
psd_ylim_orig_2 = [];
%-
% SUPPMOTOR_I = {7,10};
% SUPPMOTOR_I = {7};
% SUPPMOTOR_I = {12};
% SUPPMOTOR_I = {6};
% SUPPMOTOR_I = {};
% SUPPMOTOR_I = {3,12};
% SUPPMOTOR_I = {};
SUPPMOTOR_I = {5};
theta_4=[-1,6];
alpha_4=[-2,13];
beta_4=[-1,10];
%-
% OCCIPITAL_I = {4,10};
% OCCIPITAL_I = {5,8};
% OCCIPITAL_I = {9,11};
% OCCIPITAL_I = {10};
% OCCIPITAL_I = {7,13};
% OCCIPITAL_I = {4,6};
% OCCIPITAL_I = {5};
% OCCIPITAL_I = {7,8};
OCCIPITAL_I = {8};
theta_6=[-1,8];
alpha_6=[-1,15];
beta_6=[-1,7];
psd_ylim_fooof_6 = [];
psd_ylim_orig_6 = [];
%-
% CAUDATE_I = 9;
% CAUDATE_I = 13;
% CAUDATE_I = 14;
% CAUDATE_I = 13;
% CAUDATE_I = 9;
% CAUDATE_I = {9};
% CAUDATE_I = {8};
% CAUDATE_I = {13};
CAUDATE_I = {3};
theta_7=[-1,5.5];
alpha_7=[-1.5,9];
beta_7=[-1,7];
psd_ylim_fooof_7 = [];
psd_ylim_orig_7 = [];
%-
% CUNEUS_I = 13;
% CUNEUS_I = 3;
% CUNEUS_I = 5;
% CUNEUS_I = 5;
% CUNEUS_I = 10;
% CUNEUS_I = {5};
% CUNEUS_I = {14};
CUNEUS_I = {7};
theta_5=[-1,8];
alpha_5=[-1.5,17];
beta_5=[-1,8];
psd_ylim_fooof_5 = [];
psd_ylim_orig_5 = [];
%-
% TEMPORAL_I = {12,13};
% TEMPORAL_I = {4,3};
% TEMPORAL_I = {4,8};
% TEMPORAL_I = {3,11};
% TEMPORAL_I = {7,10};
TEMPORAL_I = {6,12};
theta_8=[-1,5];
alpha_8=[-1.5,12];
beta_8=[-1,7];
psd_ylim_fooof_8 = [];
psd_ylim_orig_8 = [];
%% (OA PAPER FIGURE CLUSTER SPEED/TERRAIN FOR EACH CLUSTER) ============ %%
%{
%## PARAMS
%-
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%-
SUBPLOT_BOTTOM = 0.2;
SUBPLOT_INIT_SHIFT = 0.05;
PLOT_STRUCT = struct('figure_position_inch',[3,3,6,9.5],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'clim',[],...
    'font_size',10,...
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
% FIG_POS = [100 300 1420 250];
%## COLORS
% g = 1;
%- colors
COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2,3,4,5,6,7,8,9,10];
COLOR_LIM_ERR = 0.05;
%% ===================================================================== %%
CLUSTERS_TO_PLOT = main_cl_inds(2:end); %valid_clusters(1:end-1); %main_cl_inds(2:end);
for k_i = 1:length(CLUSTERS_TO_PLOT)
    k = CLUSTERS_TO_PLOT(k_i);
    %## ANATOMY
    dip1 = STUDY.cluster(k).all_diplocs;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'single';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
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
    %## AXES LIMITS
    switch k
        case SENSORIMOTOR_I % sensorimotor area
            ylim_value_theta = theta_1;
            ylim_value_alpha = alpha_1;
            ylim_value_beta  = beta_1;
        case POSTERIORP_I % posterior area
            ylim_value_theta = theta_2;
            ylim_value_alpha = alpha_2;
            ylim_value_beta  = beta_2;
        case CINGULATE_I % cingulate
            ylim_value_theta = theta_3;
            ylim_value_alpha = alpha_3;
            ylim_value_beta  = beta_3;
        case SUPPMOTOR_I % supplementary motor
            ylim_value_theta = theta_4;
            ylim_value_alpha = alpha_4;
            ylim_value_beta  = beta_4;
        case CUNEUS_I
            ylim_value_theta = theta_5;
            ylim_value_alpha = alpha_5;
            ylim_value_beta  = beta_5;
        case OCCIPITAL_I % occipital
            ylim_value_theta = theta_6;
            ylim_value_alpha = alpha_6;
            ylim_value_beta  = beta_6;
        case CAUDATE_I % caudate
            ylim_value_theta = theta_7;
            ylim_value_alpha = alpha_7;
            ylim_value_beta  = beta_7;
        case TEMPORAL_I
            ylim_value_theta = theta_8;
            ylim_value_alpha = alpha_8;
            ylim_value_beta  = beta_8;
    end
    %% ================================================================= %%
    for g = 1:2
        switch g
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
        end
        
        fig = figure('color','none','renderer','Painters');
        sgtitle(atlas_name,'FontName',PLOT_STRUCT.font_name,'FontSize',14,'FontWeight','bold','Interpreter','none');
        set(fig,'Units','inches','Position',[0.5,0.5,6.5,9])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        hold on;
        %## dipole plots (terrain)
        %-
        dimw = 5816;
        dimh = 4715;
        subplot(3,4,2)
        im = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',k)]);
        im = imresize(im,0.8);
        imshow(im);
        ax = gca;
        
        im_resize1 = 0.75;
        set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[0.05,0.75,0.3*im_resize,0.25*im_resize]);
        set(ax,'Position',[0.30,0.75,0.3*im_resize1,0.25*im_resize1]); %[left bottom width height]
        %-
        subplot(3,4,3)
        im = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',k)]);
        im = imresize(im,0.75);
        imshow(im);
        ax = gca;
        im_resize2 = 0.74;
        set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[0.35,0.75,0.25*im_resize,0.25*im_resize]); 
        set(ax,'Position',[0.32+0.32*im_resize1,0.75,0.25*im_resize1,0.25*im_resize1]);  %[left bottom width height]
        %-
        dimw = 5816;
        dimh = 4715;
        subplot(3,4,4)
        im = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',k)]);
        im = imresize(im,0.75);
        imshow(im);
        ax = gca;
        
        im_resize = 0.75;
        set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[0.6,0.75,0.358*im_resize,0.25*im_resize]);
        set(ax,'Position',[0.3+0.3*im_resize+0.25*im_resize,0.75,0.358*im_resize,0.25*im_resize]); %[left bottom width height]

        %## topo plot 
        im_resize = 0.8;
%         subplot(3,3,4)
        subplot(3,4,1)
        std_topoplot_CL(STUDY,k,'together');
        colormap(linspecer); %colormap_ersp)
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'color','w')
%         fig_i.Children(1).Title.String = sprintf('CL%i N=%i',k,length(STUDY.cluster(k).sets));
        fig_i.Children(1).Title.String = sprintf('N=%i',length(STUDY.cluster(k).sets));
        fig_i.Children(1).Title.Interpreter = 'none';
        fig_i.Children(1).FontSize = 12; %PLOT_STRUCT.font_size;
        fig_i.Children(1).FontName = PLOT_STRUCT.font_name;
        fig_i.Children(1).FontWeight = 'bold';
%         fig_i.Children(1).OuterPosition = [0,0,1,1];
%         fig_i.Children(1).Position = [0.05,0.475,0.225,0.25];  %[left bottom width height]
        fig_i.Children(1).OuterPosition = [0,0,1,1];
        fig_i.Children(1).Position = [0.08,0.725,0.225*im_resize,0.25*im_resize];  %[left bottom width height]
        %%
        im_resize = 1;
        %## non-fooof psd
%         subplot(3,3,5)
        subplot(4,4,5)
        hold on;
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));
        end
        axs = [];
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end 
        % plot the aperiodic line
        for i = 1:length(spec_data_original{g}{k})
            aperiodic_fit = fooof_apfit_store{g}{k}{i}';
            dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',1,'displayname','aperiodic fit');
        end
        ax = gca;
        xlim([4 40]);ylim([-30 -5]);
        axsignif = highlight_CL(ax, fooof_freq, pcond_org{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
        ylabel('10*log_{10}(Power)');
        xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        title('PSD')
%         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[0.12,0.475,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
        hold off;
        %## fooof psd (terrain)
%         subplot(3,3,6)
        subplot(3,4,7)
        hold on;
        axs = [];
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));% input need to be a row vector
        end
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end  
        ax = gca;
        axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        xlim([4 40]);
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
%         ylabel('10*log_{10}(Power)');
        xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        legend([axs,dash],'fontsize',9,'fontname',PLOT_STRUCT.font_name);
        lg1 = legend('boxoff');
        
        set(lg1,'Position',[0.8,0.6,0.2,0.1]);
        title('Flattened PSD')
    %     set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[0.25+0.25*im_resize,0.475,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
        hold off;
        %%
        %## violin plot's theta/alpha/beta (terrain)
        horiz_shift= 0;
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            T_plot = psd_feature_table(psd_feature_table.study == num2str(g) & psd_feature_table.cluster == num2str(k),:);
%             subplot(3,3,6+i)
            if i>1
                subplot(3,4,9+i)
            else
                subplot(3,4,9+i)
            end
            hold on;
            violinplot(T_plot.(measure_name),T_plot.cond,...
                'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
                'MarkerSize',15,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin',...
                'full','width',0.3,'QuartileStyle','shadow');
            ylabel('');
            xlabel('');
            if i == 1;ylabel('10*log_{10}(Flattened PSD)');else;ylabel('');end
            if g == 2&&i==1;xlabel('speed (m/s)');end
            ax = gca;
            set(ax,'xticklabel', xtick_label_g,'XTickLabelRotation',45,'FontWeight','bold');
            switch i 
                 case 1
                     ylim(ylim_value_theta);
                 case 2
                     ylim(ylim_value_alpha);
                 case 3
                     ylim(ylim_value_beta);
            end
            % axis padded
            fig_i = get(groot,'CurrentFigure');
            box off
            title(title_plot{i},'FontWeight','bold','FontSize',12);
            T_stats_plot = psd_feature_stats(psd_feature_stats.study == num2str(g) & psd_feature_stats.cluster == num2str(k),:);
            switch i 
                case 1
                    anova_stats = T_stats_plot.theta_anova;
                    cond2_stats = T_stats_plot.theta_cond2_pval;
                    cond3_stats = T_stats_plot.theta_cond3_pval;
                    cond4_stats = T_stats_plot.theta_cond4_pval;
                    regress_sig = T_stats_plot.Th_num_pval;
                    regressline_stats = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                    R2 = T_stats_plot.Th_num_R2;
                case 2
                    anova_stats = T_stats_plot.alpha_anova;
                    cond2_stats = T_stats_plot.alpha_cond2_pval;
                    cond3_stats = T_stats_plot.alpha_cond3_pval;
                    cond4_stats = T_stats_plot.alpha_cond4_pval;
                    regress_sig = T_stats_plot.A_num_pval;
                    regressline_stats = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                    R2 = T_stats_plot.A_num_R2;
                case 3
                    anova_stats = T_stats_plot.beta_anova;
                    cond2_stats = T_stats_plot.beta_cond2_pval;
                    cond3_stats = T_stats_plot.beta_cond3_pval;
                    cond4_stats = T_stats_plot.beta_cond4_pval;
                    regress_sig = T_stats_plot.B_num_pval;
                    regressline_stats = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                    R2 = T_stats_plot.B_num_R2;
            end
            if g == 1 % terrain condition, categorical values
                 if anova_stats < 0.05
                    if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
                    if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
                    if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
                 end
            elseif g == 2 % g == 2 % speed conditions, numeric values, regression line
                if anova_stats < 0.05
                    % plot line
                    x = 0:5;
                    y = x*regressline_stats(2)*0.25 + regressline_stats(1);
                    plot(x,y,'-','color','k','linewidth',1);
                    xlim([0 5]);
                    if regress_sig > 0.01 & regress_sig < 0.05
                        text(1,gety(gca)*1.2,{['* ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',9,'FontName',PLOT_STRUCT.font_name);
                    elseif regress_sig <= 0.01 & regress_sig > 0.001 
                        text(1,gety(gca)*1.2,{['** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',9,'FontName',PLOT_STRUCT.font_name);
                    else
                        text(1,gety(gca)*1.2,{['*** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',9,'FontName',PLOT_STRUCT.font_name);
                    end
                end
            end
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
    %         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.1+horiz_shift,0.15,0.25,0.225]);  %[left bottom width height]
            hold off;
            %- iterate
            horiz_shift = horiz_shift + 0.3;
        end
        hold off;
        pause(1)
        %% ================================================================= %%
%         exportgraphics(gcf,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_des%i.pdf',k,g)],'ContentType','vector','Resolution',300)
%         exportgraphics(gcf,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_des%i.jpg',k,g)],'Resolution',300)
        exportgraphics(gcf,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_des%i.tiff',k,g)],'Resolution',1200)
        close(gcf);
    end
    
    
end
%}
%% ===================================================================== %%
%## PARAMS

%-
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
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
% FIG_POS = [100 300 1420 250];
%## COLORS
% g = 1;
%- colors
COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2,3,4,5,6,7,8,9,10];
COLOR_LIM_ERR = 0.05;
%%
CLUSTERS_TO_PLOT = main_cl_inds(2:end);
[STUDY,centroid] = std_centroid(STUDY,ALLEEG,CLUSTERS_TO_PLOT,'dipole');
%{
txt_store = cell(length(CLUSTERS_TO_PLOT),1);
for k_i = 1:length(CLUSTERS_TO_PLOT)
    k = CLUSTERS_TO_PLOT(k_i);
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
    sprintf('CENTROID: Dip Center: [%0.1f,%0.1f,%0.1f]\n\n',STUDY.cluster(k).centroid.dipole.posxyz)];
end
cellfun(@(x) disp(x),txt_store);
%}
%% ===================================================================== %%
CLUSTERS_TO_PLOT = main_cl_inds(2:end); %valid_clusters(1:end-1); %main_cl_inds(2:end);
for k_i = 1:length(CLUSTERS_TO_PLOT)
    k = CLUSTERS_TO_PLOT(k_i);
    %## ANATOMY
    dip1 = STUDY.cluster(k).all_diplocs;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'single';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
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
    %## AXES LIMITS
    switch k
        case SENSORIMOTOR_I % sensorimotor area
            ylim_value_theta = theta_1;
            ylim_value_alpha = alpha_1;
            ylim_value_beta  = beta_1;
        case POSTERIORP_I % posterior area
            ylim_value_theta = theta_2;
            ylim_value_alpha = alpha_2;
            ylim_value_beta  = beta_2;
        case CINGULATE_I % cingulate
            ylim_value_theta = theta_3;
            ylim_value_alpha = alpha_3;
            ylim_value_beta  = beta_3;
        case SUPPMOTOR_I % supplementary motor
            ylim_value_theta = theta_4;
            ylim_value_alpha = alpha_4;
            ylim_value_beta  = beta_4;
        case CUNEUS_I
            ylim_value_theta = theta_5;
            ylim_value_alpha = alpha_5;
            ylim_value_beta  = beta_5;
        case OCCIPITAL_I % occipital
            ylim_value_theta = theta_6;
            ylim_value_alpha = alpha_6;
            ylim_value_beta  = beta_6;
        case CAUDATE_I % caudate
            ylim_value_theta = theta_7;
            ylim_value_alpha = alpha_7;
            ylim_value_beta  = beta_7;
        case TEMPORAL_I
            ylim_value_theta = theta_8;
            ylim_value_alpha = alpha_8;
            ylim_value_beta  = beta_8;
    end
    %% ================================================================= %%
%     for g = 1:2
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
%         subplot(3,3,4)
        subplot(3,4,1)
        std_topoplot_CL(STUDY,k,'together');
        colormap(linspecer); %colormap_ersp)
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'color','w')
%         fig_i.Children(1).Title.String = sprintf('CL%i N=%i',k,length(STUDY.cluster(k).sets));
        fig_i.Children(1).Title.String = sprintf('N=%i',length(STUDY.cluster(k).sets));
        fig_i.Children(1).Title.Interpreter = 'none';
        fig_i.Children(1).FontSize = 12; %PLOT_STRUCT.font_size;
        fig_i.Children(1).FontName = PLOT_STRUCT.font_name;
        fig_i.Children(1).FontWeight = 'bold';
%         fig_i.Children(1).OuterPosition = [0,0,1,1];
%         fig_i.Children(1).Position = [0.05,0.475,0.225,0.25];  %[left bottom width height]
        fig_i.Children(1).OuterPosition = [0,0,1,1];
        fig_i.Children(1).Units = 'Inches';
        fig_i.Children(1).Position = [0.5,BOT_D-0.175,0.225*im_resize,0.25*im_resize];  %[left bottom width height]
        %## dipole plots (terrain)
        %-
        dpi = 1000;
        target_sz = 1.25; %inch
        target_dim = 1; %2==width, 1==height
        im1 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',k)]);
        
        if target_dim==1
            scale = target_sz/(size(im1,1)/dpi);
            im1 = imresize(im1,scale);
        elseif target_dim==2
            scale = target_sz/(size(im1,2)/dpi);
            im1 = imresize(im1,scale);
        end
        
        im2 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',k)]);
        im2(1:300,:,:) = [];
        im2(end-300:end,:,:)=[];
        if target_dim==1
            scale = target_sz/(size(im2,1)/dpi);
            im2 = imresize(im2,scale);
        elseif target_dim==2
            scale = target_sz/(size(im2,2)/dpi);
            im2 = imresize(im2,scale);
        end

        im3 = imread([cluster_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',k)]);
%         im3(1,:,:) = [];
        if target_dim==1
            scale = target_sz/(size(im3,1)/dpi);
            im3 = imresize(im3,scale);
        elseif target_dim==2
            scale = target_sz/(size(im3,2)/dpi);
            im3 = imresize(im3,scale);
        end
        
        %-
        
%         szw = 1.4;
%         szh = 1.25;
        szw1 = (size(im1,2)/dpi); %/PG_SIZE(1); %width
        szh1 = (size(im1,1)/dpi); %/PG_SIZE(2); %+0.0001; %height
        szw2 = (size(im2,2)/dpi); %/PG_SIZE(1); %width
        szh2 = (size(im2,1)/dpi); %+0.05; %/PG_SIZE(2); %+0.01;
        szw3 = (size(im3,2)/dpi); %/PG_SIZE(1);
        szh3 = (size(im3,1)/dpi); %/PG_SIZE(2);
        szw = max([szw1,szw2,szw3]);
        %-
%         subplot(3,8,3);
        axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D,BOT_D,szw,szh1],'PositionConstraint','outerposition');
        imshow(im1,'border','tight');
%         ax = gca;
%         set(ax,'Units','Inches');
%         set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[LEFT_D,BOT_D,szw,szh1]); %[left bottom width height]
%         set(ax,'PositionConstraint','outerposition');
        %-
%         subplot(3,8,5)
        axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+((szw-szw1)-(szw-szw2)),BOT_D,szw,szh2],'PositionConstraint','outerposition');
        imshow(im2,'border','tight');
%         ax = gca;
%         set(ax,'Units','Inches');
%         set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[LEFT_D+szw+GAP,BOT_D,szw,szh2]);  %[left bottom width height]
%         set(ax,'PositionConstraint','outerposition');
        %-
%         subplot(3,8,7)
        axes('Units','Inches','OuterPosition',[0,0,1,1],'Position',[LEFT_D+szw1+szw2+0.05,BOT_D,szw,szh3],'PositionConstraint','outerposition');
        imshow(im3,'border','tight');
%         ax = gca;
%         set(ax,'Units','Inches');
%         set(ax,'OuterPosition',[0,0,1,1]);
%         set(ax,'Position',[LEFT_D+szw+szw+GAP,BOT_D,szw,szh3]); %[left bottom width height]
%         set(ax,'PositionConstraint','outerposition');
       
        %%
        g = 2;
        switch g
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
        end
        im_resize = 0.5;
        %## non-fooof psd
%         subplot(3,3,5)
        subplot(4,8,9)
        hold on;
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));
        end
        axs = [];
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end 
        % plot the aperiodic line
        for i = 1:length(spec_data_original{g}{k})
            aperiodic_fit = fooof_apfit_store{g}{k}{i}';
            dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',1,'displayname','ap. fit');
        end
        ax = gca;
        xlim([4 40]);ylim([-30 -5]);
        axsignif = highlight_CL(ax, fooof_freq, pcond_org{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
        ylabel('10*log_{10}(Power)');
        xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        title('PSD')
%         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[0.08,0.575,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
        hold off;
        %## fooof psd (terrain)
%         subplot(3,3,6)
        subplot(3,8,11)
        hold on;
        axs = [];
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));% input need to be a row vector
        end
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end
        %-
        ax = gca;
        axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        xlim([4 40]);
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
%         ylabel('10*log_{10}(Power)');
        xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        %- legend
        legend([axs,dash],'FontSize',9,'FontName',PLOT_STRUCT.font_name);
        [lg1,icons,plots,txt] = legend('boxoff');
        set(lg1,'Position',[0.22+0.3*im_resize,0.60,0.2,0.1]);
        lg1.ItemTokenSize(1) = 18;
        %-
        title('Flattened PSD')
    %     set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        carry_ov = 0.12+0.3*im_resize;
        set(ax,'Position',[carry_ov,0.575,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
%         icons(2).XData = [0.05 0.1];
        hold off;
        %%
        g = 1;
        switch g
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
        end
        im_resize = 0.5;
        %## non-fooof psd
%         subplot(3,3,5)
        subplot(4,8,13)
        hold on;
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));
        end
        axs = [];
        for i = 1:length(spec_data_original{g}{k})
            data = spec_data_original{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end 
        % plot the aperiodic line
        for i = 1:length(spec_data_original{g}{k})
            aperiodic_fit = fooof_apfit_store{g}{k}{i}';
            dash = plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',1,'displayname','ap. fit');
        end
        ax = gca;
        xlim([4 40]);ylim([-30 -5]);
        axsignif = highlight_CL(ax, fooof_freq, pcond_org{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
%         ylabel('10*log_{10}(Power)');
        xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        title('PSD')
%         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[0.3+carry_ov,0.575,0.3*im_resize,0.25*im_resize]);  %[left bottom width height]
        hold off;
        %## fooof psd (terrain)
%         subplot(3,3,6)
        subplot(3,8,16)
        hold on;
        axs = [];
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                color_dark(i,:),color_light(i,:));% input need to be a row vector
        end
        for i = 1:length(fooof_diff_store{g}{k})
            data = fooof_diff_store{g}{k}{i}';
            ax = plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2,'displayname',sprintf('%s',xtick_label_g{i}));
            axs = [axs, ax];
        end  
        ax = gca;
        axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
        xlim([4 40]);
        plot([0 40],[0 0],'--','color','black');
        xlabel('Frequency(Hz)');
%         ylabel('10*log_{10}(Power)');
        xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    %     set(ax,'LineWidth',2)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
        %- legend
        legend([axs,dash],'FontSize',9,'FontName',PLOT_STRUCT.font_name);
        [lg1,icons,plots,txt] = legend('boxoff');
        set(lg1,'Position',[0.25+0.35*im_resize+carry_ov+0.3*im_resize,0.6,0.2,0.1]);
        lg1.ItemTokenSize(1) = 18;
        title('Flattened PSD')
    %     set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
        set(ax,'OuterPosition',[0,0,1,1]);
        set(ax,'Position',[0.325+carry_ov+0.3*im_resize,0.575,0.35*im_resize,0.25*im_resize]);  %[left bottom width height]
        hold off;
        %%
        im_resize= 0.5;
        g = 2;
        switch g
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
        end
        %## violin plot's theta/alpha/beta (terrain)
        horiz_shift= 0;
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            T_plot = psd_feature_table(psd_feature_table.study == num2str(g) & psd_feature_table.cluster == num2str(k),:);
%             subplot(3,3,6+i)
            if i>1
                subplot(3,8,17+i)
            else
                subplot(3,8,17+i)
            end
            hold on;
            
            
            violinplot(T_plot.(measure_name),T_plot.cond,...
                'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
                'MarkerSize',8,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin',...
                'full','width',0.3,'QuartileStyle','shadow');
            ylabel('');
            xlabel('');
            if i == 1;ylabel('10*log_{10}(Flattened PSD)');else;ylabel('');end
            if g == 2&&i==1;xlabel('speed (m/s)');end
            ax = gca;
            set(ax,'xticklabel', xtick_label_g,'XTickLabelRotation',45,'FontWeight','bold');
%             switch i 
%                  case 1
%                      ylim(ylim_value_theta);
%                  case 3
%                      ylim(ylim_value_alpha);
%                  case 3
%                      ylim(ylim_value_beta);
%             end
            ylim_eval = [round(prctile([T_plot.(measure_name)],1,'all'),1)+round(prctile([T_plot.(measure_name)],1,'all'),1)*0.8,round(prctile([T_plot.(measure_name)],99,'all')+round(prctile([T_plot.(measure_name)],99,'all'))*0.7,1)]; 
            ylim(ylim_eval)
            % axis padded
            fig_i = get(groot,'CurrentFigure');
            box off
            title(title_plot{i},'FontWeight','bold','FontSize',12);
            T_stats_plot = psd_feature_stats(psd_feature_stats.study == num2str(g) & psd_feature_stats.cluster == num2str(k),:);
            switch i 
                case 1
                    anova_stats = T_stats_plot.theta_anova;
                    cond2_stats = T_stats_plot.theta_cond2_pval;
                    cond3_stats = T_stats_plot.theta_cond3_pval;
                    cond4_stats = T_stats_plot.theta_cond4_pval;
                    regress_sig = T_stats_plot.Th_num_pval;
                    regressline_stats = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                    R2 = T_stats_plot.Th_num_R2;
                case 2
                    anova_stats = T_stats_plot.alpha_anova;
                    cond2_stats = T_stats_plot.alpha_cond2_pval;
                    cond3_stats = T_stats_plot.alpha_cond3_pval;
                    cond4_stats = T_stats_plot.alpha_cond4_pval;
                    regress_sig = T_stats_plot.A_num_pval;
                    regressline_stats = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                    R2 = T_stats_plot.A_num_R2;
                case 3
                    anova_stats = T_stats_plot.beta_anova;
                    cond2_stats = T_stats_plot.beta_cond2_pval;
                    cond3_stats = T_stats_plot.beta_cond3_pval;
                    cond4_stats = T_stats_plot.beta_cond4_pval;
                    regress_sig = T_stats_plot.B_num_pval;
                    regressline_stats = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                    R2 = T_stats_plot.B_num_R2;
            end
            if g == 1 % terrain condition, categorical values
                 if anova_stats < 0.05
                    if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
                    if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
                    if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
                 end
            elseif g == 2 % g == 2 % speed conditions, numeric values, regression line
                if anova_stats < 0.05
                    % plot line
                    x = 0:5;
                    y = x*regressline_stats(2)*0.25 + regressline_stats(1);
                    plot(x,y,'-','color','k','linewidth',1);
                    xlim([0 5]);
                    if regress_sig > 0.01 & regress_sig < 0.05
                        text(0.51,gety(gca)*1.3,{['* ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',8,'FontName',PLOT_STRUCT.font_name);
                    elseif regress_sig <= 0.01 & regress_sig > 0.001 
                        text(0.5,gety(gca)*1.3,{['** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',8,'FontName',PLOT_STRUCT.font_name);
                    else
                        text(0.5,gety(gca)*1.3,{['*** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',8,'FontName',PLOT_STRUCT.font_name);
                    end
                end
            end
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
    %         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.08+horiz_shift,0.375,0.25*im_resize,0.225*im_resize]);  %[left bottom width height]
            hold off;
            %- iterate
            horiz_shift = horiz_shift + 0.3*im_resize;
        end
        hold off;
        %%
        g = 1;
        switch g
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
        end
        %## violin plot's theta/alpha/beta (terrain)
%         horiz_shift= 0;
        for i = 1:length(measure_name_plot)
            measure_name = measure_name_plot{i};
            T_plot = psd_feature_table(psd_feature_table.study == num2str(g) & psd_feature_table.cluster == num2str(k),:);
%             subplot(3,3,6+i)
            if i>1
                subplot(3,8,21+i)
            else
                subplot(3,8,21+i)
            end
            hold on;
            violinplot(T_plot.(measure_name),T_plot.cond,...
                'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
                'MarkerSize',8,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin',...
                'full','width',0.3,'QuartileStyle','shadow');
            ylabel('');
            if i == 1
                xlabel('Terrain');
            else
                xlabel('');
            end
%             if i == 1;ylabel('10*log_{10}(Flattened PSD)');else;ylabel('');end
            if g == 2&&i==1;xlabel('speed (m/s)');end
            ax = gca;
            set(ax,'xticklabel', xtick_label_g,'XTickLabelRotation',45,'FontWeight','bold');
%             switch i 
%                  case 1
%                      ylim(ylim_value_theta);
%                  case 2
%                      ylim(ylim_value_alpha);
%                  case 3
%                      ylim(ylim_value_beta);
%             end
            ylim_eval = [round(prctile([T_plot.(measure_name)],1,'all'),1)+round(prctile([T_plot.(measure_name)],1,'all'),1)*0.8,round(prctile([T_plot.(measure_name)],99,'all')+round(prctile([T_plot.(measure_name)],99,'all'))*0.7,1)]; 
            ylim(ylim_eval)
            % axis padded
            fig_i = get(groot,'CurrentFigure');
            box off
            title(title_plot{i},'FontWeight','bold','FontSize',12);
            T_stats_plot = psd_feature_stats(psd_feature_stats.study == num2str(g) & psd_feature_stats.cluster == num2str(k),:);
            switch i 
                case 1
                    anova_stats = T_stats_plot.theta_anova;
                    cond2_stats = T_stats_plot.theta_cond2_pval;
                    cond3_stats = T_stats_plot.theta_cond3_pval;
                    cond4_stats = T_stats_plot.theta_cond4_pval;
                    regress_sig = T_stats_plot.Th_num_pval;
                    regressline_stats = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                    R2 = T_stats_plot.Th_num_R2;
                case 2
                    anova_stats = T_stats_plot.alpha_anova;
                    cond2_stats = T_stats_plot.alpha_cond2_pval;
                    cond3_stats = T_stats_plot.alpha_cond3_pval;
                    cond4_stats = T_stats_plot.alpha_cond4_pval;
                    regress_sig = T_stats_plot.A_num_pval;
                    regressline_stats = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                    R2 = T_stats_plot.A_num_R2;
                case 3
                    anova_stats = T_stats_plot.beta_anova;
                    cond2_stats = T_stats_plot.beta_cond2_pval;
                    cond3_stats = T_stats_plot.beta_cond3_pval;
                    cond4_stats = T_stats_plot.beta_cond4_pval;
                    regress_sig = T_stats_plot.B_num_pval;
                    regressline_stats = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                    R2 = T_stats_plot.B_num_R2;
            end
            if g == 1 % terrain condition, categorical values
                 if anova_stats < 0.05
                    if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
                    if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
                    if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
                 end
            elseif g == 2 % g == 2 % speed conditions, numeric values, regression line
                if anova_stats < 0.05
                    % plot line
                    x = 0:5;
                    y = x*regressline_stats(2)*0.25 + regressline_stats(1);
                    plot(x,y,'-','color','k','linewidth',1);
                    xlim([0 5]);
                    if regress_sig > 0.01 & regress_sig < 0.05
                        text(0.5,gety(gca)*1.2,{['* ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',7,'FontName',PLOT_STRUCT.font_name);
                    elseif regress_sig <= 0.01 & regress_sig > 0.001 
                        text(0.5,gety(gca)*1.2,{['** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',7,'FontName',PLOT_STRUCT.font_name);
                    else
                        text(0.5,gety(gca)*1.2,{['*** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]},'FontSize',7,'FontName',PLOT_STRUCT.font_name);
                    end
                end
            end
            set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold')
    %         set(ax,'OuterPosition',[0.1,0.35,0.15,0.3]);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[0.1+horiz_shift,0.375,0.25*im_resize,0.225*im_resize]);  %[left bottom width height]
            hold off;
            %- iterate
            horiz_shift = horiz_shift + 0.3*im_resize;
        end
        hold off;
        pause(1)
        %% ================================================================= %%
%         exportgraphics(gcf,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_des%i.pdf',k,g)],'ContentType','vector','Resolution',300)
%         exportgraphics(gcf,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_des%i.jpg',k,g)],'Resolution',300)
        exportgraphics(fig,[save_dir filesep sprintf('TopoPSD_Violins_cl%i_combo.tiff',k)],'Resolution',1000)
        close(fig);
%     end
    
    
end
