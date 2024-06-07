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
SPEED_CUTOFF = 0.1;
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
% study_dir_name = '04162024_MIM_OAN57_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- study info
SUB_GROUP_FNAME = 'all_spec';
% SUB_GROUP_FNAME = 'group_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K =13;
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
if ~exist(spec_data_dir,'dir')
    error('spec_data dir does not exist');
end
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',spec_data_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',spec_data_dir);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
%- pull out pertinent study info
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
% subject_chars = {STUDY.datasetinfo.subject};
% fPaths = {STUDY.datasetinfo.filepath};
% fNames = {STUDY.datasetinfo.filename};
try
    GROUP_CHARS = STUDY.design(1).variable(1).value;
catch
    GROUP_CHARS = [];
    fprintf('No grouping variable...\n\n');
end
CLUSTER_PICKS = main_cl_inds(2:end);
DESIGN_INDS = 1:length(STUDY.design);
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
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
% pyversion('C:\Users\jsalminen\AppData\Local\Programs\Python\Python37\python.EXE');
% pyversion('C:\Users\jsalminen\AppData\Local\Microsoft\WindowsApps\python37.EXE');
%pyversion('/cygdrive/c/Users/jsalminen/AppData/Local/Microsoft/WindowsApps/python3.8')
NOTE: windows appstore installs do not work!!
% pyversion('C:\Users\jsalminen\AppData\Local\Programs\Python\Python311\python.exe')
pyversion('c:\User\python.exe')
pe = pyenv;
%}
%%
%## RE-POP PARAMS
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
    'groupstats',ERSP_STAT_PARAMS.groupstats,...
    'method',ERSP_STAT_PARAMS.method,...
    'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
    'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
    'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
    'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
%% ================================================================== %%
% tt = readtable('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\_save\MIM_Redcap_Blood_Sppb_Grip.xlsx',...
%     'DatetimeType','text')
% tt = readtable('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\_save\MindInMotion-JacobReportEEG_DATA_LABELS_2024-03-28_1136.csv');
% delinds = zeros(size(tt,1),1);
% for i = 1:size(tt,1)
%     chk = ~isempty(tt.x12_TimeToWalk400Meters_RecordTimeThatFirstFootCrossesTheFinish{i}) ||...
%         ~isnan(tt.TotalSPPBScore_addScoresFromBalance_Gait_AndChairStand__(i)) ||...
%         ~isempty(tt.x8_BloodDraw{i}) ||...
%         ~isempty(tt.x1_HasAnyPainOrArthritisInYourHandsGottenMuchWorseRecently_{i});
%     if chk
%     else
%         delinds(i) = 1;
%     end
% end
% tt(logical(delinds),:) = [];
% new_400m = zeros(size(tt,1),1);
% for i = 1:size(tt,1)
%     tmp = datetime(tt.x12_TimeToWalk400Meters_RecordTimeThatFirstFootCrossesTheFinish(i),'Inputformat',"mm:ss");
%     tmp = minute(tmp)*60+second(tmp);
%     new_400m(i) = tmp;
% end
% tt.x400m_seconds = new_400m;
% writetable(tt,'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\MIM_Redcap_Blood_Sppb_Grip.xlsx');
%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
speed_alleeg = cell(length(STUDY.datasetinfo),2);
for i = 1:length(STUDY.datasetinfo)
    ss = STUDY.datasetinfo(i).subject;
    ind = speed_table.Var1==ss;
    if any(ind)
        speed_alleeg{i,1} = speed_table.Var1(ind);
        speed_alleeg{i,2} = double(speed_table.Var2(ind));
    else 
        fprintf('Can''t find subject %s''s speed\n',ss)
    end
end
speed_alleeg = speed_alleeg(~cellfun(@isempty,speed_alleeg(:,1)),:);
%% (TABLE) GENERATE FOOOF VALUES ======================================= %%
table_len = 0;
c_chars = nan();
g_chars = nan();
for des_i = DESIGN_INDS
    for cl_i = CLUSTER_PICKS
        s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
        for i = 1:length(STUDY.design(des_i).variable)
            if strcmp(STUDY.design(des_i).variable(i).label,'cond')
                c_chars = STUDY.design(des_i).variable(i).value;
            elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
                g_chars = STUDY.design(des_i).variable(i).value;
            end
        end
        % try
        %     g_chars = STUDY.design(des_i).variable(1).value;
        % catch
        %     g_chars = 'group';
        % end
        % c_chars = STUDY.design(des_i).variable(2).value;
        chk = true;
        try
            if isnan(g_chars)
                chk = true;
            else
                fprintf('g_chars has some weird value in it\n')
            end
        catch
            for group_i = 1:length(g_chars)
                g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
                chk = chk && sum(g_inds) > length(STUDY.design(des_i).variable(2).value);
            end
        end
        if chk
            for group_i = 1:length(g_chars)
                try
                    if isnan(g_chars)
                        g_inds = (1:length(STUDY.cluster(cl_i).sets));
                    else
                        fprintf('g_chars has some weird value in it\n')
                    end
                catch
                    g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});                   
                end
                cl_chars = s_chars(g_inds);
                for cond_i = 1:length(c_chars)
                    for subj_i = 1:length(cl_chars)
                        table_len = table_len + 1;
                    end
                end
            end
        end
    end
end
speed_ms = zeros(table_len,1);
subj_id = zeros(table_len,1);
subj_cl_ind = zeros(table_len,1);
subj_char = categorical(repmat({''},table_len,1));
comp_id = zeros(table_len,1);
design_id = categorical(zeros(table_len,1));
cluster_id = categorical(zeros(table_len,1));
group_id = zeros(table_len,1);
group_char = categorical(repmat({''},table_len,1));
% speed_double = zeros(table_len,1);
cond_id = zeros(table_len,1);
cond_char = categorical(repmat({''},table_len,1));
% speed_cat = categorical(repmat({''},table_len,1));
% terrain_cat = categorical(repmat({''},table_len,1));
aperiodic_exp = zeros(table_len,1);
aperiodic_offset = zeros(table_len,1);
central_freq = cell(table_len,1);
power = cell(table_len,1);
r_squared = zeros(table_len,1);
theta_avg_power = zeros(table_len,1);
alpha_avg_power = zeros(table_len,1);
beta_avg_power = zeros(table_len,1);
theta = cell(table_len,1);
alpha = cell(table_len,1);
beta = cell(table_len,1);
theta_center = cell(table_len,1);
alpha_center = cell(table_len,1);
beta_center = cell(table_len,1);
FOOOF_TABLE = table(speed_ms,subj_id,subj_cl_ind,subj_char,comp_id,design_id,cond_id,...
    cond_char,group_id,cluster_id,aperiodic_exp,aperiodic_offset,...
    central_freq,power,r_squared,theta_avg_power,alpha_avg_power,beta_avg_power,...
    theta,alpha,beta,theta_center,alpha_center,beta_center);
%%
% fooof_group_results_org = cell(1,length(DESIGN_INDS));
fooof_results = cell(length(DESIGN_INDS),1);
fooof_diff_store = cell(length(DESIGN_INDS),1);
fooof_apfit_store = cell(length(DESIGN_INDS),1);
spec_data_original = cell(length(DESIGN_INDS),1);
fooof_frequencies = zeros(2,1);
subj_chk = {};
cnt = 1;
for dd = 1:length(DESIGN_INDS)
    des_i = DESIGN_INDS(dd);
    for i = 1:length(CLUSTER_PICKS)
        cl_i = CLUSTER_PICKS(i);
        %## load 
        ind_cl = [STUDY.etc.mim_gen_ersp_data.clust_ind_cl] == cl_i;
        ind_des = [STUDY.etc.mim_gen_ersp_data.des_ind] == des_i;
        ind = ind_cl & ind_des;
        file_mat = STUDY.etc.mim_gen_ersp_data(ind).spec_ss_fpaths;
        tmp = par_load(file_mat,[]);     
        
        %## RUN FOOOF
        %- get subjects in cluster
        s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
        for i = 1:length(STUDY.design(des_i).variable)
            if strcmp(STUDY.design(des_i).variable(i).label,'cond')
                c_chars = STUDY.design(des_i).variable(i).value;
            elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
                g_chars = STUDY.design(des_i).variable(i).value;
            end
        end
        try
            if isnan(g_chars)
                specdata = cell(size(tmp,1),1);
                specfreqs = cell(size(tmp,1),1);
                for k = 1:size(tmp,1) % cond
                    specdata{k,1} = tmp(k,1).specdata;
                    specfreqs{k,1} = tmp(k,1).specfreqs;
                end
            else
                fprintf('g_chars has some weird value in it\n')
            end
        catch
            specdata = cell(size(tmp,2),size(tmp,1));
            specfreqs = cell(size(tmp,2),size(tmp,1));
            for j = 1:size(tmp,1) % group
                for k = 1:size(tmp,2) % cond
                    specdata{k,j} = tmp(j,k).specdata;
                    specfreqs{k,j} = tmp(j,k).specfreqs;
                end
            end
        end
        %-
        iif = (specfreqs{cond_i,group_i} < f_range(2) & specfreqs{cond_i,group_i} > f_range(1));
        fooof_frequencies = specfreqs{cond_i,group_i}(iif);
        if ~any(f_range(2) == fooof_frequencies)
            fooof_frequencies = [fooof_frequencies; f_range(2)];
        end
        if ~any(f_range(1) == fooof_frequencies)
            fooof_frequencies = [f_range(1); fooof_frequencies];
        end
        %-
        try
            if isnan(g_chars)
                chk = true;
            else
                fprintf('g_chars has some weird value in it\n')
            end
        catch
            for group_i = 1:length(g_chars)
                g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
                chk = chk && sum(g_inds) > length(STUDY.design(des_i).variable(2).value);
            end
        end
        if chk
            for group_i = 1:size(specdata,2) % in case there is young and old adult group
                try
                    if isnan(g_chars)
                        g_inds = (1:length(STUDY.cluster(cl_i).sets));
                    else
                        fprintf('g_chars has some weird value in it\n')
                    end
                catch
                    g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});                   
                end
                %-
                cl_chars = s_chars(g_inds);
                subj_inds = STUDY.cluster(cl_i).sets(g_inds);
                cl_inds = find(g_inds);
                cl_comps = STUDY.cluster(cl_i).comps(g_inds);
                cl_speeds = zeros(length(cl_chars),1);
                fprintf('%i) subjects (N=%i): %s\n',cl_i,length(cl_chars),sprintf('%s,',cl_chars{:}));
                subj_chk = [subj_chk, cl_chars];
                for j = 1:length(cl_speeds)
                    ind = cellfun(@(x) x == categorical(cl_chars(j)),speed_alleeg(:,1));
                    cl_speeds(j) = speed_alleeg{ind,2};
                end
                for cond_i = 1:size(specdata,1) % different level of terrains/speeds
                    specdata_nolog = 10.^(specdata{cond_i,group_i}/10);
                    %## MAIN FOOOF LOOP
                    return_model = true;
                    for subj_i = 1:size(specdata{cond_i,group_i},2)
                        %- run fooof
                        fooof_results{des_i}{cond_i,group_i}{subj_i} = fooof(specfreqs{cond_i,group_i}, specdata_nolog(:,subj_i), f_range, settings, return_model);
                        %- store
                        FOOOF_TABLE.speed_ms(cnt) = cl_speeds(subj_i);
                        FOOOF_TABLE.subj_id(cnt) = subj_inds(subj_i);
                        FOOOF_TABLE.subj_cl_ind(cnt) = cl_inds(subj_i);
                        FOOOF_TABLE.subj_char(cnt) = categorical(cl_chars(subj_i));
                        FOOOF_TABLE.comp_id(cnt) = cl_comps(subj_i);
                        FOOOF_TABLE.design_id(cnt) = categorical(des_i);
                        % FOOOF_TABLE.cond_id(cnt) = cond_i;
                        if any(strcmp(c_chars(cond_i),TERRAIN_VALS))
                            FOOOF_TABLE.cond_char(cnt) = categorical(c_chars(cond_i));
                            FOOOF_TABLE.cond_id(cnt) = cond_i;
                            % FOOOF_TABLE.terrain_cat(cnt) = categorical(c_chars(cond_i));
                            % FOOOF_TABLE.speed_double(cnt) = [];
                            % FOOOF_TABLE.speed_cat(cnt) = categorical([]);
                        else
                            ind = strcmp(c_chars(cond_i),SPEED_VALS(2,:));
                            FOOOF_TABLE.cond_char(cnt) = categorical(SPEED_VALS(1,ind));
                            FOOOF_TABLE.cond_id(cnt) = cond_i;
                            % ind = strcmp(c_chars(cond_i),SPEED_VALS(2,:));
                            % FOOOF_TABLE.terrain_cat(cnt) = categorical([]);
                            % FOOOF_TABLE.speed_double(cnt) = SPEED_VALS(1,ind);
                            % FOOOF_TABLE.speed_cat(cnt) = categorical(SPEED_VALS(1,ind));
                        end
                        FOOOF_TABLE.group_id(cnt) = group_i;
                        try
                            if isnan(g_chars)
                                FOOOF_TABLE.group_char(cnt) = categorical({'all'});
                            else
                                fprintf('g_chars has some weird value in it\n')
                            end
                        catch
                            FOOOF_TABLE.group_char(cnt) = categorical(g_chars(group_i));
                        end
                        FOOOF_TABLE.cluster_id(cnt) = categorical(cl_i);
                        FOOOF_TABLE.aperiodic_exp(cnt) = fooof_results{des_i}{cond_i,group_i}{subj_i}.aperiodic_params(2);
                        FOOOF_TABLE.aperiodic_offset(cnt) = fooof_results{des_i}{cond_i,group_i}{subj_i}.aperiodic_params(1);
                        FOOOF_TABLE.central_freq{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(:,1);
                        FOOOF_TABLE.power{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(:,2);
                        % try
                        %     FOOOF_TABLE.central_freq{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(1,:);
                        %     FOOOF_TABLE.power{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(2,:);
                        % catch
                        %     fprintf('No peaks generated for subject: %s...\n',cl_chars{subj_i});
                        % end
                        FOOOF_TABLE.r_squared(cnt) = fooof_results{des_i}{cond_i,group_i}{subj_i}.r_squared;
                        %- Compute average power after flatten curve
                        fooof_diff = 10*(fooof_results{des_i}{cond_i,group_i}{subj_i}.power_spectrum) - 10*(fooof_results{des_i}{cond_i,group_i}{subj_i}.ap_fit);
                        fooof_freq = fooof_results{des_i}{cond_i,group_i}{subj_i}.freqs;
                        %- store
                        FOOOF_TABLE.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                        FOOOF_TABLE.alpha_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                        FOOOF_TABLE.beta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));
                        %- data structure needs to be freq x subject
                        fooof_diff_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = fooof_diff';
                        fooof_apfit_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = 10*(fooof_results{des_i}{cond_i,group_i}{subj_i}.ap_fit);
                        %- store original spec data
                        spec_data_original{des_i}{cl_i}{cond_i,group_i} = specdata{cond_i,group_i}(specfreqs{cond_i,group_i} >= f_range(1) & specfreqs{cond_i,group_i} <= f_range(2),:);
                        %- iterate
                        cnt = cnt + 1;
                    end
                    % i_ind = i_ind + size(specdata{cond_i,group_i},2);
                end
            end
        end
    end
end
disp(length(unique(subj_chk)));
disp(length(unique(FOOOF_TABLE.subj_id)));
save([save_dir filesep 'fooof_results.mat'],'fooof_results');
%%
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
conditions = unique(FOOOF_TABLE.cond_char);
groups = unique(FOOOF_TABLE.group_char);
%-
des_i = 1;
cond_i = 1;
group_i = 1;
IM_RESIZE = 0.5;
HZ_DIM = 4;
VERTICAL_SHIFT =  0.2;
HORIZONTAL_SHIFT = 0.2;
HORIZ_START = 0.08;
VERTICAL_START = 0.75;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%##
for j = 1:length(designs)
    des_i = double(string(designs(j)));
    %-
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
    for cond_i = 1:length(xtick_label_g)
        for group_i = 1:length(groups)
            fig = figure('color','white','renderer','Painters');
            sgtitle(sprintf('Design %i, Condition %s, Group %s',des_i,xtick_label_g{cond_i},groups(group_i)),'FontName','Arial','FontSize',14,'FontWeight','bold','Interpreter','none');
            set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
            set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
            hold on;
            set(gca,AXES_DEFAULT_PROPS{:})
            
            %-
            horiz_shift = HORIZ_START;
            hz = 0;
            vert_shift = 0;
            for i = 1:length(clusters)
                hz = hz + 1;
                cl_i = double(string(clusters(i)));
                axes();
                hold on;
                %-
                fooof_psd = fooof_diff_store{des_i}{cl_i}{cond_i,group_i}';
                fooof_psd_mean = mean(fooof_psd);
                subjs = plot(fooof_freq,fooof_psd,'color',[0,0,0,0.15],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
                mean_plot = plot(fooof_freq,fooof_psd_mean,'color',color_dark(cond_i,:),'linestyle','-','linewidth',4,'displayname','orig. subj psd');
                
                %-
                ax = gca;
                plot([0 40],[0 0],'--','color','black');
                xlim([4 40]);
                ylim([-2 10]);
                xlabel('Frequency(Hz)');
                if hz ~= 1
                    ylabel('');
                else
                    ylabel('10*log_{10}(Power)');
                end
                set(ax,'FontName','Arial',...
                    'FontSize',12,...
                    'FontWeight','bold')
                xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
                set(ax,'FontName','Arial','FontSize',10,...
                    'FontWeight','bold')
                title(sprintf('CL%i',cl_i))
                set(ax,'OuterPosition',[0,0,1,1]);
                set(ax,'Position',[horiz_shift,VERTICAL_START-vert_shift,0.3*IM_RESIZE,0.25*IM_RESIZE]);  %[left bottom width height]
                if hz < HZ_DIM
                    horiz_shift = horiz_shift + HORIZONTAL_SHIFT;
                    
                else
                    vert_shift = vert_shift + VERTICAL_SHIFT;
                    horiz_shift = HORIZ_START;
                    hz = 0;
                end
            end
            hold off;
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_cond%i_group%i_cluster_fooofs_subjs_means.tiff',cl_i,des_i,cond_i,group_i)],'Resolution',600)
            close(fig);
        end
    end
end
%% (JS) DETERMINE PEAK FREQUENCIES
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        inds = find(FOOOF_TABLE.design_id == designs(i) & FOOOF_TABLE.cluster_id == clusters(j));
        for k = 1:length(inds)
            subj_i = FOOOF_TABLE.subj_id(inds(k));
            cond_i = FOOOF_TABLE.cond_id(inds(k));
            group_i = FOOOF_TABLE.group_id(inds(k));
            % tmp = fooof_diff_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i);
            % [val,ind] = findpeaks(tmp);
            % fooof_frequencies(ind)
            %-
            % figure;
            % hold on;
            % plot(fooof_frequencies,tmp)
            % scatter(fooof_frequencies(ind),val,'k^');
            % scatter(FOOOF_TABLE.central_freq{inds(k)},FOOOF_TABLE.power{inds(k)},'go')
            % hold off;
            %-
            if ~isempty(FOOOF_TABLE.central_freq{inds(k)})
                for l = 1:length(FOOOF_TABLE.central_freq{inds(k)})
                    cf = FOOOF_TABLE.central_freq{inds(k)}(l);
                    if cf > theta_band(1) & cf <= theta_band(2)
                        FOOOF_TABLE.theta{inds(k)} = [FOOOF_TABLE.theta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > alpha_band(1) & cf <= alpha_band(2)
                        FOOOF_TABLE.alpha{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > beta_band(1) & cf <= beta_band(2)
                        FOOOF_TABLE.beta{inds(k)} = [FOOOF_TABLE.beta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    end
                end
                if length(FOOOF_TABLE.theta{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.theta{inds(k)}(:,1)-6));
                    temp_power = FOOOF_TABLE.theta{inds(k)}(:,2);
                    FOOOF_TABLE.theta_center{inds(k)} = [FOOOF_TABLE.theta{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.alpha{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.alpha{inds(k)}(:,1)-10));
                    temp_power = FOOOF_TABLE.alpha{inds(k)}(:,2);
                    FOOOF_TABLE.alpha_center{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.beta{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.beta{inds(k)}(:,1)-20));
                    temp_power = FOOOF_TABLE.beta{inds(k)}(:,2);
                    FOOOF_TABLE.beta_center{inds(k)} = [FOOOF_TABLE.beta{inds(k)}(indx,1) temp_power(indx)];
                end
            end
        end
    end
end
%## SAVE DATA
% save([save_dir filesep 'fooof_results_summary.mat'],'fooof_group_results_org');
save([save_dir filesep 'fooof_diff_store.mat'],'fooof_diff_store');
save([save_dir filesep 'fooof_apfit_store.mat'],'fooof_apfit_store');
save([save_dir filesep 'spec_data_original.mat'],'spec_data_original');
%% DETERMINE PEAK FREQUENCIES
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for i = 1:length(designs)
    des_i = designs(i);
    for j = 1:length(clusters)
        cl_i = clusters(j);
        inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == cl_i);
        for k = 1:length(inds)
            if ~isempty(FOOOF_TABLE.central_freq{inds(k)})
                for l = 1:length(FOOOF_TABLE.central_freq{inds(k)})
                    cf = FOOOF_TABLE.central_freq{inds(k)}(l);
                    if cf > theta_band(1) & cf <= theta_band(2)
                        FOOOF_TABLE.theta{inds(k)} = [FOOOF_TABLE.theta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > alpha_band(1) & cf <= alpha_band(2)
                        FOOOF_TABLE.alpha{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > beta_band(1) & cf <= beta_band(2)
                        FOOOF_TABLE.beta{inds(k)} = [FOOOF_TABLE.beta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    end
                end
                if length(FOOOF_TABLE.theta{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.theta{inds(k)}(:,1)-6));
                    temp_power = FOOOF_TABLE.theta{inds(k)}(:,2);
                    FOOOF_TABLE.theta_center{inds(k)} = [FOOOF_TABLE.theta{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.alpha{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.alpha{inds(k)}(:,1)-10));
                    temp_power = FOOOF_TABLE.alpha{inds(k)}(:,2);
                    FOOOF_TABLE.alpha_center{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.beta{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.beta{inds(k)}(:,1)-20));
                    temp_power = FOOOF_TABLE.beta{inds(k)}(:,2);
                    FOOOF_TABLE.beta_center{inds(k)} = [FOOOF_TABLE.beta{inds(k)}(indx,1) temp_power(indx)];
                end
            end
        end
    end
end
%## SAVE DATA
% save([save_dir filesep 'fooof_results_summary.mat'],'fooof_group_results_org');
save([save_dir filesep 'fooof_diff_store.mat'],'fooof_diff_store');
save([save_dir filesep 'fooof_apfit_store.mat'],'fooof_apfit_store');
save([save_dir filesep 'spec_data_original.mat'],'spec_data_original');
%%
% tmp_study = STUDY;
% RE_CALC = true;
% if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') || isfield(tmp_study.cluster,'topopol') 
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
% end
% if ~isfield(tmp_study.cluster,'topo'), tmp_study.cluster(1).topo = [];end
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% for i = 1:length(designs)
%     des_i = string(designs(i));
%     for j = 1:length(clusters) % For each cluster requested
%         cl_i = double(string(clusters(j)));
%         if isempty(tmp_study.cluster(cl_i).topo) || RE_CALC
%             inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == string(cl_i));
%             sets_i = unique([FOOOF_TABLE.subj_cl_ind(inds)]);
%             tmp_study.cluster(cl_i).sets = tmp_study.cluster(cl_i).sets(sets_i);
%             tmp_study.cluster(cl_i).comps = tmp_study.cluster(cl_i).comps(sets_i);
%             tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
%             STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
%             STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
%             STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
%             STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
%             STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
%         end
%     end
% end
    
%% Create table from group results, take mean across participants ICs    
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
for i = 1:length(designs)
    des_i = designs(i);
    for j = 1:length(clusters)
        cl_i = clusters(j);
        inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == cl_i);
        for k = 1:length(inds)
            if isempty(FOOOF_TABLE.alpha{inds(k)})
                FOOOF_TABLE.alpha{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta{inds(k)})
                FOOOF_TABLE.beta{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.alpha_center{inds(k)})
                FOOOF_TABLE.alpha_center{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta_center{inds(k)})
                FOOOF_TABLE.beta_center{inds(k)} = [NaN, NaN];
            end
            [~,idx_a] = max(FOOOF_TABLE.alpha{inds(k)}(:,2));
            [~,idx_b] = max(FOOOF_TABLE.beta{inds(k)}(:,2));
        end
    end
end
FOOOF_TABLE.subj_id = categorical(FOOOF_TABLE.subj_id);
FOOOF_TABLE.cond_id = categorical(FOOOF_TABLE.cond_id);
FOOOF_TABLE.cluster_id = categorical(FOOOF_TABLE.cluster_id);
FOOOF_TABLE.design_id = categorical(FOOOF_TABLE.design_id);
FOOOF_TABLE.group_id = categorical(FOOOF_TABLE.group_id);
writetable(FOOOF_TABLE,[save_dir filesep 'fooof_spec_table.xlsx']);
save([save_dir filesep 'psd_feature_table.mat'],'FOOOF_TABLE');
%% (STATISTICS CALCULATIONS) =========================================== %% 
%## Preliminary stats on average power (substracted background)
% Create stats table
TERRAIN_DES_ID = 1;
SPEED_DES_ID = 2;
psd_feature_stats = table;
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
for i = 1:length(designs)
    des_i = designs(i);
    for j = 1:length(clusters)
        cl_i = clusters(j);
        for k = 1:length(groups)
            group_i = groups(k);
            %- extract table for cluster and design
            inds = (FOOOF_TABLE.cluster_id == cl_i & FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.group_id == group_i);
            t1 = FOOOF_TABLE(inds,:);
            %- re-categorize variables for stats
            t1.cond = categorical(t1.cond_id);
            % t1.cond = double(t1.cond_id)*0.25;
            t1.log_alpha_avg_power = log(t1.alpha_avg_power+5);
            t1.log_theta_avg_power = log(t1.theta_avg_power+5);
            %## theta
            %- test lienar model
            mod = 'theta_avg_power ~ 1 + cond + (1|subj_id)';
            lme_theta_avg_power = fitlme(t1,mod,'FitMethod','ML');
    %         lme_theta_avg_power = fitlme(t1,'theta_avg_power ~ cond + (1|speed_ms)');
            %-
            theta_anova = anova(lme_theta_avg_power);
            %- test normality
            [theta_h,theta_p] = lillietest(lme_theta_avg_power.residuals);
            % if theta_h
            %     histogram(t1.theta_avg_power)
            %     % fitglme(t1,mod,'Distribution','Poisson');
            %     newOut = fitglme(t1,mod,'Distribution','Gamma');
            %     waitforbuttonpress;
            % end
            %## alpha
            %- test linear model
            lme_alpha_avg_power = fitlme(t1,'alpha_avg_power ~ cond + (1|subj_id)');
    %         lme_alpha_avg_power = fitlme(t1,'alpha_avg_power ~ cond + (1|speed_ms)');
            %-
            alpha_anova = anova(lme_alpha_avg_power); 
            %- test normality
            [alpha_h,alpha_p] = lillietest(lme_alpha_avg_power.residuals);
            %## beta
            lme_beta_avg_power = fitlme(t1,'beta_avg_power ~ cond + (1|subj_id)');
    %         lme_beta_avg_power = fitlme(t1,'beta_avg_power ~ cond + (1|speed_ms)');
            %-
            beta_anova = anova(lme_beta_avg_power); 
            %- test normality
            [beta_h,beta_p] = lillietest(lme_beta_avg_power.residuals);
            %## LOG stats
    %         lme_log_theta_avg_power = fitlme(t1,'log_theta_avg_power ~ cond + (1|subID)');
    %         log_Th = anova(lme_log_theta_avg_power);[h,p] = lillietest(lme_log_theta_avg_power.residuals)
    %         lme_log_alpha_avg_power= fitlme(t1,'log_alpha_avg_power ~ cond + (1|subID)');
    %         log_A = anova(lme_log_alpha_avg_power);[h,p] = lillietest(lme_log_alpha_avg_power.residuals)
            %## add stats for aperiod fit/offsets
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
            temp_stats.study = des_i;
            temp_stats.cluster = cl_i;
            temp_stats.group = group_i;
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
            
           
            if des_i == num2str(SPEED_DES_ID)
                % use continuous variable
                % t1.cond = double(string(t1.cond));
                t1.cond = double(string(t1.cond_char));
                %- theta
    %             lme_theta_avg_power_num = fitlme(t1,'theta_avg_power ~ cond + (1|speed_ms)');
                lme_theta_avg_power_num = fitlme(t1,'theta_avg_power ~ cond + (1|subj_id)');
                Th_num = anova(lme_theta_avg_power_num);[theta_h,theta_p] = lillietest(lme_theta_avg_power_num.residuals);
                %- alpha
    %             lme_alpha_avg_power_num = fitlme(t1,'alpha_avg_power ~ cond + (1|speed_ms)');
                lme_alpha_avg_power_num = fitlme(t1,'alpha_avg_power ~ cond  + (1|subj_id)');
                A_num  = anova(lme_alpha_avg_power_num);[alpha_h,alpha_p] = lillietest(lme_alpha_avg_power_num.residuals);
                %- beta
    %             lme_beta_avg_power_num = fitlme(t1,'beta_avg_power ~ cond + (1|speed_ms)');
                lme_beta_avg_power_num = fitlme(t1,'beta_avg_power ~ cond + (1|subj_id)');
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
end
psd_feature_stats.study = categorical(psd_feature_stats.study);
psd_feature_stats.cluster = categorical(psd_feature_stats.cluster);
writetable(psd_feature_stats,[save_dir filesep 'psd_band_power_stats.xlsx']);
save([save_dir filesep 'psd_band_power_stats.mat'],'psd_feature_stats');
% %% VISUALIZE 
% CLUSTER_I = 3;
% DESIGN_I = 1;
% GROUP_I = 2;
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% groups = unique(FOOOF_TABLE.group_id);
% inds = psd_feature_stats.study == num2str(des_i) & psd_feature_stats.cluster == num2str(cl_i) & psd_feature_stats.group==groups(k);
% T_stats_plot = psd_feature_stats(inds,:);


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
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%- fooof_diff_store needs to be freq x subject, and condition by row
design_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
%-
pcond = cell(length(designs),1);
pgroup = cell(length(designs),1);
pinter = cell(length(designs),1);
statcond = cell(length(designs),1);
statgroup = cell(length(designs),1);
statinter = cell(length(designs),1);
cnt = 1;

for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(fooof_diff_store{des_i}{cl_i}, stats);
        %- 
        clust_t{cnt} = cl_i;
        design_t{cnt} = des_i;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter;
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond{des_i}{cl_i} = temp_pcond;
        pgroup{des_i}{cl_i} = temp_pcond;
        pinter{des_i}{cl_i} = temp_pinter;
        statcond{des_i}{cl_i} = temp_statcond;
        statgroup{des_i}{cl_i} = temp_statgroup;
        statinter{des_i}{cl_i} = temp_statinter;
        %%%%%
        for k0 = 1:length(pcond{des_i}{cl_i})
            pcond{des_i}{cl_i}{k0}(:,2) = pcond{des_i}{cl_i}{k0}(:,1)<0.05;    
            pcond_test_t{cnt} = pcond{des_i}{cl_i}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup{des_i}{cl_i})
            if ~isempty(pgroup{des_i}{cl_i}{k0})
                pgroup{des_i}{cl_i}{k0}(:,2) = pgroup{des_i}{cl_i}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup{des_i}{cl_i}{k0}(:,1)<0.05;    
            end
        end
        for k0 = 1:length(pinter{des_i}{cl_i})
            if ~isempty(pinter{des_i}{cl_i}{k0})
                pinter{des_i}{cl_i}{k0}(:,2) = pinter{des_i}{cl_i}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter{des_i}{cl_i}{k0}(:,1)<0.05;    
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
freq_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
% subj_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
% cond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
design_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%-
pcond_org = cell(length(designs),1);
pgroup_org = cell(length(designs),1);
pinter_org = cell(length(designs),1);
statcond_org = cell(length(designs),1);
statgroup_org = cell(length(designs),1);
statinter_org = cell(length(designs),1);
cnt = 1;
for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(spec_data_original{des_i}{cl_i}, stats);
        %- 
        clust_t{cnt} = cl_i;
        design_t{cnt} = des_i;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter; %{temp_pinter};
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond_org{des_i}{cl_i} = temp_pcond;
        pgroup_org{des_i}{cl_i} = temp_pcond;
        pinter_org{des_i}{cl_i} = temp_pinter;
        statcond_org{des_i}{cl_i} = temp_statcond;
        statgroup_org{des_i}{cl_i} = temp_statgroup;
        statinter_org{des_i}{cl_i} = temp_statinter;
        for k0 = 1:length(pcond_org{des_i}{cl_i})
            pcond_org{des_i}{cl_i}{k0}(:,2) = pcond_org{des_i}{cl_i}{k0}(:,1)<0.05; 
            pcond_test_t{cnt} = pcond_org{des_i}{cl_i}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup_org{des_i}{cl_i})
            if ~isempty(pgroup_org{des_i}{cl_i}{k0})
                pgroup_org{des_i}{cl_i}{k0}(:,2) = pgroup_org{des_i}{cl_i}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup_org{des_i}{cl_i}{k0}(:,1)<0.05;   
            end
        end
        for k0 = 1:length(pinter_org{des_i}{cl_i})
            if ~isempty(pinter_org{des_i}{cl_i}{k0})
                pinter_org{des_i}{cl_i}{k0}(:,2) = pinter_org{des_i}{cl_i}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter_org{des_i}{cl_i}{k0}(:,1)<0.05; %{pinter_org{g}{k}{k0}(:,1)<0.05};  
            end
        end            
        cnt=cnt+1;
    end
end
save([save_dir filesep 'org_pcond.mat'],'pcond_org');
%{
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
% tmp_study = STUDY;
% RE_CALC = true;
% if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') || isfield(tmp_study.cluster,'topopol') 
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
%     tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
% end
% if ~isfield(tmp_study.cluster,'topo'), tmp_study.cluster(1).topo = [];end
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% for i = 1:length(designs)
%     des_i = string(designs(i));
%     for j = 1:length(clusters) % For each cluster requested
%         cl_i = double(string(clusters(j)));
%         if isempty(tmp_study.cluster(cl_i).topo) || RE_CALC
%             inds = find(FOOOF_TABLE.design_id == des_i & FOOOF_TABLE.cluster_id == string(cl_i));
%             sets_i = unique([FOOOF_TABLE.subj_cl_ind(inds)]);
%             tmp_study.cluster(cl_i).sets = tmp_study.cluster(cl_i).sets(sets_i);
%             tmp_study.cluster(cl_i).comps = tmp_study.cluster(cl_i).comps(sets_i);
%             tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
%             STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
%             STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
%             STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
%             STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
%             STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
%         end
%     end
% end
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
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%-
for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        data_min = min([fooof_diff_store{des_i}{cl_i}{:}],[],'all');
        data_max = max([fooof_diff_store{des_i}{cl_i}{:}],[],'all');
        cnt = 1;
        figure('color','white');
        for k = 1:size(fooof_diff_store{des_i}{cl_i},1)
            for l = 1:size(fooof_diff_store{des_i}{cl_i},2)
                hold on;
                subplot(size(fooof_diff_store{des_i}{cl_i},1),size(fooof_diff_store{des_i}{cl_i},2),cnt)
                data = fooof_diff_store{des_i}{cl_i}{k,l};
                if des_i == num2str(TERRAIN_DES_ID)
                    if l == 1
                        tmpc = color.terrain(i,:);
                    else
                        tmpc = color.terrain(i,:)*0.5;
                    end
                else
                    if l == 2
                        tmpc = color.speed(i,:);
                    else
                        tmpc = color.speed(i,:)*0.5;
                    end
                end
                plot(fooof_freq,data,'color',tmpc);
                ylabel('log10(Power)')
                ylim([data_min data_max]);
                xlabel('Frequency(Hz)');
                % title(['Cluster ',num2str(cl_i)]);
                cnt = cnt + 1;
            end
        end
        hold off;
        drawnow;
    end
end
%% Sanity check: Plot distribution of aperiodic params (exp), central frequency, and goodness of fit
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%-
for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
        figure();
        hold on;
        set(gcf,'color','white');
        subplot(2,2,1)
        boxchart(temp_table.cond_id,temp_table.aperiodic_exp);
        ylabel('Aperodic exponent');
        subplot(2,2,3)
        boxchart(temp_table.cond_id,temp_table.r_squared);
        ylabel('R squared');
        hold off;
        drawnow;
    end
end
%% 
COLORS_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%-
des_i = 1;
cl_i = 3;
cond_i = 1;
group_i = 1;
%-
subj_ap_fits = fooof_apfit_store{des_i}{cl_i}{cond_i,group_i}';
ap_fit_mean = mean(subj_ap_fits);
subj_orig_psd = spec_data_original{des_i}{cl_i}{cond_i,group_i}';
orig_psd_mean = mean(subj_orig_psd);
fooof_psd = fooof_diff_store{des_i}{cl_i}{cond_i,group_i}';
fooof_psd_mean = mean(fooof_psd);

%## SUBJECTS AP FITS
% AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',[0.5,0.5,4,4])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
% set(gca,AXES_DEFAULT_PROPS{:})
hold on;
% mean_plot = plot(fooof_freq,orig_psd_mean,'color',COLOR_MAPS_SPEED(1,:),'linestyle','-','linewidth',3,'displayname','orig. psd mean');
% dash = plot(fooof_freq,ap_fit_mean,'color',COLOR_MAPS_SPEED(1,:),'linestyle','-.','linewidth',3,'displayname','ap. fit');
subjs = plot(fooof_freq,subj_orig_psd,'color',[COLOR_MAPS_SPEED(1,:),0.75],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
subjs_ap = plot(fooof_freq,subj_ap_fits,'color',[COLOR_MAPS_SPEED(2,:),0.5],'linestyle','-.','linewidth',2,'displayname','orig. subj psd');
hold off;
ax = gca;
xlim([4 40]);
ylim([-35 -5]);
% plot([0 40],[0 0],'--','color','black');
xlabel('Frequency(Hz)');
ylabel('10*log_{10}(Power)');
set(ax,'FontName','Arial',...
    'FontSize',14,...
    'FontWeight','bold')
%-
exportgraphics(fig,[save_dir filesep sprintf('ap_fit_subjs_tmpplot.tiff')],'Resolution',600)

%## MEAN AP FITS
% fig = figure;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',[0.5,0.5,4,4])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
% set(gca,AXES_DEFAULT_PROPS{:})
hold on;
mean_plot = plot(fooof_freq,orig_psd_mean,'color',COLOR_MAPS_SPEED(1,:),'linestyle','-','linewidth',3,'displayname','orig. psd mean');
dash = plot(fooof_freq,ap_fit_mean,'color',COLOR_MAPS_SPEED(2,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
% subjs = plot(fooof_freq,subj_orig_psd,'color',[COLOR_MAPS_SPEED(1,:),0.75],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
% subjs_ap = plot(fooof_freq,subj_ap_fits,'color',[COLOR_MAPS_SPEED(2,:),0.5],'linestyle','-.','linewidth',2,'displayname','orig. subj psd');
hold off;
ax = gca;
xlim([4 40]);
ylim([-35 -5]);
% plot([0 40],[0 0],'--','color','black');
xlabel('Frequency(Hz)');
ylabel('10*log_{10}(Power)');
set(ax,'FontName','Arial',...
    'FontSize',14,...
    'FontWeight','bold')
%-
exportgraphics(fig,[save_dir filesep sprintf('ap_fit_mean_tmpplot.tiff')],'Resolution',600)

%## FLATTEND SUBJS
% fig = figure;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',[0.5,0.5,4,4])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
% set(gca,AXES_DEFAULT_PROPS{:})
hold on;
colormap(linspecer);
% mean_plot = plot(fooof_freq,orig_psd_mean,'color',COLOR_MAPS_SPEED(1,:),'linestyle','-','linewidth',3,'displayname','orig. psd mean');
% dash = plot(fooof_freq,ap_fit_mean,'color',COLOR_MAPS_SPEED(2,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
subjs = plot(fooof_freq,fooof_psd,'linestyle','-','linewidth',2,'displayname','orig. subj psd');
% subjs_ap = plot(fooof_freq,subj_ap_fits,'color',[COLOR_MAPS_SPEED(2,:),0.5],'linestyle','-.','linewidth',2,'displayname','orig. subj psd');
hold off;
ax = gca;
xlim([4 40]);
ylim([-2 7.5]);
% plot([0 40],[0 0],'--','color','black');
xlabel('Frequency(Hz)');
ylabel('10*log_{10}(Power)');
set(ax,'FontName','Arial',...
    'FontSize',14,...
    'FontWeight','bold')
%-
exportgraphics(fig,[save_dir filesep sprintf('fooof_subjs_tmpplot.tiff')],'Resolution',600)

%## FLATTEND MEAN
% fig = figure;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',[0.5,0.5,4,4])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
% set(gca,AXES_DEFAULT_PROPS{:})
colormap(linspecer);
% mean_plot = plot(fooof_freq,orig_psd_mean,'color',COLOR_MAPS_SPEED(1,:),'linestyle','-','linewidth',3,'displayname','orig. psd mean');
% dash = plot(fooof_freq,ap_fit_mean,'color',COLOR_MAPS_SPEED(2,:),'linestyle','-.','linewidth',2,'displayname','ap. fit');
mean_plot = plot(fooof_freq,fooof_psd_mean,'linestyle','-','linewidth',4,'displayname','orig. subj psd');
subjs = plot(fooof_freq,fooof_psd,'color',[COLOR_MAPS_SPEED(2,:),0.2],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
% subjs_ap = plot(fooof_freq,subj_ap_fits,'color',[COLOR_MAPS_SPEED(2,:),0.5],'linestyle','-.','linewidth',2,'displayname','orig. subj psd');

ax = gca;

plot([0 40],[0 0],'--','color','black');
xlim([4 40]);
ylim([-2 7.5]);
xlabel('Frequency(Hz)');
ylabel('10*log_{10}(Power)');
set(ax,'FontName','Arial',...
    'FontSize',14,...
    'FontWeight','bold')
hold off;
%-
exportgraphics(fig,[save_dir filesep sprintf('fooof_subjs_tmpplot.tiff')],'Resolution',600)

%}