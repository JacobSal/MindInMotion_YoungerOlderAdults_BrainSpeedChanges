%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/spca_scripts/run_spca_g_psd_fooof.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
clear java;
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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%% (PARAMETERS) ======================================================== %%
%## hard define
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
%- study path
SPCA_DNAME = '01102025_mim_yaoa_spca_calcs';
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME =  '01192025_mim_yaoa_nopowpow_crit_speed';
ANALYSIS_DNAME = 'spca_fooof_psd_anl';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study';
%-
studies_fpath = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
spca_dir = [studies_fpath filesep sprintf('%s',SPCA_DNAME)];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb10'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5'];
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%-
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDY
%## SPCA STUDY
if ~ispc
    tmp = load('-mat',[spca_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spca_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    SPCA_STUDY = tmp.STUDY;
end
%## CLUSTER STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end

%## CONVERT GROUPS?
% fprintf('Assigning new groups...\n');
% groups_old = {'H1000',{'H2000','H3000'}};
% groups_new = {'ya','oa'};
% for i = 1:length(STUDY.datasetinfo)
%     chk = cellfun(@(x) contains(STUDY.datasetinfo(i).group,x),groups_old);
%     STUDY.datasetinfo(i).group = groups_new{chk};
% end

%## DESIGNS
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000','H2000','H3000'}}};
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},...
%             'variable1','group','values1',{'ya','oa'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
%             'variable1','group','values1',{'ya','oa'}}};
%% POP PARAMS
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');

%-
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%% ================================================================== %%
%## READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
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
%% LOAD IN TEST FILE & GRAB PARAMS
condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond});
subj_i = 35;
rest_psd = par_load(STUDY_GAIT.datasetinfo(subj_i).filepath,sprintf('gait_psd_spca.mat'));
freqs = rest_psd.icatimefopts.freqs;
%-
spca_table = par_load(cluster_k_dir,'spca_cluster_table_psd.mat');
%% TEST PLOT
%{
%-
subj_i = 1;
cluster_i = 3;
%-
tmp_t = spca_table(strcmp(spca_table.subj_c,STUDY.datasetinfo(subj_i).subject) & double(spca_table.cluster_n) == cluster_i,:);
figure;
hold on;
plot(freqs,tmp_t.psd_orig_avg_c{1},'DisplayName','original avg.');
% plot(freqs,tmp_t.psd_orig_baselined_c{1},'DisplayName','original avg. based');
% plot(freqs,tmp_t.psd_corr_based_c{1},'DisplayName','corrected based');
% plot(freqs,tmp_t.psd_corr_psc1{1},'DisplayName','pc1');
% plot(freqs,tmp_t.psd_corr_unbase_c{1},'DisplayName','corrected + rest');
% plot(freqs,tmp_t.psd_rest_c{1},'DisplayName','rest');

legend();
ylabel('10*log_{10}(PSD)');
xlabel('Frequency (Hz)');
% title(sprintf('%s: cluster %i',STUDY.datasetinfo(subj_i).subject,cluster_i))
title('Component 5');
hold off;
%}
%% (TABLE) GENERATE FOOOF VALUES ======================================= %%
try
    group_chars = STUDY.design(1).variable(1).value;
catch
    group_chars = [];
    fprintf('No grouping variable...\n\n');
end
cluster_inds = main_cl_inds(2:end);
design_inds = 1:length(STUDY.design);
table_len = 0;
c_chars = nan();
g_chars = nan();
for des_i = design_inds
    for cl_i = cluster_inds
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
speed_diff_stat = nan(table_len,1);
speed_div_stat = nan(table_len,1);
speed_diffdiv_stat = nan(table_len,1);
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
alpha1_avg_power = zeros(table_len,1);
beta1_avg_power = zeros(table_len,1);
alpha2_avg_power = zeros(table_len,1);
beta2_avg_power = zeros(table_len,1);
theta_band = cell(table_len,1);
alpha_band = cell(table_len,1);
beta_band = cell(table_len,1);
alpha1_band = cell(table_len,1);
beta1_band = cell(table_len,1);
alpha2_band = cell(table_len,1);
beta2_band = cell(table_len,1);
theta = cell(table_len,1);
alpha = cell(table_len,1);
beta = cell(table_len,1);
alpha1 = cell(table_len,1);
alpha2 = cell(table_len,1);
beta1 = cell(table_len,1);
beta2 = cell(table_len,1);
theta_center = cell(table_len,1);
alpha_center = cell(table_len,1);
beta_center = cell(table_len,1);
alpha1_center = cell(table_len,1);
alpha2_center = cell(table_len,1);
beta1_center = cell(table_len,1);
beta2_center = cell(table_len,1);

FOOOF_TABLE = table(speed_ms,subj_id,subj_cl_ind,subj_char,comp_id,design_id,cond_id,...
    cond_char,speed_diff_stat,speed_div_stat,speed_diffdiv_stat,group_id,cluster_id,aperiodic_exp,aperiodic_offset,...
    central_freq,power,r_squared,theta_avg_power,alpha_avg_power,beta_avg_power,...
    alpha1_avg_power,alpha2_avg_power,beta1_avg_power,beta2_avg_power,...
    theta_band,alpha_band,beta_band,...
    alpha1_band,alpha2_band,beta1_band,beta2_band,...
    alpha1,alpha2,beta1,beta2,...
    theta,alpha,beta,theta_center,alpha_center,beta_center,...
    alpha1_center,alpha2_center,beta1_center,beta2_center);
%%
% fooof_group_results_org = cell(1,length(DESIGN_INDS));
try
    group_chars = STUDY.design(1).variable(1).value;
catch
    group_chars = [];
    fprintf('No grouping variable...\n\n');
end
cluster_inds = main_cl_inds;
design_inds = 1:length(STUDY.design);
fooof_results = cell(length(design_inds),1);
fooof_diff_store = cell(length(design_inds),1);
fooof_apfit_store = cell(length(design_inds),1);
spec_data_original = cell(length(design_inds),1);
fooof_frequencies = zeros(2,1);
subj_chk = {};
cnt = 1;
ff = fopen([save_dir filesep 'subject_cluster_numbers.txt'],'w');
for dd = 1:length(design_inds)
    des_i = design_inds(dd);
    for cc = 1:length(cluster_inds)
        cl_i = cluster_inds(cc);
        %- get subjects in cluster
        s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
        for i = 1:length(STUDY.design(des_i).variable)
            if strcmp(STUDY.design(des_i).variable(i).label,'cond')
                c_chars = STUDY.design(des_i).variable(i).value;
            elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
                g_chars = STUDY.design(des_i).variable(i).value;
            end
        end
        
        %- get cluster  & condition indices
        if ~isempty(g_chars)
            specdata = cell(length(c_chars),length(g_chars));
            specdata_crop = cell(length(c_chars),length(g_chars));
        else
            specdata = cell(length(c_chars),1);
            specdata_crop = cell(length(c_chars),1);
        end
        
        for cond_i = 1:length(c_chars)
            for group_i = 1:length(g_chars)
                % inds_cl = cellfun(@(x) chk_cell(x,cl_i),spca_table.cluster_c);
                inds_cl = spca_table.cluster_n == cl_i;
                inds_cond = strcmp(spca_table.cond_c,c_chars{cond_i});
                inds_grp = spca_table.group_n == group_i;
                % inds_grp = cellfun(@(x) x==group_i,spca_table.group_n);
                inds = inds_cl & inds_cond & inds_grp;
                % spca_table.subj_c{inds}
                %- check
                g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
                % fprintf('CL%i) Condition: %s, Group: %s\n\t Cluster N=%i, Subjects Added=%i\n',...
                %     cl_i,c_chars{cond_i},g_chars{group_i},...
                %     length(STUDY.cluster(cl_i).sets),sum(inds));
                fprintf('CL%i) Condition: %s, Group: %s\n\t Cluster N=%i, Subjects Added=%i\n',...
                    cl_i,c_chars{cond_i},g_chars{group_i},...
                    sum(g_inds),sum(inds));
                tmp = cat(3,spca_table.psd_corr_unbase_c{inds});
                tmp = permute(tmp,[2,1,3]);
                specdata{cond_i,group_i} = squeeze(tmp);
                % specdata_crop{cond_i,group_i} = squeeze(mean(tmp(freq_crop,:),1));
                specfreqs{cond_i,group_i} = freqs;
            end
        end
        
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
                fprintf(ff,'%i) subjects (N=%i): %s\n',cl_i,length(cl_chars),sprintf('%s,',cl_chars{:}));
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
                            FOOOF_TABLE.speed_diff_stat(cnt) = nan();
                            % FOOOF_TABLE.terrain_cat(cnt) = categorical(c_chars(cond_i));
                            % FOOOF_TABLE.speed_double(cnt) = [];
                            % FOOOF_TABLE.speed_cat(cnt) = categorical([]);
                        else
                            ind = strcmp(c_chars(cond_i),SPEED_VALS(2,:));
                            FOOOF_TABLE.cond_char(cnt) = categorical(SPEED_VALS(1,ind));
                            FOOOF_TABLE.cond_id(cnt) = cond_i;
                            % FOOOF_TABLE.speed_diff_stat(cnt) = sqrt((double(string(SPEED_VALS(1,ind)))-cl_speeds(subj_i))^2)/double(string(SPEED_VALS{1,end}));
                            FOOOF_TABLE.speed_diff_stat(cnt) = (double(string(SPEED_VALS(1,ind)))-cl_speeds(subj_i))/double(string(SPEED_VALS{1,end}));
                            FOOOF_TABLE.speed_div_stat(cnt) = (double(string(SPEED_VALS(1,ind)))/cl_speeds(subj_i));
                            FOOOF_TABLE.speed_diffdiv_stat(cnt) = (double(string(SPEED_VALS(1,ind)))-cl_speeds(subj_i))/cl_speeds(subj_i);
                            
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
                        % FOOOF_TABLE.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                        FOOOF_TABLE.alpha_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha_band_lims(1) & fooof_freq < alpha_band_lims(2)));
                        FOOOF_TABLE.beta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta_band_lims(1) & fooof_freq < beta_band_lims(2)));
                        FOOOF_TABLE.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band_lims(1) & fooof_freq < theta_band_lims(2)));
                        FOOOF_TABLE.alpha1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha1_band_lims(1) & fooof_freq < alpha1_band_lims(2)));
                        FOOOF_TABLE.beta1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta1_band_lims(1) & fooof_freq < beta1_band_lims(2)));
                        FOOOF_TABLE.alpha2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha2_band_lims(1) & fooof_freq < alpha2_band_lims(2)));
                        FOOOF_TABLE.beta2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta2_band_lims(1) & fooof_freq < beta2_band_lims(2)));
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
        else
            fprintf('Error. Missing subjects for cluster %i, design %i',cl_i, des_i)
        end
    end
end
fclose(ff);
disp(length(unique(subj_chk)));
disp(length(unique(FOOOF_TABLE.subj_id)));
FOOOF_TABLE = FOOOF_TABLE(FOOOF_TABLE.design_id~=categorical(0),:);
par_save([save_dir filesep 'fooof_results.mat'],'fooof_results');
%##
% tmp = load([save_dir filesep 'fooof_results.mat']);
% fooof_results = tmp.fooof_results;
%% (APPEND RESTING STATE) ============================================== %%
%{
% fooof_save = FOOOF_TABLE;
tmp_fooof_table = FOOOF_TABLE;
subj_chars = unique(tmp_fooof_table.subj_char);
cluster_inds = unique(tmp_fooof_table.cluster_id);
g_chars = unique(tmp_fooof_table.group_char);
des_i = 3;
cond_i = 1;
%-
tmpd = dir([spca_dir filesep STUDY.datasetinfo(subj_i).subject filesep 'GAIT_EPOCHED' filesep '*' filesep sprintf('cond%s_spca_psd.mat',condition_gait{cond_i})]);
gait_epoch_subf = strsplit(tmpd.folder,filesep);
gait_epoch_subf = gait_epoch_subf{end};
spca_fpath = [spca_dir filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep gait_epoch_subf];
rest_psd = par_load(spca_fpath,sprintf('cond%s_spca_psd.mat',condition_gait{cond_i}));
freqs = rest_psd.freqs;
%-
cnt = size(tmp_fooof_table,1) + 1;
for cc = 1:length(cluster_inds)
    cl_i = double(string(cluster_inds(cc)));
    %- cluster_params
    study_subji = STUDY.cluster(cl_i).sets;
    s_chars = {STUDY.datasetinfo(study_subji).subject};
    g_chars = {STUDY.datasetinfo(study_subji).group};
    for ss = 1:length(s_chars)
        subj_i = study_subji(ss);
        cl_chars = s_chars(ss);
        %-
        subj_inds = STUDY.cluster(cl_i).sets(ss);
        cl_inds = ss; %find(g_inds);
        cl_comps = STUDY.cluster(cl_i).comps(ss);
        %-
        ind = cellfun(@(x) x == categorical(cl_chars),speed_alleeg(:,1));
        cl_speeds = speed_alleeg{ind,2};
        %- extract rest condition
        inds = strcmp(spca_table.subj_c,cl_chars) & spca_table.cluster_n == cl_i;
        tmpt = spca_table(inds,:);
        specdata = tmpt{1,'psd_rest_c'}; specdata = specdata{1};
        %- run fooof
        specdata_nolog = 10.^(specdata/10);
        fooof_results = fooof(freqs, specdata_nolog, f_range, settings, return_model);
        %- store
        tmp_fooof_table.speed_ms(cnt) = cl_speeds;
        tmp_fooof_table.subj_id(cnt) = subj_inds;
        tmp_fooof_table.subj_cl_ind(cnt) = cl_inds;
        tmp_fooof_table.subj_char(cnt) = categorical(cl_chars);
        tmp_fooof_table.comp_id(cnt) = cl_comps;
        tmp_fooof_table.design_id(cnt) = categorical(des_i);
        tmp_fooof_table.cond_char(cnt) = categorical({'rest'});
        tmp_fooof_table.cond_id(cnt) = cond_i;
        tmp_fooof_table.speed_diff_stat(cnt) = 0;
        tmp_fooof_table.speed_div_stat(cnt) = 0;
        tmp_fooof_table.speed_diffdiv_stat(cnt) = 0;       
        tmp_fooof_table.group_id(cnt) = group_i;
        try
            if isnan(g_chars)
                tmp_fooof_table.group_char(cnt) = categorical({'all'});
            else
                fprintf('g_chars has some weird value in it\n')
            end
        catch
            tmp_fooof_table.group_char(cnt) = categorical(g_chars(group_i));
        end
        tmp_fooof_table.cluster_id(cnt) = categorical(cl_i);
        tmp_fooof_table.aperiodic_exp(cnt) = fooof_results.aperiodic_params(2);
        tmp_fooof_table.aperiodic_offset(cnt) = fooof_results.aperiodic_params(1);
        tmp_fooof_table.central_freq{cnt} = fooof_results.peak_params(:,1);
        tmp_fooof_table.power{cnt} = fooof_results.peak_params(:,2);
        % try
        %     tmp_fooof_table.central_freq{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(1,:);
        %     tmp_fooof_table.power{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(2,:);
        % catch
        %     fprintf('No peaks generated for subject: %s...\n',cl_chars{subj_i});
        % end
        tmp_fooof_table.r_squared(cnt) = fooof_results.r_squared;
        %- Compute average power after flatten curve
        fooof_diff = 10*(fooof_results.power_spectrum) - 10*(fooof_results.ap_fit);
        fooof_freq = fooof_results.freqs;
        %- store
        % tmp_fooof_table.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
        tmp_fooof_table.alpha_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha_band_lims(1) & fooof_freq < alpha_band_lims(2)));
        tmp_fooof_table.beta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta_band_lims(1) & fooof_freq < beta_band_lims(2)));
        tmp_fooof_table.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band_lims(1) & fooof_freq < theta_band_lims(2)));
        tmp_fooof_table.alpha1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha1_band_lims(1) & fooof_freq < alpha1_band_lims(2)));
        tmp_fooof_table.beta1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta1_band_lims(1) & fooof_freq < beta1_band_lims(2)));
        tmp_fooof_table.alpha2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha2_band_lims(1) & fooof_freq < alpha2_band_lims(2)));
        tmp_fooof_table.beta2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta2_band_lims(1) & fooof_freq < beta2_band_lims(2)));
        %- data structure needs to be freq x subject
        % fooof_diff_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = fooof_diff';
        % fooof_apfit_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = 10*(fooof_results.ap_fit);
        %- store original spec data
        % spec_data_original{des_i}{cl_i}{cond_i,group_i} = specdata{cond_i,group_i}(specfreqs{cond_i,group_i} >= f_range(1) & specfreqs{cond_i,group_i} <= f_range(2),:);
        %- iterate
        cnt = cnt + 1;
    end
end
par_save([save_dir filesep 'rest_psd_feature_table.mat'],'tmp_fooof_table');
%}
%% (UPDATED FOOOF ALG.) ================================================ %%
% fooof_save = FOOOF_TABLE;
tmp_fooof_table = FOOOF_TABLE;
subj_chars = unique(spca_table.subj_c);
cluster_inds = main_cl_inds; 
des_i = 3;
cond_i = 1;
desdes = cat(1,STUDY.design.variable);
c_chars_d = desdes(strcmp({desdes.label},'cond'));
c_chars_d = {c_chars_d.value};
g_chars_d = desdes(strcmp({desdes.label},'group'));
g_chars_d = {g_chars_d.value};
cnt = size(tmp_fooof_table,1) + 1;
for cc = 1:length(cluster_inds)
    cl_i = double(string(cluster_inds(cc)));
    %- cluster_params
    study_subji = STUDY.cluster(cl_i).sets;
    s_chars = {STUDY.datasetinfo(study_subji).subject};
    g_chars = {STUDY.datasetinfo(study_subji).group};
    %-
    tmptc = spca_table(spca_table.cluster_n == cl_i,:);
    for ss = 1:length(s_chars)
        subj_i = study_subji(ss);
        cl_chars = s_chars(ss);
        %-
        subj_inds = STUDY.cluster(cl_i).sets(ss);
        cl_inds = ss; %find(g_inds);
        cl_comps = STUDY.cluster(cl_i).comps(ss);
        %-
        ind = cellfun(@(x) x == categorical(cl_chars),speed_alleeg(:,1));
        cl_speeds = speed_alleeg{ind,2};    
        %-
        tmptcs = tmptc(strcmp(tmptc.subj_c,cl_chars),:);
        %- extract condition data
        for cc_i = 1:length(c_chars)
            %- id design
            if cellfun(@(x) strcmp(c_chars{cc_i},c_chars))
            %- get condition data
            tmptcsc = tmptcs(strcmp(tmptcs.cond_c,c_chars{cc_i}),:);
            specdata = tmptcsc{1,'psd_corr_unbase_c'}; specdata = specdata{1};
            %- run fooof
            specdata_nolog = 10.^(specdata/10);
            fooof_results = fooof(freqs, specdata_nolog, f_range, settings, return_model);
            %- store
            tmp_fooof_table.speed_ms(cnt) = cl_speeds;
            tmp_fooof_table.subj_id(cnt) = subj_inds;
            tmp_fooof_table.subj_cl_ind(cnt) = cl_inds;
            tmp_fooof_table.subj_char(cnt) = categorical(cl_chars);
            tmp_fooof_table.comp_id(cnt) = cl_comps;
            tmp_fooof_table.design_id(cnt) = categorical(des_i);
            tmp_fooof_table.cond_char(cnt) = categorical({'rest'});
            tmp_fooof_table.cond_id(cnt) = cond_i;
            tmp_fooof_table.speed_diff_stat(cnt) = 0;
            tmp_fooof_table.speed_div_stat(cnt) = 0;
            tmp_fooof_table.speed_diffdiv_stat(cnt) = 0;       
            tmp_fooof_table.group_id(cnt) = group_i;
            try
                if isnan(g_chars)
                    tmp_fooof_table.group_char(cnt) = categorical({'all'});
                else
                    fprintf('g_chars has some weird value in it\n')
                end
            catch
                tmp_fooof_table.group_char(cnt) = categorical(g_chars(group_i));
            end
            tmp_fooof_table.cluster_id(cnt) = categorical(cl_i);
            tmp_fooof_table.aperiodic_exp(cnt) = fooof_results.aperiodic_params(2);
            tmp_fooof_table.aperiodic_offset(cnt) = fooof_results.aperiodic_params(1);
            tmp_fooof_table.central_freq{cnt} = fooof_results.peak_params(:,1);
            tmp_fooof_table.power{cnt} = fooof_results.peak_params(:,2);
            % try
            %     tmp_fooof_table.central_freq{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(1,:);
            %     tmp_fooof_table.power{cnt} = fooof_results{des_i}{cond_i,group_i}{subj_i}.peak_params(2,:);
            % catch
            %     fprintf('No peaks generated for subject: %s...\n',cl_chars{subj_i});
            % end
            tmp_fooof_table.r_squared(cnt) = fooof_results.r_squared;
            %- Compute average power after flatten curve
            fooof_diff = 10*(fooof_results.power_spectrum) - 10*(fooof_results.ap_fit);
            fooof_freq = fooof_results.freqs;
            %- store
            % tmp_fooof_table.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
            tmp_fooof_table.alpha_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha_band_lims(1) & fooof_freq < alpha_band_lims(2)));
            tmp_fooof_table.beta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta_band_lims(1) & fooof_freq < beta_band_lims(2)));
            tmp_fooof_table.theta_avg_power(cnt) = mean(fooof_diff(fooof_freq >= theta_band_lims(1) & fooof_freq < theta_band_lims(2)));
            tmp_fooof_table.alpha1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha1_band_lims(1) & fooof_freq < alpha1_band_lims(2)));
            tmp_fooof_table.beta1_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta1_band_lims(1) & fooof_freq < beta1_band_lims(2)));
            tmp_fooof_table.alpha2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= alpha2_band_lims(1) & fooof_freq < alpha2_band_lims(2)));
            tmp_fooof_table.beta2_avg_power(cnt) = mean(fooof_diff(fooof_freq >= beta2_band_lims(1) & fooof_freq < beta2_band_lims(2)));
            %- data structure needs to be freq x subject
            % fooof_diff_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = fooof_diff';
            % fooof_apfit_store{des_i}{cl_i}{cond_i,group_i}(:,subj_i) = 10*(fooof_results.ap_fit);
            %- store original spec data
            % spec_data_original{des_i}{cl_i}{cond_i,group_i} = specdata{cond_i,group_i}(specfreqs{cond_i,group_i} >= f_range(1) & specfreqs{cond_i,group_i} <= f_range(2),:);
            %- iterate
            cnt = cnt + 1;
        end
    end
end
par_save(tmp_fooof_table,[save_dir filesep 'new_psd_feature_table.mat']);

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
                    % if cf > theta_band(1) & cf <= theta_band(2)
                    %     FOOOF_TABLE.theta{inds(k)} = [FOOOF_TABLE.theta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    % elseif cf > alpha_band(1) & cf <= alpha_band(2)
                    %     FOOOF_TABLE.alpha{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    % elseif cf > beta_band(1) & cf <= beta_band(2)
                    %     FOOOF_TABLE.beta{inds(k)} = [FOOOF_TABLE.beta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    % end
                    if cf > theta_band_lims(1) & cf <= theta_band_lims(2)
                        FOOOF_TABLE.theta{inds(k)} = [FOOOF_TABLE.theta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > alpha_band_lims(1) & cf <= alpha_band_lims(2)
                        FOOOF_TABLE.alpha{inds(k)} = [FOOOF_TABLE.alpha{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > beta_band_lims(1) & cf <= beta_band_lims(2)
                        FOOOF_TABLE.beta{inds(k)} = [FOOOF_TABLE.beta{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > alpha1_band_lims(1) & cf <= alpha1_band_lims(2)
                        FOOOF_TABLE.alpha1_band{inds(k)} = [FOOOF_TABLE.alpha1_band{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > alpha2_band_lims(1) & cf <= alpha2_band_lims(2)
                        FOOOF_TABLE.alpha2_band{inds(k)} = [FOOOF_TABLE.alpha2_band{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > beta1_band_lims(1) & cf <= beta1_band_lims(2)
                        FOOOF_TABLE.beta1_band{inds(k)} = [FOOOF_TABLE.beta1_band{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    elseif cf > beta2_band_lims(1) & cf <= beta2_band_lims(2)
                        FOOOF_TABLE.beta2_band{inds(k)} = [FOOOF_TABLE.beta2_band{inds(k)}; cf, FOOOF_TABLE.power{inds(k)}(l)];
                    end
                    % alpha1_band = [8,10.5];
                    % alpha2_band = [10.5,13];
                    % beta1_band = [13,20];
                    % beta2_band = [20,30];
                end
                % if length(FOOOF_TABLE.theta{inds(k)}) > 1
                %     [~,indx] = min(abs(FOOOF_TABLE.theta{inds(k)}(:,1)-6));
                %     temp_power = FOOOF_TABLE.theta{inds(k)}(:,2);
                %     FOOOF_TABLE.theta_center{inds(k)} = [FOOOF_TABLE.theta{inds(k)}(indx,1) temp_power(indx)];
                % end
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
                if length(FOOOF_TABLE.theta{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.theta{inds(k)}(:,1)-6));
                    temp_power = FOOOF_TABLE.theta{inds(k)}(:,2);
                    FOOOF_TABLE.theta_center{inds(k)} = [FOOOF_TABLE.theta{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.alpha1_band{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.alpha1_band{inds(k)}(:,1)-10));
                    temp_power = FOOOF_TABLE.alpha1_band{inds(k)}(:,2);
                    FOOOF_TABLE.alpha1_center{inds(k)} = [FOOOF_TABLE.alpha1_band{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.alpha2_band{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.alpha2_band{inds(k)}(:,1)-20));
                    temp_power = FOOOF_TABLE.alpha2_band{inds(k)}(:,2);
                    FOOOF_TABLE.alpha2_center{inds(k)} = [FOOOF_TABLE.alpha2_band{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.beta1_band{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.beta1_band{inds(k)}(:,1)-6));
                    temp_power = FOOOF_TABLE.beta1_band{inds(k)}(:,2);
                    FOOOF_TABLE.beta1_center{inds(k)} = [FOOOF_TABLE.beta1_band{inds(k)}(indx,1) temp_power(indx)];
                end
                if length(FOOOF_TABLE.beta2_band{inds(k)}) > 1
                    [~,indx] = min(abs(FOOOF_TABLE.beta2_band{inds(k)}(:,1)-10));
                    temp_power = FOOOF_TABLE.beta2_band{inds(k)}(:,2);
                    FOOOF_TABLE.beta2_center{inds(k)} = [FOOOF_TABLE.beta2_band{inds(k)}(indx,1) temp_power(indx)];
                end
            end
        end
    end
end
%## SAVE DATA
% par_save([save_dir filesep 'fooof_results_summary.mat'],'fooof_group_results_org');
par_save(fooof_diff_store,[save_dir filesep 'fooof_diff_store.mat']);
par_save(fooof_apfit_store,[save_dir filesep 'fooof_apfit_store.mat']);
par_save(spec_data_original,[save_dir filesep 'spec_data_original.mat']);
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
            %-
            if isempty(FOOOF_TABLE.alpha1_band{inds(k)})
                FOOOF_TABLE.alpha1_band{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta1_band{inds(k)})
                FOOOF_TABLE.beta1_band{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.alpha1_center{inds(k)})
                FOOOF_TABLE.alpha1_center{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta1_center{inds(k)})
                FOOOF_TABLE.beta1_center{inds(k)} = [NaN, NaN];
            end
            %-
            if isempty(FOOOF_TABLE.alpha2_band{inds(k)})
                FOOOF_TABLE.alpha2_band{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta2_band{inds(k)})
                FOOOF_TABLE.beta2_band{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.alpha2_center{inds(k)})
                FOOOF_TABLE.alpha2_center{inds(k)} = [NaN, NaN];
            end
            if isempty(FOOOF_TABLE.beta2_center{inds(k)})
                FOOOF_TABLE.beta2_center{inds(k)} = [NaN, NaN];
            end
            [~,idx_a] = max(FOOOF_TABLE.alpha{inds(k)}(:,2));
            [~,idx_b] = max(FOOOF_TABLE.beta{inds(k)}(:,2));
            [~,idx_a] = max(FOOOF_TABLE.alpha1_band{inds(k)}(:,2));
            [~,idx_b] = max(FOOOF_TABLE.beta1_band{inds(k)}(:,2));
            [~,idx_a] = max(FOOOF_TABLE.alpha2_band{inds(k)}(:,2));
            [~,idx_b] = max(FOOOF_TABLE.beta2_band{inds(k)}(:,2));
        end
    end
end
FOOOF_TABLE.subj_id = categorical(FOOOF_TABLE.subj_id);
FOOOF_TABLE.cond_id = categorical(FOOOF_TABLE.cond_id);
FOOOF_TABLE.cluster_id = categorical(FOOOF_TABLE.cluster_id);
FOOOF_TABLE.design_id = categorical(FOOOF_TABLE.design_id);
FOOOF_TABLE.group_id = categorical(FOOOF_TABLE.group_id);
writetable(FOOOF_TABLE,[save_dir filesep 'fooof_spec_table.xlsx']);
par_save(FOOOF_TABLE,[save_dir filesep 'psd_feature_table.mat']);
%## LOAD
%{
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
%}
%% Perform time series stats on the flattened curve
%## STATS
try
    STUDY.etc = rmfield(STUDY.etc,'statistics');
end
STUDY = pop_statparams(STUDY,...
    'groupstats','off',...
    'condstats','on',...
    'method','perm',...
    'singletrials','off',...
    'mode','fieldtrip',...
    'fieldtripalpha',NaN,...
    'fieldtripmethod','montecarlo',...
    'fieldtripmcorrect','fdr',...
    'fieldtripnaccu',4000);
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
        %## Ho : all samples come from the same distribution
        %## Ha : all samples come from different distributions
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
% par_save([save_dir filesep 'fooof_psd_stats.mat'],'table_out');
par_save(pcond,[save_dir filesep 'fooof_pcond.mat']);
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
par_save(pcond_org,[save_dir filesep 'org_pcond.mat']);
