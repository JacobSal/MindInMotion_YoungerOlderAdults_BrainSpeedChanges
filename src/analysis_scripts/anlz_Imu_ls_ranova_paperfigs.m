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
    STUDY_DIR = fileparts(SCRIPT_DIR);
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

% SUBJ_PICS = {ALLEEG.subject};
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
colormap(linspecer);
% study_dir_name = '03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
% study_dir_name = '04232024_MIM_OAN57_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';

%## soft define
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
% save_dir = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'raw_data_vis'];
save_dir = [studies_fpath filesep study_dir_name filesep 'behavioral_data'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
%##
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_mri_study');
%- override with already processed ALLEEG set
SUBJ_PICS = cell(1,length(GROUP_NAMES));
SUBJ_ITERS = cell(1,length(GROUP_NAMES));
for i = 1:length(STUDY.datasetinfo)
    gi = find(strcmp(STUDY.datasetinfo(i).group,GROUP_NAMES));
    SUBJ_PICS{1,gi} = [SUBJ_PICS{1,gi}, {STUDY.datasetinfo(i).subject}];
    if isempty(SUBJ_ITERS{1,gi})
        SUBJ_ITERS{1,gi} = 1;
    else
        SUBJ_ITERS{1,gi} = [SUBJ_ITERS{1,gi}, SUBJ_ITERS{1,gi}(end)+1];
    end
end
%%
% CATEGORIES = {'YoungAdult','HF_OlderAdult','LF_OlderAdult'};
CATEGORIES = {'YoungAdult','HF_OlderAdult','LF_OlderAdult'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
N_TRIALS = 16;
%- Loop through directory
% table_imu_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
table_imu_meas = nan(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
% table_trial_imu = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_trial_imu = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
table_trialtag_imu = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
% table_subj_imu = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_imu = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
table_header_names_imu = cell(length([SUBJ_PICS{:}]),1);
% table_subj_cat_imu = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_cat_imu = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
groupid_imu = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
designid_imu = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
fname_imu = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
%-
% table_ls_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
table_ls_meas = nan(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
% table_trial_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_trial_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
table_trialtag_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
% table_subj_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
table_header_names_ls = cell(length([SUBJ_PICS{:}]),1);
% table_subj_cat_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_cat_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
groupid_ls = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
designid_ls = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
fname_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% %- Loop through directory
% table_imu_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
% table_trial_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% table_subj_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% table_header_names_imu_conds = cell(length([SUBJ_PICS{:}]),1);
% table_subj_cat_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% groupid_imu_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
% designid_imu_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
% fname_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% %-
% table_ls_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
% table_trial_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% table_subj_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% table_header_names_ls_conds = cell(length([SUBJ_PICS{:}]),1);
% table_subj_cat_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% groupid_ls_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
% designid_ls_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
% fname_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
% subj_stack = [];

%%
cnt_imu = 1;
cnt_ls = 1;
fid = fopen([save_dir filesep 'load_in.txt'],'w');
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        %## LOAD CONDITION WISE BEHAVIORS
        folder_imu_conds = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Conditions'];
        folder_ls_conds = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Conditions'];
        dir_imu_conds = dir([folder_imu_conds filesep '*.mat']);
        dir_ls_conds = dir([folder_ls_conds filesep '*.mat']);
        %## LOAD TRIAL WISE BEHAVIORS
        folder_imu = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Trials'];
        folder_ls = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Trials'];
        dir_imu = dir([folder_imu filesep 'Outcome_Measures' filesep '*.mat']);
        dir_ls = dir([folder_ls filesep 'Outcome_Measures' filesep '*.mat']);
        fprintf(fid,'Loading subject %s...\n',SUBJ_PICS{group_i}{subj_i});
        vals = zeros(length(dir_imu),2);
        %## PRINTS
        fprintf(fid,'\n');
        fprintf(fid,'%s) has %i GRF files...\n',SUBJ_PICS{group_i}{subj_i},length(dir_ls));
        fprintf(fid,'%s) has %i IMU files...\n',SUBJ_PICS{group_i}{subj_i},length(dir_imu));
        %## IMU PER TRIAL (SACRAL MEASURES ONLY) 
        for f_i = 1:length(dir_imu)
            trial_i = [];
            cat_i = [];
            des_i = [];
            group_id = [];
            %-
            tmp = load([dir_imu(f_i).folder filesep dir_imu(f_i).name]);
            tmp = tmp.myStructure;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_imu{cnt_imu} = fieldnames(tmp);
            table_imu_meas(cnt_imu,1:length(values)) = values;
            tmp = strsplit(dir_imu(f_i).name,'.');
            %=
            fprintf(fid,'IMU Assigning %s.\n',dir_imu(f_i).name);
            if any(strcmp(tmp{1},{'TM_med_2','TM_med_1'}))
                trial_i = 'med';
                trial_t = 3;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
                trial_t = 2;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
                trial_t = 1;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
                trial_t = 4;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0.25';
                trial_t = 1;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0.5';
                trial_t = 2;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0.75';
                trial_t = 3;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1.0';
                trial_t = 4;
                des_i = 2;
            else
                trial_i = tmp{1};
                des_i = nan();
            end
            %-
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
                group_id = 1;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
                group_id = 2;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
                group_id = 3;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
                group_id = 3;
            else
                cat_i = nan();
                group_id = nan();
            end
            %-
            table_subj_cat_imu(cnt_imu) = categorical({cat_i});
            table_trial_imu(cnt_imu) = categorical({trial_i}); %tmp{1};
            table_subj_imu(cnt_imu) = categorical({SUBJ_PICS{group_i}{subj_i}});
            fname_imu{cnt_imu} = dir_imu(f_i).name;
            groupid_imu(cnt_imu) = group_id;
            designid_imu(cnt_imu) = des_i;
            table_trialtag_imu(cnt_imu) = categorical(trial_t);
            cnt_imu = cnt_imu + 1;
        end
        %## LOADSOL PER TRIAL AVERAGES (GAIT)
        for f_i = 1:length(dir_ls)
            trial_i = [];
            cat_i = [];
            des_i = [];
            group_id = [];
            %-
            tmp = load([dir_ls(f_i).folder filesep dir_ls(f_i).name]);
            tmp = tmp.summarizedTrial;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_ls{cnt_ls} = fieldnames(tmp);
            table_ls_meas(cnt_ls,1:length(values)) = values;
            tmp = strsplit(dir_ls(f_i).name,'.');
            fprintf(fid,'GRF Assigning %s.\n',dir_ls(f_i).name);
            %-
            if any(strcmp(tmp{1},{'TM_med_2','TM_med_1'}))
                trial_i = 'med';
                trial_t = 3;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
                trial_t = 2;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
                trial_t = 1;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
                trial_t = 4;
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0.25';
                trial_t = 1;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0.5';
                trial_t = 2;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0.75';
                trial_t = 3;
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1.0';
                trial_t = 4;
                des_i = 2;
            else
                trial_i = tmp{1};
                des_i = nan();
            end
            %-
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
                group_id = 1;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
                group_id = 2;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
                group_id = 3;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
                group_id = 3;
            else
                cat_i = nan();
                group_id = nan();
            end
            table_trial_ls(cnt_ls) = categorical({trial_i}); %tmp{1};
            table_subj_ls(cnt_ls) = categorical({SUBJ_PICS{group_i}{subj_i}});
            table_subj_cat_ls(cnt_ls) = categorical({cat_i});
            fname_ls{cnt_ls} = dir_ls(f_i).name;
            groupid_ls(cnt_ls) = group_id;
            designid_ls(cnt_ls) = des_i;
            table_trialtag_ls(cnt_ls) = categorical(trial_t);
            cnt_ls = cnt_ls + 1;
        end
    end
end
fclose(fid);
%% IMU TABLE
%-
% rej = ~any(table_imu_meas == 0,2);
% meas_in = table_imu_meas(rej,:);
% subj_in = table_subj_imu(rej);
% trial_in = table_trial_imu(rej);
% subj_cat_in = table_subj_cat_imu(rej);
% fname_in = fname_imu(rej);
% group_id_in  = groupid_imu(rej);
% design_id_in = designid_imu(rej);
%-
rej = ~all(isnan(table_imu_meas),2);
meas_in = table_imu_meas(rej,:);
subj_in = table_subj_imu(rej);
trial_in = table_trial_imu(rej);
trial_tag_in = table_trialtag_imu(rej);
subj_cat_in = table_subj_cat_imu(rej);
fname_in = fname_imu(rej);
group_id_in  = groupid_imu(rej);
design_id_in = designid_imu(rej);
%-
table_imu_out = array2table(meas_in,'VariableNames',table_header_names_imu{1});
table_imu_out.subj_char = subj_in;
table_imu_out.trial_char = trial_in;
table_imu_out.trial_tag_in = trial_tag_in;
table_imu_out.subj_cat = subj_cat_in;
table_imu_out.file_name = fname_in;
table_imu_out.group_id = group_id_in;
table_imu_out.design_id = design_id_in;
writetable(table_imu_out,[save_dir filesep 'imu_table_out.xlsx']);
save([save_dir filesep 'imu_table_out.mat'],'table_imu_out');
%% LOADSOL TABLE
% rej = ~any(table_ls_meas == 0,2);
% meas_in = table_ls_meas(rej,:);
% subj_in = table_subj_ls(rej);
% trial_in = table_trial_ls(rej);
% subj_cat_in = table_subj_cat_ls(rej);
% fname_in = fname_ls(rej);
% group_id_in  = groupid_ls(rej);
% design_id_in = designid_ls(rej);
%-
rej = ~all(isnan(table_ls_meas),2);
meas_in = table_ls_meas(rej,:);
subj_in = table_subj_ls(rej);
trial_in = table_trial_ls(rej);
trial_tag_in = table_trialtag_ls(rej);
subj_cat_in = table_subj_cat_ls(rej);
fname_in = fname_ls(rej);
group_id_in  = groupid_ls(rej);
design_id_in = designid_ls(rej);
%-
table_ls_out = array2table(meas_in,'VariableNames',table_header_names_ls{1});
table_ls_out.subj_char = subj_in;
table_ls_out.trial_char = trial_in;
table_ls_out.trial_tag_in = trial_tag_in;
table_ls_out.subj_cat = subj_cat_in;
table_ls_out.file_name = fname_in;
table_ls_out.group_id = group_id_in;
table_ls_out.design_id = design_id_in;
writetable(table_ls_out,[save_dir filesep 'ls_table_out.xlsx']);
save([save_dir filesep 'ls_table_out.mat'],'table_ls_out')
%% LOAD IN IF AVAILABLE
table_ls_out = load([save_dir filesep 'ls_table_out.mat']);
table_ls_out = table_ls_out.table_ls_out;
table_imu_out = load([save_dir filesep 'imu_table_out.mat']);
table_imu_out = table_imu_out.table_imu_out;
%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
table_new_imu = table_imu_out;
table_new_ls = table_ls_out;
% SPEED_CUTOFF = 0.1;
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
table_new_imu.terrain_speed = zeros(size(table_new_imu,1),1);
table_new_ls.terrain_speed = zeros(size(table_new_ls,1),1);
for i = 1:size(speed_table,1)
    ss = speed_table.Var1(i);
    ss_speed = speed_table.Var2(i);
%     ss = table_new_imu.subj_char(i);
    ind1 = table_new_imu.subj_char==ss;
    ind2 = table_new_ls.subj_char==ss;
    if ss_speed < SPEED_CUTOFF
%         inds_del = table_new_imu.subj_char == ss;
%         table_new_imu(inds_del) = [];
        table_new_imu(ind1,:) = [];
        table_new_ls(ind2,:) = [];
    else
        table_new_imu.terrain_speed(ind1) = ss_speed;
        table_new_ls.terrain_speed(ind2) = ss_speed;
    end
end
%%
MasterTable = mim_read_master_sheet('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\subject_mastersheet_notes.xlsx');
tmp = fieldnames(MasterTable(:,22:38));
tmp = tmp(1:end-3);
tmp_table = MasterTable(:,22:39);
og_length = size(tmp_table,2);
tmp_table.subject_code = categorical(MasterTable.subject_code);
tmp_table.stability_score = MasterTable.flat_low_med_high_rating_of_stability;
og_flat_fn = {'flat_t1_s','flat_t2_s','flat_t3_s'};
og_low_fn = {'low_t1_s','low_t2_s','low_t3_s'};
og_med_fn = {'med_t1_s','med_t2_s','med_t3_s'};
og_high_fn = {'high_t1_s','high_t2_s','high_t3_s'};
og_avg_fn = {og_flat_fn,og_low_fn,og_med_fn,og_high_fn};
% table_new_imu.stability_rating = zeros(size(table_new_imu,1),1);
% table_new_ls.stability_rating = zeros(size(table_new_ls,1),1);
% tmp_table_imu = table_imu_out;
% tmp_table_ls = table_ls_out;
tmp_table_imu = table_new_imu;
tmp_table_ls = table_new_ls;
trials = {'flat','low','med','high'};
%##
for i = 1:size(tmp_table,1)
    %- get subject index
    ss = tmp_table.subject_code(i);
    ss_var = tmp_table.stability_score(i);
    ss_var = strsplit(ss_var{1},';');
    ss_var = ss_var(~cellfun(@isempty,ss_var));
%     ss = table_new_imu.subj_char(i);
    ind1 = tmp_table_imu.subj_char==ss;
    ind2 = tmp_table_ls.subj_char==ss;
    if any(ind1) & ~isempty(ss_var)
        for j = 1:length(trials)
            sub = strsplit(ss_var{j},',');
            dsub = double(string(sub));
            %- get subject and trial type index in imu/ls sheets
            ind11 = tmp_table_imu.subj_char==ss & tmp_table_imu.trial_char==trials{j};
            ind22 = tmp_table_ls.subj_char==ss & tmp_table_ls.trial_char==trials{j};
            if ~isempty(sub{1})
                k_imu = find(ind11);
                k_ls = find(ind22);
                for k = 1:length(k_ls)
                    if any(ind22)
                        tmp_table_ls.(sprintf('%s','stability_rating'))(k_ls(k)) = double(string(sub{k}));
                    else
                        tmp_table_ls.(sprintf('%s','stability_rating'))(k_ls(k)) = nan();
                    end
                end
                for k = 1:length(k_imu)
                    if any(ind11)
                        tmp_table_imu.(sprintf('%s','stability_rating'))(k_imu(k)) = double(string(sub{k}));
                    else
                        tmp_table_imu.(sprintf('%s','stability_rating'))(k_imu(k)) = nan();
                    end
                end
            else
                for k = 1:2
                    tmp_table_imu.(sprintf('%s','stability_rating'))(k_imu(k)) = nan();
                    tmp_table_ls.(sprintf('%s','stability_rating'))(k_ls(k)) = nan();
                end
            end
            fn = fieldnames(tmp_table);
            for f = 1:length(fn)-5
                tmp_table_imu(ind11,fn{f}) = tmp_table(tmp_table.subject_code==ss,f);
                tmp_table_ls(ind22,fn{f}) = tmp_table(tmp_table.subject_code==ss,f);
            end
            %-
            %{
            if strcmp(trials{j},'flat')
                tmp_table_imu(ind11,'og_flat_mean_ms') = mean(tmp_table_imu(ind11,og_flat_fn(:)),2)./3;
                tmp_table_ls(ind22,'og_flat_mean_ms') = mean(tmp_table_ls(ind22,og_flat_fn(:)),2)./3;
            end
            if strcmp(trials{j},'low')
                tmp_table_imu(ind11,'og_low_mean_ms') = mean(tmp_table_imu(ind11,og_low_fn(:)),2)./3;
                tmp_table_ls(ind22,'og_low_mean_ms') = mean(tmp_table_ls(ind22,og_low_fn(:)),2)./3;
            end
            if strcmp(trials{j},'med')
                tmp_table_imu(ind11,'og_med_mean_ms') = mean(tmp_table_imu(ind11,og_med_fn(:)),2)./3;
                tmp_table_ls(ind22,'og_med_mean_ms') = mean(tmp_table_ls(ind22,og_med_fn(:)),2)./3;
            end
            if strcmp(trials{j},'high')
                tmp_table_imu(ind11,'og_high_mean_ms') = mean(tmp_table_imu(ind11,og_high_fn(:)),2)./3;
                tmp_table_ls(ind22,'og_high_mean_ms') = mean(tmp_table_ls(ind22,og_high_fn(:)),2)./3;
            end
            %}
            %-
            if strcmp(trials{j},'flat')
                tmp_table_imu(ind11,'og_mean_ms') = 3./mean(tmp_table_imu(ind11,og_flat_fn(:)),2);
                tmp_table_ls(ind22,'og_mean_ms') = 3./mean(tmp_table_ls(ind22,og_flat_fn(:)),2);
            end
            if strcmp(trials{j},'low')
                tmp_table_imu(ind11,'og_mean_ms') = 3./mean(tmp_table_imu(ind11,og_low_fn(:)),2);
                tmp_table_ls(ind22,'og_mean_ms') = 3./mean(tmp_table_ls(ind22,og_low_fn(:)),2);
            end
            if strcmp(trials{j},'med')
                tmp_table_imu(ind11,'og_mean_ms') = 3./mean(tmp_table_imu(ind11,og_med_fn(:)),2);
                tmp_table_ls(ind22,'og_mean_ms') = 3./mean(tmp_table_ls(ind22,og_med_fn(:)),2);
            end
            if strcmp(trials{j},'high')
                tmp_table_imu(ind11,'og_mean_ms') = 3./mean(tmp_table_imu(ind11,og_high_fn(:)),2);
                tmp_table_ls(ind22,'og_mean_ms') = 3./mean(tmp_table_ls(ind22,og_high_fn(:)),2);
            end
            tmp_table_imu(ind11,'og_length_m') = {3};
            tmp_table_ls(ind22,'og_length_m') = {3};
        end
    end
end
%##
subjects = unique(tmp_table_imu.subj_char);
for i = 1:size(subjects)
    %- get subject index
    ss = subjects(i);
    ind1 = tmp_table_imu.subj_char==ss;
    ind2 = tmp_table_ls.subj_char==ss;
    %-
    if any(ind1) & ~isempty(ss_var)
        ind11 = find(ind1);
        ind22 = find(ind2);
        %-
        tmp_table_imu(end+1,:) = tmp_table_imu(ind11(1),:);
        fn = fieldnames(tmp_table_imu);
        for ff = 1:length(fn)-3
            try
                tmp_table_imu(end,ff) = {nan()};
            catch ME
                if strcmp(ME,'MATLAB:invalidConversion')
                    tmp_table_imu(end,ff) = categorical({''});
                end
            end
        end
        tmp_table_imu(end,:).trial_char = categorical({'og_avg'});
        tmp_table_imu(end,:).design_id = 1;
        tmp_table_imu(end,:).group_id = tmp_table_imu(ind11(1),:).group_id;
        tmp_table_imu(end,:).file_name = {''};
        tmp_table_imu(end,:).trial_tag_in = categorical(5);
        tmp_table_imu(end,:).og_mean_ms = tmp_table_imu(ind11(1),:).terrain_speed;
        %-
        tmp_table_ls(end+1,:) = tmp_table_ls(ind22(1),:);
        fn = fieldnames(tmp_table_ls);
        for ff = 1:length(fn)-3
            try
                tmp_table_ls(end,ff) = {nan()};
            catch ME
                if strcmp(ME,'MATLAB:invalidConversion')
                    tmp_table_ls(end,ff) = categorical({''});
                end
            end
        end
        tmp_table_ls(end,:).trial_char = categorical({'og_avg'});
        tmp_table_ls(end,:).design_id = 1;
        tmp_table_ls(end,:).group_id = tmp_table_ls(ind11(1),:).group_id;
        tmp_table_ls(end,:).file_name = {''};
        tmp_table_ls(end,:).trial_tag_in = categorical(5);
        tmp_table_ls(end,:).og_mean_ms = tmp_table_ls(ind11(1),:).terrain_speed;
    end
end
%##
fn = fieldnames(tmp_table);
for f = 1:length(fn)-6
    try
        tmp_table_imu(tmp_table_imu.(fn{f}) == 0,fn{f}) = {nan()};
        tmp_table_ls(tmp_table_ls.(fn{f}) == 0,fn{f}) = {nan()};
    catch ME
        if (strcmp(ME.identifier,'MATLAB:datetime:InvalidComparison'))
            fprintf('NaT detected, no value change needed...\n');
        end
    end
end
tmp_table_imu.stability_rating(tmp_table_imu.stability_rating==0) = nan();
tmp_table_ls.stability_rating(tmp_table_imu.stability_rating==0) = nan();
tmp_table_imu.og_length_m(tmp_table_imu.og_length_m==0) = nan();
tmp_table_ls.og_length_m(tmp_table_imu.og_length_m==0) = nan();
writetable(tmp_table_imu,[save_dir filesep 'imu_table_trial.xlsx']);
writetable(tmp_table_ls,[save_dir filesep 'ls_table_trial.xlsx']);
save([save_dir filesep 'imu_table_trial.mat'],'tmp_table_imu')
save([save_dir filesep 'ls_table_trial.mat'],'tmp_table_ls')
%##
table_imu_out = tmp_table_imu;
table_ls_out = tmp_table_ls;
%% average across trials
subjs = unique(table_imu_out.subj_char);
trials = unique(table_imu_out.trial_char);
table_new = [];
for i = 1:length(subjs)
    % tmpvals = table_imu_out(table_imu_out.subj_char == subjs(i),:);
    % subj_cat = unique(table_imu_out(table_imu_out.subj_char == subjs(i),:).subj_cat);
    for j = 1:length(trials)
        ind = table_imu_out.subj_char == subjs(i) & table_imu_out.trial_char == trials(j);
        tmp = table_imu_out(ind,:); %tmpvals(tmpvals.trial_char == trials(j),:);
%         tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
        if ~isempty(tmp)
            tmpvals = varfun(@mean, tmp, 'InputVariables', @isnumeric);
            tmpvals.subj_char = tmp.subj_char(1);
            tmpvals.trial_char = tmp.trial_char(1);
            tmpvals.subj_cat = tmp.subj_cat(1);
            tmpvals.file_name = tmp.file_name(1);
            tmpvals.group_id = tmp.group_id(1);
            tmpvals.design_id = tmp.design_id(1);
            tmpvals.trial_tag_in = tmp.trial_tag_in(1);
            table_new = [table_new; tmpvals];
        end
    end
end
table_new_imu = table_new;

%% average across trials
subjs = unique(table_ls_out.subj_char);
trials = unique(table_ls_out.trial_char);
table_new = [];
for i = 1:length(subjs)
    % tmpvals = table_ls_out(table_ls_out.subj_char == subjs(i),:);
    % subj_cat = unique(table_ls_out(table_ls_out.subj_char == subjs(i),:).subj_cat);
    for j = 1:length(trials)
        % tmp = tmpvals(tmpvals.trial_char == trials(j),:);
%         tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
        ind = table_ls_out.subj_char == subjs(i) & table_ls_out.trial_char == trials(j);
        tmp = table_ls_out(ind,:); %tmpvals(tmpvals.trial_char == trials(j),:);
        if ~isempty(tmp)
            tmpvals = varfun(@mean, tmp, 'InputVariables', @isnumeric);
            tmpvals.subj_char = tmp.subj_char(1);
            tmpvals.trial_char = tmp.trial_char(1);
            tmpvals.subj_cat = tmp.subj_cat(1);
            tmpvals.file_name = tmp.file_name(1);
            tmpvals.group_id = tmp.group_id(1);
            tmpvals.design_id = tmp.design_id(1);
            tmpvals.trial_tag_in = tmp.trial_tag_in(1);
            table_new = [table_new; tmpvals];
        end
    end
end
table_new_ls = table_new;
%% (SAVE) ============================================================== %%
writetable(table_new_imu,[save_dir filesep 'imu_table_meantrial.xlsx']);
writetable(table_new_ls,[save_dir filesep 'ls_table_meantrial.xlsx']);
save([save_dir filesep 'ls_table_meantrial.mat'],'table_new_ls')
save([save_dir filesep 'imu_table_meantrial.mat'],'table_new_imu')
%% (LOAD IF AVAILABLE) ================================================= %%
table_new_ls = load([save_dir filesep 'ls_table_meantrial.mat']);
table_new_ls = table_new_ls.table_new_ls;
table_new_imu = load([save_dir filesep 'imu_table_meantrial.mat']);
table_new_imu = table_new_imu.table_new_imu;
%% (STATS STRUCT) ====================================================== %%
% DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
%     'measure_tag',categorical({''}),...
%     'design_tag',categorical({''}),...
%     'mod_tag',categorical({''}),...
%     'mod_resp_terms',{''},...
%     'rnd_terms',{''},...
%     'anova_preds_terms',{''},...
%     'anova_preds_p',{[]},...
%     'anova_preds_stat',{[]},...
%     'anova_preds_df',{[]},...
%     'mod_preds_terms',{{''}},...
%     'mod_preds_p',[],...
%     'mod_preds_stat',[],...
%     'mod_preds_coeff',[],...
%     'mod_r2',[],...
%     'multi_comp_terms',{''},...
%     'multi_comp_t1_t2',[],...
%     'multi_comp_p',[],...
%     'multi_comp_coeff',[],...
%     'multi_comp_lci_uci',[],...
%     'norm_test_p',[],...
%     'norm_test_h',[],...
%     'effect_size',[],...
%     'effect_size_calc',{''});
% stats_struct = DEF_STATS_TRACK_STRUCT;
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
%% PLOT ================================================================ %%
%## VISUAL
%- colors
COLORS_MAPS_TERRAIN = linspecer(4);
YELLOW = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);YELLOW;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%-
IM_RESIZE = 1.1;
AX_H  = 0.2;
AX_W = 0.275;
AX_HORIZ_SHIFT = 0.1;
AX_INIT_HORIZ = 0.08;
AX_INIT_VERT = 0.6;
AX_VERT_SHIFT = 0.125;
%-
GROUP_LAB_FONTSIZE = 10;
XLAB_FONTSIZE = 10;
YLAB_FONTSIZE = 10;
XTICK_FONTSIZE = 7;
TITLE_FONTSIZE = 10;
FONT_SIZE_VIO = 10;
FONT_SIZE_VIO_REG = 7;
%-
FONT_NAME = 'Arial';
GROUP_LAB_FONTWEIGHT = 'bold ';
XLAB_FONTWEIGHT = 'bold';
YLAB_FONTWEIGHT = 'bold';
TITLE_FONTWEIGHT = 'bold';
%-
%-
LAB_A_YOFFSET = -0.1;
LAB_A_XOFFSET = -0.125;
LAB_B_YOFFSET = 0.12;
LAB_B_XOFFSET = -0.12;
LAB_C_YOFFSET = 0.12;
LAB_C_XOFFSET = -0.12;
LAB_D_YOFFSET = 0.12;
LAB_D_XOFFSET = -0.12;
XLABEL_OFFSET = .115;
GROUP_LAB_YOFFSET = -0.29;
FIGURE_POSITION =[0,0,6.5,9];
TITLE_XSHIFT = 0.4;
TITLE_YSHIFT = 0.975;
TITLE_BOX_SZ = [0.4,0.4];
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
    'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%## FACTORS
meas_names = {'mean_APexc_COV','mean_MLexc_COV','mean_StepDur_cov','mean_og_mean_ms'};
meas_ylabel = {'COV (%)','COV (%)','COV (%)','Walking Speed (m/s)'};
meas_titles = {{'Anteroposterior Excursion';'Coefficient of Variation'},...
    {'Mediolateral Excursion';'Coefficient of Variation'},...
    {'Step Duration';'Coefficient of Variation'},...
    {'Overground Walking';'Speeds'}};
% YLIMS = {[0,90],[0,95]};
YLIMS = {[0,55],[0,40],[0,22.5],[0,2.25]};
%## TEMPS
horiz_shift = 0;
vert_shift = 0;
des_i = 2;
cnts = 1;
%% FIGURE INIT
fig = figure('color','white','renderer','Painters');
% annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
%     'String',atlas_name,'HorizontalAlignment','center',...
%     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
%     'FontSize',14,'FontWeight','Bold','Units','normalized');
set(fig,'Units','inches','Position',FIGURE_POSITION)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
set(gca,AXES_DEFAULT_PROPS{:});
%% (SPEED) OG SPEEDS - MAIN EFFECTS ================================ %%
meas_i = 4;
des_i = 1;
conds_keep = {1};
%-
tmp = (cellfun(@(x) table_new_ls.trial_tag_in == categorical(x),conds_keep,'UniformOutput',false));
tmp = any(cat(2,tmp{:}),2);
inds = table_new_ls.design_id == des_i & tmp;
tmp_table = table_new_ls(inds,:);
tmp_table = table(categorical(string(tmp_table.trial_char)),categorical(double(tmp_table.trial_tag_in)),...
    categorical(string(tmp_table.subj_char)),...
    categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
    'VariableNames',{'cond_char','trial_tag_in','subj_char','group_id',meas_names{meas_i}});
%##
tmp_table(tmp_table.(meas_names{meas_i}) == 0,:) = [];
tmp_table(isnan(tmp_table.(meas_names{meas_i})),:) = [];
tmp_table(isinf(tmp_table.(meas_names{meas_i})),:) = [];
%##
% mod_out = sprintf('%s~1+group_id+trial_tag_in',meas_names{meas_i});
mod_out = sprintf('%s~1+group_id',meas_names{meas_i});
stats_out = fitlm(tmp_table,mod_out);
pred_terms = stats_out.CoefficientNames;
% anova_out = anova(stats_out);
[p,t,anova_out,terms] = anovan(tmp_table.(meas_names{meas_i}),{tmp_table.trial_tag_in, tmp_table.group_id},...
    'sstype',3,'varnames',{'trial_tag_in','group_id'},'model','linear','Display','off');
[comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[2],...
    'display','off','Alpha',0.05); % compaarisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
disp(stats_out)
%- test normality
[norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
%- intercept only model
altmod_out = sprintf('%s ~ 1',meas_names{meas_i});
altstats_out = fitlm(tmp_table,altmod_out);
%- alternative f2?
R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
alt_f2 = (R22-R21)/(1-R22);
%- populate struct
stats_struct(cnts).stat_test_mod = mod_out;
stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
stats_struct(cnts).design_tag = categorical(des_i);
stats_struct(cnts).mod_tag = categorical(2);
% stats_struct(cnts).resp_terms = meas_names(meas_i);
stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
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
stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Properties.RowNames';
stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
stats_struct(cnts).multi_comp_terms = gnames';
stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
[h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
stats_struct(cnts).multi_comp_p = adj_p;
stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
stats_struct(cnts).norm_test_p = norm_p;
stats_struct(cnts).norm_test_h = norm_h;
stats_struct(cnts).effect_size = alt_f2;
stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
cnts = cnts + 1;
stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
%-
tmp = t(:,7)';
tmp = tmp(~cellfun(@isempty,tmp));
anova_grp_p = tmp{3};
tmp = t(:,6)';
tmp = tmp(~cellfun(@isempty,tmp));
anova_grp_stat = tmp{3};
grp_p = [0,stats_out.Coefficients.pValue(2:3)']';
multi_comp_t1_t2 = comparisons(:,1:2);
speed_xvals = (0:5)*0.25;
if anova_grp_p > 0.01 && anova_grp_p < 0.05
    str = {sprintf('* F_{stat}=%0.3g',anova_grp_stat),'',''};
elseif anova_grp_p <= 0.01 && anova_grp_p > 0.001
    str = {sprintf('** F_{stat}=%0.3g',anova_grp_stat),'',''};
elseif anova_grp_p <= 0.001
    str = {sprintf('*** F_{stat}=%0.3g',anova_grp_stat),'',''};
else
    str = {'','',''};
end
STATS_STRUCT = struct('anova',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
    'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
    'pvals',{{}},...
    'pvals_pairs',{{}},...
    'pvals_grp',{num2cell(adj_p)},...
    'pvals_grp_pairs',{num2cell(multi_comp_t1_t2,2)},...
    'display_stats_char',true,...
    'stats_char',{str},...
    'group_order',categorical({''}),...
    'bracket_yshift_perc',0.5,...
    'bracket_y_perc',0.5,...
    'bracket_rawshifty_upper',0,...
    'bracket_rawshifty_lower',0,...
    'grp_sig_offset_x',0,...
    'grp_sig_offset_y',0);
%## PLOT =============================================================== %%
VIOLIN_PARAMS = {'width',0.1,...
    'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
    'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
    'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
%[-0.35,-0.1,0.15,0.40,0.65],
PLOT_STRUCT = struct('color_map',[COLORS_MAPS_TERRAIN;COLOR_MAPS_SPEED(2,:)*0.8],...
    'cond_labels',categorical({'flat','low','med','high','tm'}),'group_labels',categorical({'YA','OHMA','OLMA'}),...
    'cond_offsets',[-0.3375,-0.1125,0.1125,0.3375,0.5625],...
    'group_offsets',[0.125,0.475,0.812],...
    'y_label',meas_ylabel{meas_i},...
    'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
    'font_name','Arial','x_label','','do_combine_groups',false,...
    'regresslab_txt_size',FONT_SIZE_VIO_REG);
%##
ax = axes();
ax = group_violin(tmp_table,meas_names{meas_i},'trial_tag_in','group_id',...
    ax,...
    'VIOLIN_PARAMS',VIOLIN_PARAMS,...
    'PLOT_STRUCT',PLOT_STRUCT,...
    'STATS_STRUCT',STATS_STRUCT);
set(ax,'OuterPosition',[0,0,1,1]);
set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
%- set ylabel & title
% set(ax,'FontName','Arial','FontSize',XTICK_FONTSIZE,'FontWeight','normal')
% yticks(ax,'FontSize',XTICK_FONTSIZE)
yt = yticks(ax);
xlh = xlabel(ax,'','Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
pos1=get(xlh,'Position');
pos1(1,2)=pos1(1,2)+XLABEL_OFFSET;
set(xlh,'Position',pos1);
ylabel(ax,meas_ylabel{meas_i},'FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
title(ax,meas_titles{meas_i},...
    'FontSize',TITLE_FONTSIZE,'FontWeight',TITLE_FONTWEIGHT);%## LETTER
annotation('textbox',[AX_INIT_HORIZ+LAB_D_XOFFSET+(0.1/2)+horiz_shift,AX_INIT_VERT+vert_shift+LAB_B_YOFFSET+(0.1/2),.1,.1],...
    'String','A)','HorizontalAlignment','left',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',14,'FontWeight','Bold','Units','normalized');
%## SAVE
% exportgraphics(fig,[save_dir filesep 'figure_behavioral_data_paper.tiff'],'Resolution',1000);
%## SHIFT
horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
%% (SPEED) AP COV - MAIN EFFECTS ======================================= %%
meas_i = 1;
des_i = 2;
%-
inds = table_new_imu.design_id == des_i;
tmp_table = table_new_imu(inds,:);
tmp_table = table(double(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
    categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
    'VariableNames',{'cond_char','subj_char','group_id',meas_names{meas_i}});
%-
%##
tmp_table(tmp_table.(meas_names{meas_i}) == 0,:) = [];
tmp_table(isnan(tmp_table.(meas_names{meas_i})),:) = [];
tmp_table(isinf(tmp_table.(meas_names{meas_i})),:) = [];
%## RANOVA
tt_out = unstack(tmp_table,meas_names{meas_i},'cond_char');
mod_out = sprintf('x0_25,x0_5,x0_75,x1~group_id+1');
tmptmp = table([1,2,3,4]','VariableNames',{'speed'});
%- run test
rm_out = fitrm(tt_out,mod_out,'WithinDesign',tmptmp);
[ran_out,A,C,D] = ranova(rm_out);
%- test coeffs
disp(cat(1,A{:})*table2array(rm_out.Coefficients)*C)
[ran_mult] = rm_out.multcompare('group_id');
%- test normality
% [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
sphr = rm_out.mauchly(C);
sphr_corr = rm_out.epsilon(C);
%## RANOVA
%- populate struct
stats_struct(cnts).stat_test_mod = rm_out.BetweenModel;
stats_struct(cnts).measure_tag = categorical({meas_names{meas_i}});
stats_struct(cnts).measure_char = meas_names{meas_i};
stats_struct(cnts).design_tag = categorical(des_i);
stats_struct(cnts).design_num = des_i;
% stats_struct(cnts).cluster_num = cl_i;
stats_struct(cnts).mod_tag = categorical({'ranova_2'});
stats_struct(cnts).mod_resp_terms = rm_out.ResponseNames;
%- anova
stats_struct(cnts).anova_preds_terms = ran_out.Properties.RowNames';
stats_struct(cnts).anova_preds_df = ran_out.DF';
stats_struct(cnts).anova_preds_stat = ran_out.F;
%- liberal
% tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
% stats_struct(cnts).anova_preds_p = tmp(1:2)';
%- conservative
tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
stats_struct(cnts).anova_preds_p = tmp(1:2)';
%-
stats_struct(cnts).anova_preds_stat = ran_out.F;
%- model
stats_struct(cnts).mod_preds_coeff = rm_out.Coefficients; % remove from table before xlsx
stats_struct(cnts).mod_preds_terms = rm_out.BetweenFactorNames';
stats_struct(cnts).multi_comp_terms = ran_mult.Properties.VariableNames(1:2);
stats_struct(cnts).multi_comp_p = ran_mult.pValue;
stats_struct(cnts).multi_comp_t1_t2 = num2cell(ran_mult{:,1:2},2)';
stats_struct(cnts).multi_comp_lci_uci = [ran_mult.Lower,ran_mult.Upper];
% stats_struct(cnts).norm_test_p = norm_p;
% stats_struct(cnts).norm_test_h = norm_h;
cnts = cnts + 1;
stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
%## PLOT =============================================================== %%
%## RANOVA VISUALIZE
% tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
comparisons = [1,2;1,3;2,3];
[h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(ran_mult.pValue);
adj_p = [adj_p(1),adj_p(2),adj_p(4)];
anova_p = 0;
% anova_grp_p = ran_out.pValue(2);
anova_grp_p = tmp(2);
% ran_means_1 = [rm_out.Coefficients{1,:}].*[0.25,0.5,0.75,1];
% ran_means_2 = [rm_out.Coefficients{2,:}].*[0.25,0.5,0.75,1]+ran_means_1;
% ran_means_3 = [rm_out.Coefficients{3,:}].*[0.25,0.5,0.75,1]+ran_means_1;
ran_means_1 = [rm_out.Coefficients{1,:}];
ran_means_2 = [rm_out.Coefficients{2,:}]+ran_means_1;
ran_means_3 = [rm_out.Coefficients{3,:}]+ran_means_1;
if anova_grp_p > 0.01 && anova_grp_p < 0.05
    str = {sprintf('* pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.01 && anova_grp_p > 0.001
    str = {sprintf('** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.001
    str = {sprintf('*** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
else
    str = {'','',''};
end
STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(adj_p)},...
        'pvals_grp_pairs',{num2cell(comparisons,2)},...
        'regress_pval',{{0,0,0}},...
        'regress_line',{{ran_means_2,...
            ran_means_3,...
            ran_means_1}},...
        'line_type',{'means'},...
        'stats_char',{str},...
        'display_stats_char',true,...
        'bracket_yshift_perc',0.5,...
        'bracket_y_perc',0.5,...
        'bracket_rawshifty_upper',0,...
        'bracket_rawshifty_lower',0,...
        'grp_sig_offset_x',0,...
        'grp_sig_offset_y',0);

%## PLOT =============================================================== %%
inds = [stats_struct.mod_tag] == categorical(2) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
tmp_stats = stats_struct(inds);
VIOLIN_PARAMS = {'width',0.1,...
    'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
    'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
    'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
    'cond_labels',unique(tmp_table.cond_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_offsets',[0.125,0.475,0.812],...
    'y_label',meas_ylabel{meas_i},...
    'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
    'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
    'regresslab_txt_size',FONT_SIZE_VIO_REG);

%##
ax = axes();
ax = group_violin(tmp_table,meas_names{meas_i},'cond_char','group_id',...
    ax,...
    'VIOLIN_PARAMS',VIOLIN_PARAMS,...
    'PLOT_STRUCT',PLOT_STRUCT,...
    'STATS_STRUCT',STATS_STRUCT);
set(ax,'OuterPosition',[0,0,1,1]);
set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
%- set ylabel & title
% set(ax,'FontName','Arial','FontSize',XTICK_FONTSIZE,'FontWeight','normal')
% yticks(ax,'FontSize',XTICK_FONTSIZE)
yt = yticks(ax);
xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
pos1=get(xlh,'Position');
pos1(1,2)=pos1(1,2)+XLABEL_OFFSET;
set(xlh,'Position',pos1);
ylabel(ax,meas_ylabel{meas_i},'FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
title(ax,meas_titles{meas_i},...
    'FontSize',TITLE_FONTSIZE,'FontWeight',TITLE_FONTWEIGHT);
%## LETTER
annotation('textbox',[AX_INIT_HORIZ+LAB_A_XOFFSET+(0.1/2)+horiz_shift,AX_INIT_VERT+LAB_A_YOFFSET+AX_H*IM_RESIZE+(0.1/2),.1,.1],...
    'String','B)','HorizontalAlignment','left',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',14,'FontWeight','Bold','Units','normalized');
%## SHIFT
% horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
horiz_shift = 0;
vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
%% (SPEED) ML COV - MAIN EFFECTS ======================================= %%
meas_i = 2;
des_i = 2;
%-
inds = table_new_imu.design_id == des_i;
tmp_table = table_new_imu(inds,:);
tmp_table = table(double(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
    categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
    'VariableNames',{'cond_char','subj_char','group_id',meas_names{meas_i}});
%-
%##
tmp_table(tmp_table.(meas_names{meas_i}) == 0,:) = [];
tmp_table(isnan(tmp_table.(meas_names{meas_i})),:) = [];
tmp_table(isinf(tmp_table.(meas_names{meas_i})),:) = [];
%## RANOVA
tt_out = unstack(tmp_table,meas_names{meas_i},'cond_char');
mod_out = sprintf('x0_25,x0_5,x0_75,x1~group_id+1');
tmptmp = table([1,2,3,4]','VariableNames',{'speed'});
%- run test
rm_out = fitrm(tt_out,mod_out,'WithinDesign',tmptmp);
[ran_out,A,C,D] = ranova(rm_out);
%- test coeffs
disp(cat(1,A{:})*table2array(rm_out.Coefficients)*C)
[ran_mult] = rm_out.multcompare('group_id');
%- test normality
% [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
sphr = rm_out.mauchly(C);
sphr_corr = rm_out.epsilon(C);
%## RANOVA
%- populate struct
stats_struct(cnts).stat_test_mod = rm_out.BetweenModel;
stats_struct(cnts).measure_tag = categorical({meas_names{meas_i}});
stats_struct(cnts).measure_char = meas_names{meas_i};
stats_struct(cnts).design_tag = categorical(des_i);
stats_struct(cnts).design_num = des_i;
% stats_struct(cnts).cluster_num = cl_i;
stats_struct(cnts).mod_tag = categorical({'ranova_2'});
stats_struct(cnts).mod_resp_terms = rm_out.ResponseNames;
%- anova
stats_struct(cnts).anova_preds_terms = ran_out.Properties.RowNames';
stats_struct(cnts).anova_preds_df = ran_out.DF';
stats_struct(cnts).anova_preds_stat = ran_out.F';
%- liberal
% tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
% stats_struct(cnts).anova_preds_p = tmp(1:2)';
%- conservative
tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
stats_struct(cnts).anova_preds_p = tmp(1:2)';
%- model
stats_struct(cnts).mod_preds_coeff = rm_out.Coefficients; % remove from table before xlsx
stats_struct(cnts).mod_preds_terms = rm_out.BetweenFactorNames';
stats_struct(cnts).multi_comp_terms = ran_mult.Properties.VariableNames(1:2);
stats_struct(cnts).multi_comp_p = ran_mult.pValue;
stats_struct(cnts).multi_comp_t1_t2 = num2cell(ran_mult{:,1:2},2)';
stats_struct(cnts).multi_comp_lci_uci = [ran_mult.Lower,ran_mult.Upper];
% stats_struct(cnts).norm_test_p = norm_p;
% stats_struct(cnts).norm_test_h = norm_h;
cnts = cnts + 1;
stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
%## RANOVA VISUALIZE
% tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
comparisons = [1,2;1,3;2,3];
[h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(ran_mult.pValue);
adj_p = [adj_p(1),adj_p(2),adj_p(4)];
anova_p = 0;
% anova_grp_p = ran_out.pValue(2);
anova_grp_p = tmp(2);
% ran_means_1 = [rm_out.Coefficients{1,:}].*[0.25,0.5,0.75,1];
% ran_means_2 = [rm_out.Coefficients{2,:}].*[0.25,0.5,0.75,1]+ran_means_1;
% ran_means_3 = [rm_out.Coefficients{3,:}].*[0.25,0.5,0.75,1]+ran_means_1;
ran_means_1 = [rm_out.Coefficients{1,:}];
ran_means_2 = [rm_out.Coefficients{2,:}]+ran_means_1;
ran_means_3 = [rm_out.Coefficients{3,:}]+ran_means_1;
if anova_grp_p > 0.01 && anova_grp_p < 0.05
    str = {sprintf('* pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.01 && anova_grp_p > 0.001
    str = {sprintf('** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.001
    str = {sprintf('*** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
else
    str = {'','',''};
end
STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(adj_p)},...
        'pvals_grp_pairs',{num2cell(comparisons,2)},...
        'regress_pval',{{0,0,0}},...
        'regress_line',{{ran_means_2,...
            ran_means_3,...
            ran_means_1}},...
        'line_type',{'means'},...
        'stats_char',{str},...
        'display_stats_char',true,...
        'bracket_yshift_perc',0.5,...
        'bracket_y_perc',0.5,...
        'bracket_rawshifty_upper',0,...
        'bracket_rawshifty_lower',0,...
        'grp_sig_offset_x',0,...
        'grp_sig_offset_y',0);
%## PLOT =============================================================== %%
inds = [stats_struct.mod_tag] == categorical(2) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
tmp_stats = stats_struct(inds);
VIOLIN_PARAMS = {'width',0.1,...
    'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
    'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
    'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
    'cond_labels',unique(tmp_table.cond_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_offsets',[0.125,0.475,0.812],...
    'y_label',meas_ylabel{meas_i},...
    'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
    'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
    'regresslab_txt_size',FONT_SIZE_VIO_REG);

%##
ax = axes();
ax = group_violin(tmp_table,meas_names{meas_i},'cond_char','group_id',...
    ax,...
    'VIOLIN_PARAMS',VIOLIN_PARAMS,...
    'PLOT_STRUCT',PLOT_STRUCT,...
    'STATS_STRUCT',STATS_STRUCT);
set(ax,'OuterPosition',[0,0,1,1]);
set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
%- set ylabel & title
% set(ax,'FontName','Arial','FontSize',XTICK_FONTSIZE,'FontWeight','normal')
% yticks(ax,'FontSize',XTICK_FONTSIZE)
yt = yticks(ax);
xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
pos1=get(xlh,'Position');
pos1(1,2)=pos1(1,2)+XLABEL_OFFSET;
set(xlh,'Position',pos1);
ylabel(ax,meas_ylabel{meas_i},'FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
title(ax,meas_titles{meas_i},...
    'FontSize',TITLE_FONTSIZE,'FontWeight',TITLE_FONTWEIGHT);
%## LETTER
annotation('textbox',[AX_INIT_HORIZ+LAB_B_XOFFSET+(0.1/2)+horiz_shift,AX_INIT_VERT+vert_shift+LAB_B_YOFFSET+(0.1/2),.1,.1],...
    'String','C)','HorizontalAlignment','left',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',14,'FontWeight','Bold','Units','normalized');
%## SHIFT
% horiz_shift = 0;
% vert_shift = vert_shift - AX_H*IM_RESIZE - AX_VERT_SHIFT*IM_RESIZE;
horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
%% (SPEED) STEP DUR. COV - MAIN EFFECTS ================================ %%
meas_i = 3;
des_i = 2;
%-
inds = table_new_ls.design_id == des_i;
tmp_table = table_new_ls(inds,:);
tmp_table = table(double(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
    categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
    'VariableNames',{'cond_char','subj_char','group_id',meas_names{meas_i}});
% tmp_table = table(categorical(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
%     categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
%     'VariableNames',{'cond_char','subj_char','group_id',meas_names{meas_i}});
%-
%##
tmp_table(tmp_table.(meas_names{meas_i}) == 0,:) = [];
tmp_table(isnan(tmp_table.(meas_names{meas_i})),:) = [];
tmp_table(isinf(tmp_table.(meas_names{meas_i})),:) = [];
%## RANOVA
tt_out = unstack(tmp_table,meas_names{meas_i},'cond_char');
mod_out = sprintf('x0_25,x0_5,x0_75,x1~group_id+1');
tmptmp = table([1,2,3,4]','VariableNames',{'speed'});
%- run test
rm_out = fitrm(tt_out,mod_out,'WithinDesign',tmptmp);
[ran_out,A,C,D] = ranova(rm_out);
%- test coeffs
disp(cat(1,A{:})*table2array(rm_out.Coefficients)*C)
[ran_mult] = rm_out.multcompare('group_id');
%- test normality
% [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
sphr = rm_out.mauchly(C);
sphr_corr = rm_out.epsilon(C);
%## RANOVA
%- populate struct
stats_struct(cnts).stat_test_mod = rm_out.BetweenModel;
stats_struct(cnts).measure_tag = categorical({meas_names{meas_i}});
stats_struct(cnts).measure_char = meas_names{meas_i};
stats_struct(cnts).design_tag = categorical(des_i);
stats_struct(cnts).design_num = des_i;
% stats_struct(cnts).cluster_num = cl_i;
stats_struct(cnts).mod_tag = categorical({'ranova_2'});
stats_struct(cnts).mod_resp_terms = rm_out.ResponseNames;
%- anova
stats_struct(cnts).anova_preds_terms = ran_out.Properties.RowNames';
stats_struct(cnts).anova_preds_df = ran_out.DF';
stats_struct(cnts).anova_preds_stat = ran_out.F;
%- liberal
% tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
% stats_struct(cnts).anova_preds_p = tmp(1:2)';
%- conservative
tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
stats_struct(cnts).anova_preds_p = tmp(1:2)';
%-

%- model
stats_struct(cnts).mod_preds_coeff = rm_out.Coefficients; % remove from table before xlsx
stats_struct(cnts).mod_preds_terms = rm_out.BetweenFactorNames';
stats_struct(cnts).multi_comp_terms = ran_mult.Properties.VariableNames(1:2);
stats_struct(cnts).multi_comp_p = ran_mult.pValue;
stats_struct(cnts).multi_comp_t1_t2 = num2cell(ran_mult{:,1:2},2)';
stats_struct(cnts).multi_comp_lci_uci = [ran_mult.Lower,ran_mult.Upper];
stats_struct(cnts).norm_test_p = norm_p;
stats_struct(cnts).norm_test_h = norm_h;
cnts = cnts + 1;
stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
%## RANOVA VISUALIZE
% tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
comparisons = [1,2;1,3;2,3];
[h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(ran_mult.pValue);
adj_p = [adj_p(1),adj_p(2),adj_p(4)];
anova_p = 0;
% anova_grp_p = ran_out.pValue(2);
anova_grp_p = tmp(2);
% ran_means_1 = [rm_out.Coefficients{1,:}].*[0.25,0.5,0.75,1];
% ran_means_2 = [rm_out.Coefficients{2,:}].*[0.25,0.5,0.75,1]+ran_means_1;
% ran_means_3 = [rm_out.Coefficients{3,:}].*[0.25,0.5,0.75,1]+ran_means_1;
ran_means_1 = [rm_out.Coefficients{1,:}];
ran_means_2 = [rm_out.Coefficients{2,:}]+ran_means_1;
ran_means_3 = [rm_out.Coefficients{3,:}]+ran_means_1;
if anova_grp_p > 0.01 && anova_grp_p < 0.05
    str = {sprintf('* pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.01 && anova_grp_p > 0.001
    str = {sprintf('** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
elseif anova_grp_p <= 0.001
    str = {sprintf('*** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
else
    str = {'','',''};
end
STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(adj_p)},...
        'pvals_grp_pairs',{num2cell(comparisons,2)},...
        'regress_pval',{{0,0,0}},...
        'regress_line',{{ran_means_2,...
            ran_means_3,...
            ran_means_1}},...
        'line_type',{'means'},...
        'stats_char',{str},...
        'display_stats_char',true,...
        'bracket_yshift_perc',0.5,...
        'bracket_y_perc',0.5,...
        'bracket_rawshifty_upper',0,...
        'bracket_rawshifty_lower',0,...
        'grp_sig_offset_x',0,...
        'grp_sig_offset_y',0);
%## PLOT =============================================================== %%
inds = [stats_struct.mod_tag] == categorical(2) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
tmp_stats = stats_struct(inds);
VIOLIN_PARAMS = {'width',0.1,...
    'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
    'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
    'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
    'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
    'cond_labels',unique(tmp_table.cond_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
    'cond_offsets',[-0.35,-0.1,0.15,0.40],...
    'group_offsets',[0.125,0.475,0.812],...
    'y_label',meas_ylabel{meas_i},...
    'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
    'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
    'regresslab_txt_size',FONT_SIZE_VIO_REG);

%##
ax = axes();
ax = group_violin(tmp_table,meas_names{meas_i},'cond_char','group_id',...
    ax,...
    'VIOLIN_PARAMS',VIOLIN_PARAMS,...
    'PLOT_STRUCT',PLOT_STRUCT,...
    'STATS_STRUCT',STATS_STRUCT);
set(ax,'OuterPosition',[0,0,1,1]);
set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
%- set ylabel & title
% set(ax,'FontName','Arial','FontSize',XTICK_FONTSIZE,'FontWeight','normal')
% yticks(ax,'FontSize',XTICK_FONTSIZE)
yt = yticks(ax);
xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
pos1=get(xlh,'Position');
pos1(1,2)=pos1(1,2)+XLABEL_OFFSET %+0.17;
set(xlh,'Position',pos1);
ylabel(ax,meas_ylabel{meas_i},'FontSize',YLAB_FONTSIZE,'FontWeight',YLAB_FONTWEIGHT);
title(ax,meas_titles{meas_i},...
    'FontSize',TITLE_FONTSIZE,'FontWeight',TITLE_FONTWEIGHT);
%## LETTER
annotation('textbox',[AX_INIT_HORIZ+LAB_C_XOFFSET+(0.1/2)+horiz_shift,AX_INIT_VERT+vert_shift+LAB_B_YOFFSET+(0.1/2),.1,.1],...
    'String','D)','HorizontalAlignment','left',...
    'VerticalAlignment','top','LineStyle','none','FontName',FONT_NAME,...
    'FontSize',14,'FontWeight','Bold','Units','normalized');
%## SHIFT
% horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
%## SAVE
exportgraphics(fig,[save_dir filesep 'figure_behavioral_data_paper_ranova.tiff'],'Resolution',1000);

%%
save([save_dir filesep 'figure_behavioral_data_stats_ranova.mat'],'stats_struct')
stats_struct = struct2table(stats_struct);
stats_struct(isundefined(stats_struct.measure_tag),:) = [];
% stats_struct(:,'anova_preds_terms')
writetable(stats_struct,[save_dir filesep 'figure_behavioral_data_stats_ranova.xlsx'])