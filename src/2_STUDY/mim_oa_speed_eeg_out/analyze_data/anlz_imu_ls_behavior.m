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
% table_subj_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
table_header_names_ls = cell(length([SUBJ_PICS{:}]),1);
% table_subj_cat_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_cat_ls = categorical(repmat({''},length([SUBJ_PICS{:}])*N_TRIALS,1));
groupid_ls = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
designid_ls = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
fname_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
%- Loop through directory
table_imu_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
table_trial_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_imu_conds = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
groupid_imu_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
designid_imu_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
fname_imu_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
%-
table_ls_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
table_trial_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_ls_conds = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
groupid_ls_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
designid_ls_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,1);
fname_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
subj_stack = [];

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
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0.25';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0.5';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0.75';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1.0';
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
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0.25';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0.5';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0.75';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1.0';
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
            cnt_ls = cnt_ls + 1;
        end
        %{
        %## IMU TRIAL AVERAGE (SACRAL MEASURES ONLY) CONDITIONS
        for f_i = 1:length(dir_imu_conds)
            trial_i = [];
            cat_i = [];
            des_i = [];
            group_id = [];
            tmp = load([dir_imu_conds(f_i).folder filesep dir_imu_conds(f_i).name]);
            tmp = tmp.mergedOutputStruct;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_imu_conds{cnt_ls} = fieldnames(tmp);
            table_imu_meas_conds(cnt_ls,1:length(values)) = values;
            tmp = strsplit(dir_imu_conds(f_i).name,'.');
            fprintf(fid,'IMU Assigning %s.\n',dir_imu_conds(f_i).name);
            if any(strcmp(tmp{1},{'TM_med'}))
                trial_i = 'med';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low'}))
                trial_i = 'low';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat'}))
                trial_i = 'flat';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high'}))
                trial_i = 'high';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25'}))
                trial_i = '0p25';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5'}))
                trial_i = '0p5';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75'}))
                trial_i = '0p75';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0'}))
                trial_i = '1p0';
                des_i = 2;
            else
                trial_i = tmp{1};
                des_i = nan();
            end
            
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
                group_id = 1;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
                group_id = 1;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
                group_id = 1;
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
                group_id = 1;
            else
                cat_i = nan();
                group_id = nan();
            end
            table_trial_imu_conds{cnt_ls} = trial_i; %tmp{1};
            table_subj_imu_conds{cnt_ls} = SUBJ_PICS{group_i}{subj_i};
            table_subj_cat_imu_conds{cnt_ls} = cat_i;
            fname_imu_conds{cnt_ga} = dir_imu_conds(f_i).name;
            groupid_imu_conds(cnt_ga) = group_id;
            designid_imu_conds(cnt_ga) = des_i;
            cnt_ls = cnt_ls + 1;
        end
        %## LOADSOL TRIAL AVERAGES (GAIT)
        for f_i = 1:length(dir_ls_conds)
            trial_i = [];
            cat_i = [];
            des_i = [];
            group_id = [];
            tmp = load([dir_ls_conds(f_i).folder filesep dir_ls_conds(f_i).name]);
            tmp = tmp.mergedOutputStruct;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_ls_conds{cnt_ga} = fieldnames(tmp);
            table_ls_meas_conds(cnt_ga,1:length(values)) = values;
            tmp = strsplit(dir_ls_conds(f_i).name,'.');
            fprintf(fid,'GRF Assigning %s.\n',dir_ls_conds(f_i).name);
            %-
            if any(strcmp(tmp{1},{'TM_med'}))
                trial_i = 'med';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_low'}))
                trial_i = 'low';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_flat'}))
                trial_i = 'flat';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'TM_high'}))
                trial_i = 'high';
                des_i = 1;
            elseif any(strcmp(tmp{1},{'SP_0p25'}))
                trial_i = '0p25';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p5'}))
                trial_i = '0p5';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_0p75'}))
                trial_i = '0p75';
                des_i = 2;
            elseif any(strcmp(tmp{1},{'SP_1p0'}))
                trial_i = '1p0';
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
            table_subj_cat_ls_conds{cnt_ga} = cat_i;
            table_trial_ls_conds{cnt_ga} = trial_i; %tmp{1};
            table_subj_ls_conds{cnt_ga} = SUBJ_PICS{group_i}{subj_i};
            fname_ls_conds{cnt_ga} = dir_ls_conds(f_i).name;
            groupid_ls_conds(cnt_ga) = group_id;
            designid_ls_conds(cnt_ga) = des_i;
            cnt_ga = cnt_ga + 1;
        end
        %}
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
subj_cat_in = table_subj_cat_imu(rej);
fname_in = fname_imu(rej);
group_id_in  = groupid_imu(rej);
design_id_in = designid_imu(rej);
%-
table_imu_out = array2table(meas_in,'VariableNames',table_header_names_imu{1});
table_imu_out.subj_char = subj_in;
table_imu_out.trial_char = trial_in;
table_imu_out.subj_cat = subj_cat_in;
table_imu_out.file_name = fname_in;
table_imu_out.group_id = group_id_in;
table_imu_out.design_id = design_id_in;
writetable(table_imu_out,[save_dir filesep 'imu_table_out.xlsx']);
save([save_dir filesep 'imu_table_out.mat'],'table_imu_out');
%- rearrange headers
% table_imu_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
% table_imu_out = table(categorical(table_subj_vec),categorical(table_trial_vec),table_imu_meas(:,1),table_imu_meas(:,2),...
%     table_imu_meas(:,3),table_imu_meas(:,4),table_imu_meas(:,5),table_imu_meas(:,6),...
%     table_imu_meas(:,7),table_imu_meas(:,8),table_imu_meas(:,9),table_imu_meas(:,10),...
%     table_imu_meas(:,11),table_imu_meas(:,12),'VariableNames',[{'Subject'},{'TrialName'},table_header_names{1}']);
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
subj_cat_in = table_subj_cat_ls(rej);
fname_in = fname_ls(rej);
group_id_in  = groupid_ls(rej);
design_id_in = designid_ls(rej);
table_ls_out = array2table(meas_in,'VariableNames',table_header_names_ls{1});
table_ls_out.subj_char = subj_in;
table_ls_out.trial_char = trial_in;
table_ls_out.subj_cat = subj_cat_in;
table_ls_out.file_name = fname_in;
table_ls_out.group_id = group_id_in;
table_ls_out.design_id = design_id_in;
writetable(table_ls_out,[save_dir filesep 'ls_table_out.xlsx']);
save([save_dir filesep 'ls_table_out.mat'],'table_ls_out')
%- rearrange headers
% table_ls_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
%% IMU TABLE CONDITIONS (PRECALCULATED)
% rej = ~any(table_imu_meas_conds == 0,2);
% table_imu_meas_conds = table_imu_meas_conds(rej,:);
% table_subj_vec_conds = table_subj_vec_conds(rej);
% table_trial_vec_conds = table_trial_vec_conds(rej);
% table_subj_cat_vec_conds = table_subj_cat_vec_conds(rej);
% table_imu_out_conds = array2table(table_imu_meas_conds,'VariableNames',table_header_names_conds{1});
% table_imu_out_conds.SubjectName = categorical(table_subj_vec_conds);
% table_imu_out_conds.TrialName = categorical(table_trial_vec_conds);
% table_imu_out_conds.SubjectCategory = categorical(table_subj_cat_vec_conds);
% writetable(table_imu_out_conds,[save_dir filesep 'imu_table_out.xlsx']);
%- rearrange headers
% table_imu_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
% table_imu_out = table(categorical(table_subj_vec),categorical(table_trial_vec),table_imu_meas(:,1),table_imu_meas(:,2),...
%     table_imu_meas(:,3),table_imu_meas(:,4),table_imu_meas(:,5),table_imu_meas(:,6),...
%     table_imu_meas(:,7),table_imu_meas(:,8),table_imu_meas(:,9),table_imu_meas(:,10),...
%     table_imu_meas(:,11),table_imu_meas(:,12),'VariableNames',[{'Subject'},{'TrialName'},table_header_names{1}']);
%% LOADSOL TABLE CONDITIONS (PRECALCULATED)
% rej = ~any(table_ls_meas_conds == 0,2);
% table_ls_meas_conds = table_ls_meas_conds(rej,:);
% table_subj_ls_conds = table_subj_ls_conds(rej);
% table_trial_ls_conds = table_trial_ls_conds(rej);
% table_subj_cat_ls_conds = table_subj_cat_ls_conds(rej);
% table_ls_out_conds = array2table(table_ls_meas_conds,'VariableNames',table_header_names_ls_conds{1});
% table_ls_out_conds.SubjectName = categorical(table_subj_ls_conds);
% table_ls_out_conds.TrialName = categorical(table_trial_ls_conds);
% table_ls_out_conds.SubjectCategory = categorical(table_subj_cat_ls_conds);
% writetable(table_ls_out_conds,[save_dir filesep 'ls_table_out.xlsx']);
%- rearrange headers
% table_ls_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
%% LOAD IN IF AVAILABLE
table_ls_out = load([save_dir filesep 'ls_table_out.mat']);
table_ls_out = table_ls_out.table_ls_out;
table_imu_out = load([save_dir filesep 'imu_table_out.mat']);
table_imu_out = table_imu_out.table_imu_out;
%% READ IN SUBJECT STABILITY SCORES (TRIAL)
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\subject_mastersheet_notes.xlsx');
tmp = fieldnames(MasterTable(:,22:38));
tmp = tmp(1:end-3);
tmp_table = MasterTable(:,22:38);
og_length = size(tmp_table,2);
tmp_table.subject_code = categorical(MasterTable.subject_code);
tmp_table.stability_score = MasterTable.flat_low_med_high_rating_of_stability;
% table_new_imu.stability_rating = zeros(size(table_new_imu,1),1);
% table_new_ls.stability_rating = zeros(size(table_new_ls,1),1);
tmp_table_imu = table_imu_out;
tmp_table_ls = table_ls_out;
trials = {'flat','low','med','high'};
for i = 1:size(tmp_table,1)
    %- get subject index
    ss = tmp_table.subject_code(i);
    ss_var = tmp_table.stability_score(i);
    ss_var = strsplit(ss_var{1},';');
    ss_var = ss_var(~cellfun(@isempty,ss_var));
%     ss = table_new_imu.SubjectName(i);
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
            for f = 1:length(fn)-6
                tmp_table_imu(ind11,fn{f}) = tmp_table(tmp_table.subject_code==ss,f);
                tmp_table_ls(ind22,fn{f}) = tmp_table(tmp_table.subject_code==ss,f);
            end
            tmp_table_imu(ind11,'og_length_m') = {3};
            tmp_table_ls(ind22,'og_length_m') = {3};
        end
    end
end
fn = fieldnames(tmp_table);
for f = 1:length(fn)-6
    try
        tmp_table_imu(tmp_table_imu.(fn{f}) == 0,fn{f}) = {nan()};
        tmp_table_ls(tmp_table_ls.(fn{f}) == 0,fn{f}) = {nan()};
    catch ME
       switch ME.identifier
           case 'MATLAB:datetime:InvalidComparison'
               fprintf('Datetime value detected. Keeping as default.\n');
       end
    end
end
%%
tmp_table_imu.stability_rating(tmp_table_imu.stability_rating==0) = nan();
tmp_table_ls.stability_rating(tmp_table_ls.stability_rating==0) = nan();
tmp_table_imu.og_length_m(tmp_table_imu.og_length_m==0) = nan();
tmp_table_ls.og_length_m(tmp_table_ls.og_length_m==0) = nan();
writetable(tmp_table_imu,[save_dir filesep 'imu_table_trial.xlsx']);
writetable(tmp_table_ls,[save_dir filesep 'ls_table_trial.xlsx']);
save([save_dir filesep 'imu_table_trial.mat'],'tmp_table_imu')
save([save_dir filesep 'ls_table_trial.mat'],'tmp_table_ls')
table_ls_out = tmp_table_ls;
table_imu_out = tmp_table_imu;
%% LOAD IN IF AVAILABLE
table_ls_out = load([save_dir filesep 'ls_table_out.mat']);
table_ls_out = table_ls_out.table_ls_out;
table_imu_out = load([save_dir filesep 'imu_table_out.mat']);
table_imu_out = table_imu_out.table_imu_out;
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
            tmpvalsnot = varfun(@(x) x(1), tmp, 'InputVariables', @(x) ~isnumeric(x));
            tmpnames = cellfun(@(x) strsplit(x,'_'),tmpvalsnot.Properties.VariableNames,'UniformOutput',false);
            tmpnames = cellfun(@(x) strjoin(x(2:end),'_'),tmpnames,'UniformOutput',false);
            tmpvalsnot.Properties.VariableNames = tmpnames;
            tmpvals = cat(2,tmpvalsnot,tmpvals);
            % tmpvals.subj_char = tmp.subj_char(1);
            % tmpvals.trial_char = tmp.trial_char(1);
            % tmpvals.subj_cat = tmp.subj_cat(1);
            % tmpvals.file_name = tmp.file_name(1);
            % tmpvals.group_id = tmp.group_id(1);
            % tmpvals.design_id = tmp.design_id(1);
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
            tmpvalsnot = varfun(@(x) x(1), tmp, 'InputVariables', @(x) ~isnumeric(x));
            tmpnames = cellfun(@(x) strsplit(x,'_'),tmpvalsnot.Properties.VariableNames,'UniformOutput',false);
            tmpnames = cellfun(@(x) strjoin(x(2:end),'_'),tmpnames,'UniformOutput',false);
            tmpvalsnot.Properties.VariableNames = tmpnames;
            tmpvals = cat(2,tmpvalsnot,tmpvals);
            % tmpvals.subj_char = tmp.subj_char(1);
            % tmpvals.trial_char = tmp.trial_char(1);
            % tmpvals.subj_cat = tmp.subj_cat(1);
            % tmpvals.file_name = tmp.file_name(1);
            % tmpvals.group_id = tmp.group_id(1);
            % tmpvals.design_id = tmp.design_id(1);
            table_new = [table_new; tmpvals];
        end
    end
end
table_new_ls = tabl
%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
% SPEED_CUTOFF = 0.1;
SPEED_CUTOFF = 0;
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
table_new_imu.terrain_speed = zeros(size(table_new_imu,1),1);
table_new_ls.terrain_speed = zeros(size(table_new_ls,1),1);
for i = 1:size(speed_table,1)
    ss = speed_table.Var1(i);
    ss_speed = speed_table.Var2(i);
%     ss = table_new_imu.SubjectName(i);
    ind1 = table_new_imu.subj_char==ss;
    ind2 = table_new_ls.subj_char==ss;
    if ss_speed < SPEED_CUTOFF
%         inds_del = table_new_imu.SubjectName == ss;
%         table_new_imu(inds_del) = [];
        table_new_imu(ind1,:) = [];
        table_new_ls(ind2,:) = [];
    else
        table_new_imu.terrain_speed(ind1) = ss_speed;
        table_new_ls.terrain_speed(ind2) = ss_speed;
    end
end
writetable(table_new_imu,[save_dir filesep 'imu_table_meantrial.xlsx']);
writetable(table_new_ls,[save_dir filesep 'ls_table_meantrial.xlsx']);
save([save_dir filesep 'ls_table_meantrial.mat'],'table_new_ls')
save([save_dir filesep 'imu_table_meantrial.mat'],'table_new_imu')
%% READ IN SUBJECT STABILITY SCORES
%{
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\subject_mastersheet_notes.xlsx');
tmp_table = table(categorical(MasterTable.subject_code),MasterTable.flat_low_med_high_rating_of_stability);
% table_new_imu.stability_rating = zeros(size(table_new_imu,1),1);
% table_new_ls.stability_rating = zeros(size(table_new_ls,1),1);
trials = {'flat','low','med','high'};
for i = 1:size(tmp_table,1)
    ss = tmp_table.Var1(i);
    ss_var = tmp_table.Var2(i);
    ss_var = strsplit(ss_var{1},';');
    ss_var = ss_var(~cellfun(@isempty,ss_var));
%     ss = table_new_imu.SubjectName(i);
    ind1 = table_new_imu.SubjectName==ss;
    ind2 = table_new_ls.SubjectName==ss;
    if any(ind1) & ~isempty(ss_var)
        for j = 1:length(trials)
            sub = strsplit(ss_var{j},',');
            ind11 = table_new_imu.SubjectName==ss & table_new_imu.TrialName==trials{j};
            ind22 = table_new_ls.SubjectName==ss & table_new_ls.TrialName==trials{j};
            if ~isempty(sub{1})
                for k = 1:length(sub)
                    if any(ind11)
                        table_new_imu.(sprintf('%s_%i','stability_rating',k))(ind11) = double(string(sub{k}));
                        table_new_ls.(sprintf('%s_%i','stability_rating',k))(ind22) = double(string(sub{k}));
                    else
                        table_new_imu.(sprintf('%s_%i','stability_rating',k))(ind11) = nan();
                        table_new_ls.(sprintf('%s_%i','stability_rating',k))(ind22) = nan();
                    end
                end
            else
                for k = 1:2
                    table_new_imu.(sprintf('%s_%i','stability_rating',k))(ind11) = nan();
                    table_new_ls.(sprintf('%s_%i','stability_rating',k))(ind22) = nan();
                end
            end
        end
    end
end
writetable(table_new_imu,[save_dir filesep 'imu_table_meantrial.xlsx']);
writetable(table_new_ls,[save_dir filesep 'ls_table_meantrial.xlsx']);
%}
%%
%{
table_new_imu = readtable([save_dir filesep 'imu_table_meantrial.xlsx']);
table_new_ls = readtable([save_dir filesep 'ls_table_meantrial.xlsx']);
%}
%% VIOLIN PLOT IMU
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,420,420];
FIG_POSITION_GRP = [100,100,720,420];
VIOLIN_WIDTH = 0.3;
VIOLIN_WIDTH_GROUP = 0.1;
%-
% meas_names = {'nanmean_APexc_mean','nanmean_MLexc_mean','nanmean_APexc_COV','nanmean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
% meas_units = {'m','m','%','%'};
% meas_titles = {'Anteriorposterior Excursion Mean','Mediolateral Excursion Mean','Anteroposterior Excursion Coefficient of Variation','Mediolateral Excursion Coefficient of Variation'};
% meas_ylabel = {'Distance','Distance','Coefficient of Variation','Coefficient of Variation'};
% YLIMS = {[0,0.15],[0,0.3],[0,70],[0,47.5]};
% meas_names = {'nanmean_APexc_COV','nanmean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
meas_names = {'mean_APexc_COV','mean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
meas_units = {'%','%'};
meas_titles = {{'Anteroposterior Excursion';'Coefficient of Variation'},{'Mediolateral Excursion';'Coefficient of Variation'}};
% meas_titles = {{'Anteroposterior Excursion\nCoefficient of Variation'},{'Mediolateral Excursion\nCoefficient of Variation'}};

meas_ylabel = {'Coefficient of Variation (%)','Coefficient of Variation (%)'};
% YLIMS = {[0,90],[0,95]};
YLIMS = {[0,70],[0,47.5]};
%-
speed_chars = {'0p25','0p5','0p75','1p0'};
terrain_chars = {'flat','low','med','high'};
trial_names_speed = {'0.25','0.50','0.75',...
               '1.00'};
trial_names_terrain = {'flat',...
    'low','med.',...
    'high'};
if length(SUBJ_PICS) == 3
    g_cats = categorical({'YoungAdult';'HF_OlderAdult';'LF_OlderAdult'});
    group_names = {'Younger Adults','Older Adults\newlineHigh Function','Older Adults\newlineLow Function'};
elseif length(SUBJ_PICS) == 2
    g_cats = categorical({'HF_OlderAdult';'LF_OlderAdult'});
    group_names = {'Older Adults\newlineHigh Function','Older Adults\newlineLow Function'};
end
%- colors
COLORS_MAPS_TERRAIN = linspecer(size(terrain_chars,2));
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(size(speed_chars,2)*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%%
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'measure_tag',categorical({''}),...
    'design_tag',categorical({''}),...
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
%## PARAMS
IM_RESIZE = 1;
AX_H  = 0.2;
AX_W = 0.275;
AX_HORIZ_SHIFT = 0.06;
AX_INIT_VERT_VIO = 0.6;
AXES_FONT_SIZE_VIO = 10;
GROUP_LAB_FONTSIZE = AXES_FONT_SIZE_VIO;
GROUP_LAB_FONTWEIGHT = 'bold ';
XLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
YLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
XTICK_FONTSIZE = AXES_FONT_SIZE_VIO;
XLAB_FONTWEIGHT = 'bold';
YLAB_FONTWEIGHT = 'bold';
TITLE_FONTSIZE = AXES_FONT_SIZE_VIO;
TITLE_FONTWEIGHT = 'bold';
XLABEL_OFFSET = -.05;
GROUP_LAB_YOFFSET = -0.275;
%%
des_i = 2;
cnts = 1;
for meas_i = 1:length(meas_names)
    %%
    inds = table_new_imu.design_id == des_i;
    tmp_table = table_new_imu(inds,:);
    tmp_table = table(double(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
        categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
        'VariableNames',{'trial_char','subj_char','group_id',meas_names{meas_i}});
    %% (SPEED) STATISTICS GROUP & CONDITION EFFECT
    %## MAIN EFFECT GROUP & CONDITION + INTERACTION
    mod_out = sprintf('%s~1+group_id+trial_char+trial_char:group_id',meas_names{meas_i});
    stats_out = fitlm(tmp_table,mod_out);
    pred_terms = stats_out.CoefficientNames;
    % anova_out = anova(stats_out);
    [p,t,anova_out,terms] = anovan(tmp_table.(meas_names{meas_i}),{tmp_table.trial_char, tmp_table.group_id},...
        'sstype',3,'varnames',{'trial_char','group_id'},'model','interaction');
    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[1,2],...
        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
    disp(stats_out)
    %- test normality
    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
    %- get effects
    % [~,bnames,~] = stats_out.fixedEffects();
    % [~,brnames,bretable] = stats_out.randomEffects();
    %- intercept only model
    altmod_out = sprintf('%s ~ 1',meas_names{meas_i});
    altstats_out = fitlm(tmp_table,altmod_out);
    %- cohens f2
    % simp_f2 = (1-stats_out.Rsquared.Adjusted)/stats_out.Rsquared.Adjusted;
    %- see. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/ , eq. 3
    % vfull = var(stats_out.Residuals.Pearson);
    % vnull = var(altstats_out.Residuals.Pearson);
    % resd_f2 = (vnull-vfull)/vnull;
    % resd_f2 = (1-resd_f2)/resd_f2;
    %- alternative f2?
    R21 = altstats_out.SSR/altstats_out.SST;
    R22 = stats_out.SSR/stats_out.SST;
    alt_f2 = (R22-R21)/(1-R22);
    %- Convert summary to char array
    % txt = evalc('stats_out');
    % fid = fopen([save_dir filesep sprintf('%s_mdl_speed_group_int.txt',meas_names{meas_i})],'wt');
    % fprintf(fid,'%s',txt);
    % for i = 2:length(pred_terms)
    %     fprintf(fid,'%i: %s\n',i,pred_terms{i});
    % end
    %- populate struct
    stats_struct(cnts).stat_test_mod = mod_out;
    stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
    stats_struct(cnts).design_tag = categorical(des_i);
    stats_struct(cnts).mod_tag = categorical(1);
    % stats_struct(cnts).resp_terms = meas_names(meas_i);
    stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
    stats_struct(cnts).anova_preds_terms = t(:,1);
    stats_struct(cnts).anova_preds_p = t(:,7);
    stats_struct(cnts).anova_preds_stat = t(:,6);
    stats_struct(cnts).anova_preds_df = t(:,3);
    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Properties.RowNames;
    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
    stats_struct(cnts).multi_comp_terms = gnames;
    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
    stats_struct(cnts).multi_comp_p = comparisons(:,6);
    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
    stats_struct(cnts).norm_test_p = norm_p;
    stats_struct(cnts).norm_test_h = norm_h;
    stats_struct(cnts).effect_size = alt_f2;
    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
    cnts = cnts + 1;
    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;

    %% MAIN EFFECT GROUP & CONDITION (ASSUMING INTERACTIONS ARE NOT SIGNIFICANT
    mod_out = sprintf('%s~1+group_id+trial_char',meas_names{meas_i});
    stats_out = fitlm(tmp_table,mod_out);
    pred_terms = stats_out.CoefficientNames;
    % anova_out = anova(stats_out);
    [p,t,anova_out,terms] = anovan(tmp_table.(meas_names{meas_i}),{tmp_table.trial_char, tmp_table.group_id},...
        'sstype',3,'varnames',{'trial_char','group_id'},'model','linear');
    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[1,2],...
        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
    disp(stats_out)
    %- test normality
    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
    %- get effects
    % [~,bnames,~] = stats_out.fixedEffects();
    % [~,brnames,bretable] = stats_out.randomEffects();
    %- intercept only model
    altmod_out = sprintf('%s ~ 1',meas_names{meas_i});
    altstats_out = fitlm(tmp_table,altmod_out);
    %- cohens f2, this is wrong because it uses correlation coefficient (I
    %think)?
    % simp_f2 = (1-stats_out.Rsquared.Adjusted)/stats_out.Rsquared.Adjusted;
    %- see. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/ , eq. 3
    % vfull = var(stats_out.Residuals.Pearson);
    % vnull = var(altstats_out.Residuals.Pearson);
    % resd_f2 = (vnull-vfull)/vnull;
    % resd_f2 = (1-resd_f2)/resd_f2;
    %- alternative f2?
    R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
    R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
    alt_f2 = (R22-R21)/(1-R22);
    %- Convert summary to char array
    % txt = evalc('stats_out');
    % fid = fopen([save_dir filesep sprintf('%s_mdl_speed_group.txt',meas_names{meas_i})],'wt');
    % fprintf(fid,'%s',txt);
    % for i = 2:length(pred_terms)
    %     fprintf(fid,'%i: %s\n',i,pred_terms{i});
    % end
    % fclose(fid);
    %- populate struct
    stats_struct(cnts).stat_test_mod = mod_out;
    stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
    stats_struct(cnts).design_tag = categorical(des_i);
    stats_struct(cnts).mod_tag = categorical(2);
    % stats_struct(cnts).resp_terms = meas_names(meas_i);
    stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
    stats_struct(cnts).anova_preds_terms = t(:,1);
    stats_struct(cnts).anova_preds_p = t(:,7);
    stats_struct(cnts).anova_preds_stat = t(:,6);
    stats_struct(cnts).anova_preds_df = t(:,3);
    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Properties.RowNames;
    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
    stats_struct(cnts).multi_comp_terms = gnames;
    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
    stats_struct(cnts).multi_comp_p = comparisons(:,6);
    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
    stats_struct(cnts).norm_test_p = norm_p;
    stats_struct(cnts).norm_test_h = norm_h;
    stats_struct(cnts).effect_size = alt_f2;
    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
    cnts = cnts + 1;
    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
    close('all');
    %%
    horiz_shift = 0;
    vert_shift = 0;
    FIGURE_POSITION =[0,0,6.5,9];
    FONT_NAME = 'Arial';
    FONT_SIZE_VIO = 10;
    FONT_SIZE_VIO_REG = 7;
    % IM_RESIZE = 1;
    % AX_H  = 0.2;
    % AX_W = 0.275;
    AX_HORIZ_SHIFT = 0.06;
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
    %-
    inds = [stats_struct.mod_tag] == categorical(1) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
    tmp_stats = stats_struct(inds);
    VIOLIN_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
    PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
        'cond_labels',unique(tmp_table.trial_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
        'cond_offsets',[-0.35,-0.1,0.15,0.40],...
        'group_offsets',[0.125,0.475,0.812],...
        'y_label',meas_ylabel{meas_i},...
        'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
        'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
        'regresslab_txt_size',FONT_SIZE_VIO_REG);
    %-
    cond_ind = 2;
    anova_p = tmp_stats.anova_preds_p{cond_ind};
    grp_ind = 3;
    anova_grp_p = tmp_stats.anova_preds_p{grp_ind};
    intact_ind = 4;
    anova_intact_p = tmp_stats.anova_preds_p{intact_ind};
    speed_p = tmp_stats.mod_preds_p([2,5:6]);
    grp_p = tmp_stats.mod_preds_p([2,5:6]);
    speed_r2 = tmp_stats.mod_r2;
    speed_xvals = (0:5)*0.25;
    STATS_STRUCT = struct('anova',{{anova_intact_p,anova_intact_p,anova_intact_p}},...
        'anova_grp',{{anova_intact_p,anova_intact_p,anova_intact_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(grp_p)},...
        'pvals_grp_pairs',{{[1,1],[1,2],[1,3]}},...
        'regress_pval',{num2cell(speed_p)},...
        'regress_line',{{[tmp_stats.mod_preds_coeff(1), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(3), tmp_stats.mod_preds_coeff(2)+tmp_stats.mod_preds_coeff(5)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(4), tmp_stats.mod_preds_coeff(2)+tmp_stats.mod_preds_coeff(6)]}},...
        'r2_coeff',{repmat(speed_r2,3,1)},...
        'regress_xvals',speed_xvals,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'do_include_intercept',true);
    %-
    fig = figure('color','white','renderer','Painters');
    
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    set(fig,'Units','inches','Position',FIGURE_POSITION)
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:});
    ax = axes();
    ax = group_violin(tmp_table,meas_names{meas_i},'trial_char','group_id',...
        ax,...
        'VIOLIN_PARAMS',VIOLIN_PARAMS,...
        'PLOT_STRUCT',PLOT_STRUCT,...
        'STATS_STRUCT',STATS_STRUCT);
    set(ax,'OuterPosition',[0,0,1,1]);
    set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]    
    % ax = gca;
    ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
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
    exportgraphics(fig,[save_dir filesep sprintf('speed_groupwise_intact_%s.tiff',meas_names{meas_i})],'Resolution',1000);
    %%
    horiz_shift = 0;
    vert_shift = 0;
    FIGURE_POSITION =[0,0,6.5,9];
    FONT_NAME = 'Arial';
    FONT_SIZE_VIO = 10;
    FONT_SIZE_VIO_REG = 7;
    % IM_RESIZE = 1;
    % AX_H  = 0.2;
    % AX_W = 0.275;
    % AX_HORIZ_SHIFT = 0.06;
    % AX_INIT_VERT_VIO = 0.6;
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
    %-
    inds = [stats_struct.mod_tag] == categorical(2) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
    tmp_stats = stats_struct(inds);
    VIOLIN_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
    PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
        'cond_labels',unique(tmp_table.trial_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
        'cond_offsets',[-0.35,-0.1,0.15,0.40],...
        'group_offsets',[0.125,0.475,0.812],...
        'y_label',meas_ylabel{meas_i},...
        'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
        'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
        'regresslab_txt_size',FONT_SIZE_VIO_REG);
    %-
    cond_ind = 2;
    anova_p = tmp_stats.anova_preds_p{cond_ind};
    grp_ind = 3;
    anova_grp_p = tmp_stats.anova_preds_p{grp_ind};
    speed_p = tmp_stats.mod_preds_p(2);
    grp_p = [0,tmp_stats.mod_preds_p(3:4)']';
    speed_r2 = tmp_stats.mod_r2;
    speed_xvals = (0:5)*0.25;
    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(grp_p)},...
        'pvals_grp_pairs',{{[1,1],[1,2],[1,3]}},...
        'regress_pval',{{speed_p,speed_p,speed_p}},...
        'regress_line',{{[tmp_stats.mod_preds_coeff(1), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(3), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(4), tmp_stats.mod_preds_coeff(2)]}},...
        'r2_coeff',{repmat(speed_r2,3,1)},...
        'regress_xvals',speed_xvals,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'do_include_intercept',false);
    %-
    fig = figure('color','white','renderer','Painters');
    % annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
    %     'String',atlas_name,'HorizontalAlignment','center',...
    %     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
    %     'FontSize',14,'FontWeight','Bold','Units','normalized');
    set(fig,'Units','inches','Position',FIGURE_POSITION)
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:});
    ax = axes();
    ax = group_violin(tmp_table,meas_names{meas_i},'trial_char','group_id',...
        ax,...
        'VIOLIN_PARAMS',VIOLIN_PARAMS,...
        'PLOT_STRUCT',PLOT_STRUCT,...
        'STATS_STRUCT',STATS_STRUCT);
    set(ax,'OuterPosition',[0,0,1,1]);
    set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
    %- set ylabel & title
    AXES_FONT_SIZE_VIO = 10;
    GROUP_LAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    GROUP_LAB_FONTWEIGHT = 'bold ';
    XLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    YLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    XTICK_FONTSIZE = AXES_FONT_SIZE_VIO;
    XLAB_FONTWEIGHT = 'bold';
    YLAB_FONTWEIGHT = 'bold';
    TITLE_FONTSIZE = AXES_FONT_SIZE_VIO;
    TITLE_FONTWEIGHT = 'bold';
    XLABEL_OFFSET = -.05;
    GROUP_LAB_YOFFSET = -0.275;
    % ax = gca;
    ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
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
    exportgraphics(fig,[save_dir filesep sprintf('speed_groupwise_nointact_%s.tiff',meas_names{meas_i})],'Resolution',1000);
end
%% VIOLIN PLOT LOADSOL
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,360,420];
FIG_POSITION_GRP = [100,100,720,420];
VIOLIN_WIDTH = 0.3;
VIOLIN_WIDTH_GROUP = 0.1;
%-
% meas_names = {'nanmean_StepDur','nanmean_StepDur_cov','nanmean_GaitCycleDur_cov','nanmean_GaitCycleDur',};
% meas_units = {'s','%','%','s'};
% meas_titles = {'Step Duration','Step Duration Coefficient of Variation','Gait Cycle Duration Coefficient of Variation','Gait Cycle Duration'};
% meas_ylabel = {'Duration','Coefficient of Variation','Coefficient of Variation','Duration'};
% YLIMS = {[0,2],[0,27.5],[0,30],[0,4]};
% meas_names = {'nanmean_StepDur','nanmean_StepDur_cov'};
% meas_names = {'mean_StepDur','mean_StepDur_cov','mean_StanceDur','mean_SwingDur'};
meas_names = {'mean_StepDur','mean_StepDur_cov','mean_GaitCycleDur','mean_SwingDur',...
            'mean_StanceDur','mean_SingleSupport','mean_TotalDS'};
% meas_units = {'s','s','s','s','s','s'};
meas_titles = {'Step Duration','Step Duration COV','Gait Cycle Duration','Swing Duration',...
            'Stance Duration','Single Support Time','Totals Double Support Time'};
% meas_titles = {'Step Duration',{'Step Duration';'Coefficient of Variation'}};
% meas_ylabel = {'Duration','Coefficient of Variation'};
meas_ylabel = {'Duration (s)','Variation (%)','Duration (s)','Duration (s)',...
    'Duration (s)','Duration (s)','Duration (s)'};
YLIMS = {[0,2.5],[0,80],[0,3],[0,1],[0,3],[0,1.5],[0,2]};
%-
speed_chars = {'0p25','0p5','0p75','1p0'};
terrain_chars = {'flat','low','med','high'};
trial_names_speed = {'0.25','0.50','0.75',...
               '1.00'};
trial_names_terrain = {'flat',...
    'low','med.',...
    'high'};
if length(SUBJ_PICS) == 3
    g_cats = categorical({'YoungAdult';'HF_OlderAdult';'LF_OlderAdult'});
    group_names = {'Younger Adults','Older Adults\newlineHigh Function','Older Adults\newlineLow Function'};
elseif length(SUBJ_PICS) == 2
    g_cats = categorical({'HF_OlderAdult';'LF_OlderAdult'});
    group_names = {'Older Adults\newlineHigh Function','Older Adults\newlineLow Function'};
end
%-
COLORS_MAPS_TERRAIN = linspecer(size(speed_chars,2));
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(size(speed_chars,2)*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%%
des_i = 2;
% cnts = 1;
for meas_i = 1:length(meas_names)
    %%
    inds = table_new_ls.design_id == des_i;
    tmp_table = table_new_ls(inds,:);
    tmp_table = table(double(string(tmp_table.trial_char)),categorical(string(tmp_table.subj_char)),...
        categorical(string(tmp_table.group_id)),tmp_table.(meas_names{meas_i}),...
        'VariableNames',{'trial_char','subj_char','group_id',meas_names{meas_i}});
    %% (SPEED) STATISTICS GROUP & CONDITION EFFECT
    %## MAIN EFFECT GROUP & CONDITION + INTERACTION
    mod_out = sprintf('%s~1+group_id+trial_char+trial_char:group_id',meas_names{meas_i});
    stats_out = fitlm(tmp_table,mod_out);
    pred_terms = stats_out.CoefficientNames;
    % anova_out = anova(stats_out);
    [p,t,anova_out,terms] = anovan(tmp_table.(meas_names{meas_i}),{tmp_table.trial_char, tmp_table.group_id},...
        'sstype',3,'varnames',{'trial_char','group_id'},'model','interaction');
    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[1,2],...
        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
    disp(stats_out)
    %- test normality
    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
    %- get effects
    % [~,bnames,~] = stats_out.fixedEffects();
    % [~,brnames,bretable] = stats_out.randomEffects();
    %- intercept only model
    altmod_out = sprintf('%s ~ 1',meas_names{meas_i});
    altstats_out = fitlm(tmp_table,altmod_out);
    %- cohens f2
    % simp_f2 = (1-stats_out.Rsquared.Adjusted)/stats_out.Rsquared.Adjusted;
    %- see. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/ , eq. 3
    % vfull = var(stats_out.Residuals.Pearson);
    % vnull = var(altstats_out.Residuals.Pearson);
    % resd_f2 = (vnull-vfull)/vnull;
    % resd_f2 = (1-resd_f2)/resd_f2;
    %- alternative f2?
    R21 = altstats_out.SSR/altstats_out.SST;
    R22 = stats_out.SSR/stats_out.SST;
    alt_f2 = (R22-R21)/(1-R22);
    %- Convert summary to char array
    % txt = evalc('stats_out');
    % fid = fopen([save_dir filesep sprintf('%s_mdl_speed_group_int.txt',meas_names{meas_i})],'wt');
    % fprintf(fid,'%s',txt);
    % for i = 2:length(pred_terms)
    %     fprintf(fid,'%i: %s\n',i,pred_terms{i});
    % end
    %- populate struct
    stats_struct(cnts).stat_test_mod = mod_out;
    stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
    stats_struct(cnts).design_tag = categorical(des_i);
    stats_struct(cnts).mod_tag = categorical(1);
    % stats_struct(cnts).resp_terms = meas_names(meas_i);
    stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
    stats_struct(cnts).anova_preds_terms = t(:,1);
    stats_struct(cnts).anova_preds_p = t(:,7);
    stats_struct(cnts).anova_preds_stat = t(:,6);
    stats_struct(cnts).anova_preds_df = t(:,3);
    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Properties.RowNames;
    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
    stats_struct(cnts).multi_comp_terms = gnames;
    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
    stats_struct(cnts).multi_comp_p = comparisons(:,6);
    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
    stats_struct(cnts).norm_test_p = norm_p;
    stats_struct(cnts).norm_test_h = norm_h;
    stats_struct(cnts).effect_size = alt_f2;
    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
    cnts = cnts + 1;
    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;

    %% MAIN EFFECT GROUP & CONDITION (ASSUMING INTERACTIONS ARE NOT SIGNIFICANT
    mod_out = sprintf('%s~1+group_id+trial_char',meas_names{meas_i});
    stats_out = fitlm(tmp_table,mod_out);
    pred_terms = stats_out.CoefficientNames;
    % anova_out = anova(stats_out);
    [p,t,anova_out,terms] = anovan(tmp_table.(meas_names{meas_i}),{tmp_table.trial_char, tmp_table.group_id},...
        'sstype',3,'varnames',{'trial_char','group_id'},'model','linear');
    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[1,2],...
        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
    disp(stats_out)
    %- test normality
    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
    %- get effects
    % [~,bnames,~] = stats_out.fixedEffects();
    % [~,brnames,bretable] = stats_out.randomEffects();
    %- intercept only model
    altmod_out = sprintf('%s ~ 1',meas_names{meas_i});
    altstats_out = fitlm(tmp_table,altmod_out);
    %- cohens f2, this is wrong because it uses correlation coefficient (I
    %think)?
    % simp_f2 = (1-stats_out.Rsquared.Adjusted)/stats_out.Rsquared.Adjusted;
    %- see. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/ , eq. 3
    % vfull = var(stats_out.Residuals.Pearson);
    % vnull = var(altstats_out.Residuals.Pearson);
    % resd_f2 = (vnull-vfull)/vnull;
    % resd_f2 = (1-resd_f2)/resd_f2;
    %- alternative f2?
    R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
    R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
    alt_f2 = (R22-R21)/(1-R22);
    %- Convert summary to char array
    % txt = evalc('stats_out');
    % fid = fopen([save_dir filesep sprintf('%s_mdl_speed_group.txt',meas_names{meas_i})],'wt');
    % fprintf(fid,'%s',txt);
    % for i = 2:length(pred_terms)
    %     fprintf(fid,'%i: %s\n',i,pred_terms{i});
    % end
    % fclose(fid);
    %- populate struct
    stats_struct(cnts).stat_test_mod = mod_out;
    stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
    stats_struct(cnts).design_tag = categorical(des_i);
    stats_struct(cnts).mod_tag = categorical(2);
    % stats_struct(cnts).resp_terms = meas_names(meas_i);
    stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
    stats_struct(cnts).anova_preds_terms = t(:,1);
    stats_struct(cnts).anova_preds_p = t(:,7);
    stats_struct(cnts).anova_preds_stat = t(:,6);
    stats_struct(cnts).anova_preds_df = t(:,3);
    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Properties.RowNames;
    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
    stats_struct(cnts).multi_comp_terms = gnames;
    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
    stats_struct(cnts).multi_comp_p = comparisons(:,6);
    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
    stats_struct(cnts).norm_test_p = norm_p;
    stats_struct(cnts).norm_test_h = norm_h;
    stats_struct(cnts).effect_size = alt_f2;
    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
    cnts = cnts + 1;
    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
    close('all');
    %%
    horiz_shift = 0;
    vert_shift = 0;
    FIGURE_POSITION =[0,0,6.5,9];
    FONT_NAME = 'Arial';
    FONT_SIZE_VIO = 10;
    FONT_SIZE_VIO_REG = 7;
    % IM_RESIZE = 1;
    % AX_H  = 0.2;
    % AX_W = 0.275;
    AX_HORIZ_SHIFT = 0.06;
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
    %-
    inds = [stats_struct.mod_tag] == categorical(1) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
    tmp_stats = stats_struct(inds);
    VIOLIN_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
    PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
        'cond_labels',unique(tmp_table.trial_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
        'cond_offsets',[-0.35,-0.1,0.15,0.40],...
        'group_offsets',[0.125,0.475,0.812],...
        'y_label',meas_ylabel{meas_i},...
        'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
        'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
        'regresslab_txt_size',FONT_SIZE_VIO_REG);
    %-
    cond_ind = 2;
    anova_p = tmp_stats.anova_preds_p{cond_ind};
    grp_ind = 3;
    anova_grp_p = tmp_stats.anova_preds_p{grp_ind};
    intact_ind = 4;
    anova_intact_p = tmp_stats.anova_preds_p{intact_ind};
    speed_p = tmp_stats.mod_preds_p([2,5:6]);
    grp_p = tmp_stats.mod_preds_p([2,5:6]);
    speed_r2 = tmp_stats.mod_r2;
    speed_xvals = (0:5)*0.25;
    STATS_STRUCT = struct('anova',{{anova_intact_p,anova_intact_p,anova_intact_p}},...
        'anova_grp',{{anova_intact_p,anova_intact_p,anova_intact_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(grp_p)},...
        'pvals_grp_pairs',{{[1,1],[1,2],[1,3]}},...
        'regress_pval',{num2cell(speed_p)},...
        'regress_line',{{[tmp_stats.mod_preds_coeff(1), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(3), tmp_stats.mod_preds_coeff(2)+tmp_stats.mod_preds_coeff(5)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(4), tmp_stats.mod_preds_coeff(2)+tmp_stats.mod_preds_coeff(6)]}},...
        'r2_coeff',{repmat(speed_r2,3,1)},...
        'regress_xvals',speed_xvals,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'do_include_intercept',true);
    %-
    fig = figure('color','white','renderer','Painters');
    
    annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
        'String',atlas_name,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
        'FontSize',14,'FontWeight','Bold','Units','normalized');
    set(fig,'Units','inches','Position',FIGURE_POSITION)
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:});
    ax = axes();
    ax = group_violin(tmp_table,meas_names{meas_i},'trial_char','group_id',...
        ax,...
        'VIOLIN_PARAMS',VIOLIN_PARAMS,...
        'PLOT_STRUCT',PLOT_STRUCT,...
        'STATS_STRUCT',STATS_STRUCT);
    set(ax,'OuterPosition',[0,0,1,1]);
    set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
    %- set ylabel & title
    AXES_FONT_SIZE_VIO = 10;
    GROUP_LAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    GROUP_LAB_FONTWEIGHT = 'bold ';
    XLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    YLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    XTICK_FONTSIZE = AXES_FONT_SIZE_VIO;
    XLAB_FONTWEIGHT = 'bold';
    YLAB_FONTWEIGHT = 'bold';
    TITLE_FONTSIZE = AXES_FONT_SIZE_VIO;
    TITLE_FONTWEIGHT = 'bold';
    XLABEL_OFFSET = -.05;
    GROUP_LAB_YOFFSET = -0.275;
    % ax = gca;
    ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
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
    exportgraphics(fig,[save_dir filesep sprintf('speed_groupwise_intact_%s.tiff',meas_names{meas_i})],'Resolution',1000);
    %%
    horiz_shift = 0;
    vert_shift = 0;
    FIGURE_POSITION =[0,0,6.5,9];
    FONT_NAME = 'Arial';
    FONT_SIZE_VIO = 10;
    FONT_SIZE_VIO_REG = 7;
    % IM_RESIZE = 1;
    % AX_H  = 0.2;
    % AX_W = 0.275;
    AX_HORIZ_SHIFT = 0.06;
    AX_INIT_VERT_VIO = 0.6;
    TITLE_XSHIFT = 0.4;
    TITLE_YSHIFT = 0.975;
    TITLE_BOX_SZ = [0.4,0.4];
    AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
    %-
    inds = [stats_struct.mod_tag] == categorical(2) & [stats_struct.measure_tag] == categorical(meas_names(meas_i));
    tmp_stats = stats_struct(inds);
    VIOLIN_PARAMS = {'width',0.1,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
        'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
    PLOT_STRUCT = struct('color_map',COLOR_MAPS_SPEED,...
        'cond_labels',unique(tmp_table.trial_char),'group_labels',categorical({'YA','OHMA','OLMA'}),...
        'cond_offsets',[-0.35,-0.1,0.15,0.40],...
        'group_offsets',[0.125,0.475,0.812],...
        'y_label',meas_ylabel{meas_i},...
        'title',meas_names{meas_i},'font_size',FONT_SIZE_VIO,'ylim',YLIMS{meas_i},...
        'font_name','Arial','x_label','speed (m/s)','do_combine_groups',false,...
        'regresslab_txt_size',FONT_SIZE_VIO_REG);
    %-
    cond_ind = 2;
    anova_p = tmp_stats.anova_preds_p{cond_ind};
    grp_ind = 3;
    anova_grp_p = tmp_stats.anova_preds_p{grp_ind};
    speed_p = repmat(tmp_stats.mod_preds_p(2),3,1);
    grp_p = [0,tmp_stats.mod_preds_p(3:4)']';
    speed_r2 = tmp_stats.mod_r2;
    speed_xvals = (0:5)*0.25;
    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
        'pvals',{{}},...
        'pvals_pairs',{{}},...
        'pvals_grp',{num2cell(grp_p)},...
        'pvals_grp_pairs',{{[1,1],[1,2],[1,3]}},...
        'regress_pval',{num2cell(speed_p)},...
        'regress_line',{{[tmp_stats.mod_preds_coeff(1), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(3), tmp_stats.mod_preds_coeff(2)],...
            [tmp_stats.mod_preds_coeff(1)+tmp_stats.mod_preds_coeff(4), tmp_stats.mod_preds_coeff(2)]}},...
        'r2_coeff',{repmat(speed_r2,3,1)},...
        'regress_xvals',speed_xvals,...
        'subject_char',[],... % this option when filled prints removal of nan() info
        'group_order',categorical({''}),...
        'do_include_intercept',false);
    %-
    fig = figure('color','white','renderer','Painters');
    % annotation('textbox',[TITLE_XSHIFT-(TITLE_BOX_SZ(1)/2/2),TITLE_YSHIFT-(TITLE_BOX_SZ(2)/2),TITLE_BOX_SZ(1),TITLE_BOX_SZ(2)],...
    %     'String',atlas_name,'HorizontalAlignment','center',...
    %     'VerticalAlignment','middle','LineStyle','none','FontName',FONT_NAME,...
    %     'FontSize',14,'FontWeight','Bold','Units','normalized');
    set(fig,'Units','inches','Position',FIGURE_POSITION)
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    set(gca,AXES_DEFAULT_PROPS{:});
    ax = axes();
    ax = group_violin(tmp_table,meas_names{meas_i},'trial_char','group_id',...
        ax,...
        'VIOLIN_PARAMS',VIOLIN_PARAMS,...
        'PLOT_STRUCT',PLOT_STRUCT,...
        'STATS_STRUCT',STATS_STRUCT);
    set(ax,'OuterPosition',[0,0,1,1]);
    set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
    %- set ylabel & title
    AXES_FONT_SIZE_VIO = 10;
    GROUP_LAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    GROUP_LAB_FONTWEIGHT = 'bold ';
    XLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    YLAB_FONTSIZE = AXES_FONT_SIZE_VIO;
    XTICK_FONTSIZE = AXES_FONT_SIZE_VIO;
    XLAB_FONTWEIGHT = 'bold';
    YLAB_FONTWEIGHT = 'bold';
    TITLE_FONTSIZE = AXES_FONT_SIZE_VIO;
    TITLE_FONTWEIGHT = 'bold';
    XLABEL_OFFSET = -.05;
    GROUP_LAB_YOFFSET = -0.275;
    % ax = gca;
    ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
    ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
    ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
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
    exportgraphics(fig,[save_dir filesep sprintf('speed_groupwise_nointact_%s.tiff',meas_names{meas_i})],'Resolution',1000);
end