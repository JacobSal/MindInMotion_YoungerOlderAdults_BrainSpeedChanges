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
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
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
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
colormap(linspecer);
%## soft define
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
% save_dir = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'raw_data_vis'];
save_dir = [studies_fpath filesep 'mim_yaoa_mristudy'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## DATASET INFORMATION
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_mri_study');
%- override with already processed ALLEEG set
% SUBJ_PICS = cell(1,length(GROUP_NAMES));
% SUBJ_ITERS = cell(1,length(GROUP_NAMES));
% for i = 1:length(STUDY.datasetinfo)
%     gi = find(strcmp(STUDY.datasetinfo(i).group,GROUP_NAMES));
%     SUBJ_PICS{1,gi} = [SUBJ_PICS{1,gi}, {STUDY.datasetinfo(i).subject}];
%     if isempty(SUBJ_ITERS{1,gi})
%         SUBJ_ITERS{1,gi} = 1;
%     else
%         SUBJ_ITERS{1,gi} = [SUBJ_ITERS{1,gi}, SUBJ_ITERS{1,gi}(end)+1];
%     end
% end
%%
% CATEGORIES = {'YoungAdult','HF_OlderAdult','LF_OlderAdult'};
CATEGORIES = {'YoungAdult','HF_OlderAdult','LF_OlderAdult'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
N_TRIALS = 16;
%- Loop through directory
table_imu_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
table_trial_vec = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_vec = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_vec = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_ls_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
table_trial_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_ls = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
%- Loop through directory
table_imu_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
table_trial_vec_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_vec_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_conds = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_vec_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_ls_meas_conds = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
table_trial_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_ls_conds = cell(length([SUBJ_PICS{:}]),1);
table_subj_cat_ls_conds = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
subj_stack = [];
cnt_ga = 1;
cnt_ls = 1;
%%
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
        %## IMU (SACRAL MEASURES ONLY) 
        for f_i = 1:length(dir_imu)
            tmp = load([dir_imu(f_i).folder filesep dir_imu(f_i).name]);
            tmp = tmp.myStructure;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names{cnt_ls} = fieldnames(tmp);
            table_imu_meas(cnt_ls,1:length(values)) = values;
            tmp = strsplit(dir_imu(f_i).name,'.');
            fprintf(fid,'IMU Assigning %s.\n',dir_imu(f_i).name);
            if any(strcmp(tmp{1},{'TM_med_2','TM_med_1'}))
                trial_i = 'med';
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0p25';
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0p5';
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0p75';
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1p0';
            else
                trial_i = tmp{1};
            end
            table_trial_vec{cnt_ls} = trial_i; %tmp{1};
            table_subj_vec{cnt_ls} = SUBJ_PICS{group_i}{subj_i};
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
            end
            table_subj_cat_vec{cnt_ls} = cat_i;
            cnt_ls = cnt_ls + 1;
        end
        %## LOADSOL (GAIT)
        for f_i = 1:length(dir_ls)
            tmp = load([dir_ls(f_i).folder filesep dir_ls(f_i).name]);
            tmp = tmp.summarizedTrial;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_ls{cnt_ga} = fieldnames(tmp);
            table_ls_meas(cnt_ga,1:length(values)) = values;
            tmp = strsplit(dir_ls(f_i).name,'.');
            fprintf(fid,'GRF Assigning %s.\n',dir_ls(f_i).name);
            if any(strcmp(tmp{1},{'TM_med_2','TM_med_1'}))
                trial_i = 'med';
            elseif any(strcmp(tmp{1},{'TM_low_2','TM_low_1'}))
                trial_i = 'low';
            elseif any(strcmp(tmp{1},{'TM_flat_2','TM_flat_1'}))
                trial_i = 'flat';
            elseif any(strcmp(tmp{1},{'TM_high_2','TM_high_1'}))
                trial_i = 'high';
            elseif any(strcmp(tmp{1},{'SP_0p25_2','SP_0p25_1'}))
                trial_i = '0p25';
            elseif any(strcmp(tmp{1},{'SP_0p5_2','SP_0p5_1'}))
                trial_i = '0p5';
            elseif any(strcmp(tmp{1},{'SP_0p75_2','SP_0p75_1'}))
                trial_i = '0p75';
            elseif any(strcmp(tmp{1},{'SP_1p0_2','SP_1p0_1'}))
                trial_i = '1p0';
            else
                trial_i = tmp{1};
            end
            table_trial_ls{cnt_ga} = trial_i; %tmp{1};
            table_subj_ls{cnt_ga} = SUBJ_PICS{group_i}{subj_i};
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
            else
                cat_i = 'NaN';
            end
            table_subj_cat_ls{cnt_ga} = cat_i;
            cnt_ga = cnt_ga + 1;
        end
        
        %## IMU (SACRAL MEASURES ONLY) CONDITIONS
        for f_i = 1:length(dir_imu_conds)
            tmp = load([dir_imu_conds(f_i).folder filesep dir_imu_conds(f_i).name]);
            tmp = tmp.mergedOutputStruct;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_conds{cnt_ls} = fieldnames(tmp);
            table_imu_meas_conds(cnt_ls,1:length(values)) = values;
            tmp = strsplit(dir_imu_conds(f_i).name,'.');
            fprintf(fid,'IMU Assigning %s.\n',dir_imu_conds(f_i).name);
            if any(strcmp(tmp{1},{'TM_med'}))
                trial_i = 'med';
            elseif any(strcmp(tmp{1},{'TM_low'}))
                trial_i = 'low';
            elseif any(strcmp(tmp{1},{'TM_flat'}))
                trial_i = 'flat';
            elseif any(strcmp(tmp{1},{'TM_high'}))
                trial_i = 'high';
            elseif any(strcmp(tmp{1},{'SP_0p25'}))
                trial_i = '0p25';
            elseif any(strcmp(tmp{1},{'SP_0p5'}))
                trial_i = '0p5';
            elseif any(strcmp(tmp{1},{'SP_0p75'}))
                trial_i = '0p75';
            elseif any(strcmp(tmp{1},{'SP_1p0'}))
                trial_i = '1p0';
            else
                trial_i = tmp{1};
            end
            table_trial_vec_conds{cnt_ls} = trial_i; %tmp{1};
            table_subj_vec_conds{cnt_ls} = SUBJ_PICS{group_i}{subj_i};
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
            end
            table_subj_cat_vec_conds{cnt_ls} = cat_i;
            cnt_ls = cnt_ls + 1;
        end
        %## LOADSOL (GAIT) CONDITIONS
        for f_i = 1:length(dir_ls_conds)
            tmp = load([dir_ls_conds(f_i).folder filesep dir_ls_conds(f_i).name]);
            tmp = tmp.mergedOutputStruct;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_ls_conds{cnt_ga} = fieldnames(tmp);
            table_ls_meas_conds(cnt_ga,1:length(values)) = values;
            tmp = strsplit(dir_ls_conds(f_i).name,'.');
            fprintf(fid,'GRF Assigning %s.\n',dir_ls_conds(f_i).name);
            if any(strcmp(tmp{1},{'TM_med'}))
                trial_i = 'med';
            elseif any(strcmp(tmp{1},{'TM_low'}))
                trial_i = 'low';
            elseif any(strcmp(tmp{1},{'TM_flat'}))
                trial_i = 'flat';
            elseif any(strcmp(tmp{1},{'TM_high'}))
                trial_i = 'high';
            elseif any(strcmp(tmp{1},{'SP_0p25'}))
                trial_i = '0p25';
            elseif any(strcmp(tmp{1},{'SP_0p5'}))
                trial_i = '0p5';
            elseif any(strcmp(tmp{1},{'SP_0p75'}))
                trial_i = '0p75';
            elseif any(strcmp(tmp{1},{'SP_1p0'}))
                trial_i = '1p0';
            else
                trial_i = tmp{1};
            end
            table_trial_ls_conds{cnt_ga} = trial_i; %tmp{1};
            table_subj_ls_conds{cnt_ga} = SUBJ_PICS{group_i}{subj_i};
            if contains(SUBJ_PICS{group_i}{subj_i},'H1')
                cat_i = 'YoungAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H2')
                cat_i = 'HF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'H3')
                cat_i = 'LF_OlderAdult';
            elseif contains(SUBJ_PICS{group_i}{subj_i},'NH3')
                cat_i = 'LF_OlderAdult';
            else
                cat_i = 'NaN';
            end
            table_subj_cat_ls_conds{cnt_ga} = cat_i;
            cnt_ga = cnt_ga + 1;
        end
   
    end
end
fclose(fid);
%% IMU TABLE
rej = ~any(table_imu_meas == 0,2);
table_imu_meas = table_imu_meas(rej,:);
table_subj_vec = table_subj_vec(rej);
table_trial_vec = table_trial_vec(rej);
table_subj_cat_vec = table_subj_cat_vec(rej);
table_imu_out = array2table(table_imu_meas,'VariableNames',table_header_names{1});
table_imu_out.SubjectName = categorical(table_subj_vec);
table_imu_out.TrialName = categorical(table_trial_vec);
table_imu_out.SubjectCategory = categorical(table_subj_cat_vec);
writetable(table_imu_out,[save_dir filesep 'imu_table_out.xlsx']);
%- rearrange headers
% table_imu_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
% table_imu_out = table(categorical(table_subj_vec),categorical(table_trial_vec),table_imu_meas(:,1),table_imu_meas(:,2),...
%     table_imu_meas(:,3),table_imu_meas(:,4),table_imu_meas(:,5),table_imu_meas(:,6),...
%     table_imu_meas(:,7),table_imu_meas(:,8),table_imu_meas(:,9),table_imu_meas(:,10),...
%     table_imu_meas(:,11),table_imu_meas(:,12),'VariableNames',[{'Subject'},{'TrialName'},table_header_names{1}']);
%% LOADSOL TABLE
rej = ~any(table_ls_meas == 0,2);
table_ls_meas = table_ls_meas(rej,:);
table_subj_ls = table_subj_ls(rej);
table_trial_ls = table_trial_ls(rej);
table_subj_cat_ls = table_subj_cat_ls(rej);
table_ls_out = array2table(table_ls_meas,'VariableNames',table_header_names_ls{1});
table_ls_out.SubjectName = categorical(table_subj_ls);
table_ls_out.TrialName = categorical(table_trial_ls);
table_ls_out.SubjectCategory = categorical(table_subj_cat_ls);
writetable(table_ls_out,[save_dir filesep 'ls_table_out.xlsx']);
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
%% average across trials
subjs = unique(table_imu_out.SubjectName);
trials = unique(table_imu_out.TrialName);
table_new = [];
for i = 1:length(subjs)
    tmpvals = table_imu_out(table_imu_out.SubjectName == subjs(i),:);
    subj_cat = unique(table_imu_out(table_imu_out.SubjectName == subjs(i),:).SubjectCategory);
    for j = 1:length(trials)
        tmp = tmpvals(tmpvals.TrialName == trials(j),:);
%         tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
        tmp = varfun(@mean, tmp, 'InputVariables', @isnumeric);
        tmp.SubjectName = subjs(i);
        tmp.TrialName = trials(j);
        tmp.SubjectCategory = subj_cat;
        table_new = [table_new; tmp];
    end
end
table_new_imu = table_new;

%% average across trials
subjs = unique(table_ls_out.SubjectName);
trials = unique(table_ls_out.TrialName);
table_new = [];
for i = 1:length(subjs)
    tmpvals = table_ls_out(table_ls_out.SubjectName == subjs(i),:);
    subj_cat = unique(table_ls_out(table_ls_out.SubjectName == subjs(i),:).SubjectCategory);
    for j = 1:length(trials)
        tmp = tmpvals(tmpvals.TrialName == trials(j),:);
%         tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
        tmp = varfun(@mean, tmp, 'InputVariables', @isnumeric);
        tmp.SubjectName = subjs(i);
        tmp.TrialName = trials(j);
        tmp.SubjectCategory = subj_cat;
        table_new = [table_new; tmp];
    end
end
table_new_ls = table_new;

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
    ind = table_new_imu.SubjectName==ss;
    if ss_speed < SPEED_CUTOFF
%         inds_del = table_new_imu.SubjectName == ss;
%         table_new_imu(inds_del) = [];
        table_new_imu(ind,:) = [];
        table_new_ls(ind,:) = [];
    else
        table_new_imu.terrain_speed(ind) = ss_speed;
        table_new_ls.terrain_speed(ind) = ss_speed;
    end
end
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
    ind1 = tmp_table_imu.SubjectName==ss;
    ind2 = tmp_table_ls.SubjectName==ss;
    if any(ind1) & ~isempty(ss_var)
        for j = 1:length(trials)
            sub = strsplit(ss_var{j},',');
            dsub = double(string(sub));
            %- get subject and trial type index in imu/ls sheets
            ind11 = tmp_table_imu.SubjectName==ss & tmp_table_imu.TrialName==trials{j};
            ind22 = tmp_table_ls.SubjectName==ss & tmp_table_ls.TrialName==trials{j};
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
    tmp_table_imu(tmp_table_imu.(fn{f}) == 0,fn{f}) = {nan()};
    tmp_table_ls(tmp_table_ls.(fn{f}) == 0,fn{f}) = {nan()};
end
tmp_table_imu.stability_rating(tmp_table_imu.stability_rating==0) = nan();
tmp_table_ls.stability_rating(tmp_table_imu.stability_rating==0) = nan();
tmp_table_imu.og_length_m(tmp_table_imu.og_length_m==0) = nan();
tmp_table_ls.og_length_m(tmp_table_imu.og_length_m==0) = nan();
writetable(tmp_table_imu,[save_dir filesep 'imu_table_trial.xlsx']);
writetable(tmp_table_ls,[save_dir filesep 'ls_table_trial.xlsx']);
save([save_dir filesep 'imu_table_trial.mat'],'tmp_table_imu')
save([save_dir filesep 'ls_table_trial.mat'],'tmp_table_ls')
%%
%{
table_new_imu = readtable([save_dir filesep 'imu_table_meantrial.xlsx']);
table_new_ls = readtable([save_dir filesep 'ls_table_meantrial.xlsx']);
%}
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
meas_ylabel = {'Coefficient of Variation','Coefficient of Variation'};
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
% table_in = table_imu_out;
for meas_i = 1:length(meas_names)
    table_in = table_new_imu; %table_imu_out;
    meas_names{meas_i} = meas_names{meas_i};
    dat_mean = nanmean(table_in.(meas_names{meas_i}));
    dat_std = nanstd(table_in.(meas_names{meas_i}));
    in_crop = table_in.(meas_names{meas_i})>dat_mean-3*dat_std & table_in.(meas_names{meas_i})<dat_mean+3*dat_std;
    table_wi = table_in(in_crop,:);
    outliers = table_in(~in_crop,:);
    a = unstack(table_wi,meas_names{meas_i},'SubjectName');
    b = unstack(table_wi,meas_names{meas_i},'TrialName');
    b_2 = unstack(table_wi,meas_names{meas_i},'SubjectCategory');
    b_2 = b_2(:,12:end);
%     g_cats = unique(table_wi.SubjectCategory);
    cond_1 = cell(length(speed_chars),length(g_cats));
    cond_2 = cell(length(terrain_chars),length(g_cats));
    for g_i = 1:size(cond_1,2)
        for c_i = 1:size(cond_1,1)
            tmp = b_2.TrialName == speed_chars{c_i};
            g_ind = strcmp(char(g_cats(g_i)),b_2.Properties.VariableNames);
            cond_1{c_i,g_i} = b_2{tmp,g_ind};
        end
    end
    for g_i = 1:size(cond_2,2)
        for c_i = 1:size(cond_2,1)
            tmp = b_2.TrialName == terrain_chars{c_i};
            g_ind = strcmp(char(g_cats(g_i)),b_2.Properties.VariableNames);
            cond_2{c_i,g_i} = b_2{tmp,g_ind};
        end
    end
    tmp = b(:,13:end);
    b_speed = tmp(:,3:6);
    b_terrain = [tmp(:,7),tmp(:,9),tmp(:,10),tmp(:,8)];
    %% By subject plot
%     figure;
%     title(sprintf('%s Across Subjects',meas_names{meas_i}));
%     hold on;
%     violinplot(a(:,15:end),a(:,15),...
%         'ViolinColor',linspecer(size(a(:,13:end),2)),...
%         'width',VIOLIN_WIDTH,...
%         'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
%         'ShowMedian',true,'Bandwidth',range(range(a{:,15:end}))*0.2,'QuartileStyle','shadow',...
%         'HalfViolin','full','DataStyle','scatter');
%     ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
%     xlabel('Subject Code');
%     xtickangle(45)
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = [10,100,1280,320];
% %     saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',meas_names{meas_i})]);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',meas_names{meas_i})]);
    %% (SPEED) STATISTICS GROUP & CONDITION EFFECT
    %## REFORM TABLE VARIABLES
    inds = [];
    for i = 1:length(speed_chars)
        inds = [inds; find(table_new_imu.TrialName == speed_chars{i})];
    end
    tmp = table_new_imu.(meas_names{meas_i});
    meas_in = double(tmp(inds,:));
    tmp = table_new_imu.TrialName;
    [vals_tn,~,cati] = unique(tmp(inds,:));
    cat_1 = categorical(cati); % trial categories
    tmp = table_new_imu.SubjectCategory;
    cat_2 = categorical(tmp(inds,:)); % subject categories
    tmp_t = table(meas_in,cat_1,cat_2);
    %## MIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2';
    mdl_speed_mix_gc = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','normal');
    disp(mdl_speed_mix_gc)
    %- Convert summary to char array
    txt = evalc('mdl_speed_mix_gc');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## (MAIN TEST-SPEED) CONDITION NO MIXED EFFECTS
    modelspec = 'meas_in~1+cat_1';
    mdl_speed_mixc = fitglme(tmp_t,modelspec);
%     [p,T,stats,terms] = anovan(tmp_t.meas_in,{tmp_t.cat_1},'alpha',0.05,'sstype',2,'model','linear');
    [stats] = anova(mdl_speed_mixc);
    speed_mixc_f = stats{2,5};
%     comp = multcompare(stats);
    modelspec = 'meas_in~1';
    mdl_comp = fitlme(tmp_t,modelspec);
    disp(mdl_speed_mixc)
    disp(anova(mdl_speed_mixc));
    fid = fopen([save_dir filesep sprintf('%s_across_mixspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    %- converst Summary anova to char array
    txt = evalc('anova(mdl_speed_mixc)');
    fprintf(fid,'ANOVA EVAL\n');
    fprintf(fid,'%s',txt);
    fprintf(fid,'\n');
    %- Convert summary to char array
    fprintf(fid,'FITGLME EVAL\n');
    txt = evalc('mdl_speed_mixc');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    %- multiple comparisons
    for y = 1:numel(mdl_speed_mixc.CoefficientNames)
        fprintf(fid,'Contrast position %i: %s\n', y, char(mdl_speed_mixc.CoefficientNames{y}));
    end
    p_mcomp = [];
%     [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [1,0,0,0]); fprintf(fid,'multicomp intercept:intercept, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,0,0]); fprintf(fid,'multicomp intercept:cat_1_2, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,1,0]); fprintf(fid,'multicomp intercept:cat_1_3, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,0,1]); fprintf(fid,'multicomp intercept:cat_1_4, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,2,0]); fprintf(fid,'multicomp cat_1_2:cat_1_3,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,0,2]); fprintf(fid,'multicomp cat_1_2:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,1,2]); fprintf(fid,'multicomp cat_1_3:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    
    [h,crit_p,adj_ci_cvrg,adj_p_speed] = fdr_bh(p_mcomp,0.05,'pdep','no');
    fprintf(fid,'False Discovery Corrected:\n');
    fprintf(fid,'%0.4f\n',adj_p_speed); fprintf(fid,'%i\n',h); fprintf(fid,'critical fdr_p: %0.4f\n',crit_p);
    %- cohens f2 test
%     out = mes1way(tmp_t.meas_in,'eta2','group',double(tmp_t.cat_1));
%     R2 = mdl_terrain_mixc.Rsquared.Ordinary;
    R21 = mdl_comp.SSR/mdl_comp.SST;
    R22 = mdl_speed_mixc.SSR/mdl_speed_mixc.SST; %mdl_terrain_mixc.Rsquared.Ordinary; %mdl_terrain_mixc.Rsquared.Adjusted;
% 	R2 = mdl_terrain_mixc.Rsquared.Adjusted;
    cohens_f2 = (R22-R21)/(1-R22);
    fprintf(fid,'cohens f2: %0.4f\n',cohens_f2);
    fprintf('Cohens f2: %0.4f\n',cohens_f2);
    fclose(fid);
    %% (TERRAIN) STATISTICS GROUP & CONDITION EFFECT
    %## REFORM TABLE VARIABLES
    inds = [];
    for i = 1:length(terrain_chars)
        inds = [inds; find(table_new_imu.TrialName == terrain_chars{i})];
    end
    tmp = table_new_imu.(meas_names{meas_i});
    meas_in = double(tmp(inds,:));
    tmp = table_new_imu.TrialName;
    [vals_tn,~,cati] = unique(tmp(inds,:));
    cat_1 = categorical(cati); % trial categories
    tmp = table_new_imu.SubjectCategory;
    cat_2 = categorical(tmp(inds,:)); % subject categories
    tmp = table_new_imu.terrain_speed;
    coev_1 = double(tmp(inds,:));
    tmp_t = table(meas_in,cat_1,cat_2,coev_1);
    %## MIXED & MIXED-EFFECTS GROUP & COND (COVARIATE TERRAIN SPEED)
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2+(1|coev_1)';
    mdl_terrain_mix_gc_cov = fitglme(tmp_t,modelspec,'Distribution','normal');
    disp(mdl_terrain_mix_gc_cov)
    %- Convert summary to char array
    txt = evalc('mdl_terrain_mix_gc_cov');
    fid = fopen([save_dir filesep sprintf('%s_across_mixxedsubjectcatsNterrain_covarspeed_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## (MAIN TEST-TERRAIN) CONDITION MIXED EFFECTS
    modelspec = 'meas_in~1+cat_1+(1|coev_1)';
%     mdl_terrain_mixc = fitglme(tmp_t,modelspec);
    mdl_terrain_mixc = fitlme(tmp_t,modelspec);
    %-
%     modelspec = 'meas_in~1+(1|coev_1)';
    modelspec = 'meas_in~1';
    mdl_comp = fitlme(tmp_t,modelspec);
%     [p,T,stats,terms] = anovan(tmp_t.meas_in,{tmp_t.cat_1},'alpha',0.05,'sstype',2,'model','linear');
    [stats] = anova(mdl_terrain_mixc);
    terrain_mixc_f = stats{2,5};
%     comp = multcompare(stats);
    disp(mdl_terrain_mixc)
    disp(anova(mdl_terrain_mixc));
    fid = fopen([save_dir filesep sprintf('%s_across_mixterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    %- converst Summary anova to char array
    txt = evalc('anova(mdl_terrain_mixc)');
    fprintf(fid,'ANOVA EVAL\n');
    fprintf(fid,'%s',txt);
    fprintf(fid,'\n');
    %- Convert summary to char array
    fprintf(fid,'FITGLME EVAL\n');
    txt = evalc('mdl_terrain_mixc');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    %- multiple comparisons
    for y = 1:numel(mdl_terrain_mixc.CoefficientNames)
        fprintf(fid,'Contrast position %i: %s\n', y, char(mdl_terrain_mixc.CoefficientNames{y}));
    end
    p_mcomp = [];
%     [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [1,0,0,0]); fprintf(fid,'multicomp intercept:intercept, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,0,0]); fprintf(fid,'multicomp intercept:cat_1_2, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,1,0]); fprintf(fid,'multicomp intercept:cat_1_3, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,0,1]); fprintf(fid,'multicomp intercept:cat_1_4, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,2,0]); fprintf(fid,'multicomp cat_1_2:cat_1_3,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,0,2]); fprintf(fid,'multicomp cat_1_2:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,1,2]); fprintf(fid,'multicomp cat_1_3:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    
    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(p_mcomp,0.05,'pdep','no');
    fprintf(fid,'False Discovery Corrected:\n');
    fprintf(fid,'%0.4f\n',adj_p); fprintf(fid,'%i\n',h); fprintf(fid,'critical fdr_p: %0.4f\n',crit_p);
    %- cohens f2 test
%     out = mes1way(tmp_t.meas_in,'eta2','group',double(tmp_t.cat_1));
%     R2 = mdl_terrain_mixc.Rsquared.Ordinary;
    R21 = mdl_comp.SSR/mdl_comp.SST;
    R22 = mdl_terrain_mixc.SSR/mdl_terrain_mixc.SST; %mdl_terrain_mixc.Rsquared.Ordinary; %mdl_terrain_mixc.Rsquared.Adjusted;
% 	R2 = mdl_terrain_mixc.Rsquared.Adjusted;
    cohens_f2 = (R22-R21)/(1-R22);
    fprintf(fid,'cohens f2: %0.4f\n',cohens_f2);
    fprintf('Cohens f2: %0.4f\n',cohens_f2);
    fclose(fid);
    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_terrain,trial_names_terrain','width',VIOLIN_WIDTH,...
        'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',COLORS_MAPS_TERRAIN,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    %- plot outliers
%     vals_x = get(gca,'XTick');
%     for c_i = 1:length(terrain_chars)
%         idx = outliers.TrialName == terrain_chars{c_i};
%         vals_y = outliers.(meas_names{meas_i})(idx);
%         if ~isempty(vals_y)
%             in_x = repmat(vals_x(c_i),size(vals_y));
%             scatter(in_x,vals_y,'rx','jitter','on', 'jitterAmount', 0.05)
%         end
%     end
    %- plot siglines
    fig_i = get(groot,'CurrentFigure');
    ind = 1;
    if terrain_mixc_f < 0.05
        if adj_p(2) < 0.05
            sigline([1 2],'*',[],[],adj_p(2)); %mdl_terrain_c.Coefficients.pValue(2))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
    %         ind = ind+2;
        end
        if adj_p(3) < 0.05
            sigline([1 3],'*',[],[],adj_p(3)); %mdl_terrain_c.Coefficients.pValue(3))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
    %         ind = ind+2;
        end
        if adj_p(1) < 0.05
            sigline([1 4],'*',[],[],adj_p(1)); %mdl_terrain_c.Coefficients.pValue(4))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
        end
    end
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    ax = gca; %fig_i.CurrentAxes;
    set(ax, 'box', 'off')
    set(ax,'LineWidth',1)
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_terrain);
    xtickangle(30);
    %-
%     xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
%     pos1=get(xlh,'Position');
%     pos1(1,2)=pos1(1,2)-0.02;
%     set(xlh,'Position',pos1);
    %-
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
    set(ax,'FontName','Arial','FontSize',14,'FontWeight','bold')
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.tiff',meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.pdf',meas_names{meas_i})],'ContentType','vector','Resolution',300);
    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_speed,trial_names_speed,'width',VIOLIN_WIDTH,...
        'GroupOrder',b_speed.Properties.VariableNames,'ViolinColor',COLOR_MAPS_SPEED,...
        'ShowWhiskers',true,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_speed{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    %- plot outliers
%     vals_x = get(gca,'XTick');
%     for c_i = 1:length(speed_chars)
%         idx = outliers.TrialName == speed_chars{c_i};
%         vals_y = outliers.(meas_names{meas_i})(idx);
%         if ~isempty(vals_y)
%             in_x = repmat(vals_x(c_i),size(vals_y));
%             scatter(in_x,vals_y,'rx','jitter','on', 'jitterAmount', 0.05)
%         end
%     end
    %- plot siglines
    fig_i = get(groot,'CurrentFigure');
    ind = 1;
    if speed_mixc_f < 0.05
        if adj_p_speed(2) < 0.05
            sigline([1 2],'*',[],[],adj_p_speed(2)); %mdl_terrain_c.Coefficients.pValue(2))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
    %         ind = ind+2;
        end
        if adj_p_speed(3) < 0.05
            sigline([1 3],'*',[],[],adj_p_speed(3)); %mdl_terrain_c.Coefficients.pValue(3))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
    %         ind = ind+2;
        end
        if adj_p_speed(1) < 0.05
            sigline([1 4],'*',[],[],adj_p_speed(1)); %mdl_terrain_c.Coefficients.pValue(4))
            fig_i.Children.Children(ind).FontSize = 16;
            fig_i.Children.Children(ind).FontName = 'Arial';
            fig_i.Children.Children(ind+1).LineWidth = 1.25;
        end
    end
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    ax = gca; %fig_i.CurrentAxes;
    set(ax, 'box', 'off')
    set(ax,'LineWidth',1)
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_speed);
    xtickangle(30);
    %-
    xlh = xlabel('Speed (m/s)','Units','normalized');
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    %-
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    set(ax,'FontName','Arial','FontSize',14,'FontWeight','bold')
    hold off;
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.tiff',meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.pdf',meas_names{meas_i})],'ContentType','vector','Resolution',300);
    %% (PLOT) Trial & Subject Category Plot for High vs Low function OA
    vals = cat(1,cond_1{1,:});
    bandwidth = range(vals)*0.1;
    for g_i = 1:size(cond_1,2)
        n=max(cellfun(@numel,cond_1(:,g_i)));
        for c_i = 1:size(cond_1,1)
    %         tmp = cond_1{c_i,:};
            cond_1{c_i,g_i}(end+1:n,:) = nan;
        end
    end
    varargin = {'width',VIOLIN_WIDTH_GROUP,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',bandwidth,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter'};
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
%     c_maps = linspecer(size(cond_1,1));
%     offset = 0;
    cnt = 1;
    cnt_g = 1;
    xticks = [];
    xtick_labs = {};
    %- group primary, cond secondary
    save_lab = 'group_main';
    cond_offsets = [-0.3,-0.1,0.1,0.3];
    tmp_outliers = cell(size(cond_1));
    for g_i=1:size(cond_1,2)
        for c_i=1:size(cond_1,1)
            offset = cond_offsets(c_i);
%             varargin{end} = {COLORS_MAPS(c_i,:)};
            violins(cnt) = Violin(cond_1(c_i,g_i),g_i,varargin{:});
            violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
            violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
            violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
            violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
            violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
            violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
            violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
            violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
            violins(cnt).ViolinColor = {COLOR_MAPS_SPEED(c_i,:)};
            xticks = [xticks,g_i+offset];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_speed{c_i})}];
            xtick_labs = [xtick_labs,{sprintf('%s',trial_names_speed{c_i})}];
            g_ind = (char(g_cats(g_i))==outliers.SubjectCategory) & (char(speed_chars(c_i))==outliers.TrialName);
%             tmp_outliers{g_i,c_i} = outliers.(meas_names{meas_i})(g_ind);
            vals_y = outliers.(meas_names{meas_i})(g_ind);
            if ~isempty(vals_y)
                in_x = repmat(xticks(cnt),size(vals_y));
                scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
            end
            cnt = cnt+1;
        end
%         text(nanmean(xticks(cnt_g)),-8,char(group_names(g_i)),'FontSize',11);
%         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
        cnt_g = cnt_g+cnt-1;
    end
    %- set figure color, units, and size
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 12;
    %- set axes units, and size
    ax = gca; %fig_i.CurrentAxes;
    set(ax, 'box', 'off')
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',11,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.20 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
    set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',xtick_labs);
    xtickangle(30);
    xlh = xlabel('Speed (m/s)','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.12;
    set(xlh,'Position',pos1);
    %-
    %- set group labels
    if size(cond_1,2) == 2
        shift = 0;
        for g_i = 1:size(cond_1,2)
            text(0.25+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.50;
        end
    elseif size(cond_1,2) == 3
        shift = 0;
        for g_i = 1:size(cond_1,2)
            text(0.15+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.34;
        end
    else
        error('here');
    end
    %- set ylabel & title
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.hdf',save_lab,meas_names{meas_i})]);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.tiff',save_lab,meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'ContentType','vector','Resolution',300);
    %% (PLOT) Trial & Subject Category Plot for High vs Low function OA
    vals = cat(1,cond_2{1,:});
    bandwidth = range(vals)*0.1;
    for g_i = 1:size(cond_2,2)
        n=max(cellfun(@numel,cond_2(:,g_i)));
        for c_i = 1:size(cond_2,1)
    %         tmp = cond_2{c_i,:};
            cond_2{c_i,g_i}(end+1:n,:) = nan;
        end
    end
    varargin = {'width',VIOLIN_WIDTH_GROUP,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',bandwidth,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter'};
    figure;
    hold on;
    cnt = 1;
    cnt_g = 1;
    xticks = [];
    xtick_labs = {};
    %- group primary, cond secondary
    save_lab = 'group_main';
    cond_offsets = [-0.3,-0.1,0.1,0.3];
    for g_i = 1:size(cond_2,2)
        for c_i=1:size(cond_2,1)
            offset = cond_offsets(c_i);
            violins(cnt) = Violin(cond_2(c_i,g_i),g_i,varargin{:});
            violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
            violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
            violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
            violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
            violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
            violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
            violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
            violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
            violins(cnt).ViolinColor = {COLORS_MAPS_TERRAIN(c_i,:)};
            xticks = [xticks,g_i+offset];
            xtick_labs = [xtick_labs,{sprintf('%s',trial_names_terrain{c_i})}];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_terrain{c_i})}];
            g_ind = (char(g_cats(g_i))==outliers.SubjectCategory) & (char(terrain_chars(c_i))==outliers.TrialName);
%             tmp_outliers{g_i,c_i} = outliers.(meas_names{meas_i})(g_ind);
            vals_y = outliers.(meas_names{meas_i})(g_ind);
            if ~isempty(vals_y)
                in_x = repmat(xticks(cnt),size(vals_y));
                scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
            end
            cnt = cnt+1;
        end
        cnt_g = cnt_g+cnt-1;
    end
    %- set figure color, units, and size
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 12;
    %- set axes units, and size
    ax = gca; %fig_i.CurrentAxes;
    set(ax, 'box', 'off')
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',11,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.20 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
    set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',xtick_labs);
    xtickangle(30);
    xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.12;
    set(xlh,'Position',pos1);
%     xticks = get(ax,'XTick');
    %- set group labels
    if size(cond_2,2) == 2
        shift = 0;
        for g_i = 1:size(cond_2,2)
            text(0.25+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
            shift = shift+0.50;
        end
    elseif size(cond_2,2) == 3
        shift = 0;
        for g_i = 1:size(cond_2,2)
            text(0.15+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
            shift = shift+0.34;
        end
    else
        error('here');
    end
    %- set ylabel & title
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.tiff',save_lab,meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'ContentType','vector','Resolution',300);
    
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
meas_names = {'mean_StepDur','mean_GaitCycleDur','mean_SwingDur',...
            'mean_StanceDur','mean_SingleSupport','mean_TotalDS'};
meas_units = {'s','s','s','s','s','s'};
meas_titles = {'mean_StepDur','mean_GaitCycleDur','mean_SwingDur',...
            'mean_StanceDur','mean_SingleSupport','mean_TotalDS'};
% meas_titles = {'Step Duration',{'Step Duration';'Coefficient of Variation'}};
% meas_ylabel = {'Duration','Coefficient of Variation'};
meas_ylabel = {'Duration','Duration','Duration','Duration','Duration','Duration'};
YLIMS = {[0,2.5],[0,3],[0,1],[0,3],[0,1.5],[0,2]};
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
for meas_i = 1:length(meas_names)
    table_in = table_new_ls; %table_ls_out;
    dat_mean = nanmean(table_in.(meas_names{meas_i}));
    dat_std = nanstd(table_in.(meas_names{meas_i}));
    in_crop = table_in.(meas_names{meas_i})>dat_mean-3*dat_std & table_in.(meas_names{meas_i})<dat_mean+3*dat_std;
    table_wi = table_in(in_crop,:);
    outliers = table_in(~in_crop,:);
    meas_names{meas_i} = meas_names{meas_i};
    a = unstack(table_wi,meas_names{meas_i},'SubjectName');
    b = unstack(table_wi,meas_names{meas_i},'TrialName');
    b_2 = unstack(table_wi,meas_names{meas_i},'SubjectCategory');
    b_2 = b_2(:,65:end);
%     g_cats = unique(table_wi.SubjectCategory);
    cond_1 = cell(length(speed_chars),length(g_cats));
    cond_2 = cell(length(terrain_chars),length(g_cats));
    for g_i = 1:size(cond_1,2)
        for c_i = 1:size(cond_1,1)
            tmp = b_2.TrialName == speed_chars{c_i};
            g_ind = strcmp(char(g_cats(g_i)),b_2.Properties.VariableNames);
            cond_1{c_i,g_i} = b_2{tmp,g_ind};
        end
    end
    for g_i = 1:size(cond_2,2)
        for c_i = 1:size(cond_2,1)
            tmp = b_2.TrialName == terrain_chars{c_i};
            g_ind = strcmp(char(g_cats(g_i)),b_2.Properties.VariableNames);
            cond_2{c_i,g_i} = b_2{tmp,g_ind};
        end
    end
    tmp = b(:,66:end);
    b_speed = tmp(:,2:5);
    b_terrain = [tmp(:,6),tmp(:,8),tmp(:,9),tmp(:,7)];
    %% (SPEED) STATISTICS GROUP & CONDITION EFFECT
    %## REFORM TABLE VARIABLES
    inds = [];
    for i = 1:length(speed_chars)
        inds = [inds; find(table_new_ls.TrialName == speed_chars{i})];
    end
    tmp = table_new_ls.(meas_names{meas_i});
    meas_in = double(tmp(inds,:));
    tmp = table_new_ls.TrialName;
    [vals_tn,~,cati] = unique(tmp(inds,:));
    cat_1 = categorical(cati); % trial categories
    tmp = table_new_ls.SubjectCategory;
    cat_2 = categorical(tmp(inds,:)); % subject categories
    tmp_t = table(meas_in,cat_1,cat_2);
    %## MIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2';
    mdl_speed_mix_gc = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_speed_mix_gc)
    %- Convert summary to char array
    txt = evalc('mdl_speed_mix_gc');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## (MAIN TEST) CONDITION MIXED EFFECTS
    modelspec = 'meas_in~1+cat_1';
    mdl_speed_mixc = fitglme(tmp_t,modelspec);
    %-
    modelspec = 'meas_in~1';
%     mdl_terrain_mixc = fitglme(tmp_t,modelspec);
    mdl_comp = fitlme(tmp_t,modelspec);
    stats = anova(mdl_speed_mixc);
    speed_mixc_f = stats{2,5};
    disp(mdl_speed_mixc)
    disp(anova(mdl_speed_mixc));
    fid = fopen([save_dir filesep sprintf('%s_across_mixspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    %- converst Summary anova to char array
    txt = evalc('anova(mdl_speed_mixc)');
    fprintf(fid,'ANOVA EVAL\n');
    fprintf(fid,'%s',txt);
    fprintf(fid,'\n');
    %- Convert summary to char array
    fprintf(fid,'FITGLME EVAL\n');
    txt = evalc('mdl_speed_mixc');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    %- multiple comparisons
    for y = 1:numel(mdl_speed_mixc.CoefficientNames)
        fprintf(fid,'Contrast position %i: %s\n', y, char(mdl_speed_mixc.CoefficientNames{y}));
    end
    p_mcomp = [];
%     [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [1,0,0,0]); fprintf(fid,'multicomp intercept:intercept, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,0,0]); fprintf(fid,'multicomp intercept:cat_1_2, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,1,0]); fprintf(fid,'multicomp intercept:cat_1_3, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,0,1]); fprintf(fid,'multicomp intercept:cat_1_4, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,2,0]); fprintf(fid,'multicomp cat_1_2:cat_1_3,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,1,0,2]); fprintf(fid,'multicomp cat_1_2:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_speed_mixc, [0,0,1,2]); fprintf(fid,'multicomp cat_1_3:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    
    
    [h,crit_p,adj_ci_cvrg,adj_p_speed] = fdr_bh(p_mcomp,0.05,'pdep','no');
    fprintf(fid,'False Discovery Corrected:\n');
    fprintf(fid,'%0.4f\n',adj_p_speed); fprintf(fid,'%i\n',h); fprintf(fid,'critical fdr_p: %0.4f\n',crit_p);
    %- cohens f2 test
%     out = mes1way(tmp_t.meas_in,'eta2','group',double(tmp_t.cat_1));
%     R2 = mdl_terrain_mixc.Rsquared.Ordinary;
    R21 = mdl_comp.SSR/mdl_comp.SST;
    R22 = mdl_speed_mixc.SSR/mdl_speed_mixc.SST; %mdl_terrain_mixc.Rsquared.Ordinary; %mdl_terrain_mixc.Rsquared.Adjusted;
% 	R2 = mdl_terrain_mixc.Rsquared.Adjusted;
    cohens_f2 = (R22-R21)/(1-R22);
    fprintf(fid,'cohens f2: %0.4f\n',cohens_f2);
    fprintf('Cohens f2: %0.4f\n',cohens_f2);
    %-
    fclose(fid);
    %% (TERRAIN) STATISTICS GROUP & CONDITION EFFECT
    %## REFORM TABLE VARIABLES
    inds = [];
    for i = 1:length(terrain_chars)
        inds = [inds; find(table_new_ls.TrialName == terrain_chars{i})];
    end
    tmp = table_new_ls.(meas_names{meas_i});
    meas_in = double(tmp(inds,:));
    tmp = table_new_ls.TrialName;
    [vals_tn,~,cati] = unique(tmp(inds,:));
    cat_1 = categorical(cati); % trial categories
    tmp = table_new_ls.SubjectCategory;
    cat_2 = categorical(tmp(inds,:)); % subject categories
    tmp = table_new_ls.terrain_speed;
    coev_1 = double(tmp(inds,:));
    tmp_t = table(meas_in,cat_1,cat_2,coev_1);
    %## MIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2';
    mdl_terrain_mix_gc = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_terrain_mix_gc)
    %- Convert summary to char array
    txt = evalc('mdl_terrain_mix_gc');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## MIXED & MIXED-EFFECTS GROUP & COND (COVARIATE TERRAIN SPEED)
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2+(1|coev_1)';
    mdl_terrain_mix_gc_cov = fitglme(tmp_t,modelspec);
    disp(mdl_terrain_mix_gc_cov)
    %- Convert summary to char array
    txt = evalc('mdl_terrain_mix_gc_cov');
    fid = fopen([save_dir filesep sprintf('%s_across_mixxedsubjectcatsNterrain_covarspeed_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## (MAIN TEST) CONDITION MIXED EFFECTS
    modelspec = 'meas_in~1+cat_1+(1|coev_1)';
%     mdl_terrain_mixc = fitglme(tmp_t,modelspec);
    mdl_terrain_mixc = fitlme(tmp_t,modelspec);
    %-
    modelspec = 'meas_in~1+(1|coev_1)';
%     mdl_terrain_mixc = fitglme(tmp_t,modelspec);
    mdl_comp = fitlme(tmp_t,modelspec);
    stats = anova(mdl_terrain_mixc);
    terrain_mixc_f = stats{2,5};
    disp(mdl_terrain_mixc)
    disp(anova(mdl_terrain_mixc));
    fid = fopen([save_dir filesep sprintf('%s_across_mixterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    %- converst Summary anova to char array
    txt = evalc('anova(mdl_terrain_mixc)');
    fprintf(fid,'ANOVA EVAL\n');
    fprintf(fid,'%s',txt);
    fprintf(fid,'\n');
    %- Convert summary to char array
    fprintf(fid,'FITGLME EVAL\n');
    txt = evalc('mdl_terrain_mixc');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    %- multiple comparisons
    for y = 1:numel(mdl_terrain_mixc.CoefficientNames)
        fprintf(fid,'Contrast position %i: %s\n', y, char(mdl_terrain_mixc.CoefficientNames{y}));
    end
    p_mcomp = [];
%     [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [1,0,0,0]); fprintf(fid,'multicomp intercept:intercept, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,0,0]); fprintf(fid,'multicomp intercept:cat_1_2, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,1,0]); fprintf(fid,'multicomp intercept:cat_1_3, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,0,1]); fprintf(fid,'multicomp intercept:cat_1_4, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,2,0]); fprintf(fid,'multicomp cat_1_2:cat_1_3,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,1,0,2]); fprintf(fid,'multicomp cat_1_2:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    [pVal F df1 df2] = coefTest(mdl_terrain_mixc, [0,0,1,2]); fprintf(fid,'multicomp cat_1_3:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    
    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(p_mcomp,0.05,'pdep','no');
    fprintf(fid,'False Discovery Corrected:\n');
    fprintf(fid,'%0.4f\n',adj_p); fprintf(fid,'%i\n',h); fprintf(fid,'critical fdr_p: %0.4f\n',crit_p);
    %- cohens f2 test
%     out = mes1way(tmp_t.meas_in,'eta2','group',double(tmp_t.cat_1));
%     R2 = mdl_terrain_mixc.Rsquared.Ordinary;
    R21 = mdl_comp.SSR/mdl_comp.SST;
    R22 = mdl_terrain_mixc.SSR/mdl_terrain_mixc.SST; %mdl_terrain_mixc.Rsquared.Ordinary; %mdl_terrain_mixc.Rsquared.Adjusted;
% 	R2 = mdl_terrain_mixc.Rsquared.Adjusted;
    cohens_f2 = (R22-R21)/(1-R22);
    fprintf(fid,'cohens f2: %0.4f\n',cohens_f2);
    fprintf('Cohens f2: %0.4f\n',cohens_f2);
    fclose(fid);
    %% By trial plot
    figure;
    hold on;
    violinplot(b_terrain,trial_names_terrain','width',VIOLIN_WIDTH,...
        'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',COLORS_MAPS_TERRAIN,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(terrain_chars)
        idx = outliers.TrialName == terrain_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'ko','jitter','on', 'jitterAmount', 0.05)
        end
    end
    hold off;
    fig_i = get(groot,'CurrentFigure');
    ind = 1;
    if terrain_mixc_f < 0.05
    if adj_p(2) < 0.05
        sigline([1 2],'*',[],[],adj_p(2)); %mdl_terrain_c.Coefficients.pValue(2))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
%         ind = ind+2;
    end
    if adj_p(3) < 0.05
        sigline([1 3],'*',[],[],adj_p(3)); %mdl_terrain_c.Coefficients.pValue(3))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
%         ind = ind+2;
    end
    if adj_p(1) < 0.05
        sigline([1 4],'*',[],[],adj_p(1)); %mdl_terrain_c.Coefficients.pValue(4))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
    end
    end
    
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    ax = gca; %fig_i.CurrentAxes;
    set(ax, 'box', 'off')
    set(ax,'LineWidth',1)
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_terrain);
    xtickangle(30);
%     xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
%     pos1=get(xlh,'Position');
%     pos1(1,2)=pos1(1,2)-0.02;
%     set(xlh,'Position',pos1);
    %- set ylabel
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    set(ax,'FontName','Arial','FontSize',14,'FontWeight','bold')
    hold off;
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})],'Resolution',300); 
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.tiff',meas_names{meas_i})],'Resolution',900); 
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.pdf',meas_names{meas_i})],'ContentType','vector','Resolution',300,'ContentType','vector'); 
    %% By trial plot
    figure;
    hold on;
    violinplot(b_speed,trial_names_speed,'width',VIOLIN_WIDTH,...
        'GroupOrder',b_speed.Properties.VariableNames,'ViolinColor',COLOR_MAPS_SPEED,...
        'ShowWhiskers',true,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_speed{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(speed_chars)
        idx = outliers.TrialName == speed_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'ko','jitter','on', 'jitterAmount', 0.05)
        end
    end
    fig_i = get(groot,'CurrentFigure');
    ind = 1;
    if speed_mixc_f < 0.05
    if adj_p_speed(2) < 0.05
        sigline([1 2],'*',[],[],adj_p_speed(2)); %mdl_terrain_c.Coefficients.pValue(2))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
%         ind = ind+2;
    end
    if adj_p_speed(3) < 0.05
        sigline([1 3],'*',[],[],adj_p_speed(3)); %mdl_terrain_c.Coefficients.pValue(3))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
%         ind = ind+2;
    end
    if adj_p_speed(1) < 0.05
        sigline([1 4],'*',[],[],adj_p_speed(1)); %mdl_terrain_c.Coefficients.pValue(4))
        fig_i.Children.Children(ind).FontSize = 16;
        fig_i.Children.Children(ind).FontName = 'Arial';
        fig_i.Children.Children(ind+1).LineWidth = 1.25;
    end
    end
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_speed);
    xtickangle(30);
    xlh = xlabel('Speed (m/s)','Units','normalized');
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    %- set xlabel 2
%     text(0.45,-0.2,0,'Terrain Difficulty','FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
%         'Units','normalized');
    %- set ylabel
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    set(ax,'FontName','Arial','FontSize',14,'FontWeight','bold')
    set(ax, 'box', 'off')
    hold off;
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})],'Resolution',300); 
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.tiff',meas_names{meas_i})],'Resolution',900); 
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.pdf',meas_names{meas_i})],'ContentType','vector','Resolution',300); 
    
    %% (PLOT) Trial & Subject Category Plot for High vs Low function OA
    vals = cat(1,cond_1{1,:});
    bandwidth = range(vals)*0.1;
    for g_i = 1:size(cond_1,2)
        n=max(cellfun(@numel,cond_1(:,g_i)));
        for c_i = 1:size(cond_1,1)
    %         tmp = cond_1{c_i,:};
            cond_1{c_i,g_i}(end+1:n,:) = nan;
        end
    end
    varargin = {'width',VIOLIN_WIDTH_GROUP,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',bandwidth,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter'};
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    c_maps = linspecer(size(cond_1,1));
    offset = 0;
    cnt = 1;
    cnt_g = 1;
    xticks = [];
    xtick_labs = {};
%     plot_pos = [0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25];
    g1 = strsplit(char(g_cats(1)),'_');
    g2 = strsplit(char(g_cats(2)),'_');
    g_labs = [{g1{1}},{g2{1}}];
    %-
%     save_lab = 'cond_main';
%     for c_i=1:size(cond_1,1)
%         for g_i = 1:size(cond_1,2)
%             if g_i == 1
%                 offset = -0.15;
%             else
%                 offset = 0.15;
%             end
%             violins(cnt) = Violin(cond_1(c_i,g_i),c_i,varargin{:});
%             violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
%             violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
%             violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
%             violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
%             violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
%             violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
%             violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
%             violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
%             violins(cnt).ViolinColor = {c_maps(c_i,:)};
%             xticks = [xticks,c_i+offset];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_speed{c_i})}];
%             cnt = cnt+1;
%         end
%     end
    %-
    save_lab = 'group_main';
    cond_offsets = [-0.3,-0.1,0.1,0.3];
    for g_i = 1:size(cond_1,2)
        for c_i=1:size(cond_1,1)
            offset = cond_offsets(c_i);
            violins(cnt) = Violin(cond_1(c_i,g_i),g_i,varargin{:});
            violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
            violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
            violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
            violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
            violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
            violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
            violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
            violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
            violins(cnt).ViolinColor = {COLOR_MAPS_SPEED(c_i,:)};
            xticks = [xticks,g_i+offset];
            xtick_labs = [xtick_labs,{sprintf('%s',trial_names_speed{c_i})}];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_speed{c_i})}];
            g_ind = (char(g_cats(g_i))==outliers.SubjectCategory) & (char(speed_chars(c_i))==outliers.TrialName);
%             tmp_outliers{g_i,c_i} = outliers.(meas_names{meas_i})(g_ind);
            vals_y = outliers.(meas_names{meas_i})(g_ind);
            if ~isempty(vals_y)
                in_x = repmat(xticks(cnt),size(vals_y));
                scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
            end
            cnt = cnt+1;
        end
    end
    %- set figure color, units, and size
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    %- set axes units, and size
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',10,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.20 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
    set(ax,'XTick',sort(xticks));
    %-
%     shift = 0;
%     mdl_spec='Var1~1+Var2';
%     for g_i=1:size(cond_1,2)
%         tmp = cat(2,cond_1{:,g_i});
%         for subj_i = 1:size(cond_1{1,1},1)
%             if all(~isnan(tmp(subj_i,:)))
%                 y_vals = tmp(subj_i,:);
%                 x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
%                 tb = table(y_vals',x_vals');
%                 out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
% %                 p = plot(ax,out.Residuals.Raw',x_vals);
%                 p = plot(ax,out);
%                  p(end-1,1).Visible='off';
%                 p(end,1).Visible='off';
%                 p(1).Visible = 'off'; %[0,0,0,0.2];
%                 if out.Coefficients.Estimate(2) > 0
%                     p(2).Color = [0,0,0.7,0.70];
%                 else
%                     p(2).Color = [0.7,0,0,0.70];
%                 end
%                 
% %                 plot(ax,x_vals,y_vals);
%             end
%         end
%         shift = shift + length(tmp(1,:));
%     end
    set(ax,'XTickLabel',xtick_labs);
    xtickangle(30);
    xlh = xlabel('Speed (m/s)','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.13;
    set(xlh,'Position',pos1);
%     xticks = get(ax,'XTick');
    %- set group labels
    if size(cond_1,2) == 2
        shift = 0;
        for g_i = 1:size(cond_1,2)
            text(0.25+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.50;
        end
    elseif size(cond_1,2) == 3
        shift = 0;
        for g_i = 1:size(cond_1,2)
            text(0.15+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.34;
        end
    else
        error('here');
    end
    %- set ylabel & title
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})]);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.tiff',save_lab,meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'ContentType','vector','Resolution',300); 
    %% (PLOT) Trial & Subject Category Plot for High vs Low function OA
    vals = cat(1,cond_2{1,:});
    bandwidth = range(vals)*0.1;
    for g_i = 1:size(cond_2,2)
        n=max(cellfun(@numel,cond_2(:,g_i)));
        for c_i = 1:size(cond_2,1)
    %         tmp = cond_2{c_i,:};
            cond_2{c_i,g_i}(end+1:n,:) = nan;
        end
    end
    varargin = {'width',VIOLIN_WIDTH_GROUP,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',bandwidth,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter'};
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    c_maps = linspecer(size(cond_2,1));
    offset = 0;
    cnt = 1;
    cnt_g = 1;
    xticks = [];
    xtick_labs = {};
%     plot_pos = [0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25];
    g1 = strsplit(char(g_cats(1)),'_');
    g2 = strsplit(char(g_cats(2)),'_');
    g_labs = [{g1{1}},{g2{1}}];
    %- cond priority, group secondary
%     save_lab = 'cond_main';
%     for c_i=1:size(cond_2,1)
%         for g_i = 1:size(cond_2,2)
%             if g_i == 1
%                 offset = -0.15;
%             else
%                 offset = 0.15;
%             end
%             violins(cnt) = Violin(cond_2(c_i,g_i),c_i,varargin{:});
%             violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
%             violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
%             violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
%             violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
%             violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
%             violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
%             violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
%             violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
%             violins(cnt).ViolinColor = {c_maps(c_i,:)};
%             xticks = [xticks,c_i+offset];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_terrain{c_i})}];
%             cnt = cnt+1;
%         end
%     end
    %- group primary, cond secondary
    save_lab = 'group_main';
    cond_offsets = [-0.3,-0.1,0.1,0.3];
    %-
    %{
    figure;
    hold on;
    shift = 0;
    for g_i=1:size(cond_2,2)
        tmp = cat(2,cond_2{:,g_i});
        for subj_i = 1:size(cond_2{1,1},1)
            if ~isnan(tmp(subj_i,1))
                x_vals = (1:length(tmp(subj_i,:)))+shift;
                plot(x_vals,tmp(subj_i,:));
            end
        end
        shift = shift + length(tmp(1,:));
    end
    hold off;
    %}
%     vals_out = [];
%     for g_i=1:size(cond_2,2)
%         tmp = cat(2,cond_2{:,g_i});
%         vals_out = [vals_out,tmp];
%     end
    %-
    for g_i = 1:size(cond_2,2)
        for c_i=1:size(cond_2,1)
            offset = cond_offsets(c_i);
            violins(cnt) = Violin(cond_2(c_i,g_i),g_i,varargin{:});
            violins(cnt).ScatterPlot.XData = violins(cnt).ScatterPlot.XData+offset;
            violins(cnt).ViolinPlot.XData = violins(cnt).ViolinPlot.XData+offset;
            violins(cnt).WhiskerPlot.XData = violins(cnt).WhiskerPlot.XData+offset;
            violins(cnt).MedianPlot.XData = violins(cnt).MedianPlot.XData+offset;
            violins(cnt).NotchPlots(1).XData = violins(cnt).NotchPlots(1).XData+offset;
            violins(cnt).NotchPlots(2).XData = violins(cnt).NotchPlots(2).XData+offset;
            violins(cnt).MeanPlot.XData = violins(cnt).MeanPlot.XData+offset;
            violins(cnt).ViolinPlotQ.XData = violins(cnt).ViolinPlotQ.XData+offset;
            violins(cnt).ViolinColor = {COLORS_MAPS_TERRAIN(c_i,:)};
            xticks = [xticks,g_i+offset];
            xtick_labs = [xtick_labs,{sprintf('%s',trial_names_terrain{c_i})}];
%             xtick_labs = [xtick_labs,{sprintf('%s %s',g_labs{g_i},trial_names_terrain{c_i})}];
            g_ind = (char(g_cats(g_i))==outliers.SubjectCategory) & (char(terrain_chars(c_i))==outliers.TrialName);
%             tmp_outliers{g_i,c_i} = outliers.(meas_names{meas_i})(g_ind);
            vals_y = outliers.(meas_names{meas_i})(g_ind);
            if ~isempty(vals_y)
                in_x = repmat(xticks(cnt),size(vals_y));
                scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
            end
            cnt = cnt+1;
        end
%         text(nanmean(xticks(cnt_g)),-8,char(group_names(g_i)),'FontSize',11);
%         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
%         cnt_g = cnt_g+cnt-1;
    end
    %- set figure color, units, and size
    fig_i = get(groot,'CurrentFigure');
    
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 12;
    %- set axes units, and size
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',11,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.20 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
    set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',xtick_labs);
    xtickangle(30);
    %-
%     xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
%     pos1=get(xlh,'Position');
%     pos1(1,2)=pos1(1,2)-0.12;
%     set(xlh,'Position',pos1);
    %-
%     shift = 0;
%     mdl_spec='Var1~1+Var2';
%     for g_i=1:size(cond_2,2)
%         tmp = cat(2,cond_2{:,g_i});
%         for subj_i = 1:size(cond_2{1,1},1)
%             if all(~isnan(tmp(subj_i,:)))
%                 y_vals = tmp(subj_i,:);
%                 x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
%                 tb = table(y_vals',x_vals');
%                 out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
% %                 p = plot(ax,out.Residuals.Raw',x_vals);
%                 p = plot(ax,out);
%                 p(end-1,1).Visible='off';
%                 p(end,1).Visible='off';
%                 p(1).Visible = 'off'; %[0,0,0,0.2];
%                 if out.Coefficients.Estimate(2) > 0
%                     p(2).Color = [0,0,0.7,0.70];
%                 else
%                     p(2).Color = [0.7,0,0,0.70];
%                 end
% 
% %                 plot(ax,x_vals,y_vals);
%             end
%         end
%         shift = shift + length(tmp(1,:));
%     end
    %- set group labels
    if size(cond_2,2) == 2
        shift = 0;
        for g_i = 1:size(cond_2,2)
            text(0.25+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.50;
        end
    elseif size(cond_2,2) == 3
        shift = 0;
        for g_i = 1:size(cond_2,2)
            text(0.15+shift,-0.15,0,char(group_names(g_i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
                'Units','normalized');
    %         text(nanmean(xticks(cnt_g))-0.20,-17,char(group_names(g_i)),'FontSize',11);
    %         disp(xticks(cnt_g))
            shift = shift+0.34;
        end
    else
        error('here');
    end
    %- set ylabel & title
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})]);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.jpg',save_lab,meas_names{meas_i})],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.tiff',save_lab,meas_names{meas_i})],'Resolution',900);
%     exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'ContentType','vector','Resolution',300); 
    
end
%% OLD CODE SNIPETS
%{
%% STATISTITCS
    %% (SPEED) CONDITION EFFECT ONLY
    modelspec = sprintf('%s~1+TrialName',meas_names{meas_i});
    tmp_t = [];
    for i = 1:length(speed_chars)
        tmp_t = [tmp_t; table_new_ls(table_new_ls.TrialName == speed_chars{i},:)];
    end
    mdl_speed = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar',meas_names{meas_i},...
        'PredictorVar',{'TrialName','SubjectCategory'});
    disp(mdl_speed)
    %- Convert summary to char array
    txt = evalc('mdl_speed');
    fid = fopen([save_dir filesep sprintf('%s_across_speed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    fclose(fid);
    %% (TERRAIN) CONDITION EFFECT ONLY
    modelspec = sprintf('%s~1+TrialName',meas_names{meas_i});
    tmp_t = [];
    for i = 1:length(terrain_chars)
        tmp_t = [tmp_t; table_new_ls(table_new_ls.TrialName == terrain_chars{i},:)];
    end
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar',meas_names{meas_i},...
        'PredictorVar',{'TrialName','SubjectCategory'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_terrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    fclose(fid);
    %% GROUP EFFECT ONLY
    modelspec = sprintf('%s~1+SubjectCategory',meas_names{meas_i});
    tmp_t = [];
    for i = 1:length(terrain_chars)
        tmp_t = [tmp_t; table_new_ls(table_new_ls.TrialName == terrain_chars{i},:)];
    end
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar',meas_names{meas_i},...
        'PredictorVar',{'TrialName','SubjectCategory'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_subjectcats_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    fclose(fid);
    %% (TERRAIN,NONMIXED) GROUP & CONDITION EFFECT 
    modelspec = sprintf('%s~1+TrialName+SubjectCategory',meas_names{meas_i});
%     modelspec = sprintf('%s~1+TrialName+SubjectCategory+TrialName:SubjectCategory',meas_names{meas_i});
    tmp_t = [];
    for i = 1:length(terrain_chars)
        tmp_t = [tmp_t; table_new_ls(table_new_ls.TrialName == terrain_chars{i},:)];
    end
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar',meas_names{meas_i},...
        'PredictorVar',{'TrialName','SubjectCategory'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_nonmixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    fclose(fid);
    %% (TERRAIN,MIXED) GROUP & CONDITION EFFECT
    inds = [];
    modelspec = 'meas_in~(1+cat_1)+(1+cat_2)+cat_1:cat_2';
    for i = 1:length(terrain_chars)
        inds = [inds; find(table_new_ls.TrialName == terrain_chars{i})];
    end
    tmp = table_new_ls.(meas_names{meas_i});
    meas_in = double(tmp(inds,:));
    tmp = table_new_ls.TrialName;
%     cat_1 = categorical(tmp(inds,:));
    [vals_tn,~,cati] = unique(tmp(inds,:));
%     for i = 1:length(vals)
%         fprintf('%i: %s\n',i,vals(i));
%     end
    cat_1 = categorical(cati);
    tmp = table_new_ls.SubjectCategory;
    cat_2 = categorical(tmp(inds,:));
    tmp_t = table(meas_in,cat_1,cat_2);
    %-
%     X = [dummyvar(tmp_t.cat_1), dummyvar(tmp_t.cat_2)]; % DummyVarCoding -> full
%     disp(rank(X)) % 3 < size(X, 2) --> 3 < 4  --> rank deficient
%     [~, R] = qr(X, 0);
%     find(abs(diag(R)) < 1e-6)
    %-
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    %- need econometrics toolbox for this;
%     r = mdl_terrain.Coefficients;
%     R = zeros(size(mdl_terrain.CoefficientCovariance));
%     EstCov = mdl_terrain.CoefficientCovariance;
%     alpha = 0.05;
%     [h,pValue,stat,cValue] = waldtest(r,R,EstCov,alpha);
%     model2 = fitlm(tmp_t, modelspec);
    
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## (TERRAIN,MIXED) GROUP & CONDITION EFFECT
% %     modelspec = sprintf('%s~1+TrialName+SubjectCategory',meas_names{meas_i});
%     modelspec = sprintf('%s~1+TrialName*SubjectCategory',meas_names{meas_i});
%     tmp_t = [];
%     for i = 1:length(terrain_chars)
%         tmp_t = [tmp_t; table_new_ls(table_new_ls.TrialName == terrain_chars{i},:)];
%     end
%     mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar',meas_names{meas_i},...
%         'PredictorVar',{'TrialName','SubjectCategory'});
%     disp(mdl_terrain)
%     %- Convert summary to char array
%     txt = evalc('mdl_terrain');
%     fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
%     fprintf(fid,'%s',txt);
%     fclose(fid);
%}