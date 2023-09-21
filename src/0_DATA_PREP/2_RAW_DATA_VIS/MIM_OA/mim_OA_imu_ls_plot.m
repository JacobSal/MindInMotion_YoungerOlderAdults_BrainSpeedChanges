%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_proc/run_conn_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'MIM_OA_YA'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (DATASET INFORMATION) =============================================== %%
%## (MIND IN MOTION) DATASET SPECIFIC PARAMS (05/24/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H1008','H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
SUBJ_NO_MRI = {'H2010', 'H2036', 'H2041', 'H2072', 'H3018','H3120'};
SUBJ_1YA = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1032','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047'}; % JACOB,SAL (04/18/2023)
% SUBJ_2MA = {'H2017', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
%     'H2020', 'H2021', 'H2022', 'H2023',...
%     'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2037', 'H2038',...
%     'H2039', 'H2042', 'H2052', 'H2059', 'H2062', 'H2082',...
%     'H2090', 'H2095', 'H2111', 'H2117'};
% SUBJ_3MA = {'H3029','H3034','H3039','H3042','H3046',...
%     'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107',...
%     'NH3006', 'NH3007', 'NH3008', 'NH3010',...
%     'NH3021', 'NH3025', 'NH3026',...
%     'NH3030', 'NH3036',...
%     'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
%     'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
%     'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
%     'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
% SUBJ_DEBUG = {'H2117','NH3082','H3063','NH3006','NH3025','NH3114','H2007',...
%     'H3034','NH3055','H3073','NH3104','NH3051','NH3123','H3092','NH3082',...
%     'NH3056','NH3036','H3046','H3053','NH3007','H3077','H3047','NH3071'};
SUBJ_2MA = {'H2002','H2007','H2008',...
    'H2013','H2015','H2017','H2020','H2021',...
    'H2O22','H2023','H2025','H2026','H2027',...
    'H2033','H2034','H2037','H2038','H2039',...
    'H2042','H2052','H2059','H2062','H2082',...
    'H2090','H2095','H2111','H2117'};
SUBJ_3MA = {'H3024','H3029','H3034','H3039','H3042','H3046','H3047','H3053',...
    'H3063','H3072','H3073','H3077','H3092','H3103',...
    'H3107',...
    'NH3006','NH3007','NH3008','NH3010','NH3021',...
    'NH3025','NH3026','NH3028','NH3030','NH3036',...
    'NH3040','NH3041','NH3043','NH3051','NH3054',...
    'NH3055','NH3056','NH3058','NH3059','NH3066',...
    'NH3068','NH3069','NH3070','NH3071','NH3074',...
    'NH3076','NH3082','NH3086','NH3090','NH3102',...
    'NH3104','NH3105','NH3106','NH3108','NH3110',...
    'NH3112','NH3113','NH3114','NH3123','NH3128',...
    };
SUBJ_NO_MRI = {'H2010','H2012','H2018','H2036','H2041',...
    'H2072','H3018','H3120','NH3002','NH3009','NH3027','NH3129'};
SUBJ_DONT_INC = {'NH3004','NH3023'};
% (08/20/2023) JS, NH3004 has no headscan; NH3023 has no headscan clicks;
%- (YA,OA) Subject Picks 
SUBJ_PICS = {SUBJ_1YA,SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_1YA),1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- OA Subject Picks
% SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H2000''s','H3000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- test
% SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H2000''s','H3000''s'}; 
% SUBJ_ITERS = {[1,2],[1,2]};
fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
colormap(linspecer);
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
%- Subject Directory information
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'raw_data_vis'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%%
CATEGORIES = {'YoungAdult','HF_OlderAdult','LF_OlderAdult'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
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
subj_stack = [];
cnt_ga = 1;
cnt_ls = 1;
%%
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        %## LOAD CONDITION WISE BEHAVIORS
        %{
        folder_imu = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Conditions'];
        folder_ls = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Conditions'];
        dir_imu = dir([folder_imu filesep '*.mat']);
        dir_ls = dir([folder_ls filesep '*.mat']);
        %}
        %## LOAD TRIAL WISE BEHAVIORS
        folder_imu = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Trials'];
        folder_ls = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Trials'];
        dir_imu = dir([folder_imu filesep 'Outcome_Measures' filesep '*.mat']);
        dir_ls = dir([folder_ls filesep 'Outcome_Measures' filesep '*.mat']);
        fprintf('Loading subject %s...\n',SUBJ_PICS{group_i}{subj_i});
        vals = zeros(length(dir_imu),2);
        %## IMU (SACRAL MEASURES ONLY) 
        for f_i = 1:length(dir_imu)
            tmp = load([dir_imu(f_i).folder filesep dir_imu(f_i).name]);
            tmp = tmp.myStructure;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names{cnt_ls} = fieldnames(tmp);
            table_imu_meas(cnt_ls,1:length(values)) = values;
            tmp = strsplit(dir_imu(f_i).name,'.');
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
   
    end
end
%% IMU TABLE
table_imu_meas = table_imu_meas(~any(table_imu_meas == 0,2),:);
table_subj_vec = table_subj_vec(~any(table_imu_meas == 0,2));
table_trial_vec = table_trial_vec(~any(table_imu_meas == 0,2));
table_subj_cat_vec = table_subj_cat_vec(~any(table_imu_meas == 0,2));
table_imu_out = array2table(table_imu_meas,'VariableNames',table_header_names{1});
table_imu_out.SubjectName = categorical(table_subj_vec);
table_imu_out.TrialName = categorical(table_trial_vec);
table_imu_out.SubjectCategory = categorical(table_subj_cat_vec);
%- rearrange headers
% table_imu_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
% table_imu_out = table(categorical(table_subj_vec),categorical(table_trial_vec),table_imu_meas(:,1),table_imu_meas(:,2),...
%     table_imu_meas(:,3),table_imu_meas(:,4),table_imu_meas(:,5),table_imu_meas(:,6),...
%     table_imu_meas(:,7),table_imu_meas(:,8),table_imu_meas(:,9),table_imu_meas(:,10),...
%     table_imu_meas(:,11),table_imu_meas(:,12),'VariableNames',[{'Subject'},{'TrialName'},table_header_names{1}']);
%% LOADSOL TABLE
table_ls_meas = table_ls_meas(~any(table_ls_meas == 0,2),:);
table_subj_ls = table_subj_ls(~any(table_ls_meas == 0,2));
table_trial_ls = table_trial_ls(~any(table_ls_meas == 0,2));
table_subj_cat_ls = table_subj_cat_ls(~any(table_ls_meas == 0,2));
table_ls_out = array2table(table_ls_meas,'VariableNames',table_header_names_ls{1});
table_ls_out.SubjectName = categorical(table_subj_ls);
table_ls_out.TrialName = categorical(table_trial_ls);
table_ls_out.SubjectCategory = categorical(table_subj_cat_ls);
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
        tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
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
        tmp = varfun(@nanmean, tmp, 'InputVariables', @isnumeric);
        tmp.SubjectName = subjs(i);
        tmp.TrialName = trials(j);
        tmp.SubjectCategory = subj_cat;
        table_new = [table_new; tmp];
    end
end
table_new_ls = table_new;
%% VIOLIN PLOT IMU
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,420,420];
FIG_POSITION_GRP = [100,100,720,420];
VIOLIN_WIDTH = 0.3;
VIOLIN_WIDTH_GROUP = 0.1;
% meas_names = {'APexc_mean','APexc_COV','MLexc_mean','MLexc_COV','UDexc_mean',...
%     'UDexc_COV','PeakAntVel_mean','PeakAntVel_COV','PeakLatVel_mean',...
%     'PeakLatVel_COV','PeakUpDownVel_mean','PeakUpDownVel_COV'};
% meas_units = {'m','m','m','m','m',...
%     'm','m/s','m/s','m/s',...
%     'm/s','m/s','m/s'};
meas_names = {'nanmean_APexc_mean','nanmean_MLexc_mean','nanmean_APexc_COV','nanmean_MLexc_COV'}; %{'APexc_COV','MLexc_COV'};
meas_units = {'m','m','%','%'};
meas_titles = {'Anteriorposterior Excursion Mean','Mediolateral Excursion Mean','Anteriorposterior Excursion Coefficient of Variation','Mediolateral Excursion Coefficient of Variation'};
meas_ylabel = {'Distance (m)','Distance (m)','Coefficient of Variation','Coefficient of Variation'};
group_labels = {};
YLIMS = {[0,0.2],[0,0.2],[0,55],[0,40]};
% trial_names = {'SP_0p25_1','SP_0p25_2','SP_0p5_1','SP_0p5_2','SP_0p75_1',...
%                 'SP_0p75_2','SP_1p0_1','SP_1p0_2','TM_flat_1','TM_flat_2',...
%                 'TM_high_1','TM_high_2','TM_low_1','TM_low_2','TM_med_1',...
%                 'TM_med_2'};
% trial_names_speed = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
%                 '0.75 m/s','1.00 m/s','1.00 m/s'};
% trial_names_terrain = {'flat terr.','flat terr.',...
%     'low terr.','low terr.','med terr.',...
%     'med terr.','high terr.','high terr.'};
speed_chars = {'0p25','0p5','0p75','1p0'};
terrain_chars = {'flat','low','med','high'};
% trial_names_speed = {'0.25 m/s','0.50 m/s','0.75 m/s',...
%                '1.00 m/s'};
% trial_names_terrain = {'flat terr.',...
%     'low terr.','med terr.',...
%     'high terr.'};
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
COLORS_MAPS = linspecer(size(speed_chars,2));
custom_yellow = [254,223,0]/255;
COLORS_MAPS = [COLORS_MAPS(3,:);custom_yellow;COLORS_MAPS(4,:);COLORS_MAPS(2,:)];
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
    b_speed = tmp(:,2:5);
    b_terrain = [tmp(:,6),tmp(:,8),tmp(:,9),tmp(:,7)];
    %## By subject plot
%     figure;
%     title(sprintf('%s Across Subjects',meas_names{meas_i}));
%     hold on;
%     violinplot(a(:,13:end),a(:,13),...
%         'ViolinColor',linspecer(size(a(:,13:end),2)),...
%         'ShowWhiskers',false,'ShowNotches',true);
%     ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
%     xlabel('Subject Code');
%     xtickangle(45)
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
% %     saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',meas_names{meas_i})]);

    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_terrain,trial_names_terrain','width',VIOLIN_WIDTH,...
        'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',COLORS_MAPS,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(terrain_chars)
        idx = outliers.TrialName == terrain_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
        end
    end
%     ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
%     xlabel('Trial Code');
%     title(meas_titles{meas_i});
%     xtickangle(45)
%     ylim(YLIMS{meas_i});
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
%     fig_i.Children.FontSize = 13;
%     fig_i.Children.XTickLabel = trial_names_terrain;
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',12,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_terrain);
    xtickangle(30);
    xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.pdf',meas_names{meas_i})],'Resolution',300);
    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_speed,trial_names_speed,'width',VIOLIN_WIDTH,...
        'GroupOrder',b_speed.Properties.VariableNames,'ViolinColor',COLORS_MAPS,...
        'ShowWhiskers',true,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_speed{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(speed_chars)
        idx = outliers.TrialName == speed_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
        end
    end
    
    %- consider a colormap that matches for each condition across
    %trials....
%     ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
%     xlabel('Trial Code');
%     xtickangle(45);
%     title(meas_titles{meas_i});
%     ylim(YLIMS{meas_i});
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
%     fig_i.Children.FontSize = 13;
%     fig_i.Children.XTickLabel =  trial_names_speed;
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',12,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_speed);
    xtickangle(30);
    xlh = xlabel('Speed (m/s)','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.pdf',meas_names{meas_i})],'Resolution',300);
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
%     plot_pos = [0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25];
%     g1 = strsplit(char(g_cats(1)),'_');
%     g2 = strsplit(char(g_cats(2)),'_');
%     g_labs = [{g1{1}},{g2{1}}];
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
            violins(cnt).ViolinColor = {COLORS_MAPS(c_i,:)};
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
    shift = 0;
    mdl_spec='Var1~1+Var2';
    for g_i=1:size(cond_1,2)
        tmp = cat(2,cond_1{:,g_i});
        for subj_i = 1:size(cond_1{1,1},1)
            if all(~isnan(tmp(subj_i,:)))
                y_vals = tmp(subj_i,:);
                x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
                tb = table(y_vals',x_vals');
                out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
%                 p = plot(ax,out.Residuals.Raw',x_vals);
                p = plot(ax,out);
                 p(end-1,1).Visible='off';
                p(end,1).Visible='off';
                p(1).Visible = 'off'; %[0,0,0,0.2];
                p(2).Color = [0,0,0.7,0.60];
%                 plot(ax,x_vals,y_vals);
            end
        end
        shift = shift + length(tmp(1,:));
    end
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
%     saveas(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.hdf',save_lab,meas_names{meas_i})]);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'Resolution',300);
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
%     c_maps = linspecer(size(cond_2,1));
%     offset = 0;
    cnt = 1;
    cnt_g = 1;
    xticks = [];
    xtick_labs = {};
%     plot_pos = [0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25];
%     g1 = strsplit(char(g_cats(1)),'_');
%     g2 = strsplit(char(g_cats(2)),'_');
%     g_labs = [{g1{1}},{g2{1}}];
    %-
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
            violins(cnt).ViolinColor = {COLORS_MAPS(c_i,:)};
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
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',11,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.20 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
    set(ax,'XTick',sort(xticks));
    %-
    shift = 0;
    mdl_spec='Var1~1+Var2';
    for g_i=1:size(cond_2,2)
        tmp = cat(2,cond_2{:,g_i});
        for subj_i = 1:size(cond_1{1,1},1)
            if all(~isnan(tmp(subj_i,:)))
                y_vals = tmp(subj_i,:);
                x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
                tb = table(y_vals',x_vals');
                out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
%                 p = plot(ax,out.Residuals.Raw',x_vals);
                p = plot(ax,out);
                p(end-1,1).Visible='off';
                p(end,1).Visible='off';
                p(1).Visible = 'off'; %[0,0,0,0.2];
                p(2).Color = [0,0,0.7,0.60];
%                 plot(ax,x_vals,y_vals);
            end
        end
        shift = shift + length(tmp(1,:));
    end
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
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'Resolution',300);
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
    tmp_t = table(meas_in,cat_1,cat_2);
    %## MIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## NONMIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_nonmixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## CONDITION ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_terrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## GROUP ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_subjectcats_terrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
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
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## NONMIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_nonmixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## CONDITION ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_speed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## GROUP ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_subjcats_speed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
end
%% VIOLIN PLOT LOADSOL
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,360,420];
FIG_POSITION_GRP = [100,100,720,420];
VIOLIN_WIDTH = 0.3;
VIOLIN_WIDTH_GROUP = 0.1;
% meas_names = {'GaitCycleDur_L';'GaitCycleDur_R';'GaitCycleDur';'GaitCycleDur_Lcov';...
%               'GaitCycleDur_Rcov';'GaitCycleDur_cov';'GaitCycleDur_asym';'GaitCycleDur_covAsym';...
%               'StanceDur_L';'StanceDur_R';'StanceDur';'StanceDur_Lcov';'StanceDur_Rcov';...
%               'StanceDur_cov';'StanceDur_asym';'StanceDur_covAsym';'SwingDur_L';...
%               'SwingDur_R';'SwingDur';'SwingDur_Lcov';'SwingDur_Rcov';'SwingDur_cov';...
%               'SwingDur_asym';'SwingDur_covAsym';'SingleSupport_L';'SingleSupport_R';...
%               'SingleSupport';'SingleSupport_Lcov';'SingleSupport_Rcov';'SingleSupport_cov';...
%               'SingleSupport_asym';'SingleSupport_covAsym';'StepDur_L';'StepDur_R';...
%               'StepDur';'StepDur_Lcov';'StepDur_Rcov';'StepDur_cov';'StepDur_asym';...
%               'StepDur_covAsym'}; 
% %           ;'InitialDS_L';'InitialDS_R';'InitialDS';'InitialDS_Lcov';...
% %               'InitialDS_Rcov';'InitialDS_cov';'InitialDS_asym';'InitialDS_covAsym';...
% %               'TerminalDS_L';'TerminalDS_R';'TerminalDS';'TerminalDS_Lcov';'TerminalDS_Rcov';...
% %               'TerminalDS_cov';'TerminalDS_asym';'TerminalDS_covAsym';'TotalDS_L';...
% %               'TotalDS_R';'TotalDS';'TotalDS_Lcov';'TotalDS_Rcov';'TotalDS_cov';...
% %               'TotalDS_asym';'TotalDS_covAsym'};
% meas_units = {'s','s','s','s',...
%                's','s','s','s',...
%                's','s','s','s','s',...
%                's','s','s','s',...
%                's','s','s','s','s',...
%                's','s','s','s',...
%                's','s','s','s',...
%                's','s','s','s',...
%                's','s','s','s','s',...
%                's'};
meas_names = {'nanmean_StepDur','nanmean_StepDur_cov','nanmean_GaitCycleDur_cov','nanmean_GaitCycleDur',};
meas_units = {'s','%','%','s'};
meas_titles = {'Step Duration','Step Duration Coefficient of Variation','Gait Cycle Duration Coefficient of Variation','Gait Cycle Duration'};
meas_ylabel = {'Duration','Coefficient of Variation','Coefficient of Variation','Duration'};
% YLIMS = {[0.3,1.3],[0,25],[0,65],[0,4]};
YLIMS = {[0,2],[0,40],[0,65],[0,4]};
% trial_names = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
%                 '0.75 m/s','1.00 m/s','1.00 m/s','flat terr.','flat terr.',...
%                 'high terr.','high terr.','low terr.','low terr.','med terr.',...
%                 'med terr.'};
% trial_names_speed = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
%                 '0.75 m/s','1.00 m/s','1.00 m/s'};
% trial_names_terrain = {'flat terr.','flat terr.',...
%     'low terr.','low terr.','med terr.',...
%     'med terr.','high terr.','high terr.'};
speed_chars = {'0p25','0p5','0p75','1p0'};
terrain_chars = {'flat','low','med','high'};
% trial_names_speed = {'0.25 m/s','0.50 m/s','0.75 m/s',...
%                '1.00 m/s'};
% trial_names_terrain = {'flat terr.',...
%     'low terr.','med terr.',...
%     'high terr.'};
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
COLORS_MAPS = linspecer(size(speed_chars,2));
custom_yellow = [254,223,0]/255;
COLORS_MAPS = [COLORS_MAPS(3,:);custom_yellow;COLORS_MAPS(4,:);COLORS_MAPS(2,:)];
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
    tmp = b(:,65:end);
%     b_speed = tmp(:,2:9);
    b_speed = tmp(:,2:5);
%     b_terrain = [tmp(:,10:11),tmp(:,14:15),tmp(:,16:17),tmp(:,12:13)];
    b_terrain = [tmp(:,6),tmp(:,8),tmp(:,9),tmp(:,7)];
    %## By subject plot
%     figure;
%     title(sprintf('%s Across Subjects',meas_names{meas_i}));
%     hold on;
%     violinplot(a(:,64:end),...
%         'ViolinColor',linspecer(size(a(:,64:end),2)));
%     ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
%     xlabel('Subject Code');
%     xtickangle(45)
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
% %     saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',meas_names{meas_i})]);
    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
%     violinplot(b_terrain,trial_names_terrain','width',0.4,...
%         'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',linspecer(size(b_terrain,2)),...
%         'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
%         'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
%         'HalfViolin','left','DataStyle','histogram');
    violinplot(b_terrain,trial_names_terrain','width',VIOLIN_WIDTH,...
        'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',COLORS_MAPS,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(speed_chars)
        idx = outliers.TrialName == speed_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
        end
    end
    %- consider a colormap that matches for each condition across
    %trials....
%     ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}));
%     xlabel('Trial Code');
%     xtickangle(45);
%     title(meas_titles{meas_i});
%     ylim(YLIMS{meas_i});
%     hold off;
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
%     fig_i.Children.FontSize = 13;
%     fig_i.Children.XTickLabel =  trial_names_speed;
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',12,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_terrain);
    xtickangle(30);
    xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    %- set category labels
%     text(0.45,-0.2,0,'Speeds (m/s)','FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
%         'Units','normalized');
    %- set ylabel
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.pdf',meas_names{meas_i})],'Resolution',300); 
    %% By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_speed,trial_names_speed,'width',VIOLIN_WIDTH,...
        'GroupOrder',b_speed.Properties.VariableNames,'ViolinColor',COLORS_MAPS,...
        'ShowWhiskers',true,'ShowNotches',false,'ShowBox',true,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1,'QuartileStyle','shadow',...
        'HalfViolin','full','DataStyle','scatter');
    vals_x = get(gca,'XTick');
    for c_i = 1:length(terrain_chars)
        idx = outliers.TrialName == terrain_chars{c_i};
        vals_y = outliers.(meas_names{meas_i})(idx);
        if ~isempty(vals_y)
            in_x = repmat(vals_x(c_i),size(vals_y));
            scatter(in_x,vals_y,'k*','jitter','on', 'jitterAmount', 0.05)
        end
    end
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = FIG_POSITION;
%     fig_i.Children.FontSize = 13;
%     fig_i.Children.XTickLabel =  trial_names_terrain;
    %-
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Color','w')
    set(fig_i,'Units','inches','Position',[3 3 6 5])
    set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    fig_i.Children.FontSize = 13;
    ax = gca; %fig_i.CurrentAxes;
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',12,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
%     set(ax,'XTick',sort(xticks));
    set(ax,'XTickLabel',trial_names_speed);
    xtickangle(30);
    xlh = xlabel('Speed (m/s)','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.02;
    set(xlh,'Position',pos1);
    %- set xlabel 2
%     text(0.45,-0.2,0,'Terrain Difficulty','FontSize',11,'FontWeight','bold','HorizontalAlignment','center',...
%         'Units','normalized');
    %- set ylabel
    ylabel(sprintf('%s (%s)',meas_ylabel{meas_i},meas_units{meas_i}),'FontSize',12);
    title(meas_titles{meas_i});
    ylim(YLIMS{meas_i});
    hold off;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
%     saveas(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.pdf',meas_names{meas_i})],'Resolution',300); 
    
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
            violins(cnt).ViolinColor = {COLORS_MAPS(c_i,:)};
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
    shift = 0;
    mdl_spec='Var1~1+Var2';
    for g_i=1:size(cond_1,2)
        tmp = cat(2,cond_1{:,g_i});
        for subj_i = 1:size(cond_1{1,1},1)
            if all(~isnan(tmp(subj_i,:)))
                y_vals = tmp(subj_i,:);
                x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
                tb = table(y_vals',x_vals');
                out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
%                 p = plot(ax,out.Residuals.Raw',x_vals);
                p = plot(ax,out);
                 p(end-1,1).Visible='off';
                p(end,1).Visible='off';
                p(1).Visible = 'off'; %[0,0,0,0.2];
                p(2).Color = [0,0,0.7,0.60];
%                 plot(ax,x_vals,y_vals);
            end
        end
        shift = shift + length(tmp(1,:));
    end
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
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_speed_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'Resolution',300); 
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
            violins(cnt).ViolinColor = {COLORS_MAPS(c_i,:)};
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
    xlh = xlabel('Terrain Difficulty','Units','normalized','FontSize',12);
    pos1=get(xlh,'Position');
    pos1(1,2)=pos1(1,2)-0.12;
    set(xlh,'Position',pos1);
    %-
    shift = 0;
    mdl_spec='Var1~1+Var2';
    for g_i=1:size(cond_2,2)
        tmp = cat(2,cond_2{:,g_i});
        for subj_i = 1:size(cond_2{1,1},1)
            if all(~isnan(tmp(subj_i,:)))
                y_vals = tmp(subj_i,:);
                x_vals = xticks((1:length(tmp(subj_i,:)))+shift);
                tb = table(y_vals',x_vals');
                out = fitlm(x_vals,y_vals); %fitlm(tb,mdl_spec);
%                 p = plot(ax,out.Residuals.Raw',x_vals);
                p = plot(ax,out);
                 p(end-1,1).Visible='off';
                p(end,1).Visible='off';
                p(1).Visible = 'off'; %[0,0,0,0.2];
                p(2).Color = [0,0,0.7,0.60];
%                 plot(ax,x_vals,y_vals);
            end
        end
        shift = shift + length(tmp(1,:));
    end
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
    exportgraphics(fig_i,[save_dir filesep sprintf('Across_terran_TrialsSubjCat_Fig_%s_%s.pdf',save_lab,meas_names{meas_i})],'Resolution',300); 
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
    tmp_t = table(meas_in,cat_1,cat_2);
    %## MIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2+cat_1:cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## NONMIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_nonmixedsubjectcatsNterrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## CONDITION ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_terrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## GROUP ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_subjectcats_terrain_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
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
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'},'Distribution','poisson');
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_mixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## NONMIXED GROUP & COND
    modelspec = 'meas_in~1+cat_1+cat_2';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1','cat_2'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_nonmixedsubjectcatsNspeed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## CONDITION ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_speed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
    %## GROUP ONLY
    modelspec = 'meas_in~1+cat_1';
    mdl_terrain = fitglm(tmp_t,modelspec,'Intercept',true,'ResponseVar','meas_in',...
        'PredictorVar',{'cat_1'});
    disp(mdl_terrain)
    %- Convert summary to char array
    txt = evalc('mdl_terrain');
    fid = fopen([save_dir filesep sprintf('%s_across_subjcats_speed_all_mdl.txt',meas_names{meas_i})],'wt');
    fprintf(fid,'%s',txt);
    for i = 1:length(vals_tn)
        fprintf(fid,'%i: %s\n',i,vals_tn(i));
    end
    fclose(fid);
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