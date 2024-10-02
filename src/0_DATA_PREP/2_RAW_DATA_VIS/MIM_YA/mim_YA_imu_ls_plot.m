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
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
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
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
global SUBMODULES_DIR
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    DO_UNIX = false;
    PATH_EXT = 'M';
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
end
%## DEBUG: PATHROOT OVERRIDE
if DO_UNIX
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src' filesep '_test' filesep 'changL' filesep 'imu_ls_processing']; %- path to setWorkspace
%- define the direcotry for your scripts
run_dir = [source_dir filesep '1_CONVERT_LS_IMU']; %
%- cd to source directory
cd(source_dir)
%- define directory for submodules
SUBMODULES_DIR = [PATH_ROOT filesep REPO_NAME filesep 'submodules'];
%- addpath for local folder
if exist(source_dir,'dir')
    addpath(source_dir);
else
    error('''source_dir'' does not exist. Make one or fix path to use these scripts');
end
if exist(run_dir,'dir')
    addpath(run_dir);
else
    error('''run_dir'' does not exist. Make one or fix path to use these scripts');
end
%- set workspace
setWorkspace
%% PARPOOL SETUP
if DO_UNIX
    eeg_options;
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i',pp.NumWorkers);
    fprintf('Number of threads: %i',pp.NumThreads);
    %- make meta data directory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- Subject Directory Information
% PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % value to help with file looping
% OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % value to help with file looping
% if DO_UNIX
%     OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR);
% else
%     OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR,'M');
% end
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
if exist(STUDIES_DIR,'dir')
    mkdir(STUDIES_DIR);
end
if exist(SUBJINF_DIR,'dir')
    mkdir(SUBJINF_DIR);
end
%## DATASET SPECIFIC
%- MIND IN MOTION
% SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
%     'H1022','H1024','H1026','H1027','H1033','H1034'};
SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)

SUBJ_2HMA = {'H2017', 'H2010', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2036', 'H2037', 'H2038',...
    'H2039', 'H2041', 'H2042', 'H2052', 'H2059', 'H2062', 'H2072', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
% SUBJ_3HMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
%     'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120'}; % JACOB,SAL(02/23/2023)
SUBJ_3NHMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- Subject Picks
SUBJ_PICS = {SUBJ_2HMA,SUBJ_3NHMA};
GROUP_NAMES = {'H2000''s','H3000''s'};
TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
%- Subject Picks
% SUBJ_PICS = {SUBJ_YNG,SUBJ_HMA,SUBJ_HMA2,SUBJ_NHMA};
% SUBJ_PICS = {SUBJ_YNG};

% SUBJ_ITERS = {1:length(SUBJ_YNG),1:length(SUBJ_HMA),1:length(SUBJ_HMA2),1:length(SUBJ_NHMA)};
% SUBJ_ITERS = {1:length(SUBJ_YNG),[],[]};
% SUBJ_ITERS = {[26,27,28,29,30,31],[],[],[]};
%- Subject Directory Information
PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % value to help with file looping
OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % value to help with file looping
if DO_UNIX
    OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR,'dferris');
else
    OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR,'M');
end
% OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping

%% ===================================================================== %%
%## PARAMS
%- datetime override 
% dt = '20122022';
%- hard define
study_fName = sprintf('copy_study');
%- soft define
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep '%s'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%%
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
table_imu_meas = zeros(length([SUBJ_PICS{:}])*16,12); % 12 measures
table_trial_vec = cell(length([SUBJ_PICS{:}])*16,1);
table_subj_vec = cell(length([SUBJ_PICS{:}])*16,1);
table_header_names = cell(length([SUBJ_PICS{:}]),1);

table_ls_meas = zeros(length([SUBJ_PICS{:}])*16,64); % 64 measures
table_trial_ls = cell(length([SUBJ_PICS{:}])*16,1);
table_subj_ls = cell(length([SUBJ_PICS{:}])*16,1);
table_header_names_ls = cell(length([SUBJ_PICS{:}]),1);
subj_stack = [];
cnt = 1;
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        %{
        folder_imu = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Conditions'];
        folder_ls = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Conditions'];
        dir_imu = dir([folder_imu filesep '*.mat']);
        dir_ls = dir([folder_ls filesep '*.mat']);
        %}
        folder_imu = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'IMU' filesep 'Processed_Trials'];
        folder_ls = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'Loadsol' filesep 'Processed_Trials'];
        dir_imu = dir([folder_imu filesep 'Outcome_Measures' filesep '*.mat']);
        dir_ls = dir([folder_ls filesep 'Outcome_Measures' filesep '*.mat']);
        %## IMU (SACRAL MEASURES ONLY) 
        
        for f_i = 1:length(dir_imu)
            tmp = load([dir_imu(f_i).folder filesep dir_imu(f_i).name]);
            tmp = tmp.myStructure;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names{cnt} = fieldnames(tmp);
            table_imu_meas(cnt,1:length(values)) = values;
            tmp = strsplit(dir_imu(f_i).name,'.');
            table_trial_vec{cnt} = tmp{1};
            table_subj_vec{cnt} = SUBJ_PICS{group_i}{subj_i};
            cnt = cnt + 1;
        end
        
        %## LOADSOL (GAIT)
        fprintf('Loading subject %s...\n',SUBJ_PICS{group_i}{subj_i});
        for f_i = 1:length(dir_ls)
            tmp = load([dir_ls(f_i).folder filesep dir_ls(f_i).name]);
            tmp = tmp.summarizedTrial;
            values = cellfun(@(x)(tmp.(x)),fieldnames(tmp));
            table_header_names_ls{cnt} = fieldnames(tmp);
            table_ls_meas(cnt,1:length(values)) = values;
            tmp = strsplit(dir_ls(f_i).name,'.');
            table_trial_ls{cnt} = tmp{1};
            table_subj_ls{cnt} = SUBJ_PICS{group_i}{subj_i};
            cnt = cnt + 1;
        end
%         cnt = cnt + 1;
   
    end
end
%% IMU TABLE
% table_imu_meas = table_imu_meas(~any(table_imu_meas == 0,2),:);
% table_subj_vec = table_subj_vec(~any(table_imu_meas == 0,2));
% table_trial_vec = table_trial_vec(~any(table_imu_meas == 0,2));
idx = ~cellfun(@isempty,table_subj_vec);
table_imu_meas = table_imu_meas(idx,:);
table_trial_vec = table_trial_vec(idx);
table_subj_vec = table_subj_vec(idx);
table_imu_out = array2table(table_imu_meas,'VariableNames',table_header_names{1});
table_imu_out.SubjectName = categorical(table_subj_vec);
table_imu_out.TrialName = categorical(table_trial_vec);
% table_imu_out = table(categorical(table_subj_vec),categorical(table_trial_vec),table_imu_meas(:,1),table_imu_meas(:,2),...
%     table_imu_meas(:,3),table_imu_meas(:,4),table_imu_meas(:,5),table_imu_meas(:,6),...
%     table_imu_meas(:,7),table_imu_meas(:,8),table_imu_meas(:,9),table_imu_meas(:,10),...
%     table_imu_meas(:,11),table_imu_meas(:,12),'VariableNames',{'Subject','TrialName',table_header_names{1}{:}});
%% LOADSOL TABLE
% table_ls_meas = table_ls_meas(~any(table_ls_meas == 0,2),:);
% table_subj_ls = table_subj_ls(~any(table_ls_meas == 0,2));
% table_trial_ls = table_trial_ls(~any(table_ls_meas == 0,2));
idx = ~cellfun(@isempty,table_subj_ls);
table_ls_meas = table_ls_meas(idx,:);
table_trial_ls = table_trial_ls(idx);
table_subj_ls = table_subj_ls(idx);
table_ls_out = array2table(table_ls_meas,'VariableNames',table_header_names_ls{1});
table_ls_out.SubjectName = categorical(table_subj_ls);
table_ls_out.TrialName = categorical(table_trial_ls);
% table_ls_out = table(categorical(table_subj_ls),categorical(table_trial_ls),table_ls_meas(:,1),table_ls_meas(:,2),...
%     table_ls_meas(:,3),table_ls_meas(:,4),table_ls_meas(:,5),table_ls_meas(:,6),...
%     table_ls_meas(:,7),table_ls_meas(:,8),table_ls_meas(:,9),table_ls_meas(:,10),...
%     table_ls_meas(:,11),table_ls_meas(:,12),'VariableNames',{'Subject','TrialName',table_header_names_ls{1}{:}});
% table_ls_out = table(categorical(table_subj_vec),categorical(table_trial_ls),table_ls_meas(:,1),'VariableNames',{'Subject','TrialName',table_header_names_ls{1}{:}});


%% VIOLIN PLOT IMU
meas_names = {'APexc_mean','APexc_COV','MLexc_mean','MLexc_COV','UDexc_mean',...
    'UDexc_COV','PeakAntVel_mean','PeakAntVel_COV','PeakLatVel_mean',...
    'PeakLatVel_COV','PeakUpDownVel_mean','PeakUpDownVel_COV'};
% trial_names = {'SP_0p25_1','SP_0p25_2','SP_0p5_1','SP_0p5_2','SP_0p75_1',...
%                 'SP_0p75_2','SP_1p0_1','SP_1p0_2','TM_flat_1','TM_flat_2',...
%                 'TM_high_1','TM_high_2','TM_low_1','TM_low_2','TM_med_1',...
%                 'TM_med_2'};
trial_names = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
                '0.75 m/s','1.00 m/s','1.00 m/s','flat terr.','flat terr.',...
                'high terr.','high terr.','low terr.','low terr.','med terr.',...
                'med terr.'};
table_in = table_imu_out; %table_out(:,2:3);
% table_in = unstack(table_in,'TrialName','APexc_mean','GroupingVariables','Subject');
% table_in = unstack(table_in,{'Subject','TrialName'},'APexc_mean');
measure_name = 'APexc_mean';
a = unstack(table_in,measure_name,'SubjectName');
b = unstack(table_in,measure_name,'TrialName');

% table_in = table_out.(meas_names{1});
%## By subject plot
%{
figure;
title('Mean AP Excursion Across Subjects');
hold on;
violinplot(a(:,13:end))
ylabel('Mean AP Excursion IMU (m)');
xlabel('Subject Code');
hold off;
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',measure_name)]);
saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',measure_name)]);
%}
%## By trial plot
figure;
title('Mean AP Excursion Across Trials');
hold on;
violinplot([b{:,13:end}],trial_names,'width',0.4)
ylabel('Mean AP Excursion IMU (m)');
xlabel('Trial Code');
hold off;
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',measure_name)]);
saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.jpg',measure_name)]);

%% VIOLIN PLOT LOADSOL
meas_names = {'GaitCycleDur_L';'GaitCycleDur_R';'GaitCycleDur';'GaitCycleDur_Lcov';...
              'GaitCycleDur_Rcov';'GaitCycleDur_cov';'GaitCycleDur_asym';'GaitCycleDur_covAsym';...
              'StanceDur_L';'StanceDur_R';'StanceDur';'StanceDur_Lcov';'StanceDur_Rcov';...
              'StanceDur_cov';'StanceDur_asym';'StanceDur_covAsym';'SwingDur_L';...
              'SwingDur_R';'SwingDur';'SwingDur_Lcov';'SwingDur_Rcov';'SwingDur_cov';...
              'SwingDur_asym';'SwingDur_covAsym';'SingleSupport_L';'SingleSupport_R';...
              'SingleSupport';'SingleSupport_Lcov';'SingleSupport_Rcov';'SingleSupport_cov';...
              'SingleSupport_asym';'SingleSupport_covAsym';'StepDur_L';'StepDur_R';...
              'StepDur';'StepDur_Lcov';'StepDur_Rcov';'StepDur_cov';'StepDur_asym';...
              'StepDur_covAsym';'InitialDS_L';'InitialDS_R';'InitialDS';'InitialDS_Lcov';...
              'InitialDS_Rcov';'InitialDS_cov';'InitialDS_asym';'InitialDS_covAsym';...
              'TerminalDS_L';'TerminalDS_R';'TerminalDS';'TerminalDS_Lcov';'TerminalDS_Rcov';...
              'TerminalDS_cov';'TerminalDS_asym';'TerminalDS_covAsym';'TotalDS_L';...
              'TotalDS_R';'TotalDS';'TotalDS_Lcov';'TotalDS_Rcov';'TotalDS_cov';...
              'TotalDS_asym';'TotalDS_covAsym'};
trial_names = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
                '0.75 m/s','1.00 m/s','1.00 m/s','flat terr.','flat terr.',...
                'high terr.','high terr.','low terr.','low terr.','med terr.',...
                'med terr.'};
GROUP_NUM = 1;
table_in = table_ls_out; %table_out(:,2:3);
% table_in = unstack(table_in,'TrialName','APexc_mean','GroupingVariables','Subject');
% table_in = unstack(table_in,{'Subject','TrialName'},'APexc_mean');
measure_name = 'GaitCycleDur';
plot_name = measure_name;
a = unstack(table_in,measure_name,'SubjectName');
b = unstack(table_in,measure_name,'TrialName');

% table_in = table_out.(meas_names{1});
%## By subject plot
%{
figure;
title([plot_name ' Across Subjects']);
hold on;
violinplot(a(:,65:end),SUBJ_PICS{GROUP_NUM},'width',0.4)
ylabel('Time (s)');
xlabel('Subject Code');
hold off;
% axis padded
fig_i = get(groot,'CurrentFigure');
fig_i.Position = [200,200,1820,920];
saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',measure_name)]);
saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',measure_name)]);
%}
%## By trial plot
figure;
title([plot_name ' Across Subjects']);
hold on;
violinplot([b{:,65:end}],trial_names,'width',0.4)
ylabel('Time (s)');
xlabel('Trial Code');
hold off;
% axis padded
fig_i = get(groot,'CurrentFigure');
fig_i.Position = [200,200,1280,920];
saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',measure_name)]);
saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.jpg',measure_name)]);

