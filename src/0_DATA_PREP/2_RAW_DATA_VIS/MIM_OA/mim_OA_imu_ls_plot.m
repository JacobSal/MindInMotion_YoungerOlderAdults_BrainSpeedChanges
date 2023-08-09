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
SUBJ_2MA = {'H2017', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2037', 'H2038',...
    'H2039', 'H2042', 'H2052', 'H2059', 'H2062', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
SUBJ_DEBUG = {'H2117','NH3082','H3063','NH3006','NH3025','NH3114','H2007',...
    'H3034','NH3055','H3073','NH3104','NH3051','NH3123','H3092','NH3082',...
    'NH3056','NH3036','H3046','H3053','NH3007','H3077','H3047','NH3071'};
%- (YA,OA) Subject Picks 
% SUBJ_PICS = {SUBJ_1YA,SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_1YA),1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- OA Subject Picks
SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- test
% SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H2000''s','H3000''s'}; 
% SUBJ_ITERS = {[1,2],[1,2]};
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
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];%'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
N_TRIALS = 16;
%- Loop through directory
table_imu_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,12); % 12 measures
table_trial_vec = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_vec = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names = cell(length([SUBJ_PICS{:}]),1);

table_ls_meas = zeros(length([SUBJ_PICS{:}])*N_TRIALS,64); % 64 measures
table_trial_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_subj_ls = cell(length([SUBJ_PICS{:}])*N_TRIALS,1);
table_header_names_ls = cell(length([SUBJ_PICS{:}]),1);
subj_stack = [];
cnt_ga = 1;
cnt_ls = 1;
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
                trial_i = 'med';
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
                trial_i = 'med';
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
            table_trial_ls{cnt_ga} = tmp{1};
            table_subj_ls{cnt_ga} = SUBJ_PICS{group_i}{subj_i};
            cnt_ga = cnt_ga + 1;
        end
   
    end
end
%% IMU TABLE
table_imu_meas = table_imu_meas(~any(table_imu_meas == 0,2),:);
table_subj_vec = table_subj_vec(~any(table_imu_meas == 0,2));
table_trial_vec = table_trial_vec(~any(table_imu_meas == 0,2));
table_imu_out = array2table(table_imu_meas,'VariableNames',table_header_names{1});
table_imu_out.SubjectName = categorical(table_subj_vec);
table_imu_out.TrialName = categorical(table_trial_vec);
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
table_ls_out = array2table(table_ls_meas,'VariableNames',table_header_names_ls{1});
table_ls_out.SubjectName = categorical(table_subj_ls);
table_ls_out.TrialName = categorical(table_trial_ls);
%- rearrange headers
% table_ls_out = [table_ls_out(:,end-1), table_ls_out(:,end), table_ls_out(:,1:end-2)];
%% VIOLIN PLOT IMU
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,920,480];
% meas_names = {'APexc_mean','APexc_COV','MLexc_mean','MLexc_COV','UDexc_mean',...
%     'UDexc_COV','PeakAntVel_mean','PeakAntVel_COV','PeakLatVel_mean',...
%     'PeakLatVel_COV','PeakUpDownVel_mean','PeakUpDownVel_COV'};
% meas_units = {'m','m','m','m','m',...
%     'm','m/s','m/s','m/s',...
%     'm/s','m/s','m/s'};
meas_names = {'APexc_COV','MLexc_COV'};
meas_units = {'%%','%%'};
YLIMS = {[0,50],[0,50]};
% trial_names = {'SP_0p25_1','SP_0p25_2','SP_0p5_1','SP_0p5_2','SP_0p75_1',...
%                 'SP_0p75_2','SP_1p0_1','SP_1p0_2','TM_flat_1','TM_flat_2',...
%                 'TM_high_1','TM_high_2','TM_low_1','TM_low_2','TM_med_1',...
%                 'TM_med_2'};
trial_names_speed = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
                '0.75 m/s','1.00 m/s','1.00 m/s'};
trial_names_terrain = {'flat terr.','flat terr.',...
    'low terr.','low terr.','med terr.',...
    'med terr.','high terr.','high terr.'};
% table_in = table_imu_out;
for meas_i = 1:length(meas_names)
    table_in = table_imu_out;
    meas_names{meas_i} = meas_names{meas_i};
    dat_mean = mean(table_in.(meas_names{meas_i}));
    dat_std = std(table_in.(meas_names{meas_i}));
    in_crop = table_in.(meas_names{meas_i})>dat_mean-3*dat_std & table_in.(meas_names{meas_i})<dat_mean+3*dat_std;
    table_in = table_in(in_crop,:);
    a = unstack(table_in,meas_names{meas_i},'SubjectName');
    b = unstack(table_in,meas_names{meas_i},'TrialName');
    tmp = b(:,12:end);
%     b_speed = tmp(:,2:9);
    b_speed = tmp(:,2:5);
%     b_terrain = [tmp(:,10:11),tmp(:,14:15),tmp(:,16:17),tmp(:,12:13)];
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

    %## By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_terrain,trial_names_terrain','width',0.4,...
        'GroupOrder',b_terrain.Properties.VariableNames,'ViolinColor',linspecer(size(b_terrain,2)),...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1);
    ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
    xlabel('Trial Code');
    xtickangle(45)
    ylim(YLIMS{meas_i});
    hold off;
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = FIG_POSITION;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
    saveas(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    %## By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_speed,trial_names_speed,'width',0.4,...
        'GroupOrder',b_speed.Properties.VariableNames,'ViolinColor',linspecer(size(b_speed,2)),...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1);
    %- consider a colormap that matches for each condition across
    %trials....
    ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
    xlabel('Trial Code');
    xtickangle(45);
    ylim(YLIMS{meas_i});
    hold off;
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = FIG_POSITION;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
    saveas(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})]);
end
%% VIOLIN PLOT LOADSOL
% FIG_POSITION = [100,100,1480,520];
FIG_POSITION = [100,100,920,480];
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
meas_names = {'StanceDur_L','SwingDur_R','GaitCycleDur_cov','GaitCycleDur'};
meas_units = {'s','s','%%','s'};
YLIMS = {[0,1.5],[0,1.5],[0,50],[0,1.5]};
% trial_names = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
%                 '0.75 m/s','1.00 m/s','1.00 m/s','flat terr.','flat terr.',...
%                 'high terr.','high terr.','low terr.','low terr.','med terr.',...
%                 'med terr.'};
trial_names_speed = {'0.25 m/s','0.25 m/s','0.50 m/s','0.50 m/s','0.75 m/s',...
                '0.75 m/s','1.00 m/s','1.00 m/s'};
trial_names_terrain = {'flat terr.','flat terr.',...
    'low terr.','low terr.','med terr.',...
    'med terr.','high terr.','high terr.'};
% table_in = table_ls_out;
for meas_i = 1:length(meas_names)
    table_in = table_ls_out;
    meas_names{meas_i} = meas_names{meas_i};
    dat_mean = mean(table_in.(meas_names{meas_i}));
    dat_std = std(table_in.(meas_names{meas_i}));
    in_crop = table_in.(meas_names{meas_i})>dat_mean-3*dat_std & table_in.(meas_names{meas_i})<dat_mean+3*dat_std;
%     in_down = table_in.(meas_names{meas_i})<dat_mean+3*dat_std;
%     in_crop = ~(in_up+in_down);
    table_in = table_in(in_crop,:);
    meas_names{meas_i} = meas_names{meas_i};
    a = unstack(table_in,meas_names{meas_i},'SubjectName');
    b = unstack(table_in,meas_names{meas_i},'TrialName');
    tmp = b(:,64:end);
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
    %## By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_terrain,trial_names_speed,'width',0.4,'GroupOrder',b_terrain.Properties.VariableNames,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1);
    ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
    xlabel('Trial Code');
    ylim(YLIMS{meas_i})
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
    hold off;
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = FIG_POSITION;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
    saveas(fig_i,[save_dir filesep sprintf('Across_terrain_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    %## By trial plot
    figure;
%     title(sprintf('%s Across Trials',meas_names{meas_i}));
    hold on;
    violinplot(b_speed,trial_names_terrain,'width',0.4,'GroupOrder',b_speed.Properties.VariableNames,...
        'ShowWhiskers',false,'ShowNotches',false,'ShowBox',false,...
        'ShowMedian',true,'Bandwidth',range(b_terrain{:,2})*0.1);
    ylabel(sprintf('IMU %s (%s)',meas_names{meas_i},meas_units{meas_i}));
    xlabel('Trial Code');
    xtickangle(45)
    ylim(YLIMS{meas_i});
    hold off;
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = FIG_POSITION;
%     saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.fig',meas_names{meas_i})]);
    saveas(fig_i,[save_dir filesep sprintf('Across_speed_Trials_Fig_%s.jpg',meas_names{meas_i})]);
    %## STATISTICS
    
end

%% STATISTIFCS
modelspec = 'linear';
mdl = fitglm(table_ls_out,modelspec,'ResponseVar',meas_names{meas_i},'PredictorVar',{'flat','low','med','high','0p25','0p5','0p75','1p0'});
% mdl = fitglm(table_imu_out,modelspec,'ResponseVar',meas_names{meas_i},'PredictorVar',{'flat','low','med','high','0p25','0p5','0p75','1p0'});
