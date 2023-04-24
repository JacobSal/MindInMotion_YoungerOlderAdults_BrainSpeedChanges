%%% =================================================================== %%%
% IDENTIFY BRAIN-LIKE COMPONENTS
%%% =================================================================== %%%
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/1_ic_validation/MIM/run_ic_rejection_crit.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%- TIC
tic
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    DO_UNIX = false;
    PATH_EXT = 'M';
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep '1_ic_validation' filesep 'MIM'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP
if ~ispc
%     eeg_options;
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, 'option_parallel', 0,...
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
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data directory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions('option_parallel',0,'option_storedisk',1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
%- MIND IN MOTION (SUBSET (07/25/2022)
% SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
%     'H1022','H1024','H1026','H1027','H1033','H1034'};
% SUBJ_HMA = {'H2002', 'H2010', 'H2015', 'H2017', 'H2020', 'H2021', 'H2022', 'H2023',...
%     'H2025', 'H2026', 'H2034', 'H2059', 'H2062', 'H2082', 'H2095'};
% SUBJ_NMA = {'NH3008', 'NH3043', 'NH3055', 'NH3059', 'NH3069', ...
%     'NH3070', 'NH3074', 'NH3086', 'NH3090', 'NH3104', 'NH3105', 'NH3106', 'NH3112', 'NH3114'};
%- MIND IN MOTION (SUBSET (03/10/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};

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
% SUBJ_PICS = {SUBJ_2HMA,SUBJ_3HMA,SUBJ_3NHMA};
SUBJ_PICS = {SUBJ_2HMA,SUBJ_3NHMA};
GROUP_NAMES = {'H2000''s','H3000''s'};
% SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3HMA),1:length(SUBJ_3NHMA)}; % JACOB,SAL(02/23/2023)
SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3NHMA)}; 
%% ===================================================================== %%
%## PROCESSING PARAMS
%- statistics
STAT_ALPHA = 0.05;
%- datetime override
% dt = '16032023_OA_subset';
% dt = '23032023_OA_subset';
dt = '04092023_MIM_OA_subset_N85_speed_terrain';
%- hard define
% load_trials = {'0p25'};
% load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
load_trials = {'0p25','0p5','0p75','1p0'};
% load_trials = {'flat','low','med','high'}; 
%- soft define
subjinfDir = [SUBJINF_DIR filesep sprintf('%s',dt)];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'ic_rejection_crit'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(subj_list)]);
else
    POOL_SIZE = 1;
end
subj_list = [SUBJ_PICS{:}];
RUN_ALL_CRITERIA = true;
RUN_POWPOWCAT = true;
for subj_i = 53:length(subj_list)
% parfor (subj_i = 1:length(subj_list),POOL_SIZE)
    subj_str = subj_list{subj_i};
    %- Load EEG
    EEG = pop_loadset('filename',sprintf('%s_allcond_ICA_TMPEEG.set',subj_str),'filepath',[load_dir filesep subj_str filesep 'ICA']);
    out_fPath = [save_dir filesep subj_str];
    if ~exist(out_fPath,'dir')
        mkdir(out_fPath)
    end
    if RUN_ALL_CRITERIA
       crit = mim_reject_ics(EEG,out_fPath);
       close all;
    else
        tmp = load([out_fPath filesep sprintf('%s_ICRej.mat',subj_str)]);
        Output_ICRejection = tmp.Output_ICRejection;
        IC_all_brain = Output_ICRejection.IC_all_brain;
    end
    %% Last check Criteria 5: PowPow Cat Cross-Frequency Power-Power Coupling Analysis: A Useful Cross-Frequency Measure to Classify ICA-Decomposed EEG
    % It takes a long time to run. should pick only the ones that are
    % classified as 'brain'
    % %powpowcat parameters
    %{
    if RUN_POWPOWCAT
        upperFreqLimit = 100; %Frequency in Hz
        inputDataType = 2; %1, electrode data; 2, ICA time series
        methodType = 2;%1, Pearson's correlation; 2, Speaman's correlation
        numIterations = [];
        IC_powpow = find(IC_all_brain >= 8);
        fprintf('PowPowCAT parameters:\n upperFreqLimit= %i Hz\n inputDataType = ICs\n methodType= Spearman''s correlation (non-parametric)\n numIterations = %i\n',upperFreqLimit,numIterations);
        %run PowPowCAT
        % make a copy of EEG for powpowcat processing only
        if ~isempty(IC_powpow)
            EEG_powpow = EEG;
            EEG_powpow.icaact = EEG_powpow.icaact(IC_powpow,:);
            EEG_powpow = calc_PowPowCAT_CL(EEG_powpow, upperFreqLimit, inputDataType, methodType, numIterations);%CL version does not run stats
            EEG_powpow.setname = 'powpowcat';
            [ALLEEG, EEG_powpow, CURRENTSET] = eeg_store(ALLEEG, EEG_powpow , 0);
            eeglab redraw
            Plot35IcPushbutton_powpow(EEG_powpow,length(IC_powpow),IC_powpow)
            saveas(gcf,[out_fPath filesep spsubj_str,['powpowcat.fig']))
        end
    end
    %}
    %- SAVE IC_rejection 
%     EEG = pop_saveset(EEG,'filepath',EEG.filepath,'filename',sprintf('%s_IC_rej.set',EEG.subject));
end

