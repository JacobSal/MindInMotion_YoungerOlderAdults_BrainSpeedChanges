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

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM/run_a_ic_rejection_crit.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM'];
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
%% (PARAMETERS) ======================================================== %%
%## hard define
%- statistics
STAT_ALPHA = 0.05;
%- datetime override
% dt = '16032023_OA_subset';
% dt = '23032023_OA_subset';
% dt = '04092023_MIM_OA_subset_N85_speed_terrain';
% dt = '04172023_MIM_OA_subset_N85_speed_terrain_merge';
dt = '04172023_MIM_YA_fullset_N_speed_terrain_merge';
% load_trials = {'0p25'};
load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% load_trials = {'0p25','0p5','0p75','1p0'};
% load_trials = {'flat','low','med','high'};
%## soft define
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
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

