%   Project Title: Run BATCH processing for EEG filtering and artifact cleaning
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230223.0
%   Previous Version: n/a
%   Summary: 

%- run script
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/run_batch_preprocess.sh
%- run amica
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/run_singlenode_amica_run.sh
%- run amica_2nd
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/run_singlenode_amica_run_2nd.sh
%- mim dipfit
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/source_mim_mcc_dipfit_exe.sh

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
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% DEFINE SOURCE DIRECTORY & CD ======================================== %%
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '1_BATCH_PREP' filesep 'MIM_OA'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = true;
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
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('oa');
% subject_names = [SUBJ_PICS{:}];
%-
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
subject_names = [SUBJ_PICS{:}];
%% (PROCESSING PARAMS) ================================================= %%
%## hard define
%- dataset name
DATA_SET = 'MIM_dataset';
%- datetime override
% OA_PREP_FNAME = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FNAME = '07122023_OAN79_iccRX0p9_iccREMG0p3'; % JACOB,SAL(07/12/2023)
% OA_PREP_FNAME = '07142023_OAN79_iccRX0p55_iccREMG0p3_changparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(10/26/2023)
% OA_PREP_FNAME = '10302023_OAN82_iccRX0p60_iccREMG0p4_newparams'; % JACOB,SAL(10/30/2023)
% OA_PREP_FNAME = 'EMG_ANALYSIS'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p60_iccREMG0p3_newparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '11262023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(07/14/2023)
OA_PREP_FNAME = '01132024_antsnorm_iccREEG0p65_iccREMG0p4_skull0p0042';
%## soft define
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET]; % JACOB,SAL(02/23/2023)
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',OA_PREP_FNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (FIND FILES) ======================================================== %%
fnames          = cell(1,length([SUBJ_ITERS{:}]));
fpaths          = cell(1,length([SUBJ_ITERS{:}]));

for subj_i = 1:length(subject_names)
    try
        tmp = dir([save_dir filesep subject_names{subj_i} filesep 'clean' filesep '*.set']);
        fnames{subj_i} = tmp.name;
        fpaths{subj_i} = tmp.folder;

    catch e
        fprintf('Error. Subject %s has no .fdt file...\n',subject_names{subj_i});
    end
end
tmp = ~cellfun(@isempty,fnames);
fnames = fnames(tmp);
tmp = ~cellfun(@isempty,fpaths);
fpaths = fpaths(tmp);
%% CONCATENATE FLOAT DATA AND SAVE]
zeros
for subj_i = 1:length(subject_names)
    EEG = pop_loadset('filepath',fpaths{subj_i},'filename',fnames{subj_i});
    
end
