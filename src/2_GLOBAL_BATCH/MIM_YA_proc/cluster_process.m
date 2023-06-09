%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/MIM_YA_proc/run_cluster_process.sh

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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'MIM_YA_proc'];
%- cd to source directory
cd(source_dir)
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
    'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047','H1047'}; % JACOB,SAL (04/18/2023)
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
%- (OY) Subject Picks 
SUBJ_PICS = {SUBJ_1YA}; 
GROUP_NAMES = {'H1000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_1YA)}; 
%- (OA) Subject Picks 
% SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H2000''s','H3000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
% TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- eeglab_cluster.m spectral params
FREQ_LIMITS = [1,100];
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
PAD_RATIO = 2;
%- datetime override
dt = '05302023_MIM_YAN33_subset_prep_verified_gait';
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
end
%% ===================================================================== %%
%- NOTES: spec_mode, 'psd' does not work.
[STUDY,ALLEEG] = eeglab_cluster(STUDY,ALLEEG,...
                    'SPEC_MODE',SPEC_MODE,...
                    'FREQ_LIMITS',FREQ_LIMITS,...
                    'CYCLE_LIMITS',CYCLE_LIMITS,...
                    'FREQ_FAC',FREQ_FAC,...
                    'PAD_RATIO',PAD_RATIO,...
                    'DO_ERSP_CALC',false,...
                    'DO_SPEC_CALC',false,...
                    'KMEANS_CLUSTER_N',10,...
                    'PAD_RATIO',4);
%## PCA reduction algorithm
%- RECALCULATE ICAACT MATRICES
%     ALLEEG = eeg_checkset(ALLEEG,'loaddata');
%     for subj_i = 1:length(ALLEEG)
%         if isempty(ALLEEG(subj_i).icaact)
%             fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
%             ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
%         end
%     end
% [STUDY,ALLEEG,~,~] = cluster_pca_reduce(STUDY,ALLEEG);
%## ICA reduction algorithm
[STUDY,~,~] = cluster_ica_reduce(STUDY);
%## ADD ANATOMICAL LABELS
[~, atlas_cell] = add_anatomical_labels(STUDY,ALLEEG);
STUDY.etc.add_anatomical_labels = atlas_cell;
% STUDY.etc.pipe_params = params;
% [ALLEEG,MAIN_STUDY] = parfunc_rmv_subjs(tmp,MAIN_STUDY,rmv_subj);
%- Save
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                study_fName_1,save_dir,...
                                'RESAVE_DATASETS','on');