%   Project Title: Transfer major project files from long term storage
%   drive to hypercomputer storage
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230223.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/run_fem_norm_dipfit.sh

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
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep 'i_HEADMODEL' filesep '2_dipole_fit' filesep 'MIM'];
%- cd to source directory
cd(source_dir)
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
    % see. eeg_optionsbackup.m for all eeglab options.
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
    %then the number of iterations in your for loop
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (DATASET INFORMATION) =============================================== %%
%## (MIND IN MOTION) DATASET SPECIFIC PARAMS (05/24/2023)
SUBJ_1YA = {'H1002','H1004','H1007','H1009',...
    'H1010','H1011','H1012','H1013','H1017',...
    'H1018','H1019','H1020','H1022','H1024',...
    'H1025','H1026','H1027','H1029','H1030','H1031',...
    'H1032','H1033','H1034','H1035','H1036',...
    'H1037','H1038','H1039','H1041','H1042',...
    'H1044','H1045','H1046','H1047','H1048'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2002','H2007','H2008',...
    'H2013','H2015','H2017','H2020','H2021',...
    'H2022','H2023','H2025','H2026','H2027',...
    'H2033','H2034','H2037','H2038','H2039',...
    'H2042','H2052','H2059','H2062','H2082',...
    'H2090','H2095','H2111','H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3053',...
    'H3063','H3072','H3077','H3103',...
    'H3107',...
    'NH3006','NH3007','NH3008','NH3010','NH3021',...
    'NH3026','NH3030','NH3036','NH3040',...
    'NH3041','NH3043','NH3054',...
    'NH3055','NH3058','NH3059','NH3066',...
    'NH3068','NH3069','NH3070','NH3074',...
    'NH3076','NH3086','NH3090','NH3102',...
    'NH3104','NH3105','NH3106','NH3108','NH3110',...
    'NH3112','NH3113','NH3114','NH3123','NH3128',...
    };
SUBJ_SLOW_WALKERS = {'H3042','H3046','H3047','H3073',...
    'H3092','NH3025','NH3051','NH3056','NH3071','NH3082'};
SUBJ_NO_MRI = {'H2010','H2012','H2018','H2036','H2041',...
    'H2072','H3018','H3120','NH3002','NH3009','NH3027','NH3129'};
SUBJ_MISSING_COND = {'H3024','NH3028'};
% SUBJ_UNKNOWN_ERR = {'NH3108','NH3030','NH3040','NH3025'};
% (08/21/2023) JS, 
% NH3108 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3030 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3040 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3025 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% (08/22/2023) JS, NH3108 bug seems to be related to entry errors in the
% Trial_Cropping_V2_test.xlsx sheet used to remove bad time ranges
% identified during collection. (fixed)
% NH3030 bug was due to how the CropTrialCheckFunc_checkLoadsol.m
% interpreted subject characters. It would consider 'NH3030_FU' as
% 'NH3030'. Changed from 'contains' to 'strcmp' func. (fixed)
% NH3040 bug was due to an entry error in Trial_Cropping_V2_test.xlsx (fixed)
SUBJ_DONT_INC = {'NH3004','NH3023'};
% (08/20/2023) JS, NH3004 has no headscan; NH3023 has no headscan clicks;
% fprintf('Total subjects processing: %i\n',sum([length(SUBJ_2MA),length(SUBJ_3MA)]));
% fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%- (OY) Subject Picks 
% SUBJ_PICS = {SUBJ_1YA}; 
% GROUP_NAMES = {'H1000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_1YA)}; 
%- (OA&YA) Subject Picks 
SUBJ_PICS = {SUBJ_1YA,SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_1YA),1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- (OA) Subject Picks 
% SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
% GROUP_NAMES = {'H2000''s','H3000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%- (0A) DEBUG SUBSET (06/17/2023)
% SUBJ_PICS = {SUBJ_DEBUG};
% GROUP_NAMES = {'debug'}; 
% SUBJ_ITERS = {1:length(SUBJ_DEBUG)};
%- test
% SUBJ_PICS = {{'H2062','NH3040'}};
% GROUP_NAMES = {'test'}; 
% SUBJ_ITERS = {1:length(SUBJ_PICS{:})};
fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PROCESSING PARAMS) ================================================= %%
%## Hard Define
%- subject choices (Combined OA & YA)
% YA_PREP_FPATH = '04182023_YA_N37_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FPATH = '07042023_OA_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
%- MRI normalization
NORMALIZE_MRI = true;
%## Soft Define
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)

%% Store fNames and fPaths
working_dirs    = cell(1,length([SUBJ_ITERS{:}]));
% fiducial_fPaths = cell(3,length([SUBJ_ITERS{:}]));
fiducial_fPaths = cell(1,length([SUBJ_ITERS{:}]));
chanlocs_fPaths = cell(1,length([SUBJ_ITERS{:}]));
simnibs_fPaths  = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}])); 
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
dipfit_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
ants_fpaths     = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; 
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %## Generate Headmodels from MRI and Headscan
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        working_dirs{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI'];
        %- MRI fPaths (ants)
        ants_fpaths = [working_dirs{cnt} filesep ];
        %- Segmentation (SimNIBS)
        simnibs_fPaths{cnt} = [working_dirs{cnt} filesep sprintf('%s_masks_contr.nii.gz',SUBJ_PICS{group_i}{subj_i})];
        %- Fiducials (acpc_rs)
        fiducial_fPaths{1,cnt} = [working_dirs{cnt} filesep 'mri_acpc_rs.mat'];
%         fiducial_fPaths{2,cnt} = [working_dirs{cnt} filesep 'mri_acpc.mat'];
%         fiducial_fPaths{3,cnt} = [working_dirs{cnt} filesep 'ctf_fiducials.mat'];
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        fNames{cnt} = tmp.name;
        %- Generated DIPFIT fPaths
        dipfit_fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'head_model' filesep 'dipfit_struct.mat'];
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')));
        fprintf('SimNIBS Segmentation Exists: %i\n',exist(simnibs_fPaths{cnt},'file'));
        fprintf('DIPFIT Exists: %i\n',exist(dipfit_fPaths{cnt},'file'));
%         fprintf('MRI Fiducials Exist: %i\n',exist(fiducial_fPaths{1,cnt},'file') && exist(fiducial_fPaths{2,cnt},'file') && exist(fiducial_fPaths{3,cnt},'file'));
        fprintf('MRI Fiducials Exist: %i\n',exist(fiducial_fPaths{1,cnt},'file'));
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
inds1 = logical(cellfun(@(x) exist(x,'file'),dipfit_fPaths));
inds2 = logical(cellfun(@(x) exist(x,'file'),fiducial_fPaths));
inds = inds1 & inds2;
working_dirs = working_dirs(inds);
dipfit_fPaths = dipfit_fPaths(inds);
fPaths = fPaths(inds);
fNames = fNames(inds);
fiducial_fPaths = fiducial_fPaths(:,inds);
subjectNames = subjectNames(inds);
%% SET POOLSIZE
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(working_dirs)]);
else
    POOL_SIZE = 1;
end
%% LOOP THROUGH PARTICIPANTS
LOOP_VAR = 1:length(working_dirs);
parfor (subj_i = LOOP_VAR, POOL_SIZE) % (05/24/2023) JS, parfor might not
% be possible for this loop. a problem with ft_sourceplot.
% for subj_i = LOOP_VAR
    [EEG,dipfit_fem_norm] = mim_norm_dipfit(fPaths{subj_i},fNames{subj_i},...
        fiducial_fPaths{1,subj_i},dipfit_fPaths{subj_i},...
        'MRI_NORM_METHOD','ants',...
        'ANTS_FPATH',);
    par_save(dipfit_fem_norm,fPaths{subj_i},'dipfit_fem_norm.mat')
    EEG = pop_saveset(EEG,'filepath',EEG.filepath,'filename',EEG.filename); 
end
%## TIME
toc