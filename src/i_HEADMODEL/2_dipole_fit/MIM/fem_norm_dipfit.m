%   Project Title: Transfer major project files from long term storage
%   drive to hypercomputer storage
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230223.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/run_fem_dipole_fit.sh

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
    %then the number of iterations in your for loop)
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
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
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%- MIND IN MOTION (SUBSET (03/10/2023)
% SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
%             'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
%             'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
% SUBJ_MISSING_TRIAL_DATA = {'H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
%     'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
% SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
%     'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1032','H1033','H1034','H1035',...
%     'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047','H1047'}; % JACOB,SAL (04/18/2023)
SUBJ_NO_MRI = {'H2010', 'H2036', 'H2041', 'H2072', 'H3018','H3120'};
SUBJ_2HMA = {'H2017', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2037', 'H2038',...
    'H2039', 'H2042', 'H2052', 'H2059', 'H2062', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3NHMA = {'H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107',...
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
SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3NHMA)};
%- Subject Picks
% SUBJ_PICS = {SUBJ_YNG};
% SUBJ_ITERS = {1:length(SUBJ_YNG)};
% GROUP_NAMES = {'H1000s'};
%% ================================================================= %%
%## PROCESSING PARAMS
%- OA
OA_PREP_FPATH = '07042023_OA_prep_verified'; % JACOB,SAL(04/10/2023)
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
%- YA
% YA_PREP_FPATH = '04182023_YA_N37_prep_verified'; % JACOB,SAL(04/10/2023)
% OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep YA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
%% global script chain (VERSION 1)
%- datetime override
% dt = '04172023_MIM_OA_subset_N85_speed_terrain_merge';
dt = 'test';
%## PATH & TEMP STUDY NAME
%- hard define
NORMALIZE_MRI = true;
%- soft define
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
working_dirs    = cell(1,length([SUBJ_ITERS{:}]));
fiducial_fPaths = cell(3,length([SUBJ_ITERS{:}]));
chanlocs_fPaths = cell(1,length([SUBJ_ITERS{:}]));
simnibs_fPaths  = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}])); 
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; 
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %## Generate Headmodels from MRI and Headscan
        working_dirs{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI'];
        %- Segmentation (SimNIBS)
        simnibs_fPaths{cnt} = [working_dirs{cnt} filesep sprintf('%s_masks_contr.nii.gz',SUBJ_PICS{group_i}{subj_i})];
        %- Fiducials (acpc_rs)
        fiducial_fPaths{1,cnt} = [working_dirs{cnt} filesep 'mri_acpc_rs.mat'];
        fiducial_fPaths{2,cnt} = [working_dirs{cnt} filesep 'mri_acpc.mat'];
        fiducial_fPaths{3,cnt} = [working_dirs{cnt} filesep 'ctf_fiducials.mat'];
        %- Chanlocs fPaths
        chanlocs_fPaths{cnt} = [working_dirs{cnt} filesep 'CustomElectrodeLocations.txt'];
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        fNames{cnt} = tmp.name;
        %- Prints
        if ~isempty(working_dirs{cnt})
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')));
%             fprintf('MRI working directory Exists: %i\n',exist(working_dirs{cnt},'dir'));
            fprintf('SimNIBS Segmentation Exists: %i\n',exist(simnibs_fPaths{cnt},'file'));
            fprintf('Digitized Electrodes Exists: %i\n',exist(chanlocs_fPaths{cnt},'file'));
            fprintf('MRI Fiducials Exist: %i\n',exist(fiducial_fPaths{1,cnt},'file') && exist(fiducial_fPaths{2,cnt},'file') && exist(fiducial_fPaths{3,cnt},'file'));
            fprintf('vol.mat check: %i\n',exist([working_dirs{cnt} filesep 'vol.mat'],'file'))
            fprintf('elec_aligned.mat check: %i\n',exist([working_dirs{cnt} filesep 'elec_aligned.mat'],'file'))
        end
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
working_dirs = working_dirs(~cellfun(@isempty,working_dirs));
fPaths = fPaths(~cellfun(@isempty,working_dirs));
fNames = fNames(~cellfun(@isempty,working_dirs));
chanlocs_fPaths = chanlocs_fPaths(~cellfun(@isempty,working_dirs));
fiducial_fPaths = fiducial_fPaths(:,~cellfun(@isempty,working_dirs));
subjectNames = subjectNames(~cellfun(@isempty,working_dirs));
%% SET POOLSIZE
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(working_dirs)]);
else
    POOL_SIZE = 1;
end
%% LOOP THROUGH PARTICIPANTS
LOOP_VAR = 1:length(working_dirs);
% parfor (subj_i = LOOP_VAR, POOL_SIZE)
for subj_i = LOOP_VAR
    %% LOAD SUBJECT DATA
    tic
    %- load MRI
    fprintf('Load MRI info...\n');
    tmp = load(fiducial_fPaths{1,subj_i});
    mri_acpc_rs = tmp.mri_acpc_rs;
    %- assign amica_folder
    fprintf('Load EEG...\n');
%     amica_folder = [fPaths filesep fNames];
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    %- load dipfit_fem.mat
    fprintf('Load dipfit_fem.mat...\n');
    tmp = load([working_dirs{subj_i} filesep 'dipfit' filesep 'dipfit_struct.mat']);
    dipfit_fem = tmp.SAVEVAR;
    %% 
    dippos = zeros(length(dipfit_fem.component),3);
    for i=1:length(dipfit_fem.component)
        %- 
        EEG.dipfit_fem.model(i).posxyz = dipfit_fem.dip(i).pos;
        EEG.dipfit_fem.model(i).momxyz = reshape(dipfit_fem.dip(i).mom, 3, length(dipfit_fem.dip(i).mom)/3)';
        if ~isempty(dipfit_fem.dip(i).rv)
            EEG.dipfit_fem.model(i).rv     = dipfit_fem.dip(i).rv;
        else
            EEG.dipfit_fem.model(i).rv     = NaN;
        end
        %- 
        EEG.dipfit_fem.model(i).diffmap = dipfit_fem.Vmodel(i) - dipfit_fem.Vdata(i);
        EEG.dipfit_fem.model(i).sourcepot = dipfit_fem.Vmodel(i);
        EEG.dipfit_fem.model(i).datapot   = dipfit_fem.Vdata(i);

        dippos(i,:) = dipfit_fem.dip(i).pos;
    end
    %% NORMALIZE MRI
    if NORMALIZE_MRI
        cfg = [];
        cfg.nonlinear = 'yes';
        cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
        mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
        %- PLOT
        ft_sourceplot(cfg,mri_norm);
        
    end
    %% Convert the dipole location to MNI space
    dipfit_fem_pos = reshape([dipfit_fem.dip(:).pos],3,[])';
    dipfit_fem_mnipos = ft_warp_apply(mri_norm.params,ft_warp_apply(mri_norm.initial,dipfit_fem_pos), 'individual2sn');
    dipfit_fem_mni_voxinds = round(ft_warp_apply(pinv(mri_norm.transform), dipfit_fem_mnipos ));
    for i=1:length(dipfit_fem.component)   
        dipfit_fem.dip(i).mnipos = dipfit_fem_mnipos(i,:);
        dipfit_fem.dip(i).mni_voxinds = dipfit_fem_mni_voxinds(i,:);
        EEG.dipfit_fem.model(i).mnipos = dipfit_fem_mnipos(i,:);
        EEG.dipfit_fem.model(i).mni_voxinds = dipfit_fem_mni_voxinds(i,:);
        EEG.dipfit_fem.model(i).pos_old = EEG.dipfit_fem.model(i).posxyz;
        EEG.dipfit_fem.model(i).posxyz = EEG.dipfit_fem.model(i).mnipos;
    end
    %% save
    EEG.dipfit = EEG.dipfit_fem;
    dipfit_fem_model = EEG.dipfit;
    save([fPaths{cnt} filesep 'dipfit_fem_norm.mat'],'dipfit_fem_model')
    EEG = pop_saveset(EEG);  
end
%## TIME
toc