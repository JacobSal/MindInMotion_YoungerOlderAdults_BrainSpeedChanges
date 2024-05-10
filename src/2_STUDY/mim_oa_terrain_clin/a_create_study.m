%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_terrain_clin/run_a_create_study.sh

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
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
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
%- study group and saving
SESSION_NUMBER = '1';
SAVE_ALLEEG = true;
SAVE_EEG = true; %true;
OVERRIDE_DIPFIT = true;
%- epoching params
DO_SLIDING_WINDOW = false;
%* sliding window
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
%* gait
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-0.5,4.5]; %[-1,3]; %[-0.5,5]; % [-1,3] captures gait events well , [-0.5,5] captures gait events poorly
% (10/13/2023) changing from [-1,4.25] to [-0.5,4.5] to match chang's
% (10/25/2023) changing from [-0.5,4.5] to [-1,4.25] as it seems to help
% with frequency decomposition artifact during ERSP creation
% paper
% (01/23/2024) changing from [-1,4.25] to [-0.5,4.5] to match chang
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
if DO_SLIDING_WINDOW
    SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
    TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
else
    SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
    TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
end
%- datetime override
% dt = '03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name_from = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name_from = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name_from = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name_to = '04232024_MIM_OAN57_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%## soft define
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
study_fname_to = 'epoch_study';
study_fname_from = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
load_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name_from)];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name_to)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('oa_spca');
subjects_to = [SUBJ_PICS{:}];
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fname_from '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fname_from '_UNIX.study'],'filepath',load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fname_from '.study'],'filepath',load_dir);
    end
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[load_dir filesep sprintf('%s_UNIX.study',study_fname_from)]);
    MAIN_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[load_dir filesep sprintf('%s.study',study_fname_from)]);
    MAIN_STUDY = tmp.STUDY;
end
%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_STUDY.datasetinfo.filepath};
fNames = {MAIN_STUDY.datasetinfo.filename};
subjects_from = {MAIN_STUDY.datasetinfo.subject};
alleeg_fpaths = cell(length(MAIN_STUDY.datasetinfo),1);
inds = zeros(length(subjects_from),1);
for i = 1:length(subjects_to)
    if any(strcmp(subjects_to{i},subjects_from))
        inds(strcmp(subjects_to{i},subjects_from)) = 1;
    end
end
inds = find(inds);
% inds = cellfun(@(x) find(strcmp(x,subjects_from)),subjects_to);
ALLEEG = MAIN_ALLEEG(inds);
% for i = 1:length(ALLEEG)
parfor (i = 1:length(ALLEEG),SLURM_POOL_SIZE)
    fprintf('Adding subject: %s\n',ALLEEG(i).subject);
    tmp = strsplit(ALLEEG(i).filepath,filesep);
    ind = find(strcmp(tmp,study_dir_name_from));
    tmp = strjoin(tmp(ind+1:end),filesep);
    %-
    ALLEEG(i) = eeg_checkset(ALLEEG(i),'loaddata');
    if isempty(ALLEEG(i).icaact)
        fprintf('%s) Recalculating ICA activations\n',ALLEEG(i).subject);
        ALLEEG(i).icaact = (ALLEEG(i).icaweights*ALLEEG(i).icasphere)*ALLEEG(i).data(ALLEEG(i).icachansind,:);
        ALLEEG(i).icaact = reshape( ALLEEG(i).icaact, size(ALLEEG(i).icaact,1), ALLEEG(i).pnts, ALLEEG(i).trials);
    end
    ALLEEG(i).filepath = [save_dir filesep tmp];
    mkdir(ALLEEG(i).filepath);
    pop_saveset(ALLEEG(i),'filepath',ALLEEG(i).filepath,...
        'filename',ALLEEG(i).filename,...
        'savemode','twofiles');
    
end
%%
%##
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fname_to,...
                                'filename',study_fname_to,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
STUDY.etc.a_epoch_process = MAIN_STUDY.etc.a_epoch_process;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
                                        
%% Version History
%{
v2.0; (04/28/2023) JS: Splitting up the epoching, plotting, and
connectivity calculation process to avoid bugs with STUDY files.
v1.0; (11/11/2022), JS: really need to consider updating bootstrap
    algorithm with parallel computing. Taking ~ 1 day per
    condition for all subjects and the bottle neck is entirely the
    bootstrap.

    Note: validateattributes and assert functions may be helpful
    in more clearly defining function inputs.
        e.g.  DO_PHASE_RND = true;
          errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
          validationFcn = @(x) assert(islogical(x),errorMsg);
v1.0; (12/5/2022) Need to adapt this to include all conditions
    within each SUBJ structure so connectivity can be calculated
    for the ALLEEG structure rather than the EEG structure.
    *** maybe try to ditch the SUBJ strucutre entirely for this
    round?
v1.0.01132023.0 : Initializing versioning for future iterations.
%}

