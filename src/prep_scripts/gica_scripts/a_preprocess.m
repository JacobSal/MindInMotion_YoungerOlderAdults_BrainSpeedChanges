%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 

%- run script
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/run_prep_a_preprocess.sh

%- run amica
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/b_run_singlenode_amica.sh

%- mim dipfit
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/c_run_mim_mcc_dipfit.sh

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
%## Determine Working Directories
ADD_ALL_SUBMODS = true;
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    SRC_DIR = fileparts(SCRIPT_DIR);
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
% SUBJ_PICS = {{'H3046','H3047','H3073','H3077','H3092', ...
%     'NH3023','NH3025','NH3027','NH3028', ...
%     'NH3051','NH3056','NH3071','NH3082','NH3123'}};
% SUBJ_PICS = {{'NH3028'}};
%% (PROCESSING PARAMS) ================================================= %%
%## hard define
%- dataset name
DATA_SET = 'MIM_dataset';
%- datetime override
% PREPROC_DNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
PREPROC_DNAME = '02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej';
%## soft define
%- path for local data
raw_dataset_dir = [PATHS.data_dir filesep DATA_SET]; % JACOB,SAL(02/23/2023)
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
save_dir = [studies_dir filesep sprintf('%s',PREPROC_DNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
group_cats = {'H1000','H2000','H3000'};
ICA_REGEXP = '%s_cleanEEG_EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set';
%-
subj_chars = [SUBJ_PICS{:}];
eeg_fpaths = cell(1,length(subj_chars));
%- loop
for subj_i = 1:length(subj_chars)
    %## PARAMS
    fprintf('==== (%i) Loading Subject %s ====\n',subj_i,subj_chars{subj_i})
    %- eeg fpath
    eeg_fpaths{subj_i} = [raw_dataset_dir filesep subj_chars{subj_i} filesep 'EEG' filesep 'Trials']; %[ica_data_dir filesep subj_chars{subj_i} filesep 'clean'];
    tmp = dir([eeg_fpaths{subj_i} filesep '*.set']);
    if ~isempty(tmp)
        chks = cell(length(tmp),1);
        for i = 1:length(tmp)
            tt = strsplit(tmp(i).name,'.');
            chks{i} = tt{1};
        end
    end
    fprintf('conditions available: %s\n\n',strjoin(chks,','));
end
%% (PARFOR) PREPROCESS EEG
amica_cmd = cell(length(eeg_fpaths),1);
params = cell(length(eeg_fpaths),1);
parfor subj_i = 1:length(eeg_fpaths)
% for subj_i = 1:length(subj_chars) %find(strcmp(subj_chars,'NH3113'))
    fprintf('Running subject %s...\n',subj_chars{subj_i})
    %## PREP for MAIN_FUNC
    if ~exist([save_dir filesep subj_chars{subj_i}],'dir')
        mkdir([save_dir filesep subj_chars{subj_i}]);
    end
    %## RUN MAIN_FUNC
    try
        [EEG,amica_cmd{subj_i},params{subj_i}] = mim_preproc_func(subj_chars{subj_i},...
            eeg_fpaths{subj_i}, ...
            [save_dir filesep subj_chars{subj_i}], ...
            studies_dir);
        % fprintf('%s\n',amica_cmd{subj_i}{2})
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
    end
end

%% (AMICA PARAM FILE RECOVERY) 
% preprocess_pipeline = 'EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10';
% for subj_i = 1:length(eeg_fpaths)
%     fpath = [save_dir filesep subj_chars{subj_i} filesep 'clean'];
%     float_fname = sprintf('%s_cleanEEG_%s.fdt',subj_chars{subj_i},preprocess_pipeline);    
%     set_fname = sprintf('%s_cleanEEG_%s.set',subj_chars{subj_i},preprocess_pipeline);
%     %-
%     EEG = pop_loadset('filepath',fpath,'filename',set_fname);
%     %-
%     [EEG,amica_cmd] = mim_prep_hpg_amica(EEG, ...
%         [fpath filesep float_fname], ...
%         fpath, ...
%         'jsalminen@ufl.edu', ...
%         1);
% end
