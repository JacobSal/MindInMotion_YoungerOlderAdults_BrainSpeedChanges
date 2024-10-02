%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh

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
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = fileparts(STUDY_DIR);
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_mri_study');
%% (TRANSFER IMU/LS CSV DATA TO _STUDIES) ============================== %%
% EEG_TRIAL_FNAMES = {'([mM]otor\w*)_(\d)(.*)','[rR]est(.*)',...
%     '([sS][pP]_0p25_\d)(.*)','([sS][pP]_0p5_\d)(.*)','([sS][pP]_0p75_\d)(.*)','([sS][pP]_1p0_\d)(.*)',...
%     '([tT][mM]_flat_\d)(.*)','([tT][mM]_low_\d)(.*)','([tT][mM]_med_\d)(.*)','([tT][mM]_high_\d)(.*)'};
EEG_TRIAL_FNAMES = {'([mM]otor\w*)_(\d)(.set)','[rR]est(.set)',...
    '([sS][pP]_0p25_\d)(.set)','([sS][pP]_0p5_\d)(.set)','([sS][pP]_0p75_\d)(.set)','([sS][pP]_1p0_\d)(.set)',...
    '([tT][mM]_flat_\d)(.set)','([tT][mM]_low_\d)(.set)','([tT][mM]_med_\d)(.set)','([tT][mM]_high_\d)(.set)'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
%- 
chk_dir = dir([R_MIND_IN_MOTION_DIR filesep '*']);
exc_fnames = {'.DS_Store','.','..','(do not use)_H1032_OLD_2022-06-14','old do not use)H1030',...
    'SubjectFolderTemplate - Copy','TempExportTrials','Thumbs.db'};
subj_fnames = {chk_dir.name};
chk = cellfun(@(x) contains(subj_fnames,x),exc_fnames,'UniformOutput',false);
chk = any(cat(1,chk{:}));
subj_fnames = subj_fnames(~chk);
subj_fpath = chk_dir.folder;
%## CHECK EEG TRIAL DATA
def_trial_struct = struct('subj_char',{{}},...
    'speed_0p25',{{}},...
    'speed_0p5',{{}},...
    'speed_0p75',{{}},...
    'speed_1p0',{{}},...
    'terrain_flat',{{}},...
    'terrain_low',{{}},...
    'terrain_med',{{}},...
    'terrain_high',{{}},...
    'rest_trial',{{}},...
    'motorimagery_trials',{{}});
trial_struct = def_trial_struct;
cnt = 1;
ext = ['EEG' filesep 'Trials'];
speed_order = {'speed_0p25','speed_0p5','speed_0p75','speed_1p0'};
terrain_order = {'terrain_flat','terrain_low','terrain_med','terrain_high'};
for subj_i = 1:length(subj_fnames)
    subj_dir = dir([subj_fpath filesep subj_fnames{subj_i} filesep ext]);
    trial_struct(cnt).subj_char = subj_fnames{subj_i};
    %- motorimagery
    tmp = regexp({subj_dir.name},EEG_TRIAL_FNAMES{1},'match');
    trial_struct(cnt).motorimagery_trials = [tmp{~cellfun(@isempty,tmp)}];
    % trial_struct(cnt).motor_imager_trials = cell2csv_util([tmp{~cellfun(@isempty,tmp)}]);
    %- rest
    tmp = regexp({subj_dir.name},EEG_TRIAL_FNAMES{2},'match');
    trial_struct(cnt).rest_trial = [tmp{~cellfun(@isempty,tmp)}];
    %- speed
    cnts = 1;
    for match_i = 3:6
        tmp = regexp({subj_dir.name},EEG_TRIAL_FNAMES{match_i},'match');
        trial_struct(cnt).(speed_order{cnts}) = [tmp{~cellfun(@isempty,tmp)}];
        cnts = cnts + 1;
    end
    %- terrain
    cnts = 1;
    for match_i = 7:10
        tmp = regexp({subj_dir.name},EEG_TRIAL_FNAMES{match_i},'match');
        trial_struct(cnt).(terrain_order{cnts}) = [tmp{~cellfun(@isempty,tmp)}];
        cnts = cnts + 1;
    end
    cnt = cnt + 1;
    trial_struct(cnt) = def_trial_struct;
end
trial_struct = struct2table(trial_struct);
writetable(trial_struct,[PATHS.src_dir filesep '_data' filesep 'MIM_Dataset',...
    filesep '_studies' filesep 'subject_mgmt' filesep 'trial_data_script_analysis.xlsx'])