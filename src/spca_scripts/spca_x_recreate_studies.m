%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/spca_scripts/run_recreate_spca_studies.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
% clearvars;
clear all;
clear classes
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SRC_DIR);
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
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
%## EPOCH PARAMS
DEF_EPOCH_PARAMS = struct('epoch_method','timewarp',...
    'percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'baseline_time_lims',[-0.5,4.5-2],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'cond_field','cond',...
    'appx_cond_len',3*60,...
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{'rest'}},...
    'do_recalc_epoch',true);
PREPROC_DIR_NAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
study_dir_name_from = '03232024_spca_analysis_OA';
study_dir_name_to = '03232024_spca_analysis_OA';
%## soft define
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
study_fname_rest_to = 'rest_epoch_study_yaoa';
study_fname_rest_from = 'rest_epoch_study';
study_fname_gait_to = 'epoch_study_yaoa';
study_fname_gait_from = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
load_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name_from)];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name_to)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
subjects_to = [SUBJ_PICS{:}];
%% ===================================================================== %%
%## LOAD STUDY
% if ~ispc
%     tmp = load('-mat',[load_dir filesep sprintf('%s_UNIX.study',study_fname_gait_from)]);
%     MAIN_STUDY = tmp.STUDY;
% else
%     tmp = load('-mat',[load_dir filesep sprintf('%s.study',study_fname_gait_from)]);
%     MAIN_STUDY = tmp.STUDY;
% end
%% INITIALIZE PARFOR LOOP VARS
ALLEEG = cell(length(subjects_to),1);
% for i = 1:length(ALLEEG)
parfor (i = 1:length(ALLEEG),SLURM_POOL_SIZE)
    %- define tmps
    % tmp_epoch_params = DEF_EPOCH_PARAMS;
    % EEG = pop_loadset('filepath',[load_dir filesep sprintf('%s',subjects_to{i}) filesep 'GAIT_EPOCHED' filesep [tmp_epoch_params.gait_trial_chars{:}]],...
    %     'filename',sprintf('%s_%s_EPOCH_TMPEEG.set',subjects_to{i},[tmp_epoch_params.gait_trial_chars{:}]));
    % fprintf('Adding subject: %s\n',EEG.subject);
    % %-
    % EEG = eeg_checkset(EEG,'loaddata');
    % if isempty(EEG.icaact)
    %     fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    %     EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    %     EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    % end
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    fpath = [load_dir filesep sprintf('%s',subjects_to{i}) filesep 'GAIT_EPOCHED' filesep [tmp_epoch_params.gait_trial_chars{:}]];
    fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subjects_to{i},[tmp_epoch_params.gait_trial_chars{:}]);
    EEG = eeglab_reload_eeg(fpath,fname);
    ALLEEG{i} = EEG;    
end
%%
ALLEEG = util_resolve_struct(ALLEEG,cellfun(@(x) x.subject,ALLEEG,'UniformOutput',false));
%##
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fname_gait_to,...
                                'filename',study_fname_gait_to,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
% STUDY.etc = MAIN_STUDY.etc;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
          
%% ===================================================================== %%
%## LOAD STUDY
% if ~ispc
%     tmp = load('-mat',[load_dir filesep sprintf('%s_UNIX.study',study_fname_rest_from)]);
%     MAIN_STUDY = tmp.STUDY;
% else
%     tmp = load('-mat',[load_dir filesep sprintf('%s.study',study_fname_rest_from)]);
%     MAIN_STUDY = tmp.STUDY;
% end
%% INITIALIZE PARFOR LOOP VARS
ALLEEG = cell(length(subjects_to),1);
parfor (i = 1:length(ALLEEG),SLURM_POOL_SIZE)
% for i = 1:length(ALLEEG)
    % %- define tmps
    % tmp_epoch_params = DEF_EPOCH_PARAMS;
    % %- grab subjects
    % EEG = pop_loadset('filepath',[load_dir filesep sprintf('%s',subjects_to{i}) filesep 'GAIT_EPOCHED' filesep [tmp_epoch_params.rest_trial_char{:}]],...
    %     'filename',sprintf('%s_%s_EPOCH_TMPEEG.set',subjects_to{i},[tmp_epoch_params.rest_trial_char{:}]));
    % fprintf('Adding subject: %s\n',EEG.subject);
    % EEG = eeg_checkset(EEG,'loaddata');
    % if isempty(EEG.icaact)
    %     fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    %     EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    %     EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    % end
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    fpath = [load_dir filesep sprintf('%s',subjects_to{i}) filesep 'GAIT_EPOCHED' filesep [tmp_epoch_params.rest_trial_char{:}]];
    fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subjects_to{i},[tmp_epoch_params.rest_trial_char{:}]);
    EEG = eeglab_reload_eeg(fpath,fname);
    ALLEEG{i} = EEG;
end
%%
ALLEEG = util_resolve_struct(ALLEEG,cellfun(@(x) x.subject,ALLEEG,'UniformOutput',false));
%##
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fname_rest_to,...
                                'filename',study_fname_rest_to,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
% STUDY.etc = MAIN_STUDY.etc;
[~,~] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
STUDY = struct.empty;
ALLEEG = struct.empty;