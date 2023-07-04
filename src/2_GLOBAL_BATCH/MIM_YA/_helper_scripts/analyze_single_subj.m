%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_YA_proc/run_conn_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SINGLE SUBJECT PROCESSING 
%## NECESSARY HARD PARAMS
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
MNI_CHAN_1005 = [PATH_EEGLAB_BEM filesep 'elec' filesep 'standard_1005.elc'];
COORD_TRANSFORM_MNI = [0 0 0 0 0 -1.5708 1 1 1];
DIP_NUM = 1;
DIP_PLOT = 'off';
DO_DIPOLE_REJECTION = true;
%-
TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
studyDir = save_dir;
%## SOFT PARAMS
subj_i = 5;
%%
fPath = fPaths{subj_i};
fName = fNames{subj_i};
subjStr = subjectNames{subj_i};
if ~ispc
    fPath = convertPath2UNIX(fPath); 
else
    fPath = convertPath2Drive(fPath);
end
fprintf(1,'Loading Subject %s\n',subjStr)
[~,EEG,~] = eeglab_loadICA(fName,fPath);
%- override fPath & fName
EEG.subject = subjStr;
fName = sprintf('%s_allcond_ICA_TMPEEG.set',EEG.subject);
fPath = [studyDir filesep EEG.subject filesep 'ICA'];
if ~strncmp(computer,'PC',2)
    fPath = convertPath2UNIX(fPath);
else
    fPath = convertPath2Drive(fPath);
end
if ~exist(fPath,'dir')
    mkdir(fPath)
end
% %- ADD DIPFIT STRUCT TO EEG STRUCT
% fprintf('%s) dipfit available: %i\n',EEG.subject,isfield(EEG,'dipfit'));
% try
%     EEG.dipfit.coord_transform;
%     EEG.dipfit.mrifile;
%     EEG.dipfit.hdmfile;
%     EEG.dipfit.coordformat;
% catch
%     fprintf('MNI pop_dipfit_settings...\n');
%      EEG = pop_dipfit_settings(EEG,'coordformat','MNI','coord_transform',COORD_TRANSFORM_MNI,...
%             'hdmfile',MNI_VOL,'mrifile',MNI_MRI,'chanfile',MNI_CHAN_1005);
% end
% try 
%     EEG.dipfit.model;
% catch
%     %- DIPFIT (see. ft_dipolefitting())
%     fprintf('pop_multifit...\n');
%     EEG = pop_multifit(EEG,[],'dipoles',DIP_NUM,'dipplot',DIP_PLOT);
%     tmp = EEG.dipfit;
%     par_save(tmp,fPath,sprintf('%s_dipfit_mni.mat',EEG.subject));
% end
% %- Check for iclabel classifications
% if ~isfield(EEG.etc,'ic_classification')
%     EEG = iclabel(EEG, 'lite');
% end
%- adapt each EEG structure for study creation
EEG.filepath = fPath;
EEG.filename = fName;
EEG.group = char(groups{subj_i});
EEG.condition = char(conditions{subj_i});
EEG.session = char(sessions{subj_i});
% fprintf(1,'Saving Subject %s\n',EEG.subject);
% [EEG] = pop_saveset(EEG,'savemode','twofiles',...
%     'filename',fName,...
%     'filepath',fPath);
%% EPOCHING
%- Recalculate ICA Matrices && Book Keeping
EEG = eeg_checkset(EEG,'loaddata');
if isempty(EEG.icaact)
    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
end
%## PARSE TRIALS
for trial_i = 1:length(TRIAL_TYPES)
    events = mim_trial_event_override(EEG.subject,TRIAL_TYPES{trial_i},TRIAL_OVERRIDE_FPATH);
    if ~isempty(events)
        f = fields(EEG.event);
        for e = 1:length(events)
            for f_i = 1:length(f)
                if ~isfield(events{e},f{f_i})
                    events{e}.(f{f_i}) = [];
                end
            end
        end
        events = [events{:}];
        fprintf('appending events...\n');
        EEG.event = [EEG.event, events];
    end
end
%- let EEGLAB rearrange the event order
EEG = eeg_checkset(EEG,'eventconsistency');
%- parse trials
ALLEEG = mim_parse_trials(EEG,TRIAL_TYPES,EPOCH_TIME_LIMITS,...
        'SLIDING_WINDOW_ONLY',true);
%## Save EEG
for trial_i = 1:length(ALLEEG)
    if ~exist([save_dir filesep ALLEEG(trial_i).subject filesep 'EPOCHED'],'dir')
        mkdir([save_dir filesep ALLEEG(trial_i).subject filesep 'EPOCHED']);
    end
    fprintf(1,'Saving Subject %s\n',ALLEEG(trial_i).subject); 
    [ALLEEG(trial_i)] = pop_saveset(ALLEEG(trial_i),'savemode','twofiles',...
        'filename',ALLEEG(trial_i).filename,...
        'filepath',[save_dir filesep ALLEEG(trial_i).subject filesep 'EPOCHED']);
end