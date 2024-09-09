function [ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,varargin)
%MIM_CREATE_ALLEEG Summary of this function goes here
%   Detailed explanation goes here
%     MATLAB CODE DOCUMENTATION
% 
%     Function Name: mim_create_alleeg
% 
%     Description:
%     This MATLAB function, 'mim_create_alleeg', is used for creating and populating an EEGLAB STUDY structure from multiple EEG datasets. It loads EEG data files, performs Independent Component Analysis (ICA), and optionally fits dipole models for source localization. The resulting data is organized into an EEGLAB STUDY structure, which can be used for various EEG analysis tasks.
% 
%     Usage:
%     [ALLEEG] = mim_create_alleeg(fNames, fPaths, subjectNames, save_dir, varargin)
% 
%     Input Parameters:
%     - fNames: A cell array of EEG data file names (e.g., {'subject1.set', 'subject2.set'}).
%     - fPaths: A cell array of file paths to the EEG data files corresponding to 'fNames'.
%     - subjectNames: A cell array of subject names associated with the EEG data files.
%     - save_dir: The directory where the processed data and STUDY structure will be saved.
%     - varargin (Optional Parameters):
%       - conditions: A cell array specifying condition labels for each EEG data file. Default is 'tmp_cnd'.
%       - groups: A cell array specifying group labels for each EEG data file. Default is 'tmp_grp'.
%       - sessions: A cell array specifying session labels for each EEG data file. Default is 'tmp_sess'.
%       - chanlocs_fPaths: A cell array of full file paths to custom channel location files (optional).
% 
%     Output:
%     - ALLEEG: A cell array of EEG structures, one for each subject, with data and information populated.
% 
%     Important Notes:
%     - This code is designed to work with EEG data and requires the EEGLAB toolbox.
%     - The 'dipfit' plugin within EEGLAB is used for dipole fitting and source localization.
% 
%     Code Designer: Jacob Salminen
%     Date: 11/25/2022
% 
%     Functions:
%     1. fem_eeglab_dipfit: Configures dipole fitting parameters for Finite Element Model (FEM) dipfit analysis.
%     2. custom_update_chanlocs: Updates EEG channel locations using custom electrode locations.
%     3. bem_eeglab_dipfit: Fits dipole models to EEG data using Boundary Element Model (BEM) dipfit analysis.
% 
%     See individual function descriptions for more details on their functionality.

% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## DEFINE DEFAULTS
%- developer params
DO_BEM_DIPFIT = false;
DO_FEM_DIPFIT = true;
FORCE_RELOAD = false;
ICA_FNAME_REGEXP = '%s_allcond_ICA_TMPEEG.set';
%- find eeglab on path
if ~ispc
    tmp = strsplit(path,':');
else
    tmp = strsplit(path,';');
end
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
MNI_CHAN_1005 = [PATH_EEGLAB_BEM filesep 'elec' filesep 'standard_1005.elc'];
COORD_TRANSFORM_MNI = [0 0 0 0 0 -1.5708 1 1 1];
%- dipfit params
DIP_NUM = 1;
DIP_PLOT = 'off';
%- conditions, groups, sessions defs
%* CONDITIONS (e.g., ["rest","rest","rest","rest","rest"])
tmp = cell(1,length(fNames));
for i = 1:length(fNames)
    tmp{i} = 'tmp_cnd';
end
CONDITIONS = tmp;
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}. Condition label for each fPath & fName provided';
cnd_validFcn = @(x) assert((ischar(x{1}) && iscell(x) && length(x) == length(fNames) || isempty(x)), errorMsg);
%* GROUPS (e.g., ["1","1","1","1","1"])
tmp = cell(1,length(fNames));
for i = 1:length(fNames)
    tmp{i} = 'tmp_grp';
end
GROUPS = tmp;
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}. Group label for each fPath & fName provided';
grp_validFcn = @(x) assert((ischar(x{1}) && iscell(x) && length(x) == length(fNames) || isempty(x)), errorMsg);
%* SESSIONS (e.g., [1,1,1,1,1])
tmp = cell(1,length(fNames));
for i = 1:length(fNames)
    tmp{i} = 'tmp_sess';
end
SESSIONS = tmp;
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}. Session label for each fPath & fName provided';
sess_validFcn = @(x) assert((ischar(x{1}) && iscell(x) && length(x) == length(fNames) || isempty(x)), errorMsg);
%- save eeg 
DO_SAVE_ICA = false;
%- CHANLOCS_FPATHS
CHANLOCS_FPATHS = {};
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}. Full file paths for chanlocs (e.g., ''/path/to/subject/dir/chanlocs.mat'')';
cfp_validFcn = @(x) assert((isempty(x) || (iscell(x) && length(x) == length(fNames))), errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'fNames',@iscell);
addRequired(p,'fPaths',@iscell);
addRequired(p,'subjectNames',@iscell);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
addOptional(p,'conditions',CONDITIONS,cnd_validFcn);
addOptional(p,'groups',GROUPS,grp_validFcn);
addOptional(p,'sessions',SESSIONS,sess_validFcn);
%## PARAMETER
addParameter(p,'DO_SAVE_ICA',DO_SAVE_ICA,@islogical);
addParameter(p,'CHANLOCS_FPATHS',CHANLOCS_FPATHS,cfp_validFcn);
addParameter(p,'FORCE_RELOAD',FORCE_RELOAD,@islogical);
parse(p,fNames,fPaths,subjectNames,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
conditions          = p.Results.conditions;
groups              = p.Results.groups;
sessions            = p.Results.sessions;
DO_SAVE_ICA         = p.Results.DO_SAVE_ICA;
CHANLOCS_FPATHS     = p.Results.CHANLOCS_FPATHS;
FORCE_RELOAD        = p.Results.FORCE_RELOAD;
%% ===================================================================== %%
%## CREATE STUDY
fprintf(1,'==== Creating Study ====\n')
%* empty ALLEEG structure for repopulating
ALLEEG = cell(1,length(fNames)); 
% for subj_i=1:length(fNames)
parfor subj_i=1:length(fNames)
    fName = sprintf(ICA_FNAME_REGEXP,subjectNames{subj_i});
    fPath = [save_dir filesep subjectNames{subj_i} filesep 'ICA'];
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    if ~exist([fPath filesep fName],'file') || FORCE_RELOAD
        fprintf(1,'Loading Subject %s\n',subjectNames{subj_i})
        [~,EEG,~] = eeglab_loadICA(fNames{subj_i},fPaths{subj_i});
        %- override fPath & fName
        EEG.subject = subjectNames{subj_i};
        %## ADD DIPFIT STRUCT TO EEG STRUCT
        fprintf('%s) dipfit available: %i\n',EEG.subject,isfield(EEG,'dipfit'));
        %## Update EEG channel location again
        if ~isempty(CHANLOCS_FPATHS)
            fprintf('Attaching customelectrode locations...\n');
            EEG = custom_update_chanlocs(EEG,CHANLOCS_FPATHS{subj_i});
        end
        if DO_BEM_DIPFIT
            fprintf('Trying Boundary Element Model Dipfit...\n');
            try
                EEG.dipfit.model;
                EEG.dipfit.coord_transform;
                EEG.dipfit.mrifile;
                EEG.dipfit.hdmfile;
                EEG.dipfit.coordformat;
            catch e
                fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
                EEG = bem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005,DIP_NUM,DIP_PLOT)
            end
        end
        if DO_FEM_DIPFIT
            fprintf('Trying Finite Element Model Dipfit...\n');
            try
                EEG.dipfit.model;
            catch e
                fprintf('A valid dipfit model needs to be calculated for FEM dipfit analysis...\n');
                fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n',...
                     'report. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
                try
                    fprintf('Trying to load .mat file...\n');
                    out = par_load(fPaths{subj_i},'dipfit_fem_norm.mat');
                    EEG.dipfit = out;
                catch e
                    error('failed.\n CALLBACK:\n%s\n',getReport(e));
                end
            end
            [EEG] = fem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005);
        end
                 
        %## CHECK ICLABEL
        if ~isfield(EEG.etc,'ic_classification')
            EEG = iclabel(EEG, 'lite');
        end
        %## MAKE EEG STRUCTURES FOR STUDY
        EEG.filepath = fPath;
        EEG.filename = fName;
        EEG.group = char(groups{subj_i});
        EEG.condition = char(conditions{subj_i});
        EEG.session = char(sessions{subj_i});
        EEG = eeg_checkset(EEG,'eventconsistency');
        if DO_SAVE_ICA
            fprintf(1,'Saving Subject %s\n',EEG.subject);
            [EEG] = pop_saveset(EEG,'savemode','twofiles',...
                'filename',fName,...
                'filepath',fPath);
        end
    else
        EEG = pop_loadset('filepath',fPath,'filename',fName);
    end
    ALLEEG{subj_i} = EEG;
end
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(ALLEEG));
for subj_i = 1:length(ALLEEG)
    fss{subj_i} = fields(ALLEEG{subj_i})';
    disp(size(fields(ALLEEG{subj_i})'));
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(ALLEEG)
    EEG = ALLEEG{subj_i};
    fs = fields(EEG);
    %- add fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false);
    out = [out{:}];
    addFs = fsPrev(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j});
        end
    end
    ALLEEG{subj_i} = orderfields(EEG);
end
%- CONCATENATE ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
end
function [EEG] = fem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005)
    tmp = [];
    if ~isfield(EEG.dipfit,'coord_transform')
        EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1]; %COORD_TRANSFORM_MNI;
        tmp = [tmp, 'added default coord_transform; '];
    end
    if ~isfield(EEG.dipfit,'mrifile')
        EEG.dipfit.mrifile = MNI_MRI;
        tmp = [tmp, 'added default mrifile; '];
    end
    if ~isfield(EEG.dipfit,'hdmfile')
        EEG.dipfit.hdmfile = MNI_VOL;
        tmp = [tmp, 'added default hdmfile; '];
    end
    if ~isfield(EEG.dipfit,'coordformat')
        EEG.dipfit.coordformat = 'MNI';
        tmp = [tmp, 'added default coordformat; '];
    end
    if ~isfield(EEG.dipfit,'chanfile')
        EEG.dipfit.chanfile = MNI_CHAN_1005;
        tmp = [tmp, 'added default chanfile; '];
    end
    if ~isfield(EEG.dipfit,'chansel')
        EEG.dipfit.chansel = (1:EEG.nbchan);
        tmp = [tmp, 'added default chansel; '];
    end
end
function [EEG] = custom_update_chanlocs(EEG,chanlocs_fPath)
    tmp = load(chanlocs_fPath);
    chanlocs_new = tmp.chanlocs_new;
    nodatchans_new = tmp.nodatchans_new;
    %- update the EEG electrode locations
    % Be cautious that not all electrodes are EEG
    % Sanity check: if we have 120 electrodes digitized
    fprintf('Found total of %i electrodes',length(chanlocs_new));
    for p = 1:length(chanlocs_new)
        elec_idx = find(strcmpi(chanlocs_new(p).labels,{EEG.chanlocs(:).labels}));
        if ~isempty(elec_idx)
            % update all available fields
            EEG.chanlocs(elec_idx).X = chanlocs_new(p).X;
            EEG.chanlocs(elec_idx).Y = chanlocs_new(p).Y;
            EEG.chanlocs(elec_idx).Z = chanlocs_new(p).Z;
            EEG.chanlocs(elec_idx).theta = chanlocs_new(p).theta;
            EEG.chanlocs(elec_idx).radius = chanlocs_new(p).radius;
            EEG.chanlocs(elec_idx).sph_theta = chanlocs_new(p).sph_theta;
            EEG.chanlocs(elec_idx).sph_phi = chanlocs_new(p).sph_phi;
        end
    end
    % Add fiducials location 
    if isempty(EEG.chaninfo.nodatchans)
        EEG.chaninfo.nodatchans = nodatchans_new;
    end
    EEG = eeg_checkchanlocs(EEG); % check the consistency of the chanloc structure
end

function [EEG] = bem_eeglab_dipfit(EEG,COORD_TRANSFORM_MNI,MNI_MRI,MNI_VOL,MNI_CHAN_1005,DIP_NUM,DIP_PLOT)
    %## FIT DIPOLES TO HEADMODEL IF NEEDED
    fprintf('MNI pop_dipfit_settings...\n');
    %## Check Files
    tmp = [];
    if ~isfield(EEG.dipfit,'coord_transform')
        EEG.dipfit.coord_transform = COORD_TRANSFORM_MNI;
        tmp = [tmp, 'added default coord_transform; '];
    end
    if ~isfield(EEG.dipfit,'mrifile')
        EEG.dipfit.mrifile = MNI_MRI;
        tmp = [tmp, 'added default mrifile; '];
    end
    if ~isfield(EEG.dipfit,'hdmfile')
        EEG.dipfit.hdmfile = MNI_VOL;
        tmp = [tmp, 'added default hdmfile; '];
    end
    if ~isfield(EEG.dipfit,'coordformat')
        EEG.dipfit.coordformat = 'MNI';
        tmp = [tmp, 'added default coordformat; '];
    end
    if ~isfield(EEG.dipfit,'chanfile')
        EEG.dipfit.chanfile = MNI_CHAN_1005;
        tmp = [tmp, 'added default chanfile; '];
    end
    if ~isfield(EEG.dipfit,'chansel')
        EEG.dipfit.chansel = (1:EEG.nbchan);
        tmp = [tmp, 'added default chansel; '];
    end
    EEG.dipfit.comment = tmp;
   
    %- pop_multifit.m
    %- DIPFIT (see. ft_dipolefitting())
    fprintf('pop_multifit...\n');
    EEG = pop_multifit(EEG,[],'dipoles',DIP_NUM,'dipplot',DIP_PLOT);
%     dipfit = EEG.dipfit;
%     save([fPath filesep sprintf('%s_dipfit_mni.mat',EEG.subject)],'dipfit');
end
