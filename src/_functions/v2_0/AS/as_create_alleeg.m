function [ALLEEG] = as_create_alleeg(fNames,fPaths,subjectNames,save_dir,varargin)
%MIM_CREATE_ALLEEG Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## DEFINE DEFAULTS
%- developer params
DO_BEM_DIPFIT = false;
FORCE_RELOAD = true;
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
SAVE_EEG = false;
%- CHANLOCS_FPATHS
CHANLOCS_FPATHS = {};
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}. Full file paths for chanlocs (e.g., /path/to/subject/dir/chanlocs.mat)';
cfp_validFcn = @(x) assert((ischar(x{1}) && iscell(x) && length(x) == length(fNames) || isempty(x)), errorMsg);

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
addParameter(p,'SAVE_EEG',SAVE_EEG,@islogical);
addParameter(p,'CHANLOCS_FPATHS',CHANLOCS_FPATHS,cfp_validFcn);
parse(p,fNames,fPaths,subjectNames,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
conditions          = p.Results.conditions;
groups              = p.Results.groups;
sessions            = p.Results.sessions;
bool_save_eeg       = p.Results.SAVE_EEG;
chanlocs_fPaths     = p.Results.CHANLOCS_FPATHS;
%% ===================================================================== %%
%## CREATE STUDY
fprintf(1,'==== Creating Study ====\n')
%* empty ALLEEG structure for repopulating
ALLEEG = cell(1,length(fNames)); 
%## Populate ALLEEG Struct
pp = gcp('nocreate');
disp(pp);
if ~isfield(pp,'NumWorkers')
    POOL_SIZE = 1;
else
    POOL_SIZE = pp.NumWorkers;
end
%## Populate ALLEEG Struct
parfor (subj_i=1:length(fNames),POOL_SIZE)
% for subj_i=1:length(fNames)
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
        if ~isempty(chanlocs_fPaths)
            fprintf('Attaching customelectrode locations...\n');
            EEG = custom_update_chanlocs(EEG,chanlocs_fPaths{subj_i});
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
        if bool_save_eeg
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
% %## BOOKKEEPING (i.e., delete fields not similar across EEG structures)
% fsPrev = {};
% for subj_i = 1:length(ALLEEG)
%     EEG = ALLEEG{subj_i};
%     fs = fields(EEG);
%     % delete fields not present in other structs.
%     out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
%     out = [out{:}];
%     delFs = fs(~out);
%     if ~isempty(fsPrev) && any(~out)
%         for j = 1:length(delFs)
%             EEG = rmfield(EEG,delFs{j});
%             fprintf("%s) Removing fields %s",EEG.subject,delFs{j})
%         end
%     else
%         fsPrev = fs;
%     end
%     ALLEEG{subj_i} = EEG;
% end
% %- CONCATENATE ALLEEG
% ALLEEG = cellfun(@(x) [[]; x], ALLEEG);

%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(ALLEEG));
for subj_i = 1:length(ALLEEG)
    fss{subj_i} = fields(ALLEEG{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(ALLEEG)
    EEG = ALLEEG{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding %s %s\n',EEG.subject,addFs{j})
        end
    end 
    ALLEEG{subj_i} = EEG;
end
%- CONCATENATE ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
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
