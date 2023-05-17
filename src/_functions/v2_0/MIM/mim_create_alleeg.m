function [ALLEEG] = mim_create_alleeg(fNames,fPaths,subjectNames,save_dir,varargin)
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
FORCE_RELOAD = false;
POOL_SIZE = 10;
ICA_FNAME_REGEXP = '%s_allcond_ICA_TMPEEG.set';
%- find eeglab on path
tmp = strsplit(path,';');
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
%* empty STUDY structure for repopulating
fsPrev = {};
%## Populate ALLEEG Struct
parfor (subj_i=1:length(fNames),POOL_SIZE)
% for subj_i=1:length(fNames)
    subjStr = subjectNames{subj_i};
    fName = sprintf(ICA_FNAME_REGEXP,subjStr);
    fPath = [save_dir filesep subjStr filesep 'ICA'];
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    if ~exist([fPath filesep fName],'file') || FORCE_RELOAD
        fprintf(1,'Loading Subject %s\n',subjStr)
        [~,EEG,~] = eeglab_loadICA(fNames{subj_i},fPaths{subj_i});
        %- override fPath & fName
        EEG.subject = subjStr;
        %## ADD DIPFIT STRUCT TO EEG STRUCT
        fprintf('%s) dipfit available: %i\n',EEG.subject,isfield(EEG,'dipfit'));
        try
            EEG.dipfit.coord_transform;
            EEG.dipfit.mrifile;
            EEG.dipfit.hdmfile;
            EEG.dipfit.coordformat;
        catch e
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
            fprintf('MNI pop_dipfit_settings...\n');
            %{
            EEG = pop_dipfit_settings(EEG,'coordformat','MNI','coord_transform',COORD_TRANSFORM_MNI,...
                    'hdmfile',MNI_VOL,'mrifile',MNI_MRI,'chanfile',MNI_CHAN_1005);
            %}

            %## Check Filese
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
        end
        %% Update EEG channel location again %has to update chanloc for dipfit
        tmp = load(chanlocs_fPaths{subj_i});
        chanlocs_new = tmp.chanlocs_new;
%         getchanlocs_new = tmp.getchanlocs_new;
        nodatchans_new = tmp.nodatchans_new;
        % update the EEG electrode locations
        % Be cautious that not all electrodes are EEG
        % Sanity check: if we have 120 electrodes digitized
        disp(['Found total ',num2str(length(chanlocs_new)),' electrodes']);
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
        %%
        try 
            EEG.dipfit.model;
        catch e
            %- DIPFIT (see. ft_dipolefitting())
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n'],e.identifier,e.message,EEG.subject);
            fprintf('pop_multifit...\n');
            EEG = pop_multifit(EEG,[],'dipoles',DIP_NUM,'dipplot',DIP_PLOT);
            tmp = EEG.dipfit;
            par_save(tmp,fPath,sprintf('%s_dipfit_mni.mat',EEG.subject));
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
        fprintf(1,'Saving Subject %s\n',EEG.subject);
        if bool_save_eeg
            [EEG] = pop_saveset(EEG,'savemode','twofiles',...
                'filename',fName,...
                'filepath',fPath,...
                'version','6');
        end
    else
        EEG = pop_loadset('filepath',fPath,'filename',fName);
    end
    ALLEEG{subj_i} = EEG; 
end
%## BOOKKEEPING (i.e., delete fields not similar across EEG structures)
for subj_i = 1:length(ALLEEG)
    EEG = ALLEEG{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    delFs = fs(~out);
    if ~isempty(fsPrev) && any(~out)
        for j = 1:length(delFs)
            EEG = rmfield(EEG,delFs{j});
            fprintf("%s) Removing fields %s",EEG.subject,delFs{j})
        end
    else
        fsPrev = fs;
    end
    ALLEEG{subj_i} = EEG;
end
%- CONCATENATE ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
end

