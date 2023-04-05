function [SUBJSTRUCT] = subjstruct_epoch(PATHS,SUBJSTRUCT,processNum,studyName,parseVar,trialType,epochTimeLimits,varargin)
%EPOCHSUBJS Summary of this function goes here
%   IN: 
%       %## Required
%       PATHS, STRUCT
%           
%       SUBJSTRUCT, STRUCT
%           
%       processNum, NUMERIC
%           
%       studyName, CHAR
%           
%       parseVar, 
%           
%       trialType
%           
%       %## Optional
%   OUT: 
%   IMPORTANT: 
% diary './2_GLOBAL_BATCH/_diaries/subjstruct_epoch'
%## TIME
tic
validPath = (@(x) ischar(x) || isempty(x));
%## DEFINE DEFAULTS
Defaults = {3*60,...
            0.0,...
            'TrialStart',...
            'TrialEnd',...
            'type',...
            'cond',...
            []};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'parseVar',@ischar);
addRequired(p,'trialType',@ischar);
addRequired(p,'epochTimeLimits',@isnumeric);
%## OPTIONAL
addOptional(p,'path_ext',Defaults{7},validPath)
addOptional(p,'trialLength',Defaults{1},@isnumeric);
addOptional(p,'percentOverlap',Defaults{2},@isnumeric);
addOptional(p,'trialBeginStr',Defaults{3},@ischar);
addOptional(p,'trialEndStr',Defaults{4},@ischar);
addOptional(p,'eventTrialParser',Defaults{5},@ischar);
addOptional(p,'eventCondParser',Defaults{6},@ischar);
%## PARAMETER
addParameter(p,'savePath',Defaults{7},validPath)
%## SET DEFAULTS
parse(p,PATHS,SUBJSTRUCT,processNum,studyName,parseVar,trialType,epochTimeLimits,varargin{:})
savePath = p.Results.savePath;
path_ext = p.Results.path_ext;
TRIAL_LENGTH = p.Results.trialLength;
PER_OVERLAP = p.Results.percentOverlap;
TRIAL_BEGIN_STR = p.Results.trialBeginStr;
TRIAL_END_STR = p.Results.trialEndStr;
EVENT_TRIAL_PARSER = p.Results.eventTrialParser;
EVENT_COND_PARSER = p.Results.eventCondParser;
try
    if strncmp(computer,'PC',2)
        DO_UNIX = false;
    else
        DO_UNIX = true;
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end
SAVE_PROCESS_NUM = 1;
%% --------------------------------------------------------------------- %%
%## PACKAGE IMPORTS?
% (09/10/2022) was having error with ft_version (resulting from a call to
% ft_defaults) by hard coding the matlab search path the error seems to be
% fixed. No promises...
%## DIARY
diaryPath = [PATHS.path4diaries filesep studyName filesep '3_diary_subjstruct_epoch.txt'];
if ~exist([PATHS.path4diaries filesep studyName],'dir')
    mkdir([PATHS.path4diaries filesep studyName]);
end
diary(diaryPath)
%## EEGLAB
eeglab;
%% PARFOR LOOP
TMP = cell(length(SUBJSTRUCT),1);
parfor subj_i = 1:length(SUBJSTRUCT)
    SUBJ = SUBJSTRUCT(subj_i);
    %% EPOCHING LOOP
    % ==== PARAMS ==== %
    %- SUBJ
    subjStr = SUBJ.subjStr;
    %## SAVE PATHS
    %- load path for source localized EEG
    fName = SUBJ.pathsEEG.sourcelocalize(processNum).filename;
    if DO_UNIX
        fPath = convertPath2UNIX(SUBJ.pathsEEG.sourcelocalize(processNum).filepath,path_ext);
    else
        fPath = convertPath2Drive(SUBJ.pathsEEG.sourcelocalize(processNum).filepath,path_ext);
    end
    %- save path for epoched EEG
    if ~isempty(SUBJ.pathsEEG.epoched(SAVE_PROCESS_NUM).filepath)
        sP = [subjStr filesep 'Epoched' filesep parseVar filesep trialType];
        [epochSavePath,~] = get_SubPath(PATHS,savePath,sP);
    end
    % ==== END: PARAMS ==== %    
    fprintf(1,'\nEpoching EEG data for Subject %s...\n',subjStr)
    fprintf(1,'Loading Subject %s\n',subjStr);    
    %## EEG FILE LOAD HANDLER
    EEG = par_load(fPath,fName,[]);
    %## EPOCH SUBJ DATA
    fprintf(1,'Parsing EEG data trials...\n')
    [EEG] = subj_epochData(PATHS,EEG,parseVar,trialType,epochTimeLimits,...
            TRIAL_LENGTH,PER_OVERLAP,TRIAL_BEGIN_STR,TRIAL_END_STR,EVENT_TRIAL_PARSER,EVENT_COND_PARSER);
    %- Update SUBJSTRUCT for saving
    SUBJ.pathsEEG.epoched(processNum).label = sprintf('Epoch_%s-%s',parseVar,trialType)
    SUBJ.pathsEEG.epoched(processNum).filepath = epochSavePath;
    SUBJ.pathsEEG.epoched(processNum).filename = [subjStr '_EPOCH_TMPEEG.set'];
    %- Save EEG
    fprintf(1,'Saving Subject %s\n',subjStr); 
    [~] = pop_saveset(EEG ,'savemode','twofiles',...
        'filename',SUBJ.pathsEEG.epoched(processNum).filename,...
        'filepath',SUBJ.pathsEEG.epoched(processNum).filepath);
    %## assign META data to SUBJSTRUCT
    SUBJ.META.epoched.parseVar = parseVar;
    SUBJ.META.epoched.trialType = trialType;
    SUBJ.META.epoched.epochTimeLimits = epochTimeLimits;
    %## collect changes to SUBJSTRUCT
    TMP{subj_i} = SUBJ;
end
%## extract SUBJSTRUCT
SUBJSTRUCT = cellfun(@(x) [[]; x], TMP);
diary('off');
end
