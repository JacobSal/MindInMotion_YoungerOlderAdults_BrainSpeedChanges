function [SUBJSTRUCT] = subjstruct_calcDipoles(PATHS,SUBJSTRUCT,studyName,processNum,varargin)
%SUBJ_CALCDIPOLES Summary of this function goes here
%   
%   IN: 
%   OUT: 
%   IMPORTANT: 
%   (04/04/2022) in consideration of time I'm deciding to forego
%   reoptimizing this function till a later date. 

%## GLOBALS

%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {[]};
p = inputParser;
%## REQUIREDfd
addRequired(p,'PATHS',@isstruct)
addRequired(p,'SUBJSTRUCT',@isstruct)
addRequired(p,'studyName',@ischar);
addRequired(p,'processNum',@isnumeric)
%## OPTIONALS
validPath = (@(x) ischar(x) || isempty(x));
addOptional(p,'path_ext',Defaults{1},validPath);
%## PARAMETERS
parse(p, PATHS, SUBJSTRUCT, studyName, processNum, varargin{:});
PATH_EXT = p.Results.path_ext;
%## SET DEFAULTS
SAVE_PROCESS_NUM = 1;
try
    if strncmp(computer,'PC',2)
        DO_UNIX = false;
    else
        DO_UNIX = true;
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end
%% ===================================================================== %%
%## PACKAGE IMPORTS?
% (09/10/2022) was having error with ft_version (resulting from a call to
% ft_defaults) by hard coding the matlab search path the error seems to be
% fixed. No promises...
%## DIARY
diaryPath = [PATHS.path4diaries filesep studyName '1_diary_calcDipoles.txt'];
if ~exist([PATHS.path4diaries filesep studyName],'dir')
    mkdir([PATHS.path4diaries filesep studyName]);
end
diary(diaryPath)
%## EEGLAB
eeglab;
%% PARFOR/FOR LOOP
fprintf(1,'==== Starting Dipole Analysis ====\n')
TMP = cell(1,length(SUBJSTRUCT));
parfor subj_i      = 1:length(SUBJSTRUCT)
% for subj_i = 1:length(SUBJSTRUCT)
    path_ext = PATH_EXT;
    % ==== PARFOR PARAMS ==== %
    if isempty(path_ext) 
        if DO_UNIX
            path_ext = 'dferris';
        else
            path_ext = 'M';
        end
    end
    SUBJ = SUBJSTRUCT(subj_i);
    %## PATHS
    subjStr     = SUBJ.subjStr;
    fPath = SUBJ.pathsEEG.ica(processNum).filepath;
    fName = SUBJ.pathsEEG.ica(processNum).filename;
    if DO_UNIX
        fPath = convertPath2UNIX(fPath,path_ext);
    else
        fPath = convertPath2Drive(fPath,path_ext);
    end    
    if ~exist(fullfile(fPath,fName),'file')
        error('No ICA assigned for subject %s\nPath: %s',subjStr,fullfile(fPath,fName));
    end
    % ==== END: PARFOR PARAMS ==== %
    %## LOAD DATA    
    fprintf(1,'loading subject %s ICA EEG from %s\n',subjStr,fPath);
    [~,EEG,~] = subj_loadICA(SUBJ,processNum,fName,fPath);
    %## CREATE/LOAD HEADMODEL     
    fprintf(1,'Attaching files for headmodel...\n');
    [SUBJ] = my_MNI4DIPFIT(PATHS, SUBJ); 
    %## DIPFIT (SOURCE LOCALIZATION)
    fprintf(1,'Running pop_multifit()...\n');
    [EEG] = my_MNIDipFit(PATHS, SUBJ, EEG); 
    %- hardcode save params    
    fPath = SUBJ.pathsEEG.sourcelocalize(SAVE_PROCESS_NUM).filepath;
    fName = [subjStr '_TMPEEG.mat'];
    fprintf(1,'saving subject %s ICA EEG from %s\n',subjStr,fPath);
    SUBJ.pathsEEG.sourcelocalize(processNum).filepath = fPath;
    SUBJ.pathsEEG.sourcelocalize(processNum).filename = fName;
    par_save(EEG,fPath,fName);
    %## PARFOR STRUCT HANDLER
    TMP{subj_i} = SUBJ;
end
%## extract SUBJSTRUCT
SUBJSTRUCT = cellfun(@(x) [[]; x], TMP);
%## saving SUBJSTRUCT
% subj_saveStruct(SUBJSTRUCT)
diary('off');
end

