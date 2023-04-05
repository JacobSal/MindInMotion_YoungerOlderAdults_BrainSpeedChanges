function [SUBJSTRUCT,study_subDirNum] = subj_resaveStruct(SUBJSTRUCT,genNewSubDirNum,varargin)
%SUBJ_RESAVESTRUCT Summary of this function goes here
%   
%   IN:
%       startPath
% 
%   OUT: 
%   IMPORTANT: 
%

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {[],...
            'M'};
p = inputParser;
%## REQUIRED
addRequired(p,'SUBJSTRUCT',@ischar)
addRequired(p,'newPath',@ischar)
%## OPTIONAL

parse(p,SUBJSTRUCT, genNewSubDirNum, varargin{:});
%## SET DEFAULTS
subDirName = p.Results.subDirName;
path_prefix = p.Results.path_prefix;
%% ===================================================================== %%
%## MAIN PROCESS

% determine current os for SUBJSTRUCT save path
sunix = contains(inputPath,'/');
sdos = contains(inputPath,'\');
if sdos
    currSep = '\';
elseif sunix
    currSep = '/';
end
try
    if strncmp(computer,'PC',2)
        osName = 'DOS';
    else
        osName = 'UNIX';
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end

switch osName
    case 'DOS'
        dosSep = '\';
        %## reassigning savepath for SUBJSTRUCT
        tmpPth = split(savePath,currSep);
        tmpPth = join(tmpPth(1:(end-1)),dosSep); 
        tmpPth = tmpPth{1};
        if ~isempty(subDirName)
            newPath = [tmpPth dosSep subDirName dosSep];
        else
            newPath = [tmpPth dosSep];
        end
        %## generate a new subdirectory number (generic double)
        if genNewSubDirNum
            [~,study_subDirNum] = genSubDirNum(convertPath2Drive(newPath,'M'),false);
        end
        newPath =  [newPath dosSep num2str(study_subDirNum) dosSep];
        %## make directory
        mkdir(convertPath2Drive(newPath,path_prefix)) %** DRIVE
        %## append file name
        tmpPth = split(savePath,currSep);
        newSave = fullfile(newPath,tmpPth{end});
        [SUBJSTRUCT.PATHS.studyPath] = newSave;
    case 'UNIX'
        unixSep = '\';
        %## reassigning savepath for SUBJSTRUCT
        tmpPth = split(savePath,currSep);
        tmpPth = join(tmpPth(1:(end-1)),unixSep); 
        tmpPth = tmpPth{1};
        if ~isempty(subDirName)
            newPath = [tmpPth unixSep subDirName unixSep];
        else
            newPath = [tmpPth unixSep];
        end
        %## generate a new subdirectory number (generic double)
        if genNewSubDirNum
            [~,study_subDirNum] = genSubDirNum(convertPath2Drive(newPath,'M'),false);
        end
        newPath =  [newPath unixSep num2str(study_subDirNum) unixSep];
        %## make directory
        mkdir(convertPath2Drive(newPath,path_prefix)) %** DRIVE
        %## append file name
        tmpPth = split(SUBJSTRUCT(1).PATHS.studyPath,currSep);
        newSave = fullfile(newPath,tmpPth{end});
        [SUBJSTRUCT.PATHS.studyPath] = newSave;
end
end

