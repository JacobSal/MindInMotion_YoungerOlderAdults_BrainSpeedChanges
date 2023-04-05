function [SUBJSTRUCT] = subj_saveStruct(SUBJSTRUCT,varargin)
%SUBJ_RESAVESTRUCT Summary of this function goes here
%   
%   IN:
%       SUBJSTRUCT
% 
%   OUT: 
%   IMPORTANT: 
%

%## TIME
tic
%## DEFINE DEFAULTS
%- validation
validPath = (@(x) ischar(x) || isempty(x));
Defaults = {[]};
p = inputParser;
%## REQUIRED
addRequired(p,'SUBJSTRUCT',@isstruct)
%## OPTIONAL
addOptional(p,'newfPath',Defaults{1},validPath)
addOptional(p,'newfName',Defaults{1},validPath)
addOptional(p,'path_prefix',Defaults{1},validPath)
parse(p,SUBJSTRUCT, varargin{:});
newfPath = p.Results.newfPath;
newfName = p.Results.newfName;
path_prefix = p.Results.path_prefix;
%## SET DEFAULTS
%% ===================================================================== %%
%## MAIN PROCESS
% inputPath = SUBJSTRUCT.PATHS.studyPath;
if ~isempty(newfPath)
    fPath = newfPath;
else
    fPath = SUBJSTRUCT(1).SAVE.filepath;
end
if ~isempty(newfName)
    fName = newfName;
    SUBJSTRUCT(1).SAVE.filename = fName;
else
    fName = SUBJSTRUCT(1).SAVE.filename;
end
inputPath = [fPath filesep ];
SUBJSTRUCT(1).SAVE.filepath = fPath;
if ~exist(fPath,'dir')
    mkdir(fPath);
end
% determine current os for SUBJSTRUCT save path
if strncmp(computer,'PC',2)
    if isempty(path_prefix)
        path_prefix = 'M';
    end
    par_save(SUBJSTRUCT,fPath,fName,path_prefix,[])
else  % isunix
    if isempty(path_prefix)
        path_prefix = 'dferris';
    end
    par_save(SUBJSTRUCT,fPath,fName,path_prefix,[])
end

fprintf('Saved SUBJSTRUCT to %s\n',inputPath);
end

