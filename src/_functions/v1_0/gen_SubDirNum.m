function [outPath,subDirNum] = gen_SubDirNum(directory,varargin)
%GENSUBDIRNUM Summary of this function goes here
%
%   IN:
%       REQUIRED:
%           directory, CHAR
%               path to the folder where you'd like to make a sub folder
%       OPTIONAL:
%       PARAMETER:
%   OUT:
%       OUTVAR, VAR
%           
%   IMPORTANT:
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {true,...
            []};
p = inputParser;
validVal = (@(x) isnumeric(x) || isempty(x));
%## REQUIRED
addRequired(p,'directory',@ischar)
%## OPTIONAL
addOptional(p,'createDir',Defaults{1},@islogical)
%## PARAMETER
addParameter(p,'override',Defaults{2},validVal)
parse(p, directory, varargin{:});
%## SET DEFAULTS
createDir = p.Results.createDir;
override = p.Results.override;
%% ===================================================================== %%
subFolderExists = 1;
subFolderEmpty = 1;
subDirNum = 0;
while subFolderExists == 1 && subFolderEmpty == 1
    subDirNum = subDirNum + 1;
    % check if folder exists
    pathIn = fullfile(directory, num2str(subDirNum));
    subFolderExists = isfolder(pathIn);
    % check if folder is empty
    subFolderEmpty = dir(pathIn); 
    subFolderEmpty = ~isempty(subFolderEmpty(cellfun(@(x) x~=0, {subFolderEmpty.bytes})));
end
if ~isempty(override)
    subDirNum = override;
end
outPath = fullfile(directory, num2str(subDirNum));
if createDir
    mkdir(outPath);
end
end

