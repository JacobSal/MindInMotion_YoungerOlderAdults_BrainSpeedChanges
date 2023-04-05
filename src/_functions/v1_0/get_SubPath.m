function [pathOut,lcl] = get_SubPath(PATHS,basePath,subPath,varargin)
%GET_SUBPATH nice function for making subfolder paths
%   
%   IN: 
%   OUT: 
%   IMPORTANT:
% 
%## GLOBALS
%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {[]};
p = inputParser;
%## REQUIRED
verify_path = (@(x) ischar(x) || isempty(x));
addRequired(p,'PATHS',@isstruct)
addRequired(p,'basePath',verify_path)
addRequired(p,'subPath',verify_path)
%## PARAMETERS
%## OPTIONAL
verify_ext = (@(x) ischar(x) || isempty(x));
addOptional(p,'path_ext',Defaults{1},verify_ext);
%## PARSE
parse(p, PATHS, basePath, subPath, varargin{:});
%## SET DEFAULTS
path_ext = p.Results.path_ext;
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
if isempty(basePath)
    lcl = PATHS.localStorage;      
else
    lcl = basePath;
end
%- convert path to current OS
if DO_UNIX
    if isempty(path_ext)
        path_ext = 'dferris';
    end
    lcl = convertPath2UNIX(lcl,path_ext);
else
    if isempty(path_ext)
        path_ext = 'M';
    end
    lcl = convertPath2Drive(lcl,path_ext);
end
pth = fullfile(lcl,subPath);
if ~isfolder(pth)
    fprintf('Making Folder %s...',pth)
    mkdir(pth)
end
pathOut = pth;  
end

