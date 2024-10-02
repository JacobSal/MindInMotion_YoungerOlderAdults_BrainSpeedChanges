function [OUTVAR] = par_load(fPath,fName,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           fPath, CHAR
%               path to the folder where your file is held
%           fName, CHAR
%               file name & extension (e.g., 'INEEG.mat')
%       OPTIONAL:
%           path_ext, CHAR (default: [])
%               extension for operating system conversion see.
%               convertPath2Drive.m & convertPath2UNIX.m
%       PARAMETER:
%   OUT:
%       OUTVAR, VAR
%           
%   IMPORTANT:
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {[]};
p = inputParser;
%## REQUIRED
verify_fName = (@(x) ischar(x) || isempty(x));
addRequired(p,'fPath',@ischar)
addRequired(p,'fName',verify_fName)
%## OPTIONAL
verify_ext = (@(x) ischar(x) || isempty(x));
addOptional(p,'path_ext',Defaults{1},verify_ext);
%## PARAMETER
%## PARSE
parse(p, fPath, fName, varargin{:});
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
FILESEP_UNIX = '/';
FILESEP_DOS = '\';
%% ===================================================================== %%
if isempty(fName)
    sunix = contains(fPath,FILESEP_UNIX);
    sdos = contains(fPath,FILESEP_DOS);
    if sunix
        spath = strsplit(fPath,FILESEP_UNIX);
        jpath = join(spath(1:end-1)',FILESEP_UNIX);
        fName = spath{end};
        fPath = jpath{1};
    elseif sdos
        spath = strsplit(fPath,FILESEP_DOS);
        jpath = join(spath(1:end-1)',FILESEP_DOS);
        fName = spath{end};
        fPath = jpath{1};
    end
end

if DO_UNIX
    %- convert path
    pathIn = convertPath2UNIX(fPath);
    %- load
    fprintf('\nLoading %s from\n%s\n',fName,pathIn);
    tmp = load([pathIn filesep fName]); 
    %- extract
    fs = fields(tmp);
    OUTVAR = tmp.(fs{1});
else
    %- convert path
    pathIn = convertPath2Drive(fPath);
    %- load
        fprintf('\nLoading %s from\n%s\n',fName,pathIn);
        tmp = load([pathIn filesep fName]);
    %- extract
    fs = fields(tmp);
    OUTVAR = tmp.(fs{1});  
end
end

