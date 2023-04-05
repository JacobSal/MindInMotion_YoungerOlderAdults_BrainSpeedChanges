function [] = par_save(SAVEVAR,path,fName,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           SAVEVAR, variable to save.
%               STRUCT, CHAR, DICT, ARRAY you want to save.
%           path, CHAR
%               path to the folder where your file is held
%           fName, CHAR
%               file name & extension (e.g., 'INEEG.mat')
%           osName, CHAR
%               Operating system load format (options: 'DOS', 'UNIX')
%       OPTIONAL:
%           ext, CHAR
%               extension for operating system conversion see.
%               convertPath2Drive.m & convertPath2UNIX.m
%       PARAMETER:
%   OUT:
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {[]};
verify_ext = (@(x) ischar(x) || isempty(x));
p = inputParser;
%## REQUIRED

addRequired(p,'SAVEVAR')
addRequired(p,'path',@ischar)
addRequired(p,'fName',verify_ext)
%## OPTIONAL
addOptional(p,'path_ext',Defaults{1},verify_ext);
addOptional(p,'fname_ext',Defaults{1},verify_ext);
%## PARAMETER
parse(p, SAVEVAR, path, fName, varargin{:});
%## SET DEFAULTS
path_ext = p.Results.path_ext;
fname_ext = p.Results.fname_ext;
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
if DO_UNIX
    % set default extension
    if isempty(path_ext)
        path_ext = 'dferris';
    end
    % convert path to os
    pathIn = convertPath2UNIX(path,path_ext);
    if ~isempty(fname_ext)
        fnames = strsplit(fName,'.');
        fnames{1} = [fnames{1} fname_ext];
        fName = join(fnames,'.');
        fName = fName{1};
    end
    % set save path
    %- if filename is included in path
    if ~isempty(fName)
        savePath = [pathIn filesep fName];
    else
        savePath = pathIn;
    end
    % save
    save(savePath,'SAVEVAR','-v7.3');
    fprintf('\nSaving %s to\n%s\n',fName,savePath);
else
    %- set default extension
    if isempty(path_ext)
        path_ext = 'M';
    end
    %- convert path to os
    pathIn = convertPath2Drive(path,path_ext);
    %- set save path
    if ~isempty(fname_ext)
        fnames = strsplit(fName,'.');
        fnames{1} = [fnames{1} fname_ext];
        fName = join(fnames,'.');
        fName = fName{1};
    end
    %- if filename is included in path
    if ~isempty(fName)
        savePath = [pathIn filesep fName];
    else
        savePath = pathIn;
    end
    %- save
    save(savePath,'SAVEVAR','-v7.3');
    fprintf('\nSaving %s to\n%s\n',fName,savePath);
end
