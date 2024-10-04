function [path_out] = convertPath2UNIX(fPath,varargin)
%prepHPG4AMICA
%   IN: 
%       fPath, CHAR
%           
%       HPC_DRIVE_NAME, CHAR
%           
%       HPC_QOS_NAME, CHAR
%           
%   OUT: 
%       path_out
%           
%   IMPORTANT: 
%       
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/06/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%## TIME
tic
%## DEFINE DEFAULTS
%- 
HPC_DRIVE_NAME = 'blue';
errorMsg = 'Value must be CHAR. This value is the name of the STORAGE FOLDER on HiperGator.'; 
hdn_validFcn = @(x) assert(ischar(x),errorMsg);
%- 
HPC_QOS_NAME = 'dferris';
errorMsg = 'Value must be CHAR. This value is the name of the USER FOLDER on HiperGator.'; 
hqn_validFcn = @(x) assert(ischar(x),errorMsg);
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'fPath',@ischar);
%## OPTIONAL
addOptional(p,'HPC_DRIVE_NAME',HPC_DRIVE_NAME,hdn_validFcn);
addOptional(p,'HPC_QOS_NAME',HPC_QOS_NAME,hqn_validFcn);
%## PARAMETER
%## SET DEFAULTS
parse(p,fPath,varargin{:});
HPC_DRIVE_NAME = p.Results.HPC_DRIVE_NAME;
HPC_QOS_NAME = p.Results.HPC_QOS_NAME;
%% ===================================================================== %%
FILESEP_UNIX = '/'; %UNIX override
FILESEP_DOS = '\';
path_out = [];
if contains(fPath,FILESEP_UNIX) && contains(fPath,FILESEP_DOS)
    fPath = fullfile(fPath);
end
unix_sep = strfind(fPath,FILESEP_UNIX);
dos_sep = strfind(fPath,FILESEP_DOS);
%##  if there are unix seps in path && the UNIX initiation sep does not
%   exist
if ~isempty(unix_sep) && ~any(unix_sep==1)
    path_out = strsplit(fPath,FILESEP_UNIX);
    path_out = join(path_out(2:end),FILESEP_UNIX);
    if ~isempty(HPC_QOS_NAME)
        path_out = [FILESEP_UNIX HPC_DRIVE_NAME FILESEP_UNIX HPC_QOS_NAME FILESEP_UNIX path_out{1}];
    else
        path_out = [FILESEP_UNIX HPC_DRIVE_NAME FILESEP_UNIX path_out{1}];
    end
%##  if there are dos seps in path && the UNIX initiation sep does not
%   exist
elseif ~isempty(dos_sep) && ~any(dos_sep==1)
    path_out = strsplit(fPath,FILESEP_DOS);
    path_out = join(path_out(2:end),FILESEP_UNIX);
    if ~isempty(HPC_QOS_NAME)
        path_out = [FILESEP_UNIX HPC_DRIVE_NAME FILESEP_UNIX HPC_QOS_NAME FILESEP_UNIX path_out{1}];
    else
        path_out = [FILESEP_UNIX HPC_DRIVE_NAME FILESEP_UNIX path_out{1}];
    end
%## if there are dos seps in path && there is a UNIX initiation sep that
%   got converted to a dos sep
elseif ~isempty(dos_sep) && any(dos_sep==1)
    path_out = strsplit(fPath,FILESEP_DOS);
    tmp = join(path_out(2:end),FILESEP_UNIX);
    if ~strcmp(path_out{2},HPC_DRIVE_NAME)
        if ~isempty(HPC_QOS_NAME)
            path_out = [ FILESEP_UNIX HPC_QOS_NAME FILESEP_UNIX tmp{1}];
        else
            path_out = [ FILESEP_UNIX tmp{1}];
        end
    else
        path_out = [FILESEP_UNIX tmp{1}];
    end
else
    path_out = fPath;
end
end

