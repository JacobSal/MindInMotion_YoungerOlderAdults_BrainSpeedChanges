function [path_out] = convertPath2Drive(fPath,varargin)
%prepHPG4AMICA
%#ok<*STREMP>
%   IN: 
%       fPath, CHAR
%           
%       drive_letter, CHAR
%           
%       HPC_DRIVE_NAME, CHAR
%           
%       HPC_QOS_NAME, CHAR
%           
%   OUT: 
%       path_out, CHAR
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
DRIVE_LETTER = 'M';
errorMsg = 'Value must be CHAR. This value is the letter/name assigned to your storage drive (e.g., ''M'').'; 
dl_validFcn = @(x) assert(ischar(x),errorMsg);
HPC_DRIVE_NAME = 'blue';
errorMsg = 'Value must be CHAR. This value is the name of the STORAGE FOLDER on HiperGator.'; 
hdn_validFcn = @(x) assert(ischar(x),errorMsg);
HPC_QOS_NAME = 'dferris';
errorMsg = 'Value must be CHAR. This value is the name of the USER FOLDER on HiperGator.'; 
hqn_validFcn = @(x) assert(ischar(x),errorMsg);
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'fPath',@ischar);
%## OPTIONAL
addOptional(p,'DRIVE_LETTER',DRIVE_LETTER,dl_validFcn);
addOptional(p,'HPC_DRIVE_NAME',HPC_DRIVE_NAME,hdn_validFcn);
addOptional(p,'HPC_QOS_NAME',HPC_QOS_NAME,hqn_validFcn);
%## PARAMETER
%## SET DEFAULTS
parse(p,fPath,varargin{:});
DRIVE_LETTER = p.Results.DRIVE_LETTER;
HPC_DRIVE_NAME = p.Results.HPC_DRIVE_NAME;
HPC_QOS_NAME = p.Results.HPC_QOS_NAME;
%## SAVE FOLDER HANDLER (HIGHEST ORDER FUNCTION ONLY)
%% ===================================================================== %%
FILESEP_DOS = '\'; %DOS override
FILESEP_UNIX = '/';
% path_out = [];
if contains(fPath,FILESEP_UNIX) && contains(fPath,FILESEP_DOS)
    fPath = fullfile(fPath);
end
unix_sep = strfind(fPath,FILESEP_UNIX);
dos_sep = strfind(fPath,FILESEP_DOS);
%## if there are unix seps in fPath && the UNIX initiation sep exists
if ~isempty(unix_sep) && any(unix_sep==1)
    %- split path using UNIX filesep
    path_out = strsplit(fPath,FILESEP_UNIX);
    %- find indices where the HiperGator drive and account names are
    hdn_idx = cellfun(@(x) strcmp(x,HPC_DRIVE_NAME), path_out);
    hqn_idx = cellfun(@(x) strcmp(x,HPC_QOS_NAME), path_out);
    dl_idx = cellfun(@(x) strcmp(x,[DRIVE_LETTER ':']), path_out);
    %- create new path
    idx2keep = ~(hdn_idx+hqn_idx+dl_idx);
    path_out = path_out(idx2keep);
    path_out = join(path_out(2:end),FILESEP_DOS);
    path_out = [DRIVE_LETTER ':' FILESEP_DOS path_out{1}];
%## if there are dos seps in path && the UNIX initiation sep got converted to
%  dos sep
elseif ~isempty(dos_sep) && any(dos_sep==1)
    %- split path using UNIX filesep
    path_out = strsplit(fPath,FILESEP_DOS);
    path_out = path_out(2:end);
    %- find indices where the
    hdn_idx = cellfun(@(x) strcmp(x,HPC_DRIVE_NAME), path_out);
    hqn_idx = cellfun(@(x) strcmp(x,HPC_QOS_NAME), path_out);
    dl_idx = cellfun(@(x) strcmp(x,[DRIVE_LETTER ':']), path_out);
    %- create new path
    idx2keep = ~(hdn_idx+hqn_idx+dl_idx);
    path_out = path_out(idx2keep);
    path_out = join(path_out,FILESEP_DOS);
    path_out = [DRIVE_LETTER ':' FILESEP_DOS path_out{1}];
%## if UNIX path got inputted without an initiation sep
elseif ~isempty(unix_sep) && ~any(unix_sep==1)
%- split path using UNIX filesep
    path_out = strsplit(fPath,FILESEP_UNIX);
    %- find indices where the
    hdn_idx = cellfun(@(x) strcmp(x,HPC_DRIVE_NAME), path_out);
    hqn_idx = cellfun(@(x) strcmp(x,HPC_QOS_NAME), path_out);
    dl_idx = cellfun(@(x) strcmp(x,[DRIVE_LETTER ':']), path_out);
    %- create new path
    idx2keep = ~(hdn_idx+hqn_idx+dl_idx);
    path_out = path_out(idx2keep);
    path_out = join(path_out,FILESEP_DOS);
    path_out = [DRIVE_LETTER ':' FILESEP_DOS path_out{1}];
%## if UNIX path got inputted without an initiation sep
% elseif ~isempty(dos_sep) && ~any(dos_sep==1) 
%     %- split path using UNIX filesep
%     path_out = strsplit(fPath,FILESEP_DOS);
%     %- find indices where the
%     hdn_idx = cellfun(@(x) strcmp(x,HPC_DRIVE_NAME), path_out);
%     hqn_idx = cellfun(@(x) strcmp(x,HPC_QOS_NAME), path_out);
%     dl_idx = cellfun(@(x) strcmp(x,[DRIVE_LETTER ':']), path_out);
%     %- create new path
%     idx2keep = ~(hdn_idx+hqn_idx+dl_idx);
%     path_out = path_out(idx2keep);
%     path_out = join(path_out,FILESEP_DOS);
%     path_out = [DRIVE_LETTER ':' FILESEP_DOS path_out{1}];
else
    path_out = fPath;
end
end

