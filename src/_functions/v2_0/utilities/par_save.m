function [] = par_save(SAVEVAR,fPath,fName,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           SAVEVAR, variable to save.
%               STRUCT, CHAR, DICT, ARRAY you want to save.
%           fPath, CHAR
%               path to the folder where your file is held
%           fName, CHAR
%               file name & extension (e.g., 'INEEG.mat')
%       OPTIONAL:
%           fname_ext, CHAR
%               for automating the renaming of file names 
%       PARAMETER:
%   OUT:
%       NONE
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
%- fPath
errorMsg = 'Value ''fPath'' must be CHAR. The path must exist.'; 
fp_validFcn = @(x) assert(ischar(x) && exist(x,'dir'),errorMsg);
%- fname_ext
FNAME_EXT = '';
errorMsg = 'Value ''fname_ext'' must be CHAR. This value is appended to ''fName'' before the file declaration.'; 
fn_validFcn = @(x) assert(ischar(x),errorMsg);
p = inputParser;
%## REQUIRED
addRequired(p,'SAVEVAR')
addRequired(p,'fPath',fp_validFcn)
addRequired(p,'fName',@ischar)
%## OPTIONAL
addOptional(p,'fname_ext',FNAME_EXT,fn_validFcn);
%## PARAMETER
parse(p, SAVEVAR, fPath, fName, varargin{:});
%## SET DEFAULTS
fname_ext = p.Results.fname_ext;
%% ===================================================================== %%
if ~ispc
    % convert path to os
    pathIn = convertPath2UNIX(fPath);
else
    %- convert path to os
    pathIn = convertPath2Drive(fPath);
end
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
%- save
%     val = GetSize(SAVEVAR);
s = whos('SAVEVAR');
fprintf(1,'%s is %0.2g bytes\n',fName,s.bytes);
% fprintf('\nSaving %s to\n%s\n',fName,savePath);
if s.bytes >= 2e9
    fprintf('\nSaving %s using ''v7.3'' to\n%s\n',fName,savePath);
    save(savePath,'SAVEVAR','-v7.3');
else
    fprintf('\nSaving %s using ''v6'' to\n%s\n',fName,savePath);
    save(savePath,'SAVEVAR','-v6');
end
end

