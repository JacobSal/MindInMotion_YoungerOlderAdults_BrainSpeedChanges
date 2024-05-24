function [msg] = transfer_folder(from_folder,to_folder,varargin)
%TRANSFERFOLDER Summary of this function goes here
%   Detailed explanation goes here
%   
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%## TIME
tic
%## DEFINE DEFAULTS
REGEXP = '';
msg = '';
p = inputParser;
%## REQUIRED
addRequired(p,'from_folder',@ischar)
addRequired(p,'to_folder',@ischar)
%## OPTIONAL
addOptional(p,'REGEXP',REGEXP,@ischar);
%## PARAMETER
parse(p, from_folder, to_folder, varargin{:});

%## SET DEFAULTS
REGEXP = p.Results.REGEXP;
% FNAME = 'tmpEEG';
%## SAVE FOLDER HANDLER (HIGHEST ORDER FUNCTION ONLY)
%% ===================================================================== %%
%## save to new path
fprintf('Transfering files. This might take awhile...\n');
% if dirOut doesn't exist, create it.
if ~exist(to_folder,'dir')
    fprintf('Making directory: %s\n',to_folder);
    mkdir(to_folder);
end
% dont copy if float/set files already exist
dirFrom = dir([from_folder filesep REGEXP]);
if ~isempty(dirFrom)
    tmpF = [dirFrom(1).folder];
    try
        copyfile(tmpF,to_folder);
    catch
        fprintf('Not copied: %s\n',tmpF);
        return;
    end
else
    error('Folder or regexp ''%s'' does not exist!',[from_folder filesep REGEXP]);
end
fprintf('done!\n')
%## TIME
toc
end

