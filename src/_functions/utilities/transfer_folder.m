function [msg] = transfer_folder(folder_from,folder_to,file_regexp)
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
msg = '';
p = inputParser;
%## REQUIRED
addRequired(p,'folder_from',@ischar)
addRequired(p,'folder_to',@ischar)
addRequired(p,'file_regexp',@ischar)
%## PARAMETER
parse(p, folder_from, folder_to, file_regexp);
%% ===================================================================== %%
%## save to new path
fprintf('Transfering files. This might take awhile...\n');
% if dirOut doesn't exist, create it.
if ~exist(folder_to,'dir')
    fprintf('Making directory: %s\n',folder_to);
    mkdir(folder_to);
end
% dont copy if float/set files already exist
dir_from = dir([folder_from filesep file_regexp]);
fprintf('Transfering %0.2g GB of data...\n',sum([dir_from.bytes])/(1e9))
if ~isempty(dir_from)
    for fi = 1:length(dir_from)
        tmpf = [dir_from(fi).folder filesep dir_from(fi).name];
        try
            copyfile(tmpf,folder_to);
        catch e
            fprintf('Not copied: %s\n%s\n\n',tmpf,getReport(e));
            return;
        end
    end
else
    error('Folder or regexp ''%s'' does not exist!',[folder_from filesep file_regexp]);
end
fprintf('transfer_folder.m done.\n')
%## TIME
toc
end

