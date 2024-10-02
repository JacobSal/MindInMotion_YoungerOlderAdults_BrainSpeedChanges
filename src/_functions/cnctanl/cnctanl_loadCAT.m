function [EEG] = cnctanl_loadCAT(EEG,loadVar,varargin)
%CNCTANL_LOADCAT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
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
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'loadVar',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,EEG,loadVar,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS
ASSIGN_BOOTSTRAP_MEAN = false;

%% ===================================================================== %%
if ispc
    fPath = convertPath2Drive(EEG.filepath);
else
    fPath = convertPath2UNIX(EEG.filepath);
end
fName = EEG.filename;
%## Load connectivity statisticsds
switch loadVar
    case 'BootStrap'
        %- Bootstrap Data Handler
        fprintf('\n==== LOADING BOOTSTRAP MEASURES ====\n')        
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file')
            error('%s does not exist.\nRun GLOBAL_BATCH to generate bootstrap permutation test values',[fPath filesep [chk{1}, '_BootStrap.mat']]);
        else
            EEG.CAT.PConn = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
        end
        %- assign mean of bootstrap as Conn value
        if ASSIGN_BOOTSTRAP_MEAN
            EEG.CAT.Conn = stat_getDistribMean(EEG.CAT.PConn);
        end
    case 'NonzeroTest'
        %- Phase Randomization Permutation Test Data Handler
        fprintf('\n==== LOADING PHASE RANDOMIZED CONNECTIVITY MEASURES ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file') 
            error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
        else
            EEG.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
        end
        fprintf('done.\n')
        %- Nonzero Statistics Data Handler
        fprintf('\n==== LOADING NONZERO STATISTICS ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file')
            error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_NonZero.mat']]);
        else
            EEG.CAT.Stats = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
        end
        fprintf('done.\n')
    otherwise
        error('''%s'' load variable is not available',loadVar);
end
end

