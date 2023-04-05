function [EEG] = cnctanl_loadCAT(EEG,SUBJSTRUCT,processNum,subjectNum,loadVar,varargin)
%CNCTANL_LOADCAT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {true};
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'subjectNum',@isnumeric);
addRequired(p,'loadVar',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,EEG,SUBJSTRUCT,processNum,subjectNum,loadVar,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%## PARAMS
ASSIGN_BOOTSTRAP_MEAN = false;
%## MAKE DIRS
%## DETECT OS
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
fPath = SUBJSTRUCT(subjectNum).pathsEEG.cnctAnl(processNum).filepath;
fName = SUBJSTRUCT(subjectNum).pathsEEG.cnctAnl(processNum).filename;
if DO_UNIX
    fPath = convertPath2UNIX(fPath,'dferris');
else
    fPath = convertPath2Drive(fPath,'M');
end

switch loadVar
    case 'BootStrap'
        %- Bootstrap Data Handler
        fprintf('\n==== LOADING BOOTSTRAP MEASURES ====\n')        
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file')
            error('Run GLOBAL_BATCH to generate bootstrap permutation test values');
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
            error('Run GLOBAL_BATCH to generate randomized permutation test values');
        else
            EEG.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
        end
        fprintf('done.\n')
        %- Nonzero Statistics Data Handler
        fprintf('\n==== LOADING NONZERO STATISTICS ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file')
            error('Run GLOBAL_BATCH to generate nonzero connectivity statistic');
        else
            EEG.CAT.Stats = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
        end
        fprintf('done.\n')
    otherwise
        error('''%s'' load variable is not available',loadVar);
end
end

