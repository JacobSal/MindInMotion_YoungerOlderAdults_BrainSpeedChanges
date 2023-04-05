function [SUBJ] = subj_i_genConnMeas(PATHS,SUBJ,processNum,studyName,studyDir,connMethods,brainChars,brainCoords,varargin)
%GENCONNMEASSUBJS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {10};
p = inputParser; 
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJ',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
addRequired(p,'connMethods',@iscell);
addRequired(p,'brainChars',@iscell);
addRequired(p,'brainCoords',@iscell);
%## OPTIONAL
%## PARAMETER
%## PARSE
parse(p,PATHS,SUBJ,processNum,studyName,studyDir,connMethods,brainChars,brainCoords,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%## PARAMS 
REGEN_FILES = true;
DO_BOOTSTRAP = false;
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
ASSIGN_BOOTSTRAP_MEAN = false; %#ok<NASGU>
ANALYSIS_NAME = SUBJ.META.epoched.trialType;
CONN_METHODS = connMethods;
SUBJ.cnctAnl.connMethods = connMethods;
%## MAKE DIRS
mkdir(studyDir);
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

%## MAIN CONNECTIVITY PROCESSING (see. indvCnctAnl.m)
fprintf('\n==== CALCULATING CONNECTIVITY MEASURES ====\n')
%## PARAMS
components = SUBJ.cnctAnl.(ANALYSIS_NAME).components;
fPath_e = SUBJ.pathsEEG.epoched(processNum).filepath;
fName_e = SUBJ.pathsEEG.epoched(processNum).filename;
fPath = SUBJ.pathsEEG.cnctAnl(processNum).filepath;
fName = SUBJ.pathsEEG.cnctAnl(processNum).filename;
if DO_UNIX
    fPath = convertPath2UNIX(fPath,'dferris');
else
    fPath = convertPath2Drive(fPath,'M');
end
subjStr = SUBJ.subjStr;
%- Prints
fprintf('%s) Processing componets:\n',subjStr)
fprintf('%i,',components'); fprintf('\n');
%## GENERATE CONNECTIVITY MEASURES
%- check to see if connectivity is already calculated
if ~exist([fPath filesep fName],'file') || REGEN_FILES
    %- Load EEG Epoch
    EEG = pop_loadset('filename',fName_e,'filepath', fPath_e);
    %- Calculate Connectivity
    tmpEEG = cnctanl_connMeas(EEG,PATHS,components,CONN_METHODS,fPath,...
        'WindowLengthSec',WINDOW_LENGTH,'WindowStepSizeSec', WINDOW_STEP_SIZE,...
        'NewSamplingRate',NEW_SAMPLE_RATE);
else
    tmpEEG = par_load(fPath,fName);
end
%% (BOOTSTRAPPING) GROUP STATISTICS 
% (09/22/2022) JS, Might want to try and speed up bootstrap by
% addapting stat_surrogateGen.m to use parfor for bootstrapping... If
% possible? doesn't seem built well in the first place, so maybe?
if DO_BOOTSTRAP
    fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')         %#ok<UNRCH>
    chk = strsplit(fName,'.');
    %- check to see if BootStrap is already calculated
    if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') || REGEN_FILES
        tmpEEG = cnctanl_groupStats(tmpEEG,PATHS,'BootStrap');
        %- save BootStrap distribution
        tmp = tmpEEG.CAT.PConn;
        par_save(tmp,fPath,fName,[],'_BootStrap');
    else
        tmpEEG.CAT.PConn = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
    end
    %- assign mean of bootstrap as Conn value
    if ASSIGN_BOOTSTRAP_MEAN
        tmpEEG.CAT.Conn = stat_getDistribMean(tmpEEG.CAT.PConn);
    end
    %- clear bootstrap calculation
    tmpEEG.CAT.PConn = [];       
end
%- save EEG 
par_save(tmpEEG,fPath,fName,[]);
end

