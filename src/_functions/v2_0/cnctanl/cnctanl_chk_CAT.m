function [chk] = cnctanl_chk_CAT(EEG,cond_chk,fpath_chk,varargin)
%CNCTANL_CHK_CAT Summary of this function goes here
%   This is a CUSTOM function
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

%## DEFINE DEFAULTS
%*
EPOCH_FNAME_REGEXP = '%s_%s_EPOCH_TMPEEG.set';
%*
fName = strsplit(EPOCH_FNAME_REGEXP,'.'); fName = [fName{1} '.mat'];
fName = strsplit(fName,'.');
fName{1} = [fName{1} '_BootStrap'];
fName = strjoin(fName,'.');
BOOTSTRAP_FNAME_REGEXP = fName;
%*
fName = strsplit(EPOCH_FNAME_REGEXP,'.'); fName = [fName{1} '.mat'];
fName = strsplit(fName,'.');
fName{1} = [fName{1} '_PhaseRnd'];
fName = strjoin(fName,'.');
PHASERND_FNAME_REGEXP = fName;
%-
CHK_BOOTSTRAP = 0;
CHK_PHASE_RND = 0;
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'cond_chk',@iscell);
addRequired(p,'fpath_chk',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'CHK_BOOTSTRAP',CHK_BOOTSTRAP,@islogical);
addParameter(p,'CHK_PHASE_RND',CHK_PHASE_RND,@islogical);
%## PARSE
parse(p,EEG,cond_chk,fpath_chk,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER3
CHK_BOOTSTRAP = p.Results.CHK_BOOTSTRAP;
CHK_PHASE_RND = p.Results.CHK_PHASE_RND;
%- PERMS

%% ===================================================================== %%
fPath = [fpath_chk filesep sprintf('cond_%s',cond_chk)];
fName = sprintf(EPOCH_FNAME_REGEXP,EEG.subject,cond_chk);
if  ~exist([fPath filesep fName],'file')
    chk(cond_i) = 1;
end
if CHK_BOOTSTRAP
    fName = sprintf(BOOTSTRAP_FNAME_REGEXP,EEG.subject,cond_chk);
    if  ~exist([fPath filesep fName],'file')
        chk(cond_i) = 1;
    end
end
if CHK_PHASE_RND
    fName = sprintf(PHASERND_FNAME_REGEXP,EEG.subject,cond_chk);
    if  ~exist([fPath filesep fName],'file')
        chk(cond_i) = 1;
    end
end

end

