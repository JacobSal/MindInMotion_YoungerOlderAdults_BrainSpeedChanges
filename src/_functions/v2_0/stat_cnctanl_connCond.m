function [ALLEEG] = stat_cnctanl_connCond(PATHS,ALLEEG,studyDir,studyName,connMethods,components,clusters,varargin)
%STAT_CNCTANL_CONNCOND Summary of this function goes here
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
addRequired(p,'PATHS',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'studyDir',@ischar);
addRequired(p,'studyName',@ischar)
addRequired(p,'connMethods',@iscell);
addRequired(p,'components',@iscell);
addRequired(p,'clusters',@iscell);
%## OPTIONAL
%## PARAMETER
parse(p,PATHS,ALLEEG,studyDir,studyName,connMethods,components,clusters,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%## PARAMS
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
if DO_UNIX
    studyDir = convertPath2UNIX(studyDir,'dferris');
else
    studyDir = convertPath2Drive(studyDir,'M');
end
%% ===================================================================== %%
%- Check ALLEEG dataset for right format
if (length(ALLEEG) ~= 2)
    error(['ALLEEG should be length 2,',...
        'where each EEG dataset contains a single condition.\n']);
end
%- assign saveDir
saveDir = [studyDir filesep sprintf('%s-%s',ALLEEG(1).subject,ALLEEG(2).subject) filesep];
if ~exist(saveDir,'dir')
    fprintf(1,'Making Directory: %s',saveDir);
    mkdir(saveDir)
end
%## LOAD CONNECTIVITY
[ALLEEG,~] = groupStatsCnctAnl(ALLEEG,PATHS,'BtwnCond');
end