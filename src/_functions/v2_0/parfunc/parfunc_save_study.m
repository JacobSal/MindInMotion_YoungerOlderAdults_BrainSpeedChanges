function [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,varargin)
%EEGLAB_ROBUST_SAVE Summary of this function goes here
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
%% Define Defaults
%*
RESAVE_DATASETS = 'off';
%* remove file extension from file name is it exists
STUDY_FNAME = STUDY.filename;
errorMsg = 'Value must be a CHAR. STUDY filename.'; 
sfn_validFcn = @(x) assert(ischar(x),errorMsg);
%*
STUDY_FPATH = STUDY.filepath;
errorMsg = 'Value must be a CHAR. STUDY filepath or the path to which the .study file will be saved to.'; 
sfp_validFcn = @(x) assert(ischar(x),errorMsg);
%*
STUDY_COND  = [];
errorMsg = 'Value must be a CHAR. STUDY condition.'; 
sc_validFcn = @(x) assert(ischar(x)||isempty(x),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
addOptional(p,'STUDY_FNAME',STUDY_FNAME,sfn_validFcn);
addOptional(p,'STUDY_FPATH',STUDY_FPATH,sfp_validFcn);
%## PARAMETER
addParameter(p,'STUDY_COND',STUDY_COND,sc_validFcn);
addParameter(p,'RESAVE_DATASETS',RESAVE_DATASETS,@ischar);
%- parse
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
STUDY_FNAME = p.Results.STUDY_FNAME;
STUDY_FPATH = p.Results.STUDY_FPATH;
STUDY_COND = p.Results.STUDY_COND;
RESAVE_DATASETS = p.Results.RESAVE_DATASETS;
%- OPTIONAL
%- PARAMETER
%% ===================================================================== %%
%## ROBUST SAVE STUDY FOR SLURM/UNIX/MAC USE
STUDY_FNAME = chk_study_fName(STUDY_FNAME);
fprintf(1,'\n==== SAVING STUDY ''%s'' ====\n',STUDY_FNAME);
usave = sprintf('%s_UNIX.study',STUDY_FNAME);
dsave = sprintf('%s.study',STUDY_FNAME);
if ~ispc
    %- make path conversions based on pathing
    for subj_i = 1:length(STUDY.datasetinfo)
        STUDY.datasetinfo(subj_i).filepath = ALLEEG(subj_i).filepath;
        STUDY.datasetinfo(subj_i).filename = ALLEEG(subj_i).filename;
        if ~isempty(STUDY_COND)
            STUDY.datasetinfo(subj_i).condition = STUDY_COND;
        end
    end
    %- Save STUDY for UNIX
    [STUDY,ALLEEG] = pop_savestudy( STUDY, ALLEEG,...
                        'filename',usave,'filepath',STUDY_FPATH,'resavedatasets',RESAVE_DATASETS);
    tmp = STUDY;
    %- Change datasetinfo filepaths and resave for accessing via DRIVE
    for subj_i = 1:length(tmp.datasetinfo)
        tmp.datasetinfo(subj_i).filepath = convertPath2Drive(ALLEEG(subj_i).filepath);
        tmp.datasetinfo(subj_i).filename = ALLEEG(subj_i).filename;
        if ~isempty(STUDY_COND)
            tmp.datasetinfo(subj_i).condition = STUDY_COND;
        end
    end
    [~,~] = pop_savestudy( tmp, ALLEEG,...
                        'filename',dsave,'filepath',STUDY_FPATH,'resavedatasets',RESAVE_DATASETS);
else
    %- make path conversions based on pathing
    for subj_i = 1:length(STUDY.datasetinfo)
        STUDY.datasetinfo(subj_i).filepath = ALLEEG(subj_i).filepath;
        STUDY.datasetinfo(subj_i).filename = ALLEEG(subj_i).filename;
        if ~isempty(STUDY_COND)
            STUDY.datasetinfo(subj_i).condition = STUDY_COND;
        end
    end
    [STUDY,ALLEEG] = pop_savestudy( STUDY, ALLEEG,...
                        'filename',dsave,'filepath',STUDY_FPATH,'resavedatasets',RESAVE_DATASETS);
    %- save STUDY for DRIVE
    tmp = STUDY;
    for subj_i = 1:length(tmp.datasetinfo)
        tmp.datasetinfo(subj_i).filepath = convertPath2UNIX(ALLEEG(subj_i).filepath);
        tmp.datasetinfo(subj_i).filename = ALLEEG(subj_i).filename;
        if ~isempty(STUDY_COND)
            tmp.datasetinfo(subj_i).condition = STUDY_COND;
        end
    end
    %- Change datasetinfo filepaths and resave for accessing via UNIX
    [~,~] = pop_savestudy( tmp, ALLEEG,...
                        'filename',usave,'filepath',STUDY_FPATH,'resavedatasets',RESAVE_DATASETS);
end
fprintf(1,'\nDONE.\n');

%% ===================================================================== %%
%## SUBFUNCTIONS
function [ftmp] = chk_study_fName(study_fName)
    tmp_fName = strsplit(study_fName,'.');
    if length(tmp_fName)>1
        ftmp = tmp_fName{1};
    else
        ftmp = study_fName;
    end
    if regexp(ftmp,'_UNIX','once')
        ftmp = ftmp(1:(regexp(ftmp,'_UNIX','once'))-1);
    end
end
end
