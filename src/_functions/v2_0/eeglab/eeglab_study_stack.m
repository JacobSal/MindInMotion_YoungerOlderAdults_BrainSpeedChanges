function [ALL_STUDY,ALL_ALLEEG] = eeglab_study_stack(ALL_SUBJSTRUCT,trial_names,trial_iters,varargin)
%EEGLAB_STUDY_STACK Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## DEFINE DEFAULTS

%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALL_SUBJSTRUCT',@iscell);
addRequired(p,'trial_names',@iscell);
addRequired(p,'trial_iters',@isnumeric);
%## OPTIONAL
%## PARAMETER
parse(p,ALL_SUBJSTRUCT,trial_names,trial_iters,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
ALL_STUDY = cell(1,length(trial_names));
ALL_ALLEEG = cell(1,length(trial_names));
for study_i = 1:length(trial_names)    
    tt = trial_names{study_i};
    ti = trial_iters(study_i);
    tmp = ALL_SUBJSTRUCT{study_i};
    fprintf('==== Processing Condition %s ====\n',tt)
    if ~ispc
        studyDir = convertPath2UNIX(tmp(1).META.batch.study(ti).filepath,'dferris');
    else
        studyDir = convertPath2Drive(tmp(1).META.batch.study(ti).filepath,'M');
    end
    try
        studyName = tmp(1).META.batch.study(ti).filename;
    catch
        error('Error. No .study file available for %s',tt);
    end
    if ~ispc
        [ALL_STUDY{study_i}, ALL_ALLEEG{study_i}] = pop_loadstudy('filename',[studyName,'_UNIX','.study'],'filepath',studyDir);        
    else
        [ALL_STUDY{study_i}, ALL_ALLEEG{study_i}] = pop_loadstudy('filename',[studyName,'.study'],'filepath',studyDir); 
    end
end
%##
toc
end

