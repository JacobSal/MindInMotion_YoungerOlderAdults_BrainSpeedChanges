function [ALLEEG] = nj_parse_trials(EEG,trial_names,epoch_limits,varargin)
%NJ_EPOCH_DATA Summary of this function goes here
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
%- Keeping as HARD defines for now (03/08/2023)
PARSE_TYPE = '';
PERCENT_OVERLAP = 0;
EVENT_FIELD_CONDITION = 'cond';
EVENT_FIELD_TRIAL = 'type';
TRIAL_BEG_CHAR = 'TrialStart';
TRIAL_END_CHAR = 'TrialEnd';
APPROX_TRIAL_LENGTH = 3*60; % seconds
% SLIDING_WINDOW_ONLY = false;
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'trial_names',@iscell);
addRequired(p,'epoch_limits',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,@ischar);
addParameter(p,'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,@ischar);
addParameter(p,'EVENT_FIELD_TRIAL',EVENT_FIELD_TRIAL,@ischar);
addParameter(p,'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,@ischar);
addParameter(p,'TRIAL_END_CHAR',TRIAL_END_CHAR,@ischar);
addParameter(p,'APPROX_TRIAL_LENGTH',APPROX_TRIAL_LENGTH,@isnumeric);
% addParameter(p,'SLIDING_WINDOW_ONLY',SLIDING_WINDOW_ONLY,@islogical);
parse(p,EEG,trial_names,epoch_limits,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
EVENT_FIELD_CONDITION = p.Results.EVENT_FIELD_CONDITION;
EVENT_FIELD_TRIAL = p.Results.EVENT_FIELD_TRIAL;
TRIAL_BEG_CHAR = p.Results.TRIAL_BEG_CHAR;
TRIAL_END_CHAR = p.Results.TRIAL_END_CHAR;
APPROX_TRIAL_LENGTH = p.Results.APPROX_TRIAL_LENGTH;
% SLIDING_WINDOW_ONLY = p.Results.SLIDING_WINDOW_ONLY;
%- PERMS

%% ===================================================================== %%
%## STEP 2) EPOCHING
%* empty ALLEEG structure for repopulating
ALLEEG = cell(1,length(trial_names)); 
fprintf(1,'\n==== EPOCHING  ====\n');
for trial_i = 1:length(trial_names)
    fprintf(1,'\n==== %s: Processing trial %s ====\n',EEG.subject,trial_names{trial_i});
    %## EPOCH SUBJ DATA
    %- override event structure to be compatible with sliding window
    tmp_all = strcmp({EEG.event.(EVENT_FIELD_CONDITION)},trial_names{trial_i});
    trial_start = find(tmp_all,1,'first');
    trial_end = find(tmp_all,1,'last');
    EEG.event(trial_start).type = TRIAL_BEG_CHAR;
    EEG.event(trial_start).cond = trial_names{trial_i};
    EEG.event(trial_end).type = TRIAL_END_CHAR;
    EEG.event(trial_end).cond = trial_names{trial_i};
    %-
    [TMP_EEG] = eeglab_sliding_window_epoch(EEG,trial_names{trial_i},(epoch_limits(2)-epoch_limits(1)),...
                'PERCENT_OVERLAP',PERCENT_OVERLAP,...
                'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,...
                'EVENT_FIELD_TRIAL',EVENT_FIELD_TRIAL,...
                'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,...
                'TRIAL_END_CHAR',TRIAL_END_CHAR,...
                'APPROX_TRIAL_LENGTH',APPROX_TRIAL_LENGTH);
    TMP_EEG.timewarp = struct([]);
    TMP_EEG.etc.epoch.parse_var = 'sliding_window';

    %## Update EEG structure
    TMP_EEG.condition = trial_names{trial_i};
    %- assign a new name
    TMP_EEG.name            = sprintf('Epoch_%s-%s',PARSE_TYPE,trial_names{trial_i});
    TMP_EEG.filename        = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,trial_names{trial_i});
    TMP_EEG.etc.epoch.epoch_limits   = epoch_limits;    
    %- check to make sure a number isn't the first character
    chk = regexp(trial_names{trial_i},'\d');
    if any(chk)
        trial_names{trial_i} = sprintf('x%s',trial_names{trial_i});
    end
    TMP_EEG.etc.epoch.condition = trial_names{trial_i};
    %## STORE IN ALLEEG
    ALLEEG{trial_i} = TMP_EEG;
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
end

