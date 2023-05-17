function [ALLEEG] = as_parse_trials(EEG,cond_names,event_names,epoch_limits,varargin)
%MIM_PARSE_TRIALS Summary of this function goes here
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
SLIDING_WINDOW_ONLY = false;
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'cond_names',@iscell);
addRequired(p,'event_names',@iscell);
addRequired(p,'epoch_limits',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,@ischar);
addParameter(p,'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,@ischar);
addParameter(p,'EVENT_FIELD_TRIAL',EVENT_FIELD_TRIAL,@ischar);
addParameter(p,'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,@ischar);
addParameter(p,'TRIAL_END_CHAR',TRIAL_END_CHAR,@ischar);
addParameter(p,'APPROX_TRIAL_LENGTH',APPROX_TRIAL_LENGTH,@isnumeric);
addParameter(p,'SLIDING_WINDOW_ONLY',SLIDING_WINDOW_ONLY,@islogical);
parse(p,EEG,cond_names,event_names,epoch_limits,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
EVENT_FIELD_CONDITION = p.Results.EVENT_FIELD_CONDITION;
EVENT_FIELD_TRIAL = p.Results.EVENT_FIELD_TRIAL;
TRIAL_BEG_CHAR = p.Results.TRIAL_BEG_CHAR;
TRIAL_END_CHAR = p.Results.TRIAL_END_CHAR;
APPROX_TRIAL_LENGTH = p.Results.APPROX_TRIAL_LENGTH;
SLIDING_WINDOW_ONLY = p.Results.SLIDING_WINDOW_ONLY;
%- PERMS

%% ===================================================================== %%
%- empty ALLEEG structure for repopulating
ALLEEG = cell(1,length(cond_names)); 
for cond_i = 1:length(cond_names)
    %- epoch depending on number of available trials
    for event_i = 1:length(event_names)       
        fprintf(1,'\n==== Processing ''%s'' ====\n',cond_names{cond_i});    
        
        if SLIDING_WINDOW_ONLY
            %- override event structure to be compatible with sliding window
            %- SLIDING WINDOW: MIM specific function
            %- NOTE: (03/29/2023) JS, may still be buggy for amanda's data
            %will need to make changes to this function to adapt to
            %amanda's event data.
            [TMP_EEG] = eeglab_sliding_window_epoch(EEG,cond_names{trial_i},(epoch_limits(2)-epoch_limits(1)),...
                    'PERCENT_OVERLAP',PERCENT_OVERLAP,...
                    'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,...
                    'EVENT_FIELD_TRIAL',EVENT_FIELD_TRIAL,...
                    'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,...
                    'TRIAL_END_CHAR',TRIAL_END_CHAR,...
                    'APPROX_TRIAL_LENGTH',APPROX_TRIAL_LENGTH);
            TMP_EEG.timewarp = struct([]);
            TMP_EEG.etc.epoch.parse_var = 'sliding_window';
        else
            %- AS processing
            % (02/14/2023) JS, May need to make this better for amanda's data
            % seems to be dropping some epochs that shouldn't be.
            tmp = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},cond_names{cond_i});
            valid_points = [find(tmp,1,'first'),find(tmp,1,'last')];
            TMP_EEG = pop_epoch(EEG,event_names(event_i),epoch_limits,...
                'eventindices',valid_points(1):valid_points(2),'epochinfo','yes');
        end
        %## Update EEG structure
        TMP_EEG.condition = cond_names{cond_i};
        %- assign a new name
        TMP_EEG.filename        = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,cond_names{cond_i},event_names{event_i});
        TMP_EEG.etc.epoch.epoch_limits   = epoch_limits;    
        %- check to make sure a number isn't the first character
        chk = regexp(cond_names{cond_i},'\d');
        if any(chk)
            cond_names{cond_i} = sprintf('x%s',cond_names{cond_i});
        end
        TMP_EEG.etc.epoch.condition = [cond_names{cond_i} '_' event_names{event_i}]; % trial_names{trial_i};
        %## STORE IN ALLEEG
        ALLEEG{cond_i} = TMP_EEG;
    end
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);

end