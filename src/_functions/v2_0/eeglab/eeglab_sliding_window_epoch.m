function [EEG] = eeglab_sliding_window_epoch(EEG,cond_char,window_len,varargin)
%   This is a CUSTOM function
%   IN: 
%       EEG, STRUCT
%           
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
PERCENT_OVERLAP = 0;
EVENT_FIELD_CONDITION = 'cond';
EVENT_FIELD_TRIAL = 'type';
TRIAL_BEG_CHAR = 'TrialStart';
TRIAL_END_CHAR = 'TrialEnd';
APPROX_TRIAL_LENGTH = 3*60; % seconds
%-
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'cond_char',@ischar);
addRequired(p,'window_len',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,@isnumeric);
addParameter(p,'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,@ischar);
addParameter(p,'EVENT_FIELD_TRIAL',EVENT_FIELD_TRIAL,@ischar);
addParameter(p,'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,@ischar);
addParameter(p,'TRIAL_END_CHAR',TRIAL_END_CHAR,@ischar);
addParameter(p,'APPROX_TRIAL_LENGTH',APPROX_TRIAL_LENGTH,@isnumeric);
parse(p,EEG,cond_char,window_len,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
EVENT_FIELD_CONDITION = p.Results.EVENT_FIELD_CONDITION;
EVENT_FIELD_TRIAL = p.Results.EVENT_FIELD_TRIAL;
TRIAL_BEG_CHAR = p.Results.TRIAL_BEG_CHAR;
TRIAL_END_CHAR = p.Results.TRIAL_END_CHAR;
APPROX_TRIAL_LENGTH = p.Results.APPROX_TRIAL_LENGTH;

%% ================================================================= %%
%- find conditions that match input string
tmp_all = strcmp({EEG.event.(EVENT_FIELD_CONDITION)},cond_char);
tmp_start = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},TRIAL_BEG_CHAR);
tmp_end = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},TRIAL_END_CHAR);
valid_idxs = find(tmp_all & (tmp_start | tmp_end));
if rem(length(valid_idxs),2) ~= 0
        fprintf('Using all events for condition ''%s'' as 1 trial',cond_char);
        trial_start = find(tmp_all,1,'first');
        trial_end = find(tmp_all,1,'last');
        EEG.event(trial_start).type = 'tmp_start';
        EEG.event(trial_start).cond = cond_char;
        EEG.event(trial_end).type = 'tmp_end';
        EEG.event(trial_end).cond = cond_char;
        tmp_all = strcmp({EEG.event.(EVENT_FIELD_CONDITION)},cond_char);
        tmp_start = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},'tmp_start');
        tmp_end = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},'tmp_end');
        valid_idxs = find(tmp_all & (tmp_start | tmp_end));
    if rem(length(valid_idxs),2) ~= 0
        error(['Not all trial seperators are available in the EEG.event field for subject %s...\n',...
               'check EEG.event indices: [%i,%i]\n'],EEG.subject,find(tmp_all,1,'first'),find(tmp_all,1,'last'));
    end
end
spc = (abs(window_len)/2)*(1-PERCENT_OVERLAP);
%- initiate loop
trial_cnt = 1;
tmp_event = EEG.event;
tmp_trials = cell(1,length(valid_idxs)/2);
%## TRIAL APPENDING LOOP
for i = 1:2:length(valid_idxs)
    % split current event structure to put in new events
    % define constants in event struct  
    lat_1   = tmp_event(valid_idxs(i)).latency+EEG.srate*spc;
    lat_2   = lat_1+EEG.srate*APPROX_TRIAL_LENGTH;
    intervals = (lat_1:EEG.srate*spc*2:lat_2);
    events_out = cell(1,length(intervals));
    code_char = sprintf('trial_%i_cond_%s',trial_cnt,cond_char);
    dt = tmp_event(valid_idxs(i)).datetime;
    %- beginning boundary event
    events_out{1} = create_event_entry(intervals(1),...
                    window_len/2,'boundary',[],[],[]);
    %- inbetween events
    for j = 2:length(intervals)
        chk_intv = intervals(j) < (tmp_event(valid_idxs(i+1)).latency - EEG.srate*spc);
        if j == 1
            event_type = TRIAL_BEG_CHAR;
        elseif j > 1
            event_type = cond_char;
        end
        if chk_intv
            events_out{j} = create_event_entry(intervals(j),...
                    1,event_type,code_char,dt,cond_char); 
        end
    end
    %- trial end event
    events_out{j} = create_event_entry(intervals(j),...
                    1,TRIAL_END_CHAR,code_char,dt,cond_char); 
    %- ending boundary event
    events_out{j+1} = create_event_entry(intervals(j)+window_len/2,...
                    window_len/2,'boundary',[],[],[]);
    %- grab boundary events
    bds_ind = find(strcmp({EEG.event.type},'boundary'));
    bds_ind = bds_ind(bds_ind > valid_idxs(i));
    bds_ind = bds_ind(bds_ind < valid_idxs(i+1));
    bds_events = EEG.event(bds_ind);
    fs = fields(bds_events);
    fsPrev = fields(events_out{1});
    %- delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); out = [out{:}];
    delFs = fs(~out);
    if ~isempty(fsPrev) && any(~out)
        fprintf("Removing fields %s\n",delFs{:})
        for j = 1:length(delFs)
            bds_events = rmfield(bds_events,delFs{j});
        end
    end
    %- unravel events
    events_out = events_out(~cellfun(@isempty,events_out));
    split = cellfun(@(x) [[]; x], events_out);
    %- concatenate trials and wrap up loop iteration  
    tmp_trials{trial_cnt} = [split, bds_events];       
    trial_cnt = trial_cnt + 1;
end
EEG.event = [tmp_trials{:}];
%- let EEGLAB rearrange the event order
EEG = eeg_checkset(EEG,'eventconsistency');
%- epoch
EEG = pop_epoch( EEG,{cond_char},[-window_len/2,window_len/2],...
        'newname',sprintf('Merged_Datasets_%s_Epochs',EEG.subject),...
        'epochinfo','yes');

%% ===================================================================== %%
%## SUBFUNCTIONS
function [event] = create_event_entry(varargin)
    dt_tmp = datetime;
    dt_tmp.Format = 'MMddyyyy';
    event = [];
    event.latency   = varargin{1}; % necessary
    event.duration  = varargin{2}; % necessary
    event.type      = varargin{3}; % necessary (eeglab: epoching field)
    event.code      = varargin{4}; % necessary
    event.urevent   = sprintf('js_%s',dt_tmp); % necessary
    event.datetime  = varargin{5}; % necessary
    event.cond      = varargin{6}; % necessary
end
end
