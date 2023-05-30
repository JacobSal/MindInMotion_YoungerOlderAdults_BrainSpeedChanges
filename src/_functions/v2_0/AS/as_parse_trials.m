function [ALLEEG,timewarp_struct] = as_parse_trials(EEG,epoch_limits,varargin)
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

%% DEFINE DEFAULTS
%## TIME
tic
%## PARAMS
%- (ADMIN PARAMS)

%- event timewarp params
BASELINE_LATENCY_MS = 100;
% COND_FIELD2PARSE = 'bounces';
CONDLABEL_CHARS = {'competitive','cooperative','moving_serve','stationary_serve'};
BOUNCES_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; % {'Serve_Human'}
TYPE_CHARS = {'Subject_hit'}; %{'Subject_receive'};
STD_TIMEWARP = 3;
EVENTS_TIMEWARP = {'Subject_hit','Subject_receive','Subject_hit'}; %{'Subject_receive','Subject_hit','Subject_receive'};
%- sliding window params
DO_SLIDING_WINDOW = false;
REGEXP_BOUNCES = {'Human','BM'};
EVENT_FIELD_CONDITION = 'bounces';
APPROX_TRIAL_LENGTH = 3*60; % seconds
PERCENT_OVERLAP = 0;
WINDOW_LENGTH = 5;
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'epoch_limits',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,@ischar);
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@ischar);
parse(p,EEG,epoch_limits,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
WINDOW_LENGTH = p.Results.WINDOW_LENGTH;
%% ===================================================================== %%
%- empty ALLEEG structure for repopulating
if DO_SLIDING_WINDOW
    %## SLIDING WINDOW CODE
    ALLEEG = cell(1,length(REGEXP_BOUNCES)); 
    timewarp_struct = cell(1,length(REGEXP_BOUNCES));
    cnt = 1;
    %- NOTE: (03/29/2023) JS, may still be buggy for amanda's data
    %will need to make changes to this function to adapt to
    %amanda's event data.
    for i = 1:length(REGEXP_BOUNCES)
        [TMP_EEG] = as_sliding_window_epoch(EEG,REGEXP_BOUNCES{i},WINDOW_LENGTH,...
            PERCENT_OVERLAP,EVENT_FIELD_CONDITION,APPROX_TRIAL_LENGTH);
        TMP_EEG.etc.epoch_type='sliding-window';
        %## STORE IN ALLEEG
        ALLEEG{cnt} = TMP_EEG;
        timewarp_struct{cnt} = struct([]);
    end
else
    %## EVENT TIMEWARPING CODE
    ALLEEG = cell(length(TYPE_CHARS)*length(BOUNCES_CHARS),1); 
    timewarp_struct = cell(1,length(TYPE_CHARS)*length(BOUNCES_CHARS));
    cnt = 1;
    for j = 1:length(TYPE_CHARS)
        %{
        %## Extract Epochs & Remove Baseline
        %## Amanda's data already epoched so this code chunk is not needed 
        TMP_EEG = pop_epoch(EEG,TYPE_CHARS(j),epoch_limits,...
            'newname',sprintf('Merged datasets %s epochs',EEG.subject),'epochinfo', 'yes');
        TMP_EEG = eeg_checkset(TMP_EEG);
        %- Remove baseline 150ms before receive or hit
        TMP_EEG = pop_rmbase(TMP_EEG,[epoch_limits(1)*1000 epoch_limits(2)*1000-epoch_limits(2)*1000-BASELINE_LATENCY_MS] ,[]);
        TMP_EEG = eeg_checkset(TMP_EEG);
        %}
        for i = 1:length(BOUNCES_CHARS)
            fprintf(1,'\n==== Processing ''%s'' ====\n',BOUNCES_CHARS{i});
            %{
            %## Timewarp Condition
            %## Amanda's data is already timewarped so this is not needed
            [TMP_TMP_EEG] = as_timewarp_epoch(TMP_EEG,BOUNCES_CHARS{i},...
                EVENTS_TIMEWARP,STD_TIMEWARP,COND_FIELD2PARSE);
            %}
            %## select bounce events
            TMP_TMP_EEG = pop_selectevent(EEG,'bounces',BOUNCES_CHARS{i},...
                'deleteevents','off','deleteepochs','on','invertepochs','off'); 
            %## select condition events
%             TMP_TMP_EEG = pop_selectevent(EEG,'condlabel',CONDLABEL_CHARS{i},...
%                 'deleteevents','off','deleteepochs','on','invertepochs','off'); 
            %- set parameters
            TMP_TMP_EEG.etc.epoch_type = sprintf('timewarp-%s',TYPE_CHARS{j});
            TMP_TMP_EEG.condition = [BOUNCES_CHARS{i} '_' TYPE_CHARS{j}];
            TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,BOUNCES_CHARS{i},TYPE_CHARS{j});
            %## STORE IN ALLEEG
            ALLEEG{cnt} = TMP_TMP_EEG;
            timewarp_struct{cnt} = TMP_TMP_EEG.timewarp;
            cnt=cnt+1;
        end
    end
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
timewarp_struct = cellfun(@(x) [[]; x], timewarp_struct);
%##
toc
end
%% ===================================================================== %%
%## SUBFUNCTION
function [EEG] = as_timewarp_epoch(EEG,cond_char,EVENTS_TIMEWARP,STD_TIMEWARP,COND_FIELD2PARSE)

%- seconds to epoch relative to first RHS
EEG = pop_selectevent(EEG,COND_FIELD2PARSE,cond_char,'deleteevents','off','deleteepochs','on','invertepochs','off'); 
%- setup timewarp structure
timewarp = make_timewarp(EEG,EVENTS_TIMEWARP,'baselineLatency',0, ...
        'maxSTDForAbsolute',STD_TIMEWARP,...
        'maxSTDForRelative',STD_TIMEWARP);
%subject specific warpto (later use to help calc grand avg warpto across subjects)
timewarp.warpto = median(timewarp.latencies);        
goodepochs = sort([timewarp.epochs]);
%probably not needed?
EEG = eeg_checkset(EEG);   
sedi = setdiff(1:length(EEG.epoch),goodepochs);
%reject outlier hits 
EEG = pop_select( EEG,'notrial',sedi);
%- store timewarp structure in EEG
EEG.timewarp = timewarp;
end
%% SUBFUNCTION 
function [EEG] = as_sliding_window_epoch(EEG,cond_char,window_len,percent_overlap,COND_CHAR_FIELD,APPROX_TRIAL_LENGTH)
%## Looking at cooperative vs competitive vs ball_machie
%- find conditions that match input string
tmp_all = strcmp({EEG.event.(COND_CHAR_FIELD)},cond_char);
% tmp_all = contains({EEG.event.(COND_CHAR_FIELD)},'Human');
fprintf('Using all events for condition ''%s'' as 1 trial',cond_char);
trial_start = find(tmp_all,1,'first');
trial_end = find(tmp_all,1,'last');
EEG.event(trial_start).type = 'tmp_start';
EEG.event(trial_start).cond = cond_char;
EEG.event(trial_end).type = 'tmp_end';
EEG.event(trial_end).cond = cond_char;
tmp_all = strcmp({EEG.event.(COND_CHAR_FIELD)},cond_char);
tmp_start = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},'tmp_start');
tmp_end = strcmp({EEG.event.(EVENT_FIELD_TRIAL)},'tmp_end');
valid_idxs = find(tmp_all & (tmp_start | tmp_end));
%-
spc = (abs(window_len)/2)*(1-percent_overlap);
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