function [ALLEEG,timewarp_struct] = as_parse_trials(EEG,varargin)
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
EPOCH_METHOD = 'sliding_window';
STRUCT_FIELD_COND = 'bounces';
STRUCT_FIELD_EVENT = 'cond';
EVENT_CHARS = {};
EVENTS_TIMEWARP = {};
COND_CHARS = {};
REGEXP_COND = '\s';
STD_TIMEWARP = 3;
APPROX_TRIAL_LENGTH = 3*60; % seconds
PERCENT_OVERLAP = 0;
WINDOW_LENGTH = 0;
EPOCH_LENGTH_TIMELIM = [0,1];
%-
EPOCH_PARAMS = [];
EPOCH_PARAMS.epoch_method = EPOCH_METHOD;
EPOCH_PARAMS.event_chars = EVENT_CHARS;
EPOCH_PARAMS.events_timewarp = EVENTS_TIMEWARP;
EPOCH_PARAMS.cond_chars = COND_CHARS;
EPOCH_PARAMS.std_timewarp = STD_TIMEWARP;
EPOCH_PARAMS.struct_field_cond = STRUCT_FIELD_COND;
EPOCH_PARAMS.struct_field_event = STRUCT_FIELD_EVENT;
EPOCH_PARAMS.percent_overlap = PERCENT_OVERLAP;
EPOCH_PARAMS.window_length = WINDOW_LENGTH;
EPOCH_PARAMS.regexp_slidingwindow = REGEXP_COND;
EPOCH_PARAMS.approx_trial_length = APPROX_TRIAL_LENGTH;
EPOCH_PARAMS.epoch_length_timelim = EPOCH_LENGTH_TIMELIM;
%- event timewarp params
% TIMEWARP_PARAMS = {'epoch_method','timewarp',...
%     'event_chars',{},...
%     'events_timewarp',{},...
%     'cond_chars',{},...
%     'std_timewarp',3,...
%     'struct_field_cond',STRUCT_FIELD_COND,...
%     'struct_field_event',STRUCT_FIELD_EVENT};
% %- event-locked params
% EVENTLOCK_PARAMS = {'epoch_method','event_locked',...
%     'event_chars',{},...
%     'events_timewarp',{},...
%     'cond_chars',{},...
%     'std_timewarp',3,...
%     'struct_field_cond',STRUCT_FIELD_COND,...
%     'struct_field_event',STRUCT_FIELD_EVENT};
% %- sliding window params
% SLIDING_WINDOW_STRUCT = {'epoch_method','sliding_window',...
%     'regexp_bounces',{'Human','BM'},...
%     'percent_overlap',0,...
%     'window_length',5};
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'epoch_params',EPOCH_PARAMS,@isstruct);
% addParameter(p,'epoch_method',EPOCH_METHOD,@ischar);
% addParameter(p,'event_chars',EVENT_CHARS,@iscell);
% addParameter(p,'events_timewarp',EVENTS_TIMEWARP,@iscell);
% addParameter(p,'cond_chars',COND_CHARS,@iscell);
% addParameter(p,'std_timewarp',STD_TIMEWARP,@isnumeric);
% addParameter(p,'struct_field_cond',STRUCT_FIELD_COND,@ischar);
% addParameter(p,'struct_field_event',STRUCT_FIELD_EVENT,@ischar);
% addParameter(p,'percent_overlap',PERCENT_OVERLAP,@isnumeric);
% addParameter(p,'window_length',WINDOW_LENGTH,@isnumeric);
% addParameter(p,'regexp_slidingwindow',REGEXP_COND,@ischar);
% addParameter(p,'approx_trial_length',APPROX_TRIAL_LENGTH,@isnumeric);
% addParameter(p,'epoch_length_timelim',EPOCH_LENGTH_TIMELIM,@isnumeric);
parse(p,EEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
epoch_params = p.Results.epoch_params;
% event_chars = p.Results.event_chars;
% events_timewarp = p.Results.events_timewarp;
% cond_chars = p.Results.cond_chars;
% std_timewarp = p.Results.std_timewarp;
% struct_field_cond = p.Results.struct_field_cond;
% struct_field_event = p.Results.struct_field_event;
% percent_overlap = p.Results.percent_overlap;
% window_length = p.Results.window_length;
% regexp_slidingwindow = p.Results.regexp_slidingwindow;
% approx_trial_length = p.Results.approx_trial_length;
% epoch_length_timelim = p.Results.epoch_length_timelim;
EPOCH_PARAMS_FIELDS = {'epoch_method','event_chars','events_timewarp',...
    'cond_chars','std_timewarp','struct_field_cond','struct_field_event',...
    'percent_overlap','window_length','regexp_slidingwindow','approx_trial_length',...
    'epoch_length_timelim'};
for f_i = 1:length(EPOCH_PARAMS_FIELDS)
    if ~isfield(epoch_params,EPOCH_PARAMS_FIELDS{f_i})
        epoch_params.(EPOCH_PARAMS_FIELDS{f_i}) = EPOCH_PARAMS.(EPOCH_PARAMS_FIELDS{f_i});
    end
end
epoch_method = epoch_params.epoch_method;
event_chars = epoch_params.event_chars;
events_timewarp = epoch_params.events_timewarp;
cond_chars = epoch_params.cond_chars;
std_timewarp = epoch_params.std_timewarp;
struct_field_cond = epoch_params.struct_field_cond;
struct_field_event = epoch_params.struct_field_event;
percent_overlap = epoch_params.percent_overlap;
window_length = epoch_params.window_length;
regexp_slidingwindow = epoch_params.regexp_slidingwindow;
approx_trial_length = epoch_params.approx_trial_length;
epoch_length_timelim = epoch_params.epoch_length_timelim;
%% ===================================================================== %%
%- empty ALLEEG structure for repopulating
switch epoch_method
    case 'sliding_window'
        %## SLIDING WINDOW CODE
        ALLEEG = cell(1,length(regexp_slidingwindow)); 
        timewarp_struct = cell(1,length(regexp_slidingwindow));
        cnt = 1;
        %- NOTE: (03/29/2023) JS, may still be buggy for amanda's data
        %will need to make changes to this function to adapt to
        %amanda's event data.
        for i = 1:length(regexp_slidingwindow)
            [TMP_EEG] = sliding_window_epoch(EEG,regexp_slidingwindow{i},window_length,...
                percent_overlap,struct_field_cond,approx_trial_length);
            TMP_EEG.etc.epoch.type='sliding-window';
    %         TMP_EEG.condition = REGEXP_BOUNCES;
            TMP_EEG.condition = regexp_slidingwindow{i};
            %## STORE IN ALLEEG
            ALLEEG{cnt} = TMP_EEG;
            timewarp_struct{cnt} = struct([]);
        end
    case 'timewarp'
        %## EVENT TIMEWARPING CODE
        ALLEEG = cell(length(event_chars)*length(cond_chars),1); 
        timewarp_struct = cell(1,length(event_chars)*length(cond_chars));
        cnt = 1;
        for j = 1:length(event_chars)
            %{
            %## Extract Epochs & Remove Baseline
            %## Amanda's data already epoched so this code chunk is not needed 
            TMP_EEG = pop_epoch(EEG,event_chars(j),epoch_limits,...
                'newname',sprintf('Merged datasets %s epochs',EEG.subject),'epochinfo', 'yes');
            TMP_EEG = eeg_checkset(TMP_EEG);
            %- Remove baseline 150ms before receive or hit
            TMP_EEG = pop_rmbase(TMP_EEG,[epoch_limits(1)*1000 epoch_limits(2)*1000-epoch_limits(2)*1000-BASELINE_LATENCY_MS] ,[]);
            TMP_EEG = eeg_checkset(TMP_EEG);
            %}
            %{
            for i = 1:length(CONDLABEL_CHARS)
                fprintf(1,'\n==== Processing ''%s'' ====\n',CONDLABEL_CHARS{i});
                %## select condition events
                TMP_TMP_EEG = pop_selectevent(EEG,'condlabel',CONDLABEL_CHARS{i},...
                    'deleteevents','off','deleteepochs','on','invertepochs','off'); 
                %## Timewarp condition events
                [TMP_TMP_EEG] = as_timewarp_epoch(TMP_TMP_EEG,...
                    EVENTS_TIMEWARP,STD_TIMEWARP);
                %## set parameters
                TMP_TMP_EEG.etc.epoch_type = sprintf('timewarp-%s',event_chars{j});
                TMP_TMP_EEG.condition = [CONDLABEL_CHARS{i} '_' event_chars{j}];
                TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,CONDLABEL_CHARS{i},event_chars{j});
                %## STORE IN ALLEEG
                ALLEEG{cnt} = TMP_TMP_EEG;
                timewarp_struct{cnt} = TMP_TMP_EEG.timewarp;
                cnt=cnt+1;
            end
            %}
            for i = 1:length(cond_chars)
                fprintf(1,'\n==== Processing ''%s'' ====\n',cond_chars{i});
                %## Select Bounce Events
                TMP_TMP_EEG = pop_selectevent(EEG,'bounces',cond_chars{i},...
                    'deleteevents','off','deleteepochs','on','invertepochs','off');
                %- Remove baseline 150ms before receive or hit?
    %             TMP_TMP_EEG = pop_rmbase(TMP_TMP_EEG,[epoch_limits(1)*1000 epoch_limits(1)+500+150] ,[]);
    %             TMP_TMP_EEG = eeg_checkset(TMP_TMP_EEG);
                %## Timewarp Bounce events
                [TMP_TMP_EEG] = timewarp_epoch(TMP_TMP_EEG,...
                    events_timewarp,std_timewarp);
                %## set parameters
                TMP_TMP_EEG.etc.epoch_type = sprintf('timewarp-%s',event_chars{j});
                TMP_TMP_EEG.etc.epoch.condition = [cond_chars{i} '_' event_chars{j}];
    %             TMP_TMP_EEG.etc.epoch.epoch_limits = ;
                TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,cond_chars{i},event_chars{j});
                TMP_TMP_EEG.condition = [cond_chars{i} '_' event_chars{j}];
                %## STORE IN ALLEEG
                ALLEEG{cnt} = TMP_TMP_EEG;
                timewarp_struct{cnt} = TMP_TMP_EEG.timewarp;
                cnt=cnt+1;
            end
        end
    case 'event_locked'
        %## EVENT LOCKED TO STIM CODE
        ALLEEG = cell(length(event_chars)*length(cond_chars),1); 
        timewarp_struct = cell(1,length(event_chars)*length(cond_chars));
        cnt = 1;
        for i = 1:length(cond_chars)
            %## POP EVENTS OF INTEREST FROM ESTABLISHED EPOCHS
            TMP_EEG = pop_selectevent(EEG,struct_field_cond,cond_chars{i},...
                    'deleteevents','off','deleteepochs','on','invertepochs','off');
            TMP_EEG = eeg_checkset(TMP_EEG);
            for j = 1:length(event_chars)
                %## EPOCH EVENT OF INTEREST (SUB_EPOCH)
                fprintf(1,'\n==== Processing %s::''%s'' ====\n',EEG.subject,cond_chars{i});
                TMP_TMP_EEG = pop_epoch(TMP_EEG,event_chars(j),epoch_length_timelim,...
                    'newname',sprintf('Merged datasets %s epochs',EEG.subject),'epochinfo', 'yes');
                %## set parameters
                TMP_TMP_EEG.etc.epoch_type = sprintf('timelock-%s',event_chars{j});
                TMP_TMP_EEG.etc.epoch.condition = [cond_chars{i} '_' event_chars{j}];
    %             TMP_TMP_EEG.etc.epoch.epoch_length_timelim = ;
                TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,cond_chars{i},event_chars{j});
                TMP_TMP_EEG.condition = [cond_chars{i} '_' event_chars{j}];
                %## STORE IN ALLEEG
                ALLEEG{cnt} = TMP_TMP_EEG;
                timewarp_struct{cnt} = struct([]);
                cnt=cnt+1;
            end
        end
        %{
        for j = 1:length(event_chars)
            %## Extract Epochs & Remove Baseline
            %## Amanda's data already epoched so this code chunk is not needed 
            TMP_EEG = pop_epoch(EEG,event_chars(j),epoch_length_timelim,...
                'newname',sprintf('Merged datasets %s epochs',EEG.subject),'epochinfo', 'yes');
            TMP_EEG = eeg_checkset(TMP_EEG);
            %- Remove baseline 150ms before receive or hit
%             TMP_EEG = pop_rmbase(TMP_EEG,[epoch_length_timelim(1)*1000 0] ,[]);
%             TMP_EEG = eeg_checkset(TMP_EEG);
            %{
            for i = 1:length(CONDLABEL_CHARS)
                fprintf(1,'\n==== Processing ''%s'' ====\n',CONDLABEL_CHARS{i});
                %## select condition events
                TMP_TMP_EEG = pop_selectevent(EEG,'condlabel',CONDLABEL_CHARS{i},...
                    'deleteevents','off','deleteepochs','on','invertepochs','off'); 
                %## Timewarp condition events
                [TMP_TMP_EEG] = as_timewarp_epoch(TMP_TMP_EEG,...
                    events_timewarp,STD_TIMEWARP);
                %## set parameters
                TMP_TMP_EEG.etc.epoch_type = sprintf('timewarp-%s',event_chars{j});
                TMP_TMP_EEG.condition = [CONDLABEL_CHARS{i} '_' event_chars{j}];
                TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,CONDLABEL_CHARS{i},event_chars{j});
                %## STORE IN ALLEEG
                ALLEEG{cnt} = TMP_TMP_EEG;
                timewarp_struct{cnt} = TMP_TMP_EEG.timewarp;
                cnt=cnt+1;
            end
            %}
            for i = 1:length(cond_chars)
                fprintf(1,'\n==== Processing ''%s'' ====\n',cond_chars{i});
                %## Select Bounce Events
                TMP_TMP_EEG = pop_selectevent(TMP_EEG,struct_field_cond,cond_chars{i},...
                    'deleteevents','off','deleteepochs','on','invertepochs','off');
                %## set parameters
                TMP_TMP_EEG.etc.epoch_type = sprintf('timelock-%s',event_chars{j});
                TMP_TMP_EEG.etc.epoch.condition = [cond_chars{i} '_' event_chars{j}];
    %             TMP_TMP_EEG.etc.epoch.epoch_length_timelim = ;
                TMP_TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,cond_chars{i},event_chars{j});
                TMP_TMP_EEG.condition = [cond_chars{i} '_' event_chars{j}];
                %## STORE IN ALLEEG
                ALLEEG{cnt} = TMP_TMP_EEG;
                timewarp_struct{cnt} = struct([]);
                cnt=cnt+1;
            end
        end
        %}
    case 'sub_epoch'
        %## EVENT NEW EVENT LOCKED
        % (02/01/2024) JS, this does not support multiple event_chars parameter
        % right now.
        ALLEEG = cell(length(event_chars)*length(cond_chars),1); 
        timewarp_struct = cell(1,length(event_chars)*length(cond_chars));
        event_i = 1;
        cnt = 1;
        for i = 1:length(cond_chars)
            fprintf(1,'\n==== Processing %s: ''%s'' ====\n',EEG.subject,cond_chars{i});
            indsb = strcmp({EEG.event.(struct_field_cond)},cond_chars{i});
            epochb = [EEG.event(indsb).epoch];
            indse = diff([epochb,epochb(end)]);
            tmpe = epochb(logical(indse));
            epoch_keep = zeros(1,length(tmpe));
            epoch_field = sprintf('event%s',struct_field_cond);
            for j = 1:length(tmpe)
                lats = [EEG.epoch(tmpe(j)).eventlatency{:}];
                tmpi = (lats == 0);
                epoch_type = [EEG.epoch(tmpe(j)).(epoch_field){tmpi}];
                if strcmp(epoch_type,cond_chars{i})
                    epoch_keep(j) = tmpe(j);
                end
            end
            epoch_keep = epoch_keep(epoch_keep~=0);
            TMP_EEG = pop_select(EEG,'trial',epoch_keep);
            %## set parameters
            TMP_EEG.etc.epoch_type = sprintf('timelock-%s-%s',cond_chars{i},event_chars{event_i});
            TMP_EEG.etc.epoch.condition = [cond_chars{i} '_' event_chars{event_i}];
%             TMP_TMP_EEG.etc.epoch.epoch_length_timelim = ;
            TMP_EEG.filename = sprintf('%s_%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,cond_chars{i},event_chars{event_i});
            TMP_EEG.condition = [cond_chars{i} '_' event_chars{event_i}];
            %## STORE IN ALLEEG
            ALLEEG{cnt} = TMP_EEG;
            timewarp_struct{cnt} = struct([]);
            cnt = cnt + 1;
        end
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
if ~isempty(timewarp_struct{1})
    timewarp_struct = cellfun(@(x) [[], x], timewarp_struct);
end
%##
toc
end
%% ===================================================================== %%
%## SUBFUNCTION
function [EEG] = timewarp_epoch(EEG,events_timewarp,STD_TIMEWARP)
%- setup timewarp structure
timewarp = make_timewarp(EEG,events_timewarp,'baselineLatency',0, ...
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
function [EEG] = sliding_window_epoch(EEG,cond_char,window_len,percent_overlap,...
    cond_char_field,approx_trial_len)
%## Extract Trial Boundaries
spc = (abs(window_len)/2)*(1-percent_overlap);
%- find conditions that match input string
tmp_all = strcmp({EEG.event.(cond_char_field)},cond_char);
% tmp_all = contains({EEG.event.(COND_CHAR_FIELD)},'Human');
fprintf('Using all events for condition ''%s'' as 1 trial\n',cond_char);
trial_start = find(tmp_all,1,'first');
trial_end = find(tmp_all,1,'last');
EEG.event(trial_start).type = 'tmp_start';
EEG.event(trial_start).cond = cond_char;
EEG.event(trial_end).type = 'tmp_end';
EEG.event(trial_end).cond = cond_char;
tmp_all = strcmp({EEG.event.cond},cond_char);
%- option 1
tmp_start = strcmp({EEG.event.type},'tmp_start');
tmp_end = strcmp({EEG.event.type},'tmp_end');
valid_idxs = find(tmp_all & (tmp_start | tmp_end));
%- option 2
% valid_idxs = find(diff(tmp_all));
%## check for trial length and boundary events(errors occur if trial is cut up by boundary
%event insert).
% for i = 1:2:length(valid_idxs)
%     dt = EEG.event(valid_idxs(i)).latency - EEG.event(valid_idxs(i+1)).latency;
%     if dt/1000 < approx_trial_len
%         tmpt = EEG.event(valid_idxs(i)).latency+(approx_trial_len*EEG.srate);
%         tmpe = create_event_entry(tmpt,1,...
%                     'appended_tmp_end','trial_mark',dt,cond_char);
%         EEG.event = [EEG.event; ];
%     end
% end
% EEG = eeg_checkset(EEG,'eventconsistency');
%- initiate loop
trial_cnt = 1;
tmp_event = EEG.event;
tmp_trials = cell(1,length(valid_idxs)/2);
%## TRIAL APPENDING LOOP
for i = 1:2:length(valid_idxs)
    % split current event structure to put in new events
    % define constants in event struct  
    lat_1   = tmp_event(valid_idxs(i)).latency+EEG.srate*spc;
    lat_2   = lat_1+EEG.srate*approx_trial_len;
    intervals = (lat_1:EEG.srate*spc*2:lat_2);
    events_out = cell(1,length(intervals));
    code_char = sprintf('trial_%i_cond_%s',trial_cnt,cond_char);
%     dt = tmp_event(valid_idxs(i)).datetime;
    dt = datetime;
    dt.Format = 'MMddyyyy';
    %- beginning boundary event
    events_out{1} = create_event_entry(intervals(1),...
                    window_len/2,'boundary',[],[],[]);
    %- inbetween events
    for j = 2:length(intervals)
        chk_intv = intervals(j) < (tmp_event(valid_idxs(i+1)).latency - EEG.srate*spc);
        if j == 1
            event_type = 'tmp_start';
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
                    1,'tmp_end',code_char,dt,cond_char); 
    %- ending boundary event
    events_out{j+1} = create_event_entry(intervals(j)+window_len/2,...
                    window_len/2,'boundary',[],[],[]);
    %- grab boundary events
    bds_ind = find(strcmp({EEG.event.type},'boundary'));
    bds_ind = bds_ind(bds_ind > valid_idxs(i));
    bds_ind = bds_ind(bds_ind < valid_idxs(i+1));
    bds_events = EEG.event(bds_ind);
    cnt = length(events_out)+1;
    for j = 1:length(bds_events)
        events_out{cnt} = create_event_entry(bds_events(j).latency,...
                    bds_events(j).duration,bds_events(j).type,code_char,dt,cond_char);
        cnt=cnt+1;
    end
    %- unravel events
    events_out = events_out(~cellfun(@isempty,events_out));
    split = cellfun(@(x) [[]; x], events_out);
    %- concatenate trials and wrap up loop iteration  
    tmp_trials{trial_cnt} = split; %[split; bds_events];       
    trial_cnt = trial_cnt + 1;
end
EEG.event = [tmp_trials{:}];
%- let EEGLAB rearrange the event order
EEG = eeg_checkset(EEG,'eventconsistency');
%- epoch
EEG = pop_epoch( EEG,{cond_char},[-window_len/2,window_len/2],...
        'newname',sprintf('Merged_Datasets_%s_Epochs',EEG.subject),...
        'epochinfo','yes');
end
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