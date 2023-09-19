function [ALLEEG,timewarp_struct] = nj_parse_trials(EEG,epoch_method,varargin)
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
%## TIME
tic
%## DEFINE DEFAULTS
%- (ADMIN PARAMS)
COND_CHAR_FIELD = 'type';
%- sliding window params
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
APPROX_TRIAL_LENGTH = 3*60; % seconds
%- gait timewarp params
EPOCH_TIME_LIMITS = [-1,3];
TIMEWARP_EVENTS = {'RHS','LTO','LHS','RTO','RHS'};
STD_TIMEWARP = 3; %2.5;
BASELINE_TIMELIMITS = [];
%- epoching params
switch epoch_method
    case 'sliding_window_standing'
        COND_CHAR_FIELD = 'type';
        COND_CHARS = {'pre','post'};
    case 'sliding_window_gait'
        COND_CHAR_FIELD = 'cond';
        COND_CHARS = {'B1','B2','B3','P1','P2','SB1','SB2','none'};
    case 'gait_notimewarp'
        COND_CHARS = {'B1','B2','B3','P1','P2','SB1','SB2','none'};
    otherwise
        error('error. choose a valid EPOCH_METHOD');        
end
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'epoch_method',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'COND_CHARS',COND_CHARS,@iscell);
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@isnumeric);
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,@isnumeric);
addParameter(p,'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,@isnumeric);
addParameter(p,'TIMEWARP_EVENTS',TIMEWARP_EVENTS,@iscell);
addParameter(p,'STD_TIMEWARP',STD_TIMEWARP,@isnumeric);
addParameter(p,'BASELINE_TIMELIMITS',BASELINE_TIMELIMITS,(@(x) isnumeric(x) || isempty(x)));
parse(p,EEG,epoch_method,varargin{:});
%## SET DEFAULTS
%- PARAMETER
%* sliding window params
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
WINDOW_LENGTH = p.Results.WINDOW_LENGTH;
COND_CHARS = p.Results.COND_CHARS;
%* gait params
EPOCH_TIME_LIMITS = p.Results.EPOCH_TIME_LIMITS;
TIMEWARP_EVENTS = p.Results.TIMEWARP_EVENTS;
STD_TIMEWARP = p.Results.STD_TIMEWARP;
BASELINE_TIMELIMITS = p.Results.BASELINE_TIMELIMITS;
if isempty(BASELINE_TIMELIMITS)
    BASELINE_TIMELIMITS = [EPOCH_TIME_LIMITS(1),EPOCH_TIME_LIMITS(2)-1000*(1/EEG.srate)]; % time in milliseconds
end
%% ===================================================================== %%
%## STEP 2) EPOCHING
%* empty ALLEEG structure for repopulating
switch epoch_method
    case {'sliding_window_standing','sliding_window_gait'}
        %- (NJ new standing data 07/16/2023) make data continuous first to
        %avoid eeg_checksset error
        EEG = eeglab_makecontinuous(EEG);
        %- (NJ) sliding window
        ALLEEG = cell(1,length(COND_CHARS)); 
        timewarp_struct = cell(1,length(COND_CHARS));
        for i = 1:length(COND_CHARS)
            fprintf(1,'\n==== %s: Processing trial %s ====\n',EEG.subject,COND_CHARS{i});
            [TMP_EEG] = sliding_window_epoch(EEG,COND_CHARS{i},WINDOW_LENGTH,PERCENT_OVERLAP,...
                COND_CHAR_FIELD,APPROX_TRIAL_LENGTH);
            %- check to make sure a number isn't the first character
            chk = regexp(COND_CHARS{i},'\d');
            if any(chk)
                COND_CHARS{i} = sprintf('x%s',COND_CHARS{i});
            end
            TMP_EEG.etc.epoch.type = 'sliding_window';
            TMP_EEG.filename = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,COND_CHARS{i});
            TMP_EEG.etc.epoch.condition = COND_CHARS{i};
            TMP_EEG.etc.epoch.epoch_limits = WINDOW_LENGTH;
            TMP_EEG.etc.epoch.perc_overlap = PERCENT_OVERLAP;
            TMP_EEG.condition = COND_CHARS{i};
            ALLEEG{i} = TMP_EEG;
            timewarp_struct{i} = struct([]);
        end
    case 'gait_notimewarp'
        %- (NJ new standing data 07/16/2023) make data continuous first to
        %avoid eeg_checksset error
        EEG = eeglab_makecontinuous(EEG);
        %- (Noelle J) Gait No Timewarp
        ALLEEG = cell(1,length(COND_CHARS)); 
        timewarp_struct = cell(1,length(COND_CHARS));
        TMP_EEG = pop_epoch( EEG, {  'RHS'  }, EPOCH_TIME_LIMITS, 'newname', sprintf('Merged datasets %s epochs',EEG.subject), 'epochinfo', 'yes'); 
        TMP_EEG = eeg_checkset(TMP_EEG);
        TMP_EEG = pop_rmbase(TMP_EEG, BASELINE_TIMELIMITS ,[]); % Remove baseline from an epoched dataset.[-1500 2998] = baseline latency range, is it removing the mean during each gait cycle
        TMP_EEG = eeg_checkset(TMP_EEG);
        for i = 1:length(COND_CHARS)
            fprintf(1,'\n==== %s: Processing trial %s ====\n',TMP_EEG.subject,COND_CHARS{i});
            %- GAIT TIMEWARP: MIM specific function
            TMP_TMP_EEG = timewarp_epoch(TMP_EEG,COND_CHARS{i},...
                TIMEWARP_EVENTS,STD_TIMEWARP);
            TMP_TMP_EEG.etc.epoch_type = 'gait_timewarp';
            TMP_TMP_EEG.filename = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,COND_CHARS{i});
            TMP_TMP_EEG.condition = COND_CHARS{i};
%             %- check to make sure a number isn't the first character
%             chk = regexp(COND_CHARS{i},'\d');
%             if any(chk)
%                 COND_CHARS{i} = sprintf('x%s',COND_CHARS{i});
%             end
            TMP_TMP_EEG.etc.epoch.condition = COND_CHARS{i};
            TMP_TMP_EEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
            ALLEEG{i} = TMP_TMP_EEG;
            timewarp_struct{i} = TMP_TMP_EEG.timewarp;
        end
        timewarp_struct = cellfun(@(x) [[]; x], timewarp_struct);
    otherwise
        error('error. choose a valide epoch_method...');
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
end
%% ===================================================================== %%
%## SUBFUNCTION
function [EEG] = timewarp_epoch(EEG,cond_char,EVENTS_TIMEWARP,STD_TIMEWARP)
%- seconds to epoch relative to first RHS
EEG = pop_selectevent(EEG, 'cond',cond_char,'deleteevents','off','deleteepochs','on','invertepochs','off');
%- setup timewarp structure
timewarp = make_timewarp(EEG,EVENTS_TIMEWARP,'baselineLatency',0, ...
        'maxSTDForAbsolute',STD_TIMEWARP,...
        'maxSTDForRelative',STD_TIMEWARP);
%subject specific warpto (later use to help calc grand avg warpto across subjects)
timewarp.warpto = median(timewarp.latencies);        
goodepochs  = sort([timewarp.epochs]);
%probably not needed?
EEG = eeg_checkset(EEG);   
sedi = setdiff(1:length(EEG.epoch),goodepochs);
%reject outlier strides & 
EEG = pop_select( EEG,'notrial',sedi);
%- store timewarp structure in EEG
EEG.timewarp = timewarp;
end
%## SUBFUNCTION 
function [EEG] = sliding_window_epoch(EEG,cond_char,window_len,percent_overlap,...
    cond_char_field,approx_trial_len)
%## Looking at cooperative vs competitive vs ball_machie
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
tmp_start = strcmp({EEG.event.type},'tmp_start');
tmp_end = strcmp({EEG.event.type},'tmp_end');
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
%## 
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