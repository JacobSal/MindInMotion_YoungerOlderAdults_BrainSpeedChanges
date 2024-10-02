function [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,varargin)
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
%- (ADMIN PARAMS)
% COND_CHAR_FIELD = 'cond';
% %- sliding window params
% window_length = 6; % sliding window length in seconds
% PERCENT_OVERLAP = 0.0; % percent overlap between epochs
% APPROX_TRIAL_LENGTH = 3*60; % seconds
% %- event timewarp params
% EPOCH_TIME_LIMITS = [-1,3]; %[-1,3]; %[-0.5,5]; % [-1,3] captures gait events well , [-0.5,5] captures gait events poorly
% EVENTS_TIMEWARP = {'RHS','LTO','LHS','RTO','RHS'};
% STD_TIMEWARP = 3; %2.5;
% BASELINE_TIMELIMITS = [];
% if do_sliding_window
%     COND_CHARS = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
% else
%     COND_CHARS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% end

%## EPOCH PARAMS
DEF_EPOCH_PARAMS = struct('epoch_method','timewarp',...
    'percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'baseline_time_lims',[-0.5,4.5-2],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'cond_field','cond',...
    'appx_cond_len',3*60,...
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'do_recalc_epoch',true);
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'EPOCH_PARAMS',DEF_EPOCH_PARAMS,@(x) validate_struct(x,DEF_EPOCH_PARAMS))
parse(p,EEG,varargin{:});
%## SET DEFAULTS
%- PARAMETER
EPOCH_PARAMS = p.Results.EPOCH_PARAMS;
EPOCH_PARAMS = set_defaults_struct(EPOCH_PARAMS,DEF_EPOCH_PARAMS);

%% ===================================================================== %%
switch EPOCH_PARAMS.epoch_method
    case 'sliding_window'
        %- (MIND IN MOTION) sliding window
        ALLEEG = cell(1,length(EPOCH_PARAMS.slide_cond_chars)); 
        timewarp_struct = cell(1,length(EPOCH_PARAMS.slide_cond_chars));
        window_length = EPOCH_PARAMS.epoch_time_lims(2)-EPOCH_PARAMS.epoch_time_lims(1);
        for i = 1:length(EPOCH_PARAMS.slide_cond_chars)
            fprintf(1,'\n==== %s: Processing trial %s ====\n',EEG.subject,EPOCH_PARAMS.slide_cond_chars{i});
            [TMP_EEG] = sliding_window_epoch(EEG,EPOCH_PARAMS.slide_cond_chars{i},window_length,EPOCH_PARAMS.percent_overlap,...
                EPOCH_PARAMS.cond_field,EPOCH_PARAMS.appx_cond_len);
            %- check to make sure a number isn't the first character
            chk = regexp(EPOCH_PARAMS.slide_cond_chars{i},'\d');
            if any(chk)
                EPOCH_PARAMS.slide_cond_chars{i} = sprintf('x%s',EPOCH_PARAMS.slide_cond_chars{i});
            end
            TMP_EEG.etc.epoch.type = 'sliding_window';
            TMP_EEG.filename = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,EPOCH_PARAMS.slide_cond_chars{i});
            TMP_EEG.etc.epoch.condition = EPOCH_PARAMS.slide_cond_chars{i};
            TMP_EEG.etc.epoch.epoch_limits = window_length;
            TMP_EEG.etc.epoch.perc_overlap = EPOCH_PARAMS.percent_overlap;
            TMP_EEG.etc.epoch.parse_var = 'sliding_window';
            TMP_EEG.timewarp = struct([]);
            TMP_EEG.condition = EPOCH_PARAMS.slide_cond_chars{i};
            ALLEEG{i} = TMP_EEG;
            timewarp_struct{i} = struct([]);
        end
    case 'timewarp'
        %- (MIND IN MOTION) Gait Timewarp
        ALLEEG = cell(1,length(EPOCH_PARAMS.gait_trial_chars)); 
        timewarp_struct = cell(1,length(EPOCH_PARAMS.gait_trial_chars));
        TMP_EEG = pop_epoch( EEG, {  'RHS'  }, EPOCH_PARAMS.epoch_time_lims, 'newname', sprintf('Merged datasets %s epochs',EEG.subject), 'epochinfo', 'yes'); 
        TMP_EEG = eeg_checkset(TMP_EEG);
        if ~isempty(EPOCH_PARAMS.baseline_time_lims)
            TMP_EEG = pop_rmbase(TMP_EEG, EPOCH_PARAMS.baseline_time_lims ,[]); % Remove baseline from an epoched dataset.[-1500 2998] = baseline latency range, is it removing the mean during each gait cycle
        end
        TMP_EEG = eeg_checkset(TMP_EEG);
        for i = 1:length(EPOCH_PARAMS.gait_trial_chars)
            fprintf(1,'\n==== %s: Processing trial %s ====\n',TMP_EEG.subject,EPOCH_PARAMS.gait_trial_chars{i});
            %- GAIT TIMEWARP: MIM specific function
            TMP_TMP_EEG = mim_timewarp_epoch(TMP_EEG,EPOCH_PARAMS.gait_trial_chars{i},...
                EPOCH_PARAMS.tw_events,EPOCH_PARAMS.tw_stdev);
            TMP_TMP_EEG.etc.epoch_type = 'gait_timewarp';
            TMP_TMP_EEG.filename = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_TMP_EEG.subject,EPOCH_PARAMS.gait_trial_chars{i});
            TMP_TMP_EEG.condition = EPOCH_PARAMS.gait_trial_chars{i};
            %- check to make sure a number isn't the first character
            chk = regexp(EPOCH_PARAMS.gait_trial_chars{i},'\d');
            if any(chk)
                EPOCH_PARAMS.gait_trial_chars{i} = sprintf('x%s',EPOCH_PARAMS.gait_trial_chars{i});
            end
            TMP_TMP_EEG.etc.epoch.condition = EPOCH_PARAMS.gait_trial_chars{i};
            TMP_TMP_EEG.etc.epoch.epoch_limits = EPOCH_PARAMS.epoch_time_lims;
            ALLEEG{i} = TMP_TMP_EEG;
            timewarp_struct{i} = TMP_TMP_EEG.timewarp;
        end
        timewarp_struct = cellfun(@(x) [[]; x], timewarp_struct);
    otherwise
        error('Choose a valide EPOCH_PARAMS.epoch_method\n');
end
%- concatenate ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
%## TIME
toc
end
%% ===================================================================== %%
%## SUBFUNCTION
function [EEG] = mim_timewarp_epoch(EEG,cond_char,events_tw,std_tw)
    %- seconds to epoch relative to first RHS
    EEG = pop_selectevent(EEG, 'cond',cond_char,'deleteevents','off','deleteepochs','on','invertepochs','off'); 
    %- setup timewarp structure
    timewarp = make_timewarp(EEG,events_tw,'baselineLatency',0, ...
            'maxSTDForAbsolute',std_tw,...
            'maxSTDForRelative',std_tw);
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