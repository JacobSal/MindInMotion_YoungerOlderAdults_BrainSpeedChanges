function [EEG] = subj_epochData(PATHS, EEG, parseVar, condStr, epochTimeLimits, varargin)
%EPOCHSUBJDATA Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
%   NOTE (04/11/2022): IMPORTANT, not sure if this function is epoching
%   appropriately. Seems to select 2 gait cycles for 1 epoch, in some
%   cases, and is deleting a lot of epochs due to exclusion critera. Might
%   need to revist how we are defining gait events. Also there seems to be
%   substantial overlap between epoch'd gait cycles...

%## TIME
tic
%## DEFINE DEFAULTS
uniqConds = unique({EEG.event.cond});
Defaults = {3*60,...
            0.0,...
            'TrialStart',...
            'TrialEnd',...
            'type',...
            'cond'};
p = inputParser;
%## REQUIRED
validCond = (@(x) ischar(x) && any(strcmp(uniqConds,x)));
addRequired(p,'PATHS',@isstruct)
addRequired(p,'EEG',@isstruct)
addRequired(p,'parseVafr',@ischar)
addRequired(p,'condStr',validCond)
addRequired(p,'epochTimeLimits',@isnumeric);
%## OPTIONAL
addOptional(p,'trialLength',Defaults{1},@isnumeric);
addOptional(p,'percentOverlap',Defaults{2},@isnumeric);
addOptional(p,'trialBeginStr',Defaults{3},@ischar);
addOptional(p,'trialEndStr',Defaults{4},@ischar);
addOptional(p,'eventTrialParser',Defaults{5},@ischar);
addOptional(p,'eventCondParser',Defaults{6},@ischar);
%## PARAMETERS
%## PARSE
parse(p, PATHS, EEG, parseVar, condStr, epochTimeLimits, varargin{:});
%## SET DEFAULTS
TRIAL_LENGTH = p.Results.trialLength;
PER_OVERLAP = p.Results.percentOverlap;
TRIAL_BEGIN_STR = p.Results.trialBeginStr;
TRIAL_END_STR = p.Results.trialEndStr;
EVENT_TRIAL_PARSER = p.Results.eventTrialParser;
EVENT_COND_PARSER = p.Results.eventCondParser;
subjStr = EEG.subject;
%% ==================================================================== %%
% EEGLAB EPOCHING
% load specific EEG and amica results
fprintf('Epoching data for %s.\n', subjStr);    
% remove urevent (waste of space)
EEG.urevent = [];
switch parseVar
    case 'MIM_Gait_Warp'
        %## EPOCH DATA FOR GAIT CYCLES (ALL TRIALS): (RHS to RHS plus extra on ends)
        %seconds to epoch relative to first RHS
        epochTimeLimits = [-1 3]; 
        EEG = pop_epoch( EEG, {  'RHS'  }, epochTimeLimits, 'newname', sprintf('Merged datasets %s epochs',subjectString), 'epochinfo', 'yes');
        
        %## TIMEWARP GAIT CYCLES
        %setup timewarp structure
        events = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
        timewarp = make_timewarp(EEG,events,'baselineLatency',0, ...
                'maxSTDForAbsolute',3,...
                'maxSTDForRelative',3); 
        %subject specific warpto (later use to help calc grand avg warpto across subjects)
        timewarp.warpto = median(timewarp.latencies);        
        goodepochs  = sort([timewarp.epochs]);
        %probably not needed?
        EEG = eeg_checkset(EEG);   
        sedi = setdiff(1:length(EEG.epoch),goodepochs);
        %reject outlier strides
        EEG = pop_select( EEG,'notrial',sedi);
        %make sure timewarp latencies are available to us later when calculating warped ersps
        EEG.timewarp = timewarp;        
        %## GET TRIALS (if asked nicely)
%         [EEG,SUBJ] =
%         getTrialsByName(PATHS,EEG,subjStr,SUBJ,condStr,'savePath',savePath);
%         % (09/04/2022) needs updating before reinitialization.
    case 'Constant'
        %- Reassign events in EEG for constant epoching
        
        if isempty(condStr) || isempty(epochTimeLimits)
            error(['LOCALERROR: Define ''condStr'' when using ''Constant'' parseVar\n',...
                   'Make sure epochTimeLimits is also defined']);
        else
            nameOfCond = sprintf('new_%s',condStr);
        end
        %- find conditions that match input string
        tmp = strcmp({EEG.event.(EVENT_COND_PARSER)},condStr);
        tmpVal1 = strcmp({EEG.event.(EVENT_TRIAL_PARSER)},TRIAL_BEGIN_STR);
        tmpVal2 = strcmp({EEG.event.(EVENT_TRIAL_PARSER)},TRIAL_END_STR);
        validTimes = find(tmp & (tmpVal1 | tmpVal2));
        if (rem(length(validTimes),2) ~= 0)
            val = max(find(tmp));
            if (val == validTimes)
                val = val+1;
            end
            validTimes = [validTimes val];
            warning('TrialEnd not found, appending last step in condition as end of trial: index %i',val)
        end
        spc = (abs(epochTimeLimits(2)-epochTimeLimits(1))/2)*(1-PER_OVERLAP);
        %- initiate loop
        trialNum = 1;
        tmpEvent = EEG.event;
        tmpTrials = [];
        for i = 1:2:length(validTimes)
            % Loop iters
            split = [];
            % determine indices of trial events
            tmpLats =  [tmpEvent.latency];
            t1 = tmpLats >= tmpEvent(validTimes(i)).latency ;
            t2 = tmpLats <= tmpEvent(validTimes(i+1)).latency;
            % get headernames of event struct
            f = fieldnames(tmpEvent)';
            f{2,1} = {};
            % create new blank entry
            newE = struct(f{:});
            % determine start and end of trials
            idx = find(t1&t2);
            % split current event structure to put in new events
            % define constants in event struct  
            lt1 = tmpEvent(idx(1)).latency+EEG.srate*spc;
            dt1 = tmpEvent(idx(1)).datetime;
            dts1 = tmpEvent(idx(1)).datestr;
            ure1 = tmpEvent(idx(1)).urevent;
            % generate urevent numbers
            cnt = idx(1);
            while isempty(ure1)
                ure1 = tmpEvent(cnt).urevent;
                cnt = cnt + 1;
            end
            %## TRIAL APPENDING LOOP
            cnt = 1;
            intv = (lt1:EEG.srate*spc:(lt1+EEG.srate*TRIAL_LENGTH));
            % trial begin
            newE = createEventEntry(tmpEvent(validTimes(i)).latency,... % latency
                            1,... % duration
                            0,... % channel
                            [],... % bvtime
                            [],... % bvmknum
                            tmpEvent(validTimes(i)).type,... % type
                            nameOfCond,... % code
                            ure1,... % urevent
                            dt1,... % datetime
                            dts1,... % datestr
                            sprintf('%s_%s_%i',subjStr,nameOfCond,trialNum),... % trialName
                            nameOfCond); % cond
            split = [split, newE];
            % events inbetween
            for j = intv
                if j < (tmpEvent(validTimes(i+1)).latency - EEG.srate*spc)
                    newE = createEventEntry(j,... % latency
                            1,... % duration
                            0,... % channel
                            [],... % bvtime
                            [],... % bvmknum
                            nameOfCond,... % type
                            nameOfCond,... % code
                            (ure1)+cnt,... % urevent
                            dt1,... % datetime
                            dts1,... % datestr
                            sprintf('%s_%s_%i',subjStr,nameOfCond,trialNum),... % trialName
                            nameOfCond); % cond
                    split = [split, newE];
                    cnt = cnt+1;
                end
            end
            % trial End
            newE = createEventEntry(tmpEvent(validTimes(i+1)).latency,... % latency
                            1,... % duration
                            0,... % channel
                            [],... % bvtime
                            [],... % bvmknum
                            tmpEvent(validTimes(i+1)).type,... % type
                            nameOfCond,... % code
                            (ure1)+cnt,... % urevent
                            dt1,... % datetime
                            dts1,... % datestr
                            sprintf('%s_%s_%i',subjStr,nameOfCond,trialNum),... % trialName
                            nameOfCond); % cond
            split = [split, newE];
            % concatenate trials and wrap up loop iteration            
            trialNum = trialNum + 1;
            tmpTrials = horzcat(tmpTrials,split);              
        end
        % reconfigure EEG.event struct
        split1 = EEG.event(1:(validTimes(1)-1));
        split2 = EEG.event((validTimes(end)+1):end);
        fsS = fields(split1);
        fsT = fields(tmpTrials);
        if length(fsS) < length(fsT)
            compare = cellfun(@(x) any(strcmp(x,fsS)),fsT,'UniformOutput',false); compare = [compare{:}];
            delFs = fsT(~compare);
            for j = 1:length(delFs)
                fprintf('Adding field ''%s'' to EEG.event\n',delFs{j});
                vs1 = cell(size(split1));
                vs2 = cell(size(split2));
                [split1.(delFs{j})] = vs1{:};
                [split2.(delFs{j})] = vs2{:};
            end
        elseif length(fsT) < length(fsS)
            compare = cellfun(@(x) any(strcmp(x,fsT)),fsS,'UniformOutput',false); compare = [compare{:}];
            delFs = fsS(~compare);
            for j = 1:length(delFs)
                fprintf('Adding field ''%s'' to EEG.event',delFs{j});
                vs1     = cell(size(tmpTrials));
                [tmpTrials.(delFs{j})] = vs1{:};
            end
        end
        EEG.event = horzcat(split1, tmpTrials, split2);
        EEG = pop_epoch( EEG, {  nameOfCond  },epochTimeLimits,'newname', sprintf('Merged_Datasets_%s_Epochs',subjStr), 'epochinfo', 'yes');
    otherwise
        error('Please provide a epoch you''d like to study.')
end
end

%% ===================================================================== %%
function [event] = createEventEntry(varargin)
    event = [];
    event.latency   = varargin{1}; 
    event.duration  = varargin{2}; 
    event.channel   = varargin{3}; 
    event.bvtime    = varargin{4};
    event.bvmknum   = varargin{5}; 
    event.type      = varargin{6};
    event.visible   = char([]);
    event.code      = varargin{7}; 
    event.urevent   = varargin{8};
    event.datetime  = varargin{9}; 
    event.datestr   = varargin{10}; 
    event.trialName = varargin{11}; 
    event.cond      = varargin{12};
end
