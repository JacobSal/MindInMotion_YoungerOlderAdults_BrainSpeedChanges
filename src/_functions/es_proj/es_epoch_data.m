function [outputArg1,outputArg2] = ES_epoch_data(inputArg1,inputArg2)
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
%-
CONDITION_CHAR = '';
uniqConds = unique({EEG.event.cond});
cc_validFcn = @(x) assert(ischar(x) && any(strcmp(uniqConds,x)));
EPOCH_TIMES = [];
et_validFcn = @(x) assert(isnumeric(x) || isempty(x));
%-
TRIAL_LENGTH = 3*60;
tl_validFcn = @(x) assert(isnumeric(x),'TRIAL_LENGTH');
PERCENT_OVERLAP = 0.0;
po_validFcn = @(x) assert(isnumeric(x),'PERCENT_OVERLAP');
TRIAL_BEG_CHAR = 'TrialStart';
tbc_validFcn = @(x) assert(ischar(x),'TRIAL_BEG_CHAR');
TRIAL_END_CHAR = 'TrialEnd';
tec_validFcn = @(x) assert(ischar(x),'TRIAL_END_CHAR');
TRIAL_PARSE_CHAR = 'type';
tpc_validFcn = @(x) assert(ischar(x),'TRIAL_PARSE_CHAR');
CONDITION_PARSE_CHAR = 'cond';
cpc_validFcn = @(x) assert(ischar(x),'CONDITION_PARSE_CHAR');
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'epoch_type',@ischar)
%## OPTIONAL
addOptional(p,'condition_char',CONDITION_CHAR,cc_validFcn)
addOptional(p,'epoch_times',EPOCH_TIMES,et_validFcn);
%## PARAMETERS
addParameter(p,'TRIAL_LENGTH',TRIAL_LENGTH,tl_validFcn);
addParameter(p,'PERCENT_OVERLAP',PERCENT_OVERLAP,po_validFcn);
addParameter(p,'TRIAL_BEG_CHAR',TRIAL_BEG_CHAR,tbc_validFcn);
addParameter(p,'TRIAL_END_CHAR',TRIAL_END_CHAR,tec_validFcn);
addParameter(p,'TRIAL_PARSE_CHAR',TRIAL_PARSE_CHAR,tpc_validFcn);
addParameter(p,'CONDITION_PARSE_CHAR',CONDITION_PARSE_CHAR,cpc_validFcn);
%## PARSE
parse(p, EEG, epoch_type, varargin{:});
%## SET DEFAULTS
TRIAL_LENGTH = p.Results.TRIAL_LENGTH;
PERCENT_OVERLAP = p.Results.PERCENT_OVERLAP;
TRIAL_BEG_CHAR = p.Results.TRIAL_BEG_CHAR;
TRIAL_END_CHAR = p.Results.TRIAL_END_CHAR;
TRIAL_PARSE_CHAR = p.Results.TRIAL_PARSE_CHAR;
CONDITION_PARSE_CHAR = p.Results.CONDITION_PARSE_CHAR;
subjStr = EEG.subject;
%% ==================================================================== %%
% EEGLAB EPOCHING
% load specific EEG and amica results
fprintf('Epoching data for %s.\n', subjStr);    
% remove urevent (waste of space)
EEG.urevent = [];
switch epoch_type
    case 'Constant'
        %- Reassign events in EEG for constant epoching
        if isempty(condition_char) || isempty(epoch_times)
            error(['Error. Define ''condStr'' when using ''Constant'' parseVar\n',...
                   'Make sure epochTimeLimits is also defined']);
        else
            nameOfCond = sprintf('new_%s',condition_char);
        end
        %- find conditions that match input string
        tmp = strcmp({EEG.event.(CONDITION_PARSE_CHAR)},condition_char);
        tmpVal1 = strcmp({EEG.event.(TRIAL_PARSE_CHAR)},TRIAL_BEG_CHAR);
        tmpVal2 = strcmp({EEG.event.(TRIAL_PARSE_CHAR)},TRIAL_END_CHAR);
        validTimes = find(tmp & (tmpVal1 | tmpVal2));
        if (rem(length(validTimes),2) ~= 0)
            val = max(find(tmp));
            if (val == validTimes)
                val = val+1;
            end
            validTimes = [validTimes val];
            warning('TrialEnd not found, appending last step in condition as end of trial: index %i',val)
        end
        spc = (abs(epoch_times(2)-epoch_times(1))/2)*(1-PERCENT_OVERLAP);
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
        EEG = pop_epoch( EEG, {  nameOfCond  },epoch_times,'newname', sprintf('Merged_Datasets_%s_Epochs',subjStr), 'epochinfo', 'yes');
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

end

