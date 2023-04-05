function [events] = mim_trial_event_override(subj_char,trial_name,excel_fPath,varargin)
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
% Copyright (C) Ryan Downey
% Copyright (C) Chang Liu (07/14/2022)
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu (03/16/2023)
% persistent MasterTable
%## DEFINE DEFAULTS
RELOAD_SHEET = true;
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'subj_char',@ischar);
addRequired(p,'trial_name',@ischar);
addRequired(p,'excel_fPath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,subj_char,trial_name,excel_fPath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS
persistent MasterTable
%% ===================================================================== %%
if isempty(excel_fPath)
    disp('couldn''t find the excel doc so the function can''t work')
end

if isempty(MasterTable) || RELOAD_SHEET
    MasterTable = readtable(excel_fPath,'Sheet','trial_indices','Range','B2:V127'); %loads it up
    MasterTable.Properties.VariableNames(1:21) = {'sorting_num','subject_code','SP_0p5_1',...
        'SP_0p5_2',	'SP_0p25_1','SP_0p25_2','SP_0p75_1','SP_0p75_2','SP_1p0_1',	'SP_1p0_2',	...
        'TM_flat_1','TM_flat_2','TM_high_1','TM_high_2','TM_low_1',	'TM_low_2',	'TM_med_1',	'TM_med_2',	'Rest','MotorImagery_1','MotorImagery_2'};
end
SubjectMatch = contains(MasterTable{:,'subject_code'},subj_char); %look thru the subj to find match
SubjectMatchInd = find(SubjectMatch); %grab corresponding index

%## CREATE EVENTS 
if isempty(SubjectMatchInd)
    events = []; % no changes necessary
else
    table_ind = ~cellfun(@isempty,regexpi(MasterTable.Properties.VariableNames, trial_name));
    if sum(table_ind) > 1
        events = cell(1,sum(table_ind)*2);
        cnt = 1;
        sub_trial = find(table_ind);
        for ind = 1:length(sub_trial)
            res = MasterTable{SubjectMatchInd,sub_trial(ind)};
%             disp(res);
            try
                res = str2num(res{1});
                res(1);
                res(2);
            catch
                continue;
            end
            if ~isnan(res(1))
                events{cnt} = create_event_entry(res(1),...
                        1,'TrialStart','override',[],trial_name);
            end
            if ~isnan(res(2))
                events{cnt+1} = create_event_entry(res(2),...
                        1,'TrialEnd','override',[],trial_name);
            end
            cnt = cnt + 2;
        end
    else
        res = MasterTable{SubjectMatchInd,table_ind}; %result as a cell
        if isnan(res)
            events = [];
            return;
        end
        events = cell(1,2);
        res = str2num(res{1});
        events{1} = create_event_entry(res(1),...
                1,'TrialStart','override',[],[]);
        events{2} = create_event_entry(res(2),...
                1,'TrialEnd','override',[],[]);
    end
    if ~isempty(events)
        events = events(~cellfun(@isempty,events));
    end
end
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

