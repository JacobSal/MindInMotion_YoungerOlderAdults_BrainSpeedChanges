function [EEG] = eeglab_gait_timewarp_epoch(EEG,cond_char,epoch_limits,varargin)
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
EVENT_FIELD_CONDITION = 'cond';
%-
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'cond_char',@ischar);
addRequired(p,'epoch_limits',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'EVENT_FIELD_CONDITION',EVENT_FIELD_CONDITION,@ischar);
parse(p,EEG,cond_char,epoch_limits,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
EVENT_FIELD_CONDITION = p.Results.EVENT_FIELD_CONDITION;
%% ================================================================= %%
%- find trials for desired condition
tmp = strcmp({EEG.event.(EVENT_FIELD_CONDITION)},cond_char);
valid_points = [find(tmp,1,'first'),find(tmp,1,'last')];
%##  EPOCH DATA FOR GAIT CYCLES (TIME WARP)
%seconds to epoch relative to first RHS
events = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
EEG = pop_epoch( EEG, {'RHS'}, epoch_limits,...
        'eventindices',valid_points(1):valid_points(2),...
        'newname', sprintf('Merged datasets %s epochs',EEG.subject),...
        'epochinfo', 'yes');
%- setup timewarp structure
timewarp = make_timewarp(EEG,events,'baselineLatency',0, ...
        'maxSTDForAbsolute',3,...
        'maxSTDForRelative',3);
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

