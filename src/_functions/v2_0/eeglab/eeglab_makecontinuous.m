function [EEG] = eeglab_makecontinuous(EEG,varargin)
%eeglab_makecontinuous Summary of this function goes here
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
% Code Date: 07/16/2023, MATLAB 2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
%## OPTIONAL
%## PARAMETER
parse(p,EEG,varargin{:});
%## SET DEFAULTS
%- PARAMETER
%% ===================================================================== %%
if size(EEG.icaact,3)>1 
    EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),size(EEG.icaact,2)*size(EEG.icaact,3));
end
if size(EEG.data,3)>1 
     EEG.data = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
end
if EEG.trials > 1 
    EEG.trials = 1;
    EEG.pnts = size(EEG.data,2);
    tmp = (size(EEG.data,2)-1)/EEG.srate;
    EEG.xmin = 0;
    EEG.xmax = tmp;
end
if isfield(EEG.event,'epoch')
    EEG.event = rmfield(EEG.event,'epoch');
end
if length(EEG.epoch) > 1
    EEG.epoch = [];
end
EEG = eeg_checkset(EEG);
end

