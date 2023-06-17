function [fig_handle] = mim_plot_ersp(ersp_data,ers_times,ersp_freqs,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Jacob Salminen
% Code Date: 06/09/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ersp_data',@isstruct);
addRequired(p,'ersp_times',@isstruct);
addRequired(p,'ersp_freqs',@isnumeric);
%## OPTIONAL
addOptional(p,'designnumber',designnumber,@isnumeric);
%## PARAMETER
parse(p,ersp_data,ersp_times,ersp_freqs,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## ASSIGN PARAMETERS FOR READING
%-
%## TIME
toc
end

