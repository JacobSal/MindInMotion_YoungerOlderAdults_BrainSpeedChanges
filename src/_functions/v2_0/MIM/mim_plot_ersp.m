function [] = mim_plot_ersp(STUDY,ALLEEG,warping_times,save_dir,varargin)
%MIM_PLOT_ERSP Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
FREQ_RANGE = [4 100];
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'warping_times',@isnumeric);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,warping_times,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%


end

