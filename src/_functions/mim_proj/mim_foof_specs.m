function [] = mim_foof_specs(STUDY,)
%MIM_FOOF_SPECS Summary of this function goes here
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
%- custom params
colormap_ersp = linspecer; %othercolor('RdYlBu11');
% colormap_ersp = colormap_ersp(end:-1:1,:);
colormap(colormap_ersp);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'fooof_struct',@isnumeric);
addRequired(p,'cluster_load_ind',@isnumeric);
addRequired(p,'des_i',@isnumeric);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,cond_test,warping_times,cluster_ind,cluster_load_ind,des_i,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%

fooof_results{g}{cond,group}{i} = fooof(specfreqs, specdata_nolog(:,i), f_range, settings, return_model);

end

