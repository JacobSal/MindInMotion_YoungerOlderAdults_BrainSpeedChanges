function [SUBJ] = my_MNI4DIPFIT(PATHS, SUBJ, varargin)
%MAKEFEM4DIPFIT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct)
addRequired(p,'SUBJ',@isstruct)
%## OPTIONAL
parse(p, PATHS, SUBJ, varargin{:});
%## SET DEFAULTS
%## PARAMS & PATHS
%% ===================================================================== %%
warning('Using BEM MNI headmodel...');
% load standard MRI and headmodel Volume.
path2BEM    = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
% load standard elec locs
pathEA = fullfile(path2BEM,'elec');

%% Update SUBJSTRUCT
SUBJ.pathsEEG.headmodel(2).label     = 'vol';
SUBJ.pathsEEG.headmodel(2).filepath  = path2BEM;
SUBJ.pathsEEG.headmodel(2).filename  = 'standard_vol.mat';
SUBJ.pathsEEG.headmodel(3).label     = 'mri';
SUBJ.pathsEEG.headmodel(3).filepath  = path2BEM;
SUBJ.pathsEEG.headmodel(3).filename  = 'standard_mri.mat';
SUBJ.pathsEEG.headmodel(4).label     = 'elec';
SUBJ.pathsEEG.headmodel(4).filepath  = pathEA;
SUBJ.pathsEEG.headmodel(4).filename  = 'standard_1005.elc';
fprintf('Headmodel assigned to SUBJSTRUCT for %s.\n', SUBJ.subjStr)
%## time
toc
end

