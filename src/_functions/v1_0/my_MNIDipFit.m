function [EEG] = my_MNIDipFit(PATHS,SUBJ,EEG,varargin)
%FINDDIPOLES Summary of this function goes here
%   This function uses fieldtrip functions to generate a dipole fit for
%   your participant population. 
%   
%   IN: 
%   OUT: 
%   IMPORTANT: 
%   (04/04/2022) in consideration of time I'm deciding to forego
%   reoptimizing this function till a later date. 

%## GLOBALS
%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {[]};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct)
addRequired(p,'SUBJ',@isstruct)
addRequired(p,'EEG',@isstruct)
%## PARAMETERS
parse(p, PATHS, SUBJ, EEG, varargin{:});
%## SET DEFAULTS
DIP_NUM = 1;
DIP_PLOT = 'off';
%## SAVE FOLDER HANDLER
%## PARAMS
%% ===================================================================== %%
for prcsNum = 1:length(SUBJ.pathsEEG.headmodel)
    if strcmp(SUBJ.pathsEEG.headmodel(prcsNum).label, 'elec')
        pathEA  = fullfile(SUBJ.pathsEEG.headmodel(prcsNum).filepath, SUBJ.pathsEEG.headmodel(prcsNum).filename);
    end
    if strcmp(SUBJ.pathsEEG.headmodel(prcsNum).label, 'vol')
        pathHM  = fullfile(SUBJ.pathsEEG.headmodel(prcsNum).filepath, SUBJ.pathsEEG.headmodel(prcsNum).filename);
    end
    if strcmp(SUBJ.pathsEEG.headmodel(prcsNum).label, 'mri')
        pathMRI = fullfile(SUBJ.pathsEEG.headmodel(prcsNum).filepath, SUBJ.pathsEEG.headmodel(prcsNum).filename);
    end
end

%## ADD DIPFIT STRUCT TO EEG STRUCT
warning('Using EEGLAB pop_multifit.m function');
% we also convert to fieldtrip format here
fprintf('pop_dipfit_settings...\n');
EEG = pop_dipfit_settings( EEG, 'hdmfile',pathHM,'coordformat','MNI',...
    'mrifile',pathMRI,'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] );

%## DIPFIT (see. ft_dipolefitting())
fprintf('pop_multifit...\n');
EEG = pop_multifit(EEG,[],'dipoles',DIP_NUM,'dipplot',DIP_PLOT);
EEG = iclabel(EEG, 'lite');

%## TIME
toc
end
%ARCHIVE_CODE-------------------------------------------------------------%
% This code is implemented from fieldtrip, but requires a lot of
% computation time. EEGLAB implements (above) a guessing function which aids in
% predicting where the dipole 'could' be and this speeds up the dipole
% finding computation time. This could be useful if we need to implement
% dipole finding on the hypergator or on a gpu as it seems easier to
% parallelize than the EEGLAB functions.
%## PREPARE VOLUME AND SENSOR MODELS
%   This may or may not be helpful for controlling for any issues/offsets
%   between the headmodel and sensor positions
% [headModel, sens] = ft_prepare_vol_sens(headModel, ftEEG.elec, 'channel', ftEEG.elec.label(eegChan));
% ftEEG.elec = sens;

%## COMPUTE LEADFIELD (forward problem)
%   This may or may not be helpful for computation time. Seems that the
%   ft_dipolefitting function computes the leadfield at every estimated
%   dipole position than calculates the inverse (using ft_inverse_dipolefit)
% leadF = ft_compute_leadfield([0,0,0]*length(idxKeep),ftEEG.elec,headModel);

%## COMPUTE DIPOLES (inverse problem)
% LINK (dipfit_nonlinear: https://github.com/sccn/dipfit/blob/master/dipfit_nonlinear.m
% Call hieracrhcy: (forward) ft_dipolefitting => ft_compute_leadfield/ft_prepare_vol_sens => leadfield_simbio (or whatever method you choose) =>  =>
% => sb_transfer => sb_calc_vecx => sb_solve => pcg
% (inverse) ft_dipolefitting => => ft_inverse_dipolefit




