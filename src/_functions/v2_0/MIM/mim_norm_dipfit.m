function [EEG,dipfit_fem_model] = mim_norm_dipfit(eeg_fpath,eeg_fname,mri_fpath,dipfit_fpath,varargin)
%MIM_NORM_DIPFIT Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Chang Liu, 
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
NORMALIZE_MRI = true;
DO_PLOT = false;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'eeg_fpath');
addRequired(p,'eeg_fname');
addRequired(p,'mri_fpath');
addRequired(p,'dipfit_fpath');
%## OPTIONAL
%## PARAMETER
parse(p,eeg_fpath,eeg_fname,mri_fpath,dipfit_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%- load MRI
fprintf('Loading MRI info...\n');
tmp = load(mri_fpath);
mri_acpc_rs = tmp.mri_acpc_rs;
%- assign amica_folder
fprintf('Loading EEG...\n');
EEG = pop_loadset('filepath',eeg_fpath,'filename',eeg_fname);
%- load dipfit_fem.mat
fprintf('Loading dipfit_fem.mat...\n');
tmp = load(dipfit_fpath);
dipfit_fem = tmp.SAVEVAR;
clear tmp
%## Reformat Dipfit Structure
empty_dip_struct = struct('posxyz',[],'momxyz',[],'rv',NaN,'diffmap',[],'sourcepot',[],'datapot',[]);
EEG.dipfit_fem = [];
EEG.dipfit_fem.model = empty_dip_struct;
dipfit_fem_pos = zeros(length([dipfit_fem.component]),3);
for i=1:length([dipfit_fem.component])
    %- 
    if ~isempty(dipfit_fem(i).dip)
        EEG.dipfit_fem.model(i).posxyz = dipfit_fem(i).dip.pos;
        EEG.dipfit_fem.model(i).momxyz = reshape(dipfit_fem(i).dip.mom, 3, length(dipfit_fem(i).dip.mom)/3)';
        if ~isempty(dipfit_fem(i).dip.rv)
            EEG.dipfit_fem.model(i).rv     = dipfit_fem(i).dip.rv;
        else
            EEG.dipfit_fem.model(i).rv     = NaN;
        end
        %- 
        EEG.dipfit_fem.model(i).diffmap = dipfit_fem(i).Vmodel - dipfit_fem(i).Vdata;
        EEG.dipfit_fem.model(i).sourcepot = dipfit_fem(i).Vmodel;
        EEG.dipfit_fem.model(i).datapot   = dipfit_fem(i).Vdata;
        dipfit_fem_pos(i,:) = dipfit_fem(i).dip.pos;
    else
        EEG.dipfit_fem.model(i) = empty_dip_struct;
    end
end
%% NORMALIZE MRI
if NORMALIZE_MRI
    cfg = [];
    cfg.nonlinear = 'yes';
    cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
    mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
    %- PLOT
    if DO_PLOT
        ft_sourceplot(cfg,mri_norm);
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[eeg_fpath filesep sprintf('normalized_mri.fig')]);
        saveas(fig_i,[eeg_fpath filesep sprintf('normalized_mri.jpg')]);
    end
end
%% Convert the dipole location to MNI space
dipfit_fem_mnipos = ft_warp_apply(mri_norm.params,ft_warp_apply(mri_norm.initial,dipfit_fem_pos), 'individual2sn');
dipfit_fem_mni_voxinds = round(ft_warp_apply(pinv(mri_norm.transform), dipfit_fem_mnipos ));
for i=1:length([dipfit_fem.component])
    EEG.dipfit_fem.model(i).mnipos = dipfit_fem_mnipos(i,:);
    EEG.dipfit_fem.model(i).mni_voxinds = dipfit_fem_mni_voxinds(i,:);
    EEG.dipfit_fem.model(i).pos_old = EEG.dipfit_fem.model(i).posxyz;
    EEG.dipfit_fem.model(i).posxyz = EEG.dipfit_fem.model(i).mnipos;
end
%% save
EEG.dipfit = EEG.dipfit_fem;
dipfit_fem_model = EEG.dipfit;
end

