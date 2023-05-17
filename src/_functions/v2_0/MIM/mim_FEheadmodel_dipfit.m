function [error_code] = mim_FEheadmodel_dipfit(working_dir,eeg_fpath,varargin)
%MIM_FEHEADMODEL_DIPFIT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu

%## TIME
tic
%## DEFINE DEFAULTS
POOL_SIZE = 20;
error_code = 0;
%- working directory containing the vol.mat & elec_aligned.mat
errorMsg = 'Value must be CHAR. working directory containing the vol.mat & elec_aligned.mat'; 
wd_validFcn = @(x) assert(ischar(x) && exist(x,'dir') && exist([x filesep 'vol.mat'],'file') && exist([x filesep 'elec_aligned.mat'],'file'),errorMsg);
%- EEG filepath
errorMsg = 'Value must be an EEGLAB EEG STRUCT.'; 
ef_validFcn = @(x) assert(ischar(x) && exist(x,'file'),errorMsg);
fprintf('Checking working_dir (%s) for vol.mat & elec_aligned.mat\n',working_dir);
fprintf('vol.mat check: %i\n',exist([working_dir filesep 'vol.mat'],'file'))
fprintf('elec_aligned.mat check: %i\n',exist([working_dir filesep 'elec_aligned.mat'],'file'))

%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'working_dir',wd_validFcn);
addRequired(p,'eeg_fpath',ef_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,working_dir,eeg_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
fprintf('Running dipole fitting on directory: %s\n',working_dir);
%- load vol & elec_aligned
tmp = load([working_dir filesep 'vol.mat']);
vol = tmp.vol;
tmp = load([working_dir filesep 'elec_aligned.mat']);
try
    elec_aligned = tmp.elec_aligned;
catch e
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
    elec_aligned = tmp.elec_aligned_init;
end
%% PREPARE VOLUME AND SENSORS
%- prepare volume and sensors (projects sensors to scalp)
fprintf('Preparing vol.mat and elec_aligned.mat\n');
if ~exist([working_dir filesep 'headmodel_fem_tr.mat'],'file')
    %## Override Simbio path to utilize parallel processing
    %- find simbio on path
    if ~ispc
        tmp = strsplit(path,':');
    else
        tmp = strsplit(path,';');
    end
    b1 = regexp(tmp,'simbio','end');
    b2 = tmp(~cellfun(@isempty,b1));
    fprintf('Current simbio path: %s\n',b2{1});
    %- assign _override directory
    b1 = regexp(tmp,'_override','end');
    b2 = tmp(~cellfun(@isempty,b1));
    override_out = b2{1};
    path(path,[override_out filesep 'simbio'])
    fprintf('Simbio path: %s\n',[override_out filesep 'simbio'])
    %- prepare volume and sensors
    [headmodel_fem_tr, ~] = ft_prepare_vol_sens(vol, elec_aligned);
    save([working_dir filesep 'headmodel_fem_tr.mat'],'headmodel_fem_tr','-v7.3');
else
    tmp = load([working_dir filesep 'headmodel_fem_tr.mat']);
    headmodel_fem_tr = tmp.headmodel_fem_tr;
end
%% CREATE SOURCEMODEL
if ~exist([working_dir filesep 'sourcemodel.mat'],'file')
    %- 
    cfg             = [];
    cfg.resolution  = 5; % changeable
    cfg.threshold   = 0.1;
    cfg.smooth      = 5;
    cfg.headmodel   = headmodel_fem_tr;
    cfg.tight       = 'yes';
    cfg.inwardshift = 0.5; % shifts dipoles away from surfaces
    cfg.unit        = 'mm';
    sourcemodel     = ft_prepare_sourcemodel(cfg);
    save([working_dir filesep 'sourcemodel.mat'],'sourcemodel');
else
    tmp = load([working_dir filesep 'sourcemodel.mat']);
    sourcemodel = tmp.sourcemodel;
end
%- Remove this field to force average referencing of leadfield matrix
try
    elec_aligned_fem    = rmfield(elec_aligned_fem,'tra'); 
catch e
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
end

% choose the available When the forward solution is computed, the lead 
% field matrix (= channels X source points matrix) is calculated for each 
% grid point taking into account the head model and the channel positions.
% remove unused elec_aligned_fem;
if ~exist([working_dir filesep 'leadfield_fem.mat'],'file')
    cfg             = [];
    cfg.grid        = sourcemodel;
    cfg.headmodel   = headmodel_fem_tr;
    cfg.elec        = elec_aligned_fem;
    cfg.reducerank  = 'no';
    % cfg.normalize   = 'column'; 
    leadfield_fem   = ft_prepare_leadfield(cfg);% This actually didn't take that long - ?1r without parallel processing
    save([working_dir filesep 'leadfield_fem.mat'],'leadfield_fem');
else
    tmp = load([working_dir filesep 'leadfield_fem.mat']);
    leadfield_fem = tmp.leadfield_fem;
end
clear sourcemodel elec_aligned vol
%% Dipole Fitting
%## NOTES
% REQUIRED
% cfg.model           = g_model; %'moving';
% cfg.nonlinear       = g_nonlinear; %'yes';
% cfg.numdipoles      = g_numdipoles; %1; %number of dipoles, if > 1 you need to define the symmetry variable
% cfg.resolution      = g_resolution; %10; %resolution of model in units <cfg.unit>, exponetially increases computation time
% cfg.unit            = g_unit; %units; %units of headmodel 
% cfg.gridsearch      = g_gridsearch; %'yes'; %gridsearch for initial params (increases processing time, but gives stronger fit).
% cfg.dipfit.maxiter  = 200;
% cfg.reducerank      = g_reducerank; %'no';
% cfg.spmversion      = g_spmversion; %'spm12'; %'cat12'?
% cfg.vol             = vol; %the headmodel created from the MRI or MNI
% cfg.senstype        = g_senstype; %sensor type
% cfg.elec            = elec; %the channels you want to consider (keeping as is will use all EEG sensors)
% OPTIONS
% cfg.warpmni         = true;
% cfg.channel         = eegChan;
% cfg.leadfield       = leadfieldFEM;
%## COARSE FIT
tmp = strsplit(eeg_fpath,filesep);
fName = tmp{end};
fPath = strjoin(tmp(1:end-1),filesep);
%- load EEG
EEG = pop_loadset('filepath',fPath,'filename',fName);
ftEEG = eeglab2fieldtrip(EEG,'componentanalysis','dipfit');
clear EEG 
%- ft_dipolefitting
% cfg = [];
% cfg.numdipoles    =  1;
% cfg.headmodel     = headmodel_fem_tr;
% cfg.sourcemodel   = leadfield_fem;
% cfg.elec          = elec_aligned; %elec_aligned_fem;
% cfg.dipfit.metric = 'rv';
% cfg.nonlinear     = 'no';
% cfg.component     = 1:size(ftEEG.topo,2);
% % cfg.component     = 2;
% % for each component scan the whole brain with dipoles using FIELDTRIPs
% % dipolefitting function
% dipfit_fem_grid        = ft_dipolefitting(cfg,ftEEG);
% save([working_dir filesep 'dipfit_fem_grid.mat'],'dipfit_fem_grid','-v6');

%## NONLINEAR FIT
if ~exist([working_dir filesep 'dipfit'], 'dir')
    mkdir([working_dir filesep 'dipfit'])
end
sources = cell(size(ftEEG.topo,2),1);
parfor (comp_i = 1:size(ftEEG.topo,2),POOL_SIZE)
    cfg = [];
    cfg.numdipoles      =  1;
    cfg.headmodel       = headmodel_fem_tr;
    cfg.sourcemodel     = leadfield_fem;
    cfg.elec            = elec_aligned; %elec_aligned;
    cfg.dipfit.metric   = 'rv';
    cfg.nonlinear       = 'yes';
    cfg.component       = comp_i
    source        = ft_dipolefitting(cfg,ftEEG);
    sources{comp_i} = source;
    display(source)
%     par_save(source,[working_dir filesep 'dipfit'],sprintf('source_%i.mat',comp_i));
end
sources = cellfun(@(x) [[]; x], sources);
par_save(sources,[working_dir filesep 'dipfit'],'dipfit_struct.mat');
%## TIME
toc
%## EXIT
% exit(error_code);
end
