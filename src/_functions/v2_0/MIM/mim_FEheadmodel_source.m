function [error_code] = mim_FEheadmodel_source(working_dir,varargin)
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
error_code = 0;
%- working directory containing the vol.mat & elec_aligned.mat
errorMsg = 'Value must be CHAR. working directory containing the vol.mat & elec_aligned.mat'; 
wd_validFcn = @(x) assert(ischar(x) && exist(x,'dir') && exist([x filesep 'vol.mat'],'file') && exist([x filesep 'elec_aligned.mat'],'file'),errorMsg);
fprintf('Checking working_dir (%s) for vol.mat & elec_aligned.mat\n',working_dir);
fprintf('vol.mat check: %i\n',exist([working_dir filesep 'vol.mat'],'file'))
fprintf('elec_aligned.mat check: %i\n',exist([working_dir filesep 'elec_aligned.mat'],'file'))
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'working_dir',wd_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,working_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
fprintf('Running sourcemodel creation on directory: %s\n',working_dir);
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
%% PREPARE VOL & ELEC
%- prepare volume and sensors (projects sensors to scalp)
fprintf('Preparing vol.mat and elec_aligned.mat\n');
if ~exist([working_dir filesep 'headmodel_fem_tr.mat'],'file')
    %## Override Simbio path to utilize parallel processing
%     ft_hastoolbox('simbio', 1);
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
%%
if ~exist([working_dir filesep 'sourcemodel.mat'],'file')
    fprintf('Making source model...\n');
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
    fprintf('Sourcemodel already generated for working_dir: %s',working_dir);
end

%## TIME
toc
%## EXIT
% exit(error_code);
end
