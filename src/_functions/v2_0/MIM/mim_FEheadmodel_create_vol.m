function [error_code] = mim_FEheadmodel_create_vol(working_dir,subject_str,varargin)
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
DO_PLOTTING = true;
error_code = 0; %#ok<NASGU>
%- working directory containing the ctf_fiducials.mat & mri_acpc_rs.mat &
%mri_acpc.mat & CustomElectrodeLocations.txt
errorMsg = 'Value must be CHAR. working directory containing the vol.mat & elec_aligned.mat'; 
% wd_validFcn = @(x) assert(ischar(x) && exist([x filesep 'ctf_fiducials.mat'],'file') && exist([x filesep 'mri_acpc_rs.mat'],'file') && exist([x filesep 'mri_acpc.mat'],'file') && exist([x filesep 'CustomElectrodeLocations.txt'],'file'),errorMsg);
wd_validFcn = @(x) assert(ischar(x)  && exist([x filesep 'ctf_fiducials.mat'],'file')...
    && exist([x filesep sprintf('%s_masks_contr.nii.gz',subject_str)],'file')...
    && exist([x filesep 'CustomElectrodeLocations.txt'],'file')...
    && exist([x filesep 'mri_acpc.mat'],'file'),errorMsg);
fprintf('Checking working_dir (%s) for ''subject_str''_masks_contr.nii.gz & CustomElectrodeLocations.txt\n',working_dir);
fprintf('ctf_fiducials.mat check: %i\n',exist([working_dir filesep 'ctf_fiducials.mat'],'file'))
% fprintf('mri_acpc_rs.mat check: %i\n',exist([working_dir filesep 'mri_acpc_rs.mat'],'file'))
fprintf('mri_acpc.mat check: %i\n',exist([working_dir filesep 'mri_acpc.mat'],'file'))
fprintf('''subject_str''_masks_contr.nii.gz check: %i\n',exist([working_dir filesep sprintf('%s_masks_contr.nii.gz',subject_str)],'file'))
fprintf('CustomElectrodeLocations.txt check: %i\n',exist([working_dir filesep 'CustomElectrodeLocations.txt'],'file'))
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'working_dir',wd_validFcn);
addRequired(p,'subject_str',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,working_dir,subject_str,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
fprintf('Running dipole fitting on directory: %s\n',working_dir);
%- load mri in acpa & ctf coordinate systems & fiducial marks
% tmp = load([working_dir filesep 'mri_acpc_rs.mat']);
% mri_acpc_rs = tmp.mri_acpc_rs;
tmp = load([working_dir filesep 'mri_acpc.mat']);
mri_acpc = tmp.mri_acpc;
tmp = load([working_dir filesep 'ctf_fiducials.mat']);
ctf_fiducials = tmp.ctf_fiducials;
%- Load the electrodes after digitized
chanloc_scan_folder = [working_dir filesep 'CustomElectrodeLocations.txt'];
chanlocs = readtable(chanloc_scan_folder);% Same output text file from getchalocs.
chanlocs.Properties.VariableNames = {'labels','X','Y','Z'};
elec.chanpos(:,1) = [chanlocs.X];
elec.chanpos(:,2) = [chanlocs.Y];
elec.chanpos(:,3) = [chanlocs.Z];
elec.elecpos      = elec.chanpos;
elec.label(:,1)   = [chanlocs.labels]';
%%
gunzip([working_dir filesep sprintf('%s_masks_contr.nii.gz',subject_str)])
simnibs_mask = ft_read_mri([working_dir filesep sprintf('%s_masks_contr.nii',subject_str)]);
simnibs_mask.coordsys = 'acpc';

segmented = simnibs_mask;
segmented.white = simnibs_mask.anatomy == 1;
segmented.gray = simnibs_mask.anatomy == 2;
segmented.csf = (simnibs_mask.anatomy == 3 | simnibs_mask.anatomy == 8);% csf + ventricles
segmented.skull = simnibs_mask.anatomy == 4;
segmented.scalp = (simnibs_mask.anatomy == 5 | simnibs_mask.anatomy == 7);%skin and eye
segmented.air = simnibs_mask.anatomy == 6; 
segmented = rmfield(segmented,'anatomy');
seg_i_headreco = ft_datatype_segmentation(segmented,'segmentationstyle','indexed');
%% FIELDTRIP SEGMENTATION
% cfg             = [];
% cfg.spmmethod   = 'old';%new method output is weird
% cfg.units       = 'mm';
% cfg.output      = {'gray','white','csf','skull','scalp'};
% % cfg.inputfile   = [working_dir filesep 'mri_acpc_rs.mat'];
% % cfg.outputfile  = [working_dir filesep 'segmentedmri.mat'];
% segmentedmri    = ft_volumesegment(cfg, mri_acpc_rs);
%% FIELDTRIP CREATE MESH
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,seg_i_headreco);
fprintf('Saving mesh file\n');
save([working_dir filesep 'mesh.mat'],'mesh')
%% FIELDTRIP CREATE CONDUCTIVITY VOLUME (SIMBIO)
cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = zeros(1,5);
scale = 1;
% order follows mesh.tissyelabel , CAUTIOUS!!!! OMg this is not the same order as in the segmentation
cfg.conductivity(strcmp(mesh.tissuelabel,'csf')) = 1.65*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'gray')) = 0.33*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'scalp')) = 0.33*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'skull')) = 0.01*scale; %0.0042*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'white')) = 0.126*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'air')) = 2.5*10^(-14)*scale;
vol = ft_prepare_headmodel(cfg, mesh);
fprintf('Saving vol file\n');
save([working_dir filesep 'vol.mat'],'vol','-v6') 

%% FIELDTRIP CONFIRM ELECTRODE ALIGNMENT
%- Convert the fiducial position from voxel into CTF 
nas = ctf_fiducials.nas;
lpa = ctf_fiducials.lpa;
rpa = ctf_fiducials.rpa;
%- grab transformation from acpc mri
vox2head = mri_acpc.transform;
%- apply warping to get ctf coordinates of fiducials
nas = ft_warp_apply(vox2head, nas, 'homogenous');
lpa = ft_warp_apply(vox2head, lpa, 'homogenous');
rpa = ft_warp_apply(vox2head, rpa, 'homogenous');
%- create a structure similar to a template set of electrodes
fid.chanpos       = [nas; lpa; rpa];       % CTF head coordinates of fiducials
fid.label         = {'nas','lhj','rhj'};    % use the same labels as those in elec
fid.unit          = 'mm';                  % use the same units as those in mri
fid.elecpos       = fid.chanpos;           % otherwise the electroderealign cannot find elecpos
%- alignment
cfg               = [];
cfg.viewmode      = 'surface';
cfg.method        = 'fiducial';
% cfg.method        = 'interactive';%interactive doesn't work well.
cfg.headshape     = vol;
cfg.elec          = elec;                  % the electrodes we want to align
cfg.elecstyle     = {'facecolor','red'};
cfg.headmodelstyle= {'facecolor','b','edgecolor','none','facealpha',0.4};
cfg.template      = fid;                   % the template we want to align to
cfg.fiducial      = {'nas', 'lhj', 'rhj'};  % labels of fiducials in fid and in elec
elec_aligned_init = ft_electroderealign(cfg);
save([working_dir filesep 'elec_aligned_init.mat'],'elec_aligned_init') 

%## Refined Electrode Alignment
cfg               = [];
cfg.method        = 'project';
cfg.elec          = elec_aligned_init;
cfg.headshape     = vol;
elec_aligned      = ft_electroderealign(cfg);
fprintf('Saving elec_aligned file\n');
save([working_dir filesep 'elec_aligned.mat'],'elec_aligned') 
%% PLOTS
if DO_PLOTTING
    %{
    %## (PLOT 1)
    seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');

    cfg              = [];
    cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial -`д´-
    cfg.anaparameter = 'anatomy';
    cfg.funcolormap  = jet(6); % distinct color per tissue
    cfg.location     = 'center';
    cfg.atlas        = seg_i;    % the segmentation can also be used as atlas

    % check segmentation quality - It's kind of bad?!!
    ft_sourceplot(cfg, seg_i);%this plot cannot be generated...I don't know why
    saveas(gcf,[working_dir filesep 'segmentation.fig']);
    %}
    cfg              = [];
    cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
    cfg.anaparameter = 'anatomy';
    cfg.funcolormap  = linspecer(7); % distinct color per tissue
    cfg.location     = 'center';
    cfg.atlas        = seg_i_headreco;    % the segmentation can also be used as atlas
    ft_sourceplot(cfg, seg_i_headreco);
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[working_dir filesep sprintf('ft_sourceplot.fig')]);
    saveas(fig_i,[working_dir filesep sprintf('ft_sourceplot.jpg')]);
    %## (PLOT 2)
    figure;
    ft_plot_mesh(mesh, 'surfaceonly', 'yes','facecolor','b','edgecolor', 'none', 'facealpha', 0.4);
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_mesh.fig')]);
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_mesh.jpg')]);
    %## (PLOT 4) Initial alignment of electrodes using fiducial marks
    figure;
    hold on;
    ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    camlight
    ft_plot_sens(elec_aligned_init,'style','.r');
    ft_plot_sens(fid,'style','xb');%plot fiducial points
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_sens_1.fig')]);
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_sens_1.jpg')]);
    %## (PLOT 3) final alignment of electrodes after projecting to scalp
    figure;
    hold on;
    ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    camlight
    ft_plot_sens(elec_aligned ,'style','.r');
    ft_plot_sens(fid,'style','xb');%plot fiducial points
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_sens_2.fig')]);
    saveas(fig_i,[working_dir filesep sprintf('ft_plot_sens_2.jpg')]);
    
end
error_code = 1;
end
