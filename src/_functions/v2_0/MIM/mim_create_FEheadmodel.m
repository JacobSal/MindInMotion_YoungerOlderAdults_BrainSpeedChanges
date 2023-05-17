function [outputArg1,outputArg2] = mim_create_FEheadmodel(inputArg1,inputArg2)
%MIM_CREATE_FEHEADMODEL Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu

%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- cell of alleegs
errorMsg = 'Value must be a CELL ARRAY of EEG structures (EEGLAB).'; 
ta_validFcn = @(x) assert(iscell(x),errorMsg);
%- array of subject
errorMsg = 'Value must be a vector of 1''s and 0''s with length equal to length(ALLEEG).'; 
as_validFcn = @(x) assert(isnumeric(x),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'TMP_ALLEEG',ta_validFcn);
addRequired(p,'STUDY',@isstruct);
addRequired(p,'rmv_subjs',as_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,TMP_ALLEEG,STUDY,rmv_subjs,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## fieldtrip segmentation
cfg           = [];
cfg.spmmethod = 'old';%new method output is weird
cfg.output    = {'gray','white','csf','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri_acpc_rs);

% saves to current folder, instead save to patient file? for now saved to SourceLocalization folder
disp(segmentedmri)

% --- plot
seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');

cfg              = [];
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial -`д´-
cfg.anaparameter = 'anatomy';
cfg.funcolormap  = jet(6); % distinct color per tissue
cfg.location     = 'center';
cfg.atlas        = seg_i;    % the segmentation can also be used as atlas

% check segmentation quality - It's kind of bad?!!
ft_sourceplot(cfg, seg_i);%this plot cannot be generated...I don't know why
saveas(gcf,fullfile(save_headmodel_folder,'segmentation.fig'));

%% Create mesh
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';

mesh = ft_prepare_mesh(cfg,segmentedmri);

figure;
ft_plot_mesh(mesh, 'surfaceonly', 'yes','facecolor','b','edgecolor', 'none', 'facealpha', 0.4);
save(fullfile(save_headmodel_folder,'mesh.mat'),'mesh') 

%% Create headmodel
cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = zeros(1,5);
scale = 1;
% order follows mesh.tissyelabel , CAUTIOUS!!!! OMg this is not the same order as in the segmentation
cfg.conductivity(find(strcmp(mesh.tissuelabel,'csf'))) = 1.65*scale;
cfg.conductivity(find(strcmp(mesh.tissuelabel,'gray'))) = 0.33*scale;
cfg.conductivity(find(strcmp(mesh.tissuelabel,'scalp'))) = 0.33*scale;
cfg.conductivity(find(strcmp(mesh.tissuelabel,'skull'))) = 0.0042*scale;
cfg.conductivity(find(strcmp(mesh.tissuelabel,'white'))) = 0.126*scale;
cfg.conductivity(find(strcmp(mesh.tissuelabel,'air'))) = 2.5*10^(-14)*scale;

vol        = ft_prepare_headmodel(cfg, mesh);

save(fullfile(save_headmodel_folder,'vol.mat'),'vol','-v7.3') 
%     mkdir(fullfile(Mdrive_folder,subjStr,num2str(NumDir)),folder_keyword);
%         disp('Saving to M drive ...')
%         save(fullfile(Mdrive_folder,subjStr,num2str(NumDir),folder_keyword,'vol.mat'),'vol','-v7.3')

%% Check the electrode alignmnet
% Load the electrodes after digitized
clear fid
chanloc_scan_folder = fullfile(MiM_config.shareFolder,subjStr,'HeadScan','CustomElectrodeLocations.txt');
chanlocs = readtable(chanloc_scan_folder);% Same output text file from getchalocs.
chanlocs.Properties.VariableNames = {'labels','X','Y','Z'};
elec.chanpos(:,1) = [chanlocs.X];
elec.chanpos(:,2) = [chanlocs.Y];
elec.chanpos(:,3) = [chanlocs.Z];
elec.elecpos      = elec.chanpos;
elec.label(:,1)   = [chanlocs.labels]';

% Convert the fiducial position from voxel into CTF 
nas = fiducial_locs.nas;
lpa = fiducial_locs.lpa;
rpa = fiducial_locs.rpa;

vox2head = mri_acpc.transform;

nas = ft_warp_apply(vox2head, nas, 'homogenous');
lpa = ft_warp_apply(vox2head, lpa, 'homogenous');
rpa = ft_warp_apply(vox2head, rpa, 'homogenous');

% create a structure similar to a template set of electrodes
fid.chanpos       = [nas; lpa; rpa];       % CTF head coordinates of fiducials
fid.label         = {'nas','lhj','rhj'};    % use the same labels as those in elec
fid.unit          = 'mm';                  % use the same units as those in mri
fid.elecpos       = fid.chanpos;           % otherwise the electroderealign cannot find elecpos

% alignment
cfg               = [];
cfg.viewmode      = 'surface';
cfg.method        = 'fiducial';
% cfg.method        = 'interactive';%interactive doesn't work well.
cfg.headshape     = vol;
cfg.elec          = elec;                  % the electrodes we want to align
cfg.elecstyle     = {'facecolor','red'};
cfg.headmodelstyle = {'facecolor','b','edgecolor', 'none', 'facealpha', 0.4};
cfg.template      = fid;                   % the template we want to align to
cfg.fiducial      = {'nas', 'lhj', 'rhj'};  % labels of fiducials in fid and in elec
elec_aligned_init = ft_electroderealign(cfg);

figure;
hold on;
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
camlight
ft_plot_sens(elec_aligned_init,'style','.r');
ft_plot_sens(fid,'style','xb');%plot fiducial points
saveas(gcf,fullfile(save_headmodel_folder,'elec_aligned_init.fig'))

%         keyboard
cfg               = [];
cfg.method        = 'project';
cfg.elec          = elec_aligned_init;
cfg.headshape     = vol;
elec_aligned      = ft_electroderealign(cfg);

figure;
hold on;
ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
camlight
%ft_plot_sens(elec_aligned,'style','.r');
ft_plot_sens(elec_aligned ,'style','.r');
ft_plot_sens(fid,'style','xb');%plot fiducial points
saveas(gcf,fullfile(save_headmodel_folder,'elec_aligned.fig'))

save(fullfile(save_headmodel_folder,'elec_aligned.mat'),'elec_aligned') 
%         save(fullfile(Mdrive_folder,subjStr,num2str(NumDir),folder_keyword,'elec_aligned.mat'),'elec_aligned')
end

