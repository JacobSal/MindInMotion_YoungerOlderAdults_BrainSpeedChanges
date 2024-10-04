% This code aims to identify fiducials ACPC and CTF on T1 image
% Originally written by Vanessa Cruz  
% Editted by Chang Liu - 2022-02-17
%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folderd
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'MIM_OA'];

%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS 
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% read mri 
subjStr = 'H2117';%Change here for a different subject
MIM_RDRIVE = 'R:\Ferris-Lab\share\MindInMotion\Data\';
MRIfilepath = fullfile(MIM_RDRIVE,subjStr,'MRI\Raw\T1.nii');
mri = ft_read_mri(MRIfilepath);
savePath = fullfile(MIM_RDRIVE,subjStr,'MRI\Processed_fiducials'); 
mkdir(savePath)
%% check 
check_acpc = isfile(fullfile(savePath,'acpc_fiducials.mat')); 
check_ctf = isfile(fullfile(savePath,'ctf_fiducials.mat'));
if check_acpc == 1
    prompt = 'ACPC fiducials already exists, do you want to continue? Y/N : ';
    str = input(prompt,'s');
    if isempty(str) || str == 'N'
        return
    end
end
if check_ctf == 1
    prompt = 'CTF fiducials already exists, do you want to continue? Y/N : '; %default N
    str = input(prompt,'s');
    if isempty(str) || str == 'N'
        return
    end
end  
%% define acpc
cfg              = [];
cfg.coordsys     = 'acpc'; 
cfg.method       = 'interactive';
cfg.viewmode    = 'ortho';
[mri_acpc] = ft_volumerealign(cfg, mri); %DONT FORGET TO SAVE
acpc_fiducials   = mri_acpc.cfg.fiducial;

save(fullfile(savePath,'mri_acpc.mat'),'mri_acpc');
save(fullfile(savePath,'acpc_fiducials.mat'),'acpc_fiducials');

%% define ctf
cfg             = [];
cfg.coordsys    = 'ctf';
cfg.method      = 'interactive';
cfg.viewmode    = 'ortho';
[mri_ctf] = ft_volumerealign(cfg, mri_acpc); %DONT FORGET TO SAVE
ctf_fiducials = mri_ctf.cfg.fiducial;

save(fullfile(savePath,'mri_ctf.mat'),'mri_ctf');
save(fullfile(savePath,'ctf_fiducials.mat'),'ctf_fiducials');

%% realign the T1 raw image to acpc coordinate
% Save mri_acpc_rs.nii
cfg             = [];
cfg.filename    = fullfile(savePath,strcat(subjStr,'_MRI_acpc'));
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_acpc);

%% reslice data and make the scan isotropic
% While ft_volumerealign does not change the anatomical data, 
% but it adjusts the transformation matrix of the data, 
% ft_volumereslice will change the anatomical data, 
% i.e. it will arrange data in field anatomy according to the coordinate system
cfg             = [];
%         cfg.dim         = [300 300 300]; % default is [256 256 256]. increase the dim so that the head won't get cut-off at the boundary
% There is one participant needs to make the dim larger than the
% default - Chang . But I don't remember who's the participant
cfg.resolution  = [1];
mri_acpc_rs     = ft_volumereslice(cfg,mri_acpc); %can rearrange anatomical data
mri_acpc_rs.coordsys    = 'acpc';

% ---- plot 
cfg     = [];
ft_sourceplot(cfg,mri_acpc_rs);

% Save mri_acpc_rs.nii
cfg             = [];
cfg.filename    = fullfile(savePath,strcat(subjStr,'_MRI_acpc_rs'));
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_acpc_rs);

% Save data as mat
save(fullfile(savePath,'mri_acpc_rs.mat'),'mri_acpc_rs')

%% Generate plots
% plot and save acpc fiducial points
% TO DO: Add figure to put crosshairs the fiducials just marked
% plotting fiducials 
ac = mri_acpc.cfg.fiducial.ac;
pc = mri_acpc.cfg.fiducial.pc;
xz = mri_acpc.cfg.fiducial.xzpoint;

vox2head = mri_acpc.transform; %same transform as mr_realign.transformprig
ac_mri = ft_warp_apply(vox2head, ac);
pc_mri = ft_warp_apply(vox2head, pc);
xz_mri = ft_warp_apply(vox2head, xz);

cfg = [];
cfg.location = ac_mri;
ft_sourceplot(cfg, mri_acpc); %
title('Anterior Commissure location')
saveas(gcf,fullfile(savePath,'fiducial_ac.fig'))

cfg.location = pc_mri;
ft_sourceplot(cfg, mri_acpc); %
title('Posterior Commissure location')
saveas(gcf,fullfile(savePath,'fiducial_pc.fig'))

cfg.location = xz_mri;
ft_sourceplot(cfg, mri_acpc); %
title('XZ Point')
saveas(gcf,fullfile(savePath,'fiducial_xz.fig'))

% plot and save ctf fiducial points
% TO DO: Add figure to put crosshairs the fiducials just marked
% plot fiducials
lpa = mri_ctf.cfg.fiducial.lpa;
rpa = mri_ctf.cfg.fiducial.rpa;
nas = mri_ctf.cfg.fiducial.nas;

vox2head = mri_ctf.transform; %same transform as mr_realign.transformprig
lpa_mri = ft_warp_apply(vox2head, lpa);
rpa_mri = ft_warp_apply(vox2head, rpa);
nas_mri = ft_warp_apply(vox2head, nas);

cfg = [];
cfg.location = lpa_mri;
ft_sourceplot(cfg, mri_ctf); %
title('Left Pre-Auricular Location')
saveas(gcf,fullfile(savePath,'fiducial_lpa.fig'))

cfg.location = rpa_mri;
ft_sourceplot(cfg, mri_ctf); %
title('Right Pre-Auricular Location')
saveas(gcf,fullfile(savePath,'fiducial_rpa.fig'))

cfg.location = nas_mri;
ft_sourceplot(cfg, mri_ctf); %
title('Nasion Location')
saveas(gcf,fullfile(savePath,'fiducial_nas.fig'))
