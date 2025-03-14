%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/step_to_step_anlz/run_sts_b_epoch_eeg_kin.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
ADD_ALL_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
        SRC_DIR = STUDY_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa');
% subj_names = [SUBJ_PICS{:}];
SUBJ_PICS = {{'H3046','H3047','H3073','H3077','H3092', ...
    'NH3023','NH3025','NH3027',' NH3028', ...
    'NH3051','NH3056','NH3071','NH3082','NH3123'}};
subj_chars = [SUBJ_PICS{:}];
%%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
ica_data_dir = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams'; 
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
ica_data_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies' filesep ica_data_dir]; % JACOB,SAL(02/23/2023)
%%
mri_fpaths = cell(length(subj_chars),1);
mri_fnames = cell(length(subj_chars),1);
%-- masks
mask_regex = '%s_masks_contr.nii.gz';
for subj_i = 1:length(subj_chars)
    mri_fpaths{subj_i} = [PATHS.data_dir filesep DATA_SET filesep subj_chars{subj_i} filesep 'MRI'];
    mri_fnames{subj_i} = sprintf(mask_regex,subj_chars{subj_i});
end
%%
mri_fpath = 'C:\Users\jsalminen\Downloads\Adults\Head\ANTS80-84Years3T_t2w_head.nii.gz';
unzip_out = gunzip(mri_fpath);
% mri_out = ft_read_mri(unzip_out{1});
mri_out = ft_read_mri(unzip_out{1},'datatype','volume','insidestyle','logical');
source3d = ft_checkdata(mri_out, 'datatype', 'source');
%## ATLAS
% atlas_fpath = 'C:\Users\jsalminen\Downloads\Adults\Brain\Atlas\ANTS80-84Years3T_brain_atlas.nii.gz';
atlas_fpath = 'C:\Users\jsalminen\Downloads\Adults\Brain\Atlas\ANTS80-84Years3T_brain_ANTS_IXI_atlas.nii.gz';
unzip_out_atlas = gunzip(mri_fpath);
% atlas_out = ft_read_atlas(unzip_out_atlas{1},'format','afni','unit','mm'); %,'map','prob');
atlas_out = ft_read_mri(unzip_out_atlas{1});
imagesc(atlas_out.parcellation)
%%
atlas_out.coordsys = 'scanras';
mri_out.coordsys = 'scanras';
%%
% parcellation3d = rmfield(source3d, 'anatomy');
% % here we make three ROIs, each is a slab that spans 1/3rd of the whole volume
% % this first representation is "probabilistic"
% parcellation3d.roi1 = false(parcellation3d.dim); parcellation3d.roi1(1:end,1:end,1:64) = true;
% parcellation3d.roi2 = false(parcellation3d.dim); parcellation3d.roi2(1:end,1:end,65:128) = true;
% parcellation3d.roi3 = false(parcellation3d.dim); parcellation3d.roi3(1:end,1:end,129:192) = true;
% % this converts it to an "indexed" representation
% parcellation3d = ft_checkdata(parcellation3d,'type','parcellation','parcellationstyle','indexed');
% cfg = [];
% cfg.anaparameter = [];
% cfg.funparameter = 'tissue';
% cfg.funcolormap = 'jet';
% ft_sourceplot(cfg, parcellation3d);
%%
%- find eeglab on path
if ~ispc
    tmp = strsplit(path,':');
else
    tmp = strsplit(path,';');
end
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM'];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
mri_mni = load(MNI_MRI);
mri_mni = mri_mni.mri;
%%
mri_fpath = 'R:\Ferris-Lab\share\MindInMotion\Data\H2042\MRI\Raw\T1.nii';
mri_out_raw = ft_read_mri(mri_fpath);
%%
mri_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\NH3105\MRI\NH3105_masks_contr.nii.gz';
dunzip_out = gunzip(mri_fpath);
mri_out_sliced = ft_read_mri(unzip_out{1});
%%
mri_fpath='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\NH3105\MRI\mri_acpc_rs.mat';
mri_template_hires='';
tmp = load(mri_fpath);
mri_acpc_rs = tmp.mri_acpc_rs;
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
cfg.template = [];
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
%## write nifti
cfg             = [];
cfg.filename    = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\NH3105\MRI\mri_acpc_rs_mninorm.nii';
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_norm);
%## write nifti
% cfg             = [];
% cfg.filename    = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\NH3105\MRI\mri_acpc_rs.nii';
% cfg.filetype    = 'nifti';
% cfg.parameter   = 'anatomy';
% ft_volumewrite(cfg, mri_acpc_rs);
%%
cfg = [];
%-
cfg.method = 'slice';
cfg.nslices = 1;
cfg.slicerange     = [50,50,50];
cfg.slicedim     = 3;
hold on;
ft_sourceplot(cfg, mri_norm);
ft_sourceplot(cfg, mri_acpc_rs);
ft_sourceplot(cfg, mri_out_raw);
ft_sourceplot(cfg, mri_out_sliced);
ft_sourceplot(cfg, mri_mni);
hold off;
%%
x_dims = 1:mri_norm.dim(1); %50;
y_dims = 1:mri_norm.dim(2);
% z_dims = 1:mri_norm.dim(3);
z_dims = 50;
in_anat = mri_norm.anatomy(x_dims,y_dims,z_dims);
in_mask = ones(size(in_anat))*0.2;
figure();
hold on;
colormap(linspecer)
hf = imagesc(in_anat);
set(hf, 'AlphaData', in_mask)
set(hf, 'AlphaDataMapping', 'scaled')
%##
x_dims = 1:mri_mni.dim(1); %50;
y_dims = 1:mri_mni.dim(2);
% z_dims = 1:mri_norm.dim(3);
z_dims = 50;
in_anat = mri_mni.anatomy(x_dims,y_dims,z_dims);
in_anat = in_anat/255;
in_mask = ones(size(in_anat))*0.2;
colormap('gray');
hf = imagesc(in_anat);
set(hf, 'AlphaData', in_mask)
set(hf, 'AlphaDataMapping', 'scaled')
hold off;
%%
% source_in = mri_out_sliced;
% segmented = source_in;
% segmented.white = source_in.anatomy == 1;
% segmented.gray = source_in.anatomy == 2;
% segmented.csf = (source_in.anatomy == 3 | source_in.anatomy == 8);% csf + ventricles
% segmented.skull = source_in.anatomy == 4;
% segmented.scalp = (source_in.anatomy == 5 | source_in.anatomy == 7);%skin and eye
% segmented.air = source_in.anatomy == 6; 
% segmented = rmfield(segmented,'anatomy');
% seg_i_headreco = ft_datatype_segmentation(segmented,'segmentationstyle','indexed');
%%
cfg              = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
cfg.anaparameter = 'anatomy';
% cfg.funcolormap  = linspecer(6); % distinct color per tissue
cfg.location     = 'center';
cfg.renderer     = 'painter';
cfg.crosshair = 'yes';

ft_sourceplot(cfg,mri_mni);
%##
mni_temp = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm8\templates\T1.nii';
ft_read_mri(mni_temp)
ft_sourceplot(cfg, mri_mni);
%% ===================================================================== %%
subj_name='NH3105';
mri_fpath=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\mri_acpc_rs.mat',subj_name);
mri_template_hires='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
mri_tpm_hires='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_pd_tal_nlin_sym_09a.nii';
spm_default_tpm='M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm12\tpm\TPM.nii';
%-
mni_mri = ft_read_mri(mri_template_hires);
defalt_tpm = ft_read_mri(spm_default_tpm);
mni_tpm = ft_read_mri(mri_tpm_hires);
%-
tmp = load(mri_fpath);
mri_acpc_rs = tmp.mri_acpc_rs;
%-
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
cfg.template = mri_template_hires;
cfg.tpm = mri_tpm_hires;
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
%% PLOT
% GRAY = colormap('gray');
% LIN = linspecer;
sx1 = size(mri_norm.anatomy,1);
sy1 = size(mri_norm.anatomy,2);
sx2 = size(mni_mri.anatomy,1);
sy2 = size(mni_mri.anatomy,2);
diffx = abs(sx1-sx2);
diffy = abs(sy1-sy2);
%-
x_dims_1 = 1:mri_norm.dim(1); %50;
y_dims_1 = 1:mri_norm.dim(2);
% z_dims = 1:mri_norm.dim(3);
z_dims_1 = 50;
x_dims_2 = 1:mni_mri.dim(1); %50;
y_dims_2 = 1:mni_mri.dim(2);
z_dims_2 = 50;
%-
tmp = mri_norm.anatomy(x_dims_1,y_dims_1,z_dims_1);
tmp = tmp/max(max(tmp)); %255;
tmp = padarray(tmp,[diffy,diffx],1,'post');
% tmp = padarray(tmp,diffy,1,'post');
in_anat_1 = cat(3,tmp,tmp,tmp);
in_mask = tmp;
%-
tmp = mni_mri.anatomy(x_dims_2,y_dims_2,z_dims_2);
tmp = tmp/255;
in_anat_2 = cat(3,tmp,tmp,tmp);
%-
green = cat(3,zeros(size(tmp)),ones(size(tmp)),zeros(size(tmp)));

% in_mask = tmp;
%##
size(in_anat)
figure();
hold on;
imagesc(in_anat_1);
imagesc(in_anat_2);
hf = imagesc(green);
set(hf, 'AlphaData', in_mask)
hold off;
