%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh


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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = fileparts(SCRIPT_DIR);
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%%
%{
%##
for subj_i = 1:length(subj_names)
    subj_name = subj_names{subj_i};
%     mri_fpath=sprintf('R:\\Ferris-Lab\\share\\MindInMotion\\Data\\%s\\MRI\\Processed_fiducials\\',subj_name);
    mri_fpath='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset';
    mri_fpath=[mri_fpath filesep subj_name filesep 'MRI'];
    mri_fname='mri_acpc_rs.mat';
    if exist([mri_fpath filesep mri_fname],'file')
        tmp = load([mri_fpath filesep mri_fname]);
        mri_acpc_rs = tmp.mri_acpc_rs;
        %## write nifti
        cfg             = [];
        cfg.filename    = [mri_fpath filesep 'mri_acpc_rs.nii'];
        cfg.filetype    = 'nifti';
        cfg.parameter   = 'anatomy';
        ft_volumewrite(cfg, mri_acpc_rs);
    else
        fprintf('%s) doesn''t have mri_acpc_rs.mat\n',subj_name);
    end
end
%%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
for subj_i = 1:length(subj_names)
    subj_name = subj_names{subj_i};
    dipfit_fPath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'head_model'];
    dip_struct = par_load(dipfit_fPath,'dipfit_struct.mat');
    coords = zeros(length({dip_struct.dip}),3);
    chan = zeros(length({dip_struct.dip}),1);
    for i = 1:length({dip_struct.dip})
        if ~isempty(dip_struct(i).dip)
            coords(i,:) = dip_struct(i).dip.pos;
            chan(i,:) = i;
        end
    end
    coords(~all(coords,2),:) = [];
    chan(~all(chan,2),:) = [];
    x_coord = coords(:,1);
    y_coord = coords(:,2);
    z_coord = coords(:,3);
    t_in = table(x_coord,y_coord,z_coord);
%     csvwrite([dipfit_fPath filesep 'dip_pos.csv'],coords);
    writetable(t_in,[dipfit_fPath filesep 'dip_pos.csv']);
    t_in = table(chan,x_coord,y_coord,z_coord);
    par_save(t_in,dipfit_fPath,'dip_pos.mat');
end
%}
%% ===================================================================== %%
%{
%## TESTING ANTS DIP POS NORMALIZATION
%## hard define
%- datset name
THRESHOLD_RV_BRAIN = 0.15;
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
parfor (subj_i = 1:length(subj_names),floor(length(subj_names)/2))
% for subj_i = 1:length(subj_names)
    subj_name = subj_names{subj_i};
    mri_path = [DATA_DIR filesep DATA_SET filesep subj_name filesep 'MRI'];
    in_fpath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'head_model'];
    out_fpath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'clean'];
    try
        %-
        norm_pos = readtable([in_fpath filesep 'dip_pos_outf.csv']);
        norm_pos = [norm_pos{:,:}];
        norm_chan = par_load(in_fpath, 'dip_pos.mat');
        %-
        fname = dir([out_fpath filesep '*.set']);
        EEG = pop_loadset('filepath',out_fpath,'filename',fname(1).name);
        %-
    %     trans_mat = load([mri_path filesep 'ants0GenericAffine.mat']);
    %     R_mat = reshape(trans_mat.AffineTransform_double_3_3,3,4);
    %     R_mat = [R_mat; 0 0 0 1];
    %     F_mat = eye(3,3);
    %     F_mat = [F_mat, trans_mat.fixed]; F_mat = [F_mat; 0 0 0 1];
        %-
        ants_mri = ft_read_mri([mri_path filesep 'antsWarped.nii.gz']);
        voxinds = round(ft_warp_apply(pinv(ants_mri.transform), norm_pos ));
        %-
        tmp = load([in_fpath filesep 'dipfit_struct']);
        dipfit_fem = tmp.SAVEVAR;
        %## Reformat Dipfit Structure
        empty_dip_struct = struct('posxyz',[nan(),nan(),nan()],'momxyz',[nan(),nan(),nan()],'rv',nan(),'diffmap',nan(),'sourcepot',nan(),'datapot',nan());
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
                    EEG.dipfit_fem.model(i).rv     = nan();
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
        %##
        for i=1:size(norm_chan,1)
            EEG.dipfit_fem.model(norm_chan.chan(i)).mnipos = norm_pos(i,:);
            EEG.dipfit_fem.model(norm_chan.chan(i)).mni_voxinds = voxinds(i,:);
            EEG.dipfit_fem.model(norm_chan.chan(i)).pos_old = [norm_chan{i,2:end}];
            EEG.dipfit_fem.model(norm_chan.chan(i)).posxyz = norm_pos(i,:);
%             disp(EEG.dipfit_fem.model(norm_chan.chan(i)));
        end
        %## SAVE
        dipfit_fem_norm = EEG.dipfit_fem;
        EEG.dipfit = EEG.dipfit_fem;
        %## Validitiy check
        IC_RV = vertcat(EEG.dipfit.model.rv);
        IC_POSXYZ = vertcat(EEG.dipfit.model.posxyz);
        ICs_RVthreshold_keep = find(IC_RV <= THRESHOLD_RV_BRAIN & all(~isnan(IC_POSXYZ),2));

        EEG = eeg_checkset(EEG);
        pop_saveset(EEG,'filepath',EEG.filepath,'filename',EEG.filename);
        par_save(dipfit_fem_norm,out_fpath,'dipfit_fem_norm_ants.mat');
    catch e
        fprintf('%s\n',getReport(e));
    end
   
end
%}
%% ===================================================================== %%

subj_name='NH3105';
mri_fpath=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\mri_acpc_rs.mat',subj_name);
mri_fpath_ants=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\antsInverseWarped.nii.gz',subj_name);
mri_fpath_antsi=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\antsWarped.nii.gz',subj_name);
mri_template_hires='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
% mri_tpm_hires='M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_pd_tal_nlin_sym_09a.nii';
tpm1 = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_csf_tal_nlin_sym_09a.nii';
tpm2 = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_gm_tal_nlin_sym_09a.nii';
% tpm3 = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_pd_tal_nlin_sym_09a.nii';
tpm4 = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_wm_tal_nlin_sym_09a.nii';
spm_default_tpm='M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm12\tpm\TPM.nii';
%-
% hr_mni_mri = ft_read_mri(mri_template_hires);
tpm1 = ft_read_mri(tpm1);
tpm2 = ft_read_mri(tpm2);
% tpm3 = ft_read_mri(tpm3);
tpm4 = ft_read_mri(tpm4);
%-
defalt_tpm = ft_read_mri(spm_default_tpm);
tmp_anat = defalt_tpm.anatomy;
%-
tmp = load(mri_fpath);
mri_acpc_rs = tmp.mri_acpc_rs;
%-
ants_norm_mri = ft_read_mri(mri_fpath_ants);
%%
%{
%## CUSTOM TPM
custom_tpm = tpm1;
anat = cat(4,tpm2.anatomy,tpm4.anatomy,tpm1.anatomy);
% custom_tpm.gray = tpm2.anatomy;
% custom_tpm.white = tpm4.anatomy;
% custom_tpm.csf = tpm1.anatomy;
custom_tpm.anatomy = cat(4,tpm2.anatomy,tpm4.anatomy,tpm1.anatomy);
tmp_anat_cust = custom_tpm.anatomy;
%-
max(max(max(tpm2.anatomy)))
max(max(max((tmp_anat(:,:,:,1)))))
%##
mni_mri = ft_read_mri(mri_template_hires);
mni_mri.coordsys = 'acpc';
% cfg = [];
% cfg.spmversion = 'spm12';
% cfg.spmmethod = 'new';
% cfg.tpm = custom_tpm_fpath;
% cfg.coordsys = 'acpc';
% cfg.output = {'tpm'};
% segment_tpm = ft_volumesegment(cfg,mni_mri);
% cfg.output = {'gray','white','csf','skull','scalp'};
% segments_out = ft_volumesegment(cfg,ssegment_tpm);
% segments_out = ft_volumesegment(cfg,mni_mri);
%## write nifti
custom_tpm_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\custom_tpm.nii';
cfg             = [];
cfg.filename    = custom_tpm_fpath;
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, custom_tpm);
%%
cfg             = [];
cfg.coordsys    = 'acpc'; % add coordsys to mri struct so that the plot can be centered
cfg.method      = 'interactive';
[mri_new_acpc] = ft_volumerealign(cfg, mri_acpc_rs); %originaly thought not needed but needed to define coorsys for segmentation, if not already defined
cfg = [];
mri_new_acpc = ft_volumereslice(cfg,mri_new_acpc);
%}
%% ===================================================================== %%
%{
cfg              = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
cfg.anaparameter = 'anatomy';
% cfg.funcolormap  = linspecer(6); % distinct color per tissue
cfg.location     = 'center';
cfg.renderer     = 'painter';
cfg.crosshair   = 'yes';
%-
% ft_sourceplot(cfg,tpm1);
% ft_sourceplot(cfg,tpm2);
% ft_sourceplot(cfg,tpm3);
% ft_sourceplot(cfg,tpm4);
%-
defalt_tpm.anatomy = tmp_anat(:,:,:,1); %1, gm; 2, wm; 3, csf; 4, skull; 5, skin; 6, air
ft_sourceplot(cfg,defalt_tpm);
%-
% segments_out.anatomy = segments_out.csf; %1, gm; 2, wm; 3, csf; 4, skull; 5, skin; 6, air
% ft_sourceplot(cfg,segments_out);
%-
custom_tpm.anatomy = tmp_anat_cust(:,:,:,1); %1, gm; 2, wm; 3, csf; 4, skull; 5, skin; 6, air
ft_sourceplot(cfg,custom_tpm);
%-
ft_sourceplot(cfg,mri_new_acpc);
ft_sourceplot(cfg,mri_acpc_rs);
%}
%% ===================================================================== %%
% out1 = spm_load_priors8(custom_tpm_fpath)
% out2 = spm_load_priors8(custom_tpm_fpath)
%-
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
% cfg.template = mri_template_hires;
% cfg.tpm = custom_tpm_fpath;
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
% mri_norm = ft_volumenormalise(cfg,mri_new_acpc);
%% PLOT
% GRAY = colormap('gray');
% LIN = linspecer;
mri_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm8\templates\T1.nii';
mni_mri = ft_read_mri(mri_fpath);
% mri_norm = mri_new_acpc;
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
% size(in_anat_1)
figure();
title('Normalized MRI & MNI alignment');
hold on;
imagesc(in_anat_1);
imagesc(in_anat_2);
hf = imagesc(green);
set(hf, 'AlphaData', in_mask)
hold off;
%% ===================================================================== %%
subj_name = 'NH3066';
mri_fpath=sprintf('R:\\Ferris-Lab\\share\\MindInMotion\\Data\\%s\\MRI\\Processed_fiducials\\',subj_name);
mri_fname='mri_acpc_rs.mat';
tmp = load([mri_fpath filesep mri_fname]);
mri_acpc_rs = tmp.mri_acpc_rs;
do_mni_old = true;
if do_mni_old
    mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm12\toolbox\OldNorm\T1.nii';
%     mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm8\templates\T1.nii';
    custom_tpm_fpath = [];
else
    mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
    custom_tpm_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\custom_tpm.nii';
end
mni_mri = ft_read_mri(mni_fpath);
% mni_mri = ft_read_mri(mni_fpath);
% mni_mri = ft_read_mri(mri_template_hires);
%##
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
% cfg.template = mni_fpath;
% cfg.tpm = custom_tpm_fpath;
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
%## write nifti
cfg             = [];
cfg.filename    = [mri_fpath filesep 'mri_acpc_rs_mninorm.nii'];
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_norm);
%%
% mri_norm = mri_new_acpc;
sx1 = size(mri_norm.anatomy,1);
sy1 = size(mri_norm.anatomy,2);
sx2 = size(mni_mri.anatomy,1);
sy2 = size(mni_mri.anatomy,2);
diffx = sx1-sx2;
diffy = sy1-sy2;
%-
x_dims_1 = 1:mri_norm.dim(1); %50;
y_dims_1 = 1:mri_norm.dim(2);
% z_dims = 1:mri_norm.dim(3);
z_dims_1 = floor(mri_norm.dim(3)/2);
x_dims_2 = 1:mni_mri.dim(1); %50;
y_dims_2 = 1:mni_mri.dim(2);
z_dims_2 = floor(mni_mri.dim(3)/2);
%-
tmp = mri_norm.anatomy(x_dims_1,y_dims_1,z_dims_1);
tmp = tmp/max(max(tmp)); %255;
if diffx < 0
    tmp = padarray(tmp,[abs(diffx),0],1,'post');
end
if diffy < 0 
    tmp = padarray(tmp,[0,abs(diffy)],1,'post');
end
% tmp = padarray(tmp,diffy,1,'post');
in_anat_1 = cat(3,tmp,tmp,tmp);
in_mask = tmp;
%-
tmp = mni_mri.anatomy(x_dims_2,y_dims_2,z_dims_2);
tmp = tmp/255;
if diffx > 0
    tmp = padarray(tmp,[abs(diffx),0],1,'post');
end
if diffy > 0 
    tmp = padarray(tmp,[0,abs(diffy)],1,'post');
end
in_anat_2 = cat(3,tmp,tmp,tmp);
%-
green = cat(3,zeros(size(tmp)),ones(size(tmp)),zeros(size(tmp)));
% in_mask = tmp;
%##
% size(in_anat_1)
fig = figure();
title('Normalized MRI & MNI alignment');
hold on;
imagesc(in_anat_1);
imagesc(in_anat_2);
hf = imagesc(green);
set(hf, 'AlphaData', in_mask)
hold off;
if do_mni_old
    exportgraphics(fig,[mri_fpath filesep 'overlay_validation_plot_oldmnitemp.jpg']);
else
    exportgraphics(fig,[mri_fpath filesep 'overlay_validation_plot_newmnitemp.jpg']);
end
%% ===================================================================== %%
subj_name='H1020';
mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
% mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm12\toolbox\OldNorm\T1.nii';
save_dir=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\',subj_name);
mri_fpath_ants=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\antsInverseWarped.nii.gz',subj_name);
mri_fpath_antsi=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\antsWarped.nii.gz',subj_name);
ants_norm_mri = ft_read_mri(mri_fpath_ants);
mri_norm = ants_norm_mri;
mni_mri = ft_read_mri(mni_fpath);
%-
% cfg = [];
% cfg.method = 'interactive';
% cfg.coordsys = 'acpc';
% mri_norm = ft_volumerealign(cfg,mri_norm);
%-
cfg              = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
cfg.anaparameter = 'anatomy';
% cfg.funcolormap  = linspecer(6); % distinct color per tissue
cfg.location     = 'center';
cfg.renderer     = 'painter';
cfg.crosshair   = 'yes';
ft_sourceplot(cfg,mri_norm)
%-
% cfg              = [];
% cfg.method = 'ortho';
% cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
% cfg.anaparameter = 'anatomy';
% % cfg.funcolormap  = linspecer(6); % distinct color per tissue
% cfg.location     = 'center';
% cfg.renderer     = 'painter';
% cfg.crosshair   = 'yes';
% ft_sourceplot(cfg,mni_mri)
%%
subj_name='H1020';
mri_fpath=sprintf('M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\%s\\MRI\\mri_acpc_rs.mat',subj_name);
tmp = load(mri_fpath);
mri_acpc_rs = tmp.mri_acpc_rs;
%-
mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
% mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\external\spm12\toolbox\OldNorm\T1.nii';
%-
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
cfg.template = mni_fpath;
% cfg.tpm = custom_tpm_fpath;
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
%-
cfg              = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
cfg.anaparameter = 'anatomy';
% cfg.funcolormap  = linspecer(6); % distinct color per tissue
cfg.location     = 'center';
cfg.renderer     = 'painter';
cfg.crosshair   = 'yes';
ft_sourceplot(cfg,mri_norm)
%%
% NOTE: AWFUL CODE DOESNT WORK
%{
%##
sx1 = size(mri_norm.anatomy,1);
sy1 = size(mri_norm.anatomy,2);
sz1 = size(mri_norm.anatomy,3);
sx2 = size(mni_mri.anatomy,1);
sy2 = size(mni_mri.anatomy,2);
sz2 = size(mni_mri.anatomy,3);
ctrx1 = floor(sx1/2);
ctry1 = floor(sy1/2);
ctrz1 = floor(sz1/2);
ctrx2 = ceil(sx2/2);
ctry2 = ceil(sy2/2);
ctrz2 = ceil(sz2/2);
diffx = sx1-sx2;
diffy = sy1-sy2;
diffz = sz1-sz2;
%##
% size(in_anat_1)
fig = figure();
for i = 1:3
    x_dims_1 = 1:mri_norm.dim(1); %50;
    y_dims_1 = 1:mri_norm.dim(2);
    z_dims_1 = 1:mri_norm.dim(3);
    x_dims_2 = 1:mni_mri.dim(1); %50;
    y_dims_2 = 1:mni_mri.dim(2);
    z_dims_2 = 1:mni_mri.dim(3);
    align_x = false; 
    align_y = false;
    align_z = false;
    title_str = {'x','y','z'};
    %-
    switch i
        case 1
            align_y = true;
            align_z = true;
            x_dims_1 = floor(mni_mri.dim(1)/2);
            x_dims_2 = floor(mni_mri.dim(1)/2);
            im1 = squeeze(mri_norm.anatomy(x_dims_1,y_dims_1,z_dims_1));
            im1 = im1/max(max(im1));
            %-
            im2 = squeeze(mni_mri.anatomy(x_dims_2,y_dims_2,z_dims_2));
            im2 = im2/255;
            %-
            t_dim_2 = x_dims_2;
            t_dim_1 = x_dims_1;

        case 2
            align_x = true;
            align_z = true;
            y_dims_1 = floor(mni_mri.dim(2)/2);
            y_dims_2 = floor(mni_mri.dim(2)/2);
            im1 = squeeze(mri_norm.anatomy(x_dims_1,y_dims_1,z_dims_1));
            im1 = im1/max(max(im1));
            %-
            im2 = squeeze(mni_mri.anatomy(x_dims_2,y_dims_2,z_dims_2));
            im2 = im2/255;
            %-
            t_dim_2 = y_dims_2;
            t_dim_1 = y_dims_1;
        case 3
            align_x = true; 
            align_y = true;
            z_dims_1 = floor(mni_mri.dim(3)/2);
            z_dims_2 = floor(mni_mri.dim(3)/2);
            im1 = mri_norm.anatomy(x_dims_1,y_dims_1,z_dims_1);
            im1 = im1/max(max(im1));
            %-
            im2 = squeeze(mni_mri.anatomy(x_dims_2,y_dims_2,z_dims_2));
            im2 = im2/255;
            %-
            t_dim_2 = z_dims_2;
            t_dim_1 = z_dims_1;
    end
    %-
    if diffy < 0 && align_y
        disp(1)
%         pre = abs(ctry1-ctry2);
%         post = abs(pre-abs(diffy));
            pre = 0;
            post = abs(diffy);
        if ~align_x
            im1 = padarray(im1,[pre,0],1,'pre');
            im1 = padarray(im1,[post,0],1,'post');
        else
            im1 = padarray(im1,[0,pre],1,'pre');
            im1 = padarray(im1,[0,post],1,'post');
        end
    end
    if diffz < 0 && align_z
        disp(2)
%         pre = abs(ctrz1-ctrz2);
%         post = abs(pre-abs(diffz));
            pre = 0;
            post = abs(diffz);
        if ~align_x
            im1 = padarray(im1,[0,pre],1,'pre');
            im1 = padarray(im1,[0,post],1,'post');
        else
            im1 = padarray(im1,[0,pre],1,'pre');
            im1 = padarray(im1,[0,post],1,'post');
        end
    end
    if diffx < 0 && align_x
        disp(3)
%         pre = abs(ctrx1-ctrx2);
%         post = abs(pre-abs(diffy));
            pre = 0;
            post = abs(diffx);
        if ~align_x
            im1 = padarray(im1,[0,pre],1,'pre');
            im1 = padarray(im1,[0,post],1,'post');
        else
            im1 = padarray(im1,[0,pre],1,'pre');
            im1 = padarray(im1,[0,post],1,'post');
        end
    end
    if diffy > 0 && align_y
        disp(4)
%         pre = abs(ctry1-ctry2);
%         post = abs(abs(diffy)-pre);
            pre = 0;
            post = abs(diffy);
        if ~align_x
            im2 = padarray(im2,[pre,0],1,'pre');
            im2 = padarray(im2,[post,0],1,'post');
        else
            im2 = padarray(im2,[0,pre],1,'pre');
            im2 = padarray(im2,[0,post],1,'post');
        end
    end
    if diffz > 0 &&  align_z
        disp(5)
%         pre = abs(ctrz1-ctrz2);
%         post = abs(abs(diffz)-pre);
            pre = 0;
            post = abs(diffz);
        if ~align_x
            im2 = padarray(im2,[0,pre],1,'pre');
            im2 = padarray(im2,[0,post],1,'post');
        else
            im2 = padarray(im2,[0,pre],1,'pre');
            im2 = padarray(im2,[0,post],1,'post');
        end
    end
    if diffx > 0 && align_x
        disp(6)
%         pre = abs(ctrx1-ctrx2);
%         post = abs(abs(diffx)-pre);
            pre = 0;
            post = abs(diffx);
        if ~align_x
            im2 = padarray(im2,[0,pre],1,'pre');
            im2 = padarray(im2,[0,post],1,'post');
        else
            im2 = padarray(im2,[0,pre],1,'pre');
            im2 = padarray(im2,[0,post],1,'post');
        end
    end
    if align_y && align_z
        im1 = im1';
        im2 = im2';
    end
    in_anat_1 = cat(3,im1,im1,im1);
    in_mask = im1;
    in_anat_2 = cat(3,im2,im2,im2);
    %-
    green = cat(3,zeros(size(im2)),ones(size(im2)),zeros(size(im2)));
    %##
    subplot(2,2,i)
    title(sprintf('MNI_%s=%i   MRI_%s=%i',title_str{i},t_dim_2,title_str{i},t_dim_1));
    figure;
    hold on;
    imagesc(in_anat_1);
    imagesc(in_anat_2);
    hf = imagesc(green);
    set(hf, 'AlphaData', in_mask)
    hold off;
end
%             exportgraphics(fig,[mri_fpath filesep 'overlay_validation_plot_oldmnitemp.jpg']);
exportgraphics(fig,[save_dir filesep sprintf('%s_mninorm_overlay_validation_plot.jpg',dt)]);

%}
