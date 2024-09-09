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
dt = datetime;
dt.Format = 'MMddyyyy_hhmmss';
blaiottac_tpm = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\blaiottaC_tpm\BlaiottaTPM.nii';
mni_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
custom_tpm_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\custom_tpm.nii';
if ~ispc
    mni_fpath = convertPath2UNIX(mni_fpath);
    custom_tpm_fpath = convertPath2UNIX(custom_tpm_fpath);
    blaiottac_tpm = convertPath2UNIX(blaiottac_tpm);
end
ANTS_FPATH = [];
NORM_MRI_METHOD = 'fieldtrip';
DO_PLOT = true;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'eeg_fpath');
addRequired(p,'eeg_fname');
addRequired(p,'mri_fpath');
addRequired(p,'dipfit_fpath');
%## OPTIONAL
%## PARAMETER
addParameter(p,'NORM_MRI_METHOD',NORM_MRI_METHOD);
addParameter(p,'ANTS_FPATH',ANTS_FPATH);
addParameter(p,'DO_PLOT',DO_PLOT);
parse(p,eeg_fpath,eeg_fname,mri_fpath,dipfit_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
NORM_MRI_METHOD = p.Results.NORM_MRI_METHOD;
ANTS_FPATH = p.Results.ANTS_FPATH;
DO_PLOT = p.Results.DO_PLOT;
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
switch NORM_MRI_METHOD
    case 'fieldtrip'
        cfg = [];
        cfg.spmversion = 'spm12';
        cfg.nonlinear = 'yes';
%         cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
        cfg.spmmethod = 'new'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
        cfg.template = mni_fpath;
%         cfg.tpm = custom_tpm_fpath;
%         cfg.tpm = blaiottac_tpm;
        mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);
        dipfit_fem_mnipos = ft_warp_apply(mri_norm.params,ft_warp_apply(mri_norm.initial,dipfit_fem_pos), 'individual2sn');
        dipfit_fem_mni_voxinds = round(ft_warp_apply(pinv(mri_norm.transform), dipfit_fem_mnipos ));
    case 'ants'
%         mri_norm = ft_read_mri(ANTS_FPATH);
%         ants_mri_norm = ft_read_mri([ANTS_FPATH filesep 'antsWarped.nii.gz']);
%         ants_mri_norm = ft_read_mri([ANTS_FPATH filesep 'ants1InverseWarp.nii.gz']);
%         ants_mri_norm = ft_read_mri([ANTS_FPATH filesep 'ants1Warp.nii.gz']);
%         ants_mri_norm = ft_read_mri([ANTS_FPATH filesep 'antsInverseWarped.nii.gz']);
%         ants_transform = load([ANTS_FPATH filesep 'ants0GenericAffine.mat']);
%         R_mat = reshape(ants_transform.AffineTransform_double_3_3,3,4);
%         R_mat = [R_mat; 0 0 0 1];
%         dipfit_fem_mnipos = ft_warp_apply(R_mat,dipfit_fem_pos,'individual2sn');
%         dipfit_fem_mni_voxinds = round(ft_warp_apply(pinv(R_mat), dipfit_fem_mnipos ));
    case 'none'
        fprintf('Not performing normalization\n');
        DO_PLOT = false;
    otherwise
        fprintf('Not performing normalization\n');
        DO_PLOT = false;
end
%% PLOT
if DO_PLOT
%         ft_sourceplot(cfg,mri_norm);
%         fig_i = get(groot,'CurrentFigure');
%         saveas(fig_i,[eeg_fpath filesep sprintf('normalized_mri.fig')]);
%         saveas(fig_i,[eeg_fpath filesep sprintf('normalized_mri.jpg')]);
    %##
    mni_mri = ft_read_mri(mni_fpath);
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
            pre = abs(ctry1-ctry2);
            post = abs(pre-abs(diffy));
%             pre = 0;
%             post = abs(diffy);
            im1 = padarray(im1,[pre,0],1,'pre');
            im1 = padarray(im1,[post,0],1,'post');
        end
        if diffz < 0 && align_z
            pre = abs(ctrz1-ctrz2);
            post = abs(pre-abs(diffz));
%             pre = 0;
%             post = abs(diffz);
            im1 = padarray(im1,[pre,0],1,'pre');
            im1 = padarray(im1,[post,0],1,'post');
        end
        if diffx < 0 && align_x
            pre = abs(ctrx1-ctrx2);
            post = abs(pre-abs(diffy));
%             pre = 0;
%             post = abs(diffx);
            im1 = padarray(im1,[pre,0],1,'pre');
            im1 = padarray(im1,[post,0],1,'post');
        end
        if diffy > 0 && align_y
            pre = abs(ctry1-ctry2);
            post = abs(abs(diffy)-pre);
%             pre = 0;
%             post = abs(diffy);
            im2 = padarray(im2,[pre,0],1,'pre');
            im2 = padarray(im2,[post,0],1,'post');
        end
        if diffz > 0 && align_z
            pre = abs(ctrz1-ctrz2);
            post = abs(abs(diffz)-pre);
%             pre = 0;
%             post = abs(diffz);
            im2 = padarray(im2,[pre,0],1,'pre');
            im2 = padarray(im2,[post,0],1,'post');
        end
        if diffx > 0 && align_x
            pre = abs(ctrx1-ctrx2);
            post = abs(abs(diffx)-pre);
%             pre = 0;
%             post = abs(diffx);
            im2 = padarray(im2,[pre,0],1,'pre');
            im2 = padarray(im2,[post,0],1,'post');
        end
        if align_y || align_z
            im1 = im1';
            im2 = im2';
        end
        in_anat_1 = cat(3,im1,im1,im1);
        in_mask = im1;
        in_anat_2 = cat(3,im2,im2,im2);
        %-
        green = cat(3,zeros(size(im1)),ones(size(im1)),zeros(size(im1)));
        %##
        subplot(2,2,i)
        title(sprintf('MNI_%s=%i   MRI_%s=%i',title_str{i},t_dim_2,title_str{i},t_dim_1));
        hold on;
        imagesc(in_anat_1);
        imagesc(in_anat_2);
        hf = imagesc(green);
        set(hf, 'AlphaData', in_mask)
        hold off;
    end
%             exportgraphics(fig,[mri_fpath filesep 'overlay_validation_plot_oldmnitemp.jpg']);
    exportgraphics(fig,[eeg_fpath filesep sprintf('%s_mninorm_overlay_validation_plot.jpg',dt)]);
end
%% Convert the dipole location to MNI space
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

