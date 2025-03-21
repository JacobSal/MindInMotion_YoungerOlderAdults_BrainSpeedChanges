%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 

%- run script
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/prep_scripts/run_prep_g_ants_update_eeg.sh

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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
subj_chars = [SUBJ_PICS{:}];
% SUBJ_PICS = {{'NH3105'}}; %'H2090','H2095','H2111','H2117','H3029','H3034','H3039','H3053'}};
% subj_chars = [SUBJ_PICS{:}];
%% ===================================================================== %%
%## MRI TEMPLATE FOR SOURCE DEPTH
HIRES_TEMPLATE = [PATHS.data_dir filesep '_resources' filesep 'mni_icbm152_nlin_sym_09a' filesep 'mni_icbm152_t1_tal_nlin_sym_09a.nii'];
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
MNI_VOL = [fpath filesep fname '_dipplotvol.mat'];
vol = load(MNI_VOL);
try
    vol = vol.vol;
catch
    vol = vol.mesh;
end

%## DIPOLE PARAMS
THRESHOLD_RV_BRAIN = 0.15;
DEPTH_ERROR = 15; % +/-15mm (30 mm sphere)
%% (PATHS) ============================================================= %%
DATA_SET = 'MIM_dataset';
STUDY_DNAME = '02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej';
%## soft define
ica_data_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies' filesep STUDY_DNAME]; % JACOB,SAL(02/23/2023)
%%
parfor subj_i = 1:length(subj_chars)
% for subj_i = 1:length(subj_chars)
    subj_name = subj_chars{subj_i};
    mri_path = [PATHS.data_dir filesep DATA_SET filesep subj_name filesep 'MRI'];
    in_fpath = [ica_data_dir filesep subj_name filesep 'head_model'];
    out_fpath = [ica_data_dir filesep subj_name filesep 'clean'];
    try
        %- load transformed dipole pos & convert to RAS from LPS
        norm_pos = readtable([in_fpath filesep 'dip_pos_outf.csv']);
        norm_pos = [norm_pos{:,:}];
        for i = 1:size(norm_pos,1)
            norm_pos(i,:) = [-norm_pos(i,1),-norm_pos(i,2),norm_pos(i,3)];
        end
        % (05/09/2024) JS, flipping x,y coords for dipole to recorrect the
        % LPS to RAS conversion made before the transformation is made. 
        %- load channel mapping file
        norm_chan = par_load(in_fpath, 'dip_pos.mat');
        %- load EEG set file
        fname = dir([out_fpath filesep '*.set']);
        EEG = pop_loadset('filepath',out_fpath,'filename',fname(1).name);
        %-
    %     trans_mat = load([mri_path filesep 'ants0GenericAffine.mat']);
    %     R_mat = reshape(trans_mat.AffineTransform_double_3_3,3,4);
    %     R_mat = [R_mat; 0 0 0 1];
    %     F_mat = eye(3,3);
    %     F_mat = [F_mat, trans_mat.fixed]; F_mat = [F_mat; 0 0 0 1];
        %- convert dipole coordinates (cartesian space in mm) to voxel
        ants_mri = ft_read_mri([mri_path filesep 'antsWarped.nii.gz']);
        voxinds = round(ft_warp_apply(pinv(ants_mri.transform), norm_pos ));
        %- load the original dipfit structure for editing
        tmp = load([in_fpath filesep 'dipfit_struct']);
        dipfit_fem = tmp.SAVEVAR;
        %## Reformat Dipfit Structure
        empty_dip_struct = struct('posxyz',[nan(),nan(),nan()],...
            'momxyz',[nan(),nan(),nan()],...
            'rv',nan(),...
            'diffmap',nan(),...
            'sourcepot',nan(),...
            'datapot',nan(),...
            'mnipos',[],...
            'mni_voxinds',[],...
            'pos_old',[]);
        EEG.dipfit_fem = [];
        EEG.dipfit_fem.model = empty_dip_struct;
        dipfit_fem_pos = zeros(length([dipfit_fem.component]),3);
        %## LOOP THROUGH NON-NORMAL DIPOLES AND ADD TO STRUCT
        for i=1:length([dipfit_fem.component])
            %- if dipole was calculated
            if ~isempty(dipfit_fem(i).dip)
                %- update dipfit_fem struct in EEG with original fem fits
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
        %## LOOP THROUGH NORMALIZED DIPOLES AND ADD TO STRUCT
        for i=1:size(norm_chan,1)
            %- get dipole depth and reject if it isn't in brain volume
            depth = ft_sourcedepth(norm_pos(i,:), vol);
            if depth <= DEPTH_ERROR
                EEG.dipfit_fem.model(norm_chan.chan(i)).mnipos = norm_pos(i,:);
                EEG.dipfit_fem.model(norm_chan.chan(i)).mni_voxinds = voxinds(i,:);
                EEG.dipfit_fem.model(norm_chan.chan(i)).pos_old = [norm_chan{i,2:end}];
                EEG.dipfit_fem.model(norm_chan.chan(i)).posxyz = norm_pos(i,:);
            else
                fprintf('%s) dip [%0.2g,%0.2g,%0.2g] outside brain volume: %0.2f...\n',subj_name,dipfit_fem(norm_chan.chan(i)).dip.pos,depth)
                EEG.dipfit_fem.model(norm_chan.chan(i)) = empty_dip_struct;
            end
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
