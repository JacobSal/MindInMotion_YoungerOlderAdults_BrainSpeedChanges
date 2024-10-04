%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 

%- run script
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/run_g_ants_update_eeg.sh

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
%% ADD STUDY, SRC, AND WORKSPACE PATHS
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
subj_names = [SUBJ_PICS{:}];
%- datset name
DATA_SET = 'MIM_dataset';
study_dir_name = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
icadata_dir = [studies_fpath filesep study_dir_name];
%%
HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
if ~ispc
    HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
else
    HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
end
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
hires_mesh = [fpath filesep fname '_dipplotvol.mat'];
hires_mri = [fpath filesep fname '_dipplotmri.mat'];
mri = load(hires_mri);
MNI_MRI = mri.mri;
MNI_VOL = hires_mesh;
%% (DEBUG) DIPOLE CONFIRMATION 
subj_i = 15;
dip_i = 1;
disp(subj_names{subj_i});
empty_dip_struct = struct('posxyz',[nan(),nan(),nan()],...
    'momxyz',[nan(),nan(),nan()],...
    'rv',nan(),...
    'diffmap',nan(),...
    'sourcepot',nan(),...
    'datapot',nan(),...
    'mnipos',[],...
    'mni_voxinds',[],...
    'pos_old',[]);
%-
subj_name = subj_names{subj_i};
mri_path = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep subj_name filesep 'MRI'];
in_fpath = [icadata_dir filesep subj_name filesep 'head_model'];
out_fpath = [icadata_dir filesep subj_name filesep 'clean'];
%-
orig_mri = ft_read_mri([mri_path filesep 'mri_acpc_rs.nii']);
% orig_mri = ft_read_mri([mri_path filesep 'mri_acpc.mat']);
ants_mri = ft_read_mri([mri_path filesep 'antsWarped.nii.gz']);
%-
norm_pos = readtable([in_fpath filesep 'dip_pos_outf.csv']);
norm_pos = [norm_pos{:,:}];
% nonras_norm_pos = norm_pos;
for i = 1:size(norm_pos,1)
    norm_pos(i,:) = [-norm_pos(i,1),-norm_pos(i,2),norm_pos(i,3)];
end
%-
norm_pos_old = readtable([in_fpath filesep 'dip_pos_outf_oldtrans.csv']);
norm_pos_old = [norm_pos_old{:,:}];
% nonras_norm_pos = norm_pos_old;
% for i = 1:size(norm_pos_old,1)
%     norm_pos_old(i,:) = [norm_pos_old(i,1),norm_pos_old(i,2),norm_pos_old(i,3)];
% end
%-
tmp = load([in_fpath filesep 'dipfit_struct']);
dipfit_fem = tmp.SAVEVAR;
norm_chan = par_load(in_fpath, 'dip_pos.mat');
%-
DIPPLOT_STRUCT = struct('rvrange',[0,30],... % this is a value from 0 to 100 (e.g., rv = 0.15 is 15)
        'summary','off',...
        'mri',orig_mri,...
        'coordformat','MNI',...
        'transform',[],...
        'image','mri',...
        'plot','on',...
        'color',{{[0,0,1]}},...
        'view',[1,1,1],...
        'mesh','off',...
        'meshdata',MNI_VOL,...
        'axistight','off',... % display the closest MRI slice to distribution
        'gui','off',...
        'num','off',...
        'cornermri','on',...
        'drawedges','off',...
        'projimg','off',...
        'projlines','off',...
        'projwidth',1,...
        'projcol',{{[0,0,1]}},...
        'dipolesize',30,...
        'dipolelength',0,...
        'pointout','off',...
        'sphere',1,...
        'spheres','off',...
        'normlen','off',...
        'dipnames',{{}},...
        'holdon','on',...
        'camera','auto',...
        'density','off');
%%
% std_topoplot_CL(STUDY,cl_i,'together');
modout = loadmodout15(out_fpath);
fname = dir([out_fpath filesep '*.set']);
[~,EEG,~] = eeglab_loadICA(fname(1).name,out_fpath);
% X = zeros([ length(comps) size(tmp) ]);
% X = zeros(1,1,1,1);
%%
%-
inds = randi(size(norm_chan,1),5,1);
chanlocs = EEG.chanlocs(EEG.icachansind);
topo_struct = struct('grid',[],...
    'x',[],...
    'y',[]);
option = 'none';
for i = 1:size(EEG.icawinv,2)
    [~, grid, ~, Xi, Yi] = topoplot(EEG.icawinv(:,i), chanlocs, ...
                                                      'verbose', 'off',...
                                                       'electrodes', 'on' ,'style','both',...
                                                       'plotrad',0.55,'intrad',0.55,...
                                                       'noplot', 'on', 'chaninfo', EEG.chaninfo);
    topo_struct(i).grid = grid;
    topo_struct(i).x = Xi;
    topo_struct(i).y = Yi;
    % X(i,:,:,:) = grid;
    %-
end
%##
inds = 14;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
% X = squeeze(X);
fig = figure;
set(gca,AXES_DEFAULT_PROPS{:})
hold on;
for i = 1:length(inds)
    dip_i = inds(i);
    if dipfit_fem(norm_chan.chan(dip_i)).dip.rv < 0.15
        figure(fig);
        sbplot(ceil(length(inds)/5),5,i)
        Xi = topo_struct(i).x;
        Yi = topo_struct(i).y;
        toporeplot(topo_struct(norm_chan.chan(dip_i)).grid, 'style', 'both',...
            'plotrad',0.5,'intrad',0.5, 'verbose', 'off','xsurface', Xi, 'ysurface', Yi );
        colormap(linspecer); 
    end
end
hold off;
exportgraphics(fig,[mri_path filesep sprintf('scalpmaps_test.jpg')],'Resolution',200);
%%
% leadf = par_load(mri_path, 'leadfield_fem.mat');
% hmod_tr = par_load(mri_path, 'headmodel_fem_tr.mat');
%% NON-NORM MRI, ORIGNAL DIP
%- reformat to eeglab style...
% inds = randi(size(norm_chan,1),5,1);
inds = 14;
mri = [];
mri.dim = orig_mri.dim;
mri.xgrid = 1:orig_mri.dim(1);
mri.ygrid = 1:orig_mri.dim(2);
mri.zgrid = 1:orig_mri.dim(3);
mri.anatomy = uint8(orig_mri.anatomy*0.3); %uint8(tmp.anatomy); %uint8(tmp.anatomy);
mri.transform = orig_mri.transform;
mri.hdr = orig_mri.hdr;
DIPPLOT_STRUCT.mri = mri;
DIPPLOT_STRUCT.meshdata = [mri_path filesep 'vol.mat'];
DIPPLOT_STRUCT.dipolelength = 1;
DIPPLOT_STRUCT.normlen= 'on';
fig = figure;
hold on;
for i = 1:length(inds)
    dip_i = inds(i);
    if dipfit_fem(norm_chan.chan(dip_i)).dip.rv < 0.15
        dip_in = struct('posxyz',dipfit_fem(norm_chan.chan(dip_i)).dip.pos,...
            'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
            'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
        dipplot(dip_in,DIPPLOT_STRUCT);
    end
end
%-
% i = 3;
% dip_i = inds(i);
% dip_in = struct('posxyz',dipfit_fem(norm_chan.chan(dip_i)).dip.pos,...
%     'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
%     'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
% dipplot(dip_in,DIPPLOT_STRUCT);
hold off;
exportgraphics(fig,[mri_path filesep sprintf('dipplot_origdip_top.tiff')],'Resolution',1000);
view([45,0,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_origdip_coronal.tiff')],'Resolution',1000);
view([0,-45,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_origdip_sagittal.tiff')],'Resolution',1000);

%% NORM MRI, NORM DIP
ants_mri = ft_read_mri([mri_path filesep 'antsWarped.nii.gz']);
%- reformat to eeglab style...
mri = [];
mri.coordsys = 'acpc';
mri.dim = ants_mri.dim;
mri.xgrid = 1:ants_mri.dim(1);
mri.ygrid = 1:ants_mri.dim(2);
mri.zgrid = 1:ants_mri.dim(3);
mri.anatomy = uint8(ants_mri.anatomy*0.5);
mri.transform = ants_mri.transform;
mri.hdr = ants_mri.hdr;
%-
if ~exist([mri_path filesep 'antsWarped_vol.mat'],'file')
    cfg             = [];
    cfg.output      = {'scalp', 'skull', 'brain'};
    cfg.spmmethod   = 'old';
    cfg.spmversion  = 'spm12';
    segmentation    = ft_volumesegment(cfg, mri);
    cfg             = [];
    cfg.tissue      = {'scalp', 'skull', 'brain'};
    cfg.numvertices = [500,1000,1500];
    mesh            = ft_prepare_mesh(cfg, segmentation);
    cfg             = [];
    cfg.method      = 'bemcp';
    ftvol_tmp       = ft_prepare_headmodel(cfg, mesh);
    %- convert to eeglab format....
    vol = [];
    vol.ft_vol = ftvol_tmp;
    vol.bnd = ftvol_tmp.bnd(3:-1:1);
    vol.cond = ftvol_tmp.cond;
    tmp_mat = ftvol_tmp.mat;
    tmp_mat = padarray(tmp_mat,[2500,0],0,'post');
    tmp_mat = sparse(tmp_mat);
    vol.mat = tmp_mat;
    vol.type = 'bemcp';
    save([mri_path filesep 'antsWarped_vol.mat'],'vol');
end
DIPPLOT_STRUCT.meshdata = [mri_path filesep 'antsWarped_vol.mat'];
%-
DIPPLOT_STRUCT.mri = mri;
DIPPLOT_STRUCT.meshdata = [mri_path filesep 'antsWarped_vol.mat'];
fig = figure;
hold on;

for i = 1:length(inds)
    dip_i = inds(i);
    if dipfit_fem(norm_chan.chan(dip_i)).dip.rv < 0.15
        DIPPLOT_STRUCT.color = {[0,0,1]};
        dip_in = struct('posxyz',norm_pos(dip_i,:),...
            'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
            'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
        dipplot(dip_in,DIPPLOT_STRUCT);
        DIPPLOT_STRUCT.color = {[1,0,0]};
        dip_in = struct('posxyz',norm_pos_old(dip_i,:),...
            'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
            'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
        dipplot(dip_in,DIPPLOT_STRUCT);
    end
end
hold off;
exportgraphics(fig,[mri_path filesep sprintf('dipplot_normmrinormdip_top.jpg')],'Resolution',200);
view([45,0,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_normmrinormdip_coronal.jpg')],'Resolution',200);
view([0,-45,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_normmrinormdip_sagittal.jpg')],'Resolution',200);

%% MNI MRI, NORM DIP
HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
if ~ispc
    HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
else
    HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
end
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
hires_mesh = [fpath filesep fname '_dipplotvol.mat'];
hires_mri = [fpath filesep fname '_dipplotmri.mat'];
mri = load(hires_mri);
mri = mri.mri;
vol = hires_mesh;
DIPPLOT_STRUCT.mri = mri;
DIPPLOT_STRUCT.meshdata = vol;
fig = figure;
hold on;
DIPPLOT_STRUCT.color = {[0,0,1]};
for i = 1:length(inds)
    dip_i = inds(i);
    if dipfit_fem(norm_chan.chan(dip_i)).dip.rv < 0.15
        DIPPLOT_STRUCT.color = {[0,0,1]};
        dip_in = struct('posxyz',norm_pos(dip_i,:),...
            'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
            'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
        dipplot(dip_in,DIPPLOT_STRUCT);
        DIPPLOT_STRUCT.color = {[1,0,0]};
        dip_in = struct('posxyz',norm_pos_old(dip_i,:),...
            'momxyz',reshape(dipfit_fem(norm_chan.chan(dip_i)).dip.mom, 3, length(dipfit_fem(norm_chan.chan(dip_i)).dip.mom)/3)',...
            'rv',dipfit_fem(norm_chan.chan(dip_i)).dip.rv);
        dipplot(dip_in,DIPPLOT_STRUCT);
    end
end
hold off;
exportgraphics(fig,[mri_path filesep sprintf('dipplot_mnimri_top.jpg')],'Resolution',200);
view([45,0,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_mnimri_coronal.jpg')],'Resolution',200);
view([0,-45,0])
exportgraphics(fig,[mri_path filesep sprintf('dipplot_mnimri_sagittal.jpg')],'Resolution',200);