%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR%#ok<GVMIS>
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
%## Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%- datset name
DATA_SET = 'MIM_dataset';
study_dir_name = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
load_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET];
%##

% spm download: https://www.fil.ion.ucl.ac.uk/spm/software/download/spmreg.php
% cat12 download: https://andysbrainbook.readthedocs.io/en/latest/CAT12/CAT12_01_DownloadInstall.html
% see.
% https://andysbrainbook.readthedocs.io/en/latest/CAT12/CAT12_03_Preprocessing.html
% for tutorial
%%
% addpath('C:\Users\jsalminen\Documents\GitHub\spm12\spm12')
% spm fmri
% gunzip('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\H1018\MRI\antsWarped.nii.gz');
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
mri = mri.mri;
vol = hires_mesh;
%%
subj_chars = [SUBJ_PICS{:}];
files = dir([load_dir filesep '*' filesep 'MRI' filesep 'antsWarped.nii.gz']);
inds = cellfun(@(x) any(contains(x,subj_chars)),{files.folder});
files = files(inds);
% inds = randi(length(files),5,1);
inds = 1:length(files);
data_in = cell(length(inds),1);
for i = 1:length(inds)
    ind = inds(i);
    fprintf('Loading subject %s...\n',subj_chars{ind})
    tt = tic;
    tmp = strsplit(files(ind).name,'.');
    tmp = strjoin([tmp(1),tmp(2)],'.');
    if ~exist([files(ind).folder filesep tmp])
        fout = gunzip([files(ind).folder filesep files(ind).name]);
    else
        fout = {[files(ind).folder filesep tmp]};
    end
    data_in(i) = fout;
    toc(tt);
end

%##
% job = struct('data',{data_in},...
%         'globals',false);
% cat_stat_quality_measures(job);

%##
job = struct('data',{[HIRES_TEMPLATE; data_in]},...
    'globals',false,...
    'data_xml','antsWarped.xml',...
    'verb',true,...
    'new_fig',true,...
    'show_violin',true);
varargout = cat_stat_homogeneity(job);


