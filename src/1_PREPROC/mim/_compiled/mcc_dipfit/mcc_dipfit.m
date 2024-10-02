function [error_code] = mim_mcc_dipfit(working_dir,eeg_fpath,source_out_fpath,varargin)
%MIM_FEHEADMODEL_DIPFIT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/run_mim_mcc_dipfit_exe.sh
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
cat_logo();
%## TIME
startj = tic;
%## DEFINE DEFAULTS
%- define output
RV_THRESHOLD = 0.5;
% FORCE_RECREATE = 0;
error_code = 0;
%## working directory containing the ctf_fiducials.mat & mri_acpc_rs.mat &
%mri_acpc.mat & CustomElectrodeLocations.txt
errorMsg = 'Value must be CHAR. working directory containing the *_masks_contr.nii.gz, CustomElectrodeLocations.txt, mri_acpc.mat, elec_aligned.mat, & ctf_fiducials.mat'; 
wd_validFcn = @(x) assert(ischar(x)  && exist([x filesep 'ctf_fiducials.mat'],'file')...
    && ~isempty(dir([x filesep '*_masks_contr.nii.gz']))...
    && exist([x filesep 'CustomElectrodeLocations.txt'],'file')...
    && exist([x filesep 'mri_acpc.mat'],'file'),errorMsg);
%## EEG filepath
errorMsg = 'Value ''eeg_fpath'' must be PATH. ''eeg_fpath'' should point to a EEGLAB .set file'; 
ef_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%## Output for source.mat filepath
errorMsg = 'Value ''source_out_fpath'' must be PATH. ''source_out_fpath'' should point to a folder that exists';
sof_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%## Conductivity Values for volume
VOL_CONDUCTIVITIES = "[1.65,0.33,0.33,0.01,0.126,2.5*10^(-14)]";
vc_validFcn = @(x) ischar(x) && isnumeric(str2num(x));
%## (FLAG) Force recreation if needed
FORCE_RECREATE = '0';
fr_validFcn = @(x) ischar(x) && isnumeric(str2double(x));
%%
%- Checks
% fprintf(['CAT CODE\n',...
% ' _._     _,-''""`-._     \n',...
% '(,-.`._,''(       |\\`-/| \n',...
% '    `-.-'' \\ )-`( , o o) \n',...
% '          `-    \\`_`"''- \n',...
% 'Code Designer: Chang Liu, Jacob Salminen\n',...
% 'Code Date: 02/21/2024, MATLAB 2020a\n',...
% 'Copyright (C) Jacob Salminen, jsalminen@ufl.edu\n',...
% 'Copyright (C) Chang Liu, liu.chang1@ufl.edu\n']);
fprintf('Checking working_dir (%s) for ''subject_str''_masks_contr.nii.gz & CustomElectrodeLocations.txt\n',working_dir);
fprintf('ctf_fiducials.mat check: %i\n',exist([working_dir filesep 'ctf_fiducials.mat'],'file'));
fprintf('mri_acpc.mat check: %i\n',exist([working_dir filesep 'mri_acpc.mat'],'file'));
dir_segmri = dir([working_dir filesep '*_masks_contr.nii.gz']);
fprintf('''subject_str''_masks_contr.nii.gz check: %i\n',~isempty(dir_segmri));
fprintf('CustomElectrodeLocations.txt check: %i\n',exist([working_dir filesep 'CustomElectrodeLocations.txt'],'file'))
%- not required but nice to have
fprintf('elec_aligned.mat check: %i\n',exist([working_dir filesep 'elec_aligned.mat'],'file'))
fprintf('vol.mat check: %i\n',exist([working_dir filesep 'vol.mat'],'file'))
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'working_dir',wd_validFcn);
addRequired(p,'eeg_fpath',ef_validFcn);
addRequired(p,'source_out_fpath',sof_validFcn);
%## PARAMETER
addParameter(p,'VOL_CONDUCTIVITIES',VOL_CONDUCTIVITIES,vc_validFcn);
addParameter(p,'FORCE_RECREATE',FORCE_RECREATE,fr_validFcn);
parse(p,working_dir,eeg_fpath,source_out_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
VOL_CONDUCTIVITIES = p.Results.VOL_CONDUCTIVITIES;
VOL_CONDUCTIVITIES = str2num(VOL_CONDUCTIVITIES);
FORCE_RECREATE = p.Results.FORCE_RECREATE;
FORCE_RECREATE = logical(str2double(FORCE_RECREATE));
%- PARAMETER
%% Prints
fprintf('VOL_CONDUCTIVITIES: [csf, gray, scalp, skull, white, air]\n');
fprintf('VOL_CONDUCTIVITIES: [');fprintf('%0.3g, ',VOL_CONDUCTIVITIES(1:end-1));fprintf('%0.3g]\n',VOL_CONDUCTIVITIES(end));
%-
fprintf('FORCE_RECREATE: %i\n',FORCE_RECREATE);
%% ===================================================================== %%
% fid = fopen([working_dir filesep 'output.txt'],'w');
%- EEGLAB options for opening EEG data
try
    fprintf(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID') '\n']);
    fprintf(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE') '\n']);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([working_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([working_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', inf);
    disp(pPool);
catch e
    fprintf('Parallel processing failed to start');
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
end
%% ===================================================================== %%
fprintf('Running dipole fitting on directory: %s\n',working_dir);
%## Create headmodel
% out_segmri = [dir_segmri.folder filesep dir_segmri.name];
if ~exist([working_dir filesep 'vol.mat'],'file') || ~exist([working_dir filesep 'elec_aligned.mat'],'file') || FORCE_RECREATE 
    [vol,elec_aligned] = fem_create_vol(working_dir,...
        [working_dir filesep 'ctf_fiducials.mat'],...
        [dir_segmri.folder filesep dir_segmri.name],...
        [working_dir filesep 'CustomElectrodeLocations.txt'],...
        [working_dir filesep 'mri_acpc.mat'],...
        VOL_CONDUCTIVITIES);
else
    %- load vol.mat
    tmp = load([working_dir filesep 'vol.mat']);
    vol = tmp.vol;
    tmp = load([working_dir filesep 'elec_aligned.mat']);
    %- load elec_aligned.mat
    try
    elec_aligned = tmp.elec_aligned;
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
        elec_aligned = tmp.elec_aligned_init;
    end
    clear tmp
end
%% PREPARE VOLUME AND SENSORS
%- prepare volume and sensors (projects sensors to scalp)
if ~exist([working_dir filesep 'headmodel_fem_tr.mat'],'file') || FORCE_RECREATE
    starti = tic;
    fprintf('Preparing vol.mat and elec_aligned.mat...\n');
    %## Override Simbio path to utilize parallel processing
    % this is done when compiling.
    %- prepare volume and sensors
    [headmodel_fem_tr, ~] = ft_prepare_vol_sens(vol, elec_aligned);
    disp(headmodel_fem_tr);
    save([working_dir filesep 'headmodel_fem_tr.mat'],'headmodel_fem_tr','-v7.3');
    endi = toc(starti);
    fprintf('Call to ft_prepare_leadfield.m took %0.2g',endi);
else
    fprintf('Loading headmodel_fem_tr.mat...\n');
    tmp = load([working_dir filesep 'headmodel_fem_tr.mat']);
    headmodel_fem_tr = tmp.headmodel_fem_tr;
end
%% SOURCEMODEL CALCULATION
if ~exist([working_dir filesep 'sourcemodel.mat'],'file') || FORCE_RECREATE
    starti = tic;
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
    endi = toc(starti);
    fprintf('Call to ft_prepare_leadfield.m took %0.2g',endi);
else
    fprintf('Loading sourcemodel.mat...\n');
    tmp = load([working_dir filesep 'sourcemodel.mat']);
    sourcemodel = tmp.sourcemodel;
end

%- Remove this field to force average referencing of leadfield matrix
try
    elec_aligned    = rmfield(elec_aligned,'tra'); 
catch e
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
end
%% LEADFIELD CALCULATION
%##
% choose the available When the forward solution is computed, the lead 
% field matrix (= channels X source points matrix) is calculated for each 
% grid point taking into account the head model and the channel positions.
if ~exist([working_dir filesep 'leadfield_fem.mat'],'file') || FORCE_RECREATE
    starti = tic;
    fprintf('Computing Leadfield...\n');
    cfg             = [];
    cfg.grid        = sourcemodel;
    cfg.headmodel   = headmodel_fem_tr;
    cfg.elec        = elec_aligned;
    cfg.reducerank  = 'no';
    % cfg.normalize   = 'column'; 
    leadfield_fem   = ft_prepare_leadfield(cfg);% This actually didn't take that long - ?1r without parallel processing
    save([working_dir filesep 'leadfield_fem.mat'],'leadfield_fem');
    endi = toc(starti);
    fprintf('Call to ft_prepare_leadfield.m took %0.2g',endi);
else
    fprintf('Loading leadfield_fem.mat...\n');
    tmp = load([working_dir filesep 'leadfield_fem.mat']);
    leadfield_fem = tmp.leadfield_fem;
end
clear sourcemodel vol
%% Dipole Fitting
%## COARSE FIT
tmp = strsplit(eeg_fpath,filesep);
fName = tmp{end};
fPath = strjoin(tmp(1:end-1),filesep);
%- load EEG
% EEG = pop_loadset('filepath',fPath,'filename',fName);
fprintf('Loading EEG...\n');
EEG = load('-mat', [fPath filesep fName]);
fid_eeg = fopen([fPath filesep EEG.data], 'r', 'ieee-le');
data = fread(fid_eeg, [EEG.trials*EEG.pnts EEG.nbchan], 'float32')';
EEG.data = data;
fclose(fid_eeg);
%## REMOVE ICA RESULTS FROM PREVIOUS ANALYSIS
fprintf('Loading ICA...\n');
EEG = rmfield(EEG,'icaweights');
EEG.icaweights = [];
EEG = rmfield(EEG,'icawinv');
EEG.icawinv = [];
EEG = rmfield(EEG,'icaact');
EEG.icaact = [];
EEG = rmfield(EEG,'icasphere');
EEG.icasphere = [];
EEG = rmfield(EEG,'icachansind');
EEG.icachansind = []; %also remove icachanind?
%## LOAD CURRENT AMICA RESULTS FOR subjStr & subDirNum
EEG = pop_loadmodout(EEG,fPath);
EEG.dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
EEG.dipfit.mrifile = [];
EEG.dipfit.hdmfile = [working_dir filesep 'headmodel_fem_tr.mat'];
EEG.dipfit.coordformat = [];
EEG.dipfit.chanfile = [];
EEG.dipfit.chansel = (1:size(EEG.icawinv,2));
%## Convert to Fieldtrip
fprintf('Converting EEGLAB to Fieldtrip...\n');
ftEEG = eeglab2fieldtrip(EEG,'componentanalysis','dipfit');
clear EEG 
%- ft_dipolefitting
%## NONLINEAR FIT
% if ~exist(source_out_fpath, 'dir')
%     mkdir(source_out_fpath)
% end
sources = cell(size(ftEEG.topo,2),1);
starti = tic;
parfor (comp_i = 1:size(ftEEG.topo,2),POOL_SIZE)
    cfg = [];
    cfg.numdipoles      =  1;
    cfg.headmodel       = headmodel_fem_tr;
    cfg.sourcemodel     = leadfield_fem;
    cfg.elec            = elec_aligned; %elec_aligned;
    cfg.dipfit.metric   = 'rv';
    cfg.nonlinear       = 'no';
    cfg.component       = comp_i;
    source              = ft_dipolefitting(cfg,ftEEG);
    sources{comp_i}     = source;
    display(source)
end
endi = toc(starti);
fprintf('Call to ft_dipolefitting.m coarse took %0.2g',endi);
sources = cellfun(@(x) [[]; x], sources);
par_save(sources,source_out_fpath,'dipfit_struct_coarse.mat');
%- filter by residual variance by removing those > 50%
tmp = [sources.dip];
inds = find([tmp.rv] < RV_THRESHOLD);
clear tmp
%## NONLINEAR FIT
sources = cell(size(ftEEG.topo,2),1);
starti = tic;
parfor (comp_i = 1:size(ftEEG.topo,2),POOL_SIZE)
    cfg = [];
    cfg.numdipoles      =  1;
    cfg.headmodel       = headmodel_fem_tr;
    cfg.sourcemodel     = leadfield_fem;
    cfg.elec            = elec_aligned; %elec_aligned;
    cfg.dipfit.metric   = 'rv';
    cfg.nonlinear       = 'yes';
    cfg.component       = comp_i;
    if any(comp_i == inds)
        source              = ft_dipolefitting(cfg,ftEEG);
        sources{comp_i}     = source;
        display(source)
    else
        fprintf('Skipping component %i due to a high residual variance in coarse fit...\n',comp_i);
        sources{comp_i} = struct('label',cell(1,1),'dip',[],'Vdata',double(0),'Vmodel',double(0),'component',comp_i,'cfg',[]);
    end
end
endi = toc(starti);
fprintf('Call to ft_dipolefitting.m nonlinear took %0.2g',endi);
sources = cellfun(@(x) [[]; x], sources);
par_save(sources,source_out_fpath,'dipfit_struct.mat');
%## TIME
endj = toc(startj);
fprintf('Call to mim_mcc_dipfit.m took %0.2g',endj-startj);
%## EXIT
% fclose(fid);
% exit(error_code);
end
%% ===================================================================== %%
function [] = par_save(SAVEVAR,fPath,fName,varargin)
%PAR_SAVE Summary of this function goes here
%   Detailed explanation goes here
%
%   IN:
%       REQUIRED:
%           SAVEVAR, variable to save.
%               STRUCT, CHAR, DICT, ARRAY you want to save.
%           fPath, CHAR
%               path to the folder where your file is held
%           fName, CHAR
%               file name & extension (e.g., 'INEEG.mat')
%       OPTIONAL:
%           fname_ext, CHAR
%               for automating the renaming of file names 
%       PARAMETER:
%   OUT:
%       NONE
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/06/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%- fPath
errorMsg = 'Value ''fPath'' must be CHAR. The path must exist.'; 
fp_validFcn = @(x) assert(ischar(x) && exist(x,'dir'),errorMsg);
%- fname_ext
% FNAME_EXT = '';
% errorMsg = 'Value ''fname_ext'' must be CHAR. This value is appended to ''fName'' before the file declaration.'; 
% fn_validFcn = @(x) assert(ischar(x),errorMsg);
p = inputParser;
%## REQUIRED
addRequired(p,'SAVEVAR')
addRequired(p,'fPath',fp_validFcn)
addRequired(p,'fName',@ischar)
%## OPTIONAL
%## PARAMETER
parse(p, SAVEVAR, fPath, fName, varargin{:});
%## SET DEFAULTS
%% ===================================================================== %%
%- save
s = whos('SAVEVAR');
fprintf(1,'%s is %0.2g bytes\n',fName,s.bytes);
% fprintf('\nSaving %s to\n%s\n',fName,savePath);
if s.bytes >= 1.9e9
    fprintf('\nSaving %s using ''v7.3'' to\n%s\n',fName,fPath);
    save([fPath filesep fName],'SAVEVAR','-v7.3');
else
    fprintf('\nSaving %s using ''v6'' to\n%s\n',fName,fPath);
    save([fPath filesep fName],'SAVEVAR','-v6');
end
end


function [vol,elec_aligned] = fem_create_vol(output_dir,ctf_fiducial_fpath,segmri_fpath,customelec_fpath,acpcmri_fpath,varargin)
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
DO_PLOTTING = false;
% error_code = 0; %#ok<NASGU>
%- working directory containing the ctf_fiducials.mat & mri_acpc_rs.mat &
%mri_acpc.mat & CustomElectrodeLocations.txt
% errorMsg = 'Value must be CHAR. working directory containing the ctf_fiducials.mat & elec_aligned.mat'; 
% wd_validFcn = @(x) assert(ischar(x) && exist([x filesep 'ctf_fiducials.mat'],'file') && exist([x filesep 'mri_acpc_rs.mat'],'file') && exist([x filesep 'mri_acpc.mat'],'file') && exist([x filesep 'CustomElectrodeLocations.txt'],'file'),errorMsg);

validFcn_1 = @(x) assert(ischar(x) && exist(x,'file'),...
    'Value must be CHAR. ctf_fiducials.mat does not exist');
validFcn_2 = @(x) assert(ischar(x) && exist(x,'file'),...
    'Value must be CHAR. *_masks_contr.nii.gz does not exist');
validFcn_3 = @(x) assert(ischar(x) && exist(x,'file'),...
    'Value must be CHAR. CustomElectrodeLocations.txt does not exist');
validFcn_4 = @(x) assert(ischar(x) && exist(x,'file'),...
    'Value must be CHAR. mri_acpc.mat does not exist');
%-
VOL_CONDUCTIVITIES = [1.65,0.33,0.33,0.01,0.126,2.5*10^(-14)];
vc_validFcn = @(x) isnumeric(x) && length(x) == 6;
%## CHECK FPATHS
% fprintf('Checking working_dir (%s) for ''subject_str''_masks_contr.nii.gz & CustomElectrodeLocations.txt\n',working_dir);
fprintf('ctf_fiducials.mat check: %i\nPath: %s\n',exist(ctf_fiducial_fpath,'file'),ctf_fiducial_fpath);
fprintf('mri_acpc.mat check: %i\nPath: %s\n',exist(acpcmri_fpath,'file'),acpcmri_fpath);
fprintf('''subject_str''_masks_contr.nii.gz check: %i\nPath: %s\n',exist(segmri_fpath,'file'),segmri_fpath);
fprintf('CustomElectrodeLocations.txt check: %i\nPath: %s\n',exist(customelec_fpath,'file'),customelec_fpath);
fprintf('Using volume conductivities of: '); fprintf('%0.2g,',VOL_CONDUCTIVITIES(1:end-1)); fprintf('%0.2g\n',VOL_CONDUCTIVITIES(end));
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'output_dir',@(x) assert(ischar(x) && exist(x,'dir'),'output_dir does not exist'));
addRequired(p,'ctf_fiducial_fpath',validFcn_1);
addRequired(p,'segmri_fpath',validFcn_2);
addRequired(p,'customelec_fpath',validFcn_3);
addRequired(p,'acpcmri_fpath',validFcn_4);
%## OPTIONAL
addOptional(p,'VOL_CONDUCTIVITES',VOL_CONDUCTIVITIES,vc_validFcn)
%## PARAMETER
parse(p,output_dir,ctf_fiducial_fpath,segmri_fpath,customelec_fpath,acpcmri_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
VOL_CONDUCTIVITIES = p.Results.VOL_CONDUCTIVITES;
%- PARAMETER
%% ===================================================================== %%
fprintf('Output directory: %s\n',output_dir);
%- load mri in acpa & ctf coordinate systems & fiducial marks
% tmp = load([working_dir filesep 'mri_acpc_rs.mat']);
% mri_acpc_rs = tmp.mri_acpc_rs;
tmp = load(acpcmri_fpath);
mri_acpc = tmp.mri_acpc;
tmp = load(ctf_fiducial_fpath);
ctf_fiducials = tmp.ctf_fiducials;
%- Load the electrodes after digitized
chanloc_scan_folder = customelec_fpath;
chanlocs = readtable(chanloc_scan_folder);% Same output text file from getchalocs.
chanlocs.Properties.VariableNames = {'labels','X','Y','Z'};
elec.chanpos(:,1) = [chanlocs.X];
elec.chanpos(:,2) = [chanlocs.Y];
elec.chanpos(:,3) = [chanlocs.Z];
elec.elecpos      = elec.chanpos;
elec.label(:,1)   = [chanlocs.labels]';
%%
unzip_out = gunzip(segmri_fpath);
simnibs_mask = ft_read_mri(unzip_out{1});
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
save([output_dir filesep 'mesh.mat'],'mesh')
%% FIELDTRIP CREATE CONDUCTIVITY VOLUME (SIMBIO)
cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = zeros(1,5);
% scale = 1;
% order follows mesh.tissyelabel , CAUTIOUS!!!! OMg this is not the same order as in the segmentation
cfg.conductivity(strcmp(mesh.tissuelabel,'csf')) = VOL_CONDUCTIVITIES(1); %1.65*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'gray')) = VOL_CONDUCTIVITIES(2); %0.33*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'scalp')) = VOL_CONDUCTIVITIES(3); %0.33*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'skull')) = VOL_CONDUCTIVITIES(4); %0.01*scale; %0.0042*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'white')) = VOL_CONDUCTIVITIES(5); %0.126*scale;
cfg.conductivity(strcmp(mesh.tissuelabel,'air')) = VOL_CONDUCTIVITIES(6); %2.5*10^(-14)*scale;
vol = ft_prepare_headmodel(cfg, mesh);
fprintf('Saving vol file\n');
save([output_dir filesep 'vol.mat'],'vol','-v6') 

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
save([output_dir filesep 'elec_aligned_init.mat'],'elec_aligned_init') 

%## Refined Electrode Alignment
cfg               = [];
cfg.method        = 'project';
cfg.elec          = elec_aligned_init;
cfg.headshape     = vol;
elec_aligned      = ft_electroderealign(cfg);
fprintf('Saving elec_aligned file\n');
save([output_dir filesep 'elec_aligned.mat'],'elec_aligned') 
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
    saveas(fig_i,[output_dir filesep sprintf('ft_sourceplot.fig')]);
    saveas(fig_i,[output_dir filesep sprintf('ft_sourceplot.jpg')]);
    %## (PLOT 2)
    figure;
    ft_plot_mesh(mesh, 'surfaceonly', 'yes','facecolor','b','edgecolor', 'none', 'facealpha', 0.4);
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_mesh.fig')]);
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_mesh.jpg')]);
    %## (PLOT 4) Initial alignment of electrodes using fiducial marks
    figure;
    hold on;
    ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    camlight
    ft_plot_sens(elec_aligned_init,'style','.r');
    ft_plot_sens(fid,'style','xb');%plot fiducial points
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_sens_1.fig')]);
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_sens_1.jpg')]);
    %## (PLOT 3) final alignment of electrodes after projecting to scalp
    figure;
    hold on;
    ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
    camlight
    ft_plot_sens(elec_aligned ,'style','.r');
    ft_plot_sens(fid,'style','xb');%plot fiducial points
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_sens_2.fig')]);
    saveas(fig_i,[output_dir filesep sprintf('ft_plot_sens_2.jpg')]);
    
end
% error_code = 1;
end

%% (NOTES)
% (REQUIRED) ft_dipolefitting
% cfg.model           = g_model; %'moving';
% cfg.nonlinear       = g_nonlinear; %'yes';
% cfg.numdipoles      = g_numdipoles; %1; %number of dipoles, if > 1 you need to define the symmetry variable
% cfg.resolution      = g_resolution; %10; %resolution of model in units <cfg.unit>, exponetially increases computation time
% cfg.unit            = g_unit; %units; %units of headmodel 
% cfg.gridsearch      = g_gridsearch; %'yes'; %gridsearch for initial params (increases processing time, but gives stronger fit).
% cfg.dipfit.maxiter  = 200;
% cfg.reducerank      = g_reducerank; %'no';
% cfg.spmversion      = g_spmversion; %'spm12'; %'cat12'?
% cfg.vol             = vol; %the headmodel created from the MRI or MNI
% cfg.senstype        = g_senstype; %sensor type
% cfg.elec            = elec; %the channels you want to consider (keeping as is will use all EEG sensors)
% OPTIONS
% cfg.warpmni         = true;
% cfg.channel         = eegChan;
% cfg.leadfield       = leadfieldFEM;
