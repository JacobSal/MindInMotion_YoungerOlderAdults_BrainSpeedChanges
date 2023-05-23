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

%## TIME
startj = tic;
%## DEFINE DEFAULTS
%- define output
RV_THRESHOLD = 0.5;
error_code = 0;
%- working directory containing the vol.mat & elec_aligned.mat
errorMsg = 'Value ''working_dir'' must be PATH. ''working_dir'' should point to a folder. Working directory must contain a fieldtrip generated vol.mat & elec_aligned.mat'; 
wd_validFcn = @(x) assert(exist(x,'dir') && exist([x filesep 'vol.mat'],'file') && exist([x filesep 'elec_aligned.mat'],'file'),errorMsg);
fprintf('Checking working_dir (%s) for vol.mat & elec_aligned.mat\n',working_dir);
fprintf('vol.mat check: %i\n',exist([working_dir filesep 'vol.mat'],'file'))
fprintf('elec_aligned.mat check: %i\n',exist([working_dir filesep 'elec_aligned.mat'],'file'))
%- EEG filepath
errorMsg = 'Value ''eeg_fpath'' must be PATH. ''eeg_fpath'' should point to a EEGLAB .set file'; 
ef_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%- Output for source.mat filepath
errorMsg = 'Value ''source_out_fpath'' must be PATH. ''source_out_fpath'' should point to a folder that exists';
sof_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'working_dir',wd_validFcn);
addRequired(p,'eeg_fpath',ef_validFcn);
addRequired(p,'source_out_fpath',sof_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,working_dir,eeg_fpath,source_out_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
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
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', 1440);
    disp(pPool);
catch e
    fprintf('Parallel processing failed to start');
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
end
%% ===================================================================== %%
fprintf('Running dipole fitting on directory: %s\n',working_dir);
%## Load Vars
%- load vol.mat
tmp = load([working_dir filesep 'vol.mat']);
vol = tmp.vol;
%- load elec_aligned.mat
tmp = load([working_dir filesep 'elec_aligned.mat']);
try
    elec_aligned = tmp.elec_aligned;
catch e
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
    elec_aligned = tmp.elec_aligned_init;
end
clear tmp
%% PREPARE VOLUME AND SENSORS
%- prepare volume and sensors (projects sensors to scalp)
if ~exist([working_dir filesep 'headmodel_fem_tr.mat'],'file')
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
if ~exist([working_dir filesep 'sourcemodel.mat'],'file')
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
if ~exist([working_dir filesep 'leadfield_fem.mat'],'file')
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

function [] = par_save(SAVEVAR,fPath,fName)
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
% tic
%## DEFINE DEFAULTS
%- fPath
errorMsg = 'Value ''fPath'' must be CHAR. The path must exist.'; 
fp_validFcn = @(x) assert(ischar(x) && exist(x,'dir'),errorMsg);
%- fName
errorMsg = 'Value ''fname_ext'' must be CHAR. This value is appended to ''fName'' before the file declaration.'; 
fn_validFcn = @(x) assert(ischar(x),errorMsg);
p = inputParser;
%## REQUIRED
addRequired(p,'SAVEVAR')
addRequired(p,'fPath',fp_validFcn)
addRequired(p,'fName',fn_validFcn)
%## OPTIONAL
%## PARAMETER
parse(p, SAVEVAR, fPath, fName);
%## SET DEFAULTS
%% ===================================================================== %%
if ~exist(fPath, 'dir')
    mkdir(fPath)
end
save([fPath filesep fName],'SAVEVAR','-v6');
fprintf('\nSaving %s to\n%s\n',fName,fPath);
%## TIME
% toc
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