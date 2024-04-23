function [err] = mcc_epoch_process(ica_dir,study_dir,jf_fpath,varargin)
%MIM_SPCA_MCC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%       REQUIRED:
%       ica_dir, CHAR
%           directory where ICA'd EEG data is located for each subject
%           you'd like to run.
%           E.G.,
%               my_ica_folder/
%               ├─ subj_1/
%               │  ├─ clean/
%               ├─ subj_2/
%               │  ├─ clean/
%               ├─ subj_n/
%               │  ├─ clean/
%       study_dir, CHAR
%           directory where the program will store information about the
%           subjects in 'ica_dir'.
%
%       jf_fpath, CHAR
%           JSON Format is as follows:
%               {
% 	                 "EPOCH_PARAMS": {
% 		                "percent_overlap": 0,
% 		                "epoch_event_char": "RHS",
% 		                "epoch_time_lims": [-0.5,4.5],
% 		                "tw_stdev": 3,
% 		                "tw_events": ["RHS","LTO","LHS","RTO","RHS"],
% 		                "path_ext": "gait_epoched",
% 		                "gait_trial_chars": ["0p25","0p5","0p75","1p0","flat","low","med","high"],
% 		                "rest_trial_char": [],
% 		                "do_recalc_epoch": true
% 	                },
% 	                "SUBJ_CHARS": null
%               }
% 
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
%% (PARAMETERS) ======================================================== %%
%## hard define
cat_logo();
str = '';
err = 0;
ADD_CLEANING_SUBMODS = false;
ADD_DIPFIT_COMPILE_SUBMODS = false;
SAVE_SEPERATE_EEG = false;
DO_SLIDING_WINDOW = false;
study_fname_icared = 'cont_icreduced_study';
study_fname_gait = 'epoch_study';
%## EPOCH PARAMS
EPOCH_PARAMS = struct('percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%## SUBJECT CHARACTERS LOAD-IN FROM FUNCTION
[SUBJ_CHARS,GROUP_NAMES,~,~,~,~,~] = mim_dataset_information('yaoa_spca');
SUBJ_CHARS = [SUBJ_CHARS{:}];
%% LOAD JSON & ASSIGN CHECKS
%- connectivity statistics & surrogates params
if ~isempty(jf_fpath)
    fid = fopen(jf_fpath,'r');
    raw = fread(fid,inf); 
    str = char(raw');
    fclose(fid); 
    params = jsondecode(str);
    if ~isempty(params.EPOCH_PARAMS)
        tmp_epoch = params.EPOCH_PARAMS;
        varargin = [varargin, 'EPOCH_PARAMS', tmp_epoch];
    end
    if ~isempty(params.SUBJ_CHARS)
        tmp_subj_chars = params.SUBJ_CHARS;
        varargin = [varargin, 'SUBJ_CHARS', tmp_subj_chars];
    end
end
%## working directory containing subject ICA .set files
errorMsg = 'Value must be CHAR. Working directory containing all subjects to be epoched & sPCA''d.'; 
id_validFcn = @(x) assert(ischar(x)  && exist(x,'dir'),errorMsg);
errorMsg = 'Value must be CHAR. New directory to store epoched .set files and study information.'; 
sd_validFcn = @(x) assert(ischar(x),errorMsg);
errorMsg = 'Value must be CHAR. File path to a .json file containing parameters to be used.'; 
jf_validFcn = @(x) assert((ischar(x) && exist(x,'file')) || isempty(x),errorMsg);
%% (PARSER) ============================================================ %%
fprintf('ica_dir: %s\n',ica_dir);
fprintf('study_dir: %s\n',study_dir);
fprintf('jf_fpath: %s\n',jf_fpath);
fprintf('%s\n',str);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ica_dir',id_validFcn);
addRequired(p,'study_dir',sd_validFcn);
addRequired(p,'jf_fpath',jf_validFcn);
%## OPTIONAL
addParameter(p,'SUBJ_CHARS',SUBJ_CHARS,@iscell);
addParameter(p,'EPOCH_PARAMS',EPOCH_PARAMS,@(x) validate_struct(x,EPOCH_PARAMS));
%## PARAMETERS
parse(p,ica_dir,study_dir,jf_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
SUBJ_CHARS = p.Results.SUBJ_CHARS;
EPOCH_PARAMS = p.Results.EPOCH_PARAMS;
if ~exist(study_dir,'dir')
    mkdir(study_dir);
end
%% PARPOOL SETUP ======================================================= %%
[PATHS,SLURM_POOL_SIZE,ALLEEG,STUDY,CURRENTSTUDY,CURRENTSET] = set_workspace();
%- create new study directory
if ~exist(conn_save_dir,'dir')
    mkdir(conn_save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([STUDY_FPATH filesep STUDY_FNAME '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
    %## LOAD STUDY
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',STUDY_FPATH);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',STUDY_FPATH);
    end
    cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
    MAIN_STUDY.cluster = cl_struct;
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
end
%% MIM STUDY STATS
%## FIND MINIMUM TRIALS ACROSS COHORT
if DO_STANDARD_TRIALS
    subj_chars = {MAIN_STUDY.datasetinfo.subject}';
    trial_counts = cell(length(subj_chars),length(TRIAL_TYPES));
    channel_rejs = cell(length(subj_chars),1);
    sample_rej_perc = cell(length(subj_chars),1);
    for subj_i=1:length(subj_chars)
        EEG = MAIN_ALLEEG(subj_i);
        for cond_i = 1:length(TRIAL_TYPES)
            tmp = sum(strcmp({MAIN_STUDY.datasetinfo(subj_i).trialinfo.cond},TRIAL_TYPES{cond_i}));
            trial_counts{subj_i,cond_i} = tmp;
        end
        channel_rejs{subj_i,1} = sum(~EEG.etc.channel_mask(strcmp({EEG.urchanlocs.type},'EEG')));
        sample_rej_perc{subj_i,1} = (length(EEG.etc.clean_sample_mask)-sum(EEG.etc.clean_sample_mask))/length(EEG.etc.clean_sample_mask);
    end
    tbl_out = table(subj_chars,trial_counts(:,1),trial_counts(:,2),trial_counts(:,3),...
        trial_counts(:,4),trial_counts(:,5),trial_counts(:,6),trial_counts(:,7),trial_counts(:,8),channel_rejs,sample_rej_perc,'VariableNames',{'subj_chars','0p25','0p5','0p75','1p0','flat','low','med','high','channel_rejs','sample_rej_perc'});
    vals = tbl_out{:,2:9};
    trial_mins = zeros(size(vals,2),1);
    trial_maxs = zeros(size(vals,2),1);
    for i = 1:size(vals,2)
        trial_mins(i) = min([vals{:,i}]);
        trial_maxs(i) = max([vals{:,i}]);
    end
    MIN_CONN_TRIALS = min(trial_mins);
end
%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% for subj_i = LOOP_VAR
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    eeg_savefpath = [fPaths{subj_i} filesep 'conn'];
    if ~exist(eeg_savefpath,'dir')
        mkdir(eeg_savefpath)
    end
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
    end
    EEG = pop_resample(EEG,RESAMLE_FREQ);
    try
        %## RUN MAIN_FUNC
        [TMP,t_out] = cnctanl_sift_pipe(EEG,components,conn_save_dir,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'PREPDATA',DEF_PREPDATA,...
            'ESTSELMOD',DEF_ESTSELMOD,...
            'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
            'ESTFITMVAR',DEF_ESTFITMVAR,...
            'PLOTORDERCRIT',DEF_PLOTORDERCRIT);
%         for i=1:length(TMP)
%             par_save(TMP(i).CAT)
%         end
        BIG_CAT = cat(1,TMP(:).CAT);
        EEG.etc.COND_CAT = BIG_CAT;
        EEG.etc.conn_table = t_out;
        EEG.etc.conn_meta.comps_out = comps_out;
        fName = strsplit(EEG.filename,'.'); fName = [fName{1} '.mat'];
        par_save(t_out,eeg_savefpath,fName,'_conntable');
        [EEG] = pop_saveset(EEG,...
            'filepath',eeg_savefpath,'filename',EEG.filename,...
            'savemode','twofiles');
        tmp{subj_i} = EEG;
    catch e
        rmv_subj(subj_i) = 1;
        EEG.CAT = struct([]);
        tmp{subj_i} = EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
        disp(e);
    end
end
pop_editoptions('option_computeica',0);
%% SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s',MAIN_ALLEEG(cellfun(@isempty,tmp)).subject);
tmp = tmp(~cellfun(@isempty,tmp));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding %s %s\n',EEG.subject,addFs{j})
        end
    end 
    tmp{subj_i} = EEG;
end
tmp = cellfun(@(x) [[]; x], tmp);
%##
tmp = eeg_checkset(tmp,'eventconsistency');
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_SAVE,...
                                'filename',STUDY_FNAME_SAVE,...
                                'filepath',STUDY_FPATH);
%## ASSIGN PARAMETERS
STUDY.etc.a_epoch_process.epoch_chars = TRIAL_TYPES;
STUDY.etc.d_cnctanl_process.params = struct('DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'PREPDATA',DEF_PREPDATA,...
            'ESTSELMOD',DEF_ESTSELMOD,...
            'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
            'ESTFITMVAR',DEF_ESTFITMVAR,...
            'PLOTORDERCRIT',DEF_PLOTORDERCRIT,...
            'comps',comps_out);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY_FNAME_SAVE,STUDY_FPATH,...
                                            'STUDY_COND',[]);      
end
%% (SUBFUNCTIONS) ====================================================== %%
function [PATHS,SLURM_POOL_SIZE,ALLEEG,STUDY,CURRENTSTUDY,CURRENTSET] = set_workspace()
    %## TIME
    TT = tic;
    %## ================================================================= %%
    TMP_PWD = dir(['.' filesep]);
    TMP_PWD = TMP_PWD(1).folder;
    fprintf(1,'Current folder: %s\n',TMP_PWD);
    %## datetime
    dt         = datetime;
    dt.Format  = 'ddMMyyyy';
    %## PARAMS ========================================================== %%
    tmp = strsplit(TMP_PWD,filesep);
    src_ind = find(strcmp(tmp,'src'));
    if ~ispc
        %- Add source directory where setup functions are 
        src_dir = [filesep strjoin(tmp(1:src_ind),filesep)];
    else
        %- Add source directory where setup functions are
        src_dir = strjoin(tmp(1:src_ind),filesep);
    end
    if ~ispc
        %- Add submodules directory where packages are 
        submodules_dir = [filesep strjoin(tmp(1:src_ind-1),filesep) filesep 'submodules'];
    else
        %- Add submodules directory where packages are 
        submodules_dir = [strjoin(tmp(1:src_ind-1),filesep) filesep 'submodules'];
    end
    %## FUNCTIONS FOLDER
    FUNC_FPATH = [src_dir filesep '_functions' filesep 'v2_0'];
    %##
    fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n-FUNCTIONS: %s\n',src_dir,submodules_dir,FUNC_FPATH);
    %## HARDCODE PATHS STRUCT =========================================== %%
    PATHS = [];
    %- submods path
    PATHS.submods_dir = submodules_dir;
    %- src folder
    PATHS.src_dir = src_dir;
    %- _data folder
    PATHS.data_dir = fullfile(src_dir,'_data');
    %- EEGLAB folder
    PATHS.eeglab_dir = [submodules_dir filesep 'eeglab'];
    %- SIFT folder
    PATHS.eeglab_dir = [submodules_dir filesep 'SIFT'];
    %- _functions folder
    PATHS.functions_dir = FUNC_FPATH;
    %## ADDPATH for FIELDTRIP =========================================== %%
    ft_defaults;
    %- always start eeglab last.
    ALLEEG=[]; STUDY=[]; CURRENTSET=0; CURRENTSTUDY=0;
    eeglab;
    %## PARPOOL SETUP =================================================== %%
    if ~ispc
        pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
        'option_single', 1, 'option_memmapdata', 0, ...
        'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
        disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
        disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
        %## allocate slurm resources to parpool in matlab
        %- get cpu's on node and remove a few for parent script.
        SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
        %- create cluster
        pp = parcluster('local');
        %- Number of workers for processing (NOTE: this number should be higher
        %then the number of iterations in your for loop)
        fprintf('Number of workers: %i\n',pp.NumWorkers);
        fprintf('Number of threads: %i\n',pp.NumThreads);
        %- make meta data dire1ory for slurm
        mkdir([TMP_PWD filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')])
        pp.JobStorageLocation = [TMP_PWD filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')];
        %- create your p-pool (NOTE: gross!!)
        pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
    else
        pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
        'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
        'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
        SLURM_POOL_SIZE = 1;
    end
    %## TIME
    fprintf('Done with workplace setup: %0.2f',toc(TT));
end

%% (SUBFUNCTIONS) ====================================================== %%
function [p] = unix_genpath(d)
    %GENPATH Generate recursive toolbox path.
    %   P = GENPATH returns a character vector containing a path name 
    %   that includes all the folders and subfolders below MATLABROOT/toolbox, 
    %   including empty subfolders.
    %
    %   P = GENPATH(FOLDERNAME) returns a character vector containing a path 
    %   name that includes FOLDERNAME and all subfolders of FOLDERNAME, 
    %   including empty subfolders.
    %   
    %   NOTE 1: GENPATH will not exactly recreate the original MATLAB path.
    %
    %   NOTE 2: GENPATH only includes subfolders allowed on the MATLAB
    %   path.
    %
    %   See also PATH, ADDPATH, RMPATH, SAVEPATH.

    %   Copyright 1984-2018 The MathWorks, Inc.
    %------------------------------------------------------------------------------

    % String Adoption
    if nargin > 0
        d = convertStringsToChars(d);
    end

    if nargin==0
      p = genpath(fullfile(matlabroot,'toolbox'));
      if length(p) > 1, p(end) = []; end % Remove trailing pathsep
      return
    end

    % initialise variables
    classsep = '@';  % qualifier for overloaded class directories
    packagesep = '+';  % qualifier for overloaded package directories
    p = '';           % path to be returned

    % Generate path based on given root directory
    files = dir(d);
    if isempty(files)
      return
    end

    % Add d to the path even if it is empty.
    p = [p d pathsep];

    % set logical vector for subdirectory entries in d
    isdir = logical(cat(1,files.isdir));
    %
    % Recursively descend through directories which are neither
    % private nor "class" directories.
    %
    dirs = files(isdir); % select only directory entries from the current listing

    for i=1:length(dirs)
       dirname = dirs(i).name;
       if    ~strcmp( dirname,'.')          && ...
             ~strcmp( dirname,'..')         && ...
             ~strncmp( dirname,classsep,1) && ...
             ~strncmp( dirname,packagesep,1) && ...
             ~strcmp( dirname,'private')    && ...
             ~strcmp( dirname,'resources') && ...
             ~strcmp( dirname,'__archive')
          p = [p genpath([d filesep dirname])]; % recursive calling of this function.
       end
    end
end
function [ALLEEG,conn_mat_table] = cnctanl_sift_pipe(ALLEEG,conn_components,save_dir,varargin)
    %CNCTANL_SIFT_PIPE Summary of this function goes here
    %   Project Title: 
    %   Code Designer: Jacob salminen
    %
    %   Version History --> See details at the end of the script.
    %   Current Version:  v1.0.20220103.0
    %   Previous Version: n/a
    %   Summary:
    %
    %   IN:
    %       REQUIRED:
    %           EEG, Struct
    %               Can be the EEG struct or the ALLEEG struct if running multiple
    %               subjects?
    %           components, CELL
    %               these are the components/channels to which we'll fit our
    %               multivariate model. Each cell contains an array of INTS
    %               defining the componenets to select for that subject.
    %       OPTIONAL:
    %           epochTimeRange
    %       PARAMETER:
    %           
    %   OUT:
    %
    %
    % CAT CODE
    %  _._     _,-'""`-._
    % (,-.`._,'(       |\`-/|
    %     `-.-' \ )-`( , o o)
    %           `-    \`_`"'-
    % Code Designer: Jacob Salminen
    % Code Date: 12/30/2022, MATLAB 2019a
    % Copyright (C) Jacob Salminen, jsalminen@ufl.edu
    %## TIME
    tic
    %% (HIGH LEVEL VARS) =================================================== %%
    DO_BOOTSTRAP            = true;
    DO_PHASE_RND            = true;
    ASSIGN_BOOTSTRAP_MEAN   = true;
    ABSVALSQ                = true;
    SPECTRAL_DB             = true;
    WINDOW_LENGTH           = 0.5; % time (s)
    WINDOW_STEP_SIZE        = 0.025; % time (s)
    FREQS                   = (4:50); % frequency (Hz) %MAX IS 120Hz
    VERBOSITY_LEVEL         = 1;
    GUI_MODE                = 'nogui';
    N_PERMS_PHASE_RND       = 200; % number of null distribution samples to draw
    N_PERMS_BOOTSRAP        = 200;
    conn_mat_table          = [];
    %% (PARAMS) ============================================================ %%
    %## PREPARE DATA PARAMS
    DEF_PREPDATA = struct('VerbosityLevel',VERBOSITY_LEVEL,...
                 'SignalType',{{'Components'}},...
                 'VariableNames',[],...
                 'Detrend',{{'verb',VERBOSITY_LEVEL,'method',{'linear'},...
                        'piecewise',{'seglength',0.33,'stepsize',0.0825},...
                        'plot',false}},...
                 'NormalizeData',{{'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}}},...
                 'resetConfigs',true,...
                 'badsegments',[],...
                 'newtrials',[],...
                 'equalizetrials',false);
    DEF_ESTSELMOD = struct('modelingApproach',{{'Segmentation VAR',...
                            'algorithm',{{'Vieira-Morf'}},...
                            'winStartIdx',[],...
                            'winlen',WINDOW_LENGTH,...
                            'winstep',WINDOW_STEP_SIZE,...
                            'taperfcn','rectwin',...
                            'epochTimeLims',[],...
                            'prctWinToSample',100,...
                            'normalize',[],...
                            'detrend',{'method','linear'},...
                            'verb',VERBOSITY_LEVEL}},...
                    'morderRange',[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/4-1)],...
                    'downdate',true,...
                    'runPll',[],...
                    'icselector',{{'aic','hq'}},...
                    'winStartIdx',[],...
                    'epochTimeLims',[],...
                    'prctWinToSample',80,...
                    'plot',[],...
                    'verb',VERBOSITY_LEVEL);
    DEF_PLOTORDERCRIT = struct('conditions', {{}},    ...
                                'icselector', {DEF_ESTSELMOD.icselector},  ...
                                'minimizer',{{'min'}}, ...
                                'prclim', 90);
    %## Display Estimates for MVAR Validations PARAMS
    DEF_ESTDISPMVAR_CHK = struct('morder',[],...
            'winlen',WINDOW_LENGTH,'winstep',WINDOW_STEP_SIZE,'verb',VERBOSITY_LEVEL);
    %## Estimate & Validate MVAR Params
    DEF_ESTVALMVAR = struct('checkWhiteness',{{'alpha',0.05,...
                                     'statcorrection','none',...
                                     'numAcfLags',50,...
                                     'whitenessCriteria',{'Ljung-Box','ACF','Box-Pierce','Li-McLeod'},...
                                     'winStartIdx',[],...
                                     'prctWinToSample',100,...
                                     'verb',VERBOSITY_LEVEL}}, ...
                             'checkResidualVariance',{{'alpha',0.05,...
                                     'statcorrection','none',...
                                     'numAcfLags',50,...
                                     'whitenessCriteria',{{}},...
                                     'winStartIdx',[],...
                                     'prctWinToSample',100,...
                                     'verb',VERBOSITY_LEVEL}}, ...
                             'checkConsistency',{{'winStartIdx',[],...
                                     'prctWinToSample',100,...
                                     'Nr',[],...
                                     'donorm',0,...
                                     'nlags',[],...
                                     'verb',VERBOSITY_LEVEL}}, ...
                             'checkStability',{{'winStartIdx',[],...
                                     'prctWinToSample',100,...
                                     'verb',VERBOSITY_LEVEL}},     ...
                             'winStartIdx',[],      ...
                             'verb',VERBOSITY_LEVEL,...
                             'plot',false);
    %## CONNECTIVITY ESTIMATION PARAMS
    DEF_ESTFITMVAR = struct('connmethods',{'dDTF08','S'}, ...
                'absvalsq',ABSVALSQ,           ...
                'spectraldecibels',SPECTRAL_DB,   ...
                'freqs',FREQS,        ...
                'verb',VERBOSITY_LEVEL);
    %## SURROGATE STATISTICS PARAMS
    %- BOOTSTRAP
    DEF_STAT_BS_CFG = struct('mode',struct('arg_direct',1,...
                                            'nperms',N_PERMS_BOOTSRAP,...
                                            'arg_selection','Bootstrap',...
                                            'saveTrialIdx',false),...
                               'modelingApproach',[],...
                               'connectivityModeling',[],...
                               'verb',1);
    %- PHASE RANDOMIZATION
    DEF_STAT_PR_CFG = struct('mode',struct('arg_direct',1,...
                                            'nperms',N_PERMS_PHASE_RND,...
                                            'arg_selection','PhaseRand'),...
                               'modelingApproach',[],...
                               'connectivityModeling',[],...
                               'verb',1);
    %% (PARSE) ============================================================= %%
    %## VALIDATION FUNCTIONS
    sd_validfunc = (@(x) exist(x,'dir'));
    dbs_validFcn = @(x) assert(islogical(x),'Value must be (true/false). Determines whether a bootstrapped distribution will be created.');
    dpr_validFcn = @(x) assert(islogical(x),'Value must be (true/false). Determines whether a phase randomized distribution will be created.');
    %##
    p = inputParser;
    %## REQUIRED
    addRequired(p,'ALLEEG',@isstruct)
    addRequired(p,'conn_components',@isnumeric)
    addRequired(p,'save_dir',sd_validfunc)
    %## PARAMETER
    %-PREPDATA_CFG
    addParameter(p,'PREPDATA',DEF_PREPDATA,@(x) validate_struct(x,DEF_PREPDATA));
    %-
    addParameter(p,'ESTSELMOD',DEF_ESTSELMOD,@(x) validate_struct(x,DEF_ESTSELMOD));
    addParameter(p,'PLOTORDERCRIT',DEF_PLOTORDERCRIT,@(x) validate_struct(x,DEF_PLOTORDERCRIT));
    %-
    addParameter(p,'ESTVALMVAR',DEF_ESTVALMVAR,@(x) validate_struct(x,DEF_ESTVALMVAR));
    %-
    addParameter(p,'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,@(x) validate_struct(x,DEF_ESTDISPMVAR_CHK));
    addParameter(p,'ESTFITMVAR',DEF_ESTFITMVAR,@(x) validate_struct(x,DEF_ESTFITMVAR));
    addParameter(p,'STAT_BS_CFG',DEF_STAT_BS_CFG,@(x) validate_struct(x,DEF_STAT_BS_CFG));
    addParameter(p,'STAT_PR_CFG',DEF_STAT_PR_CFG,@(x) validate_struct(x,DEF_STAT_PR_CFG));
    %-
    addParameter(p,'DO_PHASE_RND',DO_PHASE_RND,dpr_validFcn);
    addParameter(p,'DO_BOOTSTRAP',DO_BOOTSTRAP,dbs_validFcn);
    parse(p,ALLEEG,conn_components,save_dir,varargin{:});
    %## SET DEFAULTS
    %- Optional
    PREPDATA       = p.Results.PREPDATA;
    ESTSELMOD       = p.Results.ESTSELMOD;
    PLOTORDERCRIT = p.Results.PLOTORDERCRIT;
    %-
    % ESTVALMVAR_CW      = p.Results.ESTVALMVAR_CW;
    % ESTVALMVAR_CRV       = p.Results.ESTVALMVAR_CRV;
    % ESTVALMVAR_CC       = p.Results.ESTVALMVAR_CC;
    % ESTVALMVAR_CS       = p.Results.ESTVALMVAR_CS;
    ESTVALMVAR = p.Results.ESTVALMVAR;
    %-
    ESTDISPMVAR_CHK       = p.Results.ESTDISPMVAR_CHK;
    ESTFITMVAR      = p.Results.ESTFITMVAR;
    STAT_BS_CFG       = p.Results.STAT_BS_CFG;
    STAT_PR_CFG       = p.Results.STAT_PR_CFG;
    DO_PHASE_RND       = p.Results.DO_PHASE_RND;
    DO_BOOTSTRAP       = p.Results.DO_BOOTSTRAP;
    %##
    PREPDATA = set_defaults_struct(PREPDATA,DEF_PREPDATA);
    ESTSELMOD = set_defaults_struct(ESTSELMOD,DEF_ESTSELMOD);
    PLOTORDERCRIT = set_defaults_struct(PLOTORDERCRIT,DEF_PLOTORDERCRIT);
    ESTVALMVAR = set_defaults_struct(ESTVALMVAR,DEF_ESTVALMVAR);
    ESTDISPMVAR_CHK = set_defaults_struct(ESTDISPMVAR_CHK,DEF_ESTDISPMVAR_CHK);
    ESTFITMVAR = set_defaults_struct(ESTFITMVAR,DEF_ESTFITMVAR);
    STAT_BS_CFG = set_defaults_struct(STAT_BS_CFG,DEF_STAT_BS_CFG);
    STAT_PR_CFG = set_defaults_struct(STAT_PR_CFG,DEF_STAT_PR_CFG);
    %##
    % if isempty(DEF_ESTSELMOD.modelingApproach.algorithm.morder)
    %     ESTSELMOD.modelingApproach.algorithm.morder = [1,ceil(ALLEEG(1).srate*WINDOW_LENGTH)/2-1];
    % end
    if isempty(ESTSELMOD.morderRange)
        ESTSELMOD.morderRange=[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/4-1)];
    end
    % if isempty(ESTSELMOD.modelingApproach) && ~isempty(ESTDISPMVAR_CHK.morder)
    %     ESTSELMOD.modelingApproach = ESTDISPMVAR_CHK.morder;
    % else
    %     ESTSELMOD.morder = ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1);
    % end
    
    %% GENERATE CONNECTIVITY MEASURES
    fprintf(1,'\n==== GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s ====\n',ALLEEG(1).subject);
    %## ASSIGN PATH FOR SIFT
    %- make sure path is in right format and make sure it exists
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
    %## Connectivity
    fprintf('%s) Processing componets:\n',ALLEEG(1).subject)
    fprintf('%i,',conn_components'); fprintf('\n');
    %- exit function if there are not enough components
    if length(conn_components) < 2
        return;
    end
    %- select components from EEG
    if length(conn_components) ~= size(ALLEEG(1).icaweights,1)
        for cond_i = 1:length(ALLEEG)
            TMP_EEG = ALLEEG(cond_i);
            TMP_EEG = pop_subcomp(TMP_EEG,sort(conn_components),0,1);
            %- Recalculate ICA Matrices && Book Keeping
            if isempty(TMP_EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',TMP_EEG.subject);
                TMP_EEG.icaact = (TMP_EEG.icaweights*TMP_EEG.icasphere)*TMP_EEG.data(TMP_EEG.icachansind,:);
                TMP_EEG.icaact = reshape(TMP_EEG.icaact,size(TMP_EEG.icaact,1),TMP_EEG.pnts,TMP_EEG.trials);
            end
            ALLEEG(cond_i) = TMP_EEG;
        end
    end
    %% (MAIN CONNECTIVITY PIPELINE) ======================================== %%
    %## STEP 3: Pre-process the data
    fprintf('===================================================\n');
    disp('PRE-PROCESSISNG DATA');
    fprintf('===================================================\n');
    %- (NOTE FROM EXAMPLE SIFT SCRIPT) No piecewise detrending based on conversation with Andreas Widmann.
    % convert list of components to cell array of strings
    comp_names = [];
    for j = 1:length(conn_components) 
        comp_names = [comp_names, {num2str(conn_components(j))}];
    end
    
    PREPDATA.VariableNames = comp_names;
    cfg = struct2args(PREPDATA);
    [ALLEEG] = pop_pre_prepData(ALLEEG,GUI_MODE,...
                 cfg{:});
    
    %% STEP 4: Identify the optimal model order
    fprintf('===================================================\n');
    disp('MODEL ORDER IDENTIFICATION');
    fprintf('===================================================\n');
    %- NOTE: Here we compute various model order selection criteria for varying model
    % orders (e.g. 1 to 30) and visualize the results
    %## OPTION 1
    % cfg = struct2args(ESTSELMOD);
    % ALLEEG(cond_i) = pop_est_selModelOrder(ALLEEG(cond_i),GUI_MODE,cfg{:});
    %## OPTION 2
    for cond_i=1:length(ALLEEG)
        %* calculate the information criteria
        [ALLEEG(cond_i).CAT.IC,cfg] = est_selModelOrder(ALLEEG(cond_i),ESTSELMOD);
        if isempty(ALLEEG(cond_i).CAT.IC)
            % use canceled
            fprintf('ERROR. Model fittings didn''t produce a viable model\n...')
            return;
        end
        if ~isempty(cfg)
            %* store the configuration structure
            ALLEEG(cond_i).CAT.configs.('est_selModelOrder') = cfg;
        end
    end
    modFuncName = ['pop_' ALLEEG(1).CAT.IC.modelFitting.modelingFuncName];
    fprintf('Running %s for subject %s...\n',modFuncName,ALLEEG(1).subject)
    ALLEEG = feval(modFuncName, ALLEEG,0,ALLEEG(1).CAT.IC.modelFitting.modelingArguments);
    %## (PLOT) ============================================================= %%
    tmp_morder = zeros(1,length(ALLEEG));
    % cfg = struct2args(PLOTORDERCRIT);
    for cond_i = 1:length(ALLEEG)
        fprintf('%s) Plotting Validations\n',ALLEEG(cond_i).subject);
        handles = vis_plotOrderCriteria(ALLEEG(cond_i).CAT.IC,PLOTORDERCRIT);
        tmp_morder(cond_i) = ceil(mean(ALLEEG(cond_i).CAT.IC.hq.popt));
        saveas(handles,[save_dir filesep sprintf('%s_%i_orderResults.fig',ALLEEG(cond_i).subject,cond_i)]);
        close(handles);
    end
    %##
    if isempty(ESTDISPMVAR_CHK.morder)
        ESTDISPMVAR_CHK.morder = ceil(mean(tmp_morder));
    end
    fprintf('\n\n');
    %% STEP 5: Fit the VAR model
    fprintf('===================================================\n');
    disp('MODEL FITTING');
    fprintf('===================================================\n');
    fprintf('\n');
    %- Here we can check that our selected parameters make sense
    % cfg = struct2args(ESTDISPMVAR_CHK);
    for cond_i = 1:length(ALLEEG)
        fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n\n',ALLEEG(cond_i).condition);
        est_dispMVARParamCheck(ALLEEG(cond_i),ESTDISPMVAR_CHK)
    end
    %- Once we have identified our optimal model order, we can fit our VAR model.
    % Note that EEG.CAT.MODEL now contains the model structure with
    % coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
    % self-evident information. 
    for cond_i = 1:length(ALLEEG)
        [ALLEEG(cond_i)] = pop_est_fitMVAR(ALLEEG(cond_i),GUI_MODE,...
                ALLEEG(cond_i).CAT.configs.est_selModelOrder.modelingApproach,...
                'ModelOrder',ESTDISPMVAR_CHK.morder);
    end
    %% STEP 6: Validate the fitted model
    fprintf('===================================================\n');
    disp('MODEL VALIDATION');
    fprintf('===================================================\n');
    % Here we assess the quality of the fit of our model w.r.t. the data. This
    % step can be slow. We can obtain statistics for residual whiteness, percent consistency, and
    % model stability...
    cfg = struct2args(ESTVALMVAR);
    [ALLEEG] = pop_est_validateMVAR(ALLEEG,GUI_MODE,...
                                cfg{:});
    
    for cond_i = 1:length(ALLEEG)
        % ... and then plot the results 
        handles = vis_plotModelValidation({ALLEEG(cond_i).CAT.VALIDATION.whitestats}, ...
                                             {ALLEEG(cond_i).CAT.VALIDATION.PCstats},         ...
                                             {ALLEEG(cond_i).CAT.VALIDATION.stabilitystats});                                        
        % If you want to save this figure you can uncomment the following lines:
        saveas(handles,[save_dir filesep sprintf('%s_%i_validationResults.fig',ALLEEG(cond_i).subject,cond_i)]);
        close(handles);
    
        % To automatically determine whether our model accurately fits the data you
        % can write a few lines as follows (replace 'acf' with desired statistic):
        if ~all(ALLEEG(cond_i).CAT.VALIDATION.whitestats.acf.w)
            fprintf(1,'WARNING: Residuals are not completely white!\nModel fit may not be valid...\n');
        end
    end
    fprintf('\n\n');
    %% STEP 7: Compute Connectivity
    %- NOTE: Next we will compute various dynamical quantities, including connectivity,
    % from the fitted VAR model.
    fprintf('===================================================\n');
    fprintf('CONNECTIVITY ESTIMATION\n');
    fprintf('===================================================\n');
    cfg = struct2args(ESTFITMVAR);
    ALLEEG = pop_est_mvarConnectivity(ALLEEG,GUI_MODE, ...
                cfg{:});
            
    disp('===================================')
    disp('DONE.')
    disp('===================================')
    fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s ====\n',ALLEEG(1).subject);          
    %% ===================================================================== %%
    %## STEP 5.a) (BOOTSTRAPPING) GROUP STATISTICS 
    % (09/22/2022), JS, Might want to try and speed up bootstrap by
    % adapting stat_surrogateGen.m to use parfor for bootstrapping... If
    % possible? doesn't seem built well in the first place, so maybe?
    % (10/27/2022), JS, Revist above note again!
    % (12/7/2022), JS, need to update this boostrapping to include ALLEEG
    if DO_BOOTSTRAP
        fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
        for cond_i=1:length(ALLEEG)
            %- calculate BootStrap distribution
            %- clear PConn
            ALLEEG(cond_i).CAT.PConn  = [];
            cfg = struct2args(STAT_BS_CFG);
            [PConn,~] = stat_surrogateGen('ALLEEG',ALLEEG(cond_i),cfg{:});
            ALLEEG(cond_i).CAT.PConn = PConn;
            %- save BootStrap distribution 
            bootstrap_dist = ALLEEG(cond_i).CAT.PConn;
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            par_save(bootstrap_dist,ALLEEG(cond_i).filepath,fName,'_BootStrap');
        end
        %- assign mean of bootstrap as Conn value
        if ASSIGN_BOOTSTRAP_MEAN
            for cond_i = 1:length(ALLEEG)
                ALLEEG(cond_i).CAT.Conn = stat_getDistribMean(ALLEEG(cond_i).CAT.PConn);
            end
        end 
        for cond_i = 1:length(ALLEEG)
            %- clear bootstrap calculation
            ALLEEG(cond_i).CAT.PConn = [];
        end    
        fprintf('\n==== DONE: CALCULATING BOOTSTRAP MEASURES ====\n')
    end
    %% ===================================================================== %%
    %## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION    
    % see. stat_surrogateGen
    % see. stat_surrogateStats
    if DO_PHASE_RND
        fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
        %## (1) ALTERNATIVE CODE 
        for cond_i=1:length(ALLEEG)
            %- Generate Phase Randomized Distribution
            fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
            %- clear PConn
            ALLEEG(cond_i).CAT.PConn  = [];
            STAT_PR_CFG.modelingApproach = ALLEEG(cond_i).CAT.configs.est_fitMVAR;
            STAT_PR_CFG.connectivityModeling = ALLEEG(cond_i).CAT.configs.est_mvarConnectivity;
            cfg = struct2args(STAT_PR_CFG);
            %- FEVAL
            [PConn,~] = stat_surrogateGen('ALLEEG',ALLEEG(cond_i),cfg{:});
            ALLEEG(cond_i).CAT.PConn = PConn;
            %- Save Phase randomized distribution
            phasernd_dist = ALLEEG(cond_i).CAT.PConn;
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            par_save(phasernd_dist,ALLEEG(cond_i).filepath,fName,'_PhaseRnd');
            fprintf('done.\n')
        end
        clear TMP_EEG
        fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    end
    
    %% ===================================================================== %%
    %## STEP 6) CONNECTIVITY MATRICES & VISUALS
    %## TABLE VARS
    t_fPaths = cell(length(ALLEEG)*length(conn_estimators),1);
    t_fNames = cell(length(ALLEEG)*length(conn_estimators),1);
    t_conn_comps = cell(length(ALLEEG)*length(conn_estimators),1);
    t_conn_meas = cell(length(ALLEEG)*length(conn_estimators),1);
    t_cond_char = cell(length(ALLEEG)*length(conn_estimators),1);
    disp(ALLEEG);
    %## Generate Connectivity Matrix
    cnt = 1;
    for conn_i = 1:length(conn_estimators)
        for cond_i = 1:length(ALLEEG)
            disp(ALLEEG(cond_i).CAT.Conn);
            t_conn_comps{cnt} = conn_components;
            t_fPaths{cnt} = ALLEEG(cond_i).filepath;
            t_fNames{cnt} = ALLEEG(cond_i).filename;
            t_conn_meas{cnt} = conn_estimators{conn_i};
            t_cond_char{cnt} = ALLEEG(cond_i).condition;
            cnt = cnt + 1;
        end
    end
    %- create table
    conn_mat_table = table(t_fPaths,t_fNames,t_conn_comps,t_conn_meas,t_cond_char);
    %% ===================================================================== %%
    %## STEP 7) WRAP UP 
    %- REMOVE FIELDS && SAVE
    for cond_i = 1:length(ALLEEG)
        %- clear STATS & PCONN for space saving.
        if isfield(ALLEEG(cond_i).CAT,'Stats')
            ALLEEG(cond_i).CAT = rmfield(ALLEEG(cond_i).CAT,'Stats');
        end
        if isfield(ALLEEG(cond_i).CAT,'PConn')
            ALLEEG(cond_i).CAT = rmfield(ALLEEG(cond_i).CAT,'PConn');
        end
    end
    fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    
    %% DONE
    fprintf('DONE: %s\n',ALLEEG(1).subject);
    %## TIME
    toc
end

%% ===================================================================== %%
function [b] = validate_struct(x,DEFAULT_STRUCT)
    b = false;
    struct_name = inputname(2);
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field names
    chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
    if ~all(chk)
        fprintf(2,'\nFields for struct do not match for %s\n',struct_name);
        return
    end
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end
%% ===================================================================== %%
function [struct_out] = set_defaults_struct(x,DEFAULT_STRUCT)
    struct_out = x;
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        if isempty(vals1{ind})
            struct_out.(fs1{ind}) = DEFAULT_STRUCT.(fs2{ind});
        end
    end
end
%% ===================================================================== %%
function [args] = struct2args(struct)
    %EEGLAB_STRUCT2ARGS Summary of this function goes here
    %   Detailed explanation goes here
    %## Define Parser
    p = inputParser;
    %## REQUIRED
    addRequired(p,'struct',@isstruct);
    parse(p,struct);
    %##
    fn = fieldnames(struct);
    sc = struct2cell(struct);
    args = cell(length(fn)+length(sc),1);
    cnt = 1;
    for i = 1:length(fn)
        args{cnt} = fn{i};
        if isempty(sc{i})
            args{cnt+1} = [];
        else
            args{cnt+1} = sc{i};
        end
        cnt = cnt + 2;
    end
end



