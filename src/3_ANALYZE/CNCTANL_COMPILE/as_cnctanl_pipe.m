% function [error_code] = as_cnctanl_pipe(ALLEEG,conn_components,conn_estimators,save_dir,varargin)
function [error_code] = as_cnctanl_pipe(JSON_Params)

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
%       JSON Format is as follows:
            % {
            %     "ALLEEG": {},

            %     "conn_components" : [],
            %     "conn_estimators" : [],
            %     "save_dir" : "",

            %     "optional" : {
            %         "DO_PHASE_RND" : null,
            %         "DO_BOOTSTRAP" : null,

            %         "FREQS" : null,
            %         "FREQ_BANDS" : null,

            %         "WINDOW_LENGTH" : null,
            %         "WINDOW_STEP_SIZE" : null,

            %         "GUI_MODE" : null,
            %         "VERBOSITY_LEVEL" : null,

            %         "ESTSELMOD_CFG" : null,
            %         "ESTVALMVAR_CW" : null,
            %         "ESTVALMVAR_CRV" : null,
            %         "ESTVALMVAR_CC" : null,
            %         "ESTVALMVAR_CS" : null,

            %         "MORDER" :null
            %     }
            % }
            
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
        'error. on working_dir %s\n',...
        'stack. %s\n'],e.identifier,e.message,working_dir,getReport(e));
end
%% (DEFINE DEFAULTS) =================================================== %%
%## HIGH LEVEL VARS
%- connectivity statistics & surrogates params
%*
fid = fopen(JSON_Params); 
raw = fread(fid,inf); 
str = char(raw');
fclose(fid); 
params= jsondecode(str);
alleeg_fpaths = params.ALLEEG;
conn_components = params.conn_components;
conn_estimators = params.conn_estimators;
%% (READ ALLEEG)
% alleeg_fpaths = jsonfile.whatever;
ALLEEG = cell(length(alleeg_fpaths),1);
% parfor subj_i = 1:length(alleeg_fpaths)
for subj_i = 1:length(alleeg_fpaths)
    tmp = strsplit(alleeg_fpaths{subj_i},filesep);
    fName = tmp{end};
    fPath = strjoin(tmp(1:end-1),filesep);
    %- load EEG
    % EEG = pop_loadset('filepath',fPath,'filename',fName);
    fprintf('Loading EEG...\n');
    EEG = load('-mat', [fPath filesep fName]);
    fid_eeg = fopen([fPath filesep EEG.data], 'r', 'ieee-le');
    data = fread(fid_eeg, [EEG.trials*EEG.pnts EEG.nbchan], 'float32')';
    EEG.data = data;
    EEG = eeg_checkset(EEG,'ica');
    fclose(fid_eeg);
    %- Recalculate ICA Matrices && Book Keeping
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    ALLEEG{subj_i} = EEG;
end
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
%%
error_code = 0;
DO_BOOTSTRAP = false;
%* boolean for phase randomization generation
DO_PHASE_RND = true;
%*
STAT_ALPHA = 0.01;
%*
% conn_components = (1:size(ALLEEG(1).icaweights,1))';
% conn_estimators = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
WINDOW_LENGTH = 0.5; % time (s)
WINDOW_STEP_SIZE = 0.025; % time (s)
FREQS = (1:200); % frequency (Hz)
VERBOSITY_LEVEL = 1;
GUI_MODE = 'nogui';
FREQ_BANDS = {FREQS;1:7;7:12;12:28;28:48;48:60};
N_PERMS_PHASE_RND = 2000; % number of null distribution samples to draw
N_PERMS_BOOTSRAP = 2000;
MORDER = [];
%% POP_PRE_PREPDATA
PREPDATA_SIGNALTYPE = {'Components'};
PREPDATA_DETRED = {'verb',VERBOSITY_LEVEL,'method',{'linear'}};
PREPDATA_NORMDATA = {'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}};
%% POP_EST_SELMODORDER
%- structs
ESTSELMOD_DETREN = struct('arg_direct',1,'method','constant','arg_selection',1);
ESTSELMOD_ALGORITHM = struct('arg_direct',1,...
    'morder',[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1)],...
    'arg_selection',{'Vieira-Morf'});
ESTSELMOD_MODELINGAPPROACH = struct('arg_direct',1,...
    'algorithm',ESTSELMOD_ALGORITHM,...
    'morder',ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1),...
    'winStartIdx',[],...
    'winlen',WINDOW_LENGTH,...
    'winstep',WINDOW_STEP_SIZE,...
    'taperfcn','rectwin',...
    'epochTimeLims',[],...
    'prctWinToSample',100,...
    'normalize',[],...
    'detrend',ESTSELMOD_DETREN,...
    'verb',VERBOSITY_LEVEL,...
    'timer',0,...
    'setArgDirectMode',1,...
    'arg_selection','Segmentation VAR');
ESTSELMOD_RUNPLL = struct('arg_direct',1,'arg_selection',0);
ESTSELMOD_PLOT = struct('arg_direct',1,'arg_selection',0);
DEF_ESTSELMOD = struct('verb',VERBOSITY_LEVEL,...
    'modelingApproach',ESTSELMOD_MODELINGAPPROACH,...
    'morderRange',[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1)],...
    'downdate',1,...
    'runPll',ESTSELMOD_RUNPLL,...
    'icselector',{{'aic','hq'}},...
    'winStartIdx',[],...
    'epochTimeLims',[],...
    'prctWinToSample',80,...
    'plot',ESTSELMOD_PLOT);
%% POP_EST_VALIDATE_MVAR
%## cells
ESTVALMVAR_CHECKWHITENESS = {'alpha',0.05 ...
    'statcorrection','none',...
    'numAcfLags',50,...
    'whitenessCriteria',{'Ljung-Box' 'ACF' 'Box-Pierce' 'Li-McLeod'} ...
    'winStartIdx',[],...
    'prctWinToSample',100,...
    'verb',VERBOSITY_LEVEL};
ESTVALMVAR_CHECKRESIDUALVARIANCE = {'alpha',0.05,...
    'statcorrection','none',...
    'numAcfLags',50,...
    'whitenessCriteria',{},...
    'winStartIdx',[],...
    'prctWinToSample',100,...
    'verb',VERBOSITY_LEVEL};
ESTVALMVAR_CHECKCONSISTENCY = {'winStartIdx',[],...
    'prctWinToSample',100,...
    'Nr',[],...
    'donorm',0,...
    'nlags',[],...
    'verb',VERBOSITY_LEVEL};
ESTVALMVAR_CHECKSTABILITY = {'winStartIdx',[],...
    'prctWinToSample',100,...
    'verb',VERBOSITY_LEVEL};
%% POP_EST_MVARCONNECTIVITY
ABSVALSQ        = false;
SPECTRAL_DB     = false;
%% STAT_SURROGATEGEN
%- BOOTSTRAP
STAT_BOOTSTRAP_CFG = [];
STAT_BOOTSTRAP_CFG.verb = 1;
STAT_BOOTSTRAP_CFG.mode = {'Bootstrap','nperms',N_PERMS_BOOTSRAP,'saveTrialIdx',false};
%- PHASE RANDOMIZATION
STAT_PHASERND_CFG = [];
STAT_PHASERND_CFG.mode.arg_direct = 1;
STAT_PHASERND_CFG.mode.nperms = N_PERMS_PHASE_RND;
STAT_PHASERND_CFG.mode.arg_selection = 'PhaseRand';
STAT_PHASERND_CFG.modelingApproach = [];%ALLEEG(trial_i).CAT.configs.est_fitMVAR;
STAT_PHASERND_CFG.connectivityModeling = [];%ALLEEG(trial_i).CAT.configs.est_mvarConnectivity;
STAT_PHASERND_CFG.verb = 1;
%% (VALIDATION OPERATIONS) ============================================= %%
%- figure savepath
% sd_validfunc = (@(x) exist(x,'dir'));
% dbs_validFcn = @(x) assert(islogical(x),'Value must be (true/false). Determines whether a bootstrapped distribution will be created.');
% dpr_validFcn = @(x) assert(islogical(x),'Value must be (true/false). Determines whether a phase randomized distribution will be created.');
% esmc_validate = @(x) assert(isempty(x) || isstruct(x),'Value must be EMPTY or STRUCT. See pop_est_selModOrder.m for more details.');
% evmc_validate = @(x) assert(isempty(x) || iscell(x),'Value must be EMPTY or CELL of CHARS & VALUES. See. pop_est_validateMVAR.m for more details.');
% mo_validate = @(x) assert(isempty(x) || isnumeric(x),'Value must be EMPTY or NUMERIC. Determines the model order for MVAR.');
%% (PARSE) ============================================================= %%


% %## REQUIRED
% addRequired(p,'ALLEEG',@isstruct)
% addRequired(p,'conn_components',@isnumeric)
% addRequired(p,'conn_estimators',@iscell)
% addRequired(p,'save_dir',sd_validfunc)
% %## PARAMETER
% addParameter(p,'DO_PHASE_RND', DO_PHASE_RND,dpr_validFcn);
% addParameter(p,'DO_BOOTSTRAP', DO_BOOTSTRAP,dbs_validFcn);
% addParameter(p,'FREQS',FREQS,@isnumeric)
% addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@isnumeric) % number of current set in previous load
% addParameter(p,'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,@isnumeric) % MIM folder overwrite in special cases
% addParameter(p,'GUI_MODE',GUI_MODE,@ischar)
% addParameter(p,'VERBOSITY_LEVEL',VERBOSITY_LEVEL,@isnumeric)
% addParameter(p,'ESTSELMOD_CFG',DEF_ESTSELMOD,esmc_validate);
% addParameter(p,'ESTVALMVAR_CW',ESTVALMVAR_CHECKWHITENESS,evmc_validate);
% addParameter(p,'ESTVALMVAR_CRV',ESTVALMVAR_CHECKRESIDUALVARIANCE,evmc_validate);
% addParameter(p,'ESTVALMVAR_CC',ESTVALMVAR_CHECKCONSISTENCY,evmc_validate);
% addParameter(p,'ESTVALMVAR_CS',ESTVALMVAR_CHECKSTABILITY,evmc_validate);
% addParameter(p,'FREQ_BANDS',FREQ_BANDS,@iscell);
% addParameter(p,'MORDER',MORDER,mo_validate);
% parse(p,ALLEEG,conn_components,conn_estimators,save_dir,varargin{:});

%## SET DEFAULTS
%- Optional
WINDOW_LENGTH       = params.optional.WINDOW_LENGTH;               % sliding window length in seconds
WINDOW_STEP_SIZE    = params.optional.WINDOW_STEP_SIZE;             % sliding window step size in seconds
GUI_MODE            = params.optional.GUI_MODE;                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = params.optional.VERBOSITY_LEVEL;              % Verbosity Level (0=no/minimal output, 2=graphical output)
FREQS               = params.optional.FREQS;
ESTSELMOD_CFG       = params.optional.ESTSELMOD_CFG;
ESTVALMVAR_CW       = params.optional.ESTVALMVAR_CW;
ESTVALMVAR_CRV      = params.optional.ESTVALMVAR_CRV;
ESTVALMVAR_CC       = params.optional.ESTVALMVAR_CC;
ESTVALMVAR_CS       = params.optional.ESTVALMVAR_CS;
MORDER              = params.optional.MORDER;
%- param mods
if isempty(VERBOSITY_LEVEL)
end
%*
if isempty(ESTSELMOD_CFG)
    ESTSELMOD_CFG = DEF_ESTSELMOD;
end
%*
if isempty(ESTVALMVAR_CW)
    ESTVALMVAR_CW = ESTVALMVAR_CHECKWHITENESS;
end
if isempty(ESTVALMVAR_CRV)
    ESTVALMVAR_CRV = ESTVALMVAR_CHECKRESIDUALVARIANCE;
end
if isempty(ESTVALMVAR_CC)
    ESTVALMVAR_CC = ESTVALMVAR_CHECKCONSISTENCY;
end
if isempty(ESTVALMVAR_CS)
    ESTVALMVAR_CS = ESTVALMVAR_CHECKSTABILITY;
end


%% ===================================================================== %%
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

%% (MAIN CONNECTIVITY PIPELINE) ========================================= %%

%## STEP 3: Pre-process the data
fprintf('===================================================\n');
disp('PRE-PROCESSISNG DATA');
fprintf('===================================================\n');
%- NOTE:
%- (NOTE FROM EXAMPLE SIFT SCRIPT) No piecewise detrending based on conversation with Andreas Widmann.
% convert list of components to cell array of strings
comp_names = [];
for j = 1:length(conn_components)
    comp_names = [comp_names, {num2str(conn_components(j))}];
end

%- ALTERNATIVE
% input_cell = {'VerbosityLevel',VERBOSITY_LEVEL,...
%              'SignalType',PREPDATA_SIGNALTYPE,...
%              'VariableNames',comp_names,...
%              'Detrend',PREPDATA_DETRED,...
%              'NormalizeData',PREPDATA_NORMDATA,...
%              'resetConfigs',true,...
%              'badsegments',[],...
%              'newtrials',[],...
%              'equalizetrials',false};
[ALLEEG] = pop_pre_prepData(ALLEEG,'nogui',...
    'VerbosityLevel',VERBOSITY_LEVEL,...
    'SignalType',PREPDATA_SIGNALTYPE,...
    'VariableNames',comp_names,...
    'Detrend',PREPDATA_DETRED,...
    'NormalizeData',PREPDATA_NORMDATA,...
    'resetConfigs',true,...
    'badsegments',[],...
    'newtrials',[],...
    'equalizetrials',false);

%% STEP 4: Identify the optimal model order
fprintf('===================================================\n');
disp('MODEL ORDER IDENTIFICATION');
fprintf('===================================================\n');
%- NOTE: Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results
% Open the model-fitting GUI for model fitting.
% Once model is fit results will be return in EEG structure
for cond_i=1:length(ALLEEG)
    %* calculate the information criteria
    [ALLEEG(cond_i).CAT.IC,cfg] = est_selModelOrder('EEG',ALLEEG(cond_i),ESTSELMOD_CFG);
    if ~isempty(cfg)
        %* store the configuration structure
        ALLEEG(cond_i).CAT.configs.('est_selModelOrder') = cfg;
    end
end
% (08/12/2023) JS, this may be buggy on HiPerGator/Hypercomputing systems due to the
% gui initiation. buggy because of feval statement: "
% modFuncName = ['pop_' ALLEEG(1).CAT.IC.modelFitting.modelingFuncName];
% ALLEEG = feval(modFuncName, ALLEEG,0,ALLEEG(1).CAT.IC.modelFitting.modelingArguments);
% (08/16/2023) JS, this function can only be ran on a single EEG dataset,
% if using ALLEEG (i.e., a set of EEGs with all conditions for a subject)
% this won't work.

%## (PLOT)
tmp_morder = zeros(1,length(ALLEEG));
for cond_i = 1:length(ALLEEG)
    fprintf('%s) Plotting Validations',ALLEEG(cond_i).subject);
    handles = vis_plotOrderCriteria(ALLEEG(cond_i).CAT.IC,'conditions', [],    ...
        'icselector', ESTSELMOD_CFG.icselector,  ...
        'minimizer', 'min', ...
        'prclim', 90);
    tmp_morder(cond_i) = ceil(mean(ALLEEG(cond_i).CAT.IC.hq.popt));
    saveas(handles,[save_dir filesep sprintf('%s_%i_orderResults.fig',ALLEEG(cond_i).subject,cond_i)]);
    close(handles);
end
if isempty(MORDER)
    model_order = ceil(mean(tmp_morder));
else
    model_order = MORDER;
end

% Finally, we can automatically select the model order which minimizes one
% of the criteria (or you can set this manually based on above figure)
% ModelOrder = ceil(mean(ModelOrder)); % (05/16/2022) Probably a better way of selecting model order

% As an alternative to using the minimum of the selection criteria over
% model order, you can find the "elbow" in the plot of model order versus
% selection criterion value. This is useful in cases where the selection
% criterion does not have a clear minimum. For example, the lines below
% plot and select the elbow location (averaged across windows) for the AIC
% criterion
%
% vis_plotOrderCriteria(EEG(1).CAT.IC,{},{},'elbow');
% ModelOrder = ceil(mean(EEG(1).CAT.IC.aic.pelbow));
fprintf('\n\n');
%% STEP 5: Fit the VAR model
fprintf('===================================================\n');
disp('MODEL FITTING');
fprintf('===================================================\n');
fprintf('\n');
%- Here we can check that our selected parameters make sense
for cond_i = 1:length(ALLEEG)
    fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n\n',ALLEEG(cond_i).condition);
    est_dispMVARParamCheck(ALLEEG(cond_i),struct('morder',model_order','winlen',WINDOW_LENGTH,'winstep',WINDOW_STEP_SIZE,'verb',VERBOSITY_LEVEL))
end
%- Once we have identified our optimal model order, we can fit our VAR model.
% Fit a model using the options specifed for model order selection.
% Note that EEG.CAT.MODEL now contains the model structure with
% coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
% self-evident information.
[ALLEEG] = pop_est_fitMVAR(ALLEEG,GUI_MODE,...
    ALLEEG(1).CAT.configs.est_selModelOrder.modelingApproach,...
    'ModelOrder',model_order);

%- Alternately, we can fit the VAR parameters using a Kalman filter (see
% doc est_fitMVARKalman for more info on arguments)
% EEG.CAT.MODEL = est_fitMVARKalman(EEG,0,'updatecoeff',0.0005,'updatemode',2,'morder',ModelOrder,'verb',2,'downsampleFactor',50);
fprintf('\n\n');
%% STEP 6: Validate the fitted model
fprintf('===================================================\n');
disp('MODEL VALIDATION');
fprintf('===================================================\n');

% Here we assess the quality of the fit of our model w.r.t. the data. This
% step can be slow. We can obtain statistics for residual whiteness, percent consistency, and
% model stability...
[ALLEEG] = pop_est_validateMVAR(ALLEEG,GUI_MODE,...
    'checkWhiteness',ESTVALMVAR_CW, ...
    'checkResidualVariance',ESTVALMVAR_CRV, ...
    'checkConsistency',ESTVALMVAR_CC, ...
    'checkStability',ESTVALMVAR_CS,     ...
    'winStartIdx',[],      ...
    'verb',VERBOSITY_LEVEL,...
    'plot',false);

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
    %
    if ~all(ALLEEG(cond_i).CAT.VALIDATION.whitestats.acf.w)
        fprintf(1,'WARNING: Residuals are not completely white!\n');
    end
end
fprintf('\n\n');

%% STEP 7: Compute Connectivity
%- NOTE: Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model. We can compute these for a range of
% frequencies (here 1-40 Hz). See 'doc est_mvarConnectivity' for a complete
% list of available connectivity and spectral estimators.

fprintf('===================================================\n');
fprintf('CONNECTIVITY ESTIMATION\n');
fprintf('===================================================\n');

ALLEEG = pop_est_mvarConnectivity(ALLEEG,GUI_MODE, ...
    'connmethods', conn_estimators, ...
    'absvalsq',ABSVALSQ,           ...
    'spectraldecibels',SPECTRAL_DB,   ...
    'freqs',FREQS,        ...
    'verb',VERBOSITY_LEVEL);

disp('===================================')
disp('DONE.')
disp('===================================')
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s ====\n',ALLEEG(1).subject);
%% ===================================================================== %%
%## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION
% see. stat_surrogateGen
% see. stat_surrogateStats
if DO_PHASE_RND
    fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    %## (1) ALTERNATIVE CODE
    for trial_i=1:length(ALLEEG)
        %- Generate Phase Randomized Distribution
        fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
        %         [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'PhaseRnd');
        %- clear PConn
        ALLEEG(trial_i).CAT.PConn  = [];
        STAT_PHASERND_CFG.modelingApproach = ALLEEG(trial_i).CAT.configs.est_fitMVAR;
        STAT_PHASERND_CFG.connectivityModeling = ALLEEG(trial_i).CAT.configs.est_mvarConnectivity;
        %- FEVAL
        [PConn,~] = feval(@stat_surrogateGen,'ALLEEG',ALLEEG(trial_i),STAT_PHASERND_CFG);
        ALLEEG(trial_i).CAT.PConn = PConn;
        %- Save Phase randomized distribution
        phasernd_dist = ALLEEG(trial_i).CAT.PConn;
        fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(phasernd_dist,ALLEEG(trial_i).filepath,fName,'_PhaseRnd');
        fprintf('done.\n')
        %- Conduct Nonzero Test (Need Phase Randomized Distribution)
        %         fprintf('\n==== NONZERO STATISTICS ====\n')
        %         [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'NonZero');
        %         %- Save Nonzero Stats
        %         nonzero_stats = ALLEEG(trial_i).CAT.Stats;
        %         fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        %         par_save(nonzero_stats,ALLEEG(trial_i).filepath,fName,'_NonZero');
        
        fprintf('done.\n')
    end
    clear TMP_EEG
    fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
end
%% ===================================================================== %%
%## STEP 5.a) (BOOTSTRAPPING) GROUP STATISTICS
% (09/22/2022), JS, Might want to try and speed up bootstrap by
% adapting stat_surrogateGen.m to use parfor for bootstrapping... If
% possible? doesn't seem built well in the first place, so maybe?
% (10/27/2022), JS, Revist above note again!
% (12/7/2022), JS, need to update this boostrapping to include ALLEEG
if DO_BOOTSTRAP
    fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
    for trial_i=1:length(ALLEEG)
        %- calculate BootStrap distribution
        %         ALLEEG(trial_i) = cnctanl_groupStats(ALLEEG(trial_i),'BootStrap');
        %- clear PConn
        ALLEEG(trial_i).CAT.PConn  = [];
        [PConn,~] = feval(@stat_surrogateGen,'ALLEEG',ALLEEG(trial_i),STAT_BOOTSTRAP_CFG);
        ALLEEG(trial_i).CAT.PConn = PConn;
        %- save BootStrap distribution
        bootsrap_dist = ALLEEG(trial_i).CAT.PConn;
        fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(bootsrap_dist,ALLEEG(trial_i).filepath,fName,'_BootStrap');
        %- BootStrap corrected statistics per condition
        %         fprintf('\n==== BETWEEN CONDITION STATISTICS ====\n')
        %         [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'BtwnCond');
        %         %- Save Nonzero Stats
        %         TMP_EEG = ALLEEG(trial_i).CAT.Stats;
        %         fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        %         par_save(TMP_EEG,ALLEEG(trial_i).filepath,fName,'_BtwnCond');
    end
    %- assign mean of bootstrap as Conn value
    if ASSIGN_BOOTSTRAP_MEAN
        for trial_i = 1:length(ALLEEG)
            ALLEEG(trial_i).CAT.Conn = stat_getDistribMean(ALLEEG(trial_i).CAT.PConn);
        end
    end
    for trial_i = 1:length(ALLEEG)
        %- clear bootstrap calculation
        ALLEEG(trial_i).CAT.PConn = [];
    end
    fprintf('\n==== DONE: CALCULATING BOOTSTRAP MEASURES ====\n')
end
%% ===================================================================== %%
%## STEP 6) CONNECTIVITY MATRICES & VISUALS
%## TABLE VARS
t_fPaths = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
t_fNames = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
% t_conn_methods = cell(length(ALLEEG),1);
t_conn_comps = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
t_conn_mats = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
t_conn_freqs = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
t_conn_sigs = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
t_conn_meas = cell(length(ALLEEG)*length(conn_estimators)*length(FREQ_BANDS),1);
disp(ALLEEG);
%## Generate Connectivity Matrix
cnt = 1;
for conn_i = 1:length(conn_estimators)
    for freq_i = 1:length(FREQ_BANDS)
        for trial_i = 1:length(ALLEEG)
            disp(ALLEEG(trial_i).CAT.Conn);
            %- calculate average connnectivity for each component pair across time
            [statMat,extract_sig] = gen_connMatrix(ALLEEG(trial_i),conn_estimators{conn_i},FREQ_BANDS{freq_i},STAT_ALPHA);
            %- store outputs into a tabularable ("table able") format.
            t_conn_comps{cnt} = conn_components;
            t_fPaths{cnt} = ALLEEG(trial_i).filepath;
            t_fNames{cnt} = ALLEEG(trial_i).filename;
            t_conn_freqs{cnt} = FREQ_BANDS{freq_i};
            t_conn_mats{cnt} = statMat;
            t_conn_sigs{cnt} = extract_sig;
            t_conn_meas{cnt} = conn_estimators{conn_i};
            cnt = cnt + 1;
        end
    end
end
%- create table
conn_mat_table = table(t_fPaths,t_fNames,t_conn_comps,t_conn_meas,t_conn_freqs,t_conn_sigs,t_conn_mats);
%% ===================================================================== %%
%## STEP 7) WRAP UP
%- REMOVE FIELDS && SAVE
for trial_i = 1:length(ALLEEG)
    %- clear STATS & PCONN for space saving.
    if isfield(ALLEEG(trial_i).CAT,'Stats')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'Stats');
    end
    if isfield(ALLEEG(trial_i).CAT,'PConn')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'PConn');
    end
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');

%% DONE
fprintf('DONE: %s\n',ALLEEG(1).subject);
%## SAVE
% ALLEEG
% conn_mat_table
% bootsrap_dist
% phasernd_dist
%## TIME
error_code = 1;
toc
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