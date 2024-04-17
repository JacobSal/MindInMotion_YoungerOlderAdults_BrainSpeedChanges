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
% DEF_PREPDATA = struct('VerbosityLevel',VERBOSITY_LEVEL,...
%              'SignalType',{{'Components'}},...
%              'VariableNames',[],...
%              'Detrend',{{'verb',VERBOSITY_LEVEL,'method',{'linear'}}},...
%              'NormalizeData',{{'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}}},...
%              'resetConfigs',true,...
%              'badsegments',[],...
%              'newtrials',[],...
%              'equalizetrials',false);
%## ESTIAMTE MODEL ORDER PARAMS
% DEF_ESTSELMOD = struct('verb',VERBOSITY_LEVEL,...
%     'modelingApproach',struct('arg_direct',1,...
%                             'algorithm',struct('arg_direct',1,...
%                                             'morder',[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1)],...
%                                             'arg_selection',{{'Vieira-Morf'}}),...
%                             'morder',ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1),...
%                             'winStartIdx',[],...
%                             'winlen',WINDOW_LENGTH,...
%                             'winstep',WINDOW_STEP_SIZE,...
%                             'taperfcn','rectwin',...
%                             'epochTimeLims',[],...
%                             'prctWinToSample',100,...
%                             'normalize',[],...
%                             'detrend',struct('arg_direct',1,...
%                                             'method','constant',...
%                                             'arg_selection',1),...
%                             'verb',VERBOSITY_LEVEL,...
%                             'timer',0,...
%                             'setArgDirectMode',1,...
%                             'arg_selection','Segmentation VAR'),...
%     'morderRange',[1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1)],...
%     'downdate',1,...
%     'runPll',struct('arg_direct',1,'arg_selection',0),...
%     'icselector',{{'aic','hq'}},...
%     'winStartIdx',[],...
%     'epochTimeLims',[],...
%     'prctWinToSample',80,...
%     'plot',struct('arg_direct',1,...
%                 'arg_selection',0));
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
% addParameter(p,'ESTVALMVAR_CW',DEF_ESTVALMVAR_CW,@(x) validate_struct(x,DEF_ESTVALMVAR_CW));
% addParameter(p,'ESTVALMVAR_CRV',DEF_ESTVALMVAR_CRV,@(x) validate_struct(x,DEF_ESTVALMVAR_CRV));
% addParameter(p,'ESTVALMVAR_CC',DEF_ESTVALMVAR_CC,@(x) validate_struct(x,DEF_ESTVALMVAR_CC));
% addParameter(p,'ESTVALMVAR_CS',DEF_ESTVALMVAR_CS,@(x) validate_struct(x,DEF_ESTVALMVAR_CS));
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


