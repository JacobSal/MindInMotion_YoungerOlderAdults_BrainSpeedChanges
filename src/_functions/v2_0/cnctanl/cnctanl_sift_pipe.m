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
topt = tic();
% com.mathworks.mwswing.MJUtilities.initJIDE;
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
    'morderRange',[10,30],...
    'downdate',true,...
    'RunInParallel',{{'profile','local',...
        'numWorkers',2}},...
    'icselector',{{'aic','hq'}},...
    'winStartIdx',[],...
    'epochTimeLims',[],...
    'prctWinToSample',80,...
    'plot',[],...
    'verb',VERBOSITY_LEVEL);
%- (04/21/2024) JS, changing runPll to on and running each subject at a
%a time to test rpocessing time increases
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
                                 'verb',0}}, ...
                         'checkResidualVariance',{{'alpha',0.05,...
                                 'statcorrection','none',...
                                 'numAcfLags',50,...
                                 'whitenessCriteria',{},...
                                 'winStartIdx',[],...
                                 'prctWinToSample',100,...
                                 'verb',0}}, ...
                         'checkConsistency',{{'winStartIdx',[],...
                                 'prctWinToSample',100,...
                                 'Nr',[],...
                                 'donorm',0,...
                                 'nlags',[],...
                                 'verb',0}}, ...
                         'checkStability',{{'winStartIdx',[],...
                                 'prctWinToSample',100,...
                                 'verb',0}},     ...
                         'prctWinToSample',70,...
                         'winStartIdx',[],      ...
                         'verb',0,...
                         'plot',false);
    %- (04/19/2024) JS, setting verbosity to 0 because when ran on
    %hipergator with 2023b the txt output creates very large (~10GB) log
    %files that are difficult to open. Perhaps only est_checkMVARConsistency is to blame... Specifically this:
    %       Warning: The specified AR model is unstable at timeindex 1.
    %       > In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('tvarsim', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/sim/tvarsim.m', 112)" style="font-weight:bold">tvarsim</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/sim/tvarsim.m',112,0)">line 112</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('est_consistency', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_consistency.m', 78)" style="font-weight:bold">est_consistency</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_consistency.m',78,0)">line 78</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('est_checkMVARConsistency', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_checkMVARConsistency.m', 135)" style="font-weight:bold">est_checkMVARConsistency</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_checkMVARConsistency.m',135,0)">line 135</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('est_validateMVAR', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_validateMVAR.m', 210)" style="font-weight:bold">est_validateMVAR</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/est/est_validateMVAR.m',210,0)">line 210</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('pop_est_validateMVAR', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/pop/pop_est_validateMVAR.m', 147)" style="font-weight:bold">pop_est_validateMVAR</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/SIFT/pop/pop_est_validateMVAR.m',147,0)">line 147</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('cnctanl_sift_pipe', '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/cnctanl/cnctanl_sift_pipe.m', 356)" style="font-weight:bold">cnctanl_sift_pipe</a> (<a href="matlab: opentoline('/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/cnctanl/cnctanl_sift_pipe.m',356,0)">line 356</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('parallel_function>make_general_channel/channel_general', '/apps/matlab/r2023b/toolbox/matlab/lang/parallel_function.m', 837)" style="font-weight:bold">parallel_function>make_general_channel/channel_general</a> (<a href="matlab: opentoline('/apps/matlab/r2023b/toolbox/matlab/lang/parallel_function.m',837,0)">line 837</a>)
    %       In <a href="matlab:matlab.internal.language.introspective.errorDocCallback('parallel.internal.parfor.cppRemoteParallelFunction', '/apps/matlab/r2023b/toolbox/parallel/cluster/+parallel/+internal/+parfor/cppRemoteParallelFunction.m', 53)" style="font-weight:bold">parallel.internal.parfor.cppRemoteParallelFunction</a> (<a href="matlab: opentoline('/apps/matlab/r2023b/toolbox/parallel/cluster/+parallel/+internal/+parfor/cppRemoteParallelFunction.m',53,0)">line 53</a>)
    %
%## CONNECTIVITY ESTIMATION PARAMS
DEF_ESTFITMVAR = struct('connmethods',{{'dDTF08','S'}}, ...
            'absvalsq',ABSVALSQ,           ...
            'spectraldecibels',SPECTRAL_DB,   ...
            'freqs',FREQS,        ...
            'verb',1);
    %- (04/20/2024) JS, setting verb to 0 to suppress memory dlg boxes and
    %just proc error. 
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
%- MAIN PROCESSING PARAMS
addParameter(p,'PREPDATA',DEF_PREPDATA,@(x) validate_struct(x,DEF_PREPDATA));
addParameter(p,'ESTSELMOD',DEF_ESTSELMOD,@(x) validate_struct(x,DEF_ESTSELMOD));
addParameter(p,'PLOTORDERCRIT',DEF_PLOTORDERCRIT,@(x) validate_struct(x,DEF_PLOTORDERCRIT));
addParameter(p,'ESTVALMVAR',DEF_ESTVALMVAR,@(x) validate_struct(x,DEF_ESTVALMVAR));
addParameter(p,'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,@(x) validate_struct(x,DEF_ESTDISPMVAR_CHK));
addParameter(p,'ESTFITMVAR',DEF_ESTFITMVAR,@(x) validate_struct(x,DEF_ESTFITMVAR));
%- SURROGATE STATS PARAMS
addParameter(p,'STAT_BS_CFG',DEF_STAT_BS_CFG,@(x) validate_struct(x,DEF_STAT_BS_CFG));
addParameter(p,'STAT_PR_CFG',DEF_STAT_PR_CFG,@(x) validate_struct(x,DEF_STAT_PR_CFG));
%- SWITCHES
addParameter(p,'DO_PHASE_RND',DO_PHASE_RND,dpr_validFcn);
addParameter(p,'DO_BOOTSTRAP',DO_BOOTSTRAP,dbs_validFcn);
parse(p,ALLEEG,conn_components,save_dir,varargin{:});
%## SET DEFAULTS
%- Optional
PREPDATA       = p.Results.PREPDATA;
ESTSELMOD       = p.Results.ESTSELMOD;
PLOTORDERCRIT = p.Results.PLOTORDERCRIT;
ESTVALMVAR = p.Results.ESTVALMVAR;
ESTDISPMVAR_CHK       = p.Results.ESTDISPMVAR_CHK;
ESTFITMVAR      = p.Results.ESTFITMVAR;
%-
STAT_BS_CFG       = p.Results.STAT_BS_CFG;
STAT_PR_CFG       = p.Results.STAT_PR_CFG;
%-
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
tt = tic();
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
toc(tt);
%% (MAIN CONNECTIVITY PIPELINE) ======================================== %%
%## STEP 3: Pre-process the data
fprintf('===================================================\nPRE-PROCESSISNG DATA\n===================================================\n');
tt = tic();
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
fprintf('%s) Pre-Processing Data Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%% STEP 4: Identify the optimal model order
fprintf('===================================================\nMODEL ORDER IDENTIFICATION\n===================================================\n');
tt = tic();
%- NOTE: Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results
%## OPTION 1
% cfg = struct2args(ESTSELMOD);
% ALLEEG(cond_i) = pop_est_selModelOrder(ALLEEG(cond_i),GUI_MODE,cfg{:});
%## OPTION 2
cfg = struct2args(ESTSELMOD);
for cond_i=1:length(ALLEEG)
    %* calculate the information criteria
    [ALLEEG(cond_i).CAT.IC,cfg] = est_selModelOrder(ALLEEG(cond_i),cfg{:});
    % [ALLEEG(cond_i).CAT.IC,cfg] = est_selModelOrder(ALLEEG(cond_i),ESTSELMOD);
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
ALLEEG = feval(modFuncName,ALLEEG,GUI_MODE,ALLEEG(1).CAT.IC.modelFitting.modelingArguments);
fprintf('%s) Model Order Identification Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%## (PLOT) ============================================================= %%
tt = tic();
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
fprintf('%s) Plotting Order Criteria Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%% STEP 5: Fit the VAR model
fprintf('===================================================\nMODEL FITTING\n===================================================\n');
fprintf('\n');
tt = tic();
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
fprintf('%s) Model Fitting Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%% STEP 6: Validate the fitted model
fprintf('===================================================\nMODEL VALIDATION\n===================================================\n');
tt = tic();
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
fprintf('%s) Model Validation Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%% STEP 7: Compute Connectivity
%- NOTE: Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model.
fprintf('===================================================\nCONNECTIVITY ESTIMATION\n===================================================\n');
tt = tic();
% cfg = struct2args(ESTFITMVAR);
% tmp = ALLEEG;
% cfg = struct2args(ESTFITMVAR);
% Conn = pop_est_mvarConnectivity(tmp,GUI_MODE, ...
%             'connmethods',{ESTFITMVAR.connmethods},...
%             'absvalsq',ESTFITMVAR.absvalsq,...
%             'spectraldecibels',ESTFITMVAR.spectraldecibels,...
%             'freqs',ESTFITMVAR.freqs,...
%             'verb',1);
cfg = struct2args(ESTFITMVAR);
for cond_i = 1:length(ALLEEG)
    [Conn,~] = est_mvarConnectivity('ALLEEG',ALLEEG(cond_i),'MODEL',ALLEEG(cond_i).CAT.MODEL,...
                cfg{:});
    if ~isempty(Conn)
        ALLEEG(cond_i).CAT.Conn = Conn; 
    end
    % clear any existing visualization GUI config files
    visFields = fieldnames(ALLEEG(cond_i).CAT.configs);
    visFields = visFields(~cellfun(@isempty,strfind(visFields,'vis_')));
    for k=1:length(visFields)
        ALLEEG(cond_i).CAT.configs.(visFields{k}) = struct([]);
    end
    
    if ~isempty(cfg)
        % store the configuration structure
        ESTFITMVAR.arg_direct = 0;
        ALLEEG(cond_i).CAT.configs.('est_mvarConnectivity') = ESTFITMVAR;
    end
end
fprintf('%s) Connectivity Estimation  Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
%% ===================================================================== %%
%## STEP 5.a) (BOOTSTRAPPING) GROUP STATISTICS 
% (09/22/2022), JS, Might want to try and speed up bootstrap by
% adapting stat_surrogateGen.m to use parfor for bootstrapping... If
% possible? doesn't seem built well in the first place, so maybe?
% (10/27/2022), JS, Revist above note again!
% (12/7/2022), JS, need to update this boostrapping to include ALLEEG
if DO_BOOTSTRAP
    fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
    tt = tic();
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
    fprintf('%s) Bootstrap Calculation Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
end
%% ===================================================================== %%
%## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION    
% see. stat_surrogateGen
% see. stat_surrogateStats
if DO_PHASE_RND
    fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    tt = tic();
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
    fprintf('%s) Phase Randomization Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
end

%% ===================================================================== %%
%## STEP 6) CONNECTIVITY MATRICES & VISUALS
%## TABLE VARS
t_fPaths = cell(length(ALLEEG)*length(ESTFITMVAR.connmethods),1);
t_fNames = cell(length(ALLEEG)*length(ESTFITMVAR.connmethods),1);
t_conn_comps = cell(length(ALLEEG)*length(ESTFITMVAR.connmethods),1);
t_conn_meas = cell(length(ALLEEG)*length(ESTFITMVAR.connmethods),1);
t_cond_char = cell(length(ALLEEG)*length(ESTFITMVAR.connmethods),1);
disp(ALLEEG);
%## Generate Connectivity Matrix
cnt = 1;
for conn_i = 1:length(ESTFITMVAR.connmethods)
    for cond_i = 1:length(ALLEEG)
        disp(ALLEEG(cond_i).CAT.Conn);
        t_conn_comps{cnt} = conn_components;
        t_fPaths{cnt} = ALLEEG(cond_i).filepath;
        t_fNames{cnt} = ALLEEG(cond_i).filename;
        t_conn_meas{cnt} = ESTFITMVAR.connmethods{conn_i};
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
%% DONE
fprintf(1,'DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s: %0.1f\n\n',ALLEEG(1).subject,toc(topt));          

end
