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

function [EEG] = cnctanl_connMeas(EEG,PATHS,components,connMeasures,varargin)
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {0.5,...
    0.025,...
    [],...
    'nogui',...
    1,...
    [],...
    []};
p = inputParser;
validPath = (@(x) ischar(x) || isempty(x));
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'PATHS',@isstruct)
addRequired(p,'components',@isnumeric)
addRequired(p,'connMeasures',@iscell)
%## OPTIONAL
addOptional(p,'savePath',Defaults{6},validPath)
%## PARAMETER
addParameter(p,'WindowLengthSec',Defaults{1}, @isnumeric) % number of current set in previous load
addParameter(p,'WindowStepSizeSec',Defaults{2}, @isnumeric) % MIM folder overwrite in special cases
addParameter(p,'NewSamplingRate',Defaults{3}, @isnumeric)
addParameter(p,'GUI_MODE',Defaults{4}, @ischar)
addParameter(p,'VERBOSITY_LEVEL',Defaults{5}, @isnumeric)
parse(p, EEG, PATHS, components, connMeasures, varargin{:});

%## SET DEFAULTS
% Defaults OR Override Values
WindowLengthSec     = p.Results.WindowLengthSec;               % sliding window length in seconds
WindowStepSizeSec   = p.Results.WindowStepSizeSec;             % sliding window step size in seconds
NewSamplingRate     = p.Results.NewSamplingRate;               % new sampling rate (if downsampling)
GUI_MODE            = p.Results.GUI_MODE;                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = p.Results.VERBOSITY_LEVEL;               % Verbosity Level (0=no/minimal output, 2=graphical output)
savePath            = p.Results.savePath;
try
    if strncmp(computer,'PC',2)
        DO_UNIX = false;
    else
        DO_UNIX = true;
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end

if DO_UNIX
    savePath = convertPath2UNIX(savePath,'dferris');
else
    savePath = convertPath2Drive(savePath,'M');
end
%% PREP SAVE DIRECTORIES
% [savePath,~] = get_SubPath(PATHS,savePath,'SIFT');

%% STEP 3: Pre-process the data
fprintf('===================================================\n');
disp('PRE-PROCESSING DATA');
fprintf('===================================================\n');

% select components from EEG
EEG = pop_subcomp(EEG, components, 0, 1);

% resample data
if ~isempty(NewSamplingRate)
    EEG = pop_resample( EEG, NewSamplingRate);
end

% convert list of components to cell array of strings
ComponentNames = [];
for j = 1:length(components)
    ComponentNames = [ComponentNames, {num2str(components(j))}];
end

% This is a SIFT toolbox function. 
% NOTE (01/14/2022): look into the varargins of the function to determine 
% what they do. (See. pre_prepData.m (SIFT toolbox) for more information 
% regarding varargins).
%% ==== PARAMS ==== %%
SEQLEN          = 0.5; % sec
STEPSIZE        = 0.025; % sec
%% ==== END: PARAMS ==== %%
cfg             = [];
cfg.verb = 1;
cfg.sigtype = [];
    cfg.sigtype.arg_direct = 0;
    cfg.sigtype.arg_selection = 'Components';
cfg.varnames = ComponentNames;
cfg.diff = [];
    cfg.diff.arg_direct = 0;
    cfg.diff.arg_selection = 0;
cfg.detrend = [];
    cfg.detrend.arg_direct = 0;
    cfg.detrend.verg = 1;
    cfg.detrend.method = {'linear'};
    cfg.detrend.piecewise = [];
        cfg.detrend.piecewise.arg_direct = 0;
        cfg.detrend.piecewise.seglength = SEQLEN; % sec
        cfg.detrend.piecewise.stepsize = STEPSIZE; % sec
        cfg.detrend.piecewise.arg_selection = 1;
    cfg.detrend.plot = 1;
    cfg.detrend.arg_selection = 1;
cfg.normalize = [];
    cfg.normalize.arg_direct = 0;
    cfg.normalize.verb = 1;
    cfg.normalize.method = {'time','ensemble'};
    cfg.normalize.arg_selection = 1;
cfg.resetConfigs  = 0;
cfg.aamp = [];
cfg.badsegments = [];
cfg.newtrials = [];
cfg.equalizetrials = 0;

%## CONSIDERATION:
% (NOTE FROM EXAMPLE SIFT SCRIPT) No piecewise detrending based on conversation with Andreas Widmann.
[EEG,cfg] = feval(@pre_prepData,'EEG',EEG,cfg);
EEG.CAT.configs.('pre_prepData') = cfg;
fprintf('\n\n');

%% STEP 4: Identify the optimal model order
fprintf('===================================================\n');
disp('MODEL ORDER IDENTIFICATION');
fprintf('===================================================\n');
% Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results

%## PARAMS
MORDER                  = [1 30];
INFC_SELECTOR           = {'sbc' 'aic' 'fpe' 'hq'};
APPRCH_ALG              = 'Vieira-Morf';
APPRCH_TAPER            = 'rectwin';
APPRCH_PRCT_WIN2SAMP    = 100;
PRCT_WIN2SAMP           = 80;
SEQLEN                  = 0.5; % sec
STEPSIZE                = 0.025; % sec
%END
% cfg = arg_tovals(arg_report('rich',@est_selMModelOrder,[{'EEG',EEG},varargin]),false);
cfg = [];
cfg.verb = VERBOSITY_LEVEL;
cfg.modelingApproach = [];
    cfg.modelingApproach.arg_direct = 0;
    cfg.modelingApproach.algorithm = [];
        cfg.modelingApproach.algorithm.arg_direct = 0;
        cfg.modelingApproach.algorithm.morder = 10;
        cfg.modelingApproach.algorithm.arg_selection = APPRCH_ALG;
    cfg.modelingApproach.morder = 10;
    cfg.modelingApproach.winStartIdx =  [];
    cfg.modelingApproach.winlen = SEQLEN;
    cfg.modelingApproach.winstep = STEPSIZE;
    cfg.modelingApproach.taperfcn = APPRCH_TAPER;
    cfg.modelingApproach.epochTimeLims = [];
    cfg.modelingApproach.prctWinToSample = APPRCH_PRCT_WIN2SAMP;
    cfg.modelingApproach.normalize = [];
    cfg.modelingApproach.detrend = [];
        cfg.modelingApproach.detrend.arg_direct = 0;
        cfg.modelingApproach.detrend.method = 'constant';
        cfg.modelingApproach.detrend.arg_selection = 1;
    cfg.modelingApproach.verb = VERBOSITY_LEVEL;
    cfg.modelingApproach.timer = 0;
    cfg.modelingApproach.setArgDirectMode = 1;
    cfg.modelingApproach.arg_selection = 'Segmentation VAR';
cfg.morderRange = MORDER;
cfg.downdate = 1;
cfg.runPll = [];
    cfg.runPll.arg_direct = 0;
    cfg.runPll.arg_selection = 0;
cfg.icselector = INFC_SELECTOR;
cfg.winStartIdx = [];
cfg.epochTimeLims = [];
cfg.prctWinToSample = PRCT_WIN2SAMP;
cfg.plot = [];
    cfg.plot.arg_direct = 0;
    cfg.plot.arg_selection = 0;
% pop_est_selModelOrder
[IC,cfg] = feval(@est_selModelOrder,'EEG',EEG,cfg);
EEG.CAT.IC = IC;
EEG.CAT.configs.('est_selModelOrder') = cfg;
% To plot the results, use this:
handles = vis_plotOrderCriteria(EEG.CAT.IC,'conditions', [],    ...
                                            'icselector', {'sbc','aic','fpe','hq'},  ...
                                            'minimizer', 'min', ...
                                            'prclim', 90);
ModelOrder = ceil(mean(EEG.CAT.IC.hq.popt));

% If you want to save this figure you can uncomment the following lines:
saveas(handles,[savePath filesep sprintf('%s_orderResults.fig',EEG.subject)]);
close(handles);

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

% Here we can check that our selected parameters make sense
fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n\n',EEG.condition);

est_dispMVARParamCheck(EEG,struct('morder',ModelOrder','winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',VERBOSITY_LEVEL))

% Once we have identified our optimal model order, we can fit our VAR model.

% Fit a model using the options specifed for model order selection (STEP 4)
[EEG] = pop_est_fitMVAR(EEG,GUI_MODE, ...
        EEG.CAT.configs.est_selModelOrder.modelingApproach, ...
        'ModelOrder',ModelOrder);

% Note that EEG.CAT.MODEL now contains the model structure with
% coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
% self-evident information

% Alternately, we can fit the VAR parameters using a Kalman filter (see
% doc est_fitMVARKalman for more info on arguments)
%
% EEG.CAT.MODEL = est_fitMVARKalman(EEG,0,'updatecoeff',0.0005,'updatemode',2,'morder',ModelOrder,'verb',2,'downsampleFactor',50);
fprintf('\n\n');
%% STEP 6: Validate the fitted model
fprintf('===================================================\n');
disp('MODEL VALIDATION');
fprintf('===================================================\n');

% Here we assess the quality of the fit of our model w.r.t. the data. This
% step can be slow.

%## PARAMS
cfg                     = [];
CHKW_ALPHA              = 0.05;
CHKW_CORR               = 'none';
CHKW_LAGS               = 50;
CHKW_CRIT               = {'Ljung-Box' 'ACF' 'Box-Pierce' 'Li-McLeod'};
CHKW_PRCT_WIN2SAMP      = 100;

CHKR_ALPHA              = 0.05;
CHKR_CORR               = 'none';
CHKR_LAGS               = 50;
CHKR_CRIT               = {};
CHKR_PRCT_WIN2SAMP      = 100;

CHKC_PRCT_WIN2SAMP      = 100;
CHKC_DONORM             = 0;
CHKC_VERB               = 0;

CHKS_PRCT_WIN2SAMP      = 100;
%END

% We can obtain statistics for residual whiteness, percent consistency, and
% model stability ...
[EEG] = pop_est_validateMVAR(EEG,GUI_MODE,...
                            'checkWhiteness', ...
                                {'alpha' CHKW_ALPHA ...
                                 'statcorrection' CHKW_CORR ...
                                 'numAcfLags' CHKW_LAGS         ...
                                 'whitenessCriteria' CHKW_CRIT ...
                                 'winStartIdx' [] ...
                                 'prctWinToSample' CHKW_PRCT_WIN2SAMP  ...
                                 'verb' 0}, ...
                             'checkResidualVariance',...
                                {'alpha' CHKR_ALPHA ...
                                 'statcorrection' CHKR_CORR ...
                                 'numAcfLags' CHKR_LAGS    ...
                                 'whitenessCriteria' CHKR_CRIT  ...
                                 'winStartIdx' []        ...
                                 'prctWinToSample' CHKR_PRCT_WIN2SAMP   ...
                                 'verb' 0}, ...
                             'checkConsistency',    ...
                                {'winStartIdx' []   ...
                                 'prctWinToSample' CHKC_PRCT_WIN2SAMP ...
                                 'Nr' []                ...
                                 'donorm' CHKC_DONORM         ...
                                 'nlags' []         ...
                                 'verb' CHKC_VERB}, ...
                             'checkStability',  ...
                                {'winStartIdx' []   ...
                                 'prctWinToSample' CHKS_PRCT_WIN2SAMP ...
                                 'verb' 0},     ...
                             'prctWinToSample',70,  ...
                             'winStartIdx',[],      ...
                             'verb',VERBOSITY_LEVEL,...
                             'plot',false);

% ... and then plot the results 
handles = vis_plotModelValidation({EEG.CAT.VALIDATION.whitestats}, ...
                                     {EEG.CAT.VALIDATION.PCstats},         ...
                                     {EEG.CAT.VALIDATION.stabilitystats});                                        
% If you want to save this figure you can uncomment the following lines:
saveas(handles,[savePath filesep sprintf('%s_validationResults.fig',EEG.subject)]);
close(handles);

% To automatically determine whether our model accurately fits the data you
% can write a few lines as follows (replace 'acf' with desired statistic):
%
if ~all(EEG.CAT.VALIDATION.whitestats.acf.w)
    msgbox('Residuals are not completely white!');
end

fprintf('\n\n');

%% STEP 7: Compute Connectivity
fprintf('===================================================\n');
disp('CONNECTIVITY ESTIMATION');
fprintf('===================================================\n');

% Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model. We can compute these for a range of
% frequencies (here 1-40 Hz). See 'doc est_mvarConnectivity' for a complete
% list of available connectivity and spectral estimators.

%## PARAMS
cfg             = [];
% CONN_METHOD      = {'mCoh' 'dDTF08' 'ffDTF' 'S'};
ABSVALSQ        = true;
SPECTRAL_DB     = true;
FREQS           = [1:100];
%END

EEG = pop_est_mvarConnectivity(EEG,GUI_MODE, ...
            'connmethods', connMeasures, ...
            'absvalsq',ABSVALSQ,           ...
            'spectraldecibels',SPECTRAL_DB,   ...
            'freqs',FREQS,        ...
            'verb',VERBOSITY_LEVEL);

disp('done')
disp('===================================')
%## TIME
toc
%% SUBFUNCTIONS

%% DOCUMENTATION PER FUNCTION 
% ----------------------------------------------------------------------- %
%## FUNCTION: PRE_PREPDATA
% cfg             = [];
% cfg.verb = 1;
% cfg.sigtype = [];
%     cfg.sigtype.arg_direct = 0;
%     cfg.sigtype.arg_selection = 'Components';
% cfg.varnames = {};
% cfg.diff = [];
%     cfg.diff.arg_direct = 0;
%     cfg.diff.arg_selection = 0;
% cfg.detrend = [];
%     cfg.detrend.arg_direct = 0;
%     cfg.detrend.verg = 1;
%     cfg.detrend.method = {'linear'};
%     cfg.detrend.piecewise = [];
%         cfg.detrend.piecewise.arg_direct = 0;
%         cfg.detrend.piecewise.seglength = SEQLEN; % sec
%         cfg.detrend.piecewise.stepsize = STEPSIZE; % sec
%         cfg.detrend.piecewise.arg_selection = 1;
%     cfg.detrend.plot = 1;
%     cfg.detrend.arg_selection = 1;
% cfg.normalize = [];
%     cfg.normalize.arg_direct = 0;
%     cfg.normalize.verb = 1;
%     cfg.normalize.method = {'time','ensemble'};
%     cfg.normalize.arg_selection = 1;
% cfg.resetConfigs  = 0;
% cfg.aamp = [];
% cfg.badsegments = [];
% cfg.newtrials = [];
% cfg.equalizetrials = 0;
% ----------------------------------------------------------------------- %
%## FUNCTION: EST_SELMODELORDER
% cfg = [];
% cfg.verb = VERBOSITY_LEVEL;
% cfg.modelingApproach = [];
%     cfg.modelingApproach.arg_direct = 0;
%     cfg.modelingApproach.algorithm = [];
%         cfg.modelingApproach.algorithm.arg_direct = 0;
%         cfg.modelingApproach.algorithm.morder = 10;
%         cfg.modelingApproach.algorithm.arg_selection = APPRCH_ALG;
%     cfg.modelingApproach.morder = 10;
%     cfg.modelingApproach.winStartIdx =  [];
%     cfg.modelingApproach.winlen = SEQLEN;
%     cfg.modelingApproach.winstep = STEPSIZE;
%     cfg.modelingApproach.taperfcn = APPRCH_TAPER;
%     cfg.modelingApproach.epochTimeLims = [];
%     cfg.modelingApproach.prctWinToSample = APPRCH_PRCT_WIN2SAMP;
%     cfg.modelingApproach.normalize = [];
%     cfg.modelingApproach.detrend = [];
%         cfg.modelingApproach.detrend.arg_direct = 0;
%         cfg.modelingApproach.detrend.method = 'constant';
%         cfg.modelingApproach.detrend.arg_selection = 1;
%     cfg.modelingApproach.verb = VERBOSITY_LEVEL;
%     cfg.modelingApproach.timer = 0;
%     cfg.modelingApproach.setArgDirectMode = 1;
%     cfg.modelingApproach.arg_selection = 'Segmentation VAR';
% cfg.morderRange = MORDER;
% cfg.downdate = 1;
% cfg.runPll = [];
%     cfg.runPll.arg_direct = 0;
%     cfg.runPll.arg_selection = 0;
% cfg.icselector = INFC_SELECTOR;
% cfg.winStartIdx = [];
% cfg.epochTimeLims = [];
% cfg.prctWinToSample = PRCT_WIN2SAMP;
% cfg.plot = [];
%     cfg.plot.arg_direct = 0;
%     cfg.plot.arg_selection = 0;
% ----------------------------------------------------------------------- %
%## FUNCTION: 
% ----------------------------------------------------------------------- %
%## FUNCTION: 
% ----------------------------------------------------------------------- %
%## FUNCTION: 