function [ALLEEG] = cnctanl_conn_meas(ALLEEG,comps,conn_measures,varargin)
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
%## DEFINE DEFAULTS

FREQS = (1:60);
WINDOW_LENGTH = 0.5; % time (s)
WINDOW_STEP_SIZE = 0.025; % time (s)
SAMPLING_RATE  = []; % empty (i.e., []) or integer
VERBOSITY_LEVEL = 1;
GUI_MODE = 'nogui';
%- figure savepath
sd_validfunc = (@(x) exist(x,'dir'));
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct)
addRequired(p,'comps',@isnumeric)
addRequired(p,'conn_measures',@iscell)
addRequired(p,'save_dir',sd_validfunc)
%## OPTIONAL

%## PARAMETER
addParameter(p,'FREQS',FREQS, @isnumeric) 
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH, @isnumeric) % number of current set in previous load
addParameter(p,'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE, @isnumeric) % MIM folder overwrite in special cases
addParameter(p,'SAMPLING_RATE',SAMPLING_RATE, @isnumeric)
addParameter(p,'GUI_MODE',GUI_MODE, @ischar)
addParameter(p,'VERBOSITY_LEVEL',VERBOSITY_LEVEL, @isnumeric)

parse(p, ALLEEG, comps, conn_measures, varargin{:});

%## SET DEFAULTS
%- Optional
save_dir     = p.Results.save_dir;
%- parameters
WINDOW_LENGTH    = p.Results.WINDOW_LENGTH;               % sliding window length in seconds
WINDOW_STEP_SIZE = p.Results.WINDOW_STEP_SIZE;             % sliding window step size in seconds
SAMPLING_RATE    = p.Results.SAMPLING_RATE;               % new sampling rate (if downsampling)
GUI_MODE         = p.Results.GUI_MODE;                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL  = p.Results.VERBOSITY_LEVEL;              % Verbosity Level (0=no/minimal output, 2=graphical output)
FREQS            = p.Results.FREQS;
%% ===================================================================== %%
%## STEP 3: Pre-process the data
fprintf('===================================================\n');
disp('PRE-PROCESSING DATA');
fprintf('===================================================\n');

% resample data
if ~isempty(SAMPLING_RATE)
    ALLEEG = pop_resample(ALLEEG, SAMPLING_RATE);
end

% convert list of components to cell array of strings
comp_names = [];
for j = 1:length(comps) 
    comp_names = [comp_names, {num2str(comps(j))}];
end

%%
% This is a SIFT toolbox function. 
% NOTE (01/14/2022): look into the varargins of the function to determine 
% what they do. (See. pre_prepData.m (SIFT toolbox) for more information 
% regarding varargins).
%##
% g             = [];
% %- set verbosity and assing components
% g.verb = 1;
% g.sigtype = [];
%     g.sigtype.arg_direct = 1;
%     g.sigtype.arg_selection = 'Components';
% g.varnames = ComponentNames;
% g.diff = [];
%     g.diff.arg_direct = 0;
%     g.diff.arg_selection = 0;
% %- detrend using linear non-piecewise
% % least-squares fit method
% g.detrend = [];
%     g.detrend.arg_direct = 1;
%     g.detrend.verg = 1;
%     g.detrend.method = {'linear'};
%     g.detrend.piecewise = [];
%         g.detrend.piecewise.arg_direct = 0;
%         g.detrend.piecewise.seglength = WINDOW_LENGTH; % sec
%         g.detrend.piecewise.stepsize = WINDOW_STEP_SIZE; % sec
%         g.detrend.piecewise.arg_selection = 1;
%     g.detrend.plot = 1;
%     g.detrend.arg_selection = 1;
% %- normalize across time and ensemble (trials)
% g.normalize = [];
%     g.normalize.arg_direct = 1;
%     g.normalize.verb = 1;
%     g.normalize.method = {'time','ensemble'};
%     g.normalize.arg_selection = 1;
% %- other bonus options
% g.resetConfigs  = 0;
% g.aamp = [];
% g.badsegments = [];
% g.newtrials = [];
% g.equalizetrials = 0;
% 
% %## CONSIDERATION:
% % (NOTE FROM EXAMPLE SIFT SCRIPT) No piecewise detrending based on conversation with Andreas Widmann.
% [EEG,g] = feval(@pre_prepData,'EEG',EEG,g);
% EEG.CAT.configs.('pre_prepData') = g;
% fprintf('\n\n');

%- ALTERNATIVE
[ALLEEG] = pop_pre_prepData(ALLEEG,'nogui',...
             'VerbosityLevel',VERBOSITY_LEVEL,...
             'SignalType',{'Components'},...
             'VariableNames',comp_names,...
             'Detrend',{'verb',VERBOSITY_LEVEL,'method',{'linear'}},...
             'NormalizeData',{'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}},...
             'resetConfigs',true,...
             'badsegments',[],...
             'newtrials',[],...
             'equalizetrials',false);

%% STEP 4: Identify the optimal model order
fprintf('===================================================\n');
disp('MODEL ORDER IDENTIFICATION');
fprintf('===================================================\n');
% Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results

%## PARAMS
MORDER_RANGE            = [1, ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1)];
MORDER                  = ceil((ALLEEG(1).srate*WINDOW_LENGTH)/2-1);
INFC_SELECTOR           = {'aic' 'hq'};
APPRCH_ALG              = {'Vieira-Morf'};
APPRCH_TAPER            = 'rectwin';
APPRCH_PRCT_WIN2SAMP    = 100;
PRCT_WIN2SAMP           = 80;
%END
% cfg = arg_tovals(arg_report('rich',@est_selMModelOrder,[{'EEG',EEG},varargin]),false);
cfg = [];
cfg.verb = VERBOSITY_LEVEL;
cfg.modelingApproach = [];
    cfg.modelingApproach.arg_direct = 1;
    cfg.modelingApproach.algorithm = [];
        cfg.modelingApproach.algorithm.arg_direct = 1;
        cfg.modelingApproach.algorithm.morder = MORDER;
        cfg.modelingApproach.algorithm.arg_selection = APPRCH_ALG;
    cfg.modelingApproach.morder = MORDER;
    cfg.modelingApproach.winStartIdx =  [];
    cfg.modelingApproach.winlen = WINDOW_LENGTH;
    cfg.modelingApproach.winstep = WINDOW_STEP_SIZE;
    cfg.modelingApproach.taperfcn = APPRCH_TAPER;
    cfg.modelingApproach.epochTimeLims = [];
    cfg.modelingApproach.prctWinToSample = APPRCH_PRCT_WIN2SAMP;
    cfg.modelingApproach.normalize = [];
    cfg.modelingApproach.detrend = [];
        cfg.modelingApproach.detrend.arg_direct = 1;
        cfg.modelingApproach.detrend.method = 'constant';
        cfg.modelingApproach.detrend.arg_selection = 1;
    cfg.modelingApproach.verb = VERBOSITY_LEVEL;
    cfg.modelingApproach.timer = 0;
    cfg.modelingApproach.setArgDirectMode = 1;
    cfg.modelingApproach.arg_selection = 'Segmentation VAR';
cfg.morderRange = MORDER_RANGE;
cfg.downdate = 1;
cfg.runPll = [];
    cfg.runPll.arg_direct = 1;
    cfg.runPll.arg_selection = 0;
cfg.icselector = INFC_SELECTOR;
cfg.winStartIdx = [];
cfg.epochTimeLims = [];
cfg.prctWinToSample = PRCT_WIN2SAMP;
cfg.plot = [];
    cfg.plot.arg_direct = 1;
    cfg.plot.arg_selection = 0;
% pop_est_selModelOrder
% for cond_i=1:length(ALLEEG)
%     [IC,cfg] = feval(@est_selModelOrder,'EEG',ALLEEG(cond_i),cfg);
%     ALLEEG(cond_i).CAT.IC = IC;
%     ALLEEG(cond_i).CAT.configs.('est_selModelOrder') = cfg;
% end
for cond_i=1:length(ALLEEG)
    % calculate the information criteria
    ALLEEG(cond_i).CAT.IC = est_selModelOrder('EEG',ALLEEG(cond_i),cfg);
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cond_i).CAT.configs.('est_selModelOrder') = cfg;
    end
end
% Open the model-fitting GUI for model fitting. 
% Once model is fit results will be return in EEG structure

%- ALTERNATIVE
%** NOTE: this may be buggy on HiPerGator/Hypercomputing systems due to the
% gui initiation. buggy because of feval statement: "
% modFuncName = ['pop_' ALLEEG(1).CAT.IC.modelFitting.modelingFuncName];
% ALLEEG = feval(modFuncName, ALLEEG,0,ALLEEG(1).CAT.IC.modelFitting.modelingArguments);

% ALLEEG = pop_est_selModelOrder(ALLEEG,'nogui', ...
%         'modelingApproach',         ...
%             {'Segmentation VAR',     ...
%                 'algorithm', APPRCH_ALG ... % Was 'Vieira-Morf'. 'Group Lasso (DAL/SCSA)' is recommended but slow. 'ARfit' is the simplest. (04/18/2021 Makoto)
%                 'winStartIdx', []    ...
%                 'winlen',  WINDOW_LENGTH    ...
%                 'winstep', WINDOW_STEP_SIZE  ...
%                 'taperfcn', APPRCH_TAPER    ... % For long sliding window like 20 s, rectangle window is suitable. (04/18/2021 Makoto)
%                 'epochTimeLims', []      ...
%                 'prctWinToSample', APPRCH_PRCT_WIN2SAMP   ...
%                 'normalize', []          ...
%                 'detrend', {'method', 'linear'} ...
%                 'verb', VERBOSITY_LEVEL},      ...
%         'morderRange',MORDER ,  ... % Typically below 10 using the elbow method. (04/18/2021 Makoto)
%         'downdate',true,        ...
%         'runPll',[],            ...
%         'icselector',INFC_SELECTOR,  ...
%         'winStartIdx',[],       ...
%         'epochTimeLims',[],     ...
%         'prctWinToSample',PRCT_WIN2SAMP,   ...
%         'plot', [], ...
%         'verb',VERBOSITY_LEVEL);

% To plot the results, use this:
tmp_morder = zeros(1,length(ALLEEG));
for cond_i = 1:length(ALLEEG)
    fprintf('%s) Plotting Validations',ALLEEG(cond_i).subject);
    handles = vis_plotOrderCriteria(ALLEEG(cond_i).CAT.IC,'conditions', [],    ...
                                            'icselector', INFC_SELECTOR,  ...
                                            'minimizer', 'min', ...
                                            'prclim', 90);
    tmp_morder(cond_i) = ceil(mean(ALLEEG(cond_i).CAT.IC.hq.popt));
    saveas(handles,[save_dir filesep sprintf('%s_%i_orderResults.fig',ALLEEG(cond_i).subject,cond_i)]);
    close(handles);
end
ModelOrder = ceil(mean(tmp_morder));

% If you want to save this figure you can uncomment the following lines:


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
for cond_i = 1:length(ALLEEG)
    fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n\n',ALLEEG(cond_i).condition);
    est_dispMVARParamCheck(ALLEEG(cond_i),struct('morder',ModelOrder','winlen',WINDOW_LENGTH,'winstep',WINDOW_STEP_SIZE,'verb',VERBOSITY_LEVEL))
end
% Once we have identified our optimal model order, we can fit our VAR model.

% Fit a model using the options specifed for model order selection (STEP 4)
[ALLEEG] = pop_est_fitMVAR(ALLEEG,GUI_MODE,...
        ALLEEG(1).CAT.configs.est_selModelOrder.modelingApproach,...
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
%- Whiteness Check
CHKW_ALPHA              = 0.05;
CHKW_CORR               = 'none';
CHKW_LAGS               = 50;
CHKW_CRIT               = {'Ljung-Box' 'ACF' 'Box-Pierce' 'Li-McLeod'};
CHKW_PRCT_WIN2SAMP      = 100;
%- Residuals Check
CHKR_ALPHA              = 0.05;
CHKR_CORR               = 'none';
CHKR_LAGS               = 50;
CHKR_CRIT               = {};
CHKR_PRCT_WIN2SAMP      = 100;
%- Consistency Check
CHKC_PRCT_WIN2SAMP      = 100;
CHKC_DONORM             = 0;
CHKC_VERB               = 0;
%- Stability Check
CHKS_PRCT_WIN2SAMP      = 100;
%END

% We can obtain statistics for residual whiteness, percent consistency, and
% model stability ...
[ALLEEG] = pop_est_validateMVAR(ALLEEG,GUI_MODE,...
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
fprintf('===================================================\n');
disp('CONNECTIVITY ESTIMATION');
fprintf('===================================================\n');

% Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model. We can compute these for a range of
% frequencies (here 1-40 Hz). See 'doc est_mvarConnectivity' for a complete
% list of available connectivity and spectral estimators.

%## PARAMS
ABSVALSQ        = false;
SPECTRAL_DB     = false;
%ENDn 

ALLEEG = pop_est_mvarConnectivity(ALLEEG,GUI_MODE, ...
            'connmethods', conn_measures, ...
            'absvalsq',ABSVALSQ,           ...
            'spectraldecibels',SPECTRAL_DB,   ...
            'freqs',FREQS,        ...
            'verb',VERBOSITY_LEVEL);

disp('done')
disp('===================================')

%## TIME
toc
end
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