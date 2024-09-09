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

function [EEG] = bsmart_cnctanl_connMeas(EEG,components,connMeasures,varargin)
%## TIME
tic
%## DEFINE DEFAULTS
%- optionals
fig_savepath = [];
%- parameters
SLIDE_WIN_LEN          = 0.5; % sec
SLIDE_STEP_SIZE        = 0.025; % sec
WINDOW_LENGTH = SLIDE_WIN_LEN*EEG.srate; %ceil(EEG.xmax-EEG.xmin)*EEG.srate;
MAX_MORDER = 30;
SAMPLING_RATE = [];
%- Detrending & Normalization
DETREND_E = false;
DETREND_T = true;
NORM_E = true;
NORM_T = true;
%## PARSE
p = inputParser;
validPath = (@(x) ischar(x) || isempty(x));
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'components',@isnumeric)
addRequired(p,'connMeasures',@iscell)
%## OPTIONAL
addOptional(p,'fig_savepath',fig_savepath,validPath)
%## PARAMETER
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@isnumeric) % number of current set in previous load
addParameter(p,'MAX_MORDER',MAX_MORDER,@isnumeric) % number of current set in previous load
addParameter(p,'SLIDE_WIN_LEN',SLIDE_WIN_LEN,@isnumeric) % number of current set in previous load
addParameter(p,'SLIDE_STEP_SIZE',SLIDE_STEP_SIZE,@isnumeric) % number of current set in previous load
addParameter(p,'SAMPLING_RATE',SAMPLING_RATE, @isnumeric)

parse(p, EEG, components, connMeasures, varargin{:});

%## SET DEFAULTS
%- Optionals
%- Parameters
WINDOW_LENGTH       = p.Results.WINDOW_LENGTH;              
MAX_MORDER          = p.Results.MAX_MORDER;             
SLIDE_WIN_LEN       = p.Results.SLIDE_WIN_LEN;               
SLIDE_STEP_SIZE     = p.Results.SLIDE_STEP_SIZE;             
fig_savepath        = p.Results.fig_savepath;
SAMPLING_RATE    = p.Results.SAMPLING_RATE;
%- modifications
% convert to OS
if strncmp(computer,'PC',2)
    fig_savepath = convertPath2UNIX(fig_savepath,'dferris');
else
    fig_savepath = convertPath2Drive(fig_savepath,'M');
end
% Check to see if integers

%% STEP 1: Pre-process the data
fprintf('===================================================\n');
fprintf('%10s PRE-PROCESSING DATA\n','');
fprintf('===================================================\n');

%- select components from EEG
EEG = pop_subcomp(EEG, components, 0, 1);

% resample data
if ~isempty(SAMPLING_RATE)
    EEG = pop_resample( EEG, SAMPLING_RATE);
end

%- convert list of components to cell array of strings
component_names = cell(1,length(components));
for j = 1:length(components)
    component_names{j} = num2str(components(j));
end

%% CONVERT EEGLAB DATA TO BSMART FORMAT
%- (Components) convert EEG to data for BSMART 
% data = time x channels/comps x trials
data = zeros(size(EEG.icaweights,1),size(EEG.data,2),size(EEG.data,3));
for dim = 1:size(EEG.data,3)
    data(:,:,dim) = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:,dim);
end
data = reshape(data,size(data,2),size(data,1),size(data,3));
writedat('dataset.bin',data)
%% STEP 2: Detrend and Normalize data
fprintf('===================================================\n');
fprintf('%10s MODEL DETRENDING & NORMALIZATION\n','');
fprintf('===================================================\n');
%## data Detrending
% NOTE: may want to consider adding linear least-squares to coinside with
% SIFT methods. For now using g.detrend.method = {'constant'};
%- randomly select a channel/component to observe
chan_i = randi(size(data,2));
%- subtract ensemble mean
if DETREND_E
    fprintf('Removing ensemble mean\n');
    figure;
    title('Remove Ensemble Mean');
    hold on;
    plot(data(:,chan_i,1),'DisplayName','Before')
    [data] = pre_sube(data);
    plot(data(:,chan_i,1),'DisplayName','After')
    legend()
    hold off;
end
%- subtract temporal mean
if DETREND_T
    fprintf('Removing temporal mean\n');
    figure;
    title('Remove Temporal Mean');
    hold on;
    plot(data(:,chan_i,1),'DisplayName','Before')
    [data] = pre_subt(data);
    plot(data(:,chan_i,1),'DisplayName','After')
    legend()
    hold off;
end
%## data Normalization
%- subtract temporal mean and divide by standard deviation
if NORM_T
    fprintf('Removing temporal mean and divide by standard deviation\n');
    figure;
    title('Normalize Across Time');
    hold on;
    plot(data(:,chan_i,1),'DisplayName','Before')
    [data] = pre_subt_divs(data);
    plot(data(:,chan_i,1),'DisplayName','After')
    legend()
    hold off;
end
%- subtract ensemble mean and divide by standard deviation
if NORM_E
    fprintf('Removing ensemble mean and divide by standard deviation\n');
    figure;
    title('Normalize Across Ensemble');
    hold on;
    plot(data(:,chan_i,1),'DisplayName','Before')
    [data] = pre_sube_divs(data);
    plot(data(:,chan_i,1),'DisplayName','After')
    legend()
    hold off;
end

fprintf('\n\n');
%## Hualou, Liang Granger-Geweke Causality (adaptation of
%bsmart-toolbox.org
[M,S] = compute_allnpCGCvar3(X,fs,fRes);
[M] = computeCGCvar3(X,fs,para)
%% BSMART Toolbox implementation
%##WARNING: FOR CODE BELOW, must recreate .exe's before code is usable
%{
%% STEP 3: Identify the optimal model order
fprintf('===================================================\n');
fprintf('%10s MODEL ORDER IDENTIFICATION \n','');
fprintf('===================================================\n');

%- Open GUI
%## data Model Order Selection using Information Criterion
% AIC = aic_test(data,SLIDE_WIN_LEN,MAX_MORDER);
% fprintf('Akaike information criterion: %i\n',AIC)
% mod_morder = min(AIC); % this may need fixing
% mod_morder = 10;
fprintf('\n\n');
%% STEP 4: Fit the VAR model
fprintf('===================================================\n');
fprintf('%10s MODEL FITTING\n','');
fprintf('===================================================\n');
% mod_start_samp  = 1;
% mod_end_samp    = size(data,1)-1;
%- moving window multivariate autoregressive model fitting
% [path_mod_AR, path_mod_COV] = mov_mul_model(data,mod_morder,mod_start_samp,mod_end_samp,SLIDE_WIN_LEN);
fprintf('\n\n');
%% STEP 5: Validate the fitted model
fprintf('===================================================\n');
fprintf('%10s MODEL VALIDATION\n','');
fprintf('===================================================\n');

%## data Whiteness test (need model order)
% mod_resid = whiteness_test(data,SLIDE_WIN_LEN,mod_morder);

%## model consistency test (need fitted model coefficients/covariances)
% mod_ratio = consistencytest(data,mod_arcoeff,mod_arnoise);

%## stability test
% mod_burnin = []; % not sure what this is
% mod_le = lyap_batch(mod_arcoeff,mod_arnoicse,mod_burnin);

%## self power/auto power
% [power] = autopower(A,Ve,100,200);
%% STEP 6: Compute Connectivity
fprintf('===================================================\n');
fprintf('%10s CONNECTIVITY ESTIMATION\n','');
fprintf('===================================================\n');

%## Coherence metric
% [coherence] = paircoherence(aredat,arndat,n,fs); %coherence()?

%## Power metric
% [power] = autopower(aredat,arndat,n,fs); %spectrum_analysis()?

%## Granger metric
% [Fx2y,Fy2x] = mov_mul_ga(dat,startp,endp,window,order,fs,freq);
%}
end
%% SUBFUNCTIONS

%% DOCUMENTATION PER FUNCTION 
