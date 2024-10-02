function [EEG,cfg] = cnctanl_groupStats(EEG,stat,varargin)
%GROUPSTATSCNCTANL Summary of this function goes here
%   Detailed explanation goes here
%   

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {'nogui',...
    1,...
    []};
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'stat',@ischar)
%## OPTIONAL
%## PARAMETER
parse(p, EEG, stat, varargin{:});



%% next we compute p-values and confidence intervals
% (CHOOSE ONE OF THE FOLLOWING)
switch stat
    case 'BootStrap'
        N_PERMS = 200;
        %% OPTIONAL STEP 8: Compute Statistics (This step is slow)
        % see. stat_surrogateGen.m
        fprintf('\n===================================================\n');
        disp('BOOTSTRAP STATISTICS')
        fprintf('===================================================\n');
        % number of bootstrap samples to draw
        EEG.CAT.PConn  = [];
        cfg = [];
        cfg.verb = 1;
        cfg.mode = {'Bootstrap','nperms',N_PERMS,'saveTrialIdx',false};
        [PConn, cfg] = feval(@stat_surrogateGen,'ALLEEG',EEG,cfg);
        EEG.CAT.PConn = PConn;
        % Bootstrap distributions are now stored in EEG.CAT.PConn
        % This is a structure containing bootstrap estimates of the connectivity,
        % etc., stored in matrices of size:
        % [num_vars x num_vars x num_freqs x num_times x num_samples]

        % we can also replace connectivity estimate with bootstrap estimate
        
    case 'PhaseRnd'
        
        %% NOTE: alternately, we can also obtain the phase-randomized null distributions for each condition
        % Create surrogate data by randomizing phases.
        % Multiply each fourier amplitude by e^{iw)
        % where w is a random phase chosen in [0 2pi]
        % (c.f. Theiler, et al 1997)
        % To generate hermitian phase distributions, we extract
        % the phase of a random matrix. This ensures the surrogate
        % spectrum is conjugate symmetric
        fprintf('\n===================================================\n');
        disp('PHASE RANDOMIZED BOOTSRAP')
        fprintf('===================================================\n');
        %- clear PConn
        EEG.CAT.PConn  = [];
        %- 
        N_PERMS = 200; % number of null distribution samples to draw
        %- set configuration structure
        cfg = [];
        cfg.mode.arg_direct = 1;
        cfg.mode.nperms = N_PERMS;
        cfg.mode.arg_selection = 'PhaseRand';
        %-
        cfg.modelingApproach = EEG.CAT.configs.est_fitMVAR;
        cfg.connectivityModeling = EEG.CAT.configs.est_mvarConnectivity;
        cfg.verb = 1;
        %- FEVAL
        [PConn, cfg] = feval(@stat_surrogateGen,'ALLEEG',EEG,cfg);
        EEG.CAT.PConn = PConn;
        
        % Phase randomized distributions are now stored in EEG.CAT.PConn
        % This is a structure containing estimates of the connectivity etc., under 
        % the null hypothesis of no connectivity (random phase).
        % The distributions are stored in matrices of size:
        % [num_vars x num_vars x num_freqs x num_times x num_samples]

    case 'BtwnCond'
        %% 1) Between-condition test:
        %     For conditions A and B, the null hypothesis is either
        %     A(i,j)<=B(i,j), for a one-sided test, or
        %     A(i,j)=B(i,j), for a two-sided test
        %     A p-value for rejection of the null hypothesis can be
        %     obtained by taking the difference of the distributions
        %     computing the probability
        %     that a sample from the difference distribution is non-zero
        if length(EEG)<2
            error('You need two datasets to compute between-condition statistics')
        end
        % Note this function will return a new EEG dataset with the condition
        % differences (Set A - Set B) in the order specified in datasetOrder
        fprintf('\n===================================================\n');
        disp('Between Condition Test')
        fprintf('===================================================\n');
        for i = 1:length(EEG)
            EEG(i).CAT.Stats = [];
        end
        cfg = [];
        cfg.statTest = 'Hab';
%         cfg.statTest.tail = [];
%         cfg.connmethods = {'dDTF08'}; % (default: 'all')
        cfg.verb = 1;
        [Stats, ~, cfg] = feval(@stat_surrogateStats,'ALLEEG',EEG,cfg);
        for i = 1:length(EEG)
            EEG(i).CAT.Stats = Stats;
        end
        %(08/07/2022),JS: convert to base stat_ functions
        
        % Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
    case 'DevFromBase'
        %% 2) Devation from baseline test
        %     For conditions A, the null hypothesis is
        %     C(i,j)=baseline_mean(C). This is a two-sided test.
        %     A p-value for rejection of the null hypothesis can be
        %     obtained by obtaining the distribution of the difference from
        %     baseline mean and computing the probability
        %     that a sample from this distribution is non-zero

        %(08/07/2022),JS: convert to base stat_ functions
        
        % Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
    case 'NonZero'
        %% 3) Test for non-zero connectivity
        %     We are testing with respect to a phase-randomized null
        %     distribution. A p-value for rejection of the null hypothesis
        %     can be obtained by computing the probability that the
        %     observed connectivity is a random sample from the null distribution
        fprintf('\n===================================================\n');
        disp('NonZero Test')
        fprintf('===================================================\n');
        %-
        EEG.CAT.Stats = [];
        %- 
        STAT_ALPHA = 0.05; 
        STAT_MCORRECTION = 'fdr';
        STAT_TAIL = 'right';
        STAT_TESTMETHOD = 'quantile';
        %-
        cfg = [];
        cfg.statTest =  {'Hnull',...
                         'testMethod',STAT_TESTMETHOD,...
                         'tail',STAT_TAIL,...
                         'alpha',STAT_ALPHA,...
                         'mcorrection',STAT_MCORRECTION,...
                         'statcondargs',{}};
%         cfg.connmethods = {'dDTF08'}; % (default: 'all')
        cfg.verb = 1;
        [Stats, ~, cfg] = feval(@stat_surrogateStats,'ALLEEG',EEG,cfg);
        EEG.CAT.Stats = Stats;
%         EEG.CAT.Conn = ConnNew;
        % Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
    case 'AnlStat'
        %% OPTIONAL STEP 8b: Compute analytic statistics
        % This computes analytic alpha-significance thresholds, p-values, and confidence
        % intervals for select connectivity estimators (RPDC, nPDC).
        % These are asymptotic estimators and may not be accurate for small sample
        % sizes. However, they are very fast and usually a reasonable estimate.
        
        %(08/07/2022),JS: convert to base stat_ functions

        % Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
end
%% SUBFUNCTIONS
end

