function [ALLEEG] = cnctanl_nonzero_test(EEG,varargin)
%CNCTANL_NONZERO_TEST Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
%## OPTIONAL
%## PARAMETER
parse(p, EEG, varargin{:});
%%
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
for cond_i = 1:lenght(EEG.etc.cond_files)
    EEG.CAT = EEG.etc.COND_CAT(cond_i);
    %- Phase Randomization Permutation Test Data Handler
    fprintf('\n==== LOADING PHASE RANDOMIZED CONNECTIVITY MEASURES ====\n')
    if ispc
        fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
    else
        fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
    end
    fName = EEG.etc.cond_files(cond_i).fName;
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file') 
        error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
    else
        EEG.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
    end
    fprintf('done.\n')
    %- Nonzero Statistics Data Handler
    fprintf('\n==== LOADING NONZERO STATISTICS ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file')
        error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_NonZero.mat']]);
    else
        EEG.CAT.Stats = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
    end
    fprintf('done.\n')
    ALLEEG{cond_i} = EEG;
end
ALLEEG = cellfun(@(x) [[],x],ALLEEG);
%%
%{
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
for cond_i = 1:length(EEG.etc.cond_files)
    if ispc
        fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
    else
        fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
    end
    fName = EEG.etc.cond_files(cond_i).fName;
    ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
    %- Phase Randomization Permutation Test Data Handler
    fprintf('\n==== LOADING PHASE RANDOMIZED CONNECTIVITY MEASURES ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file') 
        error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
    else
        ALLEEG{cond_i}.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
    end
    fprintf('done.\n')
    %% 3) Test for non-zero connectivity
    %     We are testing with respect to a phase-randomized null
    %     distribution. A p-value for rejection of the null hypothesis
    %     can be obtained by computing the probability that the
    %     observed connectivity is a random sample from the null distribution
    fprintf('\n===================================================\n');
    disp('NonZero Test')
    fprintf('===================================================\n');
    %-
    ALLEEG{cond_i}.CAT.Stats = [];
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
    [Stats, ~, cfg] = feval(@stat_surrogateStats,'ALLEEG',ALLEEG{cond_i},cfg);
    ALLEEG{cond_i}.CAT.Stats = Stats;

    %         EEG.CAT.Conn = ConnNew;
    % Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
    % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
    
end
ALLEEG = cellfun(@(x) [[],x],ALLEEG);
%}

