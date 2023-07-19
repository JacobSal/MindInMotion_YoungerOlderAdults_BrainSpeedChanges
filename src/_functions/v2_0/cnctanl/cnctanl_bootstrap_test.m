function [ALLEEG,Stats,nonzero_stats] = cnctanl_bootstrap_test(EEG,varargin)
%CNCTANL_BOOTSTRAP_TEST Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
ALPHA = 0.05;
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
%## OPTIONAL
%## PARAMETER
parse(p, EEG, varargin{:});
%%
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
nonzero_stats = cell(length(EEG.etc.COND_CAT),1);
for cond_i = 1:length(EEG.etc.cond_files)
    if ispc
        fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
    else
        fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
    end
    fName = EEG.etc.cond_files(cond_i).fName;
    ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
    ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
    %- Phase Randomization Permutation Test Data Handler
    fprintf('\n==== LOADING BOOTSTRAPPED CONNECTIVITY MEASURES ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') 
        error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
    else
        ALLEEG{cond_i}.CAT.PConn = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
    end
    fprintf('done.\n')
    %- Nonzero Statistics Data Handler
    fprintf('\n==== LOADING NONZERO STATISTICS ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file')
        error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_NonZero.mat']]);
    else
        nonzero_stats{cond_i} = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
    end
end
ALLEEG = cellfun(@(x) [[],x],ALLEEG);
%## 1) Between-condition test:
%     For conditions A and B, the null hypothesis is either
%     A(i,j)<=B(i,j), for a one-sided test, or
%     A(i,j)=B(i,j), for a two-sided test
%     A p-value for rejection of the null hypothesis can be
%     obtained by taking the difference of the distributions
%     computing the probability
%     that a sample from the difference distribution is non-zero
if length(ALLEEG)<2
    error('You need two datasets to compute between-condition statistics')
end
% Note this function will return a new EEG dataset with the condition
% differences (Set A - Set B) in the order specified in datasetOrder
fprintf('\n===================================================\n');
disp('Between Condition Test')
fprintf('===================================================\n');
for cond_i = 1:length(ALLEEG)
    ALLEEG(cond_i).CAT.Stats = [];
end
cfg = [];
cfg.statTest = 'Hab';
%         cfg.statTest.tail = [];
%         cfg.connmethods = {'dDTF08'}; % (default: 'all')
cfg.verb = 1;
[Stats, ~, cfg] = feval(@stat_surrogateStats,'ALLEEG',ALLEEG,cfg);
%## Significance Masking
%{
connmethods = ALLEEG(1).CAT.configs.est_mvarConnectivity.connmethods;
for cond_i = 1:length(ALLEEG)
    ALLEEG(cond_i).CAT.Stats = Stats;
    for conn_i = 1:length(connmethods)
        %## DEBUG
        %{
        i = 2;
        j = 6;
        bs_tmp_p = squeeze(ALLEEG(cond_i).CAT.Stats.(connmethods{conn_i}).pval(i,j,:,:));
        nz_tmp_p = squeeze(nonzero_stats{cond_i}.(connmethods{conn_i}).pval(i,j,:,:));
        nz_tmp_th = squeeze(nonzero_stats{cond_i}.(connmethods{conn_i}).thresh(i,j,:,:));
        bs_tmp_ci = squeeze(ALLEEG(cond_i).CAT.Stats.(connmethods{conn_i}).ci(i,j,:,:,:));
        %}
        vals = ALLEEG(cond_i).CAT.Conn.(connmethods{conn_i});
        tmp_1 = ALLEEG(cond_i).CAT.Stats.(connmethods{conn_i}).pval(:,:,:,:)<ALPHA;
        tmp_2 = nonzero_stats{cond_i}.(connmethods{conn_i}).pval<ALPHA;
        vals = vals.*tmp_1.*tmp_2;
        ALLEEG(cond_i).CAT.masked_conn.(connmethods{conn_i}) = vals;
        ALLEEG(cond_i).CAT.masked_conn.alpha = ALPHA;
%         EEG.CAT.Stats.(connmethods{conn_i}).pval(:,:,:,:) = tmp;
    end
end
%}
%(08/07/2022),JS: convert to base stat_ functions

% Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
% The dimensionality is [num_vars x num_vars x num_freqs x num_times]
end

