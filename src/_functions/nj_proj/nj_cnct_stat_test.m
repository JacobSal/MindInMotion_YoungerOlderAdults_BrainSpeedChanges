function [ALLEEG,stats_store,cond_tests,nonzero_stats,boot_avgs] = as_cnct_stat_test(EEG,varargin)
%CNCTANL_BOOTSTRAP_TEST Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
ALPHA = 0.05;
CONN_METHODS = {};
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
%## OPTIONAL
%## PARAMETER
addParameter(p,'CONN_METHODS',CONN_METHODS,@iscell)
addParameter(p,'ALPHA',ALPHA,@isnumeric)
parse(p, EEG, varargin{:});
CONN_METHODS = p.Results.CONN_METHODS;
ALPHA = p.Results.ALPHA;
%%
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
nonzero_stats = cell(length(EEG.etc.COND_CAT),1);
parfor cond_i = 1:length(EEG.etc.cond_files)
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
    ALLEEG{cond_i}.condition = sprintf('cond_%i',cond_i);
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
stats_store = cell(length(ALLEEG),length(ALLEEG));
cond_tests = cell(length(ALLEEG),length(ALLEEG));
boot_avgs = cell(length(ALLEEG),length(ALLEEG));
done = [];
for cond_i = 1:length(ALLEEG)
    for cond_j = 1:length(ALLEEG)
        if any((cond_j==done)) %|| cond_i==cond_j
            continue;
        end
        tmp = [ALLEEG(cond_i);ALLEEG(cond_j)];
        %- Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
        [stats,boot_avg,~] = stat_surrogateStats('ALLEEG',tmp,...
            'statTest',{'Hab',...
                'tail','both',... %[left|right|one|both]
                'testMethod','quantile',...
                'computeci',true,...
                'alpha',ALPHA,...
                'mcorrection','fdr',...
                'statcondargs',{'mode','perm'}},...
            'connmethods',CONN_METHODS,...
            'VerbosityLevel',1);
        stats_store{cond_i,cond_j} = stats;
        cond_tests{cond_i,cond_j} = {ALLEEG(cond_i).condition,ALLEEG(cond_j).condition};
        boot_avgs{cond_i,cond_j} = boot_avg;
    end
    done = [done cond_i];
end
for cond_i = 1:length(ALLEEG)
    ALLEEG(cond_i).CAT.PConn = [];
end
%## TIME
toc

end

