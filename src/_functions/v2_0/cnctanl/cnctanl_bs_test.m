function [ALLEEG,bootstrap_masekdconn,boot_avgs] = cnctanl_bs_test(EEG,varargin)
%CNCTANL_BOOTSTRAP_TEST Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
ALPHA = 0.05;
CONN_METHODS = {};
SAVE_DIR = EEG.filepath;
DISPLAYNAMES = {};
%-
% CLIM = [0,0.005];
% FREQSCALE = 'log';
% DO_NONZEROD_BOOTSTRAP = true;
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
%## OPTIONAL
%## PARAMETER
addParameter(p,'SAVE_DIR',SAVE_DIR,@ischar)
addParameter(p,'CONN_METHODS',CONN_METHODS,@iscell)
addParameter(p,'ALPHA',ALPHA,@isnumeric)
addParameter(p,'DISPLAYNAMES',DISPLAYNAMES,@iscell);
parse(p, EEG, varargin{:});
CONN_METHODS = p.Results.CONN_METHODS;
ALPHA = p.Results.ALPHA;
SAVE_DIR = p.Results.SAVE_DIR;
DISPLAYNAMES = p.Results.DISPLAYNAMES;
%% ===================================================================== %%
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
bootstrap_pconn = cell(length(EEG.etc.COND_CAT),1);
bootstrap_masekdconn = cell(length(EEG.etc.COND_CAT),length(EEG.etc.COND_CAT),length(CONN_METHODS));
% parfor cond_i = 1:length(EEG.etc.cond_files)
for cond_i = 1:length(EEG.etc.cond_files)
    if ispc
        fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
    else
        fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
    end
    fName = EEG.etc.cond_files(cond_i).fName;
    ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
    ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
    %- Bootstrap Test Data Handler
    fprintf('\n==== LOADING BOOTSTRAPPED CONNECTIVITY MEASURES ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') 
        error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
    else
        bootstrap_pconn{cond_i} = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
    end
    fprintf('done.\n')
    %- condition override
%     ALLEEG{cond_i}.condition = sprintf('cond_%i',cond_i);
end
ALLEEG = cellfun(@(x) [[],x],ALLEEG);
%## 0) Make unthresheld plot
% for cond_i = 1:length(ALLEEG)
%     Stats = [];
%     for meth_i = 1:length(CONN_METHODS)
%         tf_plot(ALLEEG(cond_i),{ALLEEG(cond_i).condition},Stats,ALLEEG(cond_i).CAT.Conn,CONN_METHODS{meth_i},'nonzero','nothresh',DISPLAYNAMES,SAVE_DIR);
%     end
% end
%## 2) Between-condition test:
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
    ALLEEG(cond_i).CAT.PConn = bootstrap_pconn{cond_i};
end
boot_stats = cell(length(ALLEEG),length(ALLEEG));
cond_tests = cell(length(ALLEEG),length(ALLEEG));
boot_avgs = cell(length(ALLEEG),length(ALLEEG));
done = [];
for cond_i = 1:length(ALLEEG)
    for cond_j = 1:length(ALLEEG)
        if any((cond_j==done)) || cond_i==cond_j
            continue;
        end
        tmp = [ALLEEG(cond_i);ALLEEG(cond_j)];
        %- Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
        [Stats,boot_avg,~] = stat_surrogateStats('ALLEEG',tmp,...
            'statTest',{'Hab',...
                'tail','both',... %[left|right|one|both]
                'testMethod','quantile',...
                'computeci',true,...
                'alpha',ALPHA,...
                'mcorrection','fdr',... %[fdr|bonferonni|numvars]
                'statcondargs',{'mode','perm'}},...
            'connmethods',CONN_METHODS,...
            'VerbosityLevel',1);
        % (08/03/2023) JS, changing mcorrection from fdr to bonferoni to be
        % a little less agressive on removing false positives.
        % (08/03/2023) JS, changing back to fdr.
        boot_stats{cond_i,cond_j} = Stats;
        cond_tests{cond_i,cond_j} = {ALLEEG(cond_i).condition,ALLEEG(cond_j).condition};
        boot_avgs{cond_i,cond_j} = boot_avg;
        for meth_i = 1:length(CONN_METHODS)
            [conn_masked] = tf_plot(tmp,{ALLEEG(cond_i).condition,ALLEEG(cond_j).condition},...
                Stats,boot_avg,CONN_METHODS{meth_i},'bootstrap','thresh',DISPLAYNAMES,SAVE_DIR);
            bootstrap_masekdconn{cond_i,cond_j,meth_i} = conn_masked;
        end
    end
    done = [done cond_i];
end
for cond_i = 1:length(ALLEEG)
    ALLEEG(cond_i).CAT.PConn = [];
    ALLEEG(cond_i).CAT.Stats = [];
end
% Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
% The dimensionality is [num_vars x num_vars x num_freqs x num_times]
end

function [conn_store] = tf_plot(TMP_EEG,condOrder,stats,conn,conn_method,...
    stat_name,fig_tag,display_names,save_dir)
CLIM = [0,0.005];
ALPHA = 0.05;
PLOT_CI = false;
FREQSCALE = 'log';
SAVE_FIG = true;
dims = [size(conn.(conn_method)),2];
conn_store = zeros(dims);
for eeg_i = 1:length(TMP_EEG)
    [~,~,~,new_conn] = jsedit_vis_TimeFreqGrid('ALLEEG',TMP_EEG,'Conn',conn,...
        'plotCondDiff',{'condOrder',condOrder},...
        'stats',stats,...
        'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
        'msubset','all',...
        'MatrixLayout',{'Full','estimator',conn_method,'clim',CLIM},...
        'thresholding',{'Statistics','plotci',PLOT_CI,'sigthreshmethod','pval','alpha',ALPHA},...
        'transform','linear',...
        'freqscale',FREQSCALE,... 
        'NodeLabels',display_names,...
        'events',{{0,'r',':',2}},...
        'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
        'FrequencyMarkerColor',[0,0,0],...
        'backgroundColor',[1,1,1],...
        'textColor',[0,0,0],...
        'linecolor',[0,0,0],...
        'patchcolor',[0,0,0],...
        'axesFontSize',10,...
        'topoplot','Topoplot'); %,'estimator',CONN_METHODS
    %- plot edits
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
    %-
    conn_store(:,:,:,:,eeg_i) = new_conn;
    %- save figures
    if SAVE_FIG
        if length(condOrder) > 1
%             saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s-%s.fig',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1},condOrder{2})]);
            saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s-%s.jpg',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1},condOrder{2})]);
        else
%             saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s.fig',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1})]);
            saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s.jpg',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1})]);
        end
        close(fig_i)
    end
end
end