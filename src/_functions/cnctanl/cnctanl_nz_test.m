function [ALLEEG,phasernd_maskedconn] = cnctanl_nz_test(EEG,varargin)
%CNCTANL_BOOTSTRAP_TEST Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
ALPHA = 0.05;
CONN_METHODS = {};
SAVE_DIR = EEG.filepath;
DISPLAYNAMES = {};
DO_PLOT = true;
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
addParameter(p,'DO_PLOT',DO_PLOT,@islogical);
parse(p, EEG, varargin{:});
CONN_METHODS = p.Results.CONN_METHODS;
ALPHA = p.Results.ALPHA;
DO_PLOT = p.Results.DO_PLOT;
SAVE_DIR = p.Results.SAVE_DIR;
DISPLAYNAMES = p.Results.DISPLAYNAMES;
%% ===================================================================== %%
ALLEEG = cell(length(EEG.etc.COND_CAT),1);
nonzero_stats = cell(length(EEG.etc.COND_CAT),1);
phasernd_pconn = cell(length(EEG.etc.COND_CAT),1);
phasernd_maskedconn = cell(length(EEG.etc.COND_CAT),length(CONN_METHODS));
parfor cond_i = 1:length(EEG.etc.cond_files)
    if ispc
        fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
    else
        fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
    end
    fName = EEG.etc.cond_files(cond_i).fName;
    ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
    ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
    %- PhaseRnd test
    fprintf('\n==== LOADING PHASERAND STATISTICS ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file')
        error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
    else
        phasernd_pconn{cond_i} = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
    end
    %- condition override
    ALLEEG{cond_i}.condition = sprintf('cond_%i',cond_i);
end
ALLEEG = cellfun(@(x) [[],x],ALLEEG);
%## 0) Make unthresheld plot
if DO_PLOT
    for cond_i = 1:length(ALLEEG)
        Stats = [];
        for meth_i = 1:length(CONN_METHODS)
            tf_plot(ALLEEG(cond_i),{ALLEEG(cond_i).condition},Stats,ALLEEG(cond_i).CAT.Conn,CONN_METHODS{meth_i},'nonzero','nothresh',DISPLAYNAMES,SAVE_DIR);
        end
    end
end
%## 1) Nonzero test for phase randomized data
for cond_i = 1:length(ALLEEG)
    %% 3) Test for non-zero connectivity
    %     We are testing with respect to a phase-randomized null
    %     distribution. A p-value for rejection of the null hypothesis
    %     can be obtained by computing the probability that the
    %     observed connectivity is a random sample from the null distribution
    fprintf('\n===================================================\n');
    disp('NonZero Test')
    fprintf('===================================================\n');
    %-
    ALLEEG(cond_i).CAT.PConn = phasernd_pconn{cond_i};
    %-
    ALLEEG(cond_i).CAT.Stats = [];
    %-
    [Stats,~,~] = stat_surrogateStats('ALLEEG',ALLEEG(cond_i),...
                     'statTest',{'Hnull',...
                        'tail','both',... %[left|right|one|both]
                        'testMethod','quantile',...
                        'computeci',true,...
                        'alpha',ALPHA,...
                        'mcorrection','fdr',...
                        'statcondargs',{'mode','perm'}},...
                    'connmethods',CONN_METHODS,...
                    'VerbosityLevel',1);
    ALLEEG(cond_i).CAT.Stats = Stats;
    nonzero_stats{cond_i} = Stats;
    for meth_i = 1:length(CONN_METHODS)
        if DO_PLOT
            [conn_masked] = tf_plot(ALLEEG(cond_i),{ALLEEG(cond_i).condition},...
                Stats,ALLEEG(cond_i).CAT.Conn,CONN_METHODS{meth_i},'nonzero',...
                'thresh',DISPLAYNAMES,SAVE_DIR);
            if isfield(ALLEEG(cond_i).CAT,'Stats')
                stat_thresh = ALLEEG(cond_i).CAT.Stats.(CONN_METHODS{meth_i}).pval < ALPHA;
            end
            %- connectivity extraction
    %         conn_masked_2  = squeeze(ALLEEG(cond_i).CAT.Conn.(CONN_METHODS{meth_i})(:,:,:,:)).*squeeze(stat_thresh(:,:,:,:));
    %         disp(squeeze(conn_masked_2(1,3,:,:)));
    %         disp(squeeze(conn_masked(1,3,:,:)));
            phasernd_maskedconn{cond_i,meth_i} = conn_masked; %conn_masked
        end
    end
    close all
end
for cond_i = 1:length(ALLEEG)
    ALLEEG(cond_i).CAT.PConn = [];
    ALLEEG(cond_i).CAT.Stats = [];
end
% Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
% The dimensionality is [num_vars x num_vars x num_freqs x num_times]
end
%% ===================================================================== %%
function [conn_store] = tf_plot(TMP_EEG,condOrder,stats,conn,conn_method,...
            stat_name,fig_tag,display_names,save_dir)
    CLIM = [0,0.005];
    ALPHA = 0.05;
    PLOT_CI = false;
    FREQSCALE = 'log';
    SAVE_FIG = true;
%     dims = [size(conn.(conn_method)),2];
%     conn_store = zeros(dims);
    [~,~,~,conn_store] = jsedit_vis_TimeFreqGrid('ALLEEG',TMP_EEG,'Conn',conn,...
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
%     conn_store(:,:,:,:) = new_conn;
    %- save figures
    if SAVE_FIG
        if length(condOrder) > 1
    %             saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s-%s.fig',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1},condOrder{2})]);
            saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s-%s.jpg',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(1).condition,condOrder{1},condOrder{2})]);
        else
    %             saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s.fig',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(eeg_i).condition,condOrder{1})]);
            saveas(fig_i,[save_dir filesep sprintf('%s_%s_%s_valsshown-%s_TimeFreqChart_%s.jpg',fig_tag,TMP_EEG(1).subject,stat_name,TMP_EEG(1).condition,condOrder{1})]);
        end
    %         close(fig_i)
    end
end