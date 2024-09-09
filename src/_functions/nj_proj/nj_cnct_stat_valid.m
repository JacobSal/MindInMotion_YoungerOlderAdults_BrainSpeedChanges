function [conn_store_subj] = nj_cnct_stat_valid(EEG,STUDY,cluster_iters,cluster_names,component_vec,save_dir,varargin)
    %AS_CNCT_STAT_VALID Summary of this function goes here
    %   Detailed explanation goes here
    %## TIME
    tic
    %## DEFINE DEFAULTS
    SAVE_FIG = true;
    ALPHA = 0.05;
    CONN_METHODS = {};
    %- PLOT DEFAULTS
    %* color limits
    CLIM = [0,0.005];
    %* frequency scaling
    FREQSCALE = 'log';
    PLOT_CI = true;
    %## 
    p = inputParser;
    %## REQUIRED
    addRequired(p,'EEG',@isstruct)
    addRequired(p,'STUDY',@isstruct);
    addRequired(p,'cluster_iters',@isnumeric)
    addRequired(p,'cluster_names',@iscell)
    addRequired(p,'component_vec',@isnumeric);
    addRequired(p,'save_dir',@ischar);
    %## OPTIONAL
    %## PARAMETER
    addParameter(p,'ALPHA',ALPHA,@isnumeric)
    addParameter(p,'CONN_METHODS',CONN_METHODS,@iscell);
    parse(p,EEG,STUDY,cluster_iters,cluster_names,component_vec,save_dir,varargin{:});
    ALPHA = p.Results.ALPHA;
    CONN_METHODS = p.Results.CONN_METHODS;
    %##
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
    %% ===================================================================== %%
    %## Set Variables
    STAT_CHAR = {'bootstrap','nonzeros'};
    stat_dim = 2;
    eeg_dim = 2;
    conn_dim =  length(EEG.etc.COND_CAT);
    freq_dim = length(EEG.etc.COND_CAT(1).Conn.freqs);
    time_dim = length(EEG.etc.COND_CAT(1).Conn.winCenterTimes);
    conn_store_subj = nan(length(STUDY.cluster),length(STUDY.cluster),freq_dim,time_dim,conn_dim,conn_dim,stat_dim,eeg_dim);
    %## Gather Statistic Tests
    [EEG,stats_bootstrap,cond_bootstrap,stats_nonzero] = nj_cnct_stat_test(EEG,...
        'CONN_METHODS',CONN_METHODS,...
        'ALPHA',ALPHA);
    %## Plot Validations & Store Values
    for stat_i = 1:2
        done = [];
        for cond_i = 1:length(EEG)
            for cond_j = 1:length(EEG)
                if any((cond_j == done))
                    continue;
                end
                fprintf('Plotting condition test %s-%s\n',cond_bootstrap{cond_i,cond_j}{1},cond_bootstrap{cond_i,cond_j}{2})
                for eeg_i = 1:length(EEG)
                    if stat_i == 1
                        EEG(eeg_i).CAT.Stats = stats_bootstrap{cond_i,cond_j};
                    else
                        EEG(eeg_i).CAT.Stats = stats_nonzero{cond_i};
                    end
                end
                %## LOOP MEAT
                %- generate display names based on CLUSTER_ASSIGNMENTS
                orig_comps = squeeze(component_vec);
                [tmpcl,idxcl] = sort(orig_comps);
                display_names = cell(length(orig_comps),1);
                for i = 1:length(orig_comps)
                    if orig_comps(i) ~= 0
                        display_names{i} = sprintf('%s_ic%i',cluster_names{(i == cluster_iters)},orig_comps(i));
                    end
                end
                display_names = display_names(idxcl);
                display_names = display_names(~cellfun(@isempty,display_names));
                clust_inds = idxcl(tmpcl~=0);
                %## Perform Condition T-tests with Bootstrapped Distributions
                tmp = [EEG(cond_i);EEG(cond_j)];
                for eeg_i = 1:length(tmp)
                    [~,~,~,ConnMatrix] = jsedit_vis_TimeFreqGrid('ALLEEG',tmp,'Conn',tmp(eeg_i).CAT.Conn,...
                        'plotCondDiff',{'condOrder',cond_bootstrap{cond_i,cond_j}},...
                        'stats',tmp(eeg_i).CAT.Stats,...
                        'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
                        'msubset','all',...
                        'MatrixLayout',{'Full','estimator',CONN_METHODS{1},'clim',CLIM},...
                        'thresholding',{'Statistics','plotci',PLOT_CI,'sigthreshmethod','pval','alpha',ALPHA},...
                        'pcontour',{true},... %'contourcolor',[0,0,0]
                        'transform',FREQSCALE,...
                        'freqscale',FREQSCALE,... %'yord',{'4','7','12','28','48','60'},...
                        'events',{{0,'r',':',2}},...
                        'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
                        'FrequencyMarkerColor',[0,0,0],...
                        'backgroundColor',[1,1,1],...
                        'textColor',[0,0,0],...
                        'linecolor',[0,0,0],...
                        'patchcolor',[0,0,0],...
                        'NodeLabels',display_names,...
                        'axesFontSize',10,...
                        'topoplot','Topoplot'); %,'estimator',CONN_METHODS
                    %- plot edits
                    fig_i = get(groot,'CurrentFigure');
                    set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
    %                 for i = 1:length(display_names)
    %                     set(fig_i.Children(i),'YTick',{0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307})
    %                     set(fig_i.Children(i),'YTickLabel',{'1','7','18','29','47','76'})
    %                 end
                    %- turn zeros into 'nan'
                    tmp_conn = ConnMatrix;
                    tmp_conn(tmp_conn == 0) = nan();
                    %- store stats
                    for j = 1:size(ConnMatrix,1)
                        for k = 1:size(ConnMatrix,2)
                            conn_store_subj(clust_inds(j),clust_inds(k),:,:,cond_i,cond_j,stat_i,eeg_i) = squeeze(tmp_conn(j,k,:,:));
                        end
                    end
                    %- save figures
                    if SAVE_FIG
                        saveas(fig_i,[save_dir filesep sprintf('%s_valsshown-%s_TimeFreqChart_%s-%s.fig',STAT_CHAR{stat_i},tmp(eeg_i).condition,cond_bootstrap{cond_i,cond_j}{1},cond_bootstrap{cond_i,cond_j}{2})]);
                        saveas(fig_i,[save_dir filesep sprintf('%s_valsshown-%s_TimeFreqChart_%s-%s.jpg',STAT_CHAR{stat_i},tmp(eeg_i).condition,cond_bootstrap{cond_i,cond_j}{1},cond_bootstrap{cond_i,cond_j}{2})]);
                        close(fig_i)
                    end
                end
            end
            done = [done cond_i];
        end
    end
    end
    
    