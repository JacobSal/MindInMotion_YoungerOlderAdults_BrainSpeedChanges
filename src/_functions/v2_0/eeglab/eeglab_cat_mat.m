function [conn_mat] = eeglab_cat_mat(STUDY,ALLEEG,varargin)
%EEGLAB_CAT_MAT Summary of this function goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%% Define Defaults
CLUSTER_ITERS = 2:length(STUDY.cluster);
CLUSTER_ASSIGNMENTS = cell(1,length(CLUSTER_ITERS));
for i = 2:length(STUDY.cluster)
    CLUSTER_ASSIGNMENTS{i-1} = num2str(i);
end

%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'cluster_components',@isnumeric);
%## OPTIONAL
addOptional(p,'cluster_iters',CLUSTER_ITERS,@iscell);
addOptional(p,'cluster_assignments',CLUSTER_ASSIGNMENTS,@iscell);
%## PARAMETER
parse(p,STUDY, ALLEEG, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
mat_nan = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),length(MAIN_ALLEEG));
mat_nan_stdv = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),length(MAIN_ALLEEG));
mat_out_nan = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),length(MAIN_ALLEEG),length(load_trials));
mat_out = cell(1,length(load_trials));
mat_out_stdv = cell(1,length(load_trials));
%- extract component array
comps_store = zeros(length(MAIN_ALLEEG),length(MAIN_STUDY.cluster));
for clus_i = 2:length(MAIN_STUDY.cluster)
    sets_i = MAIN_STUDY.cluster(clus_i).sets;
    for j = 1:length(sets_i)
        comps_store(sets_i(j),clus_i) = MAIN_STUDY.cluster(clus_i).comps(j);
    end
end
%- extract connectivity arrays
for cond_i = 1:length(load_trials)
    clusterNames = cell(1,length(CLUSTER_ASSIGNMENTS));
    %- loop through each condition   
    name_trial = load_trials{cond_i};
    %- loop through subjects in each condition
    for subj_i = 1:length(MAIN_ALLEEG)
        chk = (comps_store(subj_i,1:end)>0);
        if any(chk)
            statMat = MAIN_ALLEEG(subj_i).etc.js_processing(cond_i).META.(conn_meas).connStatMatrix;
            comps = MAIN_ALLEEG(subj_i).etc.js_processing(cond_i).META.connComponents;
            meanMat = squeeze(statMat(1,:,:));
%             stdvMat = squeeze(statMat(2,:,:));
%             medMat = squeeze(statMat(3,:,:));
            clust_idx = zeros(1,length(comps));
            comp_idx = zeros(1,length(comps));
            for i = 1:length(comps)
%                 chk = find(comps(i) == comps_store(1:end,subj_i));
                chk = find(comps(i) == comps_store(subj_i,1:end));
                if ~isempty(chk)
                    clust_idx(i) = chk;
                    comp_idx(i) = i;
                end
            end
            clust_idx = clust_idx(clust_idx ~= 0);
            comp_idx = comp_idx(comp_idx ~= 0);
            if ~isempty(comp_idx)
                val_in = meanMat(comp_idx,comp_idx);
                val_in(isnan(val_in)) = 0;
                mat_nan(clust_idx,clust_idx,subj_i) = val_in;
            else
                continue;
            end                    
        else
            continue;
        end
    end
    %- store
    mat_out_nan(:,:,:,cond_i) = mat_nan;
    mat_out{cond_i} = mat_nan; % nanmean(mat_nan,3);
    mat_out_stdv{cond_i} = mat_nan_stdv;
    %- plot
    cnt = 1;
    for i = 1:length(MAIN_STUDY.cluster)
        N = length(MAIN_STUDY.cluster(i).sets);
%         clusterNames{i} = sprintf('(N=%i) Cluster %i',N,i);
        if any(i == CLUSTER_ITERS)
            %- assign name if available
            idx = find(i == CLUSTER_ITERS);
            clusterNames{idx} = sprintf('(N=%i) %s',N,CLUSTER_ASSIGNMENTS{(i == CLUSTER_ITERS)});
            cnt = cnt + 1;            
        end
    end
    
    %## PLOT
    %- assign matrix
%     mat_nan(isnan(mat_nan)) = 0;
%     tmp = mean(mat_nan,3);
    tmp = nanmean(mat_nan,3);
    I = eye(size(tmp));
    I = (I == 0);
    tmp = tmp.*I;
    tmp(tmp == 0) = nan();
%     tmp = log(tmp);
    %- delte unused clusters
    tmp = tmp(CLUSTER_ITERS,:);
    tmp = tmp(:,CLUSTER_ITERS);
    tmp_mat_nan = mat_nan;
    tmp_mat_nan = squeeze(tmp_mat_nan(CLUSTER_ITERS,:,:));
    tmp_mat_nan = squeeze(tmp_mat_nan(:,CLUSTER_ITERS,:));
    %- plot
    figure;
    hnd = heatmap(tmp,'Colormap',jet); %,'CellLabelColor', 'None');
    title(sprintf('''%s'' Connectivity mean across clusters',name_trial));
    hnd.YDisplayLabels = clusterNames;
    hnd.XDisplayLabels = clusterNames;
%         hnd.ColorLimits = [min(min(tmp))-0.2,max(max(tmp))+0.2]; %for log scale
%     hnd.ColorLimits = [-5.5,-2];
    hnd.ColorLimits = [0,0.075];
    hnd.GridVisible = 'off';
    hnd.CellLabelFormat = '%0.1g';
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep sprintf('%s_meanMat_%s.fig','CLUSTERs',name_trial)]);
    saveas(fig_i,[save_dir filesep sprintf('%s_meanMat_%s.jpg','CLUSTERs',name_trial)]);
end
end

