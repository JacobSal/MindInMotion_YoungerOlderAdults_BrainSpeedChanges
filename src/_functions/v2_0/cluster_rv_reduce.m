function [STUDY,comps_store,comps_rej] = cluster_rv_reduce(STUDY,ALLEEG,varargin)
%CLUSTER_ICA_REDUCE Summary of this function goes here
%   Detailed explanation goes here
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

%## TIME
tic
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%% ===================================================================== %%
comps_store = zeros(length(STUDY.cluster),length(STUDY.datasetinfo));
comps_rej = cell(length(STUDY.cluster),length(STUDY.datasetinfo));
fprintf('==== Choosing best component for each cluster ====\n');
for cluster_i = 2:length(STUDY.cluster)    
    %- subset sets and comps
    tmpsets = unique(STUDY.cluster(cluster_i).sets);
    tmps = [];
    tmpc = [];
    for i = 1:length(tmpsets)
        subj_i = tmpsets(i);
        idx = (STUDY.cluster(cluster_i).sets == subj_i);
        comps_clust = STUDY.cluster(cluster_i).comps(idx);
        %- choosing component based on minimum dipole fit residual variance
        if ~isempty(comps_clust) && length(comps_clust) > 1
            [val,idx] = min([ALLEEG(subj_i).dipfit.model(comps_clust).rv]); %min(comps_clust);
            chc = comps_clust(idx);
            %- prints % stores
            fprintf('Cluster %i) Subject %i''s choice component (%i) rv: %i\n',cluster_i,subj_i,chc,val);
            comps_store(cluster_i,subj_i) = chc;
            comps_rej{cluster_i,subj_i} = comps_clust(comps_clust ~= chc);
            fprintf('Cluster %i) Subject %i''s outlier components:',cluster_i,subj_i); fprintf('%i, ',comps_rej{cluster_i,subj_i}); fprintf('\n');
            tmps = [tmps, repmat(subj_i,1,length([comps_rej{cluster_i,subj_i}]))];
            tmpc = [tmpc, comps_rej{cluster_i,subj_i}];
        else
            comps_store(cluster_i,subj_i) = comps_clust;
        end
    end
    %- remove zeros
    tmp = comps_store(cluster_i,:);
    out = tmp(tmp ~= 0);
    %- assign
    STUDY.cluster(cluster_i).sets = tmpsets;
    STUDY.cluster(cluster_i).comps = out;
    %-
    STUDY.cluster(end+1).sets = tmps;
    STUDY.cluster(end).comps = tmpc;
    STUDY.cluster(end).name = sprintf('Outlier clust_%i',cluster_i);
    STUDY.cluster(end).parent = STUDY.cluster(cluster_i).parent;
    STUDY.cluster(end).algorithm = {'minimum residual variance of subject components'};
end
%- parentcluster alterations
all_sets = [];
all_comps = [];
for clust_i = 2:length(STUDY.cluster)
    all_sets = [all_sets, STUDY.cluster(clust_i).sets];
    all_comps = [all_comps, STUDY.cluster(clust_i).comps];
end
STUDY.cluster(1).comps = all_comps;
STUDY.cluster(1).sets = all_sets;
fprintf('done.\n');
%## TIME
toc
end
