function [STUDY,comps_store] = cluster_ica_reduce(STUDY,varargin)
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
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%% ===================================================================== %%
comps_store = zeros(length(STUDY.cluster),length(STUDY.datasetinfo));
fprintf('==== Choosing best component for each cluster ====\n');
for cluster_i = 2:length(STUDY.cluster)    
    %- subset sets and comps
    tmpsets = unique(STUDY.cluster(cluster_i).sets);
    for subj_i = tmpsets
        idx = (STUDY.cluster(cluster_i).sets == subj_i);
        comps_clust = STUDY.cluster(cluster_i).comps(idx);
        %- just choosing a component based on the AMICA ICA sorting algorithm
        if ~isempty(comps_clust) && length(comps_clust) > 1
            comps_store(cluster_i,subj_i) = min(comps_clust);
            fprintf('Cluster %i) Subject %i: using component %i...\n',cluster_i,subj_i,min(comps_clust));
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

