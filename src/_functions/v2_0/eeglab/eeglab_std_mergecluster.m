function [merge_cluster_struct] = eeglab_std_mergecluster(cluster_struct_1,cluster_struct_2,varargin)
%EEGLAB_STD_MERGECLUSTER Summary of this function goes here
%   Detailed explanation goes here%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (04/21/2023)
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'cluster_struct_1',@isstruct);
addRequired(p,'cluster_struct_2',@isstruct);
%## OPTIONAL
%## PARAMETER
parse(p,cluster_struct_1,cluster_struct_2,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%- get unique subjects for each
unq_sbj1 = unique(cluster_struct_1(1).sets);
unq_sbj2 = unique(cluster_struct_2(1).sets);
%-
rmv_subjs = setdiff(unq_sbj1,unq_sbj2);
tmp_subjs = zeros(length(unique(STUDY.cluster(1).sets)),2);
tmp_subjs(:,1) = unique(STUDY.cluster(1).sets);
rmv_subjs = STUDY.etc.rmvd_subj.inds;
% set_inds = unique(tmpS.cluster(1).sets);
iter = 1;
for subj_i = 1:length(tmp_subjs)
    if any(subj_i == tmp_rmv_subjs)
        continue;
    else
        tmp_subjs(subj_i,2) = iter;
        iter = iter + 1;
    end
end
fprintf('==== removing subjects from cluster indices ====\n');
%- create an unscrambling array for removing subjects
iter = 1;
for subj_i = 1:length(tmp_subjs)
    if any(subj_i == tmp_rmv_subjs)
        continue;
    else
        tmp_subjs(subj_i,2) = iter;
        iter = iter + 1;
    end
end
for subj_i = 1:length(tmp_rmv_subjs)
    for cluster_i = 2:length(STUDY.cluster)
        inds = (STUDY.cluster(cluster_i).sets == tmp_rmv_subjs(subj_i));
        if any(inds)
            STUDY.cluster(cluster_i).sets(inds) = [];
            STUDY.cluster(cluster_i).comps(inds) = [];
        else
            continue;
        end
    end
end
for cluster_i = 2:length(STUDY.cluster)
    for comp_i = 1:length(STUDY.cluster(cluster_i).sets)
        inds = (tmp_subjs(:,1) == STUDY.cluster(cluster_i).sets(comp_i));
        if any(inds)
            STUDY.cluster(cluster_i).sets(comp_i) = tmp_subjs(inds,2);
        else
            continue;
        end
    end
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
end

