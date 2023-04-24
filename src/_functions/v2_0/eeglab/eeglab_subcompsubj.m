function [STUDY,ALLEEG] = eeglab_subcompsubj(STUDY,ALLEEG,arr_comps,arr_subjs,varargin)
%EEGLAB_SUBCOMP Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- cell array of components
errorMsg = 'Value must be a cell array with length equal to length(ALLEEG). Each cell contains a vector of indices corresponding to components you want to remove.'; 
ac_validFcn = @(x) assert(iscell(x) && length(x) == length(ALLEEG),errorMsg);
%- array of subject
errorMsg = 'Value must be a vector of 1''s and 0''s with length equal to length(ALLEEG).'; 
as_validFcn = @(x) assert(isnumeric(x) && length(x) == length(ALLEEG),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'STUDY',@isstruct);
addRequired(p,'arr_comps',ac_validFcn);
addRequired(p,'arr_subjs',as_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,arr_comps,arr_subjs,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## REMOVE COMPS
chars_subjs = ALLEEG(logical(arr_subjs)).subject;
for subj_i = 1:length(ALLEEG)
    if length(arr_comps{subj_i}) < 3
        warning('error. %s. pop_subcomp() requires 3 or more components to remain in the EEG structure. Change arr_comps.',ALLEEG(subj_i).subject);
        continue;
    end
    ALLEEG(subj_i) = pop_subcomp(ALLEEG(subj_i),arr_comps{subj_i},0,0);
end
ALLEEG = ALLEEG(~logical(arr_subjs));
% fprintf('\n==== Making Study Modifications ====\n')
[tmp, ALLEEG] = std_editset([],ALLEEG,...
                                'name',study_name);
tmp.cluster = STUDY.cluster;
tmp.design = STUDY.design(~logical(arr_subjs));
tmp.etc.removed_compsubj.subj_inds = find(arr_subjs);
tmp.etc.removed_compsubj.comps = arr_comps;
tmp.etc.removed_compsubj.cluster = STUDY.cluster;
tmp.etc.removed_compsubj.subj_chars = chars_subjs;
%- 
tmp_subjs = zeros(length(unique(tmp.cluster(1).sets)),2);
tmp_subjs(:,1) = unique(tmp.cluster(1).sets);
tmp_rmv_subjs = tmp.etc.removed_compsubj.subj_inds;
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
fprintf('==== removing cluster indices ====\n');
for subj_i = 1:length(tmp_rmv_subjs)
    for cluster_i = 2:length(tmp.cluster)
        inds = (tmp.cluster(cluster_i).sets == tmp_rmv_subjs(subj_i));
        if any(inds)
            tmp.cluster(cluster_i).sets(inds) = [];
            tmp.cluster(cluster_i).comps(inds) = [];
        else
            continue;
        end
    end
end
for cluster_i = 2:length(tmp.cluster)
    for comp_i = 1:length(tmp.cluster(cluster_i).sets)
        inds = (tmp_subjs(:,1) == tmp.cluster(cluster_i).sets(comp_i));
        if any(inds)
            tmp.cluster(cluster_i).sets(comp_i) = tmp_subjs(inds,2);
        else
            continue;
        end
    end
end
%- parentcluster alterations
all_sets = [];
all_comps = [];
for clust_i = 2:length(tmp.cluster)
    all_sets = [all_sets, tmp.cluster(clust_i).sets];
    all_comps = [all_comps, tmp.cluster(clust_i).comps];
end
tmp.cluster(1).comps = all_comps;
tmp.cluster(1).sets = all_sets;


fprintf('==== Reorganizing dipfit indices ====\n');
for subj_i = 1:length(ALLEEG)
    dipfit_i = ALLEEG.dipfit.model;
    vals = [];
    %- extract cluster number and associated component number
    for clust_i = 2:length(tmp.cluster)
        idx = (tmp.cluster(clust_i).sets == subj_i);
        if any(idx)
            vals = [vals; clust_i, find(idx), tmp.cluster(clust_i).comps(idx)];
        end
    end
    %- 
    [~,idx] = sort(vals);
    for val_i = 1:size(vals,1)
        tmp.cluster(vals(idx(val_i,3),1)).comps(vals(idx(val_i,3),2)) = val_i;
    end
    tmp.datasetinfo(subj_i).comps = idx(:,3)';
end
all_sets = [];
all_comps = [];
for clust_i = 2:length(tmp.cluster)
    all_sets = [all_sets, tmp.cluster(clust_i).sets];
    all_comps = [all_comps, tmp.cluster(clust_i).comps];
end
tmp.cluster(1).comps = all_comps;
tmp.cluster(1).sets = all_sets;

%- might not be needed
% tmp.datasetinfo = STUDY.datasetinfo(~logical(arr_subjs)); 
%- study modifications
tmp.etc.rmvd_subj.inds = arr_subjs;
tmp.etc.rmvd_subj.chars = chars_subjs;
end

%{
arr_comps = cell(length(MAIN_ALLEEG));
arr_subjs = zeros(length(MAIN_ALLEEG),1);
rmv_subj_inds = [1,6,10];
for i = 1:length(MAIN_ALLEEG)
    arr_comps{i} = (1:size(MAIN_ALLEEG(i).icaweights,1));
    if any(i == rmv_subj_inds)
        arr_subjs(i) = 1;
        continue;
    end
end
%}