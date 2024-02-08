function [STUDY,ALLEEG] = as_create_study(ALLEEG,cluster_info_fpath,study_fName,study_fPath,varargin)
%MIM_CREATE_STUDY Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
AMANDA_CLUSTER_ITERS = [1,3:12];
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'cluster_info_fpath',@ischar);
addRequired(p,'study_fName',@ischar);
addRequired(p,'study_fPath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,ALLEEG,cluster_info_fpath,study_fName,study_fPath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
% DIP_NUM = 1;
% DIP_PLOT = 'off';
%## MAKE DIRS
if ~exist(study_fPath,'dir')
    mkdir(study_fPath);
end
%% ===================================================================== %%
% pp = gcp('nocreate');
% disp(pp);
% if ~isfield(pp,'NumWorkers')
%     POOL_SIZE = 1;
% else
%     POOL_SIZE = pp.NumWorkers;
% end
%## LOAD AMANDAS CLUSTER INFO
fprintf('\nLoading Amandas Cluster Information...\n');
tmp = load(cluster_info_fpath);
tmp = tmp.cluster_info;
tmp_cluster = tmp; %tmp(AMANDA_CLUSTER_ITERS);
%- extract component array
comps_out = zeros(length(tmp_cluster),length(ALLEEG));
compList = [];
setList  = [];
for clus_i = 2:length(tmp_cluster)
    sets_i = tmp_cluster(clus_i).sets;
    for j = 1:length(sets_i)
        comps_out(clus_i,sets_i(j)) = tmp_cluster(clus_i).comps(j);
        compList = [compList tmp_cluster(clus_i).comps(j)];
        setList = [setList repmat(sets_i(j),1,length(tmp_cluster(clus_i).comps(j)))];
    end
end
tmp_cluster(1).ursets = tmp_cluster(1).sets;
tmp_cluster(1).urcomps = tmp_cluster(1).comps;
tmp_cluster(1).sets = setList;
tmp_cluster(1).comps = compList;
%% REMOVE COMPS?
%{
tmp_rmv_subjs = zeros(1,length(ALLEEG));
for subj_i = 1:length(ALLEEG)
    tmp = tmp_cluster(1).comps(tmp_cluster(1).sets == subj_i);
    if max(tmp) == size(ALLEEG(subj_i).icaweights,1)
        continue;
    elseif length(tmp) < 3 
        tmp_rmv_subjs(subj_i) = 1;
        continue;
    else
        ALLEEG(subj_i) = pop_subcomp(ALLEEG(subj_i),tmp,0,1);
    end
end
ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
[tmp_cluster] = rmv_subj_inds(tmp_cluster,find(tmp_rmv_subjs));
%}
%% CREATE STUDY
% initiailize study
fprintf('\n==== Making Study Modifications ====\n')
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName,...
                                'filename',study_fName,...
                                'filepath',study_fPath);
%## (AS) DIPOLE STRUCT
STUDY.urcluster = tmp_cluster;
STUDY.cluster = tmp_cluster;
% STUDY.etc.rmvd_subj.inds = tmp_rmv_subjs;
%## SAVE
% parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
% parfor subj_i = 1:length(ALLEEG)
%     ALLEEG(subj_i).etc.full_setfile.filename = ALLEEG(subj_i).filename;
%     ALLEEG(subj_i).etc.full_setfile.filepath = ALLEEG(subj_i).filepath;
%     ALLEEG(subj_i).filename = sprintf('%s_%s_ICA_TMPEEG',ALLEEG(subj_i).subject,'reducedcomps');
%     ALLEEG(subj_i) = pop_saveset(ALLEEG(subj_i),'filename',ALLEEG(subj_i).filename,'filepath',ALLEEG(subj_i).filepath);
% end
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
end
%% SUBFUNCTIONS
function [tmp_cluster] = rmv_subj_inds(tmp_cluster,tmp_rmv_subjs)
tmp_subjs = zeros(length(unique(tmp_cluster(1).sets)),2);
tmp_subjs(:,1) = unique(tmp_cluster(1).sets);
% set_inds = unique(tmpS.cluster(1).sets);
fprintf('==== Reorganizing Cluster Information ====\n');
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
%- remove subjects that got rejected (timewarping problem, <3 brain
%comps,...)
for subj_i = 1:length(tmp_rmv_subjs)
    for cluster_i = 2:length(tmp_cluster)
        inds = (tmp_cluster(cluster_i).sets == tmp_rmv_subjs(subj_i));
        if any(inds)
            tmp_cluster(cluster_i).sets(inds) = [];
            tmp_cluster(cluster_i).comps(inds) = [];
        else
            continue;
        end
    end
end
%- renumber subjects to align cluster struct to ALLEEG struct
for cluster_i = 2:length(tmp_cluster)
    for comp_i = 1:length(tmp_cluster(cluster_i).sets)
        inds = (tmp_subjs(:,1) == tmp_cluster(cluster_i).sets(comp_i));
        if any(inds)
            tmp_cluster(cluster_i).sets(comp_i) = tmp_subjs(inds,2);
        else
            continue;
        end
    end
end
%- parentcluster recomputation
all_sets = [];
all_comps = [];
for clust_i = 2:length(tmp_cluster)
    all_sets = [all_sets, tmp_cluster(clust_i).sets];
    all_comps = [all_comps, tmp_cluster(clust_i).comps];
end
tmp_cluster(1).comps = all_comps;
tmp_cluster(1).sets = all_sets;
end
