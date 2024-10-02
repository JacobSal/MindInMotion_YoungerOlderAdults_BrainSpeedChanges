function [TMP_ALLEEG,STUDY] = parfunc_rmv_subjs(TMP_ALLEEG,STUDY,rmv_subjs,varargin)
%EEGLAB_SUBCOMP Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- cell of alleegs
errorMsg = 'Value must be a CELL ARRAY of EEG structures (EEGLAB).'; 
ta_validFcn = @(x) assert(iscell(x),errorMsg);
%- array of subject
errorMsg = 'Value must be a vector of 1''s and 0''s with length equal to length(ALLEEG).'; 
as_validFcn = @(x) assert(isnumeric(x),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'TMP_ALLEEG',ta_validFcn);
addRequired(p,'STUDY',@isstruct);
addRequired(p,'rmv_subjs',as_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,TMP_ALLEEG,STUDY,rmv_subjs,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
fprintf('==== Reformatting Study ====\n');
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(TMP_ALLEEG));
for subj_i = 1:length(TMP_ALLEEG)
    fss{subj_i} = fields(TMP_ALLEEG{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(TMP_ALLEEG)
    EEG = TMP_ALLEEG{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
        end
    end 
    TMP_ALLEEG{subj_i} = EEG;
end
TMP_ALLEEG = cellfun(@(x) [[]; x], TMP_ALLEEG);
%## CREATE NEW STUDY STRUCTURED
fprintf('==== Removing Subjects From Study ====\n');
[tmp_study, TMP_ALLEEG] = std_editset(STUDY,TMP_ALLEEG,...
                            'addchannellabels','off',...
                            'commands',{'remove',find(rmv_subjs)});
tmp_study = std_maketrialinfo(tmp_study, TMP_ALLEEG);
%% ===================================================================== %%
%## Study Modifications
%- NOTE: it seems that when you use 'commands' option in std_editset the
%STUDY.cluster strucutre will be cleared every time you save :( . Need to
%override later on in the pipeline or make some weird edits to STUDY
%struct... for a later date (04/28/2023).
% disp(MAIN_STUDY.cluster); %## DEBUG
% disp(tmp_study.cluster); %## DEBUG
%- using 'commands' in std_editset sets STUDY.cluster = [] & STUDY.changrp = [];
tmp_study.etc.changrp_save = STUDY.changrp;
tmp_study.etc.cluster_save = STUDY.cluster;
%- removed subject information
tmp_study.etc.rmvd_subj.inds = find(rmv_subjs);
tmp_study.etc.rmvd_subj.chars = {TMP_ALLEEG(rmv_subjs).subject};
%## 
tmp_subjs = zeros(length(unique(tmp_study.cluster(1).sets)),2);
tmp_subjs(:,1) = unique(tmp_study.cluster(1).sets);
tmp_rmv_subjs = tmp_study.etc.rmvd_subj.inds;
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
    for cluster_i = 2:length(tmp_study.cluster)
        inds = (tmp_study.cluster(cluster_i).sets == tmp_rmv_subjs(subj_i));
        if any(inds)
            tmp_study.cluster(cluster_i).sets(inds) = [];
            tmp_study.cluster(cluster_i).comps(inds) = [];
        else
            continue;
        end
    end
end
%- renumber subjects to align cluster struct to ALLEEG struct
for cluster_i = 2:length(tmp_study.cluster)
    for comp_i = 1:length(tmp_study.cluster(cluster_i).sets)
        inds = (tmp_subjs(:,1) == tmp_study.cluster(cluster_i).sets(comp_i));
        if any(inds)
            tmp_study.cluster(cluster_i).sets(comp_i) = tmp_subjs(inds,2);
        else
            continue;
        end
    end
end
%- parentcluster recomputation
all_sets = [];
all_comps = [];
for clust_i = 2:length(tmp_study.cluster)
    all_sets = [all_sets, tmp_study.cluster(clust_i).sets];
    all_comps = [all_comps, tmp_study.cluster(clust_i).comps];
end
tmp_study.cluster(1).comps = all_comps;
tmp_study.cluster(1).sets = all_sets;
%- outputs
STUDY = tmp_study;
%## TOC
toc
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