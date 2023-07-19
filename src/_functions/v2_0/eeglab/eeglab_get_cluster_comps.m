function [comps_out,main_cl_inds,outlier_cl_inds,valid_clusters] = eeglab_get_cluster_comps(STUDY,varargin)
%EEGLAB_GET_CLUSTER_COMPS Summary of this function goes here
%   Function generates a matrix (NxM, N = number of clusters, M = number of
%   subjects) of DOUBLES. Each element in the matrix is the component
%   number for each i'th subject in STUDY.datasetinfo/ALLEEG.
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## (PARSER) DEFINE DEFAULTS

%## (PARSER) Define Parser
p = inputParser;
%## (PARSER) REQUIRED
addRequired(p,'STUDY',@isstruct)
%## (PARSER) OPTIONAL
%## (PARSER) PARAMETER
parse(p,STUDY,varargin{:});
%## (PARSER) SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%## (PARSER) MAKE DIRS
%% ===================================================================== %%
%## Extract components for each cluster & subject
%- PARAMS
comps_out = zeros(length(STUDY.cluster),length(STUDY.datasetinfo));
main_cl_inds = zeros(1,length(STUDY.cluster));
outlier_cl_inds = zeros(1,length(STUDY.cluster));
%- loop through all clusters but the parent cluster
main_cl_inds(1) = 1;
for clus_i = 2:length(STUDY.cluster)
    %- if cluster name contains the word 'Outlier' or 'Reject' then the cluster is an
    %outlier
    chk = isempty(regexpi(STUDY.cluster(clus_i).name,'outlier')) && isempty(regexpi(STUDY.cluster(clus_i).name,'reject'));
    if chk
        sets_i = STUDY.cluster(clus_i).sets;
        main_cl_inds(clus_i) = 1;
        for j = 1:length(sets_i)
            comps_out(clus_i,sets_i(j)) = STUDY.cluster(clus_i).comps(j);
        end
    else
        outlier_cl_inds(clus_i) = 1;
    end
end
tmp = cellfun(@length,cellfun(@unique,{STUDY.cluster.sets},'UniformOutput',false),'UniformOutput',false);
tmp = cell2mat(tmp);
valid_clusters = find(tmp(3:end)>=0.5*(length(STUDY.subject)))+2;
%- remove all outlier cluster components
comps_out = comps_out(logical(main_cl_inds),:);
main_cl_inds = find(main_cl_inds);
outlier_cl_inds = find(outlier_cl_inds);
if all(comps_out==0)
    error('STUDY cluster information not generated');
end
end

