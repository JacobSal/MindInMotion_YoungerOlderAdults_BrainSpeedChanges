function [cluster_struct, centroid] = evaluate_cluster(STUDY, ALLEEG, cluster_struct, evaluation_method,varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu    
%% Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@istruct);
addRequired(p,'cluster_struct',@ischar);
addRequired(p,'evaluation_method',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,cluster_struct,evaluation_method,varargin{:});
%## SET DEFAULTS
%% ===================================================================== %%
switch evaluation_method
    case 'min_rv'
        cluster_struct = cluster_struct;
        for k = 3:length(cluster_struct)
            clear subj_in_cluster min_rv comp_ind_store
            subj_in_cluster = unique(cluster_struct(k).sets);%Subjects in this cluster
            for j = 1:length(unique(cluster_struct(k).sets))    
                comp_ind = cluster_struct(k).comps(find(cluster_struct(k).sets == subj_in_cluster(j)));
                if ~isempty(comp_ind)
                    comp_ind_rv = [ALLEEG(subj_in_cluster(j)).dipfit.model(comp_ind).rv];
                    [~,min_rv_comp_ind] = min(comp_ind_rv);
                    comp_ind_store(1,j) = comp_ind(min_rv_comp_ind);
                    min_rv(1,j) = min(comp_ind_rv);
                end
            end
            % update cluster_update
            cluster_struct(k).sets = subj_in_cluster;
            cluster_struct(k).comps = comp_ind_store;
        end
        cluster_struct = compute_dipole(ALLEEG, cluster_struct);
    case 'min_ic'
       cluster_struct = cluster_struct;
        for k = 3:length(cluster_struct)
            clear subj_in_cluster min_ic comp_ind_store comp_ind_ic
            subj_in_cluster = unique(cluster_struct(k).sets);%Subjects in this cluster
            for j = 1:length(unique(cluster_struct(k).sets))    
                comp_ind = cluster_struct(k).comps(find(cluster_struct(k).sets == subj_in_cluster(j)));
                if ~isempty(comp_ind)
                    comp_ind_ic(1,j) = [comp_ind(1)];
                end
            end
            % update cluster_update
            cluster_struct(k).sets = subj_in_cluster;
            cluster_struct(k).comps = comp_ind_ic;
        end
        cluster_struct = compute_dipole(ALLEEG, cluster_struct);
    case 'pca_decomp'
        STUDY.cluster = cluster_struct;
        [STUDY,ALLEEG,comps_out,outliers] = cluster_pca_reduce(STUDY,ALLEEG);
    otherwise
        error('Choose a valid cluster evaluation option...\n');
end

end
%% ===================================================================== %%
function [cluster_solution] = compute_dipole(ALLEEG, cluster_solution)
% Adapted from the bemobile pipeline -ish
% Compute centroid, rv from the given solution and STUDY

% compute for all clusters except parentclust and outlier clust
clsind = 3:length(cluster_solution); 

for clust = 1:length(clsind)
    clear all_diplocs
    clear residual_variances
    max_r = 0;
    len = length(cluster_solution(clsind(clust)).comps);
    tmppos = [ 0 0 0 ];
    tmpmom = [ 0 0 0 ];
    tmprv = 0;
    ndip = 0;
    for k = 1:len 
        comp  = cluster_solution(clsind(clust)).comps(k);
        abset = cluster_solution(clsind(clust)).sets(1,k);
        if ~isfield(ALLEEG(abset), 'dipfit')
           warndlg2(['No dipole information available in dataset ' num2str(abset) ], 'Aborting compute centroid dipole');
           return;
        end
        if ~isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
            ndip   = ndip +1;
            posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
            momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
            if size(posxyz,1) == 2
                if all(posxyz(2,:) == [ 0 0 0 ])
                    posxyz(2,:) = [];
                    momxyz(2,:) = [];
                end;
            end;
            tmppos = tmppos + mean(posxyz,1);
            tmpmom = tmpmom + mean(momxyz,1);
            tmprv = tmprv + ALLEEG(abset).dipfit.model(comp).rv;
            if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
               if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                   load('-mat', ALLEEG(abset).dipfit.hdmfile);
                   max_r = max(max_r, max(vol.r));
               else % old version of dipfit
                   max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
               end
            end
            all_diplocs(k,:) = mean(posxyz,1);
        end
    end
    centroid{clust}.dipole.posxyz =  tmppos/ndip;
    centroid{clust}.dipole.momxyz =  tmpmom/ndip;
    centroid{clust}.dipole.rv =  tmprv/ndip;
    if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical') & (~isfield(ALLEEG(abset).dipfit, 'hdmfile')) %old dipfit
        centroid{clust}.dipole.maxr = max_r;
    end
    cluster_solution(clsind(clust)).dipole = centroid{clust}.dipole;
    cluster_solution(clsind(clust)).all_diplocs = all_diplocs;
    
    % compute spread (sum of squared deviation from centroid)
    squared_deviations = 0;
    for IC = 1:size(all_diplocs,1)

		% Pythagoras in 3D
        dist_this_IC = sqrt((all_diplocs(IC,1) - centroid{clust}.dipole.posxyz(1))^2 +...
							(all_diplocs(IC,2) - centroid{clust}.dipole.posxyz(2))^2 +...
							(all_diplocs(IC,3) - centroid{clust}.dipole.posxyz(3))^2);
		% square and add
        squared_deviations = squared_deviations + dist_this_IC^2;

    end
    
    cluster_solution(clsind(clust)).squared_deviations = squared_deviations;
    
    % store residual variances
    
    for ICind = 1:length(cluster_solution(clsind(clust)).comps)
        thisICdatasets = cluster_solution(clsind(clust)).sets(:,ICind);
        IC = cluster_solution(clsind(clust)).comps(ICind);           
        
        RV_temp = 0;
        for dataSet = 1:length(thisICdatasets)
            RV_temp(dataSet) = ALLEEG(thisICdatasets(dataSet)).dipfit.model(IC).rv;
        end
        
        assert(sum(diff(RV_temp))==0,'RVs of the same IC seem to be different between conditions!')
        residual_variances(ICind) = RV_temp(1);
        
    end
    median_rv = median(residual_variances);
    mean_rv = mean(residual_variances);
    
    cluster_solution(clsind(clust)).residual_variances = residual_variances;
    cluster_solution(clsind(clust)).median_rv = median_rv;
    cluster_solution(clsind(clust)).mean_rv = mean_rv;
    
end

end