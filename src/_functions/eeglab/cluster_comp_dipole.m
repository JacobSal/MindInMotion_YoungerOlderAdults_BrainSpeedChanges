function [cluster_solution] = cluster_comp_dipole(ALLEEG, cluster_solution)
% Adapted from the bemobile pipeline -ish
% Compute centroid, rv from the given solution and STUDY
% compute for all clusters except parentclust and outlier clust
clsind = 1:length(cluster_solution); 
centroid = cell(length(clsind),1);
for clust = 1:length(clsind)
%     clear all_diplocs
%     clear residual_variances
    max_r = 0;
    len = length(cluster_solution(clsind(clust)).comps);
    tmppos = [ 0 0 0 ];
    tmpmom = [ 0 0 0 ];
    tmprv = 0;
    ndip = 0;
    all_diplocs = zeros(len,3);
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
                end
            end
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
    residual_variances = zeros(length(cluster_solution(clsind(clust)).comps),1);
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