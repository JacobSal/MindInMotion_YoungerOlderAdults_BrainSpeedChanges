function [STUDY,ALLEEG,comps_out,outliers] = cluster_pca_reduce(STUDY,ALLEEG,varargin)
%CLUSTER_PCA_REDUCE Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT
% NOTES:
%       U*S*V' = A where U,S,V are results from pca(), and A is the
%       input matrix;
%       ==> (I*I*V')' = (inv(U*S)*A)' ==> V = (inv(U*S)*A)'
%
%       (1) U*S can contain negative values which may change the
%       interpretation of the results
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
PCA_CUTOFF = 0.95; % units percent variance
DIST_CUTOFF = 30; % units mm
DO_MASKING = false;
%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct)
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETERS
%-
%% ===================================================================== %%
[STUDY, ~] = std_centroid(STUDY,ALLEEG,1:length(STUDY.cluster),'dipole');
%- extract centroid locations
dipfit_roi = [STUDY.cluster(1:end).centroid];
dipfit_roi = [dipfit_roi.dipole];
dipfit_roi = cat(1,dipfit_roi.posxyz);
%- loop through clusters
% pca_out = zeros(length(STUDY.cluster),length(ALLEEG),max([ALLEEG.pnts]));
comps_out = zeros(length(STUDY.cluster),length(ALLEEG));
outliers = cell(length(STUDY.cluster),length(ALLEEG));
for cluster_i = 2:length(STUDY.cluster)
    fprintf('==== Cluster %i ====\n',cluster_i);
    %- subset sets and comps
    sets_clust = unique(STUDY.cluster(cluster_i).sets);
    out = zeros(1,length(sets_clust));
%     clust_center = std_centroid(MAIN_STUDY,MAIN_ALLEEG,cluster_i,'dipole');
    for subj_i = sets_clust
        if isempty(ALLEEG(subj_i).icaact)
            error('Error. %s has an empty EEG.icaact matrix',ALLEEG(subj_i).subject);
        end
        idx = logical(STUDY.cluster(cluster_i).sets == subj_i);
        comps_clust = STUDY.cluster(cluster_i).comps(idx);
        % comps_clust = comps_clust(randperm(length(comps_clust))); % randomize order?
        %- just choosing a component based on the AMICA ICA sorting algorithm
        if ~isempty(comps_clust) && length(comps_clust) > 1
            fprintf('%s) Replacing components %s in cluster %i\n',ALLEEG(subj_i).subject,sprintf('%i,',comps_clust),cluster_i)
            %- change dipfit values 
            dipfit_out = ALLEEG(subj_i).dipfit.model(comps_clust);
            %- loop through fields
            FIELD_VALUE = 'posxyz';
            if isfield(dipfit_out,FIELD_VALUE)
                %- extract field values & calculate center
                tmp = vertcat(dipfit_out.(FIELD_VALUE));
                ctr = sum(tmp)/size(tmp,1); %dipfit_roi(cluster_i,:); %sum(tmp)/size(tmp,1);
                %## Calculate distances
                dist_to_ctr = zeros(size(tmp));
                dist_to_clust = zeros(size(tmp));
                %- get dist to center of cluster centroid
                for i = 1:length(dipfit_out)
                    dist_to_clust(i,:) = sqrt((tmp(i,:) - dipfit_roi(cluster_i,:)).^2);
                    fprintf('%7s%iIC''s distance to cluster centroid for dipfit.model.%s: %0.1f\n','',...
                            comps_clust(i),FIELD_VALUE,sum(dist_to_clust(i,:),2));
                end
                %- get dist to center of selected components
                for i = 1:length(dipfit_out)
                    dist_to_ctr(i,:) = sqrt((ctr - tmp(i,:)).^2);
                    fprintf('%7s%iIC''s distance to subject component''s center for dipfit.model.%s: %0.1f\n','',...
                            comps_clust(i),FIELD_VALUE,sum(dist_to_ctr(i,:),2));
                end
                if all(sum(dist_to_ctr,2) < DIST_CUTOFF)
                    %- pca
                    [U, S, V, pcas_sum, perc_var_ic, perc_var_pc] = eeglab_pca(ALLEEG(subj_i).icaact(comps_clust,1:end),PCA_CUTOFF);
                    %- scale dist to center
                    tmp_shift = perc_var_ic.*sum(dist_to_ctr,2);
%                     tmp_shift = sign_var.*perc_var_pc.*dist_to_ctr;
                    %- shift points
                    tmp_new = tmp+tmp_shift;
                    %- calculate new center
                    ctr_new = sum(tmp_new,1)/size(tmp_new,1);
                    dipfit_out(1).(FIELD_VALUE) = ctr_new;
                    fprintf('%7sdipfit.model.%s updated value for IC %i: (%0.1f,%0.1f,%0.1f)\n','',...
                                                FIELD_VALUE,comps_clust(1),ctr_new);
                    fprintf('%7snew dipfit.model.%s is %0.1f from centroid, and %0.1f from previous position\n','',...
                                                FIELD_VALUE,sqrt(sum((dipfit_roi(cluster_i,:)-ctr_new).^2)),...
                                                sqrt(sum((ALLEEG(subj_i).dipfit.model(comps_clust(1)).(FIELD_VALUE)-ctr_new).^2)));
                    %- assign dipole fits
%                     ALLEEG(subj_i).dipfit.model(comps_clust) = dipfit_out;
                    out(subj_i) = comps_clust(1);
                    comps_out(cluster_i,subj_i) = comps_clust(1);
                elseif any(sum(dist_to_ctr,2) < DIST_CUTOFF)
                    fprintf(2,'%7sWarning: Distances are big from centroid\n','');
                    sub_comps_clust = comps_clust(sum(dist_to_ctr,2) < DIST_CUTOFF);
                    [U, S, V, pcas_sum, perc_var_ic, perc_var_pc] = eeglab_pca(ALLEEG(subj_i).icaact(sub_comps_clust,1:end),PCA_CUTOFF);
                    %- change dipfit values 
                    dipfit_out = ALLEEG(subj_i).dipfit.model(sub_comps_clust);
                    %- extract field values & calculate center
                    tmp = vertcat(dipfit_out.(FIELD_VALUE));
                    ctr = sum(tmp,1)/size(tmp,1); %dipfit_roi(cluster_i,:); %sum(tmp)/size(tmp,1);
                    %## Calculate distances
                    dist_to_ctr = zeros(size(tmp));
                    dist_to_clust = zeros(size(tmp));
                    %- get dist to center of cluster centroid
                    for i = 1:length(dipfit_out)
                        dist_to_clust(i,:) = sqrt((tmp(i,:) - dipfit_roi(cluster_i,:)).^2);
                        fprintf('%7s%iIC''s distance to cluster centroid for dipfit.model.%s: %0.1f\n','',...
                                sub_comps_clust(i),FIELD_VALUE,sum(dist_to_clust(i,:),2));
                    end
                    %- get dist to center of selected components
                    for i = 1:length(dipfit_out)
                        dist_to_ctr(i,:) = sqrt((ctr - tmp(i,:)).^2);
                        fprintf('%7s%iIC''s distance to subject component''s center for dipfit.model.%s: %0.1f\n','',...
                                sub_comps_clust(i),FIELD_VALUE,sum(dist_to_ctr(i,:),2));
                    end
                    %- scale dist to center
                    tmp_shift = perc_var_ic.*sum(dist_to_ctr,2);
%                     tmp_shift = sign_var.*perc_var_pc.*dist_to_ctr;
                    %- shift points
                    tmp_new = tmp+tmp_shift;
                    %- calculate new center
                    ctr_new = sum(tmp_new,1)/size(tmp_new,1);
                    dipfit_out(1).(FIELD_VALUE) = ctr_new;
                    fprintf('%7sdipfit.model.%s updated value for IC %i: (%0.1f,%0.1f,%0.1f)\n','',...
                                                FIELD_VALUE,sub_comps_clust(1),ctr_new);
                    fprintf('%7snew dipfit.model.%s is %0.1f from centroid, and %0.1f from previous position\n','',...
                                                FIELD_VALUE,sqrt(sum((dipfit_roi(cluster_i,:)-ctr_new).^2)),...
                                                sqrt(sum((ALLEEG(subj_i).dipfit.model(sub_comps_clust(1)).(FIELD_VALUE)-ctr_new).^2)));
                    %- assign dipole fits
%                     ALLEEG(subj_i).dipfit.model(sub_comps_clust) = dipfit_out;
                    out(subj_i) = comps_clust(1);
                    comps_out(cluster_i,subj_i) = comps_clust(1);
                    comps_clust = sub_comps_clust;
                else
                    out(subj_i) = comps_clust(1);
                    comps_out(cluster_i,subj_i) = comps_clust(1);
                    outliers{cluster_i,subj_i} = comps_clust(2:end);
                    continue;
                end
            end
            rmv_idx = zeros(1,size(ALLEEG(subj_i).icaact,1));
            %- determine components to remove from ica
            for i = 1:length(comps_clust)
                if i ~= 1
                    rmv_idx(comps_clust(i)) = 1;
                end
            end
            %- replace desired component with PCA &assign principle component 
            % that explains the most variance of the independent components.
            % NOTE: all definitions of 'pc_out' are viable uses.
            pc_out = sign(diag(U)).*S*V'; % the idea here is to only use the principal component space 
            % pc_out = U*S*V';
            %- (1) sum pca's that belong to subject to create new
            %ica-activation
            ALLEEG(subj_i).icaact(comps_clust(1),:) = sum(pc_out(pcas_sum,:),1);
            %- (2) use the first pca 
            ALLEEG(subj_i).icaact(comps_clust(1),:) = pc_out(1,:); 
            %- remove other components within the cluster from ica.
            %{
            ori_icawinv = size(ALLEEG(subj_i).icawinv);
            ori_icasphere = size(ALLEEG(subj_i).icasphere);
            ori_icaact = size(ALLEEG(subj_i).icaact);
            ALLEEG(subj_i).icasphere = ALLEEG(subj_i).icasphere(~rmv_idx,:);
            ALLEEG(subj_i).icawinv = ALLEEG(subj_i).icawinv(:,~rmv_idx);
            ALLEEG(subj_i).icaact = ALLEEG(subj_i).icaact(~rmv_idx,:);
            fprintf('%8sicawinv size change: (%i,%i)\n','',ori_icawinv-size(ALLEEG(subj_i).icawinv));
            fprintf('%8sicasphere size change: (%i,%i)\n','',ori_icasphere-size(ALLEEG(subj_i).icasphere));
            fprintf('%8sicaact size change: (%i,%i)\n','',ori_icaact-size(ALLEEG(subj_i).icaact));
            %}
            %## Update ica 
            %- NOTE: ica_weights = EEG.icaact/EEG.data/EEG.icasphere
            %- NOTE: ica_winv = EEG.data/EEG.icaact
            ica_weights = (ALLEEG(subj_i).icaact/ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:))/ALLEEG(subj_i).icasphere;
            ica_winv = (ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:)/ALLEEG(subj_i).icaact);
%             tmp = sum(sqrt((ALLEEG(subj_i).icaweights(~rmv_idx,~rmv_idx)-ica_weights).^2),[1,2]);
            tmp = sum(sqrt((ALLEEG(subj_i).icaweights-ica_weights).^2),[1,2]);
            fprintf('%7ssum(sqrt((icaweights_new-icaweights_old).^2)) = %0.3f\n','',tmp);
            %## Update ALLEEG & comps_out
            ALLEEG(subj_i).icaweights = ica_weights;
            ALLEEG(subj_i).icawinv = ica_winv;
            %- Replace the first IC in the components to cluster with PCA
            comps_out(cluster_i,subj_i) = comps_clust(1);
             %- assign component for subject and cluster
            out(subj_i) = comps_clust(1);
            TMPS = STUDY.cluster(cluster_i).comps(idx);
            outliers{cluster_i,subj_i} = TMPS(TMPS ~= comps_clust(1));
        else
            fprintf('%s) Using component %i\n',ALLEEG(subj_i).subject,comps_clust);
            %- assign component for subject and cluster
            out(subj_i) = comps_clust;
            comps_out(cluster_i,subj_i) = comps_clust;
        end        
    end
    %- remove zeros
    out = out(out ~= 0);
    %- assign
    STUDY.cluster(cluster_i).sets = sets_clust;
    STUDY.cluster(cluster_i).comps = out;
    tmps = [];
    tmpc = [];
    for subj_i = 1:length(ALLEEG)
        tmps = [tmps, repmat(subj_i,1,length([outliers{cluster_i,subj_i}]))];
        tmpc = [tmpc, outliers{cluster_i,subj_i}];
    end
    STUDY.cluster(end+1).sets = tmps;
    STUDY.cluster(end).comps = tmpc;
    STUDY.cluster(end).name = sprintf('Outlier clust_%i',cluster_i);
    STUDY.cluster(end).parent = STUDY.cluster(cluster_i).parent;
    STUDY.cluster(end).algorithm = {'principal component analysis reduction of sifted brain components'};
end
fprintf('done.\n');
toc
end
%%  SUBFUNCTIONS
function [U, S, V, pcas_sum, perc_var_ic, perc_var_pc] = eeglab_pca(alleeg_icaact,pca_cutoff)
    %- use pca
    %- perform pca using randomization (nPerms = 200) w/ no rank
    %reduction. Note: may need to reduce rank to 1 here as it could
    %be assumed that the components clustered together have shared
    %information between them.
    [U,S,V] = pca(squeeze(alleeg_icaact),size(alleeg_icaact,1),200);
    %- use PCA coefficients to determine PC with greatest contributions
    perc_var_pc = diag(S)/sum(diag(S)); %test_S/total_var; %diag(S)/sum(diag(S))
    %- use the PCA coefficients to determine IC with greatest contribution
    perc_var_ic = diag(inv(U*S))/sum(diag(abs(inv(U*S))));
    %- loop
    tmpsum = 0;
    pcas_sum = [];
    for i = 1:length(perc_var_pc)
        fprintf('%7sPrincipal component %i explains %0.2g percent-variance\n','',i,perc_var_pc(i)*100);
        if tmpsum < pca_cutoff
            tmpsum = tmpsum + perc_var_pc(i);
            pcas_sum = [pcas_sum, i];
        end
    end
end
%% DEBUG PHRASES
%## validation plots
%{
figure;
hold on;
scatter3([tmp(:,1);ctr(1)],[tmp(:,2);ctr(1)],[tmp(:,3);ctr(1)],'DisplayName','original');
scatter3([tmp_new(:,1);ctr_new(1)],[tmp_new(:,2);ctr_new(1)],[tmp_new(:,3);ctr_new(1)],'DisplayName','shifted');
view([45,-45,45])
hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(sprintf('Subject %s',ALLEEG(subj_i).subject));
legend;
%}
%{
tt = 1000:5000;
ic_f = ALLEEG(subj_i).icaact(comps_clust,1:end);
% pc_f = sign(diag(U)).*S*V';
pc_f = U*S*V';
% pc_f = S*V';
pc_f = sum(pc_f(1:2,:),1);
figure;
hold on;
plot(ALLEEG(subj_i).times(tt),ic_f(1,tt),'DisplayName','ICA');                    
plot(ALLEEG(subj_i).times(tt),pc_f(1,tt),'DisplayName','PCA');
hold off;
xlabel('time');
ylabel('voltage (uV)');
legend();
%}
%% LOG
%- (01/13/2023), JS NOTE: need to update weights and sphere for ica
