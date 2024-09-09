function [STUDY,ALLEEG,new_cluster] = cluster_pca_reduce_simple(STUDY,ALLEEG,varargin)
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
%       (2) Make sure to resave EEG set files under a new name because this
%       code will edit the underlying ICA decomposition.
%           e.g., 
%            for subj_i = 1:length(TMP_ALLEEG)
%                EEG = eeg_checkset(TMP_ALLEEG(subj_i),'loaddata');
%                if isempty(EEG.icaact)
%                    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
%                    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
%                    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
%                end
%                TMP_ALLEEG(subj_i) = EEG;
%            end
%            [TMP_STUDY,TMP_ALLEEG,~,~] = cluster_pca_reduce(TMP_STUDY,TMP_ALLEEG);
%            for subj_i = 1:length(TMP_ALLEEG)
%                TMP_ALLEEG(subj_i).filename = sprintf('%s_pcareduced_comps.set',TMP_ALLEEG(subj_i).subject);
%                TMP_ALLEEG(subj_i) = pop_saveset(TMP_ALLEEG(subj_i),...
%                    'filename',TMP_ALLEEG(subj_i).filename,...
%                    'filepath',TMP_ALLEEG(subj_i).filepath);
%            end
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
PSD_RANDI_SAMPLES = 200;
DEFAULT_PCA_STRUCT = struct('pca_proc_type','continuous',...
    'variance_cutoff',80,... % percent variance
    'plot_psd',true,...
    'plot_topo',true,...
    'distance_cutoff',30); % mm
%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct)
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PCA_STRUCT',DEFAULT_PCA_STRUCT,@(x) validate_struct(x,DEFAULT_PCA_STRUCT));
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETERS
PCA_STRUCT = p.Results.PCA_STRUCT;
PCA_STRUCT = set_defaults_struct(PCA_STRUCT,DEFAULT_PCA_STRUCT);
%-
%% ===================================================================== %%
[STUDY, ~] = std_centroid(STUDY,ALLEEG,1:length(STUDY.cluster),'dipole');
%- extract centroid locations
dipfit_roi = [STUDY.cluster(1:end).centroid];
dipfit_roi = [dipfit_roi.dipole];
dipfit_roi = cat(1,dipfit_roi.posxyz);
%- loop through clusters
% pca_out = zeros(length(STUDY.cluster),length(ALLEEG),max([ALLEEG.pnts]));
tmp_cluster = STUDY.cluster;
comps_out = zeros(length(tmp_cluster),length(ALLEEG));
outliers = cell(length(tmp_cluster),length(ALLEEG));
clust_inds = find(~cellfun(@(x)contains(x,'parent','IgnoreCase',true) || contains(x,'outlier','IgnoreCase',true),{tmp_cluster.name}));
non_clust_inds = setdiff(2:length(tmp_cluster),clust_inds);
clust_inf_store = cell(length(clust_inds),1); % 1:subjs to pca, 2:inds to pca, 3:subjs not pca, 4:inds not pca
%##
for i = 1:length(clust_inds)
    cluster_i = clust_inds(i);
    fprintf('==== Cluster %i ====\n',cluster_i);
    %- subset sets and comps
    sets_clust = unique(tmp_cluster(cluster_i).sets);
    tmp_cell_store = cell(length(sets_clust),5); % 1:subjs to pca, 2:inds to pca, 3:subjs not pca, 4:inds not pca, 5: outliers
%     clust_center = std_centroid(MAIN_STUDY,MAIN_ALLEEG,cluster_i,'dipole');
    for k = 1:length(sets_clust)
        subj_i = sets_clust(k);
        idx = logical(tmp_cluster(cluster_i).sets == subj_i);
        comps_clust = tmp_cluster(cluster_i).comps(idx);
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
                for j = 1:length(dipfit_out)
                    dist_to_clust(j,:) = sqrt((tmp(j,:) - dipfit_roi(cluster_i,:)).^2);
                    fprintf('%7s%iIC''s distance to cluster centroid for dipfit.model.%s: %0.1f\n','',...
                            comps_clust(j),FIELD_VALUE,sum(dist_to_clust(j,:),2));
                end
                %- get dist to center of selected components
                for j = 1:length(dipfit_out)
                    dist_to_ctr(j,:) = sqrt((ctr - tmp(j,:)).^2);
                    fprintf('%7s%iIC''s distance to subject component''s center for dipfit.model.%s: %0.1f\n','',...
                            comps_clust(j),FIELD_VALUE,sum(dist_to_ctr(j,:),2));
                end
                if all(sum(dist_to_ctr,2) < PCA_STRUCT.distance_cutoff)
                    tmp_cell_store{k,2} = comps_clust;
                    tmp_cell_store{k,1} = subj_i;
                    
                elseif any(sum(dist_to_ctr,2) < PCA_STRUCT.distance_cutoff)
                    fprintf(2,'%7sWarning: Distances are big from centroid\n','');
                    sub_comps_clust = comps_clust(sum(dist_to_ctr,2) < PCA_STRUCT.distance_cutoff);
                    if length(sub_comps_clust) > 1 && length(sub_comps_clust)==length(comps_clust)
                        tmp_cell_store{k,2} = sub_comps_clust;
                        tmp_cell_store{k,1} = subj_i;
                    elseif length(sub_comps_clust) > 1 && length(sub_comps_clust)~=length(comps_clust)
                        tmp_cell_store{k,2} = sub_comps_clust;
                        tmp_cell_store{k,1} = subj_i;
                        tmp_cell_store{k,5} = comps_clust(sum(dist_to_ctr,2) > PCA_STRUCT.distance_cutoff);
                    else
                        tmp_cell_store{k,4} = sub_comps_clust;
                        tmp_cell_store{k,3} = subj_i;
                    end
                else
                    fprintf('%s) Using component %i\n',ALLEEG(subj_i).subject,comps_clust(1));
                    tmp_cell_store{k,4} = comps_clust(1);
                    tmp_cell_store{k,3} = subj_i;
                end
            end
        else
            fprintf('%s) Using component %i\n',ALLEEG(subj_i).subject,comps_clust);
            %- assign component for subject and cluster
            tmp_cell_store{k,4} = comps_clust;
            tmp_cell_store{k,3} = subj_i;
        end
    end
    clust_inf_store{i} = tmp_cell_store;
end

%%
new_cluster = struct('name',tmp_cluster(1).name,...
    'sets',tmp_cluster(1).sets,...
    'comps',tmp_cluster(1).comps,...
    'parent',tmp_cluster(1).parent,...
    'child',{tmp_cluster(1).child},...
    'preclust',tmp_cluster(1).preclust,...
    'algorithm','custom pca algorithm');
tmp_clust_store = cell(length(tmp_cluster),2);
for subj_i = 1:length(STUDY.datasetinfo)
    subj_char = STUDY.datasetinfo(subj_i).subject;
    fprintf('Processing Subject %s...\n',subj_char);
    inds = zeros(length(clust_inf_store),1);
    for ii = 1:length(clust_inf_store)
        tmp = clust_inf_store{ii}(:,1);
        inds(ii) = any(cellfun(@(x) ~isempty(x) && x==subj_i,tmp));
    end
    inds_inf = find(inds);
    inds_clust = clust_inds(inds_inf);
    for i = 1:length(inds_inf)
        cluster_i = inds_clust(i);
        tmp = clust_inf_store{inds_inf(i)}(:,1);
        ind = cellfun(@(x) ~isempty(x) && x==subj_i,tmp);
        %-
        if isempty(ALLEEG(subj_i).icaact)
            error('Error. %s has an empty EEG.icaact matrix',ALLEEG(subj_i).subject);
        end
        % idx = logical(STUDY.cluster(cluster_i).sets == subj_i);
        % comps_clust = STUDY.cluster(cluster_i).comps(idx);
        comps_inf = sort(clust_inf_store{inds_inf(i)}{ind,2});
        %## MAIN PCA PROCESSING
        switch PCA_STRUCT.pca_proc_type
            case 'continuous'
                %## PCA ALGORITHM
                CC = cov(squeeze(ALLEEG(subj_i).icaact(comps_inf,1:end))');
                %- obtain eigenvalues
                [v, d] = eig(CC);
                %- find component with greatest eigenvalue
                diag_d = diag(d);
                [~, sort_ix] = sort(diag_d,'descend');
                V = v(:,sort_ix);
                Ai = V.';
                %- variance estimates for component extraction
                var_sum = 0;
                comps_keep = 0;
                for j = 1:length(sort_ix)
                    fprintf('%7sComponent %i accounts %0.2f%% variance\n','',comps_inf(sort_ix(j)),(diag_d(sort_ix(j))/sum(diag_d))*100);
                    if var_sum < PCA_STRUCT.variance_cutoff
                        var_sum = var_sum + (diag_d(sort_ix(j))/sum(diag_d))*100;
                        comps_keep = comps_keep + 1;
                    end
                end
                %## PCA ALGORITHM
                pca_rmv = [];
                if length(sort_ix) > comps_keep
                    pca_rmv = comps_inf(sort_ix(comps_keep+1:end));
                end
                comps_keep = sort(comps_inf(sort_ix(1:comps_keep)));
                CC = cov(squeeze(ALLEEG(subj_i).icaact(comps_keep,1:end))');
                %- obtain eigenvalues
                [v, d] = eig(CC);
                %- find component with greatest eigenvalue
                diag_d = diag(d);
                [~, sort_ix] = sort(diag_d,'descend');
                V = v(:,sort_ix);
                Ai = V.';
                for j = 1:length(sort_ix)
                    fprintf('%7sComponent %i accounts %0.2f%% variance\n','',comps_keep(sort_ix(j)),(diag_d(sort_ix(j))/sum(diag_d))*100);
                end
                % pnts = size(ALLEEG(subj_i).icaact(:,:),2);
                % r = xcorr(ALLEEG(subj_i).icaact(comps_inf(sort_ix(1)),1:end),ALLEEG(subj_i).icaact(comps_inf(sort_ix(2)),1:end),pnts-1,'normalized');
                % mean(r)
                % v = conv(ALLEEG(subj_i).icaact(comps_inf(sort_ix(1)),1:end),ALLEEG(subj_i).icaact(comps_inf(sort_ix(2)),1:end));
                %- cross correlation
                % x = ALLEEG(subj_i).icaact(comps_inf(sort_ix(1)),1:end)';
                % y = ALLEEG(subj_i).icaact(comps_inf(sort_ix(2)),1:end)';
                % nx = numel(x);
                % ny = numel(y);
                % m = max(nx,ny);
                % m2 = findTransformLength(m);
                % mxl = m-1;
                % X = fft(x,m2,1);
                % Y = fft(y,m2,1);
                % if isreal(x) && isreal(y)
                %     c1 = ifft(X.*conj(Y),[],1,'symmetric');
                % else
                %     c1 = ifft(X.*conj(Y),[],1);
                % end
                % % Keep only the lags we want and move negative lags before positive
                % % lags.
                % c = [c1(m2 - mxl + (1:mxl)); c1(1:mxl+1)];
                % c = c./m;
                %- convert components to PC space
                % ica_psc = squeeze(ALLEEG(subj_i).icaact(comps_inf(sort_ix),1:end))';
                ica_psc = squeeze(ALLEEG(subj_i).icaact(comps_keep,1:end))';
                t_psc = ica_psc(:,:) * V;
                ica_act = (t_psc(:,1) * Ai(1,:))';
                %- spectrums of pcs
                if PCA_STRUCT.plot_psd
                    sub_samp = randi(size(ica_act,2)-ALLEEG(subj_i).srate*60*6,1,1);
                    fprintf('%7sComputing sample PSD...\n','');
                    frames = (ALLEEG(subj_i).srate);
                    [psd_out, f] = spectopo(ica_act(:,sub_samp:sub_samp+(ALLEEG(subj_i).srate*60*6)-1), frames, ALLEEG(1).srate,...
                        'plot','off',...
                        'boundaries',[],...
                        'chanlocs',ALLEEG(subj_i).chanlocs,...
                        'freqfac',ceil(ALLEEG(subj_i).srate/(pi^2)),...
                        'winsize',ceil((ALLEEG(subj_i).srate)/pi),...
                        'overlap',[],...
                        'verbose','off',...
                        'wintype','hamming');
                    [psd_out_old, f] = spectopo(ALLEEG(subj_i).icaact(comps_keep(sort_ix),sub_samp:sub_samp+(ALLEEG(subj_i).srate*60*6)-1), frames, ALLEEG(1).srate,...
                        'plot','off',...
                        'boundaries',[],...
                        'chanlocs',ALLEEG(subj_i).chanlocs,...
                        'freqfac',ceil(ALLEEG(subj_i).srate/(pi^2)),...
                        'winsize',ceil((ALLEEG(subj_i).srate)/pi),...
                        'overlap',[],...
                        'verbose','off',...
                        'wintype','hamming');
                    % ALLEEG(subj_i).icaact(comps_inf,1:end)
                    %## PLOT
                    freq_inds = find(f >= 3 & f <= 100);
                    fig = figure;
                    hold on;
                    for tmp_i = 1:length(comps_keep)
                        plot(f(freq_inds),psd_out(tmp_i,freq_inds),'DisplayName',sprintf('comp %i',comps_keep(sort_ix(tmp_i))));
                        plot(f(freq_inds),psd_out_old(tmp_i,freq_inds),'DisplayName',sprintf('old comp %i',comps_keep(sort_ix(tmp_i))),'LineStyle','-.');
                    end
                    hold off;
                    legend();
                    ylim([min(psd_out,[],'all'),max(psd_out,[],'all')])
                    exportgraphics(fig,[ALLEEG(subj_i).filepath filesep sprintf('%s_cl%i_psd_valid.jpg',ALLEEG(subj_i).subject,cluster_i)])
                end
                %## ASSIGN NEW ACTIVATIONS
                if length(size(ALLEEG(subj_i).icaact)) == 3
                    epoch_act = zeros(length(comps_keep),ALLEEG(subj_i).pnts,ALLEEG(subj_i).trials);
                    cnt = 1;
                    for trial_i = 1:ALLEEG(subj_i).trials
                        epoch_act(:,:,trial_i) = ica_act(:,cnt:cnt+ALLEEG(subj_i).pnts-1);
                        cnt = cnt + ALLEEG(subj_i).pnts;
                    end
                    ALLEEG(subj_i).icaact(comps_keep(sort_ix),:,:) = epoch_act;
                elseif length(size(ALLEEG(subj_i).icaact)) == 2
                    ALLEEG(subj_i).icaact(comps_keep(sort_ix),:) = ica_act;
                else
                    error('Unable to determine the dimensions of EEG.icaact\n');
                end
            case 'epoch'
                %## PCA ALGORITHM
                CC = cov(squeeze(ALLEEG(subj_i).icaact(comps_inf,1:end))');
                %- obtain eigenvalues
                [v, d] = eig(CC);
                %- find component with greatest eigenvalue
                diag_d = diag(d);
                [~, sort_ix] = sort(diag_d,'descend');
                V = v(:,sort_ix);
                Ai = V.';
                %- variance estimates for component extraction
                var_sum = 0;
                comps_keep = 0;
                for j = 1:length(sort_ix)
                    fprintf('%7sComponent %i accounts %0.1f%% variance\n','',comps_inf(sort_ix(j)),(diag_d(sort_ix(j))/sum(diag_d))*100);
                    if var_sum < PCA_STRUCT.variance_cutoff
                        var_sum = var_sum + (diag_d(sort_ix(j))/sum(diag_d))*100;
                        comps_keep = comps_keep + 1;
                    end
                end
                %## PCA ALGORITHM
                pca_rmv = [];
                 if length(sort_ix) > comps_keep
                    pca_rmv = comps_inf(sort_ix(comps_keep+1:end));
                end
                comps_keep = sort(comps_inf(sort_ix(1:comps_keep)));
                CC = cov(squeeze(ALLEEG(subj_i).icaact(comps_keep,1:end))');
                %- obtain eigenvalues
                [v, d] = eig(CC);
                %- find component with greatest eigenvalue
                diag_d = diag(d);
                [~, sort_ix] = sort(diag_d,'descend');
                V = v(:,sort_ix);
                Ai = V.';
                for j = 1:length(sort_ix)
                    fprintf('%7sComponent %i accounts %0.1f%% variance\n','',comps_keep(sort_ix(j)),(diag_d(sort_ix(j))/sum(diag_d))*100);
                end
                % ica_psc = permute(ALLEEG(subj_i).icaact(comps_inf(sort_ix),1:end,:),[3,2,1]);
                ica_psc = permute(ALLEEG(subj_i).icaact(comps_keep,1:end,:),[3,2,1]);
                ica_out = zeros(size(ica_psc));
                for trial_i = 1:size(ica_psc,1)
                    t_psc = squeeze(ica_psc(trial_i,:,:)) * V;
                    ica_out(trial_i,:,:) = t_psc(:,1) * Ai(1,:);
                end
                ica_act = permute(ica_out,[3,2,1]);
                %- spectrums of pcs
                if PCA_STRUCT.plot_psd
                    sub_samp = randi(size(ica_act,3),PSD_RANDI_SAMPLES,1);
                    psd_out = cell(size(ica_act,3),1);
                    psd_out_old = cell(size(ica_act,3),1);
                    fprintf('%7sComputing sample PSD...\n','');
                    for trial_i = 1:length(sub_samp)
                        [tmp, f] = spectopo(ica_act(:,:,sub_samp(trial_i)), size(ica_act,2), ALLEEG(1).srate,...
                            'plot','off',...
                            'boundaries',[],...
                            'chanlocs',ALLEEG(subj_i).chanlocs,...
                            'freqfac',ceil(ALLEEG(subj_i).srate/(pi^2)),...
                            'winsize',ceil((ALLEEG(subj_i).srate)/pi),...
                            'overlap',[],...
                            'verbose','off',...
                            'wintype','hamming');
                        [tmp_old, f] = spectopo(ALLEEG(subj_i).icaact(comps_keep(sort_ix),:,trial_i), size(ALLEEG(subj_i).icaact,2), ALLEEG(1).srate,...
                            'plot','off',...
                            'boundaries',[],...
                            'chanlocs',ALLEEG(subj_i).chanlocs,...
                            'freqfac',ceil(ALLEEG(subj_i).srate/(pi^2)),...
                            'winsize',ceil((ALLEEG(subj_i).srate)/pi),...
                            'overlap',[],...
                            'verbose','off',...
                            'wintype','hamming');
                        psd_out{trial_i} = tmp;
                        psd_out_old{trial_i} = tmp_old;
                        if trial_i == 1
                            fprintf('%7s','');
                        end
                        if mod(trial_i,10) == 0, fprintf('%d ', trial_i); end
                    end
                    fprintf('\n');
                    psd_out = cat(3,psd_out{:});
                    psd_out = mean(psd_out,3);
                    psd_out_old = cat(3,psd_out_old{:});
                    psd_out_old = mean(psd_out_old,3);
                    %## PLOT
                    freq_inds = find(f >= 3 & f <= 100);
                    fig = figure;
                    hold on;
                    for tmp_i = 1:length(comps_keep)
                        plot(f(freq_inds),psd_out(tmp_i,freq_inds),'DisplayName',sprintf('component %i',comps_keep(sort_ix(tmp_i))));
                        plot(f(freq_inds),psd_out_old(tmp_i,freq_inds),'DisplayName',sprintf('old comp %i',comps_keep(sort_ix(tmp_i))),'LineStyle','-.');
                    
                    end
                    hold off;
                    legend();
                    ylim([min(psd_out,[],'all'),max(psd_out,[],'all')])
                    exportgraphics(fig,[ALLEEG(subj_i).filepath filesep sprintf('%s_cl%i_psd_valid.jpg',ALLEEG(subj_i).subject,cluster_i)])
                end
                %## ASSIGN NEW ACTIVATIONS
                ALLEEG(subj_i).icaact(comps_keep,:,:) = ica_act;
            otherwise
                error('A valid PCA_STRUCT.pca_proc_type was not specified. Choose [continuous/epoch]');
        end
        %## TOPOPLOTS (BEFORE PCA DIMENSIONALITY REDUCTION)
        if PCA_STRUCT.plot_topo
            chanlocs = ALLEEG(subj_i).chanlocs(ALLEEG(subj_i).icachansind);
            topo_struct = struct('grid',[],...
                'x',[],...
                'y',[]);
            for tmp_i = 1:length(comps_keep)
                chan_i = comps_keep(tmp_i);
                [~, grid, ~, Xi, Yi] = topoplot(ALLEEG(subj_i).icawinv(:,chan_i), chanlocs, ...
                                                                  'verbose', 'off',...
                                                                   'electrodes', 'on' ,'style','both',...
                                                                   'plotrad',0.55,'intrad',0.55,...
                                                                   'noplot', 'on', 'chaninfo', ALLEEG(subj_i).chaninfo);
                topo_struct(tmp_i).grid = grid;
                topo_struct(tmp_i).x = Xi;
                topo_struct(tmp_i).y = Yi;
            end
            %## PLOT
            AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
            fig = figure;
            set(gca,AXES_DEFAULT_PROPS{:})
            hold on;
            for tmp_i = 1:length(comps_keep)
                dip_i = comps_keep(sort_ix(tmp_i));
                if ALLEEG(subj_i).dipfit.model(dip_i).rv < 0.15
                    figure(fig);
                    sbplot(ceil(length(comps_keep)/5),5,tmp_i)
                    Xi = topo_struct(tmp_i).x;
                    Yi = topo_struct(tmp_i).y;
                    toporeplot(topo_struct(tmp_i).grid, 'style', 'both',...
                        'plotrad',0.5,'intrad',0.5, 'verbose', 'off',...
                        'xsurface', Xi, 'ysurface', Yi );
                    title(sprintf('(RV: %0.3f)\nComponent %i',ALLEEG(subj_i).dipfit.model(dip_i).rv,dip_i));
                    colormap(linspecer); 
                end
            end
            hold off;
            exportgraphics(fig,[ALLEEG(subj_i).filepath filesep sprintf('%s_cl%i_topo_valid.jpg',ALLEEG(subj_i).subject,cluster_i)])
        end
        
        %## CLUSTER ASSIGNMENTS
        if length(new_cluster) < cluster_i
            new_cluster(cluster_i).sets = [];
            new_cluster(cluster_i).comps = [];
        end
        new_cluster(cluster_i).sets = [new_cluster(cluster_i).sets, subj_i];
        new_cluster(cluster_i).comps = [new_cluster(cluster_i).sets, comps_keep(sort_ix(1))];
        tmp_clust_store{cluster_i,1} = [tmp_clust_store{cluster_i,1}, repmat(subj_i,1,length(sort_ix(2:end)))];
        tmp_clust_store{cluster_i,2} = [tmp_clust_store{cluster_i,2}, comps_keep(sort_ix(2:end))];
        outliers = clust_inf_store{inds_inf(i)}{ind,5};
        if ~isempty(outliers)
            tmp_clust_store{cluster_i,1} = [tmp_clust_store{cluster_i,1}, repmat(subj_i,1,length(outliers))];
            tmp_clust_store{cluster_i,2} = [tmp_clust_store{cluster_i,2}, outliers];
        end
        if ~isempty(pca_rmv)
            tmp_clust_store{cluster_i,1} = [tmp_clust_store{cluster_i,1}, repmat(subj_i,1,length(pca_rmv))];
            tmp_clust_store{cluster_i,2} = [tmp_clust_store{cluster_i,2}, pca_rmv];
        end
    end
end
%% RECALCULATE ICA INFORMATION
%- NOTE: ica_weights = EEG.icaact/EEG.data/EEG.icasphere
%- NOTE: ica_winv = EEG.data/EEG.icaact
for subj_i = 1:length(ALLEEG)
    ica_weights = (ALLEEG(subj_i).icaact(:,:)/ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:))/ALLEEG(subj_i).icasphere;
    ica_winv = (ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:)/ALLEEG(subj_i).icaact(:,:));
    tmp = sum(sqrt((ALLEEG(subj_i).icaweights-ica_weights).^2),[1,2]);
    fprintf('%7ssum(sqrt((icaweights_new-icaweights_old).^2)) = %0.3f\n','',tmp);
    %## Update ALLEEG & comps_out
    ALLEEG(subj_i).icaweights = ica_weights;
    ALLEEG(subj_i).icawinv = ica_winv;
end
%## ASSIGN NEW ACTIVATIONS
% if length(size(ALLEEG(subj_i).icaact)) == 3
%     epoch_act = zeros(length(comps_inf),ALLEEG(subj_i).pnts,ALLEEG(subj_i).trials);
%     cnt = 1;
%     for trial_i = 1:ALLEEG(subj_i).trials
%         epoch_act(:,:,trial_i) = ica_act(:,cnt:cnt+ALLEEG(subj_i).pnts);
%         cnt = cnt + ALLEEG(subj_i).pnts+1;
%     end
%     ALLEEG(subj_i).icaact(comps_inf(sort_ix),:,:) = epoch_act;
% elseif length(size(ALLEEG(subj_i).icaact)) == 2
%     ALLEEG(subj_i).icaact(comps_inf(sort_ix),:) = ica_act;
% else
%     error('Unable to determine the dimensions of EEG.icaact\n');
% end
 %## TOPOPLOTS (BEFORE PCA DIMENSIONALITY REDUCTION)
if PCA_STRUCT.plot_topo
    for subj_i = 1:length(ALLEEG)
        inds = zeros(length(clust_inf_store),1);
        for ii = 1:length(clust_inf_store)
            tmp = clust_inf_store{ii}(:,1);
            inds(ii) = any(cellfun(@(x) ~isempty(x) && x==subj_i,tmp));
        end
        inds_inf = find(inds);
        inds_clust = clust_inds(inds_inf);
        for i = 1:length(inds_inf)
            cluster_i = inds_clust(i);
            tmp = clust_inf_store{inds_inf(i)}(:,1);
            ind = cellfun(@(x) ~isempty(x) && x==subj_i,tmp);
            comps_inf = clust_inf_store{inds_inf(i)}{ind,2};
            chanlocs = ALLEEG(subj_i).chanlocs(ALLEEG(subj_i).icachansind);
            topo_struct = struct('grid',[],...
                'x',[],...
                'y',[]);
            for tmp_i = 1:length(comps_inf)
                try
                    chan_i = comps_inf(tmp_i);
                    [~, grid, ~, Xi, Yi] = topoplot(ALLEEG(subj_i).icawinv(:,chan_i), chanlocs, ...
                                                                      'verbose', 'off',...
                                                                       'electrodes', 'on' ,'style','both',...
                                                                       'plotrad',0.55,'intrad',0.55,...
                                                                       'noplot', 'on', 'chaninfo', ALLEEG(subj_i).chaninfo);
                    topo_struct(tmp_i).grid = grid;
                    topo_struct(tmp_i).x = Xi;
                    topo_struct(tmp_i).y = Yi;
                catch e
                    fprintf('Error on subject %s component %i\n\n%s\n',ALLEEG(subj_i).subject,chan_i,getReport(e));
                end
            end
            %## PLOT
            AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
            fig = figure;
            set(gca,AXES_DEFAULT_PROPS{:})
            hold on;
            for tmp_i = 1:length(comps_inf)
                if ALLEEG(subj_i).dipfit.model(comps_inf(tmp_i)).rv < 0.15
                    figure(fig);
                    sbplot(ceil(length(comps_inf)/5),5,tmp_i)
                    Xi = topo_struct(tmp_i).x;
                    Yi = topo_struct(tmp_i).y;
                    toporeplot(topo_struct(tmp_i).grid, 'style', 'both',...
                        'plotrad',0.5,'intrad',0.5, 'verbose', 'off',...
                        'xsurface', Xi, 'ysurface', Yi );
                    title(sprintf('(RV: %0.3f)\nComponent %i',ALLEEG(subj_i).dipfit.model(comps_inf(tmp_i)).rv,comps_inf(tmp_i)));
                    colormap(linspecer); 
                end
            end
            hold off;
            exportgraphics(fig,[ALLEEG(subj_i).filepath filesep sprintf('%s_cl%i_topo_valid_new.jpg',ALLEEG(subj_i).subject,cluster_i)])
        end
    end
end
%## OUTLIER CLUSTER (MIM, CHANG)
for ii = 1:length(non_clust_inds)
    ff = fields(new_cluster);
    for f = 1:length(fields(new_cluster))
        new_cluster(non_clust_inds(ii)).(ff{f}) = tmp_cluster(non_clust_inds(ii)).(ff{f});
    end
end
%## FINAL ASSIGNMENT OF OUTLIER CLUSTERS
for ii = 1:length(tmp_clust_store)
    cluster_i = clust_inds(ii);
    new_cluster(cluster_i).name = sprintf('clust_%i',cluster_i);
    new_cluster(cluster_i).parent = tmp_cluster(cluster_i).parent;
    new_cluster(cluster_i).algorithm = [tmp_cluster(cluster_i).algorithm, sprintf('cluster_pca_reduce')];
    %-
    new_cluster(end+1).sets = tmp_clust_store{cluster_i,1};
    new_cluster(end).comps = tmp_clust_store{cluster_i,2};
    new_cluster(end).name = sprintf('Outlier clust_%i',cluster_i);
    new_cluster(end).parent = new_cluster(cluster_i).parent;
    new_cluster(end).algorithm = [tmp_cluster(cluster_i).algorithm, sprintf('cluster_pca_reduce')];
end
STUDY.cluster = new_cluster;
end
%%
function m = findTransformLength(m)
    m = 2*m;
    while true
        r = m;
        for p = [2 3 5 7]
            while (r > 1) && (mod(r, p) == 0)
                r = r / p;
            end
        end
        if r == 1
            break;
        end
        m = m + 1;
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
%## COVARIANCE METHOD PCA VALIDATION
%{
j = 1;
%-
tmp = ica_act(j,:);
%-
tmp = (ica_psc(:,j) * Ai(j,:))';
tmp = tmp(1,:);
%-
tmp = sum(ica_psc(:,j) * Ai(j,:),2)';
%## DECOMPOSE SIGNAL
for tmp_j = 1:length(sort_ix)
    for tmp_k = 1:length(sort_ix)
        pc1 = (ica_psc(:,tmp_j) * Ai(tmp_k,:))'; pc1 = pc1(1,:);
        pcr = (ica_psc(:,tmp_j) * Ai(tmp_k,:))'; pcr = sum(pcr(2:end,:));
        figure;
        hold on;
        plot(pc1(1,1:500),'color','r');
        plot(pcr(1,1:500),'color','b');
        for tmp_i = 1:length(sort_ix)
            plot(ALLEEG(subj_i).icaact(comps_inf(sort_ix(tmp_i)),1:500),'color',[0,1,0,0.7])
        end
        title(sprintf('tmp_j: %i, tmp_k: %i',tmp_j,tmp_k));
        hold off;
    end
end
%## RECREATE ORIGINAL SIGNAL
for tmp_k = 1:length(sort_ix)
    pc1 = (ica_psc(:,:) * Ai(:,tmp_k))';
    figure;
    hold on;
    plot(pc1(1,1:500),'color','r');
    for tmp_i = 1:length(sort_ix)
        plot(ALLEEG(subj_i).icaact(comps_inf(sort_ix(tmp_i)),1:500),'color',[0,1,0,0.7])
    end
    title(sprintf('tmp_j: %i, tmp_k: %i',tmp_j,tmp_k));
    hold off;
end
%}
%% LOG
%- (01/13/2023), JS NOTE: need to update weights and sphere for ica
