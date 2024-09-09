function [cluster_solutions,fig] = mim_cluster_process(STUDY,ALLEEG,varargin)
%MIM_CLUSTER_PROCESS Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are CLUSTER_STRUCT.saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% code adapted from bemobile pipeline & EEGLAB
%## TIME
tic;
%## DEFINE DEFAULTS
fig = [];
%- clustering struct
DEF_CLUSTER_STRUCT = struct('algorithm','kmeans',...
    'clust_k_num',12,... %(10-12 seems ideal, see. Makoto's PCA test)
    'clust_k_evals',(5:25),...
    'clust_k_maxiter',10000,...
    'clust_k_replicates',1000,...
    'clust_k_repeat_iters',50,...
    'clust_k_repeat_std',3,... % (inf | INT), INT uses robust_kmeans_CL
    'clust_k_robust_maxiter',5,...
    'clust_k_empty_action','drop',...
    'do_eval_clusters',true);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'CLUSTER_STRUCT',DEF_CLUSTER_STRUCT,@(x) validate_struct(x,DEF_CLUSTER_STRUCT));
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- PARAMETER
CLUSTER_STRUCT = p.Results.CLUSTER_STRUCT;
CLUSTER_STRUCT = set_defaults_struct(CLUSTER_STRUCT,DEF_CLUSTER_STRUCT);
%% ===================================================================== %%
if CLUSTER_STRUCT.do_eval_clusters
    cnt = 1;
    cluster_idx = zeros(length([STUDY.datasetinfo.comps]),length(CLUSTER_STRUCT.clust_k_evals));
    %## cluster using std_preclust weightings for each K desired (deafault here is kmeans with no outlier rejection).
    for i = 1:length(CLUSTER_STRUCT.clust_k_evals)
        %- temp variables
        params = CLUSTER_STRUCT;
        params.clust_k_num = CLUSTER_STRUCT.clust_k_evals(i);
        %- calculate clusters
        [~, clust_out] = cluster_dips(STUDY,ALLEEG,params);
        cluster_idx(1:length(clust_out.IDX),cnt) = clust_out.IDX;
        cnt = cnt+1;
    end
    % cluster_idx = cluster_idx(:,all(cluster_idx ~= 0));
    %##
    fig = plot_k_evals(STUDY,cluster_idx,STUDY.filepath);
end
%##
% (08/20/2023) JS, to force kmeans to kmeans_robust you must set an
% outliers std
CLUSTER_STRUCT.clust_k_repeat_std = CLUSTER_STRUCT.clust_k_repeat_std;
fprintf('Clustering %d times...\n', CLUSTER_STRUCT.clust_k_repeat_iters)
cluster_solutions = cell(length(CLUSTER_STRUCT.clust_k_num),1);
% parfor i = 1:length(CLUSTER_STRUCT.clust_k_evals)
for i = 1:length(CLUSTER_STRUCT.clust_k_num)
    cl_i = CLUSTER_STRUCT.clust_k_num(i);
    fprintf('Running clustering for k=%i',cl_i);
    %- temp variables
    params = CLUSTER_STRUCT;
    params.clust_k_num = cl_i;
    solutions = cell(CLUSTER_STRUCT.clust_k_repeat_iters,1);
    %- iterate
    for iter = 1:CLUSTER_STRUCT.clust_k_repeat_iters
        % start time
        tic
        fprintf('Iteration %d of %d \n',iter,CLUSTER_STRUCT.clust_k_repeat_iters)
        % clustering according to parameters - Note by Chang Liu, basically
        % doing the robust_kmeans cluster in pop_clust.m, I prefer not
        % overwrite the original STUDY
        [TMP_STUDY,~] = cluster_dips(STUDY,ALLEEG,params); 
        % store info in output struct
        solutions{iter} = TMP_STUDY.cluster;
        % stop checking the time and plot estimated time of arrival
        lastduration = toc; 
        eta = lastduration*(CLUSTER_STRUCT.clust_k_repeat_iters-iter);
        fprintf('Last duration: %1.2f s \n',round(lastduration,2))
        fprintf('ETA: %d h, %d min \n', floor(eta/3600), round(((eta/3600)-floor(eta/3600))*60))
    end
    %- store clustering info
    param_struct = struct('preclustparams',TMP_STUDY.etc.preclust,...
        'outlier_sigma',CLUSTER_STRUCT.clust_k_repeat_std,...
        'cluster_params',CLUSTER_STRUCT);
    tmp_clust_sol = struct('parameters',param_struct,...
        'solutions',{solutions});
    cluster_solutions{i} = tmp_clust_sol;
end
end
%% (SUBFUNCTION) ======================================================= %%
function [STUDY,cluster_outcome] = cluster_dips(STUDY,ALLEEG,params_in)
    %MIM_CLUSTERING_ANL Summary of this function goes here
    % CAT CODE
    %  _._     _,-'""`-._
    % (,-.`._,'(       |\`-/|
    %     `-.-' \ )-`( , o o)
    %           `-    \`_`"'-
    % Code Designers: Chang Liu, Jacob Salminen
    % Code Date: 07/17/2023, MATLAB 2020b
    % Written by Chang - 2023-4-15 to run plotERSP with parallel processing
    % output are CLUSTER_STRUCT.saved as mat
    % Copyright (C) Chang Liu,
    % Copyright (C) Jacob Salminen, jsalminen@ufl.edu
    % code adapted from bemobile pipeline & EEGLAB
    %## TIME
    tic;
    %## CHECKS
    %- 
    if params_in.clust_k_num < 2
        error('params_in.clust_k_num must be > 2 not %i',params_in.clust_k_num);
    end
    %- check which kmeans is currently in use
    flagstats = strcmp(regexp(which('kmeans'), '(?<=[\\/]toolbox[\\/])[^\\/]+', 'match', 'once'),'stats');
    if ~flagstats
        kmeansPath = fileparts(which('kmeans'));
        rmpath(kmeansPath);
        addpath(kmeansPath);
        fprintf('Using kmeans on path: %s',kmeansPath);
    end
    %##
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = (2:length(STUDY.cluster));
    else
        for index = 2:length(STUDY.cluster)
            if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) && ~strncmpi('Notclust',STUDY.cluster(index).name,8)
                rmindex = [rmindex index ];
            end
        end      
    end
    if ~isempty(rmindex)
        fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
        STUDY.cluster(rmindex)          = [];
        STUDY.cluster(clustlevel).child = [];
        if clustlevel == 1 && length(STUDY.cluster) > 1
            STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
        end
    end
    %- cluster using assigned algorithm
    clustdata = STUDY.etc.preclust.preclustdata;
    switch lower(params_in.algorithm)
        case {'kmeans','kmeanscluster'}
            if params_in.clust_k_repeat_std == Inf
                if strcmpi(params_in.algorithm, 'kmeans')
                    % Change to use Makoto's parameter                
                    [idx,C,sumd,D] = kmeans(clustdata,params_in.clust_k_num,...
                        'EmptyAction','singleton',...
                        'maxiter',params_in.clust_k_maxiter,...
                        'replicates',params_in.clust_k_replicates);
                else
                    [idx,C,sumd,D] = kmeanscluster(clustdata,params_in.clust_k_num);
                end
                [STUDY] = std_createclust(STUDY,ALLEEG,...
                    'clusterind',idx,...
                    'algorithm',{'Kmeans',params_in.clust_k_num});
            else
                % (09/08/2024) JS, using modified robust_means (author,
                % Chang Liu) instead of original.
                % [IDX,C,sumd,D,outliers_out] = robust_kmeans(clustdata,params_in.clust_k_num,params_in.clust_k_repeat_std,ROBUST_KMEANS_MAXITER,params_in.algorithm);
                [idx,C,sumd,D,outliers_out] = robust_kmeans_override(clustdata,...
                    params_in.clust_k_num,...
                    params_in.clust_k_repeat_std,...
                    params_in.clust_k_robust_maxiter,...
                    params_in.algorithm,...
                    params_in.clust_k_maxiter,...
                    params_in.clust_k_replicates,...
                    params_in.clust_k_empty_action);
                [STUDY] = std_createclust(STUDY,ALLEEG,...
                    'clusterind',idx,...
                    'algorithm',{'robust_kmeans', params_in.clust_k_num});
            end
        case 'neural network'
            [idx,C] = neural_net(clustdata,params_in.clust_k_num);
            [STUDY] = std_createclust(STUDY,ALLEEG,...
                'clusterind',idx,...
                'algorithm',{'Neural Network', params_in.clust_k_num});
        case 'affinity propogation'           
            [idx,C,sumd] = std_apcluster(clustdata,'maxits',params_in.maxiter);
            [STUDY]      = std_createclust(STUDY,ALLEEG,...
                'clusterind',idx,...
                'algorithm',{'Affinity Propagation',size(C,1)});
    otherwise
            disp('pop_clust: unknown params_in.algorithm %s return',params_in.algorithm);
            return
    end
    %- store clustering outcome        
    cluster_outcome.IDX = idx;
    cluster_outcome.C = C;
    cluster_outcome.sumd = sumd;
    cluster_outcome.D = D;
    cluster_outcome.outliers = [];
    if exist('outliers_out','var')
        cluster_outcome.outliers = outliers_out;
    end
    %##
    toc
end
%% (SUBFUNCTION) ======================================================= %%
function [fig] = plot_k_evals(STUDY,cluster_idx,save_dir)
    %% (STEP 1) CALCULATE CLUSTER SOLUTIONS
    %- find the optimal number of clusters without generating outliers
    %- (07/12/2023) CL, 30 is the maximum number clusters typically
        % used for EEG (see Makoto) % I think 5 is a reasonable lower bound
        % updated to do repeated clustering clust_CL is a function
        % that is similar as pop_clust that does k_means cluster 
        % but not generate cluster struct in the STUDY; clust_CL 
        % used kmeans to cluster for 1000 times
    %- Calculate evaulation criteria for K (kmeans alg)
    eva1 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx+1, 'silhouette'); % this is the method to find optimal number of clusters (confirmed by EEGlab mailing list)
    eva2 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx+1, 'CalinskiHarabasz');
    eva3 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx+1, 'DaviesBouldin');
    %- Calculate evaulation criteria
    fig = figure();
    %-
    subplot(1,3,1)
    plot(eva1.InspectedK,eva1.CriterionValues,'-o');hold on;
    plot(eva1.OptimalK,eva1.CriterionValues(eva1.InspectedK == eva1.OptimalK),'o','color','r');
    ylim([0,max(eva1.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('Silhouette');
    %-
    subplot(1,3,2)
    plot(eva2.InspectedK,eva2.CriterionValues,'-o');hold on;
    plot(eva2.OptimalK,eva2.CriterionValues(eva2.InspectedK == eva2.OptimalK),'o','color','r');
    ylim([0,max(eva2.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('CalinskiHarabasz');
    %-
    subplot(1,3,3)
    plot(eva3.InspectedK,eva3.CriterionValues,'-o');hold on;
    plot(eva3.OptimalK,eva3.CriterionValues(eva3.InspectedK == eva3.OptimalK),'o','color','r');
    ylim([0,max(eva3.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('DaviesBouldin');
    saveas(gcf,[save_dir filesep 'cluster_kmeans_eval.fig'])
    saveas(gcf,[save_dir filesep 'cluster_kmeans_eval.jpg'])
end
%% (SUBFUNCTION) ======================================================= %%
% Editted by Chang Liu - 2023-01-25 to increase the number of iterations of
% kmeans
% (09/08/2024) JS, edited to increase functionality of parameters
function  [IDX,C,sumd,D,outliers] = robust_kmeans_override(data,N,STD,max_iter,method,k_max_iter,k_reps,k_empty_action)
    %#ok<*ASGLU>     
    % data - pre-clustering data matrix.
    % N - number of wanted clusters.
    % k_max_iter = 10000; - iterations to cluster locally
    % k_reps = 1000; - replicates of local clustering
    % k_empty_action = 'drop'; - what to do with outliers
    % method = 'kmeans'; - algorithm to use.
    flag  = 1;
    not_outliers = 1:size(data,1);
    old_outliers = [];
    if strcmpi(method, 'kmeans')
        [IDX,C,sumd,D] = kmeans(data,N,...
            'replicates',k_reps,...
            'maxiter',k_max_iter,...
            'emptyaction',k_empty_action); % Cluster using K-means algorithm % CL: editted
    else
        [IDX,C,sumd,D] = kmeanscluster(data,N); % Cluster using K-means algorithm
    end   
    if STD >= 2 % STD for returned outlier
        rSTD = STD - 1;
    else
        rSTD = STD;
    end
    loop = 0;
    while flag
        loop =  loop + 1;
	    std_all = [];
        ref_D = 0;
	    for k = 1:N
            tmp = ['cls' num2str(k) ' = find(IDX=='  num2str(k) ')''; ' ]; %find the component indices belonging to each cluster (cls1 = ...).
            eval(tmp);
            tmp = ['std' num2str(k) ' = std(D(cls'  num2str(k) ' ,' num2str(k) ')); ' ]; %compute the std of each cluster
            eval(tmp);
            std_all = [std_all ['std' num2str(k)  '  ']];
            tmp = [ 'ref_D = ' num2str(ref_D) ' + mean(D(cls'  num2str(k) ' ,' num2str(k) '));' ];
            eval(tmp);
	    end
	    std_all = [ '[ ' std_all ' ]' ];
        std_all = eval(std_all);
        
	    % Find the outliers
        % Outlier definition - its distance from its cluster center is bigger
        % than STD times the std of the cluster, as long as the distance is bigger
        % than the mean distance times STD (avoid problems where all points turn to be outliers).
	    outliers = [];
        ref_D = ref_D/N;
	    for k = 1:N
            tmp = ['cls' num2str(k) '(find(D(find(IDX=='  num2str(k) ')'' , ' num2str(k) ') > ' num2str(STD)  '*std' num2str(k) ')); ' ];
            optionalO = eval(tmp);
            Oind = find(D(optionalO,k) >  ref_D*STD);
            outliers = [outliers optionalO(Oind)];
	    end
        if isempty(outliers) || (loop == max_iter)
            flag = 0;
        end
        l = length(old_outliers);
        returned_outliers = [];
     
        
        for k = 1:l
            tmp = sum((C-ones(N,1)*data(old_outliers(k),:)).^2,2)'; % Find the distance of each former outlier to the current cluster
            if isempty(find(tmp <= std_all*rSTD))  %Check if the outlier is still an outlier (far from each cluster center more than STD-1 times its std).
                returned_outliers = [returned_outliers old_outliers(k)];
            end
        end
        outliers = not_outliers(outliers);
        outliers = [outliers returned_outliers ];
	    tmp = ones(1,size(data,1));
	    tmp(outliers) = 0;
	    not_outliers = (find(tmp==1));
        
        if strcmpi(method, 'kmeans')
            [IDX,C,sumd,D] = kmeans(data(not_outliers,:),N,...
                'replicates',k_reps,...
                'emptyaction',k_empty_action);
        else
            [IDX,C,sumd,D] = kmeanscluster(data(not_outliers,:),N);
        end
        old_outliers = outliers;
        old_IDX = zeros(size(data,1),1);
        old_IDX(sort(not_outliers)) = IDX;
    end
    IDX = old_IDX;
end