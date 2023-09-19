function [all_solutions] = mim_cluster_process(STUDY,ALLEEG,save_dir,varargin)
%MIM_CLUSTER_PROCESS Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are CLUSTER_PARAMS.saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% code adapted from bemobile pipeline & EEGLAB
%## TIME
tic;
%## DEFINE DEFAULTS
%- high level parameters
KS_TO_TEST = (5:25);
MAK_MAXITER = 10000;
MAK_REPLICATE = 1000;
MAX_REPEATED_ITERATIONS = 50;
REPEATED_CLUSTERING_STD = 3;
ROBUST_KMEANS_MAXITER = 5;
%- clustering parameters
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',STUDY.filename,...
    'filepath',STUDY.filepath,...
    'outliers',inf(),...
    'maxiter',200);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'CLUSTER_PARAMS',CLUSTER_PARAMS,@isstruct);
addParameter(p,'REPEATED_CLUSTERING_STD',REPEATED_CLUSTERING_STD,@isnumeric);
addParameter(p,'MAX_REPEATED_ITERATIONS',MAX_REPEATED_ITERATIONS,@isnumeric);
addParameter(p,'KS_TO_TEST',KS_TO_TEST,@isnumeric);
parse(p,STUDY,ALLEEG,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
%## MAKE DIRS
CLUSTER_PARAMS = p.Results.CLUSTER_PARAMS;
REPEATED_CLUSTERING_STD = p.Results.REPEATED_CLUSTERING_STD;
MAX_REPEATED_ITERATIONS = p.Results.MAX_REPEATED_ITERATIONS;
KS_TO_TEST = p.Results.KS_TO_TEST;
%% ===================================================================== %%
cnt = 1;
% KS_TO_TEST = K_RANGE(1):K_RANGE(2);
% cluster_outcomes = cell(K_RANGE(2)-K_RANGE(1),1);
cluster_idx = zeros(length([STUDY.datasetinfo.comps]),length(KS_TO_TEST));
%## cluster using std_preclust weightings for each K desired (deafault here is kmeans with no outlier rejection).
for clust_i = KS_TO_TEST
    %- change cluster k
    CLUSTER_PARAMS.clust_num = clust_i;
    %- calculate clusters
    [~, clust_out] = cluster_dips(STUDY,ALLEEG,CLUSTER_PARAMS,MAK_MAXITER,MAK_REPLICATE,ROBUST_KMEANS_MAXITER);
    cluster_idx(1:length(clust_out.IDX),cnt) = clust_out.IDX;
    cnt = cnt+1;
end
cluster_idx = cluster_idx(:,all(cluster_idx ~= 0));
%##
plot_k_evals(STUDY,cluster_idx,save_dir);
%##
% (08/20/2023) JS, to force kmeans to kmeans_robust you must set an
% outliers std
CLUSTER_PARAMS.outliers = REPEATED_CLUSTERING_STD;
fprintf('Clustering %d times...\n', MAX_REPEATED_ITERATIONS)
all_solutions = cell(length(KS_TO_TEST),1);
parfor i = 1:length(KS_TO_TEST)
% for i = 1:length(cluster_ks)
    clust_i = KS_TO_TEST(i);
    %- temp variables
    params = CLUSTER_PARAMS;
    params.clust_num = clust_i;
    solutions = cell(MAX_REPEATED_ITERATIONS,1);
    %- iterate
    for iter = 1:MAX_REPEATED_ITERATIONS
        % start time
        tic
        fprintf('Iteration %d of %d \n',iter,MAX_REPEATED_ITERATIONS)
        % clustering according to parameters - Note by Chang Liu, basically
        % doing the robust_kmeans cluster in pop_clust.m, I prefer not
        % overwrite the original STUDY
        [TMP_STUDY,~] = cluster_dips(STUDY,ALLEEG,params,MAK_MAXITER,MAK_REPLICATE,ROBUST_KMEANS_MAXITER); 
        % store info in output struct
        solutions{iter} = TMP_STUDY.cluster;
        % stop checking the time and plot estimated time of arrival
        lastduration = toc; 
        eta = lastduration*(MAX_REPEATED_ITERATIONS-iter);
        fprintf('Last duration: %1.2f s \n',round(lastduration,2))
        fprintf('ETA: %d h, %d min \n', floor(eta/3600), round(((eta/3600)-floor(eta/3600))*60))
    end
    %- store clustering info
    param_struct = struct('preclustparams',TMP_STUDY.etc.preclust,...
        'outlier_sigma',CLUSTER_PARAMS.outliers,...
        'n_clust',clust_i,...
        'n_iters',MAX_REPEATED_ITERATIONS);
    tmp_clust_sol = struct('parameters',param_struct,...
        'solutions',{solutions});
    all_solutions{i} = tmp_clust_sol;
end
end
%% (SUBFUNCTION) ======================================================= %%
function [STUDY,cluster_outcome] = cluster_dips(STUDY,ALLEEG,CLUSTER_PARAMS,MAK_MAXITER,MAK_REPLICATE,ROBUST_KMEANS_MAXITER)
    %MIM_CLUSTERING_ANL Summary of this function goes here
    % CAT CODE
    %  _._     _,-'""`-._
    % (,-.`._,'(       |\`-/|
    %     `-.-' \ )-`( , o o)
    %           `-    \`_`"'-
    % Code Designers: Chang Liu, Jacob Salminen
    % Code Date: 07/17/2023, MATLAB 2020b
    % Written by Chang - 2023-4-15 to run plotERSP with parallel processing
    % output are CLUSTER_PARAMS.saved as mat
    % Copyright (C) Chang Liu,
    % Copyright (C) Jacob Salminen, jsalminen@ufl.edu
    % code adapted from bemobile pipeline & EEGLAB
    %## TIME
    tic;
    
    %## CHECKS
    %- 
    if CLUSTER_PARAMS.clust_num < 2
        error('CLUSTER_PARAMS.clust_num must be >2 not %i',CLUSTER_PARAMS.clust_num);
    end
    %- check which kmeans is currently in use
    flagstats = strcmp(regexp(which('kmeans'), '(?<=[\\/]toolbox[\\/])[^\\/]+', 'match', 'once'),'stats');
    if ~flagstats
        kmeansPath = fileparts(which('kmeans'));
        rmpath(kmeansPath);
        addpath(kmeansPath);
    end
    %##
    %- 
    rmindex    = [];
    clustlevel = STUDY.etc.preclust.clustlevel;
    nameclustbase = STUDY.cluster(clustlevel).name;
    if clustlevel == 1
        rmindex = [2:length(STUDY.cluster)];
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
    switch lower(CLUSTER_PARAMS.algorithm)
        case {'kmeans','kmeanscluster'}
            if CLUSTER_PARAMS.outliers == Inf
                if strcmpi(CLUSTER_PARAMS.algorithm, 'kmeans')
                    % Change to use Makoto's parameter                
                    [IDX,C,sumd,D] = kmeans(clustdata,CLUSTER_PARAMS.clust_num,...
                        'emptyaction','singleton',...
                        'maxiter',MAK_MAXITER,...
                        'replicate',MAK_REPLICATE);
                else
                    [IDX,C,sumd,D] = kmeanscluster(clustdata,CLUSTER_PARAMS.clust_num);
                end
                [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Kmeans', CLUSTER_PARAMS.clust_num});
            else
                [IDX,C,sumd,D,outliers_out] = robust_kmeans(clustdata,CLUSTER_PARAMS.clust_num,CLUSTER_PARAMS.outliers,ROBUST_KMEANS_MAXITER,CLUSTER_PARAMS.algorithm);
                [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'robust_kmeans', CLUSTER_PARAMS.clust_num});
            end
        case 'neural network'
            [IDX,C] = neural_net(clustdata,CLUSTER_PARAMS.clust_num);
            [STUDY] = std_createclust(STUDY,ALLEEG,...
                'clusterind',IDX,...
                'algorithm',{'Neural Network', CLUSTER_PARAMS.clust_num});
        case 'affinity propogation'           
            [IDX,C,sumd] = std_apcluster(clustdata,'maxits',CLUSTER_PARAMS.maxiter);
            [STUDY]      = std_createclust(STUDY,ALLEEG,...
                'clusterind',IDX,...
                'algorithm',{'Affinity Propagation',size(C,1)});
    otherwise
            disp('pop_clust: unknown CLUSTER_PARAMS.algorithm %s return',CLUSTER_PARAMS.algorithm);
            return
    end
    %- store clustering outcome        
    cluster_outcome.IDX = IDX;
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
        % used for EEG (see Makyoto) % I think 5 is a reasonable lower bound
        % updated to do repeated clustering clust_CL is a function
        % that is similar as pop_clust that does k_means cluster 
        % but not generate cluster struct in the STUDY; clust_CL 
        % used kmeans to cluster for 1000 times
    %- Calculate evaulation criteria for K (kmeans alg)
    eva1 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'silhouette'); % this is the method to find optimal number of clusters (confirmed by EEGlab mailing list)
    eva2 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'CalinskiHarabasz');
    eva3 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'DaviesBouldin');
    %- Calculate evaulation criteria
    fig = figure();
    subplot(1,3,1)
    plot(eva1.InspectedK,eva1.CriterionValues,'-o');hold on;
    plot(eva1.OptimalK,eva1.CriterionValues(eva1.InspectedK == eva1.OptimalK),'o','color','r');
    ylim([0,max(eva1.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('Silhouette');
    subplot(1,3,2)
    plot(eva2.InspectedK,eva2.CriterionValues,'-o');hold on;
    plot(eva2.OptimalK,eva2.CriterionValues(eva2.InspectedK == eva2.OptimalK),'o','color','r');
    ylim([0,max(eva2.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('CalinskiHarabasz');
    subplot(1,3,3)
    plot(eva3.InspectedK,eva3.CriterionValues,'-o');hold on;
    plot(eva3.OptimalK,eva3.CriterionValues(eva3.InspectedK == eva3.OptimalK),'o','color','r');
    ylim([0,max(eva3.CriterionValues)+0.05])
    xlabel('Clusters');ylabel('DaviesBouldin');
    saveas(gcf,[save_dir filesep 'cluster_kmeans_eval.fig'])
    saveas(gcf,[save_dir filesep 'cluster_kmeans_eval.jpg'])
end