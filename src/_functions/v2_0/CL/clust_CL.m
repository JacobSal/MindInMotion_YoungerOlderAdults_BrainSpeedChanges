% k means cluster without using the pop_clust
% This code is adapted from pop_clust from EEGLAB

% Chang Liu 2021-11-26

% Check that path to the stats toolbox comes first (conflict with Fieldtrip)
function [STUDY,ClusterOutcome] = clust_CL(STUDY, ALLEEG, varargin)

flagstats = strcmp(regexp(which('kmeans'), '(?<=[\\/]toolbox[\\/])[^\\/]+', 'match', 'once'),'stats');
if ~flagstats
    kmeansPath = fileparts(which('kmeans'));
    rmpath(kmeansPath);
    addpath(kmeansPath);
end

rmindex    = [];
clustlevel = STUDY.etc.preclust.clustlevel;
nameclustbase = STUDY.cluster(clustlevel).name;
if clustlevel == 1
    rmindex = [2:length(STUDY.cluster)];
else
    for index = 2:length(STUDY.cluster)
        if strcmpi(STUDY.cluster(index).parent{1}, nameclustbase) && ~strncmpi('Notclust',STUDY.cluster(index).name,8)
            rmindex = [ rmindex index ];
        end
    end;        
end
if ~isempty(rmindex)
    fprintf('Removing child clusters of ''%s''...\n', nameclustbase);
    STUDY.cluster(rmindex)          = [];
    STUDY.cluster(clustlevel).child = [];
    if clustlevel == 1 && length(STUDY.cluster) > 1
        STUDY.cluster(1).child = { STUDY.cluster(2).name }; % "Notclust" cluster
    end
end

%default values
algorithm = 'kmeans';
clus_num  = 20;
save     = 'off';
filename = STUDY.filename;
filepath = STUDY.filepath;
outliers = Inf; % default std is Inf - no outliers
maxiter  = 200;
    
if mod(length(varargin),2) ~= 0
    error('pop_clust(): input variables must be specified in pairs: keywords, values');
end

for k = 1:2:length(varargin)
    switch(varargin{k})
        case 'algorithm'
            algorithm  = varargin{k+1};
        case 'clus_num'
            clus_num   = varargin{k+1};
        case 'outliers'
            outliers   =  varargin{k+1};                
        case 'save'
            save       = varargin{k+1};
        case 'filename' 
            filename   = varargin{k+1};
        case 'filepath'
            filepath   = varargin{k+1};
        case 'maxiter' 
            maxiter    = varargin{k+1};
    end
end    
if clus_num < 2
    clus_num = 2;
end

clustdata = STUDY.etc.preclust.preclustdata;
switch lower(algorithm)
    case { 'kmeans' 'kmeanscluster' }
        if outliers == Inf
            if strcmpi(algorithm, 'kmeans')
%                 [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'replicates',10,'emptyaction','drop');
                % Change to use Makoto's parameter                
                [IDX,C,sumd,D] = kmeans(clustdata,clus_num,'emptyaction', 'singleton', 'maxiter', 10000, 'replicate', 1000);
                
            else
                [IDX,C,sumd,D] = kmeanscluster(clustdata,clus_num);
            end
%             [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Kmeans', clus_num});
        else
            [IDX,C,sumd,D,outliers] = robust_kmeans(clustdata,clus_num,outliers,5, algorithm);
%             [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'robust_kmeans', clus_num});
        end
    case 'neural network'
        [IDX,C] = neural_net(clustdata,clus_num);
        [STUDY] = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm',  {'Neural Network', clus_num});
    case 'Affinity Propagation'           
         [IDX,C,sumd] = std_apcluster(clustdata,'maxits',maxiter);
         [STUDY]      = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Affinity Propagation',size(C,1)});
    otherwise
        disp('pop_clust: unknown algorithm return');
        return
end
        % If save updated STUDY to disk
% if strcmpi(save,'on')
%     if (~isempty(STUDY.filename)) && (~isempty(STUDY.filepath))
%         STUDY = pop_savestudy(STUDY, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
%     else
%         STUDY = pop_savestudy(STUDY);
%     end
% end      
ClusterOutcome.IDX = IDX;
ClusterOutcome.C = C;
ClusterOutcome.sumd = sumd;
ClusterOutcome.D = D;
if exist('outliers')
    ClusterOutcome.outliers = outliers;
end
end

