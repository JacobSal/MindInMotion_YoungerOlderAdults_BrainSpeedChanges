% Create study for cluster ICs. This code only works for cluster without
% using ERSP. Precompute ERSP needed to be done on Hipergator
% Chang Liu - 2021-11-23 - V1

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%   NJacobsen notes
%   When timewarping data, save values as EEG.timewarp = timewarp;
%   EEG.timewarp.medianlatency = median(timewarp.latencies(:,:));%Warping to the median latency of my 5 events
%   By default, std_ersp will use the median of all subject's
%   timewarp.latencies(:,:) as 'timewarpms' unless individual subject 
%   warpto is indiciated using 'timewarpms', 'subject tw matrix'
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230417.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM/run_b_ic_clustering.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (PARAMETERS) ======================================================== %%
%## Hard Define
%- timewarping params
TW = 1; %Time warping option
preclust = 1; %1 - precluster study; 0 - no preclustering
clust = 1; %1 - precluster study; 0 - no preclustering
precomp_nonERSPs = 0; %1 - pulls up precompute gui; 0 - no gui
precompute_ERSP = 1; %1 - precompute ersps; 0 - don't precompute ersps
precompute_spec_scalp = 1;
erspComp='light'; %'light' - quicker computation; 'full' - with usual parameters (takes longer)
showClusterPlotGUI=1; %1 - show cluster plot gui at end; 0 - don't show it (and clear study)
groupmedian_timewarpms = 1 ;
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 5;
clustering_weights.scalp = 5;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
cluster_alg = 'kmeans';
do_multivariate_data = 1;
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
%- datetime override
dt = '04092023_MIM_OA_subset_N85_speed_terrain';
%- hard define
% load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
load_trials = {'0p25','0p5','0p75','1p0'};
% load_trials = {'flat','low','med','high'}; 
%## soft define
subjinfDir = [SUBJINF_DIR filesep sprintf('%s',dt)];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'ic_rejection_crit'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
% --------------- SETUP PATH ----------------------
Process_Imagined = 0;
folder_name = 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10'
Study_folder = folder_name;
filepath = 'STUDY-preprocess-HY_202212';
studyname = [folder_name,'.study'];
MiM_HY_config_params;
brain_score = 8;
load_folder = fullfile('M:\liu.chang1\',filepath, folder_name);
save_study_folder = fullfile('M:\liu.chang1\',filepath,folder_name,'BATCH-3-Epoch'); 
save_epoch_folder = fullfile('M:\liu.chang1\',filepath,folder_name,'BATCH-3-Epoch');    
load_IC_rej = fullfile('M:\liu.chang1\',filepath,folder_name,'Summary');    
    
% mkdir(fullfile(save_study_folder,num2str(subDirStudy)));
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(load_trials)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
study_fName = sprintf('%s_MIM_study',[load_trials{:}]);
if ~exist([load_dir filesep study_fName '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);
    end
end
all_subjStr = {ALLEEG.subject};
%%
if ~load_study
    EEG=[]; ALLEEG=[]; CURRENTSET=[]; ALLCOM=[]; CURRENTSTUDY = []; STUDY = [];

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    indx = 1;
    disp('Check do not use eeglab on R drive');
    % Create Study
    for i = [1:23]
        i
        subjStr = all_subjStr{i};
        if exist(fullfile(load_folder,subjStr,[subjStr,'_Epoched.set']))
            [STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name',studyname,...,
                'commands',...
                {{'index',indx,'load',fullfile(load_folder,subjStr,[subjStr,'_Epoched.set']),...
                'subject',subjStr},{'inbrain' 'on' 'dipselect' 0.15}},'updatedat','on');
            subjToAnalyze{indx} = subjStr;
            indx = indx + 1;
        end
    end
    [STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    eeglab redraw

    %Make STUDY design
    [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'subjselect', subjToAnalyze,'variable1','cond','values1', {'flat','low','med','high'});
    [STUDY] = std_makedesign(STUDY, ALLEEG, 2, 'subjselect', subjToAnalyze,'variable1','cond','values1', {'0p25','0p5','0p75','1p0'});

    % Save STUDY
    [STUDY ALLEEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',save_study_folder);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
end

%% Precompute Measures
% 1. Precompute ERSP: Note that mod_std_precomp_v10_2_5_5a no longer works
% with eeglab2021. It seems that precompute ERSP takes forever on local
% computer. Needs to use hipergator for parallel processing
% Notes: cycle = wavelet cycle. 

% ------------------------------------------ 
% The following code will DO STATSSSSS!!!
%{
if precompute_ERSP
    if TW == 0 % CL: Don't really understand this part yet
        switch erspComp
            case 'full'
                tic
                [STUDY ALLEEG] = std_precomp_CL_eeglab2021(STUDY, ALLEEG, 'components','ersp','on','itc','on','erspparams',{'cycles',[3 0.5],'alpha',0.05, 'padratio',2,'savetrials','off','baseline',nan}, 'recompute','on');
                toc
            case 'light'
                tic
                %Parameters for quicker computation (still works quite well)
                [STUDY ALLEEG] =std_precomp_CL_eeglab2021(STUDY, ALLEEG, 'components','ersp','on','itc','off','erspparams',{'cycles',[3 0.9],'freqs',[3 256],'nfreqs',100,'alpha',0.05,'freqscale','log','savetrials','off','baseline',nan}, 'recompute','on');
                toc
            otherwise
                error('Incorrect case for erspComp!');
        end
    elseif TW == 1 %CL Also don't understand this part
        %Run ERSPs (timewarped)
        if groupmedian_timewarpms == 1 
            warps=zeros(length(ALLEEG),length(ALLEEG(1,1).timewarp.warpto));%NJ; stored in ALLEEG.timewarp.warpto
            for i=1:length(ALLEEG)
                warps(i,:)=ALLEEG(1,i).timewarp.warpto; %NJ; stored in ALLEEG.timewarp.warpto
            end%         warpingvalues=median(warps);
            roundNear=50; %round numbers to the closest multiple of this value
            warpingvalues = round(median(warps)/roundNear)*roundNear;
        elseif groupmedian_timewarpms == 2 %use subject specific warpto values
            warpingvalues = 'subject tw matrix';
        end
        
        switch erspComp
            case 'full'
                tic
                [STUDY ALLEEG] = std_precomp_CL_eeglab2021(STUDY, ALLEEG, 'components','ersp','on','itc','on',...
                    'erspparams',{'cycles',[3 0.5],'alpha',0.05, 'padratio',2,'savetrials','off','baseline',nan,'timewarp','subject tw matrix','timewarpms', warpingvalues},...
                    'recompute','on');
                toc
            case 'light'
                tic
                %Parameters suggested by Makoto
                [STUDY ALLEEG] = std_precomp_CL_eeglab2021(STUDY, ALLEEG(1), 'components','ersp','on','itc','off','recompute','on',...
                    'erspparams',{'cycles',[3 0.9],'freqs',[3 256],'nfreqs',100,'alpha',0.05,'freqscale','log','savetrials','off','baseline',nan,'timewarp','subject tw matrix',...
                    'timewarpms', warpingvalues});
                toc
            otherwise
                error('Incorrect case for erspComp!');
        end
    end
end


%-----------------------------------------------
% This following code does not run stats for precompute ERSP
if TW == 1 %CL Also don't understand this part
    %Run ERSPs (timewarped)
    if groupmedian_timewarpms == 1 
        warps=zeros(length(ALLEEG),length(ALLEEG(1,1).timewarp.warpto));%NJ; stored in ALLEEG.timewarp.warpto
        for i=1:length(ALLEEG)
            warps(i,:)=ALLEEG(1,i).timewarp.warpto; %NJ; stored in ALLEEG.timewarp.warpto
        end%         warpingvalues=median(warps);
        roundNear=50; %round numbers to the closest multiple of this value
        warpingvalues = round(median(warps)/roundNear)*roundNear;
    elseif groupmedian_timewarpms == 2 %use subject specific warpto values
        warpingvalues = 'subject tw matrix';
    end
    % --------- copied/modified from Ryan's code
    parpool(2)
    for set_i = 1%:length(subjToAnalyze) %loop through all subjects in study
        CurrentSubjStr = subjToAnalyze(set_i).name;
%         icaFileName = [CurrentSubjStr '.icatimef'];
        
        tempStudy = STUDY; %make temp fake study to prevent errors later when we send in just one subject at a time
        tempStudy.datasetinfo = tempStudy.datasetinfo(set_i); %just current subject
        tempStudy.datasetinfo.index = 1; %change index to prevent indexing error in eeglab later (e.g. if we temporarily grab the third subject, they used to have index 3 but since they are the only subject now, they need to be index 1
        tempStudy.subject = tempStudy.subject(set_i);
        erspTic = tic;
%         [~, ~] = std_precomp(tempStudy, ALLEEG(set_i), 'components','savetrials','on','recompute',recomputeParam,'ersp','on',...
%             'erspparams',{'cycles',[3 0.8] ,'nfreqs',100,'ntimesout',200,'timewarp',ALLEEG(set_i).timewarp.latencies,'timewarpms', grandAvgWarpTo},'itc','off'); %ERSP
%         disp(['Successfully calculated ERSP (component) for ' CurrentSubjStr ' !!']);
      
        [~, ~] = std_precomp(tempStudy, ALLEEG(set_i), 'components', 'savetrials','on','recompute','on','ersp','on',...
        'erspparams',{'parallel','on','cycles',[3 0.8],'nfreqs',100,'ntimesout',200,'baseline',nan,'timewarp',ALLEEG(set_i).timewarp.latencies,'timewarpms', warpingvalues},'itc','off'); %ERSP
        disp(['Successfully calculated ERSP (component) for ' CurrentSubjStr ' !!']);
        
        erspToc = toc(erspTic);
        disp(['It took ',num2str(erspToc),' seconds to calculate an ERSP for this subject']);
        
    end

end
%}
%-----------------------------------------------
% 2. Precompute Spectrum
% if precompute_spec_scalp
%     freqrange = [1 256];
%     [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components',...
%     'scalp','on',...
%     'spec','off','specparams',{'freqrange' freqrange 'specmode' 'fft' 'logtrials' 'off'});
% end
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components',...
'scalp','on');
%% Alternatively, use the precomputeMeasures function

%% perform preclustering, add weights to the params used for clustering
command = '[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1';

if isfield(clustering_weights,'dipoles') && clustering_weights.dipoles ~= 0
    command = [command ',{''dipoles'',''weight'',clustering_weights.dipoles}'];
end
if isfield(clustering_weights,'scalp') && clustering_weights.scalp ~= 0
    command = [command ',{''scalp'',''npca'',10,''weight'',clustering_weights.scalp,''abso'',1}'];
end
if isfield(clustering_weights,'spec') && clustering_weights.spec ~= 0
    command = [command ',{''spec'',''npca'',10,''weight'',clustering_weights.spec,''freqrange'',[3 45] }'];
end
if isfield(clustering_weights,'erp') && clustering_weights.erp ~= 0
    command = [command ',{''erp'',''npca'',10,''weight'',clustering_weights.erp,''timewindow'',timewindow,''erpfilter'',''''}'];
end
if isfield(clustering_weights,'ersp') && clustering_weights.ersp ~= 0
    command = [command ',{''ersp'',''npca'',10,''freqrange'',freqrange,''timewindow'',timewindow,''norm'',1,''weight'',clustering_weights.ersp}'];
end
command = [command ');'];

eval(command)

% store essential info in STUDY struct for later reading
freqrange = [3 45];
STUDY.etc.clustering.preclustparams.clustering_weights = clustering_weights;
STUDY.etc.clustering.preclustparams.freqrange = freqrange;
% STUDY.etc.clustering.preclustparams.timewindow = timewindow;

%% Perform another pass to remove bad ICs
for i = 1:length(STUDY.datasetinfo)
    subjStr = STUDY.datasetinfo(i).subject;
    comp_study = STUDY.datasetinfo(i).comps;
    % Remove additional components based on the PowpowCat xlsx
    PowpowCat_Rej = MiM_config.PowpowCat_Rej;
    badIC = str2num(PowpowCat_Rej.IC_bad(find(strcmp(MiM_config.PowpowCat_Rej.SubjStr,subjStr))))
    comp_valid_all = setdiff(comp_study,badIC);
    STUDY.datasetinfo(i).comps = comp_valid_all ;
end
%% Perform Cluster - partly adapt from bemobil_repeated_clustering

% if precomp_nonERSPs==1
%     [STUDY ALLEEG] = pop_precomp(STUDY, ALLEEG,'components');
% end
% 
% %Precluster components
% if preclust == 1
%     [STUDY ALLEEG] = pop_preclust(STUDY,ALLEEG);
% end

% First compute the average ICs across participants
if load_study
    evaluate_clust = 1;
    lookup_atlas = 1;
    evaluate_method = 'min_rv';
    save_clutering_evaluate = [clustering_method,'_',evaluate_method];

    for n = 1:length(STUDY.datasetinfo)
        numIC(n) = size(STUDY.datasetinfo(n).comps,2);
    end
    fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
    mean_IC_allSub = floor(mean(numIC));
    
    %Cluster components
    filename_clustering_solutions = fullfile(save_study_folder, 'clustering_solution.mat');
    n_iterations = 1;
    n_clust = 10;
    outlier_sigma = 3;
    if clust 
        switch cluster_alg
            case 'affinity' % I am not sure how this works
            %   [STUDY] = pop_clust(STUDY,ALLEEG);
                [IDX,C,sumd] = std_apcluster(STUDY.etc.preclust.preclustdata,'maxits',200);
                [STUDY]      = std_createclust(STUDY, ALLEEG, 'clusterind', IDX, 'algorithm', {'Affinity Propagation',size(C,1)});
            case 'kmeans' 
                % Step 1: find the optimal number of clusters without
                % generating outliers
                clear cluster_idx
                p = 1;
                for num_clust = 5:mean_IC_allSub %30 is the maximum number clusters typically used for EEG (see Makyoto) % I think 5 is a reasonable lower bound
                    % updated to do repeated clustering
                    % clust_CL is a function that is similar as pop_clust
                    % that does k_means cluster but not generate cluster
                    % struct in the STUDY; clust_CL used kmeans to cluster
                    % for 1000 times
                    [~, ClusterOutcome{num_clust}] = clust_CL(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  num_clust , 'outliers',  Inf );
                    cluster_idx(:,p) = ClusterOutcome{num_clust}.IDX;
                    p = p+1;
                end
                eva1 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'silhouette') % this is the method to find optimal number of clusters (confirmed by EEGlab mailing list)
                eva2 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'CalinskiHarabasz');
                eva3 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'DaviesBouldin');

                figure();
                subplot(1,3,1)
                plot(eva1.InspectedK,eva1.CriterionValues,'-o');hold on;
                plot(eva1.OptimalK,eva1.CriterionValues(eva1.InspectedK == eva1.OptimalK),'o','color','r');
                xlabel('Clusters');ylabel('Silhouette');
                subplot(1,3,2)
                plot(eva2.InspectedK,eva2.CriterionValues,'-o');hold on;
                plot(eva2.OptimalK,eva2.CriterionValues(eva2.InspectedK == eva2.OptimalK),'o','color','r');
                xlabel('Clusters');ylabel('CalinskiHarabasz');
                subplot(1,3,3)
                plot(eva3.InspectedK,eva3.CriterionValues,'-o');hold on;
                plot(eva3.OptimalK,eva3.CriterionValues(eva3.InspectedK == eva3.OptimalK),'o','color','r');
                xlabel('Clusters');ylabel('DaviesBouldin');
                saveas(gcf, fullfile(save_study_folder,' Cluster Eva.fig'))

                % Step 2  save clustering_solutions for selected number of
                % clusters
                for num_clust = 5:mean_IC_allSub
                    clear clustering_solution cluster_update
                    % Note: the clustering solutions are not exactly the same but not super different
                    clustering_solutions = repeated_clustering(STUDY,ALLEEG, n_iterations, num_clust, outlier_sigma, STUDY.etc.clustering.preclustparams);
%                   clustering_solutions2 = repeated_clustering(STUDY,ALLEEG, n_iterations, num_clust, outlier_sigma, STUDY.etc.clustering.preclustparams);
                    mkdir(fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(num_clust)))
                    save(fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(num_clust),['clustering_solutions',num2str(num_clust),'.mat']),'clustering_solutions');
                    if evaluate_clust
                        switch evaluate_method
                            case 'min_rv'
                                [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
                                %% Look up cluster centroid Brodmann area
                                if lookup_atlas
                                    atlas = ft_read_atlas('D:\Dropbox (UFL)\Postdoc_Code\func\fieldtrip-20210910\template\atlas\aal\ROI_MNI_V4.nii');
                                    for i = 3:length(cluster_update)
                                        cfg            = [];
                                        cfg.roi        = cluster_update(i).dipole.posxyz;
                                        cfg.output     = 'multiple';
                                        cfg.atlas      = atlas;
                                        cfg.inputcoord = 'mni';
                                        cfg.sphere = 1;
                                        labels = ft_volumelookup(cfg, atlas);
                                        [~, indx] = max(labels.count);
                                        Atlas_name{i,1} = labels.name(indx)
                                        cluster_update(i).analabel = labels.name(indx);
                                    end
                                end
                                %     [STUDY] = pop_savestudy(STUDY, EEG, 'filename',studyname,'filepath',save_study_folder);
                                mkdir(fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(num_clust),evaluate_method));
                                save(fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(num_clust),evaluate_method,['cluster_update',num2str(num_clust),'.mat']),'cluster_update');

                            case 'average'
                        end

                    end
                end
%                     for solution = 1:(length(fields(clustering_solutions))-1)
%                         STUDY.cluster = clustering_solutions.(['solution_' num2str(solution)]);  % Study.cluster include the last iteration result   
%                         % find dipole locations and centroids
%     %                     STUDY.clusterOutcome{solution} = bemobil_dipoles(STUDY,ALLEEG);
%                     end
         end
   end
end

% keyboard
% [STUDY ALLEEG] = pop_savestudy(STUDY, EEG, 'filename',studyname,'filepath',save_study_folder);

%% Identify clusters with more than half of the participants
%keyboard
for pick_cluster = 5:mean_IC_allSub
    close all
    outputdir = fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(pick_cluster),evaluate_method);
    load(fullfile(outputdir,['cluster_update',num2str(pick_cluster),'.mat']))
    STUDY.cluster = cluster_update;
    numSubj_per_cluster = cellfun(@length,cellfun(@unique, {STUDY.cluster.sets},'UniformOutput', false),'UniformOutput', false);
    numSubj_per_cluster = cell2mat(numSubj_per_cluster);
    valid_cluster = find(numSubj_per_cluster(3:end)>=0.5*(length(STUDY.subject)))+2;
    %% Update 7/20/2022 - Plot clustering results for AHA proposal
    % Plot scalp topographs which also need to be averaged? 
    if ~isfield(STUDY.cluster,'topo'), STUDY.cluster(1).topo = []; end
    for clus = 3:length(STUDY.cluster) % For each cluster requested
        if isempty(STUDY.cluster(clus).topo)
            STUDY = std_readtopoclust_CL(STUDY,ALLEEG, clus);% Using this custom modified code to allow taking average within participant for each cluster
        end
    end
    figure;
    std_topoplot_CL(STUDY,3:length(STUDY.cluster),'together');
    set(gcf,'position',[16 582 1340 751],'color','w')
    saveas(gcf,fullfile(outputdir,'Cluster_topo_plot_averaged.fig'));

    % Plot all dipole fit locations
    std_dipplot(STUDY,ALLEEG,'clusters','all','figure','off');
    set(gcf,'position',[16 582 1340 751],'color','w')

    % Plot dipole fit locations after averaged within participants
    std_dipplot_CL(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster),'figure','off','mode','together_averaged');
    set(gcf,'position',[16 582 1340 751],'color','w')
    saveas(gcf,fullfile(outputdir,'Cluster_dipole_plot_averaged.fig'));

    STUDY.etc.dipparams.centrline = 'off';
    std_dipplot_CL(STUDY,ALLEEG,'clusters',valid_cluster,'figure','off','mode','together_averaged_only','spheres','off','projlines','off');
    set(gcf,'position',[16 582 500 300],'color','w')
    camzoom(1.2^2);
    saveas(gcf,fullfile(outputdir,'Cluster_dipole_plot_allaveraged_only.fig'));
    % Plot dipole fit locations after averaged within participants and
    % different clusters have different colors
    STUDY.etc.dipparams.centrline = 'off';
    std_dipplot_CL(STUDY,ALLEEG,'clusters',valid_cluster,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
    set(gcf,'position',[16 582 500 300],'color','w')
    camzoom(1.2^2);
    saveas(gcf,fullfile(outputdir,'Cluster_dipole_plot_allaveraged.fig'));
    % Save STUDY final time
    % [STUDY ~] = pop_savestudy(STUDY, EEG, 'filename',[studyname,'_allInfo'],'filepath',save_study_folder);
    % CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
end
%% Plot clustering results
%{
figure;
std_topoplot(STUDY,ALLEEG,'clusters','all','figure','off');
set(gcf,'position',[16 582 1340 751],'color','w')
saveas(gcf, fullfile(save_study_folder, 'Cluster_topoplot.fig'));

figure;
std_dipplot(STUDY,ALLEEG,'clusters','all','figure','off');
set(gcf,'position',[16 582 1340 751],'color','w')
saveas(gcf, fullfile(save_study_folder, 'Cluster_dipplot.fig'));

% plot spectrum
figure;
STUDY.currentdesign = 1; % 1: level of terrains
[STUDY, specdata, specfreqs] = std_specplot(STUDY,ALLEEG,'clusters',3:10,'comps','all','subject','');
set(gcf,'position',[16 582 1340 751],'color','w')

figure;
STUDY.currentdesign = 2; % 1: level of speed trials?
std_specplot(STUDY,ALLEEG,'clusters',3:10,'comps','all','subject','','specparams',{'freqrange' freqrange 'specmode' 'fft' 'logtrials' 'off'});
set(gcf,'position',[16 582 1340 751],'color','w')

% ---- plot all dipoles but each cluster has different name
[STUDY] = std_dipplot(STUDY, ALLEEG, 'clusters','all','groups','on','mode','multicolor'); 
saveas(gcf,fullfile(save_study_folder,'all_cluster_scalpmap.fig'));

if showClusterPlotGUI==1
    [STUDY] = pop_clustedit(STUDY,ALLEEG); %pulls up cluster visualization interface
else
    clear EEG;  clear EEG_orig; clear ALLEEG; clear CURRENTSET; clear STUDY;
end

%}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    