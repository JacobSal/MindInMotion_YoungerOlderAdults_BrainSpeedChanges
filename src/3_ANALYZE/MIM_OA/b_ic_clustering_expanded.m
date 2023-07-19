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

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_b_ic_clustering_expanded.sh

%{
%## RESTORE MATLABs
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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_OA'];
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
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
RECOMPUTE_SPEC = false;
RECOMPUTE_ERSP = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; %false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
STAT_ALPHA = 0.05;
%- spec generation parameters
% FREQ_LIMITS = [1,70];
SPEC_MODE = 'psd'; %'fft'; %'psd'; 
%- spec plotting parameters
SPEC_FREQLIMITS = [3,70];
SPEC_NACCU = 2000;
% SPEC_XLIM = [1,100];
LOG_TRIALS = 'on';
SPEC_YLIM = [-30,-8];
SPEC_CONDSTATS = 'on';        % ['on'|'off]
SPEC_STAT_METHOD = 'perm';    % ['param'|'perm'|'bootstrap']
SPEC_ALPHA = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
SPEC_MCORRECT = 'fdr'; % 'cluster'; %fdr
% (07/16/2023) JS, updating to fdr as per CL YA paper
SPEC_GROUPSTATS = 'off';
SPEC_STAT_MODE = 'fieldtrip'; 
SPEC_STAT_FTMETHOD = 'permutation';
SPEC_SINGLETRIALS = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'
%- ERSP gen params
CYCLE_LIMITS = [3,0.8];
FREQ_FAC = 4;
%- ERSP PARAMS
% NOTE: see. statcondfieldtrip.m or std_stat.m
COND_EVENT_CHAR = 'cond';
% STAT_SINGLETRIALS = 'on'; %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
STAT_MODE = 'fieldtrip'; 
FIELDTRIP_METHOD = 'permutation'; %['permutation'|'parametric']
% TIMEWARP_NTIMES = [];
% ERSP_PAIRED = { 'off' 'off' };
ERSP_CAXIS = [-2,2]; %[-1.5,1.5];
ERSP_FREQLIMITS = [3,60];
ERSP_MCORRECT = 'fdr'; %'cluster'; %'fdr';
% (07/16/2023) JS, updating to fdr as per CL YA paper
ERSP_STAT_METHOD = 'bootstrap'; %'perm'; % ['param'|'perm'|'bootstrap']
% (07/16/2023) JS, updating to bootstrap as per CL YA paper
ERSP_ALPHA = 0.05; % [NaN|alpha], Significance threshold (0<alpha<<1)
ERSP_CONDSTATS = 'on'; % ['on'|'off]
ERSP_GROUPSTATS = 'off';
ERSP_SUBBASELINE = 'off'; %['on'|'off'];
ERSP_SINGLETRIALS = 'off'; %['on'|'off'], % this should be off for GAIT cycle analysis. 
ERSP_NACCU = 2000;
%- clustering parameters
% DO_SPEC_CALC = true;
MIN_ICS = 5;
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
%- iterative clustering parameters
n_iterations = 50;
outlier_sigma = 3;
%- datetime override
% dt = '05182023_MIM_OA_subset_N85_oldpipe';
% dt = '05192023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06122023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06282023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
dt = '07152023_MIM_OAN79_subset_prep_verified_gait';
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
end
% all_subjStr = {ALLEEG.subject};
%## POOL LOADING
POOL_SIZE = 5;
% fprintf('Loading STUDY & ALLEEG\n');
% if exist('SLURM_POOL_SIZE','var')
%     POOL_SIZE = min([SLURM_POOL_SIZE,length(ALLEEG)]);
% else
%     POOL_SIZE = 1;
% end
%% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
%     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
grandAvgWarpTo = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
% b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)-1000*(1/ALLEEG(1).srate)];
b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)];
% ERSP_CROP_TIMES=[grandAvgWarpTo(1)+abs(ALLEEG(1).etc.epoch.epoch_limits(1))*1000, grandAvgWarpTo(5)];
ERSP_CROP_TIMES=[grandAvgWarpTo(1), grandAvgWarpTo(5)];
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',b_lims(1),b_lims(2));
disp(grandAvgWarpTo);
%% (SET PARAMS)
STUDY = pop_statparams(STUDY,'condstats',ERSP_CONDSTATS,...
        'groupstats',ERSP_GROUPSTATS,...
        'method',ERSP_STAT_METHOD,...
        'singletrials',ERSP_SINGLETRIALS,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
        'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_SUBBASELINE,...
      'ersplim',ERSP_CAXIS,'freqrange',ERSP_FREQLIMITS);
tmp_group_orig = cell(length(ALLEEG),1);
tmp_group_unif = cell(length(ALLEEG),1);
for subj_i = 1:length(ALLEEG)
    tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
    tmp_group_unif{subj_i} = 'Older Adults';
end
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',SPEC_FREQLIMITS,'logtrials',LOG_TRIALS});
    end
end
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
if RECOMPUTE_ERSP
    disp(['Grand average (across all subj) warp to: ',num2str(grandAvgWarpTo)]);
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        %- determine timewarping parameters
         if DO_TIMEWARP
            timewarp_param = EEG.timewarp.latencies;
            timewarpms_param = grandAvgWarpTo;
         else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',b_lims,...
                    'commonbase','on','trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
end
EEG = [];
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
%- store essential info in STUDY struct for later reading
freqrange = [3 45];
STUDY.etc.clustering.preclustparams.clustering_weights = clustering_weights;
STUDY.etc.clustering.preclustparams.freqrange = freqrange;
%% (STEP 1) CALCULATE CLUSTER SOLUTIONS
%- NOTE: partly adapt from bemobil_repeated_clustering
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
% mean_IC_allSub = 6;
%## Cluster Components
%- (Step 1) find the optimal number of clusters without generating outliers
%- (07/12/2023) CL, 30 is the maximum number clusters typically
    % used for EEG (see Makyoto) % I think 5 is a reasonable lower bound
    % updated to do repeated clustering clust_CL is a function
    % that is similar as pop_clust that does k_means cluster 
    % but not generate cluster struct in the STUDY; clust_CL 
    % used kmeans to cluster for 1000 times
p = 1;
ClusterOutcome = cell(mean_IC_allSub-MIN_ICS,1);
cluster_idx = zeros(length([STUDY.datasetinfo.comps]),mean_IC_allSub);
for clust_i = MIN_ICS:mean_IC_allSub 
    [~, ClusterOutcome{clust_i}] = clust_CL(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  clust_i , 'outliers',  Inf );
    cluster_idx(1:length(ClusterOutcome{clust_i}.IDX),p) = ClusterOutcome{clust_i}.IDX;
    p = p+1;
end
cluster_idx = cluster_idx(:,all(cluster_idx ~= 0));
%- Calculate evaulation criteria for K (kmeans alg)
eva1 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'silhouette'); % this is the method to find optimal number of clusters (confirmed by EEGlab mailing list)
eva2 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'CalinskiHarabasz');
eva3 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'DaviesBouldin');
%- Calculate evaulation criteria
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
saveas(gcf, [save_dir filesep 'cluster_kmeans_eval.fig'])
saveas(gcf, [save_dir filesep 'cluster_kmeans_eval.jpg'])
%% (Step 2) CALCULATE REPEATED CLUSTERED SOLUTIONS
%- NOTE: the clustering solutions are not exactly the same but not super different
parfor (clust_i = MIN_ICS:mean_IC_allSub,POOL_SIZE)
    clustering_solutions = repeated_clustering(STUDY,ALLEEG, n_iterations, clust_i, outlier_sigma, STUDY.etc.clustering.preclustparams);
    cluster_dir = [save_dir filesep clustering_method filesep num2str(clust_i)];
    if ~exist(cluster_dir,'dir')
        mkdir(cluster_dir);
    end
    par_save(clustering_solutions,cluster_dir,sprintf('clustering_solutions_%i.mat',clust_i));
    %## Use RV to Choose Components
    [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
    %## Look up cluster centroid Brodmann area
    atlas = ft_read_atlas([FIELDTRIP_PATH filesep 'template' filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii']);
    atlas_name = cell(size(cluster_update,1),2);
    for i = 3:length(cluster_update)
        cfg            = [];
        cfg.roi        = cluster_update(i).dipole.posxyz;
        cfg.output     = 'multiple';
        cfg.atlas      = atlas;
        cfg.inputcoord = 'mni';
        cfg.sphere = 1;
        labels = ft_volumelookup(cfg, atlas);
        [~, indx] = max(labels.count);
        cluster_update(i).analabel = labels.name(indx);
        if ~isempty(labels)
            [val, indx] = max(labels.count);
            if strcmp(labels.name(indx),'no_label_found')
                sub_indx = find(labels.count ~= 0 & labels.count < val);
                if isempty(sub_indx)
                    atlas_name{i,1} = ['CLs ',num2str(i)];
                    atlas_name{i,2} = labels.name{indx};
                    continue;
                end
                atlas_name{i,1} = ['CLs ',num2str(i)];
                tmp = sprintf('(N=%i) %s',labels.count(sub_indx(1)),labels.name{sub_indx(1)});
                for ic_i = sub_indx(2:end)
                    tmp = [tmp,' & ', sprintf('(N=%i) %s',labels.count(ic_i),labels.name{ic_i})];
                end
                atlas_name{i,2} = tmp;
            else
                atlas_name{i,1} = ['CLs ',num2str(i)];
                atlas_name{i,2} = labels.name{indx};
            end
        else
            warning('No label found for cluster %i',i);
        end
    end
    fprintf('Cluster \t\t Label\n');
    fprintf('________________________\n');
    for i = 3:length(cluster_update)
        label = cellstr(atlas_name{i,2});
        cl =  cellstr(atlas_name{i,1});
        fprintf('%s\t\t%s\n',cl{:},label{:})
    end
    par_save(cluster_update,cluster_dir,sprintf('cluster_update_%i.mat',clust_i));
end

%% (STEP 3) PLOT ERSPs, DIPOLEs, & PSDs
%## LOCAL SWITCHES
% ERSP_SINGLETRIALS = 'on';
XAXIS_LABEL = 'ms';
COLORAXIS_LABEL = 'dB';
ERSP_FREQSCALE = 'log'; % 'native', 'LOG'
ERSP_CHANLOCS = struct('labels', {});
RUN_CLUSTERS = (MIN_ICS:mean_IC_allSub);
RUN_CLUSTERS = RUN_CLUSTERS(RUN_CLUSTERS>=12);
for clust_i = RUN_CLUSTERS
    cluster_dir = [save_dir filesep clustering_method filesep num2str(clust_i)];
    cluster_update = par_load(cluster_dir,sprintf('cluster_update_%i.mat',clust_i));
    STUDY.cluster = cluster_update;
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cluster] = eeglab_get_cluster_comps(STUDY);
    fprintf('Clusters with more than 50%% of subjects:'); fprintf('%i,',valid_cluster(1:end-1)); fprintf('%i',valid_cluster(end)); fprintf('\n');
    fprintf('Main cluster numbers:'); fprintf('%i,',main_cl_inds(1:end-1)); fprintf('%i',main_cl_inds(end)); fprintf('\n');
    %## Update 7/20/2022 - Plot clustering results for AHA proposal
    % Plot scalp topographs which also need to be averaged? 
    if ~isfield(STUDY.cluster,'topo') 
        STUDY.cluster(1).topo = [];
    end
    for i = main_cl_inds % For each cluster requested
        if isempty(STUDY.cluster(i).topo)
            % Using this custom modified code to allow taking average within participant for each cluster
            STUDY = std_readtopoclust_CL(STUDY,ALLEEG, i);
        end
    end
    %% ASSIGN GROUP NAME & STUDY STATS
    %- get non-outlier clusters
%     [~,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
    %- Plot topographies 
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_unif{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
    end
    figure;
    std_topoplot_CL(STUDY,valid_cluster,'together');
    set(gcf,'position',[16 582 1340 751],'color','w')
%     saveas(gcf,fullfile(cluster_dir,'Cluster_topo_plot_averaged.fig'));
    saveas(gcf,fullfile(cluster_dir,'Cluster_topo_plot_averaged.jpg'));
    %- Plot all dipole fit locations
    std_dipplot(STUDY,ALLEEG,'clusters','all','figure','off');
    set(gcf,'position',[16 582 1340 751],'color','w')
    %- Plot dipole fit locations after averaged within participants
    std_dipplot_CL(STUDY,ALLEEG,'clusters',valid_cluster,'figure','off','mode','together_averaged');
    set(gcf,'position',[16 582 1340 751],'color','w')
%     saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_averaged.fig'));
    saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_averaged.jpg'));
    %- Plot dipole clusters 
    STUDY.etc.dipparams.centrline = 'off';
    std_dipplot_CL(STUDY,ALLEEG,'clusters',valid_cluster,'figure','off','mode','together_averaged_only','spheres','off','projlines','off');
    set(gcf,'position',[16 582 500 300],'color','w')
    camzoom(1.2^2);
%     saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_allaveraged_only.fig'));
    saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_allaveraged_only.jpg'));
    %- Plot dipole fit locations after averaged within participants and
    % different clusters have different colors
    STUDY.etc.dipparams.centrline = 'off';
    std_dipplot_CL(STUDY,ALLEEG,'clusters',valid_cluster,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
    set(gcf,'position',[16 582 500 300],'color','w')
    camzoom(1.2^2);
%     saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_allaveraged.fig'));
    saveas(gcf,fullfile(cluster_dir,'Cluster_dipole_plot_allaveraged.jpg'));
    %- Spec plot
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_orig{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
    end
    %     [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
    %             'subjselect', {ALLEEG.subject},...
    %             'variable1','group',...
    %             'values1',unique({ALLEEG.group}),...
    %             'variable2',COND_EVENT_CHAR,...
    %             'values2', COND_DESIGNS{des_i});
    [STUDY] = std_makedesign(STUDY, ALLEEG, 1,...
            'subjselect', {ALLEEG.subject},...
            'variable1','group',...
            'values1',unique({ALLEEG.group}));
    fprintf('==== Making Spectogram Plots ====\n');
    std_specplot(STUDY,ALLEEG,'clusters',valid_cluster,...
        'freqrange',SPEC_FREQLIMITS);
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1480 920];
%     saveas(fig_i,[cluster_dir filesep sprintf('allSpecPlot.fig')]);
    saveas(fig_i,[cluster_dir filesep sprintf('allSpecPlot.jpg')]);
    %- close all figures
    close all
    %## ersp plot per cluster per condition
    %- params (06/10/2023) JS, removing timerange loading, consolidating single
    %trials stats variable
    %% ASSIGN GROUP NAME & STUDY STATS
    STUDY = pop_statparams(STUDY,'condstats',ERSP_CONDSTATS,...
            'groupstats',ERSP_GROUPSTATS,...
            'method',ERSP_STAT_METHOD,...
            'singletrials',ERSP_SINGLETRIALS,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
            'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
    STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_SUBBASELINE,...
          'ersplim',ERSP_CAXIS,'freqrange',ERSP_FREQLIMITS);
    for des_i = 1:length(COND_DESIGNS)
        fprintf('==== Making Study Design ====\n');
        %## group & condition study
        %{
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_orig{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1','group',...
            'values1',unique({ALLEEG.group}),...
            'variable2',COND_EVENT_CHAR,...
            'values2', COND_DESIGNS{des_i});
        %}
        %## group only differences
        %{
        group_only_flag = true;
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_orig{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1','group',...
            'values1',unique({ALLEEG.group}));
        %}
        %## conditions across all older subjects
        
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_unif{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1',COND_EVENT_CHAR,...
            'values1', COND_DESIGNS{des_i});
       
        parfor (j = 1:length(valid_cluster),POOL_SIZE)
            cluster_i = valid_cluster(j);
            %% (ERSP PLOT 2) CUSTOM
            %- load ERSP data using std_readdat
            %* (06/04/2023) JS, could try and load a subject at a time,
            % baseline, then stack them together to plot)
            %* (06/10/2023) JS, removing timerange loading, need to custom
            %prune later
            [~, ersp_data, ersp_times, ersp_freqs, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
                            'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
                            'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
                            'design',des_i);
            %- across condition plot
            titles_ersp = std_figtitle('threshold', STAT_ALPHA, 'mcorrect', ERSP_MCORRECT,...
                                    'condstat', ERSP_CONDSTATS,...
                                    'statistics', STAT_MODE,... %STAT_MODE,... %'fieldtrip cluster',... %ERSP_STAT_METHOD,...
                                    'condnames', COND_DESIGNS{des_i},...
                                    'clustname', STUDY.cluster(cluster_i).name,...
                                    'subject', [], 'datatype', 'ersp', 'plotmode', 'normal', ...
                                    'effect', 'main');

            allersp = ersp_data;
            alltimes = ersp_times;
            allfreqs = ersp_freqs;
            %## subtract condition based mean ersp
            % (05/31/2023): JS I don't believe this to generate a good result.
            % Problem with logging the output after baseline subtraction.
            % (06/02/2023): JS I'm trying this again, going to try and extract
            % each epoch using event structure and sort into conditions.
%             stacked_cond_ersp = [allersp{:}];
            %## Cropping Indicies
            crop_inds = (alltimes>=ERSP_CROP_TIMES(1) & alltimes<=ERSP_CROP_TIMES(2));
            alltimes = alltimes(crop_inds);
            for cond_i = 1:length(allersp)
                allersp{cond_i} = allersp{cond_i}(:,crop_inds,:);
            end
            %## Across Design Baselineing
            %{
            tmp_allersp = [allersp{:}];
            tmp_allersp = reshape(tmp_allersp,size(allersp{1},1),size(allersp{1},2),size(allersp{1},3)*length(allersp));
            tmp = mean(tmp_allersp,3);
            tmp = mean(tmp,2);
            tmp_std = std(tmp,[],2);
            tmp = repmat(tmp,1,size(tmp_allersp,2)); 
            tmp_std = repmat(tmp_std,1,size(tmp_allersp,2));
            for cond_i = 1:length(allersp)
                allersp{cond_i} = allersp{cond_i} - tmp;
            end
            %}
            %## Within Each Condition Baselineing
            for cond_i = 1:length(allersp)
                %- interesting way to baseline using mean across trials and
                %time
                tmp = mean(allersp{cond_i},3);
    %             tmp_std = std(allersp{cond_i},[],3);
                tmp = mean(tmp,2);
                tmp_std = std(tmp,[],2);
                tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
                tmp_std = repmat(tmp_std,1,size(allersp{cond_i},2));
                %- interesting way to baseline using median across trials and
                %time
    %             tmp = median(allersp{cond_i},3);
    %             tmp = median(tmp,2);
    %             tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
                %- interesting way to baseline using median across trials and
                %time
    %             tmp = mean(allersp{cond_i},3);
    %             tmp = median(tmp,2);
    %             tmp = repmat(tmp,1,size(allersp{cond_i},2));
                %- subtract out baseline from all trials.
                allersp{cond_i} = allersp{cond_i} - tmp;
            end
            %## MULTI SUBJECT ANALYSIS
            %{
            [pcond_ersp,pgroup_ersp,pinter_ersp,stat_cond,stat_group,stat_inter] = std_stat(allersp,...
                    'condstats', ERSP_CONDSTATS,...
                    'groupstats',ERSP_GROUPSTATS,...
                    'method',ERSP_STAT_METHOD,...
                    'naccu',ERSP_NACCU,...
                    'alpha',ERSP_ALPHA,...
                    'mcorrect',ERSP_MCORRECT,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
                    'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
            %- plot
            std_plottf(alltimes, allfreqs, allersp,...
                           'datatype', 'ersp', ...
                           'groupstats', pgroup_ersp,...
                           'condstats', pcond_ersp,...
                           'interstats', pinter_ersp,...
                           'plotmode', 'normal',...
                           'titles', titles_ersp,...
                           'caxis',ERSP_CAXIS,...
                           'freqscale',ERSP_FREQSCALE,...
                           'ersplim',[],...
                           'effect','main',...
                           'threshold',ERSP_ALPHA,...
                           'maskdata','off','averagemode','ave');
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'Position',[100 100 1620 300]);
            for i = [2,4,5,6]
                tmp = fig_i.Children(i);
                set(tmp,'YTick',log([4.01,8,13,30,50,99.4843])); 
                set(tmp,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
                ylabel(tmp,sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
                xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
                set(tmp,'XTick',grandAvgWarpTo)
                xx = num2cell(grandAvgWarpTo);
                xx = cellfun(@num2str,xx,'UniformOutput',false);
                set(tmp,'XTickLabel',{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'});
                xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
                xline(tmp,grandAvgWarpTo(1),'k--');
                xline(tmp,grandAvgWarpTo(2),'k--');
                xline(tmp,grandAvgWarpTo(3),'k--');
                xline(tmp,grandAvgWarpTo(4),'k--');
                xline(tmp,grandAvgWarpTo(5),'k--');
            end
            fig_i.Name = sprintf('Cluster %i) ERSP averages across Conditions',cluster_i);
%             saveas(fig_i,[cluster_dir filesep sprintf('ersp_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
            saveas(fig_i,[cluster_dir filesep sprintf('ersp_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
            close(fig_i)
            %}
            %## SINGLE SUBJECT LOOP
            for subj_i = 1:size(allersp{1},1)
                tmp_ersp = cell(size(allersp));
                for cond_i = 1:length(allersp)
                    tmp_ersp{cond_i} = allersp{cond_i}(:,:,subj_i);
                end
                %- plot
                std_plottf(alltimes, allfreqs, tmp_ersp,...
                               'datatype', 'ersp', ...
                               'plotmode', 'normal',...
                               'titles', titles_ersp,...
                               'caxis',ERSP_CAXIS,...
                               'freqscale',ERSP_FREQSCALE,...
                               'ersplim',[],...
                               'maskdata','off','averagemode','rms');
                fig_i = get(groot,'CurrentFigure');
                set(fig_i,'Position',[100 100 1620 300]);
                for i = [2,4,5,6]
                    tmp = fig_i.Children(i);
                    set(tmp,'YTick',log([4.01,8,13,30,50,99.4843])); 
                    set(tmp,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
                    ylabel(tmp,sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
                    xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
                    set(tmp,'XTick',grandAvgWarpTo)
                    xx = num2cell(grandAvgWarpTo);
                    xx = cellfun(@num2str,xx,'UniformOutput',false);
                    set(tmp,'XTickLabel',{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'});
                    xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
                    xline(tmp,grandAvgWarpTo(1),'k--');
                    xline(tmp,grandAvgWarpTo(2),'k--');
                    xline(tmp,grandAvgWarpTo(3),'k--');
                    xline(tmp,grandAvgWarpTo(4),'k--');
                    xline(tmp,grandAvgWarpTo(5),'k--');
                end
                fig_i.Name = sprintf('Cluster %i) ERSP averages across Conditions',cluster_i);
    %             saveas(fig_i,[cluster_dir filesep sprintf('ersp_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
                saveas(fig_i,[cluster_dir filesep sprintf('%s_ersp_%i_%s.jpg',ALLEEG(subj_i).subject,cluster_i,[COND_DESIGNS{des_i}{:}])]);
                close(fig_i)
            end
            
        end
        if group_only_flag
            break
        end
    end
    STUDY = pop_statparams(STUDY, 'condstats', SPEC_CONDSTATS,...
                        'method',SPEC_STAT_METHOD,...
                        'singletrials',SPEC_SINGLETRIALS,'mode',SPEC_STAT_MODE,'fieldtripalpha',SPEC_ALPHA,...
                        'fieldtripmethod',SPEC_STAT_FTMETHOD,'fieldtripmcorrect',SPEC_MCORRECT,'fieldtripnaccu',SPEC_NACCU);
    %## MAKEDESIGN & PLOT
    for des_i = 1:length(COND_DESIGNS)
        %-%## group & condition study
        %{
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_orig{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1','group',...
            'values1',unique({ALLEEG.group}),...
            'variable2',COND_EVENT_CHAR,...
            'values2', COND_DESIGNS{des_i});
        %}
        %## group only differences
        group_only_flag = true;
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_orig{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1','group',...
            'values1',unique({ALLEEG.group}));
        %## conditions across all older subjects
        %{
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_unif{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
        end
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1',COND_EVENT_CHAR,...
            'values1', COND_DESIGNS{des_i});
        %}
        parfor (j = 1:length(valid_cluster),POOL_SIZE)
            cluster_i = valid_cluster(j);
            %-
        %     [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
        %             'clusters',cluster_i,'comps','all','subject','','freqrange', FREQ_LIMITS,'subtractsubjectmean','on');
            %-
            [~, ~, ~, ~, ~, ~] = std_specplot(STUDY, ALLEEG,...
                    'clusters',cluster_i,'freqrange',SPEC_FREQLIMITS,'plotmode','condensed',...
                    'plotconditions','together','ylim',SPEC_YLIM); %'plotgroups','together'
            fig_i = get(groot,'CurrentFigure');
            %- set figure line colors
            cc = linspecer(length(COND_DESIGNS{des_i}));
            iter = 1;
            for i = 1:length(fig_i.Children(2).Children)
                set(fig_i.Children(2).Children(i),'LineWidth',1.5);
                set(fig_i.Children(2).Children(i),'Color',horzcat(cc(iter,:),0.6));
                if iter == size(cc,1)
                    iter = 1;
                else
                    iter = iter + 1;
                end                
            end
            %- tight layout
            % (06/22/2023) JS, seems like changing the Position property also
            % erases the statistics bar. Not sure how to fix?
    %         set(fig_i.Children(2),'Position',[0.26,0.26,0.54,0.51]); %DEFAULT
    %         set(fig_i.Children(3),'Position',[0.26,0.2345,0.54,0.0255]); %DEFAULT
            set(fig_i.Children(2),'Position',[0.15,0.15,0.8,0.8]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
            set(fig_i.Children(3),'Position',[0.15,0.15-0.0255,0.8,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
            set(fig_i.Children(1),'Location','southeast') %reset Legend
            drawnow;
            %- save
%             saveas(fig_i,[cluster_dir filesep sprintf('powerspec_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
            saveas(fig_i,[cluster_dir filesep sprintf('powerspec_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
            close(fig_i)
        end
        if group_only_flag
            break
        end
    end
end