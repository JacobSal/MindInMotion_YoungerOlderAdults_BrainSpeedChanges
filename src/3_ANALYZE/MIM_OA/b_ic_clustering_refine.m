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

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_b_ic_clustering_refine.sh

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
FORCE_RECALC_SPEC = true;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result? 
SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
% NOTE: see. statcondfieldtrip.m or std_stat.m
COND_EVENT_CHAR = 'cond';
%- clustering parameters
MIN_ICS_SUBJ = [5,6,7]; %[2,3,4,5,6,7,8]; % iterative clustering
% K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
CLUSTER_SWEEP_VALS = [14,19]; %[10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
% DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = true;
% DO_K_SWEEPING = false;
% (08/21/2023) JS, this currenlty doesn't do anything but take up more
% memory.
REPEATED_CLUSTERING_STD = 3;
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',[],...
    'filepath',[],...
    'outliers',inf(),...
    'maxiter',200);
%- 
%- custom params
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = 'dipole_1'; %['dipole_',num2str(clustering_weights.dipoles),...
%     '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
%     '_spec_',num2str(clustering_weights.spec)];
STD_PRECLUST_COMMAND = {'dipoles','weight',clustering_weights.dipoles};
%- iterative clustering parameters
% n_iterations = 50;
% outlier_sigma = 3;
%- datetime override
% dt = '05182023_MIM_OA_subset_N85_oldpipe';
% dt = '05192023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06122023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06282023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07152023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
dt = '10052023_MIM_OAN70_noslowwalkers';
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
CLUSTER_PARAMS.filename = STUDY.filename;
CLUSTER_PARAMS.filepath = STUDY.filepath;
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
% STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
%         'groupstats',ERSP_STAT_PARAMS.groupstats,...
%         'method',ERSP_STAT_PARAMS.method,...
%         'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
%         'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
tmp_group_orig = cell(length(ALLEEG),1);
tmp_group_unif = cell(length(ALLEEG),1);
for subj_i = 1:length(ALLEEG)
    tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
    tmp_group_unif{subj_i} = 'Older Adults';
end
%-
%{
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = tmp_group_orig{subj_i};
    STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
end
%}
%- NOTE: partly adapt from bemobil_repeated_clustering
numIC = zeros(length(STUDY.datasetinfo),1);
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
tmp = strsplit(ALLEEG(1).filename,'.');
spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
    fprintf('Calculating Spectograms...\n');
    %- override variables for the stats
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_orig{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        for in_i = 1:length(STUDY.datasetinfo(subj_i).trialinfo)
            STUDY.datasetinfo(subj_i).trialinfo(in_i).group = tmp_group_orig{subj_i};
        end
    end
    parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/2))
        EEG = ALLEEG(subj_i);
        TMP_STUDY = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %- overrride datasetinfo to trick std_precomp to run.
        TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
        fprintf('SUBJECT: %s\n',TMP_STUDY.datasetinfo.subject);
        fprintf('GROUP: %s\n',TMP_STUDY.datasetinfo.group);
        disp(STUDY.datasetinfo(subj_i));
        TMP_STUDY.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(TMP_STUDY, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_PARAMS.specmode,'freqfac',SPEC_PARAMS.freqfac,...
                    'freqrange',SPEC_PARAMS.freqrange,'logtrials',SPEC_PARAMS.logtrials});
    end
end
%% (STEP 4) CALCULATE CLUSTER SOLUTIONS AFTER PRUNING SUBJECTS WITH LOW NUMBER OF ICs
%## GET EEGLAB PATH
tmp = strsplit(path,';');
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
%-
STUDY_COND_DESI = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
%%
if DO_K_ICPRUNE 
%     for i = 1:length(MIN_ICS_SUBJ)
    parfor (i = 1:length(MIN_ICS_SUBJ),length(MIN_ICS_SUBJ))
        tmp_dir = [save_dir filesep sprintf('icrej_%i',MIN_ICS_SUBJ(i))];
        if ~exist(tmp_dir,'dir')
            mkdir(tmp_dir)
        end
        ics_subjs = cellfun(@(x) size(x,2),{ALLEEG.icawinv});
        idx = ics_subjs > MIN_ICS_SUBJ(i);
        TMP_ALLEEG = ALLEEG(idx);
        %##
        [TMP_STUDY, TMP_ALLEEG] = std_editset([],TMP_ALLEEG,...
                                        'updatedat','off',...
                                        'savedat','off',...
                                        'name',sprintf('temp_study_rejics%i',MIN_ICS_SUBJ(i)),...
                                        'filename',sprintf('temp_study_rejics%i.study',MIN_ICS_SUBJ(i)),...
                                        'filepath',tmp_dir);
        [TMP_STUDY,TMP_ALLEEG] = std_checkset(TMP_STUDY,TMP_ALLEEG);
        [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                            TMP_STUDY.filename,TMP_STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
        %##
        [TMP_STUDY,TMP_ALLEEG] = std_preclust(TMP_STUDY,TMP_ALLEEG,1,STD_PRECLUST_COMMAND);
        %- remove components inside the brain
        [TMP_STUDY, TMP_ALLEEG] = std_editset(TMP_STUDY,TMP_ALLEEG,...
                                        'commands',{'inbrain','on'});                  
        %- store essential info in STUDY struct for later reading
        freqrange = [4 60];
        TMP_STUDY.etc.clustering.preclustparams.clustering_weights = clustering_weights;
        TMP_STUDY.etc.clustering.preclustparams.freqrange = freqrange;
        all_solutions = mim_cluster_process(TMP_STUDY,TMP_ALLEEG,tmp_dir,...
            'CLUSTER_PARAMS',CLUSTER_PARAMS,...
            'REPEATED_CLUSTERING_STD',REPEATED_CLUSTERING_STD,...
            'MAX_REPEATED_ITERATIONS',MAX_REPEATED_ITERATIONS,...
            'KS_TO_TEST',CLUSTER_SWEEP_VALS);
        %## (Step 2) CALCULATE REPEATED CLUSTERED SOLUTIONS
        %- NOTE: the clustering solutions are not exactly the same but not super different
    %     cluster_ks = K_RANGE(1):K_RANGE(2);
    %     parfor (j = 1:length(cluster_ks),length(cluster_ks))
        TMP_STUDY.etc.cluster_vars = [];
        for j = 1:length(CLUSTER_SWEEP_VALS)
            %-
            clust_i = CLUSTER_SWEEP_VALS(j);
            cluster_dir = [tmp_dir filesep num2str(clust_i)];
            if ~exist(cluster_dir,'dir')
                mkdir(cluster_dir)
            end
            %## Use RV to Choose Components
            TMP_STUDY.cluster = all_solutions{j}.solutions{1};
        %     [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
            [TMP_STUDY,~,~] = cluster_rv_reduce(TMP_STUDY,TMP_ALLEEG);
            cluster_update = TMP_STUDY.cluster;
            %- get cluster centroid and residual variance
            cluster_update = cluster_comp_dipole(TMP_ALLEEG, cluster_update);
            TMP_STUDY.cluster = cluster_update;
            %## Look up cluster centroid Brodmann area
            [~,atlas_names,~] = add_anatomical_labels(TMP_STUDY,TMP_ALLEEG);
            for k = 1:length(cluster_update)
                cluster_update(k).analabel = atlas_names{k,2};
            end
            par_save(cluster_update,cluster_dir,sprintf('cl_inf_%i.mat',clust_i));
            TMP_STUDY.etc.cluster_vars(j).fpath = [cluster_dir filesep sprintf('cl_inf_%i.mat',clust_i)];
            TMP_STUDY.etc.cluster_vars(j).inf = {'type','kmeans','k',clust_i,'params',CLUSTER_PARAMS};
        end
        %## (PLOT) Looking for dipoles outside of brain.
        %{
        vol = load(MNI_VOL);
        vol = vol.vol;
        vol_bounds = zeros(3,2,length(vol.bnd));
        rmv_dip = [];
        fig = figure('color','w');
%         ft_plot_mesh(vol.bnd(1));
%         ft_plot_mesh(vol.bnd(2));
        ft_plot_mesh(vol.bnd(3));
        hold on;
        for c = 2:length(TMP_STUDY.cluster)
            for d = 1:size(TMP_STUDY.cluster(c).all_diplocs,1)
                for b = 1:length(vol.bnd)
                    dist = vol.bnd(b).pnt-repmat(TMP_STUDY.cluster(c).all_diplocs(d,:),size(vol.bnd(b).pnt,1),1);
                    sig = sign(dist);
                    pos = sig > 0;
                    neg = sig < 0;
                    depth = ft_sourcedepth(TMP_STUDY.cluster(c).all_diplocs(d,:), vol);
%                     disp(depth);
                    if depth > 0 && b%any(all(pos,1)) || any(all(neg,1))
                        rmv_dip = [rmv_dip;[c,d,b]];
                        plot3(TMP_STUDY.cluster(c).all_diplocs(d,1),TMP_STUDY.cluster(c).all_diplocs(d,2),TMP_STUDY.cluster(c).all_diplocs(d,3),'*-');
                    end
                end
            end
        end
        hold off;
        drawnow;
        saveas(fig,[tmp_dir filesep sprintf('ics_out_of_brain.fig')]);
        %}
        %- save
        [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                            TMP_STUDY.filename,TMP_STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
    end
end
%%
if DO_K_ICPRUNE
    parfor (i = 1:length(MIN_ICS_SUBJ),length(MIN_ICS_SUBJ))
%     for i = 1:length(MIN_ICS_SUBJ)
        study_fName = sprintf('temp_study_rejics%i',MIN_ICS_SUBJ(i));
        tmp_dir = [save_dir filesep sprintf('icrej_%i',MIN_ICS_SUBJ(i))];
        if ~exist(tmp_dir,'dir')
            mkdir(tmp_dir)
        end
        if ~exist([tmp_dir filesep study_fName '.study'],'file')
            error('ERROR. study file does not exist');
        else
            if ~ispc
                [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',tmp_dir);
            else
                [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',tmp_dir);
            end
        end
        STUDY_GROUP_DESI = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',unique({TMP_STUDY.datasetinfo.group})},...
            {'subjselect',{},...
            'variable2','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',unique({TMP_STUDY.datasetinfo.group})}};
        %## grab subjects for study designs
        tmp_group_orig = cell(length(TMP_ALLEEG),1);
        tmp_group_unif = cell(length(TMP_ALLEEG),1);
        for subj_i = 1:length(TMP_ALLEEG)
            tmp_group_orig{subj_i} = TMP_ALLEEG(subj_i).group;
            tmp_group_unif{subj_i} = 'Older Adults';
        end
        %## ersp plot per cluster per condition
        TMP_STUDY = pop_statparams(TMP_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
                'groupstats',ERSP_STAT_PARAMS.groupstats,...
                'method',ERSP_STAT_PARAMS.method,...
                'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
                'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
                'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
                'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
        TMP_STUDY = pop_specparams(TMP_STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
                'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
                'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
        for j = 1:length(CLUSTER_SWEEP_VALS)
            %-
            clust_i = CLUSTER_SWEEP_VALS(j);
            cluster_dir = [tmp_dir filesep num2str(clust_i)];
            cluster_update = par_load(TMP_STUDY.etc.cluster_vars(j).fpath,[]); %par_load(cluster_dir,sprintf('cluster_inf_%i.mat',clust_i));
            %## Create Plots
            TMP_STUDY.cluster = cluster_update;
             %- get inds
            [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(TMP_STUDY);
            %- clusters to plot
            CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
            %## PLOT spectrum based information
            %{
            %- ASSIGN DESIGNS
            %* combine groups
            for subj_i = 1:length(TMP_ALLEEG)
                TMP_ALLEEG(subj_i).group = tmp_group_unif{subj_i};
                TMP_STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
            end
            %* assign studies
            TMP_STUDY.cache = [];
            for des_i = 1:length(STUDY_COND_DESI)
                [TMP_STUDY] = std_makedesign(TMP_STUDY,TMP_ALLEEG,des_i,STUDY_COND_DESI{des_i}{:});
            end
            %- SPECTRUM PLOTS
            mim_custom_spec_plots(TMP_STUDY,TMP_ALLEEG,cluster_dir,...
                'CLUSTERS_TO_PLOT',CLUSTER_PICKS);
            %}
            %- ASSIGN DESIGNS
            
            %* combine groups
            for subj_i = 1:length(TMP_ALLEEG)
                TMP_ALLEEG(subj_i).group = tmp_group_orig{subj_i};
                TMP_STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
            end
            %* assign studies
            TMP_STUDY.cache = [];
            for des_i = 1:length(STUDY_GROUP_DESI)
                [TMP_STUDY] = std_makedesign(TMP_STUDY,TMP_ALLEEG,des_i,STUDY_GROUP_DESI{des_i}{:});
            end
            %- SPECTRUM PLOTS
            mim_custom_spec_plots(TMP_STUDY,TMP_ALLEEG,cluster_dir,...
                'CLUSTERS_TO_PLOT',CLUSTER_PICKS);
            %## PLOT CLUSTER BASE INFORMATION
            %- CLUSTER DIPOLES, TOPOS
            mim_gen_cluster_figs(TMP_STUDY,TMP_ALLEEG,cluster_dir,...
                'CLUSTERS_TO_PLOT',CLUSTER_PICKS);
            %- close all figures
            close all
        end
    end
end
%%
%{
parfor (i = 1:length(MIN_ICS_SUBJ),length(MIN_ICS_SUBJ))
    study_fName = sprintf('temp_study_rejics%i',MIN_ICS_SUBJ(i));
    tmp_dir = [save_dir filesep sprintf('icrej_%i',MIN_ICS_SUBJ(i))];
    if ~exist(tmp_dir,'dir')
        mkdir(tmp_dir)
    end
    if ~exist([tmp_dir filesep study_fName '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',tmp_dir);
        else
            [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',tmp_dir);
        end
    end
    %## grab subjects for study designs
    tmp_group_orig = cell(length(TMP_ALLEEG),1);
    tmp_group_unif = cell(length(TMP_ALLEEG),1);
    for subj_i = 1:length(TMP_ALLEEG)
        tmp_group_orig{subj_i} = TMP_ALLEEG(subj_i).group;
        tmp_group_unif{subj_i} = 'Older Adults';
    end
    %## ersp plot per cluster per condition
    TMP_STUDY = pop_statparams(TMP_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
            'groupstats',ERSP_STAT_PARAMS.groupstats,...
            'method',ERSP_STAT_PARAMS.method,...
            'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
            'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
            'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
            'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    TMP_STUDY = pop_specparams(TMP_STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
            'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
            'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
    %##
    for j = 1:length(CLUSTER_SWEEP_VALS)
        %-
        clust_i = CLUSTER_SWEEP_VALS(j);
        cluster_dir = [tmp_dir filesep num2str(clust_i)];
        cluster_update = par_load(TMP_STUDY.etc.cluster_vars(j).fpath,[]); %par_load(cluster_dir,sprintf('cluster_inf_%i.mat',clust_i));
        %## Create Plots
        TMP_STUDY.cluster = cluster_update;
        [comps_out,main_cl_inds,~,valid_cluster,~,nonzero_cluster] = eeglab_get_cluster_comps(TMP_STUDY);
        CLUSTER_PICKS = nonzero_cluster; %valid_cluster; %main_cl_inds(2:end);
        fprintf('Clusters with more than 50%% of subjects:'); fprintf('%i,',valid_cluster(1:end-1)); fprintf('%i',valid_cluster(end)); fprintf('\n');
        fprintf('Main cluster numbers:'); fprintf('%i,',main_cl_inds(1:end-1)); fprintf('%i',main_cl_inds(end)); fprintf('\n');
        cl_names = {TMP_STUDY.cluster(CLUSTER_PICKS).name};
        for k = 1:length(CLUSTER_PICKS)
            cl_to_plot = CLUSTER_PICKS(k);
            fprintf('\n(k=%i, Cluster=%i) Plotting for Min ICs Per Subject of of Cluster %i\n',clust_i,cl_to_plot,MIN_ICS_SUBJ(i));
            %-
            % Plot scalp topographs which also need to be averaged? 
            if ~isfield(TMP_STUDY.cluster,'topo') 
                TMP_STUDY.cluster(1).topo = [];
            end
            for c = 1:length(cl_to_plot) % For each cluster requested
                clus_i = cl_to_plot(c);
                disp(clus_i)
                if isempty(TMP_STUDY.cluster(clus_i).topo)
                    % Using this custom modified code to allow taking average within participant for each cluster
                    TMP_STUDY = std_readtopoclust_CL(TMP_STUDY,TMP_ALLEEG,clus_i);
                end
            end
            %## SUBJECT SPECIFIC PER CLUSTER
            subj_inds = TMP_STUDY.cluster(cl_to_plot).sets;
            cl_comps = TMP_STUDY.cluster(cl_to_plot).comps;
            for s_i = 1:length(subj_inds)
                subj_ind = subj_inds(s_i);
                comp_i = cl_comps(s_i);
                subj_char = TMP_STUDY.datasetinfo(subj_ind).subject; %TMP_STUDY.subject{subj_ind};
                subj_save_dir = [cluster_dir filesep sprintf('%i',CLUSTER_PICKS(k))];
                if ~exist(subj_save_dir,'dir')
                    mkdir(subj_save_dir);
                end
                fprintf('Making Plots for Subject %s\n...',subj_char);
                %- (TOPOPLOT) 
            %             set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
            %                 figure;
                std_topoplot(TMP_STUDY,TMP_ALLEEG,'clusters',cl_to_plot,'comps',s_i);
                fig_i = get(groot,'CurrentFigure');
                set(fig_i,'position',[16 100 500 350],'color','w');
                drawnow;
                for c = 2:length(fig_i.Children)
            %                 set(fig_i.Children(c).Title,'Interpreter','none');
                    fig_i.Children(c).Title.Interpreter = 'none';
                    fig_i.Children(c).TitleFontSizeMultiplier = 1.4;
                end
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_topo_ic%i.jpg',subj_char,comp_i)]);
                %- (DIPOLE) Plot dipole clusters 
                TMP_STUDY.etc.dipparams.centrline = 'off';
            %                     std_dipplot_CL(TMP_STUDY,TMP_ALLEEG,'clusters',cl_to_plot,'comps',s_i,...
            %                         'figure','off','mode','apart','spheres','off','projlines','off');
                figure;
                tmp = linspecer(2);
                options = {'projlines','off',...
                    'axistight','off',...
                    'projimg','off',...
                    'spheres','off',...
                    'dipolelength',0,...
                    'density','off',...
                    'gui','off',...
                    'cornermri','on',...
                    'mri',TMP_ALLEEG(subj_ind).dipfit.mrifile,...
                    'coordformat',TMP_ALLEEG(subj_ind).dipfit.coordformat,...
                    'color',{tmp(1,:),tmp(2,:)},...
                    'meshdata',TMP_ALLEEG(subj_ind).dipfit.hdmfile};
                dip1 = TMP_STUDY.cluster(cl_to_plot).dipole;
                dip2 = [];
                dip2.posxyz = TMP_STUDY.cluster(cl_to_plot).all_diplocs(s_i,:);
                dip2.momxyz = [0,0,0];
                dip2.rv = TMP_STUDY.cluster(cl_to_plot).residual_variances(s_i);
                dipplot([dip1,dip2],options{:});
                fig_i = get(groot,'CurrentFigure');
                set(fig_i,'position',[16 582 300 350],'color','w')
                set(fig_i, 'DefaultAxesTickLabelInterpreter', 'none')
                camzoom(1.2^2);
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_top_ic%i.jpg',subj_char,comp_i)]);
                view([45,0,0])
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_coronal_ic%i.jpg',subj_char,comp_i)]);
                view([0,-45,0])
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_sagittal_ic%i.jpg',subj_char,comp_i)]);
                %- (SPEC) Spec plot conds for des_i and all groups
                fprintf('Plotting Spectograms for Conditions...\n');
                for ii = 1:length(TMP_ALLEEG)
                    TMP_ALLEEG(ii).group = tmp_group_unif{ii};
                    TMP_STUDY.datasetinfo(ii).group = tmp_group_unif{ii};
                end
                for des_i = 1:length(COND_DESIGNS)
                    [TMP_STUDY] = std_makedesign(TMP_STUDY,TMP_ALLEEG,des_i,...
                            'subjselect', {subj_char},...
                            'variable1',COND_EVENT_CHAR,...
                            'values1',COND_DESIGNS{des_i});
                    std_specplot(TMP_STUDY,TMP_ALLEEG,'clusters',cl_to_plot,'comps',s_i,...
                        'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed','design',des_i);
                    fig_i = get(groot,'CurrentFigure');
                    fig_i.Position = [16 582 420 360];
                    %- set figure line colors
                    cc = linspecer(length(COND_DESIGNS{des_i}));
                    iter = 1;
                    for d = 1:length(fig_i.Children(2).Children)
                        %- pane 1
                        set(fig_i.Children(2).Children(d),'LineWidth',1.5);
                        set(fig_i.Children(2).Children(d),'Color',horzcat(cc(iter,:),0.6));

                        if iter == size(cc,1)
                            iter = 1;
                        else
                            iter = iter + 1;
                        end                
                    end
                    set(fig_i.Children(2),'FontSize',13)
            %                     set(fig_i.Children(3),'FontSize',13)
                    set(fig_i.Children(2),'Position',[0.20,0.20,0.7,0.7]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
            %                     set(fig_i.Children(3),'Position',[0.20,0.20-0.0255,0.7,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
                    set(fig_i.Children(1),'Location','northeast') %reset Legend
                    drawnow;
                    saveas(fig_i,[subj_save_dir filesep sprintf('%s_psd_des%i_ic%i.jpg',subj_char,des_i,comp_i)]);
                end
                close all
            end
        end
    end
          
end
%}