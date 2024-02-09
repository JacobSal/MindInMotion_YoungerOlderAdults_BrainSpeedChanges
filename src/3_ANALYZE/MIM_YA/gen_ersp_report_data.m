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

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_ersp_report_data.sh

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
FORCE_RECALC_SPEC = false;
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
DISTS_FROM_MEDIAN = [30]; %[10,20,30,40,50,60,70]; % distance in units (mm)
K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
CLUSTER_SWEEP_VALS = [10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = false;
DO_K_SWEEPING = false;
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
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
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
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
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
%- NOTE: partly adapt from bemobil_repeated_clustering
numIC = zeros(length(STUDY.datasetinfo),1);
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
%% ===================================================================== %%
%## ADMIN SET
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
atlases = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'brainnetome' filesep 'BNA_MPM_thr25_1.25mm.nii'],...
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat'],...
    [ATLAS_PATH filesep 'vtpm' filesep 'vtpm.mat'],...
    [ATLAS_PATH filesep 'yeo' filesep 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'],...
    [ATLAS_PATH filesep 'brainweb' filesep 'brainweb_discrete.mat']}; % also a discrete version of this

SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\dipole_1_scalp_0_ersp_0_spec_0';
%- convert SUB_DIR
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true,true,true,true};
CLUSTER_K_PICKS = [14,14,14,14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics5','temp_study_rejics3','temp_study_rejics7'};
CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics5' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics3' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics7' filesep '14']};
CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat','cluster_update_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6'],...
    [SUB_DIR filesep 'subjrejs_minics5'],...
    [SUB_DIR filesep 'subjrejs_minics3'],...
    [SUB_DIR filesep 'subjrejs_minics7']};
%%
% parfor (k_i = 1:length(CLUSTER_STUDY_DIRS),length(CLUSTER_STUDY_DIRS))
for k_i = 1:length(CLUSTER_STUDY_DIRS)
    %- convert cluster directory
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    %- convert study directory
    if ~ispc
        cluster_study_dir = convertPath2UNIX(CLUSTER_STUDY_DIRS{k_i});
    else
        cluster_study_dir = convertPath2Drive(CLUSTER_STUDY_DIRS{k_i});
    end
    if ~exist([cluster_study_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '_UNIX.study'],'filepath',cluster_study_dir);
        else
            [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '.study'],'filepath',cluster_study_dir);
        end
    end
    %## grab subjects for study designs
    tmp_group_orig = cell(length(TMP_ALLEEG),1);
    tmp_group_unif = cell(length(TMP_ALLEEG),1);
    for subj_i = 1:length(TMP_ALLEEG)
        tmp_group_orig{subj_i} = TMP_ALLEEG(subj_i).group;
        tmp_group_unif{subj_i} = 'Older Adults';
    end
    %- make subjects a uniform group
    for subj_j = 1:length(TMP_ALLEEG)
        TMP_ALLEEG(subj_j).group = tmp_group_unif{subj_j};
        TMP_STUDY.datasetinfo(subj_j).group = tmp_group_unif{subj_j};
    end
    [TMP_STUDY] = std_makedesign(TMP_STUDY,TMP_ALLEEG,1,'subjselect',{TMP_ALLEEG.subject},...
        'variable1','cond','values1',{'flat','low','med','high'},...
        'variable2','group','values2',{});
    [TMP_STUDY] = std_makedesign(TMP_STUDY,TMP_ALLEEG,2,'subjselect',{TMP_ALLEEG.subject},...
        'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
        'variable2','group','values2',{});
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
    %## LOAD cluster information
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    TMP_STUDY.cluster = cluster_update;
    [~,main_cl_inds,~,valid_cluster,~,nonzero_cluster] = eeglab_get_cluster_comps(TMP_STUDY);
    %## PLOT cluster based information
    mim_gen_cluster_figs(TMP_STUDY,TMP_ALLEEG,CLUSTER_DIRS{k_i},...
        'CLUSTERS_TO_PLOT',main_cl_inds);
    %% SUBJECT SPECIFIC PER CLUSTER
    for k = 1:length(main_cl_inds)
        clust_i = main_cl_inds(k);        
        subj_inds = TMP_STUDY.cluster(clust_i).sets;
        cl_comps = TMP_STUDY.cluster(clust_i).comps;
        %-
        for s_i = 1:length(subj_inds)
            subj_k = subj_inds(s_i);
            comp_i = cl_comps(s_i);
            subj_char = TMP_STUDY.datasetinfo(subj_k).subject;
            fprintf('\n(Cluster=%i) Plotting For Subject %s\n',clust_i,subj_char);
            subj_save_dir = [cluster_dir filesep sprintf('%i',clust_i)];
            if ~exist(subj_save_dir,'dir')
                mkdir(subj_save_dir);
            end
            %- (TOPOPLOT)
            std_topoplot(TMP_STUDY,TMP_ALLEEG,'clusters',clust_i,'comps',s_i);
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'position',[16 100 500 350],'color','w');
            drawnow;
            for c = 2:length(fig_i.Children)
                fig_i.Children(c).Title.Interpreter = 'none';
                fig_i.Children(c).TitleFontSizeMultiplier = 1.4;
            end
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_topo_ic%i.jpg',subj_char,comp_i)]);
            %- (DIPOLE) Plot dipole clusters 
            TMP_STUDY.etc.dipparams.centrline = 'off';
            figure;
            tmp = linspecer(2);
            options = {'projlines','off',...
                'axistight','off',...
                'projimg','off',...
                'spheres','off',...
                'dipolelength',0.1,...
                'density','off',...
                'holdon','on',...
                'gui','off',...
                'mri',TMP_ALLEEG(subj_k).dipfit.mrifile,...
                'coordformat',TMP_ALLEEG(subj_k).dipfit.coordformat,...
                'color',{tmp(1,:),tmp(2,:)},...
                'meshdata',TMP_ALLEEG(subj_k).dipfit.hdmfile};
            dip1 = TMP_STUDY.cluster(clust_i).dipole;
            tmp = TMP_ALLEEG(subj_k).dipfit.model(comp_i);
            dip2 = [];
            dip2.posxyz = tmp.posxyz;
            dip2.momxyz = tmp.momxyz;
            dip2.rv = tmp.rv;
            %- plot dipole
            dipplot([dip1,dip2],options{:});
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'position',[16 582 300 350],'color','w')
            set(fig_i, 'DefaultAxesTickLabelInterpreter', 'none')
            camzoom(1.2^2);
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_top_cl%i_ic%i.jpg',subj_char,clust_i,comp_i)]);
            view([45,0,0])
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_coronal_cl%i_ic%i.jpg',subj_char,clust_i,comp_i)]);
            view([0,-45,0])
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_sagittal_cl%i_ic%i.jpg',subj_char,clust_i,comp_i)]);
            %## Find Anatomy
            labels = cell(length(atlases),3);
            for atlas_i = 1:length(atlases)
                atlas = ft_read_atlas(atlases{atlas_i});
                cfg              = [];
                cfg.roi        = dip1.posxyz;
                cfg.output     = 'multiple';
                cfg.atlas      = atlas;
                cfg.inputcoord = 'mni';
                %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
                cfg.sphere = 3;
                label_i = ft_volumelookup(cfg, atlas);
                if ~isempty(label_i)
                    % disp(labels.count(labels.count ~= 0))
                    [val, indx] = max(label_i.count);
                    if strcmp(label_i.name(indx),'no_label_found')
                        sub_indx = find(label_i.count ~= 0 & label_i.count < val);
                        if isempty(sub_indx)
                            atlas_name = label_i.name{indx};
                        end
                    else
                        atlas_name = label_i.name{indx};
                    end
                end
                labels{atlas_i,1} = atlas_name;
                cfg.roi        = dip2.posxyz;
                label_i = ft_volumelookup(cfg, atlas);
                if ~isempty(label_i)
                    % disp(labels.count(labels.count ~= 0))
                    [val, indx] = max(label_i.count);
                    if strcmp(label_i.name(indx),'no_label_found')
                        sub_indx = find(label_i.count ~= 0 & label_i.count < val);
                        if isempty(sub_indx)
                            atlas_name = label_i.name{indx};
                        end
                    else
                        atlas_name = label_i.name{indx};
                    end
                end
                labels{atlas_i,2} = atlas_name;
                tmp = strsplit(atlases{atlas_i},filesep);
                labels{atlas_i,3} = tmp{end};
            end
            % Convert cell to a table and use first row as variable names
            T = cell2table(labels,'VariableNames',{'cluster_centroid','subject_dipole','atlas'});
            % Write the table to a CSV file
            writetable(T,[subj_save_dir filesep sprintf('%s_atlasinf_ic%i.csv',subj_char,comp_i)])
            %- (SPEC) Spec plot conds for des_i and all groups
            fprintf('Plotting Spectograms for Conditions...\n');
            for des_i = 1:length(TMP_STUDY.design)
                std_specplot(TMP_STUDY,TMP_ALLEEG,'clusters',clust_i,'comps',s_i,...
                    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed','design',des_i);
                fig_i = get(groot,'CurrentFigure');
                fig_i.Position = [16 582 420 360];
                %- set figure line colors
                cc = linspecer(length(TMP_STUDY.design(des_i).variable.value));
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
                set(fig_i.Children(2),'FontSize',13);
                set(fig_i.Children(2),'Position',[0.20,0.20,0.7,0.7])
                set(fig_i.Children(1),'Location','northeast') %reset Legend
                drawnow;
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_psd_des%i_cl%i_ic%i.jpg',subj_char,des_i,clust_i,comp_i)]);
            end
            close all
        end
    end
end