%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_ersp_plots_opt.sh

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
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
%     SLURM_POOL_SIZE = 2;
%     pp = parcluster('local');
%     pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
end
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = true;
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
DO_SUBJ_PLOTS = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
do_multivariate_data = 1;
evaluate_method = 'min_rv';
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',... %'fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result?
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% (08/03/2023) JS, changing fieldtripnaccu to 10,000; changing
% fieldtripmcorrect to cluster; method to perm; (these are the parameters
% set inside CL's PlotAndSaveERSP_CL_V3.m...
% pipeline although this doesn't align with her YA manuscript methods?
% (08/06/2023) JS, changing fieldtripnaccu to 2000 again and mcorrect to fdr...
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
% dt = '10252023_MIM_OAN70_noslowwalkers_gait_powpow0p25';
dt = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
%## Soft Define
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%- (EDIT!) convert SUB_DIR
SUB_DIR = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [12];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_5' filesep '12']};
CLUSTER_FILES = {'cl_inf_12.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
SUB_GROUP_FNAME = []; %'H3000'; %[]; %'H2000';
% SUB_GROUP_FNAME = 'H3000'; %[]; %'H2000';
SUB_GROUP_FNAME_REGEX = []; %'H3000''s'; %[]; %'H2000''s';
% SUB_GROUP_FNAME_REGEX = 'H3000''s'; %[]; %'H2000''s';
CLUSTER_CLIM_MATCH = [];
%% (STEP 2) PLOT
%##
for k_i = 1:length(CLUSTER_K_PICKS)
% parfor (k_i = 1:length(CLUSTER_K_PICKS),length(CLUSTER_K_PICKS))
    fprintf('Loading Cluster K=%i',CLUSTER_K_PICKS(k_i));
    %## Loop Params
    clust_i = CLUSTER_K_PICKS(k_i);
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    if ~isempty(SUB_GROUP_FNAME_REGEX)
        spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
        plot_store_dir = [cluster_dir filesep 'plots_out' filesep SUB_GROUP_FNAME];
    else
        spec_data_dir = [cluster_dir filesep 'spec_data'];
        plot_store_dir = [cluster_dir filesep 'plots_out'];
    end
    if ~exist(spec_data_dir,'dir')
        error('spec_data dir does not exist');
    end
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    
    %## Load Study
    if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '_UNIX.study'],'filepath',spec_data_dir);
        else
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '.study'],'filepath',spec_data_dir);
        end
    end
    
    %## CALCULATE GRANDAVERAGE WARPTOs
    for subj_i = 1:length(ALLEEG)
        %- assign percondition timewarping
        ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    end
    allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    for subj_i = 1:length(ALLEEG)
        allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
    end
    averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
    %## (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
    TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
    ERSP_TIMERANGE=[averaged_warpto_events(1), averaged_warpto_events(end)];
    STUDY.etc.averaged_warpto_events = averaged_warpto_events;
    fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
    disp(averaged_warpto_events);
    %## RE-POP PARAMS
    STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
          'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_TIMERANGE);
    %## Cluster Update
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
    %## PLOT cluster based information
    mim_gen_cluster_figs(STUDY,ALLEEG,CLUSTER_DIRS{k_i},...
        'CLUSTERS_TO_PLOT',main_cl_inds);
    %## generate iter pairs
    cl_inds = [STUDY.etc.mim_gen_ersp_data.clust_ind_cl];
    des_inds = [STUDY.etc.mim_gen_ersp_data.des_ind];
    tmp = [cl_inds', des_inds'];
    iter_pairs = [];
    for c_i = 1:length(CLUSTER_PICKS)
        inds = CLUSTER_PICKS(c_i)==tmp(:,1);
        iter_pairs = cat(1,iter_pairs,tmp(inds,:));
    end
    %## PARFOR
    parfor (cnt = 1:size(iter_pairs,1),size(iter_pairs,1))
        %- INITIATE ITERS
        TMP_STUDY = STUDY;
        cl_inds = [TMP_STUDY.etc.mim_gen_ersp_data.clust_ind_cl];
        des_inds = [TMP_STUDY.etc.mim_gen_ersp_data.des_ind];
        des_i = iter_pairs(cnt,2);
        cluster_i = iter_pairs(cnt,1);
        cond_test = TMP_STUDY.design(des_i).variable(1).value;
        cluster_load_ind = find(logical(cl_inds == cluster_i) & logical(des_inds == des_i));
        %- PRINT OUTS
        fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
        TMP_STUDY.currentdesign = des_i;
        fprintf('Current design: %i\n',TMP_STUDY.currentdesign);
        fprintf('Statistics Parameters:\n');
        disp(TMP_STUDY.etc.statistics)
        fprintf('Statistics Fieldtrip Parameters:\n');
        disp(TMP_STUDY.etc.statistics.fieldtrip)
        fprintf('Statistics EEGLAB Parameters:\n');
        disp(TMP_STUDY.etc.statistics.eeglab)
        fprintf('ERSP Parameters:\n');
        disp(TMP_STUDY.etc.erspparams)
        %- defaults
        allersp = {};
        alltimes = [];
        allfreqs = [];
        pcond = {};
        pgroup = {};
        pinter = {};%## RUN PLOTTING
        fprintf('Plotting Cluster %i for design %i\n',cluster_i,des_i);
        mim_custom_ersp_plots(TMP_STUDY,cond_test,averaged_warpto_events,...
            cluster_i,cluster_load_ind,des_i,plot_store_dir,...
            'DO_SUBJ_PLOTS',DO_SUBJ_PLOTS,...
            'CLUSTER_CLIM_MATCH',CLUSTER_CLIM_MATCH,...
            'ALLERSP',allersp,...
            'ALLTIMES',alltimes,...
            'ALLFREQS',allfreqs,...
            'PCOND',pcond,...
            'PGROUP',pgroup,...
            'PINTER',pinter)
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}