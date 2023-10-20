%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_ersp_plots.sh

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
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
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
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
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
    'plot_freqrange',[4,60]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
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
%## ADMIN SET
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
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
% LOAD_DIFFERENT_STUDY = {true,true};
% CLUSTER_K_PICKS = [14,14];
% CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics5'};
% CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14'],...
%     [SUB_DIR filesep 'subjrejs_minics5' filesep '14']};
% CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat'};
% CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6'],...
%     [SUB_DIR filesep 'subjrejs_minics5']};
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6'};
CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14']};
CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
%% (STEP 2) PLOT
%##
% for k_i = 1:length(CLUSTER_K_PICKS)
parfor (k_i = 1:length(CLUSTER_K_PICKS),length(CLUSTER_K_PICKS))
    fprintf('Loading Cluster K=%i',CLUSTER_K_PICKS(k_i));
    %## Loop Params
    clust_i = CLUSTER_K_PICKS(k_i);
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    plot_store_dir = [cluster_dir filesep 'plots_out'];
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        error('error. path %s doesn''t exist',spec_data_dir);
    end
    %## Load Study
    % (08/03/2023) JS, this can be optimized in the future by only loding
    % in the file strucutre and maintaining a STUDY file with needed info
%     if LOAD_DIFFERENT_STUDY{k_i}
%         %- Create STUDY & ALLEEG structs
%         if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
%             error('error. study file does not exist');
%         else
%             if ~ispc
%                 tmp = load('-mat',[spec_data_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
%                 STUDY = tmp.STUDY;
%             else
%                 tmp = load('-mat',[spec_data_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
%                 STUDY = tmp.STUDY;
%             end
%         end
%     else
%         error('error. Define a valid study path...');
%     end
    %##
    if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [STUDY,ALLEEG] = pop_loadstudy('filename',sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i}),'filepath',spec_data_dir);
        else
            [STUDY,ALLEEG] = pop_loadstudy('filename',sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i}),'filepath',spec_data_dir);
        end
    end
    %## grab subjects for study designs
    tmp_group_orig = cell(length(ALLEEG),1);
    tmp_group_unif = cell(length(ALLEEG),1);
    for subj_i = 1:length(ALLEEG)
        tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
        tmp_group_unif{subj_i} = 'Older Adults';
    end
    %## RE-POP PARAMS
    ERSP_CROP_TIMES = STUDY.etc.erspparams.timerange;
    STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
          'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_CROP_TIMES);
    SPEC_PARAMS.subtractsubjectmean = 'on';
    STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
        'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
        'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
    
    %## Cluster Update
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
    %## PLOT cluster based information
%     mim_gen_cluster_figs(STUDY,TMP_ALLEEG,CLUSTER_DIRS{k_i},...
%         'CLUSTERS_TO_PLOT',main_cl_inds);
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_unif{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
    end
    %% Loop Through Designs
    for des_i = 1:length(STUDY.design)
        cond_test = STUDY.design(des_i).variable(1).value;
        fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
        STUDY.currentdesign = des_i;
        fprintf('Current design: %i\n',STUDY.currentdesign);
        fprintf('Statistics Parameters:\n');
        disp(STUDY.etc.statistics)
        fprintf('Statistics Fieldtrip Parameters:\n');
        disp(STUDY.etc.statistics.fieldtrip)
        fprintf('Statistics EEGLAB Parameters:\n');
        disp(STUDY.etc.statistics.eeglab)
        fprintf('ERSP Parameters:\n');
        disp(STUDY.etc.erspparams)
        cl_inds = [STUDY.etc.mim_gen_ersp_data.clust_ind_cl];
        des_inds = [STUDY.etc.mim_gen_ersp_data.des_ind];
%         des_cls = TMP_STUDY.etc.mim_gen_ersp_data.clust_ind_cl([TMP_STUDY.etc.mim_gen_ersp_data.des_ind] == des_i);
%         parfor (j = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
        for j = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(j);
            cluster_load_ind = find(logical(cl_inds == cluster_i) & logical(des_inds == des_i));
            %- (SPEC) Spec plot conds for des_i and all groups
            fprintf('Plotting Spectograms for Conditions...\n');
            std_specplot(STUDY,ALLEEG,'clusters',cluster_i,...
                'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed','design',des_i);
            fig_i = get(groot,'CurrentFigure');
            fig_i.Position = [16 582 420 360];
            %- set figure line colors
            cc = linspecer(length(STUDY.design(des_i).variable.value));
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
            set(fig_i.Children(3),'FontSize',13)
            set(fig_i.Children(2),'Position',[0.20,0.20,0.7,0.7]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
            set(fig_i.Children(3),'Position',[0.20,0.20-0.0255,0.7,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
            set(fig_i.Children(1),'Location','northeast') %reset Legend
            drawnow;
            saveas(fig_i,[cluster_dir filesep sprintf('allSpecPlot_des%i_cl%i.jpg',des_i,cluster_i)]);
        end
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}