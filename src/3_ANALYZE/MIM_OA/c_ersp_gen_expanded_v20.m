%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_c_ersp_gen_expanded_v20.sh

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
FORCE_RECALC_SPEC = false;
FORCE_RECALC_ERSP = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
DO_SUBJ_PLOTS = true;
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
clustering_method = 'dipole_1'; 
% ['dipole_',num2str(clustering_weights.dipoles),...
%     '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
%     '_spec_',num2str(clustering_weights.spec)];
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
% SPEC_PARAMS = struct('freqrange',[1,200],...
%     'subject','',...
%     'specmode','psd',...
%     'freqfac',4,...
%     'logtrials','on',...
%     'comps','all',...
%     'plot_freqrange',[4,60]);
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
    'freqrange',[1,200]);
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
%- convert SUB_DIR
SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster';
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true,true};
CLUSTER_K_PICKS = [14,14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_6' filesep '14'],...
    [SUB_DIR filesep 'icrej_5' filesep '14']};
CLUSTER_FILES = {'cl_inf_14.mat','cl_inf_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_6'],...
    [SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
%% (STEP 1) GENERATE ERSP & SPEC DATA FOR-EACH DESIGN & CLUSTER
%## NOTE: This Loop ABSOLUTELY CAN NOT be ran in parallel at this point.
for k_i = 1:length(CLUSTER_K_PICKS)
    %## TEMPORARIES
    parameters = []; %#ok<NASGU>
    tmp_group_orig = cell(length(ALLEEG),1);
    tmp_group_unif = cell(length(ALLEEG),1);
    %## LOAD STUDY
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
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '_UNIX.study'],'filepath',cluster_study_dir);
        else
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '.study'],'filepath',cluster_study_dir);
        end
    end
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        mkdir(spec_data_dir);
    end
    %% CALCULATE GRANDAVERAGE WARPTOs
    for subj_i = 1:length(ALLEEG)
        %- assign percondition timewarping
        ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    %     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    end
    allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    % allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    for subj_i = 1:length(ALLEEG)
        allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
    end
    % grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
    averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
    %% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
    TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
    ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)];
    STUDY.etc.averaged_warpto_events = averaged_warpto_events;
    fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
    disp(averaged_warpto_events);
    %## ersp plot per cluster per condition
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
    
    %% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
    tmp = strsplit(ALLEEG(1).filename,'.');
    spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
    topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
    if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
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
            TMP_STUDY.datasetinfo = TMP_STUDY.datasetinfo(subj_i);
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
    %% (PRECOMPUTE MEASURES) COMPUTE ERSPs
    icatimf_f = [ALLEEG(1).filepath filesep sprintf('%s.icatimef',ALLEEG(1).subject)];
    if ~exist(icatimf_f,'file') || FORCE_RECALC_ERSP
        disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
        parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/3))
            EEG = ALLEEG(subj_i);
            TMP_STUDY = STUDY;
            EEG = eeg_checkset(EEG,'loaddata');
            if isempty(EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            end
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
            %- overrride datasetinfo to trick std_precomp to run.
            TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
            TMP_STUDY.datasetinfo(1).index = 1;
            %- determine timewarping parameters
             if DO_TIMEWARP
                timewarp_param = EEG.timewarp.latencies;
                timewarpms_param = averaged_warpto_events;
             else
                 timewarp_param = [];
                 timewarpms_param = [];
            end
            %-
            if DO_BASELINE_CORRECTION
                % Baseline correction
                [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                        'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                        'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                        'trialbase','off','basenorm','on'}); %ERSP
            else
                % No baseline correction
                [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
                        'baseline',nan(),'timewarp',timewarp_param,...
                        'timewarpms',timewarpms_param}); %ERSP
            end
        end
    end
    %- load .icatimef load-in parameters
    tmp = load(icatimf_f,'-mat');
    parameters = tmp.parameters;
    %     warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
    parameters{find(strcmp(parameters,'baseline'))+1} = [averaged_warpto_events(1),averaged_warpto_events(end)];
    parameters = [parameters, {'trialbase'}, {'off'}];
    cellArray = parameters(1,2:2:length(parameters));
    fields = parameters(1,1:2:length(parameters)-1);
    parameters = cell2struct(cellArray,fields,2);
    %% MAKE DESIGNS
    %## NOTE (07/18/2023) JS, scripts/functions adapted from CL, AS, NJ scripts:
    % PlotAndSaveERSP_CL_V3.m, Stats_Plotting_ERSPs_local_allCondBaseline_V3.m,
    % Figures_Plotting_ERSPs_local_allCondBaseline_V3.m, plotERSP_parfor.m
    STUDY.cache = [];
    for subj_i = 1:length(ALLEEG)
        tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
        tmp_group_unif{subj_i} = 'Older Adults';
    end
    %## Make Study Designs
    %- combine groups
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_unif{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
    end
    %- assign studies
    STUDY.cache = [];
    for des_i = 1:length(STUDY_DESI_PARAMS)
        [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
    end
    %## Get Cluster Information
    clust_i = CLUSTER_K_PICKS(k_i);
    fprintf('Making data for K=%i...\n',clust_i);
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        mkdir(spec_data_dir)
    end
    %- assign cluster information
    STUDY.cluster = cluster_update;
    %- get inds
    [~,~,~,~,~,nonzero_cl_inds] = eeglab_get_cluster_comps(STUDY);
    CLUSTER_PICKS = nonzero_cl_inds; %1:length(TMP_STUDY.cluster);
    %- stores
    cnt = 1;
    cnt2 = 1;
    spec_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    spec_ss_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    ersp_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    ersp_norm_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    ersp_normcb_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    ersp_singtrial_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    load_ind_cl = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    clust_ind_cl = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    des_ind = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    cnt_ind = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
    %- loop designs
    for des_i = 1:length(STUDY.design)
        store_1 = cell(length(CLUSTER_PICKS),1);
        store_2 = cell(length(CLUSTER_PICKS),1);
        store_3 = cell(length(CLUSTER_PICKS),1);
        store_4 = cell(length(CLUSTER_PICKS),1);
        store_5 = cell(length(CLUSTER_PICKS),1);
        store_6 = cell(length(CLUSTER_PICKS),1);
        store_7 = cell(length(CLUSTER_PICKS),1);
        store_8 = cell(length(CLUSTER_PICKS),1);
        store_9 = cell(length(CLUSTER_PICKS),1);
        store_10 = cell(length(CLUSTER_PICKS),1);
        STUDY.currentdesign = des_i;
        design_char = [];
        for i = 1:length(STUDY.design(des_i).variable)
            if i == 1
                design_char = [STUDY.design(des_i).variable(1).value{:}];
            else
                design_char = [design_char '_' STUDY.design(des_i).variable(i).value{:}];
            end
        end
        %## (PARFOR) Compute Specs & ERSPs
        parfor (i = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
%         for i = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(i);
            %- generate power spectrums for each cluster and condition
            [~,spec_savef,spec_subjcorr_savef] = mim_gen_spec_data(STUDY,ALLEEG,...
                    averaged_warpto_events,des_i,cluster_i,design_char,spec_data_dir);
            store_1{i} = spec_savef;
            store_2{i} = spec_subjcorr_savef;
            %## DEBUG
            %{
            ersp_savef = [];
            ersp_subbase_savef = [];
            ersp_subbase_combase_savef = [];
            ersp_singletrial_subbase_savef = [];
            %}
            %- generate event related spectral perturbations for each cluster and condition
            [~,ersp_savef,ersp_subbase_savef,ersp_subbase_combase_savef,ersp_singletrial_subbase_savef] = mim_gen_ersp_data(STUDY,ALLEEG,averaged_warpto_events,...
                    parameters,des_i,cluster_i,design_char,spec_data_dir,...
                    'STAT_PARAMS',ERSP_STAT_PARAMS);
            store_3{i} = ersp_savef;
            store_4{i} = ersp_subbase_savef;
            store_5{i} = ersp_subbase_combase_savef;
            store_6{i} = ersp_singletrial_subbase_savef;
            store_7{i} = i;
            store_8{i} = cluster_i;
            store_9{i} = des_i;
        end
        spec_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_1;
        spec_ss_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_2;
        ersp_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1)  = store_3;
        ersp_norm_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1)  = store_4;
        ersp_normcb_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1)  = store_5;
        ersp_singtrial_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_6;
        load_ind_cl(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_7;
        clust_ind_cl(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_8;
        des_ind(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_9;
        for k = cnt:cnt+length(CLUSTER_PICKS)-1
            cnt_ind{k,1} = k;
        end
        cnt = cnt + length(CLUSTER_PICKS);
    end
    %## Store Paths in Struct
    gen_data_struct = struct('ersp_fpaths',ersp_fpaths,...
        'spec_fpaths',spec_fpaths,...
        'spec_ss_fpaths',spec_ss_fpaths,...
        'ersp_norm_fpaths',ersp_norm_fpaths,...
        'ersp_normcb_fpaths',ersp_normcb_fpaths,...
        'ersp_singtrial_fpaths',ersp_singtrial_fpaths,...
        'load_ind_cl',load_ind_cl,...
        'clust_ind_cl',clust_ind_cl,...
        'des_ind',des_ind,...
        'cnt_ind',cnt_ind);
    STUDY.etc.mim_gen_ersp_data = gen_data_struct;
    par_save(gen_data_struct,spec_data_dir,'spec_data_struct.mat');
    [~,~] = parfunc_save_study(STUDY,ALLEEG,...
                                    STUDY.filename,spec_data_dir,...
                                    'RESAVE_DATASETS','off');
    %% CLUSTER DIAGNOSTIC PLOTS && SUBJECT SPECIFIC PER CLUSTER
    %{
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
    %## PLOT cluster based information
%     mim_gen_cluster_figs(STUDY,ALLEEG,cluster_dir,...
%         'CLUSTERS_TO_PLOT',main_cl_inds);
    %##
    STUDY.etc.dipparams.centrline = 'off';
    parfor (k = 1:length(main_cl_inds),length(main_cl_inds))
        clust_i = main_cl_inds(k);        
        subj_inds = STUDY.cluster(clust_i).sets;
        cl_comps = STUDY.cluster(clust_i).comps;
        %-
        for s_i = 1:length(subj_inds)
            subj_k = subj_inds(s_i);
            comp_i = cl_comps(s_i);
            subj_char = STUDY.datasetinfo(subj_k).subject;
            fprintf('\n(Cluster=%i) Plotting For Subject %s\n',clust_i,subj_char);
            subj_save_dir = [cluster_dir filesep sprintf('%i',clust_i)];
            if ~exist(subj_save_dir,'dir')
                mkdir(subj_save_dir);
            end
            %- (TOPOPLOT)
            std_topoplot(STUDY,ALLEEG,'clusters',clust_i,'comps',s_i);
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'position',[16 100 500 350],'color','w');
            drawnow;
            for c = 2:length(fig_i.Children)
                fig_i.Children(c).Title.Interpreter = 'none';
                fig_i.Children(c).TitleFontSizeMultiplier = 1.4;
            end
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_topo_ic%i.jpg',subj_char,comp_i)]);
            %- (DIPOLE) Plot dipole clusters
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
                'mri',ALLEEG(subj_k).dipfit.mrifile,...
                'coordformat',ALLEEG(subj_k).dipfit.coordformat,...
                'color',{tmp(1,:),tmp(2,:)},...
                'meshdata',ALLEEG(subj_k).dipfit.hdmfile};
            dip1 = STUDY.cluster(clust_i).dipole;
            tmp = ALLEEG(subj_k).dipfit.model(comp_i);
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
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_top_ic%i.jpg',subj_char,comp_i)]);
            view([45,0,0])
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_coronal_ic%i.jpg',subj_char,comp_i)]);
            view([0,-45,0])
            saveas(fig_i,[subj_save_dir filesep sprintf('%s_dip_sagittal_ic%i.jpg',subj_char,comp_i)]);
            %## Find Anatomy
            labels = cell(length(ATLAS_FPATHS),3);
            for atlas_i = 1:length(ATLAS_FPATHS)
                atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
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
                tmp = strsplit(ATLAS_FPATHS{atlas_i},filesep);
                labels{atlas_i,3} = tmp{end};
            end
            % Convert cell to a table and use first row as variable names
            T = cell2table(labels,'VariableNames',{'cluster_centroid','subject_dipole','atlas'});
            % Write the table to a CSV file
            writetable(T,[subj_save_dir filesep sprintf('%s_atlasinf_ic%i.csv',subj_char,comp_i)])
            %- (SPEC) Spec plot conds for des_i and all groups
            fprintf('Plotting Spectograms for Conditions...\n');
            for des_i = 1:length(STUDY.design)
                std_specplot(STUDY,ALLEEG,'clusters',clust_i,'comps',s_i,...
                    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed','design',des_i);
                fig_i = get(groot,'CurrentFigure');
%                 fig_i.Position = [16 582 420 360];
                set(fig_i,'position',[16 582 420 360],'color','w')
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
                set(fig_i.Children(2),'FontSize',13);
                set(fig_i.Children(2),'Position',[0.20,0.20,0.7,0.7])
                set(fig_i.Children(1),'Location','northeast') %reset Legend
                drawnow;
                saveas(fig_i,[subj_save_dir filesep sprintf('%s_psd_des%i_ic%i.jpg',subj_char,des_i,comp_i)]);
            end
            close all            
        end
    end
    %}
    
end
%% (STEP 2) PLOT
%##
%{
clear('STUDY');
clear('ALLEEG');
clear('TMP_ALLEEG');
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
    if LOAD_DIFFERENT_STUDY{k_i}
        %- Create STUDY & ALLEEG structs
        if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
            error('error. study file does not exist');
        else
            if ~ispc
                tmp = load('-mat',[spec_data_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
                STUDY = tmp.STUDY;
            else
                tmp = load('-mat',[spec_data_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
                STUDY = tmp.STUDY;
            end
        end
    else
        error('error. Define a valid study path...');
    end
    %## RE-POP PARAMS
    STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
          'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_CROP_TIMES);
    %## Cluster Update
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    STUDY.cluster = cluster_update;
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
    %## PLOT cluster based information
%     mim_gen_cluster_figs(STUDY,ALLEEG,CLUSTER_DIRS{k_i},...
%         'CLUSTERS_TO_PLOT',main_cl_inds);
    %## Loop Through Designs
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
%             cluster_load_ind = cluster_i;
%             cluster_load_ind = TMP_STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).cluster_n(1,cluster_i);
            %- defaults
            allersp = {};
            alltimes = [];
            allfreqs = [];
            pcond = {};
            pgroup = {};
            pinter = {};
            %- LOAD OPTION 1
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data = par_load(fpath,fname);
    %         allersp = ersp_data.allerspdata;
    %         alltimes = ersp_data.alltimes;
    %         allfreqs = ersp_data.allfreqs;
            %- LOAD OPTION 2
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_norm_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data_norm = par_load(fpath,fname);
    %         allersp = ersp_data_norm.allerspdata;
    %         alltimes = ersp_data_norm.alltimes;
    %         allfreqs = ersp_data_norm.allfreqs;
    %         pcond = ersp_data_norm.pcond;
            %- LOAD OPTION 3
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_normcb_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data_normcb = par_load(fpath,fname);
    %         allersp = ersp_data_normcb.allerspdata;
    %         alltimes = ersp_data_normcb.alltimes;
    %         allfreqs = ersp_data_normcb.allfreqs;
            %- LOAD OPTION 4
    %         [~, allersp, alltimes, allfreqs, ~, ~] = std_readdata(STUDY,ALLEEG,...
    %             'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
    %             'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
    %             'design',des_i);
    %         %* get stats
    %         [pcond_ersp,pgroup_ersp,pinter_ersp,~,~,~] = std_stat(allersp,...
    %                     'condstats', ERSP_CONDSTATS,...
    %                     'groupstats',ERSP_GROUPSTATS,...
    %                     'method',ERSP_STAT_METHOD,...
    %                     'naccu',ERSP_NACCU,...
    %                     'alpha',ERSP_ALPHA,...
    %                     'mcorrect',ERSP_MCORRECT,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
    %                     'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
            %## RUN PLOTTING
            fprintf('Plotting Cluster %i for design %i\n',cluster_i,des_i);
            mim_custom_ersp_plots(STUDY,cond_test,averaged_warpto_events,...
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
end
%}
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}