%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_b_ersp_data_gen.sh

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
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','on',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
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
    'freqrange',[1,200]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
% dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
% dt = '10052023_MIM_OAN70_noslowwalkers_gait';
% dt = '10252023_MIM_OAN70_noslowwalkers_gait_powpow0p20';
% dt = '10302023_MIM_OAN70_noslowwalkers_gait_iccREMG0p4_powpow0p1';
% dt = '10302023_MIM_OAN70_newnormalize_iccREMG0p4_powpow0p1';
dt = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
%## Soft Define
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
study_fName_1 = 'epoch_study';
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
% hpg_compatability = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip\copy_private';
% if ~ispc
%     hpg_compatability = convertPath2UNIX(hpg_compatability);
% end
% addpath(hpg_compatability);
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
%- (EDIT!) convert SUB_DIR
SUB_DIR = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [12];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_5' filesep '12']};
CLUSTER_FILES = {'cl_inf_12.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
SUB_GROUP_FNAME = ['test']; %'H3000'; %'H2000';
SUB_GROUP_FNAME_REGEX = []; %'H3000''s'; %'H2000''s';
% SUB_GROUP_FNAME = 'H3000'; %'H2000';
% SUB_GROUP_FNAME_REGEX = 'H3000''s'; %'H2000''s';
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H2000''s','H3000''s'}}};
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable1','cond','values1',{'flat','low','med','high'},...
%             'variable2','group','values2',{}},...
%             {'subjselect',{},...
%             'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
%             'variable2','group','values2',{}}};
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
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '_UNIX.study'],'filepath',cluster_study_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '.study'],'filepath',cluster_study_dir);
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
    SPEC_PARAMS.subtractsubjectmean = 'off';
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
%     inds = find(cellfun(@(x) strcmp(x,SUB_GROUP_FNAME_REGEX),{STUDY.datasetinfo.group}));
%     fprintf('Running subjects:'); cellfun(@(x) fprintf('%s,',x),{STUDY.datasetinfo(inds).subject}); fprintf('\n');
%     for des_i = 1:length(STUDY_DESI_PARAMS)
%         STUDY_DESI_PARAMS{des_i}{10} = {SUB_GROUP_FNAME_REGEX};
%     end
    %## Make Study Designs
    STUDY.cache = [];
    STUDY.design = [];
    for des_i = 1:length(STUDY_DESI_PARAMS)
        [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
%         STUDY.design(des_i).variable(1).pairing ='off';
%         STUDY.design(des_i).variable(2).pairing ='off';        
    end
    %## Get Cluster Information
    clust_i = CLUSTER_K_PICKS(k_i);
    fprintf('Making data for K=%i...\n',clust_i);
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    if ~isempty(SUB_GROUP_FNAME)
        spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
    else
        spec_data_dir = [cluster_dir filesep 'spec_data'];
    end
    if ~exist(spec_data_dir,'dir')
        mkdir(spec_data_dir)
    end
    %- assign cluster information
    STUDY.cluster = cluster_update;
    %- get inds
    [~,main_cl_inds,~,~,~,nonzero_cl_inds] = eeglab_get_cluster_comps(STUDY);
    CLUSTER_PICKS = main_cl_inds(2:end);
    %- stores
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
    %## generate iter pairs
    cl_inds = repmat(CLUSTER_PICKS,1,length(STUDY.design));
    iter_pairs = zeros(size(cl_inds,2),2);
    c_inds = 1:length(CLUSTER_PICKS);
    for d_i = 1:length(STUDY.design)
        iter_pairs(c_inds,1:2) = [cl_inds(c_inds)',repmat(d_i,length(c_inds),1)];
        c_inds = c_inds + length(c_inds);
    end
    %% (PARFOR) Compute Specs & ERSPs
%     parfor (cnt = 1:size(iter_pairs,1),size(iter_pairs,1))
    for cnt = 1:size(iter_pairs,1)
        TMP_STUDY = STUDY;
        TMP_STUDY.cache = [];
        cluster_i = iter_pairs(cnt,1);
        des_i = iter_pairs(cnt,2);
        TMP_STUDY.currentdesign = des_i;
        design_char = [];
        for i = 1:length(TMP_STUDY.design(des_i).variable)
            if i == 1
                design_char = [TMP_STUDY.design(des_i).variable(1).value{:}];
            else
                design_char = [design_char '_' TMP_STUDY.design(des_i).variable(i).value{:}];
            end
        end
        %- generate power spectrums for each cluster and condition
        [~,spec_savef,spec_subjcorr_savef] = mim_gen_spec_data(TMP_STUDY,ALLEEG,...
                averaged_warpto_events,des_i,cluster_i,design_char,spec_data_dir);
        %- generate event related spectral perturbations for each cluster and condition
        [~,ersp_savef,ersp_subbase_savef,ersp_subbase_combase_savef,ersp_singletrial_subbase_savef] = mim_gen_ersp_data(TMP_STUDY,ALLEEG,averaged_warpto_events,...
                parameters,des_i,cluster_i,design_char,spec_data_dir,...
                'STAT_PARAMS',ERSP_STAT_PARAMS);
        %- store outputs
        spec_fpaths{cnt,1} = spec_savef;
        spec_ss_fpaths{cnt,1} = spec_subjcorr_savef;
        ersp_fpaths{cnt,1}  = ersp_savef;
        ersp_norm_fpaths{cnt,1}  = ersp_subbase_savef;
        ersp_normcb_fpaths{cnt,1}  = ersp_subbase_combase_savef;
        ersp_singtrial_fpaths{cnt,1} = ersp_singletrial_subbase_savef;
        load_ind_cl{cnt,1} = i;
        clust_ind_cl{cnt,1} = cluster_i;
        des_ind{cnt,1} = des_i;
        cnt_ind{cnt,1} = cnt;
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
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}