%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh


%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = fileparts(SCRIPT_DIR);
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(STUDY_DIR);
fprintf(1,'Current folder: %s\n',STUDY_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
FORCE_RECALC_ERSP = true;
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
% dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
% dt = '10052023_MIM_OAN70_noslowwalkers_gait';
% dt = '10252023_MIM_OAN70_noslowwalkers_gait_powpow0p20';
% dt = '10302023_MIM_OAN70_noslowwalkers_gait_iccREMG0p4_powpow0p1';
% dt = '10302023_MIM_OAN70_newnormalize_iccREMG0p4_powpow0p1';
% dt = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
% dt = '11302023_MIM_OAN70_antsnormalize_iccREMG0p3_powpow0p1';
% dt = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
dt = '01232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3';
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
SUB_GROUP_FNAME = [];
SUB_GROUP_FNAME_REGEX = []; 
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},...
%             'variable1','group','values1',{'H2000''s','H3000''s'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
%             'variable1','group','values1',{'H2000''s','H3000''s'}}};
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
    %- assign studies
    if ~isempty(SUB_GROUP_FNAME_REGEX)
        inds = find(cellfun(@(x) strcmp(x,SUB_GROUP_FNAME_REGEX),{STUDY.datasetinfo.group}));
        fprintf('Running subjects:'); cellfun(@(x) fprintf('%s,',x),{STUDY.datasetinfo(inds).subject}); fprintf('\n');
%         for des_i = 1:length(STUDY_DESI_PARAMS)
%             STUDY_DESI_PARAMS{des_i}{2} = {inds}; %{STUDY.datasetinfo(inds).subject};
%         end
        for des_i = 1:length(STUDY_DESI_PARAMS)
            STUDY_DESI_PARAMS{des_i}{10} = {SUB_GROUP_FNAME_REGEX};
        end
    else
        %- combine groups
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_unif{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
        end
    end
    STUDY.cache = [];
    for des_i = 1:length(STUDY_DESI_PARAMS)
        [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
    end
    %## Get Cluster Information
    clust_i = CLUSTER_K_PICKS(k_i);
    fprintf('Making data for K=%i...\n',clust_i);
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    if ~isempty(SUB_GROUP_FNAME_REGEX)
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
    CLUSTER_PICKS = main_cl_inds(2:end); % nonzero_cl_inds; %1:length(TMP_STUDY.cluster);
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
    
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
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
            close all           
        end
    end
    %-
end
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
    CLUSTER_PICKS = main_cl_inds(2:end); 
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