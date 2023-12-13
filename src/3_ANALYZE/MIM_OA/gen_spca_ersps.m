
%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_ersps.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
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
%- define the directory to the src folderd
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
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
WAVELET_STRUCT = struct('t',[0,1/500],...
    'f',(4:100),...
    'fc',1,...
    'FWHM_tc',3,...
    'squared','n');
SPCA_PARAMS = struct('analysis_type','component',...
    'event_char','RHS',...
    'epoch_min_max',[1,4.25],...
    'n_resamples',100,...
    'timewarp_events',{{'RHS','LHS','LTO','RTO'}},...
    'condition_base','rest',...
    'condition_gait',{{'flat','low','med','high','0p25','0p5','0p75','1p0'}});
%- datetime override
% dt = '10252023_MIM_OAN70_noslowwalkers_gait_powpow0p25';
% dt = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
dt = '11302023_MIM_OAN70_antsnormalize_iccREMG0p3_powpow0p1';
% dt = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%## Soft Define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
STUDY_DIR = [STUDIES_DIR filesep sprintf('%s',dt)];
%% ===================================================================== %%
%- (EDIT!) convert SUB_DIR
SUB_DIR = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET

CLUSTER_MAT_FNAME = 'cl_inf_12.mat';
CLUSTER_FPATH = [SUB_DIR filesep 'icrej_5' filesep '12'];
% CLUSTER_STUDY_FPATH = [SUB_DIR filesep 'icrej_5'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
%%
fprintf('\nLoading cluster data from %s\n',CLUSTER_FPATH);
spec_data_dir = [CLUSTER_FPATH filesep 'spec_data'];
plot_store_dir = [CLUSTER_FPATH filesep 'plots_out'];
if ~exist(spec_data_dir,'dir')
    error('spec_data dir does not exist');
end
if ~exist(plot_store_dir,'dir')
    mkdir(plot_store_dir);
end
%## Load Study
if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAME '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '_UNIX.study'],'filepath',spec_data_dir);
    else
        [STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '.study'],'filepath',spec_data_dir);
    end
end
%%
subject_chars = {MAIN_ALLEEG.subject};
clear MAIN_ALLEEG
%%
% for subj_i = 1:length(subject_chars)
parfor (subj_i = 1:length(subject_chars),floor(length(subject_chars)/3))
    %## GET REST CONDITION FROM ORIGINAL ICA
    fpath = [STUDIES_DIR filesep OA_PREP_FPATH filesep subject_chars{subj_i} filesep 'clean'];
    fname = dir([fpath filesep '*.set']);
%     EEG_full = pop_loadset('filepath',fpath,'filename',fname(1).name);
    [~,EEG_full,~] = eeglab_loadICA(fname(1).name,fpath);
    %- remove extremely noise channels (preprocessing)?
    % (12/9/2023) this step may not be necessary since we already performed
    % reject pre-ICA
%     chan_rmv = {MAIN_ALLEEG(subj_i).chaninfo.removedchans.labels};
%     ind_rmv = {MAIN_ALLEEG(subj_i).chaninfo.removedchans.urchan};
%     type_rmv = {MAIN_ALLEEG(subj_i).chaninfo.removedchans.type};
%     inds = strcmp(type_rmv,'EEG');
%     ind_rmv = ind_rmv(inds); ind_rmv = [ind_rmv{:}];
%     EEG_full = pop_select(EEG_full, 'nochannel', chan_rmv);
    %- remove brain components (ic rejection for brain comps)
    % (12/9/2023) JS, this may not be needed, or even detrimental, for the
    % success of the spca alg to remove muscle. Perhaps right flow would be
    % to (1) preproc (2) ic reject (2.5) cluster ics (3) generate ERSP (4) use EEG from after
    % step 1 to generate "clean" ERSPs (5) use step 2 EEG to reject
    % non-brain ics and only calculate averages for brain areas assigned in
    % step 2.5
%     chan_ind_keep = MAIN_ALLEEG(subj_i).etc.urreject.ic_keep;
%     EEG_full = pop_select(EEG_full, 'channel', chan_ind_keep);    
%     EEG_full = pop_subcomp(EEG_full, chan_ind_keep, 0, 1);
    %-
%     EEG_full = eeg_checkset(EEG_full,'loaddata');
%     if isempty(EEG_full.icaact)
%         fprintf('%s) Recalculating ICA activations\n',EEG_full.subject);
%         EEG_full.icaact = (EEG_full.icaweights*EEG_full.icasphere)*EEG_full.data(EEG_full.icachansind,:);
%         EEG_full.icaact = reshape( EEG_full.icaact, size(EEG_full.icaact,1), EEG_full.pnts, EEG_full.trials);
%     end
    %- get rest indicies
    inds1 = logical(strcmp({EEG_full.event.cond}, SPCA_PARAMS.condition_base));
    inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
    val_inds = find(inds1 & ~inds2);
    FROM = [EEG_full.event(val_inds(1)).latency];
    TO = [EEG_full.event(val_inds(end)).latency];
    fprintf('%s) Rest length is %0.2fs\n',subject_chars{subj_i},(TO-FROM)/1000);
    EEG_BASE = pop_select(EEG_full, 'point', [FROM; TO]');
    EEG_BASE = eeg_checkset(EEG_BASE,'loaddata');
    if isempty(EEG_BASE.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG_BASE.subject);
        EEG_BASE.icaact = (EEG_BASE.icaweights*EEG_BASE.icasphere)*EEG_BASE.data(EEG_BASE.icachansind,:);
        EEG_BASE.icaact = reshape( EEG_BASE.icaact, size(EEG_BASE.icaact,1), EEG_BASE.pnts, EEG_BASE.trials);
    end
    %##
    
    %- get gait EEG
    ALLEEG = cell(length(SPCA_PARAMS.condition_gait),1);
    for cond_i = 1:length(SPCA_PARAMS.condition_gait)
        inds1 = logical(strcmp({EEG_full.event.cond}, SPCA_PARAMS.condition_gait{cond_i}));
        inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
        val_inds = find(inds1 & ~inds2);
        FROM = [EEG_full.event(val_inds(1)).latency];
        TO = [EEG_full.event(val_inds(end)).latency];
        ALLEEG{cond_i} = pop_select(EEG_full, 'point', [FROM; TO]');
        % print
        fprintf('\n%s) Condition %s''s length is %0.2fs\n',subject_chars{subj_i},...
            SPCA_PARAMS.condition_gait{cond_i},(TO-FROM)/1000);
    end
    ALLEEG =  cellfun(@(x) [[]; x], ALLEEG);
    ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
    %- Recalculate ICA Matrices && Book Keeping
    ALLEEG = eeg_checkset(ALLEEG,'loaddata');
    if isempty(ALLEEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',ALLEEG.subject);
        ALLEEG.icaact = (ALLEEG.icaweights*ALLEEG.icasphere)*ALLEEG.data(ALLEEG.icachansind,:);
        ALLEEG.icaact = reshape( ALLEEG.icaact, size(ALLEEG.icaact,1), ALLEEG.pnts, ALLEEG.trials);
    end
    %- DEBUG
    %{
    ALLEEG = eeglab_pop_subcomp(ALLEEG,(1:4),0);
    EEG_BASE = eeglab_pop_subcomp(EEG_BASE,(1:4),0);
    %- Recalculate ICA Matrices && Book Keeping
    ALLEEG = eeg_checkset(ALLEEG,'loaddata');
    if isempty(ALLEEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',ALLEEG.subject);
        ALLEEG.icaact = (ALLEEG.icaweights*ALLEEG.icasphere)*ALLEEG.data(ALLEEG.icachansind,:);
        ALLEEG.icaact = reshape( ALLEEG.icaact, size(ALLEEG.icaact,1), ALLEEG.pnts, ALLEEG.trials);
    end
    %- Recalculate ICA Matrices && Book Keeping
    EEG_BASE = eeg_checkset(EEG_BASE,'loaddata');
    if isempty(EEG_BASE.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG_BASE.subject);
        EEG_BASE.icaact = (EEG_BASE.icaweights*EEG_BASE.icasphere)*EEG_BASE.data(EEG_BASE.icachansind,:);
        EEG_BASE.icaact = reshape( EEG_BASE.icaact, size(EEG_BASE.icaact,1), EEG_BASE.pnts, EEG_BASE.trials);
    end
    %}
    %## SPCA DENOISE ALGORITHM
    fprintf('\nRunning sPCA Algorithm...\n')
    tt = tic;
    ALLEEG.data = [];
    EEG_BASE.data = [];
    [gait_avg,ERSP,GPM,~,output_struct] = spca_time_freq_decomp(ALLEEG,EEG_BASE,...
        'SPCA_PARAMS',SPCA_PARAMS,...
        'WAVELET_STRUCT',WAVELET_STRUCT);
    [ERSP_corr, GPM_corr, PSC1, ~,COEFFs] = specPCAdenoising(ERSP);
    delete ALLEEG EEG_BASE
    fprintf('time: %0.2f',toc(tt));
    %## SAVE PCA INFORMATION
    gait_ersp_struct = [];
    gait_ersp_struct.ID         = EEG_gait.subject;
    gait_ersp_struct.Noise_cov  = output_struct.baseline_cov;% noise cov for kernel computation
    gait_ersp_struct.F_Rest     = output_struct.baseline_ersp;
    gait_ersp_struct.TF         = gait_avg;
    gait_ersp_struct.ERSP_uncor = ERSP;
    gait_ersp_struct.GPM_uncor  = GPM;
    gait_ersp_struct.ERSP       = ERSP_corr;
    gait_ersp_struct.GPM        = GPM_corr;
    gait_ersp_struct.PSC1       = PSC1;
    gait_ersp_struct.numStrides = output_struct.cycle_cnt;
    gait_ersp_struct.numValidStrides = output_struct.valid_cycle_cnt;
    gait_ersp_struct.chanlocs   = EEG_gait.chanlocs;
    par_save(gait_ersp_struct,fpath,'gait_ersp_spca.mat');
    %## PLOT
    figure(); set(gcf, 'position', [0 0 600 500]);
    plot(WAVELET_STRUCT.f, squeeze(output_struct.baseline_ersp)', 'k-');
    ylabel('Amplitude (\muV)'), ylabel('Frequency (Hz)');
    grid on; box off
    title('Baseline ERSP (rest)');
    exportgraphics(fig,[fpath filesep 'allcond_baseline_avgs.jpg']);
    %##
    CHAN_INT = randi(size(ERSP,2),1);
    fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
    exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_corr.jpg',CHAN_INT)]);
    %##
    fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
    exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_corr.jpg',CHAN_INT)]);
    %##
    fig = plot_txf(squeeze(ERSP(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
    exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_orig.jpg',CHAN_INT)]);
    %##
    fig = plot_txf(squeeze(GPM(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
    exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_orig.jpg',CHAN_INT)]);
    %##
    fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
    exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
    
    %% DENOISE
    fprintf('\nRunning Denoising with sPCA\n');
    %- get gait EEG
    ALLEEG = cell(length(SPCA_PARAMS.condition_gait),1);
    %- calculate baseline mean
%     [tf_rest,noise_cov] = txf_baseline(EEG_BASE,SPCA_PARAMS.analysis_type,WAVELET_STRUCT);
%     base_mean(1,:,:) = squeeze(mean(abs(tf_rest(10+1:end-10,:,:)),1));
    %- use already calculated baselinen mean
%     base_mean(1,:,:) = squeeze(mean(abs(output_struct.baseline_ersp(10+1:end-10,:,:)),1));
    base_mean = output_struct.baseline_ersp;
    CHAN_INT = randi(size(ERSP,2),1);
    for cond_i = 1:length(SPCA_PARAMS.condition_gait)
        inds1 = logical(strcmp({EEG_full.event.cond}, SPCA_PARAMS.condition_gait{cond_i}));
        inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
        val_inds = find(inds1 & ~inds2);
        FROM = [EEG_full.event(val_inds(1)).latency];
        TO = [EEG_full.event(val_inds(end)).latency];
        EEG_gait = pop_select(EEG_full, 'point', [FROM; TO]');
        %- print
        fprintf('\n%s) Condition %s''s length is %0.2fs\n',subject_chars{subj_i},...
            SPCA_PARAMS.condition_gait{cond_i},(TO-FROM)/1000);
        %## DEBUG
%         EEG_gait = eeglab_pop_subcomp(EEG_gait,(1:4),0);
        %-
        if isempty(EEG_gait.icaact)
            EEG_gait = eeg_checkset(EEG_gait,'loaddata');
            fprintf('%s) Recalculating ICA activations\n',EEG_gait.subject);
            EEG_gait.icaact = (EEG_gait.icaweights*EEG_gait.icasphere)*EEG_gait.data(EEG_gait.icachansind,:);
            EEG_gait.icaact = reshape( EEG_gait.icaact, size(EEG_gait.icaact,1), EEG_gait.pnts, EEG_gait.trials);
        end
        EEG_gait.data = [];
        data = permute(EEG_gait.icaact, [2,1]); % pnts x chans! --> BS way?
        n_comps = size(data, 2);
        n_freqs = length(WAVELET_STRUCT.f);
        %-
        hs_min_max = SPCA_PARAMS.epoch_min_max*EEG_full.srate;
        data = bsxfun(@minus, data, mean(data,2));
        %- time frequency transform
        [TF, morlet_params] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared);
        TF = abs(TF);
        %-
        idx_hs = find(strcmp({EEG_gait.event.type}, SPCA_PARAMS.event_char));
        %## METHOD 1
        %{
        gait_ersp_out = zeros(output_struct.cycle_cnt,SPCA_PARAMS.n_resamples,n_comps,n_freqs); %strides/trials x pnts x chans x freqs
        gait_gpm_out = zeros(output_struct.cycle_cnt,SPCA_PARAMS.n_resamples,n_comps,n_freqs); %strides/trials x pnts x chans x freqs
        %}
        %## METHOD 2
        gait_tf = zeros(length(idx_hs)-1,SPCA_PARAMS.n_resamples,n_comps*n_freqs); %strides/trials x pnts x chans x freqs
        %- step counter, increased for each valid step
        iter = 1;
        cnt = 1;
        %- resample each stride to the same legth (100 pnts)
        for cycle_cnt = 1:length(idx_hs)-1
            %- find first and last sample of stride
            cycle_edge = round([EEG_gait.event(idx_hs(cycle_cnt)).latency,...
                EEG_gait.event(idx_hs(cycle_cnt+1)).latency-1]); % first and last frame of gait cycle
            %- labels of all events within this cycle
            cycle_event = {EEG_gait.event([idx_hs(cycle_cnt):idx_hs(cycle_cnt+1)]).type};
            %- only keep labels of gait events to check their order:
            cycle_gaitEvent = cycle_event(contains(cycle_event,SPCA_PARAMS.timewarp_events));
            %-
            if hs_min_max(1) <= cycle_edge(2)-cycle_edge(1) &&... % check time until next HS
                    cycle_edge(2)-cycle_edge(1) <= hs_min_max(2) &&...
                    all(ismember(SPCA_PARAMS.timewarp_events,cycle_gaitEvent)) % oder of gait events correct
                %## METHOD 1
                %{
                tf_cycle = TF_new(cycle_edge(1):cycle_edge(2),:,:); % extract data
                tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*N_FREQS); % reshape to be able to use the resample function, skip but resample over different dimension?
                gait_tf = resample(tf_cycle,SPCA_PARAMS.n_resamples,cycle_edge(2)-cycle_edge(1)+1,0);
                gait_tf = reshape(gait_tf,SPCA_PARAMS.n_resamples,n_comps,n_freqs);
                gait_tf = 20*bsxfun(@minus,log10(gait_tf), log10(base_mean));
                [ERSP_corr, GPM_corr, PSC1,~,~] = specPCAdenoising(gait_tf,COEFFs);
                gait_ersp_out(cnt,:,:,:) = ERSP_corr;
                gait_gpm_out(cnt,:,:,:) = GPM_corr;
                %##
                fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_trial%i_ersp_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i},cnt)]);
                %##
                fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_trial%i_gpm_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i},cnt)]);
                %}
                %## METHOD 2
                tf_cycle = TF(cycle_edge(1):cycle_edge(2),:,:); % extract data
                tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*n_freqs); % reshape to be able to use the resample function, skip but resample over different dimension?
                gait_tf(cnt,:,:) = resample(tf_cycle,SPCA_PARAMS.n_resamples,cycle_edge(2)-cycle_edge(1)+1,0); % resample and store
                cnt = cnt+1;
            end
        end
        %## METHOD 1
        %{
        gait_avg = squeeze(mean(gait_ersp_out)); 
        gait_gpm_avg = squeeze(mean(gait_gpm_out));
        par_save(gait_avg,fpath,sprintf('cond%s_spca_ersp_corr.mat',SPCA_PARAMS.condition_gait{cond_i}));
        par_save(gait_gpm_avg,fpath,sprintf('cond%s_spca_gpm_corr.mat',SPCA_PARAMS.condition_gait{cond_i}));
        par_save(gait_ersp_out,fpath,sprintf('cond%s_spca_ersp_alltrials.mat',SPCA_PARAMS.condition_gait{cond_i}));
        par_save(gait_gpm_out,fpath,sprintf('cond%s_spca_gpm_alltrials.mat',SPCA_PARAMS.condition_gait{cond_i}));
        %}
        %## METHOD 2
        gait_tf = reshape(gait_tf,size(gait_tf,1),SPCA_PARAMS.n_resamples,n_comps,n_freqs); % reshape to trials x pnts x chans x freqs
        gait_avg = squeeze(mean(gait_tf)); % average over trials
        %-
        ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_mean));
        [ERSP_corr, GPM_corr, PSC1,~,~] = specPCAdenoising(ERDS,COEFFs);
        fprintf('\n%s) Plotting validations...\n',subject_chars{subj_i});
        disp([num2str(round(cnt/cycle_cnt*100)) '% of the gait cycles are valid'])
        par_save(ERSP_corr,fpath,sprintf('cond%s_spca_ersp.mat',SPCA_PARAMS.condition_gait{cond_i}));
        %##
        fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
        title(sprintf('%s) ERSP corrected',SPCA_PARAMS.condition_gait{cond_i}));
        exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
        %##
        fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
        title(sprintf('%s) GPM corrected',SPCA_PARAMS.condition_gait{cond_i}));
        exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_gpm_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
        %##
        fig = plot_txf(squeeze(ERDS(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
        title(sprintf('%s) ERSP original',SPCA_PARAMS.condition_gait{cond_i}));
        exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
        %##
%         fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
%         title(sprintf('%s) PSC1 original',SPCA_PARAMS.condition_gait{cond_i}));
%         exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_psc1_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);

    end
end
%%