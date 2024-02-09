
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
    'f',(1:200),...
    'fc',1,...
    'FWHM_tc',3,...
    'squared','n',...
    'data_type','double');
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
% dt = '11302023_MIM_OAN70_antsnormalize_iccREMG0p3_powpow0p1';
dt = '01122024_spca_analysis_subsample';
% dt = 'spca_analysis';
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%## Soft Define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
STUDY_DIR = [STUDIES_DIR filesep sprintf('%s',dt)];
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa');
subject_chars = {};
for i = 1:length(SUBJ_PICS)
    subject_chars = [subject_chars, SUBJ_PICS{i}];
end
%%
% for subj_i = 1:length(subject_chars)
parfor (subj_i = 1:length(subject_chars),floor(length(subject_chars)/3))
% for subj_i = find(strcmp(subject_chars,'NH3108'))
    %## GET REST CONDITION FROM ORIGINAL ICA
    save_dir_tmp = [STUDY_DIR filesep subject_chars{subj_i}];
    if ~exist(save_dir_tmp,'dir')
        mkdir(save_dir_tmp);
    end
    chk = ~cellfun(@(x) exist([save_dir_tmp filesep sprintf('cond%s_spca_ersp.mat',x)],'file'),SPCA_PARAMS.condition_gait);
    if any(chk)
        try
            base_mean = [];
            n_comps = [];

            spca_params_tmp = SPCA_PARAMS;
            EEG_GAIT = []; %#ok<NASGU>
            fpath = [STUDIES_DIR filesep OA_PREP_FPATH filesep subject_chars{subj_i} filesep 'clean'];
            fname = dir([fpath filesep '*.set']);
            [~,EEG_full,~] = eeglab_loadICA(fname(1).name,fpath);
            %## get rest indicies & txf
            inds1 = logical(strcmp({EEG_full.event.cond}, spca_params_tmp.condition_base));
            inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
            val_inds = find(inds1); %find(inds1 & ~inds2);
            FROM = [EEG_full.event(val_inds(1)).latency];
            TO = [EEG_full.event(val_inds(end)).latency];
            fprintf('%s) Rest length is %0.2fs\n',subject_chars{subj_i},(TO-FROM)/1000);
            EEG_BASE = pop_select(EEG_full, 'point', [FROM; TO]');

%             pop_saveset(EEG_BASE,'filepath',save_dir_tmp,...
%                 'filename',sprintf('sitting_rest.set'),...
%                 'savemode','twofiles');
            %- time frequency transform
            fprintf('\nComputing baseline time-frequency decomposition...\n');
            [tf,~] = eeg_txf_decomp(EEG_BASE,'component','WAVELET_STRUCT',WAVELET_STRUCT);
            %- average over time, keep magnitude (not power -> would amplify outliers)
            base_mean(1,:,:) = squeeze(mean(abs(tf(10+1:end-10,:,:)),1));
            %- clear vars
            tf = double.empty;
            EEG_BASE = struct.empty;
            %## get gait EEG
            EEG_GAIT = cell(length(spca_params_tmp.condition_gait),1);
            for cond_i = 1:length(spca_params_tmp.condition_gait)
                inds1 = logical(strcmp({EEG_full.event.cond}, spca_params_tmp.condition_gait{cond_i}));
                inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
                val_inds = find(inds1 & ~inds2);
                FROM = [EEG_full.event(val_inds(1)).latency];
                TO = [EEG_full.event(val_inds(end)).latency];
                EEG_GAIT{cond_i} = pop_select(EEG_full, 'point', [FROM; TO]');
                % print
                fprintf('\n%s) Condition %s''s length is %0.2fs\n',subject_chars{subj_i},...
                    spca_params_tmp.condition_gait{cond_i},(TO-FROM)/1000);
            end
            EEG_GAIT =  cellfun(@(x) [[]; x], EEG_GAIT);
            EEG_GAIT = pop_mergeset(EEG_GAIT,1:length(EEG_GAIT),1);
            %- DEBUG
            %{
            EEG_full = eeglab_pop_subcomp(ALLEEG,(1:4),0);
            ALLEEG = eeglab_pop_subcomp(ALLEEG,(1:4),0);
            EEG_BASE = eeglab_pop_subcomp(EEG_BASE,(1:4),0);
            %- Recalculate ICA Matrices && Book Keeping
            EEG_full = eeg_checkset(EEG_full,'loaddata');
            if isempty(EEG_full.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG_full.subject);
                EEG_full.icaact = (EEG_full.icaweights*EEG_full.icasphere)*EEG_full.data(EEG_full.icachansind,:);
                EEG_full.icaact = reshape( EEG_full.icaact, size(EEG_full.icaact,1), EEG_full.pnts, EEG_full.trials);
            end
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
            [gait_avg,ERSP,GPM,~,output_struct] = prep_gait_spca(EEG_GAIT,base_mean,...
                'SPCA_PARAMS',spca_params_tmp,...
                'WAVELET_STRUCT',WAVELET_STRUCT);
            %- sPCA Algorithm
            [ERSP_corr, GPM_corr, PSC1, ~,COEFFS] = specPCAdenoising(ERSP);
            fprintf('time: %0.2f',toc(tt));
            %## SAVE PCA INFORMATION
            gait_ersp_struct = [];
            gait_ersp_struct.ID         = EEG_GAIT.subject;
            gait_ersp_struct.Noise_cov  = output_struct.baseline_cov;% noise cov for kernel computation
            gait_ersp_struct.F_Rest     = base_mean;
            gait_ersp_struct.TF         = gait_avg;
            gait_ersp_struct.ERSP_uncor = ERSP;
            gait_ersp_struct.GPM_uncor  = GPM;
            gait_ersp_struct.ERSP       = ERSP_corr;
            gait_ersp_struct.GPM        = GPM_corr;
            gait_ersp_struct.PSC1       = PSC1;
            gait_ersp_struct.numStrides = output_struct.cycle_cnt;
            gait_ersp_struct.numValidStrides = output_struct.valid_cycle_cnt;
            gait_ersp_struct.chanlocs   = EEG_GAIT.chanlocs;
            par_save(gait_ersp_struct,save_dir_tmp,'gait_ersp_spca.mat');
            %- clear vars
            gait_ersp_struct = struct.empty;
            EEG_GAIT = struct.empty;
            %## PLOT
            fig = figure(); set(gcf, 'position', [0 0 600 500]);
            plot(WAVELET_STRUCT.f, squeeze(base_mean)', 'k-');
            ylabel('Amplitude (\muV)'), ylabel('Frequency (Hz)');
            grid on; box off
            title('Baseline ERSP (rest)');
            exportgraphics(fig,[save_dir_tmp filesep 'allcond_baseline_avgs.jpg']);
            %##
            CHAN_INT = randi(size(ERSP,2),1);
            fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
            exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_allcond_ersp_corr.jpg',CHAN_INT)]);
            close(fig);
            %##
            fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
            exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_allcond_gpm_corr.jpg',CHAN_INT)]);
            close(fig);
            %##
            fig = plot_txf(squeeze(ERSP(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
            exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_allcond_ersp_orig.jpg',CHAN_INT)]);
            close(fig);
            %##
            fig = plot_txf(squeeze(GPM(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
            exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_allcond_gpm_orig.jpg',CHAN_INT)]);
            close(fig);
            %##
            fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
            exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
            close(fig);
            %% DENOISE
            fprintf('\nRunning Denoising with sPCA\n');
            %- get gait EEG
        %     ALLEEG = cell(length(SPCA_PARAMS.condition_gait),1);
            %- calculate baseline mean
        %     [tf_rest,noise_cov] = txf_baseline(EEG_BASE,SPCA_PARAMS.analysis_type,WAVELET_STRUCT);
        %     base_mean(1,:,:) = squeeze(mean(abs(tf_rest(10+1:end-10,:,:)),1));
            %- use already calculated baselinen mean
        %     base_mean(1,:,:) = squeeze(mean(abs(output_struct.baseline_ersp(10+1:end-10,:,:)),1));
        %     base_mean = output_struct.baseline_ersp;
            for cond_i = 1:length(spca_params_tmp.condition_gait)
                inds1 = logical(strcmp({EEG_full.event.cond}, spca_params_tmp.condition_gait{cond_i}));
                inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
                val_inds = find(inds1 & ~inds2);
                FROM = [EEG_full.event(val_inds(1)).latency];
                TO = [EEG_full.event(val_inds(end)).latency];
                EEG_gait = pop_select(EEG_full, 'point', [FROM; TO]');
                %- print
                fprintf('\n%s) Condition %s''s length is %0.2fs\n',subject_chars{subj_i},...
                    spca_params_tmp.condition_gait{cond_i},(TO-FROM)/1000);
                [ERSP_corr,GPM_corr,output_struct] = apply_spca_cond(EEG_gait,base_mean,COEFFS);
                struct_out = struct('ersp_corr',ERSP_corr,...
                    'gpm_corr',GPM_corr,...
                    'coeffs',COEFFS,...
                    'apply_spca_cond',output_struct);
                EEG_gait = struct.empty;
                par_save(struct_out,save_dir_tmp,sprintf('cond%s_spca_ersp.mat',SPCA_PARAMS.condition_gait{cond_i}));

                %##
                fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                title(sprintf('%s) ERSP corrected',SPCA_PARAMS.condition_gait{cond_i}));
                if ~isempty(save_dir_tmp)
                    exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_cond%s_ersp_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
                end
                close(fig);
                %##
                fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                title(sprintf('%s) GPM corrected',SPCA_PARAMS.condition_gait{cond_i}));
                if ~isempty(save_dir_tmp)
                    exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_cond%s_gpm_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
                end
                close(fig);
                %##
                fig = plot_txf(squeeze(output_struct.erds_orig(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                title(sprintf('%s) ERSP original',SPCA_PARAMS.condition_gait{cond_i}));
                if ~isempty(save_dir_tmp)
                    exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
                end
                close(fig);
                %##
                fig = plot_txf(squeeze(output_struct.gpm_orig(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                title(sprintf('%s) GPM original',SPCA_PARAMS.condition_gait{cond_i}));
                if ~isempty(save_dir_tmp)
                    exportgraphics(fig,[save_dir_tmp filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
                end
                close(fig);
                %##
                % fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
                % title(sprintf('%s) PSC1 original',SPCA_PARAMS.condition_gait{cond_i}));
                % exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_psc1_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            end
        catch e
            fprintf('\nError occured on subject %s\n%s\n',subject_chars{subj_i},getReport(e));
        end
    else
        fprintf('Non-timewarped SPCA Ersps already %s', subject_chars{subj_i})
    end
end
%%