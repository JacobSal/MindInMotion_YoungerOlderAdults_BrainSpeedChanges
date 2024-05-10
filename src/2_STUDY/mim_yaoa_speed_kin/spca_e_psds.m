%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_spca_e_psds.sh

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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%% Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name & 
DATA_SET = 'MIM_dataset';
PREPROC_DIR_NAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% study_dir = '01122024_spca_analysis';
% study_dir = '02082024_spca_analysis_fix';
STUDY_DIR_NAME = '03232024_spca_analysis_OA';
STUDY_FNAME_REST = 'rest_epoch_study';
STUDY_FNAME_GAIT = 'epoch_study';
%- study group and saving
SAVE_ALLEEG = false;
% %- epoching params
% DO_SLIDING_WINDOW = false;
% %* sliding window
% WINDOW_LENGTH = 6; % sliding window length in seconds
% PERCENT_OVERLAP = 0.0; % percent overlap between epochs
% %* gai
% EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
% STD_TIMEWARP = 3;
% EPOCH_TIME_LIMITS = [-0.5,4.5];
% TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
% if DO_SLIDING_WINDOW
%     SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
%     TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
% else
%     SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
%     TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% end
%## EPOCH PARAMS
DEF_EPOCH_PARAMS = struct('epoch_method','timewarp',...
    'percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'baseline_time_lims',[-0.5,4.5-2],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'cond_field','cond',...
    'appx_cond_len',3*60,...
    'slide_cond_chars',{{}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%- compute measures for spectrum and ersp
ICLABEL_EYE_CUTOFF = 0.75;
FORCE_RECALC_PSD = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false;
DO_RECREATE_SPCA = true;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
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
SPCA_PARAMS = struct('analysis_type','component',...
    'event_char','RHS',...
    'epoch_min_max',[1,4.25],...
    'n_resamples',100,...
    'timewarp_events',{{'RHS','LHS','LTO','RTO'}},...
    'condition_base','rest',...
    'condition_gait',{{'flat','low','med','high','0p25','0p5','0p75','1p0'}});
%% (PATHS)
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',STUDY_DIR_NAME)];
ica_dir = [STUDIES_DIR filesep PREPROC_DIR_NAME];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '.study'],'filepath',save_dir);
% end
% if ~ispc
%     [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '.study'],'filepath',save_dir);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    STUDY = tmp.STUDY;
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
end
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
%## ersp plot per cluster per condition
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);

SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%%
%## ersp plot per cluster per condition
STUDY_REST = pop_statparams(STUDY_REST,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);

SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_REST = pop_specparams(STUDY_REST,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
subject_chars = {STUDY_REST.datasetinfo.subject};
% for subj_i = 1:length(ALLEEG)
parfor subj_i = 1:length(STUDY_REST.datasetinfo)
    tt = tic;
    base_txf_mean = [];
    % EEG = ALLEEG(subj_i);
    % EEG = pop_loadset('filepath',[save_dir filesep 'GAIT_EPOCHED' filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}]],...
    %     'filename',sprintf('%s_%s',[DEF_EPOCH_PARAMS.gait_trial_chars{:}]))
    EEG = pop_loadset('filepath',STUDY.datasetinfo(subj_i).filepath,...
        'filename',STUDY.datasetinfo(subj_i).filename);
    fprintf('Running Subject %s\n',EEG.subject);
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    %##
    if DO_RECREATE_SPCA || ~all(cellfun(@(x) exist([EEG.filepath filesep sprintf('cond%s_spca_psd.mat',x)],'file'),DEF_EPOCH_PARAMS.gait_trial_chars))
        %## GENERATE ERSPS
        icaspec_f = [EEG.filepath filesep sprintf('%s.icaspec',EEG.subject)];
        if ~exist(icaspec_f,'file') || FORCE_RECALC_PSD % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            TMP_STUDY = STUDY;
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
            fprintf('Done calculating PSDs for %s\n',EEG.subject);
        else
            fprintf('PSDs already calculated for %s\n',EEG.subject);
        end

        %##
        % EEG = ALLEEG_REST(subj_i);
        EEG = pop_loadset('filepath',STUDY_REST.datasetinfo(subj_i).filepath,...
            'filename',STUDY_REST.datasetinfo(subj_i).filename);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        EEG = eeglab_makecontinuous(EEG);
        icaspec_f = [EEG.filepath filesep sprintf('%s.icaspec',EEG.subject)];
        if ~exist(icaspec_f,'file') || FORCE_RECALC_PSD % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            TMP_STUDY = STUDY_REST;
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
            fprintf('Done calculating PSDs for %s\n',EEG.subject);
        else
            fprintf('PSDs already calculated for %s\n',EEG.subject);
        end

        try
            %## (LOAD REST) ================================================ %%
            epoched_fPath = strsplit(EEG.filepath,filesep);
            epoched_fPath = strjoin(epoched_fPath(1:end-1),filesep);
            fpath = [epoched_fPath filesep 'rest'];
    %         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,'rest');
            icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
            %- load .icaspec load-in parameters
            fprintf('Loading Resting Frequency Data...\n');
            tmp = load(icaspec_f,'-mat');
    %         parameters = tmp.parameters;
            %- reshape data [pnts x chans x freq]
            fn = fieldnames(tmp);
            inds = find(contains(fieldnames(tmp),'comp'));
            inds_not = find(~contains(fieldnames(tmp),'comp'));
            tmp_out = [];
            for i = 1:length(inds_not)
                tmp_out.(fn{inds_not(i)}) = tmp.(fn{inds_not(i)});
            end
            test = tmp.(fn{inds(1)});
    %         rest_txf = zeros(size(test,3),size(test,2),size(test,1),length(inds),'single');
            rest_psdf = zeros(size(test,1),size(test,2),length(inds),'double');
            for i = 1:length(inds)
                rest_psdf(:,:,i) = tmp.(fn{inds(i)}); % freq x epoch x chan
            end
            % rest_psdf = squeeze(mean(mean(abs(rest_psdf),3),2));
            rest_psdf = permute(rest_psdf,[2,1,3]);
            %- average over time, keep magnitude (not power -> would amplify outliers)
            base_txf_mean = mean(rest_psdf);
            base_txf_mean = permute(base_txf_mean,[1,3,2]);

            % base_txf_mean = zeros(1,size(rest_psdf,1),size(rest_psdf,2));
            % base_txf_mean(1,:,:) = rest_psdf; % format this as 'pnts' x chans x freq
    %         base_txf_mean(1,:,:) = squeeze(mean(abs(rest_txf(10+1:end-10,:,:)),1));
            %- clear rest_txf
            % psd_out = struct('psd_fec',rest_psdf)
            rest_psdf = double.empty;
            par_save(tmp,fpath,sprintf('rest_avg_psdf.mat'));
            %## (LOAD GAIT TIMEWARP) ======================================= %%
            fprintf('Loading Gait Frequency Data...\n');
            fpath = [epoched_fPath filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}]];
    %         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[TRIAL_TYPES{:}]);
            icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
            %- load .icaspec load-in parameters
            tmp = load(icaspec_f,'-mat');
    %         parameters = tmp.parameters;
            %## (APPLY SPCA) =============================================== %%
            fprintf('Running SPCA on All Conditions\n');
            %- sPCA Algorithm
            [psd,psd_avg,output_struct] = apply_spca_cond_psd(tmp,base_txf_mean);
            [psd_corr,based_psd_corr,PSC1,~,COEFFS] = specPCAdenoising(psd);
            
            %## SAVE PCA INFORMATION
            gait_ersp_struct = [];
            gait_ersp_struct.ID         = EEG.subject;
            gait_ersp_struct.Noise_cov  = [];% noise cov for kernel computation
            gait_ersp_struct.F_Rest     = output_struct.baseline_ersp;
            gait_ersp_struct.psd_avg         = psd_avg;
            gait_ersp_struct.psd = psd;
            gait_ersp_struct.psd_corr       = psd_corr;
            gait_ersp_struct.based_psd_corr        = based_psd_corr;
            gait_ersp_struct.PSC1       = PSC1;
            gait_ersp_struct.chanlocs   = EEG.chanlocs;
            gait_ersp_struct.icatimefopts = tmp_out;
            % gait_ersp_struct.warptimes  = averaged_warpto_events;
            % gait_ersp_struct.ntimes     = TIMEWARP_NTIMES;
            par_save(gait_ersp_struct,fpath,'gait_psd_spca.mat');
            %## PLOT
            fig = figure(); set(gcf, 'position', [0 0 600 500]);
    %         plot(tmp.freqs, squeeze(output_struct.baseline_ersp)', 'k-');
            plot(tmp.freqs, squeeze(base_txf_mean)', 'k-');
            ylabel('Amplitude (\muV)');
            xlabel('Frequency (Hz)');
            grid on; box off
            title('Baseline ERSP (rest)');
            exportgraphics(fig,[fpath filesep 'allcond_baseline_avgs_psds.jpg']);

            %## VALIDATION PLOTS
            %{
            CHAN_INT = randi(size(ERSP,2),1);
            fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],tmp.freqs);
            exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_corr.jpg',CHAN_INT)]);
            %-
            fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],tmp.freqs);
            exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_corr.jpg',CHAN_INT)]);
            %-
            fig = plot_txf(squeeze(ERSP(:,CHAN_INT,:)),[],tmp.freqs);
            exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_orig.jpg',CHAN_INT)]);
            %-
            fig = plot_txf(squeeze(GPM(:,CHAN_INT,:)),[],tmp.freqs);
            exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_orig.jpg',CHAN_INT)]);
            %-
            fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],tmp.freqs);
            exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
            %}
            %% DENOISE
            fprintf('\nRunning Denoising with sPCA\n');
            conds = unique({tmp.trialinfo.cond});
            for cond_i = 1:length(conds)
                COND_STR = conds{cond_i};
                [psd_corr,psd_avg,output_struct] = apply_spca_cond_psd(tmp,base_txf_mean,...
                    'COEFFS',COEFFS,...
                    'COND_STR',COND_STR);
                struct_out = struct('psd_corr',psd_corr,...
                    'freqs',tmp.freqs,...
                    'pc1',psd_avg,...
                    'coeffs',COEFFS,...
                    'apply_spca_cond',output_struct);
                par_save(struct_out,fpath,sprintf('cond%s_spca_psd.mat',conds{cond_i}));

                %## VALIDATION PLOTS
                %{
                fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],tmp.freqs);
                title(sprintf('%s) ERSP corrected',conds{cond_i}));
                if ~isempty(fpath)
                    exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_corr.jpg',CHAN_INT,conds{cond_i})]);
                end
                %-
                fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],tmp.freqs);
                title(sprintf('%s) GPM corrected',conds{cond_i}));
                if ~isempty(fpath)
                    exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_gpm_corr.jpg',CHAN_INT,conds{cond_i})]);
                end
                %-
                fig = plot_txf(squeeze(output_struct.erds_orig(:,CHAN_INT,:)),[],tmp.freqs);
                title(sprintf('%s) ERSP original',conds{cond_i}));
                if ~isempty(fpath)
                    exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,conds{cond_i})]);
                end
                %-
                fig = plot_txf(squeeze(output_struct.gpm_orig(:,CHAN_INT,:)),[],tmp.freqs);
                title(sprintf('%s) GPM original',conds{cond_i}));
                if ~isempty(fpath)
                    exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_gpm_orig.jpg',CHAN_INT,conds{cond_i})]);
                end
                %##
                % fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],tmp.freqs);
                % title(sprintf('%s) PSC1 original',SPCA_PARAMS.condition_gait{cond_i}));
                % exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_psc1_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
                %}
            end
            %## CLEANUP AND FREE UP DRIVE SPACE
            fpath = [epoched_fPath filesep [DEF_EPOCH_PARAMS.gait_trial_chars{:}]];
            icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
            delete(icaspec_f);
            fpath = [epoched_fPath filesep 'rest'];
            icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
            delete(icaspec_f);
            fprintf('Subject %s done. time: %0.2f',EEG.subject, toc(tt));
        catch e
            fprintf('\nError occured on subject %s\n%s\n',subject_chars{subj_i},getReport(e));
        end
    end
end