%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_spca_b_tw_timefreq.sh

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
        STUDY_DIR = SCRIPT_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = SCRIPT_DIR;
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(SCRIPT_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name & 
DATA_SET = 'MIM_dataset';
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% OA_PREP_FPATH = '01132024_antsnorm_iccREEG0p65_iccREMG0p4_skull0p0042';
% study_dir = '01122024_spca_analysis';
study_dir = '03232024_spca_analysis_OA';
% study_fname_rest = 'rest_epoch_study';
% study_fname_gait = 'epoch_study';
study_fname_rest = 'rest_epoch_study_ya';
study_fname_gait = 'epoch_study_ya';
%- study group and saving
SAVE_ALLEEG = false;
%- epoching params
DO_SLIDING_WINDOW = false;
%* sliding window
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
%* gai
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-0.5,4.5];
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
if DO_SLIDING_WINDOW
    SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
    TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
else
    SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
    TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
end
%- compute measures for spectrum and ersp
ICLABEL_EYE_CUTOFF = 0.75;
FORCE_RECALC_ERSP = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false;
DO_RECREATE_SPCA = false;
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
DATA_DIR = [PATHS.src_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [STUDIES_DIR filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '_UNIX.study'],'filepath',save_dir);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '.study'],'filepath',save_dir);
end
if ~ispc
    [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '_UNIX.study'],'filepath',save_dir);
else
    [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '.study'],'filepath',save_dir);
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
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
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
%%
%## ersp plot per cluster per condition
STUDY_REST = pop_statparams(STUDY_REST,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY_REST = pop_erspparams(STUDY_REST,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_REST = pop_specparams(STUDY_REST,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
subject_chars = {ALLEEG.subject};
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
for subj_i = 1:length(ALLEEG)
    
    tt = tic;
    base_txf_mean = [];
    EEG = ALLEEG(subj_i);
    fprintf('Running Subject %s\n',EEG.subject);
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    %##
    if DO_RECREATE_SPCA || ~all(cellfun(@(x) exist([EEG.filepath filesep sprintf('cond%s_spca_ersp.mat',x)],'file'),TRIAL_TYPES))
        %## GENERATE ERSPS
        icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
        if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            TMP_STUDY = STUDY;
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
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                        'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                        'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                        'trialbase','off','basenorm','on'}); %ERSP
            else
                % No baseline correction
                [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
                        'baseline',nan(),'timewarp',timewarp_param,...
                        'timewarpms',timewarpms_param}); %ERSP
            end
            fprintf('Done calculating timewarped ERSPs for %s\n',EEG.subject);
        else
            fprintf('Timewarped ERSPs already calculated for %s\n',EEG.subject);
        end

        %##
        EEG = ALLEEG_REST(subj_i);
        icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
        if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            TMP_STUDY = STUDY_REST;
            EEG = eeg_checkset(EEG,'loaddata');
            if isempty(EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
            end
            %- overrride datasetinfo to trick std_precomp to run.
            TMP_STUDY.datasetinfo = STUDY_REST.datasetinfo(subj_i);
            TMP_STUDY.datasetinfo(1).index = 1;
            %-
            if DO_BASELINE_CORRECTION
                % Baseline correction
                [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                        'ntimesout',TIMEWARP_NTIMES,...
                        'trialbase','off','basenorm','on'}); %ERSP
            else
                % No baseline correction
                [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES}); %ERSP
            end
            fprintf('Done calculating timewarped ERSPs for %s\n',EEG.subject);
        else
            fprintf('Timewarped ERSPs already calculated for %s\n',EEG.subject);
        end

        try
            %## (LOAD REST) ================================================ %%
            epoched_fPath = strsplit(EEG.filepath,filesep);
            epoched_fPath = strjoin(epoched_fPath(1:end-1),filesep);
            fpath = [epoched_fPath filesep 'rest'];
    %         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,'rest');
            icatimef_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
            %- load .icatimef load-in parameters
            fprintf('Loading Resting Time-Frequency Data...\n');
            tmp = load(icatimef_f,'-mat');
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
            rest_txf = zeros(size(test,1),size(test,2),size(test,3),length(inds),'double');
            for i = 1:length(inds)
                rest_txf(:,:,:,i) = tmp.(fn{inds(i)}); % freq x time x epoch x chan
            end
            rest_txf = squeeze(mean(mean(abs(rest_txf),3),2));
            rest_txf = permute(rest_txf,[2,1]);
            %- average over time, keep magnitude (not power -> would amplify outliers)
            base_txf_mean = zeros(1,size(rest_txf,1),size(rest_txf,2));
            base_txf_mean(1,:,:) = rest_txf; % format this as 'pnts' x chans x freq
    %         base_txf_mean(1,:,:) = squeeze(mean(abs(rest_txf(10+1:end-10,:,:)),1));
            %- clear rest_txf
            rest_txf = double.empty;
            par_save(base_txf_mean,fpath,sprintf('rest_avg_txf.mat'))
            fprintf('done.\n\n');
            %## (LOAD GAIT TIMEWARP) ======================================= %%
            fprintf('Loading Gait Time-Frequency Data...\n');
            fpath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
    %         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[TRIAL_TYPES{:}]);
            icatimef_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
            %- load .icatimef load-in parameters
            tmp = load(icatimef_f,'-mat');
    %         parameters = tmp.parameters;
            %## (APPLY SPCA) =============================================== %%
            fprintf('Running SPCA on All Conditions\n');
            %- sPCA Algorithm
            [ERSP,GPM,gait_avg,output_struct] = apply_spca_cond_timewarp(tmp,base_txf_mean);
            [ERSP_corr,GPM_corr,PSC1,~,COEFFS] = specPCAdenoising(ERSP);
            
            %## SAVE PCA INFORMATION
            gait_ersp_struct = [];
            gait_ersp_struct.ID         = EEG.subject;
            gait_ersp_struct.Noise_cov  = [];% noise cov for kernel computation
            gait_ersp_struct.F_Rest     = output_struct.baseline_ersp;
            gait_ersp_struct.TF         = gait_avg;
            gait_ersp_struct.ERSP_uncor = ERSP;
            gait_ersp_struct.GPM_uncor  = GPM;
            gait_ersp_struct.ERSP       = ERSP_corr;
            gait_ersp_struct.GPM        = GPM_corr;
            gait_ersp_struct.PSC1       = PSC1;
            gait_ersp_struct.chanlocs   = EEG.chanlocs;
            gait_ersp_struct.icatimefopts = tmp_out;
            gait_ersp_struct.warptimes  = averaged_warpto_events;
            gait_ersp_struct.ntimes     = TIMEWARP_NTIMES;
            par_save(gait_ersp_struct,fpath,'gait_ersp_spca.mat');
            %## PLOT
            fig = figure(); set(gcf, 'position', [0 0 600 500]);
    %         plot(tmp.freqs, squeeze(output_struct.baseline_ersp)', 'k-');
            plot(tmp.freqs, squeeze(base_txf_mean)', 'k-');
            ylabel('Amplitude (\muV)');
            xlabel('Frequency (Hz)');
            grid on; box off
            title('Baseline ERSP (rest)');
            exportgraphics(fig,[fpath filesep 'allcond_baseline_avgs.jpg']);

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
                [ERSP_corr,GPM_corr,gait_avg,output_struct] = apply_spca_cond_timewarp(tmp,base_txf_mean,...
                    'COEFFS',COEFFS,...
                    'COND_STR',COND_STR);
                struct_out = struct('ersp_corr',ERSP_corr,...
                    'gpm_corr',GPM_corr,...
                    'times',tmp.times,...
                    'freqs',tmp.freqs,...
                    'pc1',gait_avg,...
                    'coeffs',COEFFS,...
                    'apply_spca_cond',output_struct);
                par_save(struct_out,fpath,sprintf('cond%s_spca_ersp.mat',conds{cond_i}));

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
            fpath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
            icatimef_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
            delete(icatimef_f);
            fpath = [epoched_fPath filesep 'rest'];
            icatimef_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
            delete(icatimef_f);
            fprintf('Subject %s done. time: %0.2f',EEG.subject, toc(tt));
        catch e
            fprintf('\nError occured on subject %s\n%s\n',subject_chars{subj_i},getReport(e));
        end
    end
end
