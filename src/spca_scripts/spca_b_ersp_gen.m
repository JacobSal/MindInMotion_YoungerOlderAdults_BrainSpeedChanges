%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/spca_scripts/run_spca_b_ersp_gen.sh

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
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SRC_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    SRC_DIR = fileparts(SCRIPT_DIR);    
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- Study Name
STUDY_DIR_FNAME = '03232024_spca_analysis_OA';
%- Subject Directory information
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_FNAME_REST = 'rest_epoch_study_yaoa';
STUDY_FNAME_GAIT = 'gait_epoch_study_yaoa';
%-
studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
ica_data_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies' filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
study_fpath = [studies_dir filesep sprintf('%s',STUDY_DIR_FNAME)];
%- create new study directory
if ~exist(study_fpath,'dir')
    mkdir(study_fpath);
end
%- study group and saving
% TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
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
    'slide_cond_chars',{{'rest'}},...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'do_recalc_epoch',true);
%- compute measures for spectrum and ersp
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false;
DO_RECREATE_SPCA = true; %false;
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
%% ===================================================================== %%
% if ~ispc
%     [STUDY_GAIT,ALLEEG_GAIT] = pop_loadstudy('filename',[STUDY_FNAME_GAIT '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY_GAIT,ALLEEG_GAIT] = pop_loadstudy('filename',[STUDY_FNAME_GAIT '.study'],'filepath',save_dir);
% end
% if ~ispc
%     [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[STUDY_FNAME_REST '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[STUDY_FNAME_REST '.study'],'filepath',save_dir);
% end
if ~ispc
    tmp = load('-mat',[study_fpath filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
else
    tmp = load('-mat',[study_fpath filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
end
if ~ispc
    tmp = load('-mat',[study_fpath filesep sprintf('%s_UNIX.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
else
    tmp = load('-mat',[study_fpath filesep sprintf('%s.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
end
subj_chars = {STUDY_GAIT.datasetinfo.subject};
%% CALCULATE GRANDAVERAGE WARPTOs
TIMEWARP_NTIMES = floor(EEG.srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
EEG = load([STUDY_GAIT.datasetinfo(1).filepath filesep STUDY_GAIT.datasetinfo(1).filename],'-mat');
warptos = nan(length(STUDY_GAIT.datasetinfo),length(EEG.etc.timewarp_by_cond(1).warpto));
%- subject based timings
parfor subj_i = 1:length(STUDY_GAIT.datasetinfo)    
    tmp_study = STUDY_GAIT;
    fprintf('Loading subject %s\n',tmp_study.datasetinfo(subj_i).subject)
    %-
    % EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj).filepath, ...
    %     'filename',tmp_study.datasetinfo(subj).filename);
    %-
    EEG = load([tmp_study.datasetinfo(subj_i).filepath filesep tmp_study.datasetinfo(subj_i).filename],'-mat');
    %- assign percondition timewarping
    warptos(subj_i,:) = nanmedian(cat(1,EEG.etc.timewarp_by_cond.warpto));
end
averaged_warpto_events = floor(nanmean(warptos)); % tends to be longer? (e.g., [0,262,706,982,1415])
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
STUDY_GAIT.etc.averaged_warpto_events = averaged_warpto_events;
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);

%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION 
%## ersp plot per cluster per condition
STUDY_GAIT = pop_statparams(STUDY_GAIT,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY_GAIT = pop_erspparams(STUDY_GAIT,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_CROP_TIMES);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_GAIT = pop_specparams(STUDY_GAIT,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
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
subject_chars = {STUDY_GAIT.datasetinfo.subject};
for subj_i = 1:length(subject_chars)
    %## PARAM COPIES
    tmp_study_gait = STUDY_GAIT;
    tmp_study_rest = STUDY_REST;
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    %## CHECK DATA EXISTENCE
    % epoched_fPath = [save_dir filesep subj_chars{subj_i}];
    % epoched_fPath = [study_fpath filesep subj_chars{subj_i} filesep 'step_to_step'];
    %- gait
    gait_fpath = tmp_study_gait.datasetinfo(subj_i).filepath; %[epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
    % gait_fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subj_chars{subj_i},[tmp_epoch_params.gait_trial_chars{:}]);
    if ~exist(gait_fpath,'dir')
        mkdir(gait_fpath)
    end
    %- rest
    rest_fpath = tmp_study_rest.datasetinfo(subj_i).filepath; %[epoched_fPath filesep tmp_epoch_params.slide_cond_chars{:}];
    % rest_fname = sprintf('%s_%s_EPOCH_TMPEEG.set',subj_chars{subj_i},'rest');
    if ~exist(rest_fpath,'dir')
        mkdir(rest_fpath)
    end
    %##
    chk = ~all(cellfun(@(x) exist([gait_fpath filesep sprintf('cond%s_spca_ersp.mat',x)],'file'),tmp_epoch_params.gait_trial_chars));
    if DO_RECREATE_SPCA || chk
        tt = tic;
        %## GAIT ERSP
        EEG = pop_loadset('filepath',tmp_study_gait.datasetinfo(subj_i).filepath, ...
            'filename',tmp_study_gait.datasetinfo(subj_i).filename);
        tmp_dat = EEG.data;
        fprintf('Running Subject %s\n',EEG.subject);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
        if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            %- overrride datasetinfo to trick std_precomp to run.
            tmp_study_gait.datasetinfo = STUDY_GAIT.datasetinfo(subj_i);
            tmp_study_gait.datasetinfo(1).index = 1;
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
                [~, ~] = std_precomp(tmp_study_gait,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                        'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                        'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                        'trialbase','off','basenorm','on'}); %ERSP
            else
                % No baseline correction
                [~, ~] = std_precomp(tmp_study_gait,EEG,'components','savetrials','on',...
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
        EEG.data = tmp_dat;

        %## RESTING ERSP
        EEG = pop_loadset('filepath',tmp_study_rest.datasetinfo(subj_i).filepath, ...
            'filename',tmp_study_rest.datasetinfo(subj_i).filename);
        tmp_dat = EEG.data;
        fprintf('Running Subject %s\n',EEG.subject);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
        if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(EEG.subject,FINISHED_ADULTS))
            EEG = eeg_checkset(EEG,'loaddata');
            if isempty(EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
            end
            %- overrride datasetinfo to trick std_precomp to run.
            tmp_study_rest.datasetinfo = STUDY_REST.datasetinfo(subj_i);
            tmp_study_rest.datasetinfo(1).index = 1;
            %-
            if DO_BASELINE_CORRECTION
                % Baseline correction
                [~, ~] = std_precomp(tmp_study_rest,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                        'ntimesout',TIMEWARP_NTIMES,...
                        'trialbase','off','basenorm','on'}); %ERSP
            else
                % No baseline correction
                [~, ~] = std_precomp(tmp_study_rest,EEG,'components','savetrials','on',...
                        'recompute','on','ersp','on','itc','off',...
                        'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES}); %ERSP
            end
            fprintf('Done calculating timewarped ERSPs for %s\n',EEG.subject);
        else
            fprintf('Timewarped ERSPs already calculated for %s\n',EEG.subject);
        end
        EEG.data = tmp_dat;
        %## (SPCA PROCEDURE) =========================================== %%
        %- Rest data
        epoched_fPath = strsplit(EEG.filepath,filesep);
        epoched_fPath = strjoin(epoched_fPath(1:end-1),filesep);
        icatimef_f = [rest_fpath filesep sprintf('%s.icatimef',EEG.subject)];
        %- load .icatimef load-in parameters
        fprintf('Loading Resting Time-Frequency Data...\n');
        tf_mat_struct = load(icatimef_f,'-mat');
        [ersp_tf,freqs,trialinfo,ersp_nolog_tf] = spca_get_ersp_dat(tf_mat_struct);
        %- reshape data [pnts x chans x freq]
        % fn = fieldnames(tmp);
        % inds = find(contains(fieldnames(tmp),'comp'));
        % 
        % test = tmp.(fn{inds(1)});
        % rest_txf = zeros(size(test,1),size(test,2),size(test,3),length(inds),'double');
        % for i = 1:length(inds)
        %     rest_txf(:,:,:,i) = tmp.(fn{inds(i)}); % freq x time x epoch x chan
        % end
        %-- average over time, keep magnitude (not power -> would amplify outliers)
        ersp_nolog_tf = squeeze(mean(mean(abs(ersp_nolog_tf),3),2));
        ersp_nolog_tf = permute(ersp_nolog_tf,[2,1]);
        
        base_txf_mean = zeros(1,size(ersp_nolog_tf,1),size(ersp_nolog_tf,2));
        base_txf_mean(1,:,:) = ersp_nolog_tf; % format this as 'pnts' x chans x freq
        %- clear rest_txf
        ersp_nolog_tf = double.empty;
        par_save(base_txf_mean,rest_fpath,sprintf('rest_avg_txf.mat'))
        fprintf('done.\n\n');
        
        %## (LOAD GAIT TIMEWARP) ======================================= %%
        fprintf('Loading Gait Time-Frequency Data...\n');
        icatimef_f = [gait_fpath filesep sprintf('%s.icatimef',EEG.subject)];
        %- load .icatimef load-in parameters
        tf_mat_struct = load(icatimef_f,'-mat');

        %## (APPLY SPCA) =============================================== %%
        fprintf('Running SPCA on All Conditions\n');
        %- get the ERSP data from the mat file
        % [ERSP,GPM,gait_avg,output_struct] = apply_spca_cond_timewarp(tf_mat_struct,base_txf_mean);
        [ersp_tf,freqs,trialinfo,ersp_nolog_tf] = spca_get_ersp_dat(tf_mat_struct);
        %- average across trials
        gait_avg = squeeze(mean(abs(ersp_nolog_tf),3)); % average over trials and take abs() to ensure amplitude calc
        gait_avg = permute(gait_avg,[2,3,1]);
        %- baseline to resting
        ERSP = 20*bsxfun(@minus,log10(gait_avg), log10(base_txf_mean)); % decibel calculation
        %- further baseline to mean across time
        % GPM = bsxfun(@minus,ERSP,mean(ERSP));
        %- spca alg
        [ERSP,GPM,PSC1,~,COEFFS] = spca_denoising_alg(ERSP);

        %## SAVE PCA INFORMATION
        gait_ersp_struct = [];
        gait_ersp_struct.ID         = EEG.subject;
        gait_ersp_struct.Noise_cov  = [];% noise cov for kernel computation
        gait_ersp_struct.F_Rest     = output_struct.baseline_ersp;
        gait_ersp_struct.TF         = gait_avg;
        gait_ersp_struct.ERSP_uncor = ERSP;
        gait_ersp_struct.GPM_uncor  = GPM;
        gait_ersp_struct.ERSP       = ERSP;
        gait_ersp_struct.GPM        = GPM;
        gait_ersp_struct.PSC1       = PSC1;
        gait_ersp_struct.chanlocs   = EEG.chanlocs;
        gait_ersp_struct.icatimefopts = tmp_out;
        gait_ersp_struct.warptimes  = averaged_warpto_events;
        gait_ersp_struct.ntimes     = TIMEWARP_NTIMES;
        par_save(gait_ersp_struct,gait_fpath,'gait_ersp_spca.mat');

        %## PLOT
        fig = figure(); set(gcf, 'position', [0 0 600 500]);
        plot(tmp.freqs, squeeze(base_txf_mean)', 'k-');
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on; box off
        title('Baseline ERSP (rest)');
        exportgraphics(fig,[gait_fpath filesep 'allcond_baseline_avgs.jpg']);

        %## DENOISE BY CONDITION ======================================= %%
        fprintf('\nRunning Denoising with sPCA\n');
        conds = unique({trialinfo.cond});
        for cond_i = 1:length(conds)
            condition_char = conds{cond_i};
            %- get the ERSP data from the mat file
            [ersp_tf,freqs,~,ersp_nolog_tf] = spca_get_ersp_dat(tf_mat_struct,condition_char);
            %- baseline to resting
            ERSP = 20*bsxfun(@minus,log10(ersp_nolog_tf), log10(base_txf_mean)); % decibel calculation
            %- further baseline to mean across time
            % GPM = bsxfun(@minus,ERSP,mean(ERSP));
            %- spca alg
            [ERSP,GPM,PSC1,~,~] = spca_denoising_alg(ERSP,COEFFS);
            % [ERSP_corr,GPM_corr,gait_avg,output_struct] = apply_spca_cond_timewarp(tf_mat_struct,base_txf_mean,...
            %     'COEFFS',COEFFS,...
            %     'COND_STR',condition_char);
            struct_out = struct('ersp_corr',ERSP,...
                'gpm_corr',GPM,...
                'times',tf_mat_struct.times,...
                'freqs',tf_mat_struct.freqs,...
                'pc1',gait_avg,...
                'coeffs',COEFFS,...
                'apply_spca_cond',output_struct);
            par_save(struct_out,gait_fpath,sprintf('cond%s_spca_ersp.mat',conds{cond_i}));
            %## VALIDATION PLOTS
            %{
            fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) ERSP corrected',conds{cond_i}));
            if ~isempty(gait_fpath)
                exportgraphics(fig,[gait_fpath filesep sprintf('chan%i_cond%s_ersp_corr.jpg',CHAN_INT,conds{cond_i})]);
            end
            %-
            fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) GPM corrected',conds{cond_i}));
            if ~isempty(gait_fpath)
                exportgraphics(fig,[gait_fpath filesep sprintf('chan%i_cond%s_gpm_corr.jpg',CHAN_INT,conds{cond_i})]);
            end
            %-
            fig = plot_txf(squeeze(output_struct.erds_orig(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) ERSP original',conds{cond_i}));
            if ~isempty(gait_fpath)
                exportgraphics(fig,[gait_fpath filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,conds{cond_i})]);
            end
            %-
            fig = plot_txf(squeeze(output_struct.gpm_orig(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) GPM original',conds{cond_i}));
            if ~isempty(gait_fpath)
                exportgraphics(fig,[gait_fpath filesep sprintf('chan%i_cond%s_gpm_orig.jpg',CHAN_INT,conds{cond_i})]);
            end
            %##
            % fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],tmp.freqs);
            % title(sprintf('%s) PSC1 original',SPCA_PARAMS.condition_gait{cond_i}));
            % exportgraphics(fig,[gait_fpath filesep sprintf('chan%i_cond%s_psc1_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            %}
        end
        %## CLEANUP AND FREE UP DRIVE SPACE
        % fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
        icatimef_f = [gait_fpath filesep sprintf('%s.icatimef',EEG.subject)];
        delete(icatimef_f);
        % fpath = [epoched_fPath filesep 'rest'];
        icatimef_f = [rest_fpath filesep sprintf('%s.icatimef',EEG.subject)];
        delete(icatimef_f);
        fprintf('Subject %s done. time: %0.2f',EEG.subject, toc(tt));
    else
        fprintf('\nSubject %s already has sPCA ERSP data created...',subject_chars{subj_i})
    end
end
