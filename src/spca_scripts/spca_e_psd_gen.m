%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/spca_scripts/run_spca_e_psd_gen.sh

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
ADD_ALL_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
        SRC_DIR = STUDY_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
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
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','high','low','med'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%- compute measures for spectrum and ersp
FORCE_RECALC_PSD = false;
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
SPCA_PARAMS = struct('analysis_type','component',...
    'event_char','RHS',...
    'epoch_min_max',[1,4.25],...
    'n_resamples',100,...
    'timewarp_events',{{'RHS','LHS','LTO','RTO'}},...
    'condition_base','rest',...
    'condition_gait',{{'0p25','0p5','0p75','1p0','flat','high','low','med'}});
%(03/05/2025) JS, for speed analysis removing the terrain conditions for
%better calcs.
%% (PATHS)
%- datset name & 
DATA_SET = 'MIM_dataset';
ICA_DIR_FNAME = '02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej';
STUDY_DNAME = '02202025_mim_yaoa_spca_calcs';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study_all';
STUDY_FNAME_REST = 'spca_rest_slide_study_all';
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
ica_data_dir = [studies_dir filesep ICA_DIR_FNAME]; % JACOB,SAL(02/23/2023)
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
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
    'freqrange',SPEC_PARAMS.plot_freqrange, ...
    'plotmode','condensed',...
    'plotconditions','together', ...
    'ylim',SPEC_PARAMS.plot_ylim, ...
    'plotgroups','together');
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
    'freqrange',SPEC_PARAMS.plot_freqrange, ...
    'plotmode','condensed',...
    'plotconditions','together', ...
    'ylim',SPEC_PARAMS.plot_ylim, ...
    'plotgroups','together');
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
subj_chars = {STUDY_REST.datasetinfo.subject};
parfor subj_i = 1:length(subj_chars)
    %-
    %%
    fprintf('Running Subject %s\n',subj_chars{subj_i}); 
    tt = tic;    
    %## PARAMETERS
    base_txf_mean = [];
    tmp_study_rest = STUDY_REST;
    tmp_study = STUDY;
    tmp_spec_params = SPEC_PARAMS;
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    
    %## DOUBLE CHECK SAME SUBJECT FROM STUDY
    subj_ii = find(strcmp(subj_chars{subj_i},{tmp_study.datasetinfo.subject}));
    
    %## CALCULATE GAIT PSD
    %- load EEG
    EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj_ii).filepath,...
        'filename',tmp_study.datasetinfo(subj_ii).filename);       
    %-
    icaspec_f = [EEG.filepath filesep sprintf('%s.icaspec',EEG.subject)];
    if ~exist(icaspec_f,'file') || FORCE_RECALC_PSD
        tmp_study_calc = tmp_study;
        %- load data
        tmp_fdt_fpath = EEG.data;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %- overrride datasetinfo to trick std_precomp to run.
        tmp_study_calc.datasetinfo = tmp_study_calc.datasetinfo(subj_ii);
        tmp_study_calc.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp_study_calc, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',tmp_spec_params.specmode,'freqfac',tmp_spec_params.freqfac,...
                    'freqrange',tmp_spec_params.freqrange,'logtrials',tmp_spec_params.logtrials});
        fprintf('Done calculating PSDs for %s\n',EEG.subject);
        %- clear data
        EEG.data = tmp_fdt_fpath;
    else
        fprintf('PSDs already calculated for %s\n',EEG.subject);
    end    

    %## CALCULATE REST PSD
    %- load EEG
    EEG = pop_loadset('filepath',tmp_study_rest.datasetinfo(subj_i).filepath,...
        'filename',tmp_study_rest.datasetinfo(subj_i).filename);    
    %- 
    icaspec_f = [EEG.filepath filesep sprintf('%s.icaspec',EEG.subject)];
    if ~exist(icaspec_f,'file') || FORCE_RECALC_PSD
        tmp_study_calc = tmp_study_rest;
        %- load data
        tmp_fdt_fpath = EEG.data;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        % EEG = eeglab_makecontinuous(EEG);
        %(01/14/2025) JS, is this really necessary?
        %- overrride datasetinfo to trick std_precomp to run.
        tmp_study_calc.datasetinfo = tmp_study_calc.datasetinfo(subj_i);
        tmp_study_calc.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp_study_calc, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',tmp_spec_params.specmode,'freqfac',tmp_spec_params.freqfac,...
                    'freqrange',tmp_spec_params.freqrange,'logtrials',tmp_spec_params.logtrials});
        fprintf('Done calculating PSDs for %s\n',EEG.subject);
        %- clear data
        EEG.data = tmp_fdt_fpath;
    else
        fprintf('PSDs already calculated for %s\n',EEG.subject);
    end    

    %## SPCA PROCESS
    try
        %## (LOAD REST) ================================================ %%
        epoched_fPath = strsplit(EEG.filepath,filesep);
        epoched_fPath = strjoin(epoched_fPath(1:end-1),filesep);
        fpath = [epoched_fPath filesep 'rest'];
        icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        %-- load
        [psd_log_f,etc_inf] = spca_get_psd_dat(icaspec_f,[]); % freq x epoch x chan
        fprintf('Number of epochs: %i\n',size(psd_log_f,2))
        rest_psdf = permute(psd_log_f,[2,1,3]);
        %- average over time, keep magnitude (not power -> would amplify outliers)
        base_txf_mean = mean(rest_psdf);
        base_txf_mean = permute(base_txf_mean,[1,3,2]); % 1 x chan x freq
        %- clear rest_txf
        par_save(rest_psdf,fpath,sprintf('rest_avg_psdf.mat'));
        rest_psdf = double.empty;

        %## (LOAD GAIT TIMEWARP) ======================================= %%
        fprintf('Loading Gait Frequency Data...\n');
        fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
        icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        %-- load .icaspec load-in parameters
        [psd_log_f,etc_inf] = spca_get_psd_dat(icaspec_f,[]); % freq x epoch x chan
        %-- average across trials & baseline to rest
        psd_log_f = permute(psd_log_f,[2,1,3]); % epoch x freq x chan
        psd_avg = mean(psd_log_f); 
        psd_avg = permute(psd_avg,[1,3,2]);
        psd_baselined = bsxfun(@minus,psd_avg,base_txf_mean);
        %-- plot
        fprintf('Running SPCA on All Conditions\n');
        psd_in = squeeze(psd_baselined(:,:,:));
        fig = figure();
        set(gcf, 'position', [0 0 600 500]);
        ax = axes();
        plot(etc_inf.freqs,psd_in,'-k')
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on;
        box off;
        exportgraphics(fig,[fpath filesep 'allcond_avg_psds_uncorr.jpg']);

        %## (APPLY SPCA) =============================================== %%

        %-- apply spca
        [psd_corr_based,based_psd_corr,psd_corr_psc1,psd_psc,coeffs] = specPCAdenoising(psd_baselined);
        
        %## PLOT 
        psd_in = squeeze(psd_corr_based(:,:,:));
        fig = figure();
        set(gcf, 'position', [0 0 600 500]);
        ax = axes();
        plot(etc_inf.freqs,psd_in,'-k')
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on;
        box off;
        exportgraphics(fig,[fpath filesep 'allcond_avg_psds_corr.jpg']);

        %## SAVE PCA INFORMATION
        psd_struct = struct('ID', {EEG.subject}, ...
            'noise_cov', [], ... % noise covariance for kernel computation
            'psd_orig_avg', psd_avg, ...
            'psd_orig_based', psd_baselined, ...
            'psd_corr', psd_corr_based, ...
            'based_psd_corr', based_psd_corr, ...
            'psdc1', psd_corr_psc1, ...
            'chanlocs', EEG.chanlocs, ...
            'etc_inf', etc_inf ...
        );
        par_save(psd_struct,fpath,'gait_psd_spca.mat');

        %## PLOT
        fig = figure(); set(gcf, 'position', [0 0 600 500]);
        plot(etc_inf.freqs, squeeze(base_txf_mean)', 'k-');
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on; box off
        title('Baseline ERSP (rest)');
        exportgraphics(fig,[fpath filesep 'allcond_baseline_avgs_psds.jpg']);

        %% DENOISE BY CONDITION ======================================== %%
        fprintf('\nRunning Denoising with sPCA\n');
        conds = unique({etc_inf.trialinfo.cond});
        fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
        icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        psd_store = zeros(length(conds),size(psd_corr_based,2),size(psd_corr_based,3));
        for cond_i = 1:length(conds)
            %## (APPLY SPCA) =============================================== %%
            fprintf('Running SPCA on %s\n',conds{cond_i});
            [psd_log_f,etc_inf] = spca_get_psd_dat(icaspec_f,conds{cond_i});
            %-- average across trials
            psd_avg = permute(psd_log_f,[2,1,3]); % epoch x freq x chan
            psd_avg = mean(psd_avg);
            psd_avg = permute(psd_avg,[1,3,2]); % epoch x chan x freq
            psd_baselined = bsxfun(@minus,psd_avg,base_txf_mean);
            %-- apply spca
            [psd_corr_based,based_psd_corr,psd_corr_psc1,psd_psc,~] = specPCAdenoising(psd_baselined,coeffs);
            struct_out = struct('condition',conds{cond_i},...
                'psd_corr_based',psd_corr_based,...
                'psd_orig_avg',psd_avg,...
                'baseline_psd',base_txf_mean,...
                'psd_corr_psc1',psd_corr_psc1,...
                'psd_orig_baselined',psd_baselined,...
                'freqs',etc_inf.freqs,...
                'coeffs',coeffs);
            par_save(struct_out,fpath,sprintf('cond%s_spca_psd.mat',conds{cond_i}));
            %-- plotting var
            psd_store(cond_i,:,:) = psd_corr_based;
        end

        %## VALIDATION PLOT COND-BY-COND
        AXES_DEFAULT_PROPS = {'box','off', ...
            'xtick',[], ...
            'ytick',[],...
            'ztick',[], ...
            'xcolor',[1,1,1], ...
            'ycolor',[1,1,1]};
        [~,co_i] = max(EEG.etc.ic_classification.ICLabel.classifications(:,1));
        co_i = co_i - sum(EEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.75);
        % co_i = 10;
        cmaps = linspecer(length(conds));
        fig = figure();
        set(gcf, 'position', [0 0 920 920]);
        set(gca,AXES_DEFAULT_PROPS{:})
        hold on;
        ax = axes();
        lis = [];
        LINE_STRUCT = struct('line_alpha',0.75, ...
            'err_alpha',0.3, ...
            'do_err_shading',false)
        for cond_i = 1:length(conds)
            psd_in = squeeze(psd_store(cond_i,co_i,:));
            LINE_STRUCT.line_color = cmaps(cond_i,:);
            LINE_STRUCT.err_color = cmaps(cond_i,:);
            %--            
            [~,~,li] = plot_psd(ax,psd_in,etc_inf.freqs, ...
                'LINE_STRUCT',LINE_STRUCT);
            li.DisplayName = conds{cond_i};
            lis = [lis,li];
        end
        legend(lis)
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on;
        exportgraphics(fig,[fpath,sprintf('allcond_comp%i_spca_psd.jpg',co_i)]);

        %## STRIDE-BY-STRIDE IMPLEMENT
        fprintf('Running stride-by-stride SPCA\n');
        [psd_log_f,etc_inf] = spca_get_psd_dat(icaspec_f,[]);
        psd_in = permute(psd_log_f,[2,3,1]);
        psd_store = zeros(size(psd_in));
        % psd_store_b = zeros(size(psd_in));
        tt = tic();
        for t_i = 1:size(psd_in,1)
            if mod(t_i,10) == 0
                fprintf('.');
            end
            psd_baselined = bsxfun(@minus,psd_in(t_i,:,:),base_txf_mean);
            %-- apply spca
            [psd_corr_based,~,~,~,~] = specPCAdenoising(psd_baselined,coeffs);
            psd_store(t_i,:,:) = psd_corr_based;
            % psd_store_b(t_i,:,:) = psd_baselined;
        end
        fprintf('done. %0.2f min',toc(tt)/(60))
        struct_out = struct('psd_corr_based',psd_store,...
            'baseline_psd',base_txf_mean,...
            'psd_corr_psc1',[],...
            'etc_inf',etc_inf,...
            'coeffs',coeffs);
        par_save(struct_out,fpath,sprintf('allcond_sbs_spca_psd.mat'));        
        
        %## VALIDATION PLOT STEP-BY-STEP
        % co_i = 1;
        cmaps = linspecer(length(conds));
        fig = figure();
        set(gcf, 'position', [0 0 920 920]);
        set(gca,AXES_DEFAULT_PROPS{:})
        hold on;
        ax = axes();
        lis = [];
        LINE_STRUCT = struct.empty;
        for cond_i = 1:length(conds)
            t_ii = strcmp(conds{cond_i},{etc_inf.trialinfo.cond});
            psd_in = squeeze(psd_store(t_ii,co_i,:))';
            LINE_STRUCT.line_color = cmaps(cond_i,:);
            LINE_STRUCT.err_color = cmaps(cond_i,:);
            LINE_STRUCT.line_alpha = 0.75;
            LINE_STRUCT.err_alpha = 0.3;
            %--            
            [~,~,li] = plot_psd(ax,psd_in,etc_inf.freqs, ...
                'LINE_STRUCT',LINE_STRUCT);
            li.DisplayName = conds{cond_i};
            lis = [lis,li];
        end
        legend(lis)
        ylabel('Amplitude (\muV)');
        xlabel('Frequency (Hz)');
        grid on;
        exportgraphics(fig,[fpath,sprintf('allcond_comp%i_sbs_spca_psd.jpg',co_i)]);
        
        %## CLEANUP AND FREE UP DRIVE SPACE
        % fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
        % icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        % delete(icaspec_f);
        % fpath = [epoched_fPath filesep 'rest'];
        % icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        % delete(icaspec_f);
        %(01/14/2025) JS, unsure if necessary given small memory impact of
        %icaspec files (~500 MB)
        fprintf('Subject %s done. time: %0.2f',EEG.subject, toc(tt));
    catch e
        fprintf('\nError occured on subject %s\n%s\n',subj_chars{subj_i},getReport(e));
    end
end