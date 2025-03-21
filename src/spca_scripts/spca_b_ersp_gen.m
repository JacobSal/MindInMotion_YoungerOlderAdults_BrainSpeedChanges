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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');

%% (PARAMETERS) ======================================================== %%
%- compute measures for spectrum and ersp
FORCE_RECALC_ERSP = true;
DO_RECREATE_SPCA = false; %false;
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
    'do_recalc_epoch',true);% opengl('dsave', 'software') % might be needed to plot dipole plots?
%- compute measures for spectrum and ersp
FORCE_RECALC_PSD = true;
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
    'condition_gait',{{'0p25','0p5','0p75','1p0','flat','high','low','med'}});
%(03/05/2025) JS, for speed analysis removing the terrain conditions for
%better calcs.
%% (PATHS) ========================================================== %%
%- datset name
DATA_SET = 'MIM_dataset';
%- Study Name
STUDY_DNAME = '02202025_mim_yaoa_spca_calcs';
% STUDY_DNAME = 'dummy_study';
%- Subject Directory information
ICA_DIR_FNAME = '02212025_YAOAN117_iccR0p65_iccREMG0p4_chanrej_samprej';
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
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[save_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
else
    tmp = load('-mat',[save_dir filesep sprintf('%s.study',STUDY_FNAME_REST)]);
    STUDY_REST = tmp.STUDY;
end
%% CALCULATE GRANDAVERAGE WARPTOs
TIMEWARP_NTIMES = floor(500/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
averaged_warpto_events = STUDY_GAIT.etc.a_epoch_process.avg_warpto_events;
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
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
skip_b1 = {'H1010','H1012','H1017','H1019','H1007','H1011','H1009','H1004','H1013','H1018', ...
    'H1035','H2042','H1031','H1026','H1047','H1039','H2023','H1034','H2018_FU','H2012_FU'};
skip_b2 = {'H2034','H1038','H1033','H1024','H1030','H1045','H2017','H2039','H2022','H2008', ...
'H2027','H1022','H1029','H1032','H1037','H1042','H2015','H2038','H2021','H2002','H2026', ...
'H1020','H1041','H1027','H1036','H3103','H2013','H2037','H1048','H2020','H3120','NH3030','H2025'};
skip_subjs = [skip_b1,skip_b2];
subj_chars = {STUDY_GAIT.datasetinfo.subject};
% for subj_i = 1:length(subj_chars)
parfor subj_i = 1:length(STUDY_GAIT.datasetinfo)
    %%
    %## PARAM COPIES
    tmp_study_gait = STUDY_GAIT;
    tmp_study_rest = STUDY_REST;
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    subj_char = tmp_study_gait.datasetinfo(subj_i).subject;

    %## DOUBLE CHECK SAME SUBJECT FROM STUDY
    subj_ii = find(strcmp(subj_char,{tmp_study_rest.datasetinfo.subject}));    

    if ~any(strcmp(subj_char,skip_subjs))
        %## CHECK DATA EXISTENCE
        %- gait
        gait_fpath = tmp_study_gait.datasetinfo(subj_i).filepath;
        if ~exist(gait_fpath,'dir')
            mkdir(gait_fpath)
        end
        %- rest
        rest_fpath = tmp_study_rest.datasetinfo(subj_ii).filepath;
        if ~exist(rest_fpath,'dir')
            mkdir(rest_fpath)
        end
        tt = tic;
        try
            %## GAIT ERSP
            EEG = pop_loadset('filepath',tmp_study_gait.datasetinfo(subj_i).filepath, ...
                'filename',tmp_study_gait.datasetinfo(subj_i).filename);
            tmpf = strsplit(tmp_study_gait.datasetinfo(subj_i).filename,'.');
            tmpf{2} = 'fdt';
            icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',subj_char)];
            %-- calculate
            if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(subj_char,FINISHED_ADULTS))
                tmp_dat = strjoin(tmpf,'.');
                fprintf('Running Subject %s\n',subj_char);
                EEG = eeg_checkset(EEG,'loaddata');
                if isempty(EEG.icaact)
                    fprintf('%s) Recalculating ICA activations\n',subj_char);
                    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
                end
                tmps = tmp_study_gait;
                %-- overrride datasetinfo to trick std_precomp to run.
                tmps.datasetinfo = tmps.datasetinfo(subj_i);
                tmps.datasetinfo(1).index = 1;
                %-- determine timewarping parameters
                timewarp_param = EEG.timewarp.latencies;
                timewarpms_param = averaged_warpto_events;
                %-- no baseline correction
                [~, ~] = std_precomp(tmps,EEG,'components', ...
                    'savetrials','on', ...
                    'recompute','on', ...
                    'ersp','on', ...
                    'itc','off', ...
                    'erspparams',{'parallel','on', ...
                        'cycles',ERSP_PARAMS.cycles, ...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))), ...
                        'ntimesout',TIMEWARP_NTIMES, ...
                        'baseline',nan(), ...
                        'timewarp',timewarp_param, ...
                        'timewarpms',timewarpms_param}); %ERSP
                fprintf('Done calculating timewarped ERSPs for %s\n',subj_char);
                % EEG.data = tmp_dat;
            else
                fprintf('Timewarped ERSPs already calculated for %s\n',subj_char);
            end
        catch e
            fprintf(['error. identifier: %s\n',...
                    'error. %s\n',...
                    'error. on subject %s\n',...
                    'stack. %s\n'],e.identifier,e.message,tmp_study_gait.datasetinfo(subj_i).subject,getReport(e));
        end
        EEG.data = [];
        
        %## RESTING ERSP
        try
            EEG = pop_loadset('filepath',tmp_study_rest.datasetinfo(subj_ii).filepath, ...
                'filename',tmp_study_rest.datasetinfo(subj_ii).filename);
            tmpf = strsplit(tmp_study_rest.datasetinfo(subj_ii).filename,'.');
            tmpf{2} = 'fdt';
            icatimef_f = [EEG.filepath filesep sprintf('%s.icatimef',subj_char)];
            %-- calculate
            if ~exist(icatimef_f,'file') || FORCE_RECALC_ERSP % || any(strcmp(subj_char,FINISHED_ADULTS))
                tmp_dat = strjoin(tmpf,'.'); 
                fprintf('Running Subject %s\n',subj_char);
                EEG = eeg_checkset(EEG,'loaddata');
                if isempty(EEG.icaact)
                    fprintf('%s) Recalculating ICA activations\n',subj_char);
                    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
                end
                tmps = tmp_study_rest;
                %- overrride datasetinfo to trick std_precomp to run.
                tmps.datasetinfo = tmps.datasetinfo(subj_ii);
                tmps.datasetinfo(1).index = 1;
                %-- no baseline correction
                [~, ~] = std_precomp(tmps,EEG,'components', ...
                    'savetrials','on',...
                    'recompute','on', ...
                    'ersp','on', ...
                    'itc','off',...
                    'erspparams',{'parallel','on', ...
                        'cycles',ERSP_PARAMS.cycles,...
                        'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))), ...
                        'ntimesout',TIMEWARP_NTIMES}); %ERSP
                fprintf('Done calculating timewarped ERSPs for %s\n',subj_char);
                % EEG.data = tmp_dat;
            else
                fprintf('Timewarped ERSPs already calculated for %s\n',subj_char);
            end
        catch e
            fprintf(['error. identifier: %s\n',...
                    'error. %s\n',...
                    'error. on subject %s\n',...
                    'stack. %s\n'],e.identifier,e.message,tmp_study_rest.datasetinfo(subj_ii).subject,getReport(e));
        end
        EEG.data = [];
    
        %## (SPCA PROCEDURE) =========================================== %%
        try
            %## LOAD REST DATA
            icatimef_f = [rest_fpath filesep sprintf('%s.icatimef',subj_char)];
            %-- load .icatimef load-in parameters
            fprintf('Loading Resting Time-Frequency Data...\n');
            MEAN_STRUCT = struct('do_return',true, ...
                'fcn',@(x) squeeze(mean(mean(abs(x),3),2)));
            [ersp_nolog_tf,etc_inf] = spca_get_ersp_dat(icatimef_f,[], ...
                'MEAN_STRUCT',MEAN_STRUCT);
            ersp_nolog_tf = squeeze(ersp_nolog_tf);
            ersp_nolog_tf = permute(ersp_nolog_tf,[2,1]);
            base_txf_mean = zeros(1,size(ersp_nolog_tf,1),size(ersp_nolog_tf,2));
            base_txf_mean(1,:,:) = ersp_nolog_tf; % format this as 'pnts' x chans x freq
            %-- clear rest_txf
            ersp_nolog_tf = double.empty;
            par_save(base_txf_mean,rest_fpath,sprintf('rest_avg_txf.mat'))
            fprintf('done.\n\n');

            %## (LOAD GAIT TIMEWARP) ======================================= %%
            fprintf('Loading Gait Time-Frequency Data...\n');
            icatimef_f = [gait_fpath filesep sprintf('%s.icatimef',subj_char)];
            MEAN_STRUCT = struct('do_return',true, ...
                'fcn',@(x) squeeze(mean(abs(x),3)));
            [ersp_nolog_tf,etc_inf] = spca_get_ersp_dat(icatimef_f,[], ...
                'MEAN_STRUCT',MEAN_STRUCT);
            % ersp_nolog_tf = squeeze(ersp_nolog_tf); 
            ersp_nolog_tf = permute(ersp_nolog_tf,[2,3,1]); % pnts x chans x freqs
            %--
            freqs = etc_inf.freqs;
        
            %## (APPLY SPCA) =============================================== %%
            fprintf('Running SPCA on All Conditions\n');
            %-- baseline to resting
            ERSP = 20*bsxfun(@minus,log10(ersp_nolog_tf), log10(base_txf_mean)); % decibel calculation
            %-- further baseline to mean across time
            GPM = bsxfun(@minus,ERSP,mean(ERSP));
            %-- spca alg
            [ERSP_corr,GPM_corr,PSC1,~,COEFFS] = spca_denoising_alg(ERSP);
        
            %## SAVE PCA INFORMATION
            gait_ersp_struct = struct('ID',subj_char, ...
                'noise_cov',COEFFS, ...
                'f_rest',base_txf_mean, ...
                'tf',[], ...
                'ERSP_uncor',ERSP, ...
                'GPM_uncor',GPM, ...
                'ERSP',ERSP_corr, ...
                'GPM',GPM_corr, ...
                'PSC1',PSC1, ...
                'chanlocs',EEG.chanlocs, ...
                'icatimefopts',etc_inf, ...
                'warptimes',{averaged_warpto_events}, ...
                'ntimes',TIMEWARP_NTIMES);
            par_save(gait_ersp_struct,gait_fpath,'gait_ersp_spca.mat');
        
            %## PLOT
            fig = figure();
            set(gcf, 'position', [0 0 600 500]);
            plot(freqs, squeeze(20*log10(base_txf_mean))', 'k-');
            ylabel('Amplitude (\muV)');
            xlabel('Frequency (Hz)');
            grid on; box off
            title('Baseline ERSP (rest)');
            exportgraphics(fig,[gait_fpath filesep 'allcond_baseline_avgs.jpg']);
        
            %## DENOISE BY CONDITION ======================================= %%
            fprintf('\nRunning Denoising with sPCA\n');
            conds = unique({etc_inf.trialinfo.cond});
            for cond_i = 1:length(conds)
                cond_char = conds{cond_i};
                %--
                icatimef_f = [gait_fpath filesep sprintf('%s.icatimef',subj_char)];
                MEAN_STRUCT = struct('do_return',true, ...
                    'fcn',@(x) squeeze(mean(abs(x),3)));
                [ersp_nolog_tf,etc_inf] = spca_get_ersp_dat(icatimef_f,cond_char, ...
                    'MEAN_STRUCT',MEAN_STRUCT);
                ersp_nolog_tf = squeeze(ersp_nolog_tf); 
                ersp_nolog_tf = permute(ersp_nolog_tf,[2,3,1]); % pnts x chans x freqs
                %-- baseline to resting
                ERSP = 20*bsxfun(@minus,log10(ersp_nolog_tf), log10(base_txf_mean)); % decibel calculation
                %-- further baseline to mean across time
                GPM = bsxfun(@minus,ERSP,mean(ERSP));            
                %-- spca alg
                [ERSP_corr,GPM_corr,PSC1,~,~] = spca_denoising_alg(ERSP,COEFFS);
                struct_out = struct('ersp_corr',ERSP_corr,...
                    'gpm_corr',GPM_corr,...
                    'ersp_uncorr',ERSP,...
                    'gpm_uncorr',GPM,...
                    'times',etc_inf.times,...
                    'freqs',etc_inf.freqs,...
                    'pc1',PSC1,...
                    'coeffs',COEFFS,...
                    'apply_spca_cond',etc_inf);
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
            fprintf('%s) Deleting gait data...\n',subj_char)
            icatimef_f = [gait_fpath filesep sprintf('%s.icatimef',subj_char)];
            delete(icatimef_f);
            % fpath = [epoched_fPath filesep 'rest'];
            fprintf('%s) Deleting rest data...\n',subj_char)
            icatimef_f = [rest_fpath filesep sprintf('%s.icatimef',subj_char)];
            delete(icatimef_f);
            fprintf('Subject %s done. time: %0.2f',subj_char, toc(tt));
        catch e
            fprintf(['error. identifier: %s\n',...
                    'error. %s\n',...
                    'error. on subject %s\n',...
                    'stack. %s\n'],e.identifier,e.message,subj_char,getReport(e));
        end
    else
        fprintf('%s) Subject already complete.\n',subj_char);
    end
end
