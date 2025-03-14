%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/spca_scripts/run_spca_sbs_e_psd_gen.sh

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
    'condition_gait',{{'flat','low','med','high','0p25','0p5','0p75','1p0'}});
%% (PATHS)
%- datset name & 
DATA_SET = 'MIM_dataset';
ICA_DIR_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_DNAME = '01102025_mim_yaoa_spca_calcs';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study';
STUDY_FNAME_REST = 'spca_rest_slide_study';
%-
studies_dir = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
ica_data_dir = [studies_dir filesep ICA_DIR_FNAME];
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
subj_chars = {STUDY_REST.datasetinfo.subject};
parfor subj_i = 1:length(STUDY_REST.datasetinfo)
    %-
    fprintf('Running Subject %s\n',subj_chars{subj_i}); 
    tt = tic;    
    %## PARAMETERS
    base_txf_mean = [];
    tmp_study_rest = STUDY_REST;
    tmp_study = STUDY;
    tmp_spec_params = SPEC_PARAMS;
    tmp_epoch_params = DEF_EPOCH_PARAMS;
    
    %## CALCULATE GAIT PSD
    %- load EEG
    EEG = pop_loadset('filepath',tmp_study.datasetinfo(subj_i).filepath,...
        'filename',tmp_study.datasetinfo(subj_i).filename);       
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
        %- load .icaspec load-in parameters
        fprintf('Loading Resting Frequency Data...\n');
        tmp = load(icaspec_f,'-mat');
        [rest_psdf,freqs,trialinfo,psd_nolog_f] = spca_get_psd_dat(tmp)
        % %- reshape data [pnts x chans x freq]
        % fn = fieldnames(tmp);
        % inds = find(contains(fieldnames(tmp),'comp'));
        % inds_not = find(~contains(fieldnames(tmp),'comp'));
        % tmp_out = [];
        % for i = 1:length(inds_not)
        %     tmp_out.(fn{inds_not(i)}) = tmp.(fn{inds_not(i)});
        % end
        % test = tmp.(fn{inds(1)});
        % rest_psdf = zeros(size(test,1),size(test,2),length(inds),'double');
        % for i = 1:length(inds)
        %     rest_psdf(:,:,i) = tmp.(fn{inds(i)}); % freq x epoch x chan
        % end
        % rest_psdf = permute(rest_psdf,[2,1,3]);
        %- average over time, keep magnitude (not power -> would amplify outliers)
        base_txf_mean = mean(rest_psdf);
        base_txf_mean = permute(base_txf_mean,[1,3,2]);
        %- clear rest_txf
        rest_psdf = double.empty;
        par_save(tmp,fpath,sprintf('rest_avg_psdf.mat'));
        %## (LOAD GAIT TIMEWARP) ======================================= %%
        fprintf('Loading Gait Frequency Data...\n');
        fpath = [epoched_fPath filesep [tmp_epoch_params.gait_trial_chars{:}]];
        icaspec_f = [fpath filesep sprintf('%s.icaspec',EEG.subject)];
        %- load .icaspec load-in parameters
        tmp = load(icaspec_f,'-mat');
        %## (APPLY SPCA) =============================================== %%
        fprintf('Running SPCA on All Conditions\n');
        % [psd_baselined,psd_avg,output_struct] = apply_spca_cond_psd(tmp,base_txf_mean);
        [psd_f,freqs,trialinfo,psd_nolog_f] = spca_get_psd_dat(tmp)
        %- average across trials & baseline to rest
        psd_f = permute(psd_f,[2,1,3]);
        psd_avg = mean(psd_f); 
        psd_avg = permute(psd_avg,[1,3,2]);
        psd_baselined = bsxfun(@minus,psd_avg,base_txf_mean);
        %- apply spca
        % [psd_corr_based,based_psd_corr,psd_corr_psc1,psd_psc,coeffs] = specPCAdenoising(psd_baselined);
        [ersp_tf,freqs,trialinfo,ersp_nolog_tf] = spca_get_ersp_dat(psd_baselined);
        %## SAVE PCA INFORMATION
        gait_ersp_struct = [];
        gait_ersp_struct.ID             = EEG.subject;
        gait_ersp_struct.noise_cov      = []; % noise cov for kernel computation
        gait_ersp_struct.psd_orig_avg   = psd_avg;
        gait_ersp_struct.psd_orig_based = psd_baselined;
        gait_ersp_struct.psd_corr       = psd_corr_based;
        gait_ersp_struct.based_psd_corr = based_psd_corr;
        gait_ersp_struct.psdc1          = psd_corr_psc1;
        gait_ersp_struct.chanlocs       = EEG.chanlocs;
        gait_ersp_struct.icatimefopts   = tmp_out;
        par_save(gait_ersp_struct,fpath,'gait_psd_spca.mat');
        %## PLOT
        fig = figure(); set(gcf, 'position', [0 0 600 500]);
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
        [psd_f,freqs,~,~] = spca_get_psd_dat(tmp);
        psd_f = permute(psd_f,[2,3,1]) % epoch x chan x freq
        psd_corr_based = zeros(size(psd_f));
        psd_avg = psd_f;
        psd_corr_psc1 = zeros(size(psd_f));
        psd_baselined = zeros(size(psd_f));
        based_psd_corr = zeros(size(psd_f));
        % for cond_i = 1:length(conds)
        for s_i = 1:size(psd_f,2)
            %## STRIDE-BY-STRIDE IMPLEMENT OF SPCA
            psd_avg = psd_f(s_i,:,:);
            psd_baselined = bsxfun(@minus,psd_avg,base_txf_mean);
            %- apply spca
            [psd_cb,psd_bc,psd_cp1,psd_psc,~] = specPCAdenoising(psd_baselined,coeffs);
            psd_corr_based(s_i,:,:) = psd_cb;
            based_psd_corr(s_i,:,:) = psd_bc;
            psd_corr_psc1(s_i,:,:) = psd_cp1;
        end
        struct_out = struct('condition','all',...
            'psd_corr_based',psd_corr_based,...
            'psd_orig_avg',psd_avg,...
            'baseline_psd',base_txf_mean,...
            'psd_corr_psc1',psd_corr_psc1,...
            'psd_orig_baselined',psd_baselined,...
            'trialinfo',trialinfo,...
            'freqs',freqs,...
            'coeffs',coeffs);
        par_save(struct_out,fpath,sprintf('psd_sbs_struct.mat'));
        
        
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