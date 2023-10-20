%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/AS/run_c_ersp_plotting_new.sh

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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'AS'];
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
DATA_SET = 'AS_dataset';
%- compute measures for spectrum and ersp
RECOMPUTE_SPEC = false;
RECOMPUTE_ERSP = false;
DO_TIMEWARP = false;
DO_BASELINE_CORRECTION = false; %false;
%- statistics & conditions
EVENTS_TIMEWARP = {'Subject_hit','Subject_receive','Subject_hit'};
% COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
COND_CHARS =  {'2Bounce_Human','2Bounce_BM'};
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
test_1 = {'1Bounce_Human','2Bounce_Human'};
test_2 = {'2Bounce_Human','2Bounce_BM'};
test_3 = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'};
COND_EVENT_CHAR = 'bounces';
baseline_rng = [0,150];
COND_DESIGNS = {test_1,test_2,test_3};
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
% (10/20/2023) borrowing parameters from MIM study
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
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
dt = '08182023_bounces_1h_2h_2bm_JS';
%## Soft Define
%- combinations of events and conditions
EVENT_COND_COMBOS = cell(length(COND_CHARS)*length(EVENT_CHARS),1);
cnt = 1;
for cond_i = 1:length(COND_CHARS)
    for event_i = 1:length(EVENT_CHARS)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',COND_CHARS{cond_i},EVENT_CHARS{event_i});
        cnt = cnt + 1;
    end
end
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
study_load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
fprintf('Loading STUDY & ALLEEG\n');
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(EVENT_COND_COMBOS)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([study_load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',study_load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',study_load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cluster] = eeglab_get_cluster_comps(STUDY);
end
%% CALCULATE GRANDAVERAGE WARPTO
%{
%## (OPT 1)
tmp_warpto = cell(length(ALLEEG),1);
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).timewarp.warpto = [];
    ALLEEG(subj_i).timewarp.latencies = [];
    ALLEEG(subj_i).timewarp.epochs = [];
    ALLEEG(subj_i).timewarp.eventSequence = [];
    tmp = zeros(length(ALLEEG(subj_i).epoch),3);
    for epoch_i = 1:length(ALLEEG(subj_i).epoch)
         tmp(epoch_i,:) = [ALLEEG(subj_i).epoch(epoch_i).eventtimewarplatencies_1{1}, ALLEEG(subj_i).epoch(epoch_i).eventtimewarplatencies_2{1}, ALLEEG(subj_i).epoch(epoch_i).eventtimewarplatencies_3{1}];
    end
    tmp_warpto{subj_i} = tmp;
    ALLEEG(subj_i).timewarp.latencies = tmp;
    ALLEEG(subj_i).timewarp.epochs = 1:length(ALLEEG(subj_i).epoch);
    ALLEEG(subj_i).timewarp.warpto = median(tmp);
    ALLEEG(subj_i).timewarp.eventSequence = EVENTS_TIMEWARP;
end
%## (OPT 1)
% for subj_i = 1:length(ALLEEG)
%     %- assign percondition timewarping
%     ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
% %     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
% end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
% allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
grandAvgWarpTo = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
%##
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
% b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)-1000*(1/ALLEEG(1).srate)];
b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(end)];
% ERSP_CROP_TIMES=[grandAvgWarpTo(1)+abs(ALLEEG(1).etc.epoch.epoch_limits(1))*1000, grandAvgWarpTo(5)];
ERSP_CROP_TIMES=[grandAvgWarpTo(1), grandAvgWarpTo(end)];
fprintf('Using timewarp limits: [%0.1g,%0.1f]\n',b_lims(1),b_lims(2));

%}
%## ersp plot per cluster per condition
%%
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange); %,'timerange',ERSP_CROP_TIMES);
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
    TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
%     disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
%     parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/3))
    for subj_i = 1:length(ALLEEG)
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
%% (ERSP PLOT 1) EEGLAB IMPLEMENT: ACROSS CONDITIONS
% group_names = unique({ALLEEG.group});
% ersp_event_names = ALLEEG(1).timewarp.eventSequence;
% for des_i = 1:length(COND_DESIGNS)
%     fprintf('==== Making Study Design ====\n');
%     [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
%         'subjselect', {ALLEEG.subject},...
%         'variable1','cond',...
%         'values1', COND_DESIGNS{des_i});
%     for cluster_i = 2:length(STUDY.cluster)
%         fprintf('\n==== Loading Cluster %i from STUDY %s ====\n',cluster_i, STUDY.name)
%         fprintf('Cluster %i has subjects: ',cluster_i); fprintf('%i,',STUDY.cluster(cluster_i).sets); fprintf('\n');
%         fprintf('with components: '); fprintf('%i,',STUDY.cluster(cluster_i).comps); fprintf('\n');
%         %- load data into matrix using std_readdata (can load by
%         %components or clusters (using clusters here).
%         [~] = std_erspplot(STUDY, ALLEEG,...
%                             'clusters',cluster_i);
%         fig_i = get(groot,'CurrentFigure');
%         fig_i.Name = sprintf('Cluster %i) condition %s',cluster_i,[COND_DESIGNS{des_i}{:}]);
%         saveas(fig_i,[save_dir filesep sprintf('cond_%s_ersp_%i.fig',[COND_DESIGNS{des_i}{:}],cluster_i)]);
%         saveas(fig_i,[save_dir filesep sprintf('cond_%s_ersp_%i.jpg',[COND_DESIGNS{des_i}{:}],cluster_i)]);
%     end
% end
%% (CUSTOM) COMPARE ERSP ACROSS CONDITION
%{
group_names = unique({ALLEEG.group});
ersp_event_names = ALLEEG(1).timewarp.eventSequence;
%}
%## LOCAL SWITCHES
% ERSP_SINGLETRIALS = 'on';
XAXIS_LABEL = 'ms';
COLORAXIS_LABEL = 'dB';
ERSP_FREQSCALE = 'log'; % 'native', 'LOG'
ERSP_CHANLOCS = struct('labels', {});
%%
for des_i = 1:length(COND_DESIGNS)
    fprintf('==== Making Study Design ====\n');
    [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
        'subjselect', {ALLEEG.subject},...
        'variable1','cond',...
        'values1', COND_DESIGNS{des_i});
    for cluster_i = main_cl_inds
        %% (ERSP PLOT 2) CUSTOM
        %- ERSP DATA
%         [STUDY,data_ersp,times_ersp,freqs_ersp,events_ersp,params_ersp] = mim_read_ersp(STUDY,cluster_i,...
%             des_i,b_lims,ERSP_FREQLIMITS,ERSP_SINGLETRIALS);
        %- load ERSP data using std_readdat
        %* (06/04/2023) JS, could try and load a subject at a time,
        % baseline, then stack them together to plot)
        %* (06/10/2023) JS, removing timerange loading, need to custom
        %prune later
%         [STUDY, data_ersp, times_ersp, freqs_ersp, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
%                         'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
%                         'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
%                         'design',des_i,'timerange',b_lims);
        [STUDY, ersp_data, ersp_times, ersp_freqs, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
                        'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
                        'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
                        'design',des_i);
        %- across condition plot
        titles_ersp = std_figtitle('threshold', STAT_ALPHA, 'mcorrect', ERSP_MCORRECT,...
                                'condstat', ERSP_CONDSTATS,...
                                'statistics', 'Fieldtrip Montecarlo w/ Cluster',... %STAT_MODE,... %'fieldtrip cluster',... %ERSP_STAT_METHOD,...
                                'condnames', COND_DESIGNS{des_i},...
                                'clustname', STUDY.cluster(cluster_i).name,...
                                'subject', [], 'datatype', 'ersp', 'plotmode', 'normal', ...
                                'effect', 'main');
        
        allersp = ersp_data;
        alltimes = ersp_times;
        allfreqs = ersp_freqs; 
        %## subtract condition based mean ersp
        % (05/31/2023): JS I don't believe this to generate a good result.
        % Problem with logging the output after baseline subtraction.
        % (06/02/2023): JS I'm trying this again, going to try and extract
        % each epoch using event structure and sort into conditions.
        for cond_i = 1:length(allersp)
            %- interesting way to baseline using mean across trials and
            %time
            tmp = mean(allersp{cond_i},3);
%             tmp_std = std(allersp{cond_i},[],3);
            tmp = mean(tmp,2);
            tmp_std = std(tmp,[],2);
            tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
            tmp_std = repmat(tmp_std,1,size(allersp{cond_i},2));
            allersp{cond_i} = allersp{cond_i} - tmp;
        end
        %- compute statistics
        %* (06/10/2023) JS, making this fieldtrip method only
%         [pcond_ersp,pgroup_ersp,pinter_ersp,stat_cond,stat_group,stat_inter] = std_stat(allersp,...
%                                         'condstats', ERSP_CONDSTATS,...
%                                         'groupstats',ERSP_GROUPSTATS,...
%                                         'method',ERSP_STAT_METHOD,...
%                                         'naccu',ERSP_NACCU,...
%                                         'alpha',ERSP_ALPHA,...
%                                         'mcorrect',ERSP_MCORRECT,'mode','fieldtrip');
        %## Cropping Indicies
        crop_inds = (alltimes>=ERSP_CROP_TIMES(1) & alltimes<=ERSP_CROP_TIMES(2));
        alltimes = alltimes(crop_inds);
        for cond_i = 1:length(allersp)
            allersp{cond_i} = allersp{cond_i}(:,crop_inds,:);
        end
        [pcond_ersp,pgroup_ersp,pinter_ersp,stat_cond,stat_group,stat_inter] = std_stat(allersp,...
                'condstats', ERSP_CONDSTATS,...
                'groupstats',ERSP_GROUPSTATS,...
                'method',ERSP_STAT_METHOD,...
                'naccu',ERSP_NACCU,...
                'alpha',ERSP_ALPHA,...
                'mcorrect',ERSP_MCORRECT,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
                'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
        %- plot
        std_plottf(alltimes, allfreqs, allersp,...
                           'datatype', 'ersp', ...
                           'groupstats', stat_group,...
                           'condstats', stat_cond,...
                           'interstats', stat_inter,...
                           'plotmode', 'normal',...
                           'unitx',XAXIS_LABEL,...
                           'titles', titles_ersp,...
                           'events', {},... %ersp_event_names,...
                           'unitcolor', COLORAXIS_LABEL,...
                           'chanlocs', ERSP_CHANLOCS,... %ALLEEG(1).chanlocs.labels,...
                           'caxis',ERSP_CAXIS,...
                           'freqscale',ERSP_FREQSCALE,...
                           'ersplim',[],...
                           'threshold',ERSP_ALPHA,'effect','main',...
                           'maskdata','off','averagemode','rms');
        fig_i = get(groot,'CurrentFigure');
        fig_i.Name = sprintf('Cluster %i) ERSP averages across Conditions',cluster_i);
%         for i = 3:length(fig_i.Children)
%         end
%         fig_i.Children(4).YScale = 'log';
        saveas(fig_i,[save_dir filesep sprintf('custom2_Cond_ersp_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('custom2_Cond_ersp_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        
    end
end
%% (CL) CUSTOM PLOTS
if ~exist([save_dir filesep 'custom_cl_plot'],'dir')
    mkdir([save_dir filesep 'custom_cl_plot']);
end
SUB_FREQ_LIMS = [4,60];
ALPHA = 0.05;
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);
for des_i = 1:length(COND_DESIGNS)
    fprintf('==== Making Study Design ====\n');
    %## conditions across all 
    [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
        'subjselect', {ALLEEG.subject},...
        'variable1',COND_EVENT_CHAR,...
        'values1',COND_DESIGNS{des_i});
    for j = 1:length(valid_cluster)
%     parfor (j = 1:length(valid_cluster),POOL_SIZE)
        cluster_i = valid_cluster(j);
        [~, ersp_data, ersp_times, ersp_freqs, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
                        'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
                        'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
                        'design',des_i);
        %- across condition plot
%         titles_ersp = std_figtitle('threshold', STAT_ALPHA, 'mcorrect', ERSP_MCORRECT,...
%                                 'condstat', ERSP_CONDSTATS,...
%                                 'statistics', STAT_MODE,... %STAT_MODE,... %'fieldtrip cluster',... %ERSP_STAT_METHOD,...
%                                 'condnames', COND_DESIGNS{des_i},...
%                                 'clustname', STUDY.cluster(cluster_i).analabel{1},...
%                                 'subject', [], 'datatype', 'ersp', 'plotmode', 'normal', ...
%                                 'effect', 'main');
        allersp = ersp_data;
        alltimes = ersp_times;
        allfreqs = ersp_freqs;
        %## BASELINE
        baseidx = alltimes>=baseline_rng(1) & alltimes>=baseline_rng(2);
        %- 
        freqidx = find(allfreqs>=SUB_FREQ_LIMS(1) & allfreqs<=SUB_FREQ_LIMS(2));
        subjs_cl = unique(STUDY.cluster(cluster_i).sets); %Subjects in this cluster
        mean_subj = cell(length(allersp),1);
        for cond_i = 1:length(allersp)
            mean_subj{cond_i,1}(:,:,:) = zeros(size(allersp{cond_i},1),size(allersp{cond_i},2),length(subjs_cl));
        end
        based_ersps = allersp;        
        p = 1;
        sub = unique(STUDY.cluster(cluster_i).sets);
        for n = 1:length(unique(STUDY.cluster(cluster_i).sets))    
            comp_ind = STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
            % comp not using
            if cluster_i == 1 %bad comp
                for cond_i = 1:length(COND_DESIGNS{des_i})
                    based_ersps{cond_i}(:,:,p) = nan(size(allersp{cond_i},1),size(allersp{cond_i},2),1);
                    mean_subj{cond_i}(:,:,n) = nanmean(based_ersps{cond_i}(:,:,p:p + length(comp_ind)-1),3);     
                end
            else
                for cond_i = 1:length(COND_DESIGNS{des_i})
                    mean_subj{cond_i}(:,:,n) = nanmean( based_ersps{cond_i}(:,:,p:p + length(comp_ind)-1),3);
                end
            end
            p = p+length(comp_ind);
        end
        alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign).variable(1).value, ...
            'clustname', sprintf('CL_%i',cluster_i));
        % reorganize allerspdata
        cluster_allcomp_ersp_crop = cell(length(mean_subj),1);
        cluster_allcomp_ersp = cell(length(mean_subj),1);
        for cond_i = 1: length(mean_subj)
            erspdata = mean_subj{cond_i}(:,:,:);
            baseidx = find(alltimes>=grandAvgWarpTo(1) & alltimes<=grandAvgWarpTo(end));                
            baseline_allcomp = mean(erspdata(:,baseidx,:),2); % mean power for each person
            baseline = mean(baseline_allcomp,3);%mean power across participant
            cluster_allcomp_ersp{cond_i,1} = mean_subj{cond_i}(:,:,:)-repmat(baseline_allcomp,1,length(alltimes));% subtract baseline for each person
        %             cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline,1,length(alltimes)),3);
            cluster_allcomp_ersp_crop{cond_i,1} = cluster_allcomp_ersp{cond_i,1}(:,baseidx);
        end
        climMat = [min(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all') max(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all')];
        climMat_max = max(abs(climMat));
        %## (PLOT) Paper Figure for YA paper and IEEE NER - significance masked ERSP for high terrain
        freqidx = find(allfreqs>=SUB_FREQ_LIMS(1) & allfreqs<=SUB_FREQ_LIMS(2));
        figure('color','white','position',[200 200 700 150],'renderer','Painters');
        for cond_i = 1:length(allersp)
            if ~isnan(ALPHA)
                curr_ersp_temp = cluster_allcomp_ersp{cond_i,1}(freqidx,baseidx,:);% this is already sub baseline
                curr_ersp_temp_mean = mean(curr_ersp_temp,3);
                surro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
                for n = 1:2000
                    bootLatency = randi(size(curr_ersp_temp,2),[size(curr_ersp_temp,2),1]); %random time sample
                    bootFreq = 1:size(curr_ersp_temp,1);
                    bootIc = 1:size(curr_ersp_temp,3); 
                    tmpSurro = mean(curr_ersp_temp(bootFreq,bootLatency,bootIc),3);
                    surro(:,:,n) = tmpSurro; %save 2000 iterations of surrogates 
                end
                bootSurro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
                for n = 1:2000
                    bootIdx  = randi(2000,[size(curr_ersp_temp,3),1]);
                    tmpSurro = mean(surro(:,:,bootIdx),3);
                    bootSurro(:,:,n) = tmpSurro;
                end
                pvalMap = stat_surrogate_pvals(bootSurro,curr_ersp_temp_mean,'both');
                pvalMap(pvalMap>1)=1; 
                [p_masked, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalMap,0.05,'pdep',1);
                % debri removal
                [labelMap,uniqueLabelNum] = bwlabeln(p_masked);
                tmpDisp = sort(labelMap(:),'descend');
                [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
                sortOccurrence = sort(occurrence,'descend');
                threshold = 1000;
                threshOccurrence = occurrence;
                threshIdx = find(threshOccurrence<threshold);
                kMask = ismember(labelMap,idx(threshIdx));
                finalMask = p_masked-kMask;

                clust_ersp = curr_ersp_temp_mean; 
                clust_maskedersp = clust_ersp; 
                clust_maskedersp(~finalMask) = 0;
            else
                clust_ersp = mean(cluster_allcomp_ersp_mean{cond_i},3);
                clust_maskedersp = clust_ersp;
            end   

            subplot(1,length(allersp),cond_i)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(clust_maskedersp))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(clust_maskedersp ~=0 ) = 0; %0 is significant? 1 is not? 
            contourf(alltimes(baseidx), allfreqs(freqidx), clust_ersp,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes(baseidx),allfreqs(freqidx),clust_maskedersp,'AlphaData',faceAlpha_mask);
            %- add vertical line
            xline(gca,grandAvgWarpTo(1),'k--');
            xline(gca,grandAvgWarpTo(2),'k--');
            xline(gca,grandAvgWarpTo(3),'k--');
%             vline([grandAvgWarpTo(2) grandAvgWarpTo(3) grandAvgWarpTo(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',[-climMat_max, climMat_max],'xlim',[grandAvgWarpTo(1) grandAvgWarpTo(end)],...
                'ydir','norm','ylim',[allfreqs(1) SUB_FREQ_LIMS(2)],'yscale','log')
            if SUB_FREQ_LIMS(2) == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
            elseif SUB_FREQ_LIMS(2) == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
            end
            if cond_i == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
            title(sprintf('Condition %s',COND_DESIGNS{des_i}{cond_i}));
            set(gca,'xtick',grandAvgWarpTo,'xticklabel',{'rec','hit','rec'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
        end
        hp4 = get(subplot(1,length(allersp),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.008  hp4(4)-0.071]);
        c.Limits = [-climMat_max, climMat_max];
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.2;
        hL.Position(2) = .13;
        c.Position(2) = .10;
        c.Position(4) = .5; 
        saveas(gcf,[save_dir filesep 'custom_cl_plot' filesep sprintf('allerspdata_within_cl%i_stats_lim%i.fig',cluster_i,SUB_FREQ_LIMS(2))]);
        saveas(gcf,[save_dir filesep 'custom_cl_plot' filesep sprintf('allerspdata_within_cl%i_stats_lim%i.pdf',cluster_i,SUB_FREQ_LIMS(2))]);
        %## subject plots
        freqidx = find(allfreqs>=SUB_FREQ_LIMS(1) & allfreqs<=SUB_FREQ_LIMS(2));
        for i = 1:length(STUDY.cluster(cluster_i).comps)
            ic = STUDY.cluster(cluster_i).comps(i);
            sub = STUDY.datasetinfo(STUDY.cluster(cluster_i).sets(i)).subject; 

            baseidx = find(alltimes>=grandAvgWarpTo(1) & alltimes<=grandAvgWarpTo(5));
            erspdata = allersp{1}(:,:,i);
            baseline = mean(erspdata(:,baseidx,:),2);
            curr_ersp = erspdata(:,:,:)-repmat(baseline,1,length(alltimes));
            curr_ersp = mean(curr_ersp,3);
            curr_maskedersp = curr_ersp;
            %- plot
            figure('renderer','Painters');
            tftopo(curr_maskedersp,alltimes,allfreqs,'limits',... 
                [grandAvgWarpTo(1) grandAvgWarpTo(end) nan nan nan nan],...
                'vert',grandAvgWarpTo(1:5),'logfreq','native');
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
            title(strcat({'Cluster '},num2str(cluster_i)));
            cbar('vert');
            xline(gca,grandAvgWarpTo(1),'k--');
            xline(gca,grandAvgWarpTo(2),'k--');
            xline(gca,grandAvgWarpTo(3),'k--');
            xline(gca,grandAvgWarpTo(4),'k--');
            xline(gca,grandAvgWarpTo(5),'k--');
            alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, 'clustname', STUDY.cluster(cluster_i).name,...
                'subject', sub, 'compnames', num2str(ic));
            %- reorganize allerspdata
            single_comp_ersp = cell(length(allersp),1);
            single_comp_ersp_crop = cell(length(allersp),1);
            for cond_i = 1: length(allersp)
                erspdata = allersp{cond_i}(:,:,i);
                baseidx = find(alltimes>=grandAvgWarpTo(1) & alltimes<=grandAvgWarpTo(5));                
                baseline = mean(erspdata(:,baseidx,:),2);                
                single_comp_ersp{cond_i,1} = mean(allersp{cond_i}(:,:,i)-repmat(baseline,1,length(alltimes)),3);
                single_comp_ersp_crop{cond_i,1} = single_comp_ersp{cond_i,1}(:,baseidx);
            end
            std_plottf(alltimes(baseidx),allfreqs, single_comp_ersp_crop, 'datatype','ersp', 'plotmode','normal','titles',alltitles)
            %- save
            saveas(gcf,[save_dir filesep 'custom_cl_plot' filesep sprintf('%s_within_eeglab_%i_stats_lim%i.pdf',sub,cluster_i,SUB_FREQ_LIMS(2))]);
        end
        close all
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}