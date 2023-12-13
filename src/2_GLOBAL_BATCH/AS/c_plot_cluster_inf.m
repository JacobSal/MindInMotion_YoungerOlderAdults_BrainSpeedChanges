%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/AS/run_c_plot_cluster_inf.sh

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
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
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
%% DEFINE SOURCE DIRECTORY & CD ======================================== %%
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'AS'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP
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
%% ===================================================================== %%
%## Data PATHS
fprintf('\nData PATHS\n');
%## DATASET SPECIFIC
SUBJ_PICS = {{'02','03','04','05','09','11','15','16','18','19','21','22',...
            '23','24','25','27','28','29','30','31','32','33','35','36','38'}};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
%% (PARAMETERS) ======================================================== %%
fprintf('\nData Processing Parameters\n');
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
% COND_DESIGNS = {test_1,test_2,test_3};
COND_DESIGNS = {test_2};
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
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
cluster_info_fpath = [STUDIES_DIR filesep 'as_cluster_info' filesep 'cluster_info.mat'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(EVENT_COND_COMBOS)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
%     [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
% STUDY.cluster = STUDY.etc.cluster_save;
% disp(STUDY.cluster);
% disp(STUDY.urcluster);
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
%% MAKEDESIGN & PLOT
%- get inds
[~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
%- clusters to plot
CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
%## PLOT cluster based information
mim_gen_cluster_figs(STUDY,ALLEEG,save_dir,...
    'CLUSTERS_TO_PLOT',CLUSTER_PICKS);
mim_custom_spec_plots(STUDY,ALLEEG,save_dir,...
    'CLUSTERS_TO_PLOT',CLUSTER_PICKS,...
    'SPEC_PARAMS',SPEC_PARAMS);
% for des_i = 1:length(COND_DESIGNS)
%     %- (06/22/2023) JS, AS dataset requires variable1 to be changed;
%     %'condlabel' includes:['competitive','cooperative','moving_serve','stationary_serve'
% %     [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
% %             'subjselect', {ALLEEG.subject},...
% %             'variable1','condlabel',...
% %             'values1', COND_DESIGNS{des_i});
%     %- (06/22/2023) JS, AS dataset requires variable 1 to be changed
%     %'condlabel' includes:['1Bounce_Human','2Bounce_BM','2Bounce_Human','Serve_Human'
%     [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
%         'subjselect', {ALLEEG.subject},...
%         'variable1','bounces',...
%         'values1', COND_DESIGNS{des_i});
%     for cluster_i = main_cl_inds
%         %-
%     %     [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
%     %             'clusters',cluster_i,'comps','all','subject','','freqrange', FREQ_LIMITS,'subtractsubjectmean','on');
%         %-
%         [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
%                 'clusters',cluster_i,'freqrange',SPEC_FREQLIMITS,'plotmode','condensed',...
%                 'plotconditions','together','ylim',SPEC_YLIM); %'plotgroups','together'
%         fig_i = get(groot,'CurrentFigure');
%         %- set figure line colors
%         cc = linspecer(length(COND_DESIGNS{des_i}));
%         iter = 1;
%         for i = 1:length(fig_i.Children(2).Children)
%             set(fig_i.Children(2).Children(i),'LineWidth',1.5);
%             set(fig_i.Children(2).Children(i),'Color',horzcat(cc(iter,:),0.6));
%             if iter == length(cc)
%                 iter = 1;
%             else
%                 iter = iter + 1;
%             end                
%         end
%         %- tight layout
%         % (06/22/2023) JS, seems like changing the Position property also
%         % erases the statistics bar. Not sure how to fix?
% %         set(fig_i.Children(2),'Position',[0.26,0.26,0.54,0.51]); %DEFAULT
% %         set(fig_i.Children(3),'Position',[0.26,0.2345,0.54,0.0255]); %DEFAULT
%         set(fig_i.Children(2),'Position',[0.15,0.15,0.8,0.8]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
%         set(fig_i.Children(3),'Position',[0.15,0.15-0.0255,0.8,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
%         set(fig_i.Children(1),'Location','southeast') %reset Legend
%         drawnow;
%         %- save
%         saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
%         saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
%         close(fig_i)
%     end
% end
%% Version History
%{
v1.0; (11/11/2022), JS: really need to consider updating bootstrap
    algorithm with parallel computing. Taking ~ 1 day per
    condition for all subjects and the bottle neck is entirely the
    bootstrap.

    Note: validateattributes and assert functions may be helpful
    in more clearly defining function inputs.
        e.g.  DO_PHASE_RND = true;
          errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
          validationFcn = @(x) assert(islogical(x),errorMsg);
v1.0; (12/5/2022) Need to adapt this to include all conditions
    within each SUBJ structure so connectivity can be calculated
    for the ALLEEG structure rather than the EEG structure.
    *** maybe try to ditch the SUBJ strucutre entirely for this
    round?
v1.0.01132023.0 : Initializing versioning for future iterations.
%}

