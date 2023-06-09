%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_test/3_paper_MIM_HOA/run_alleeg_spectral_run.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    DO_UNIX = false;
    PATH_EXT = 'M';
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
end
%## DEBUG: PATHROOT OVERRIDE
if DO_UNIX
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '_test' filesep '3_paper_MIM_HOA'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
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
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
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
%% ================================================================= %%
fprintf('Assigning Paths\n');
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
% SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
%             'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
%             'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
% SUBJ_MISSING_TRIAL_DATA = {'H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
%     'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
SUBJ_2HMA = {'H2017', 'H2010', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2036', 'H2037', 'H2038',...
    'H2039', 'H2041', 'H2042', 'H2052', 'H2059', 'H2062', 'H2072', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3NHMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- Subject Picks
SUBJ_PICS = {SUBJ_2HMA,SUBJ_3NHMA};
GROUP_NAMES = {'H2000''s','H3000''s'};
SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3NHMA)};
%% ===================================================================== %%
fprintf('Assigning Params\n');
%## PROCESSING PARAMS
%- compute measures for spectrum and ersp
RECOMPUTE_SPEC = false;
RECOMPUTE_ERSP = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = true; %false;
FREQ_LIMITS = (1:100);
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
STAT_ALPHA = 0.05;
%- SPEC params
SPEC_XLIM = [1,100];
SPEC_YLIM = [-30,-5]; 
SPEC_MCORRECT = 'cluster'; %'fdr';
SPEC_STAT_METHOD = 'perm'; %'parametric';
SPEC_ALPHA = 0.05; %nan();
SPEC_CONDSTATS = 'on';
SPEC_GROUPSTATS = 'off';
SPEC_NACCU = 10000;
%- ERSP PARAMS
% condstats ='on';        % ['on'|'off]
% statsMethod ='perm';    % ['param'|'perm'|'bootstrap']
% Alpha = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
% mcorrect = 'cluster'; %fdr
% groupstats = 'off';
% mode = 'fieldtrip';
% singletrials = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
TIMEWARP_NTIMES = 200; % making this too big can cause overlap between gait cyles?
% ERSP_CAXIS = [-5,10]; %[-1.5,1.5];
ERSP_CAXIS = [-1.5,1.5];
ERSP_FREQRANGE = [1,70];
DO_SUBBASELINE = false;
ERSP_MCORRECT = 'cluster'; %'fdr';
ERSP_STAT_METHOD = 'perm';
ERSP_ALPHA = 0.05; %nan();
ERSP_CONDSTATS = 'on';
ERSP_GROUPSTATS = 'off';
ERSP_NACCU = 10000;
% PAD_RATIO = 2;
% SPEC_EPOCH_LIMS = [-1,3.4]; %[-2,2];
%- datetime override
dt = '05012023_MIM_OA_subset_N85_speed_terrain_merge';
%- hard define
load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
study_fName_3 = sprintf('%s_EPOCH_study',[load_trials{:}]);
%- soft define
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
fprintf('Loading STUDY & ALLEEG\n');
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(load_trials)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_3 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_3 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_3 '.study'],'filepath',load_dir);
    end
end
%% combine groups?
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = 'Older Adults';
    STUDY.datasetinfo(subj_i).group = 'Older Adults';
end
%% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- extract per condition timewarp from event struct
%     new_warp = struct('latencies',[],'epochs',[],'eventSequence',[],'warpto',[]);
%     new_warp.latencies = cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.latencies);
%     new_warp.epochs = 1:length(new_warp.latencies);
%     new_warp.eventSequence = ALLEEG(subj_i).etc.timewarp_by_cond(1).eventSequence;
%     new_warp.warpto = ALLEEG(subj_i).timewarp.warpto; %median(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto)); % could use mean here, using median for now
    all_warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto)); % could use mean here, using median for now
    %- reassign per condition timewarping
    ALLEEG(subj_i).timewarp.warpto = all_warpto;
end
% allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
all_warpto = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    all_warpto(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = median(allWarpTo);
grandAvgWarpTo = nanmedian(all_warpto);
% grandAvgWarpTo = mean(allWarpTo);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
%     [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG,...
%                         'components',...
%                         'recompute','on',...
%                         'spec','on',...
%                         'scalp','on',...
%                         'savetrials','on',...
%                         'specparams',...
%                         {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
%                         'freqrange',FREQ_LIMITS,'logtrials','off'});
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',FREQ_LIMITS,'logtrials','off'});
    end
end
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
if RECOMPUTE_ERSP
    disp(['Grand average (across all subj) warp to: ',num2str(grandAvgWarpTo)]);
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        %- determine timewarping parameters
         if DO_TIMEWARP
            timewarp_param = EEG.timewarp.latencies;
            timewarpms_param = grandAvgWarpTo;
         else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length(FREQ_LIMITS),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',[grandAvgWarpTo(1) grandAvgWarpTo(end)],...
                    'commonbase','on','trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length(FREQ_LIMITS),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
end
%% COMPARE SPECTRUM ACROSS CONDITION
group_names = unique({ALLEEG.group});
ersp_event_names = ALLEEG(1).timewarp.eventSequence;
for des_i = 1:length(COND_DESIGNS)
    fprintf('==== Making Study Design ====\n');
    [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
        'subjselect', {ALLEEG.subject},...
        'variable1','cond',...
        'values1', COND_DESIGNS{des_i});
    for cluster_i = 2:length(STUDY.cluster)
        fprintf('\n==== Loading Cluster %i from STUDY %s ====\n',cluster_i, STUDY.name)
        fprintf('Cluster %i has subjects: ',cluster_i); fprintf('%i,',STUDY.cluster(cluster_i).sets); fprintf('\n');
        fprintf('with components: '); fprintf('%i,',STUDY.cluster(cluster_i).comps); fprintf('\n');
        %- load data into matrix using std_readdata (can load by
        %components or clusters (using clusters here).
        %% (ERSP PLOT 1) EEGLAB IMPLEMENT
        %- ersp plot per cluster per condition
        STUDY = pop_statparams(STUDY, 'condstats', 'on',...
                'groupstats','off',...
                'method','perm',...
                'singletrials','off','mode','fieldtrip','fieldtripalpha',0.05,...
                'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',10000);
        t_rng = [grandAvgWarpTo(1) grandAvgWarpTo(5)];
%         rng = [0 1320];
        STUDY = pop_erspparams(STUDY, 'subbaseline','off',...
            'timerange',t_rng, 'ersplim',ERSP_CAXIS,'freqrange',ERSP_FREQRANGE);  % 'subbaseline' - ['on'|'off'] subtract the same baseline across conditions for ERSP     
        [~] = std_erspplot(STUDY, ALLEEG,...
                            'clusters',cluster_i);
        fig_i = get(groot,'CurrentFigure');
        fig_i.Name = sprintf('Cluster %i) condition %s',cluster_i,[COND_DESIGNS{des_i}{:}]);
        saveas(fig_i,[save_dir filesep sprintf('cond_%s_ersp_%i.fig',[COND_DESIGNS{des_i}{:}],cluster_i)]);
        saveas(fig_i,[save_dir filesep sprintf('cond_%s_ersp_%i.jpg',[COND_DESIGNS{des_i}{:}],cluster_i)]);
        %% (ERSP PLOT 2) CUSTOM
        %{
        %- ERSP DATA
        [STUDY, data_ersp, times_ersp, freqs_ersp, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
                        'clusters',cluster_i,'datatype','ersp','freqrange',FREQ_LIMITS,'design',des_i);
        %- combine groups...
        if size(data_ersp,2) > 1
            tmp_ersp = cell(size(data_ersp,1),1);
            for i = 1:size(data_ersp,1)
                for j = 1:size(data_ersp,2)
                    tmp_ersp{i} = cat(3,tmp_ersp{i},data_ersp{i,j});
                end
            end
            data_ersp = tmp_ersp;
        end
        %- across condition plot
        titles_ersp = std_figtitle('threshold', STAT_ALPHA, 'mcorrect', ERSP_MCORRECT,...
                                'condstat', ERSP_CONDSTATS,...
                                'statistics', ERSP_STAT_METHOD,...
                                'condnames', COND_DESIGNS{des_i},...
                                'clustname', STUDY.cluster(cluster_i).name,...
                                'subject', [], 'datatype', 'ersp', 'plotmode', 'normal', ...
                                'effect', 'main');
        
        allersp = data_ersp;
        alltimes = times_ersp;
        allfreqs = freqs_ersp;
        %- subtract out average ersp
        if DO_SUBBASELINE
            params_ersp.baseline = [grandAvgWarpTo(1) grandAvgWarpTo(5)];
            params_ersp.basenorm = 'on';
            params_ersp.trialbase = 'off';
            params_ersp.singletrials = 'off';
            params_ersp.commonbase   = 'on';
            [allersp,basesamples,basevals] = newtimefbaseln(allersp, alltimes, params_ersp);
        else
            params_ersp.singletrials = 'off';
            allersp = cellfun(@(x)newtimefbaseln(x, alltimes, params_ersp), allersp, 'uniformoutput', false);
        end
        %- convert ersp to log decibel scale
        allersp = cellfun(@(x)10*log10(x), allersp, 'uniformoutput', false);
        %- compute statistics
        [pcond_ersp,pgroup_ersp,pinter_ersp,~,~,~] = std_stat(allersp,...
                                        'condstats', ERSP_CONDSTATS,...
                                        'groupstats',ERSP_GROUPSTATS,...
                                        'method',ERSP_STAT_METHOD,...
                                        'naccu',ERSP_NACCU,...
                                        'alpha',ERSP_ALPHA,...
                                        'mcorrect',ERSP_MCORRECT,'mode','fieldtrip');
        std_plottf(alltimes, allfreqs, allersp,...
                           'datatype', 'ersp', ...
                           'groupstats', pgroup_ersp,...
                           'condstats', pcond_ersp,...
                           'interstats', pinter_ersp,...
                           'plotmode', 'normal',...
                           'titles', titles_ersp,...
                           'events', {},... %ersp_event_names,...
                           'unitcolor', 'dB',...
                           'chanlocs', struct([]),... %ALLEEG(1).chanlocs,...
                           'caxis',ERSP_CAXIS,...
                           'ersplim',[],'threshold',nan(),'effect','main','maskdata','off','averagemode','rms');
        fig_i = get(groot,'CurrentFigure');
        fig_i.Name = sprintf('Cluster %i) ERSP averages across Conditions',cluster_i);
%         saveas(fig_i,[save_dir filesep sprintf('btwnCond_ersp_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('btwnCond_ersp_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        %}
        %% (SPECTOGRAM PLOT 1) EEGLAB IMPLEMENT
        %{
        STUDY = pop_statparams(STUDY, 'condstats', 'on',...
                'groupstats','off',...
                'method','perm',...
                'singletrials','off','mode','fieldtrip','fieldtripalpha',0.05,...
                'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',10000);
        [STUDY, ~, ~, ~, ~, ~] = std_specplot(STUDY, ALLEEG,...
                'clusters',cluster_i,'freqrange',ERSP_FREQRANGE);
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        %}
        %% (SPECTOGRAM PLOT 2) CUSTOM
        %{
        %- SPECTRUM DATA
        [STUDY, data_spec, freqs_spec, ~, ~, ~] = std_readdata(STUDY,ALLEEG,...
                        'clusters',cluster_i,'datatype','spec','freqrange',FREQ_LIMITS,'design',des_i); 
        %- combine groups...
        if size(data_spec,2) > 1
            tmp_spec = cell(size(data_spec,1),1);
            for i = 1:size(data_spec,1)
                for j = 1:size(data_spec,2)
                    tmp_spec{i} = cat(3,tmp_spec{i},data_spec{i,j});
                end
            end
            data_spec = tmp_spec;
        end
        %## Preprocessing data across subjects
        disp('running stats');
        % naccu > 200 for 0.01 alpha and naccu > 2000 for 0.001 alpha
        [pcond,pgroup,pinter,statcond,statgroup,statinter] = std_stat(data_spec, 'condstats', SPEC_CONDSTATS,...
                                    'groupstats',SPEC_GROUPSTATS,'method',SPEC_STAT_METHOD,'naccu',SPEC_NACCU,...
                                    'alpha',SPEC_ALPHA,'mcorrect',SPEC_MCORRECT,'mode','eeglab');
        %## generate mean and standard deviation across group for component and condition
        color_g1 = {'r','b','g','c','y'};
        color_g2 = {'--r','--b','--g','--c','--y'};
        figure;
        title(sprintf('Cluster %i) Mean and Standard Dev Across Conditions',cluster_i));
        hold on;
        for group_i = 1:size(data_spec)
            tmpdata = data_spec{cond_i};
            for cond_i = 1:size(data_spec)
                y = mean(tmpdata,2); % your mean vector;
                x = freqs_spec;
                std_dev = std(tmpdata,1,2);
                curve2 = (y + std_dev);
                curve1 = (y - std_dev);
                x2 = [x', fliplr(x')];
                inBetween = [curve1',fliplr(curve2')];
        %         h = fill(x2, inBetween, sprintf('%s',color{cond_i}),'DisplayName',TRIAL_TYPES{cond_i});
                h.FaceAlpha = 0.35;
                plot(x, y,sprintf('%s',color{cond_i}),'LineWidth',2,'DisplayName',COND_DESIGNS{des_i}{cond_i});
            end
        end
    %     plot(x, 5*(newpcond < STAT_ALPHA),'LineWidth',2,'DisplayName',sprintf('pval < %0.02f',STAT_ALPHA));
        xlabel('Frequency (Hz)')
        ylabel('Log Power 10*log10(uV^2)');
        legend()
        xlim(SPEC_XLIM)
        ylim(SPEC_YLIM)
        hold off;
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_dir filesep sprintf('btwnCond_psd_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('btwnCond_psd_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        %}
    end
end
%%
%- Dipole Plot
% std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',CLUSTER_ITERS,'mode','multicolor');
% std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters','all','mode','multicolor');
% view([45,0,0])
% view([0,-45,0])
% view([0,0,45])
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}