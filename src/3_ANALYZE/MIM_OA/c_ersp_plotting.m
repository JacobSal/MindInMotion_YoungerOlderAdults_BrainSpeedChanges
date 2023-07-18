%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_c_ersp_plotting.sh

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
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
RECOMPUTE_SPEC = false;
RECOMPUTE_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; %false;
% DO_CUSTOM_SUBBASELINE = true;
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
SPEC_FREQLIMITS = [1,70];
% SPEC_MCORRECT = 'cluster'; %'fdr';
% SPEC_STAT_METHOD = 'perm'; %'parametric';
% SPEC_ALPHA = 0.05; %nan();
% SPEC_CONDSTATS = 'on';
% SPEC_GROUPSTATS = 'off';
% SPEC_NACCU = 2000;
%- ERSP PARAMS
% STAT_SINGLETRIALS = 'on'; %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
STAT_MODE = 'fieldtrip'; %[
FIELDTRIP_METHOD = 'montecarlo';
% TIMEWARP_NTIMES = [];
ERSP_CAXIS = [-2,2]; %[-1.5,1.5];
ERSP_FREQLIMITS = [3,60];
ERSP_MCORRECT = 'cluster'; %'fdr';
ERSP_STAT_METHOD = 'perm'; % ['param'|'perm'|'bootstrap']
ERSP_ALPHA = 0.05; % [NaN|alpha], Significance threshold (0<alpha<<1)
ERSP_CONDSTATS = 'on'; % ['on'|'off]
ERSP_GROUPSTATS = 'off';
ERSP_SUBBASELINE = 'off'; %['on'|'off'];
ERSP_SINGLETRIALS = 'off'; %['on'|'off'], % this should be off for GAIT cycle analysis. 
ERSP_NACCU = 2000;
%- datetime override
% dt = '05012023_MIM_OA_subset_N85_speed_terrain_merge';
% dt = '05192023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06122023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06282023_MIM_OAN79_subset_prep_verified_gait';
dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
%## Soft Define
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
fprintf('Loading STUDY & ALLEEG\n');
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(TRIAL_TYPES)*4]);
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
    %## load chang's algorithmic clustering
    %* cluster parameters
    pick_cluster = 15;
    clustering_weights.dipoles = 1;
    clustering_weights.scalp = 0;
    clustering_weights.ersp = 0;
    clustering_weights.spec = 0;
    cluster_alg = 'kmeans';
    do_multivariate_data = 1;
    evaluate_method = 'min_rv';
    clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
        '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
        '_spec_',num2str(clustering_weights.spec)];
    %* load cluster information
    cluster_load_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
    outputdir = [cluster_load_dir filesep clustering_method,...
        filesep num2str(pick_cluster)];
    tmp = load([outputdir filesep sprintf('cluster_update_%i.mat',pick_cluster)]);
    cluster_update = tmp.cluster_update;
    STUDY.cluster = cluster_update;
    %- get inds
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
%% combine groups?
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = 'Older Adults';
    STUDY.datasetinfo(subj_i).group = 'Older Adults';
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
grandAvgWarpTo = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
% b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)-1000*(1/ALLEEG(1).srate)];
b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)];
% ERSP_CROP_TIMES=[grandAvgWarpTo(1)+abs(ALLEEG(1).etc.epoch.epoch_limits(1))*1000, grandAvgWarpTo(5)];
ERSP_CROP_TIMES=[grandAvgWarpTo(1), grandAvgWarpTo(5)];
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',b_lims(1),b_lims(2));
disp(grandAvgWarpTo);
%## ersp plot per cluster per condition
%- params (06/09/2023) JS
% STUDY = pop_statparams(STUDY,'condstats',ERSP_CONDSTATS,...
%         'groupstats',ERSP_GROUPSTATS,...
%         'method',ERSP_STAT_METHOD,...
%         'singletrials',STAT_SINGLETRIALS,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
%         'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
% STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_SUBBASELINE,...
%       'timerange',b_lims,'ersplim',ERSP_CAXIS,'freqrange',ERSP_FREQLIMITS);
%- params (06/10/2023) JS, removing timerange loading, consolidating single
%trials stats variable
STUDY = pop_statparams(STUDY,'condstats',ERSP_CONDSTATS,...
        'groupstats',ERSP_GROUPSTATS,...
        'method',ERSP_STAT_METHOD,...
        'singletrials',ERSP_SINGLETRIALS,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
        'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_SUBBASELINE,...
      'ersplim',ERSP_CAXIS,'freqrange',ERSP_FREQLIMITS);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
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
                    'freqrange',SPEC_FREQLIMITS,'logtrials','on'});
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
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',b_lims,...
                    'commonbase','on','trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
end
EEG = [];
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
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = tmp_group_unif{subj_i};
    STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
end
for des_i = 1:length(COND_DESIGNS)
    fprintf('==== Making Study Design ====\n');
    [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
        'subjselect', {ALLEEG.subject},...
        'variable1','cond',...
        'values1', COND_DESIGNS{des_i});
    for cluster_i = 2:length(STUDY.cluster)
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
                                'statistics', STAT_MODE,... %STAT_MODE,... %'fieldtrip cluster',... %ERSP_STAT_METHOD,...
                                'condnames', COND_DESIGNS{des_i},...
                                'clustname', STUDY.cluster(cluster_i).name,...
                                'subject', [], 'datatype', 'ersp', 'plotmode', 'normal', ...
                                'effect', 'main');
        
        allersp = ersp_data;
        alltimes = ersp_times;
        allfreqs = ersp_freqs;
        %## subtract out average ersp (SINGLTRIALS MUST BE OFF)
        % (06/02/2023), JS this doesnt seem to sufficiently baseline
        % conditions, moving toward a more analytic approach
        %{
        %- Run#:(baseline,basenorm,trialbase,singletrials,commonbase);"RESULT NOTES"
        % Run1:(set,'on','off','off','on');""
        % Run2: (default,'on','off','off','on');"this seems to work well, make sure NOT to log output"
        params_ersp.baseline = [grandAvgWarpTo(1) grandAvgWarpTo(5)-1000*(1/ALLEEG(1).srate)];
        params_ersp.basenorm = 'on';
        params_ersp.trialbase = 'off';
        params_ersp.singletrials = 'off';
        params_ersp.commonbase   = 'on';
        [allersp,basesamples,basevals] = newtimefbaseln(allersp, alltimes, params_ersp);
        %}
        %## subtract out average ersp (SINGLETRIALS MUST BE ON)
        % (06/09/2023) JS, tried turning singletrials off and it definetly
        % changes the result. Raises the question: should be average across
        % subjects then baseline, or baseline across trials (i.e., gait
        % cycles) then average for each condition? It would, seems like
        % std_readdata takes the mean ersp for each subject when
        % singletrials is 'off';
        % (06/09/2023) JS, In revisting this I'm not sure if altering the
        % parameters really baselines the allersp variable. Going to try
        % and see if this changes the ersp I get when I use my custom
        % baselining as well. RESULT: not really?
%         params_ersp.powerbase = nan();
%         params_ersp.baseline = ERSP_CROP_TIMES;
%         params_ersp.basenorm = 'off';
%         params_ersp.commonbase   = 'off';
%         params_ersp.singletrials = 'on';
%         params_ersp.trialbase = 'off';
%         allersp = cellfun(@(x)newtimefbaseln(x, alltimes, params_ersp), allersp,...
%                 'uniformoutput', false);
%         [allersp,basln,mbase] = newtimefbaseln(allersp, alltimes, params_ersp);
        
        %## subtract aggregate baseline across all conds
%         trial_n = 0;
%         trial_lower = 1;
%         trial_upper = 0;
%         %- extract
% %         tmp_ersp = zeros(size(allersp{1},1),size(allersp{1},2),size(allersp{1},3)*length(allersp));
% %         for i = 1:length(allersp)
% %             trial_upper = trial_upper + size(allersp{i},3);
% %             tmp_ersp(:,:,trial_lower:trial_upper) = allersp{i};
% %             trial_lower = trial_upper+1;
% %         end
%         %- extract
%         for i = 1:length(allersp)
%             trial_n = trial_n + size(allersp{i},3);
%         end
%         tmp_ersp = [allersp{:}];
%         tmp_ersp = reshape(tmp_ersp,size(allersp{1},1),size(allersp{1},2),trial_n);
%         %- average
%         tmp = mean(tmp_ersp,3);
%         tmp = mean(tmp,2);
% %         tmp = mean(tmp_ersp,3);
% %         tmp = median(tmp,2);
%         tmp = repmat(tmp,1,size(allersp{1},2)); 
        %## subtract condition based mean ersp
        % (05/31/2023): JS I don't believe this to generate a good result.
        % Problem with logging the output after baseline subtraction.
        % (06/02/2023): JS I'm trying this again, going to try and extract
        % each epoch using event structure and sort into conditions.
        
        %## Cropping Indicies
        crop_inds = (alltimes>=ERSP_CROP_TIMES(1) & alltimes<=ERSP_CROP_TIMES(2));
        alltimes = alltimes(crop_inds);
        for cond_i = 1:length(allersp)
            allersp{cond_i} = allersp{cond_i}(:,crop_inds,:);
        end
        %## Across Design Baselineing
        %{
        tmp_allersp = [allersp{:}];
        tmp_allersp = reshape(tmp_allersp,size(allersp{1},1),size(allersp{1},2),size(allersp{1},3)*length(allersp));
        tmp = mean(tmp_allersp,3);
        tmp = mean(tmp,2);
        tmp_std = std(tmp,[],2);
        tmp = repmat(tmp,1,size(tmp_allersp,2)); 
        tmp_std = repmat(tmp_std,1,size(tmp_allersp,2));
        for cond_i = 1:length(allersp)
            allersp{cond_i} = allersp{cond_i} - tmp;
        end
        %}
        %## Within Each Condition Baselineing
        for cond_i = 1:length(allersp)
            %- interesting way to baseline using mean across trials and
            %time
            tmp = mean(allersp{cond_i},3);
%             tmp_std = std(allersp{cond_i},[],3);
            tmp = mean(tmp,2);
            tmp_std = std(tmp,[],2);
            tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
            tmp_std = repmat(tmp_std,1,size(allersp{cond_i},2));
            %- interesting way to baseline using median across trials and
            %time
%             tmp = median(allersp{cond_i},3);
%             tmp = median(tmp,2);
%             tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
            %- interesting way to baseline using median across trials and
            %time
%             tmp = mean(allersp{cond_i},3);
%             tmp = median(tmp,2);
%             tmp = repmat(tmp,1,size(allersp{cond_i},2)); 
            %- baseline using mean across condition trials
%             tmp = mean(allersp{cond_i},3);
            %- baselin using median across condition trials
%             tmp = median(allersp{cond_i},3);
            %- baseline using mean/median singular value
            % (06/09/2023) JS, this doesn't seem to produce a
            % clean/interpretable result.
%             tmp = mean(allersp{cond_i},3);
%             tmp = mean(tmp,1);
%             tmp = mean(tmp,2);
%             tmp = repmat(tmp,size(allersp{cond_i},1),size(allersp{cond_i},2)); 
            %- baseline using mean/median singular value
%             tmp = mean(allersp{cond_i},3);
%             tmp = mean(tmp,1);
%             tmp = mean(tmp,2);
%             tmp = repmat(tmp,size(allersp{cond_i},1),size(allersp{cond_i},2)); 
            %- subtract out baseline from all trials.
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
                           'titles', titles_ersp,...
                           'caxis',ERSP_CAXIS,...
                           'freqscale',ERSP_FREQSCALE,...
                           'ersplim',[],...
                           'effect','main',...
                           'threshold',ERSP_ALPHA,...
                           'maskdata','off','averagemode','ave');
%         std_plottf(alltimes, allfreqs, allersp,...
%                            'datatype', 'ersp', ...
%                            'groupstats', stat_group,...
%                            'condstats', stat_cond,...
%                            'interstats', stat_inter,...
%                            'plotmode', 'normal',...
%                            'unitx',XAXIS_LABEL,...
%                            'titles', titles_ersp,...
%                            'events', {},... %ersp_event_names,...
%                            'unitcolor', COLORAXIS_LABEL,...
%                            'chanlocs', ERSP_CHANLOCS,... %ALLEEG(1).chanlocs.labels,...
%                            'caxis',ERSP_CAXIS,...
%                            'freqscale',ERSP_FREQSCALE,...
%                            'ersplim',[],...
%                            'threshold',ERSP_ALPHA,'effect','main',...
%                            'maskdata','off','averagemode','rms');
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[100 100 1620 300]);
%         grandAvgWarpTo = [0,262,706,981,1414];
        for i = [2,4,5,6]
%             axes(fig_i.Children(2));
            tmp = fig_i.Children(i);
            set(tmp,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(tmp,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(tmp,sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
            xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
            set(tmp,'XTick',grandAvgWarpTo)
            xx = num2cell(grandAvgWarpTo);
            xx = cellfun(@num2str,xx,'UniformOutput',false);
%             set(tmp,'XTickLabel',xx)
            set(tmp,'XTickLabel',{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'});
            xlabel(tmp,'Time (ms)','Fontsize',16,'fontweight','bold');
            xline(tmp,grandAvgWarpTo(1),'k--');
            xline(tmp,grandAvgWarpTo(2),'k--');
            xline(tmp,grandAvgWarpTo(3),'k--');
            xline(tmp,grandAvgWarpTo(4),'k--');
            xline(tmp,grandAvgWarpTo(5),'k--');
        end
        fig_i.Name = sprintf('Cluster %i) ERSP averages across Conditions',cluster_i);
        saveas(fig_i,[save_dir filesep sprintf('custom2_Cond_ersp_%i_%s.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('custom2_Cond_ersp_%i_%s.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        
        %{
        %## temporary save
        tmp = fig_i.FileName;
        tmp = strsplit(tmp,'.');
        tmp{2} = 'jpg';
        tmp = strjoin(tmp,'.');
        saveas(fig_i,tmp);
        %}
        
    end
end
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