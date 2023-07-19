%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_c_ersp_gen_expanded.sh

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
RECOMPUTE_ERSP = false;
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
% dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
dt = '07152023_MIM_OAN79_subset_prep_verified_gait';
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
    pick_cluster = 21;
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
    cluster_update = par_load(outputdir,sprintf('cluster_update_%i.mat',pick_cluster));
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
b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)];
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
%% (STEP 1) GENERATE VARIOUS ERSP & SPEC W/ DIFFERENT BASELINES & NORMALIZATIONS
% NOTE (07/18/2023) JS, scripts/functions adapted from CL scripts:
% PlotAndSaveERSP_CL_V3.m, Stats_Plotting_ERSPs_local_allCondBaseline_V3.m,
% Figures_Plotting_ERSPs_local_allCondBaseline_V3.m, plotERSP_parfor.m
STUDY.cache = [];
[STUDY] = std_makedesign(STUDY,ALLEEG,1,'subjselect',{ALLEEG.subject},...
    'variable1','cond','values1',{'flat','low','med','high'});
[STUDY] = std_makedesign(STUDY,ALLEEG,2,'subjselect',{ALLEEG.subject},...
    'variable1','cond','values1',{'0p25','0p5','0p75','1p0'});
%- load .icatimef load-in parameters
icatimf_fname = sprintf('%s.icatimef',ALLEEG(1).subject);
tmp = load([ALLEEG(1).filepath filesep icatimf_fname],'-mat');
parameters = tmp.parameters;
warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
baseline_range = [warpingvalues(1) warpingvalues(end)];
% b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)];
ersp_load_params = parameters;
ersp_load_params{find(strcmp(parameters,'baseline'))+1} = baseline_range;
ersp_load_params{length(parameters)+1} = 'trialbase';
% NOTE:
% (07/18/2023) CL, I don't think it works for trial 'on' condition somehow.
% (07/18/2023) JS, Responding to 07/18/2023 CL comment: there is a switch
% inside the ersp baseline function newtimfbaseln that controls whether
% trialbase works or not. see docs/ERSP_baselining_approaches.docx
ersp_load_params{length(parameters)+2} = 'off';
cellArray = {ersp_load_params{1,2:2:length(ersp_load_params)}};
fields = {ersp_load_params{1,1:2:length(ersp_load_params)-1}};
ersp_load_params = cell2struct(cellArray,fields,2);
%- generate various baselined ersp 
mim_gen_ersp_data(STUDY,ALLEEG,grandAvgWarpTo,...
    ersp_load_params,save_dir);
%% (STEP 2) GENERATE STATS DATA

%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}