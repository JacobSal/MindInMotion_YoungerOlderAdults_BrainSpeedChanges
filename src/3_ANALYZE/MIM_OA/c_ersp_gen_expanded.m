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
%     SLURM_POOL_SIZE = 2;
%     pp = parcluster('local');
%     pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
end
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
do_multivariate_data = 1;
evaluate_method = 'min_rv';
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result?
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% (08/03/2023) JS, changing fieldtripnaccu to 10,000; changing
% fieldtripmcorrect to cluster; method to perm; (these are the parameters
% set inside CL's PlotAndSaveERSP_CL_V3.m...
% pipeline although this doesn't align with her YA manuscript methods?
% (08/06/2023) JS, changing fieldtripnaccu to 2000 again and mcorrect to fdr...
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
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
%{
[~,~] = parfunc_save_study(STUDY,ALLEEG,...
                            STUDY.filename,load_dir,...
                            'RESAVE_DATASETS','on');
%}
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    
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
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)];
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
        TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
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
    disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
    parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/3))
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
            timewarpms_param = averaged_warpto_events;
         else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                    'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                    'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
end
%- load .icatimef load-in parameters
tmp = load(icatimf_f,'-mat');
parameters = tmp.parameters;
%     warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
parameters{find(strcmp(parameters,'baseline'))+1} = [averaged_warpto_events(1),averaged_warpto_events(end)];
parameters = [parameters, {'trialbase'}, {'off'}];

%-
% myErspParams_trialbasefull = parameters;
% myErspParams_trialbasefull{find(strcmp(parameters,'trialbase'))+1} = 'full';
% % myErspParams_trialbasefull{length(myErspParams_trialbasefull)+1} = 'basenorm';
% % myErspParams_trialbasefull{length(myErspParams_trialbasefull)+2} = 'on';

cellArray = parameters(1,2:2:length(parameters));
fields = parameters(1,1:2:length(parameters)-1);
parameters = cell2struct(cellArray,fields,2);

% cellArray= {myErspParams_trialbasefull{1,2:2:length(myErspParams_trialbasefull)}};
% fields = {myErspParams_trialbasefull{1,1:2:length(myErspParams_trialbasefull)-1}};
% myErspParams_trialbasefull = cell2struct(cellArray, fields,2);
% NOTE:
% (07/18/2023) CL, I don't think it works for trial 'on' condition somehow.
% (07/18/2023) JS, Responding to 07/18/2023 CL comment: there is a switch
% inside the ersp baseline function newtimfbaseln that controls whether
% trialbase works or not. see docs/ERSP_baselining_approaches.docx. so
% a later step might imploy this freature...
%% MAKE DESIGNS
%## NOTE (07/18/2023) JS, scripts/functions adapted from CL, AS, NJ scripts:
% PlotAndSaveERSP_CL_V3.m, Stats_Plotting_ERSPs_local_allCondBaseline_V3.m,
% Figures_Plotting_ERSPs_local_allCondBaseline_V3.m, plotERSP_parfor.m
STUDY.cache = [];
% tmp_group_orig = cell(length(ALLEEG),1);
% tmp_group_unif = cell(length(ALLEEG),1);
% for subj_i = 1:length(ALLEEG)
%     tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
%     tmp_group_unif{subj_i} = 'Older Adults';
% end
%## combine groups
% for subj_i = 1:length(ALLEEG)
%     ALLEEG(subj_i).group = tmp_group_unif{subj_i};
%     STUDY.datasetinfo(subj_i).group = tmp_group_unif{subj_i};
% end
[STUDY] = std_makedesign(STUDY,ALLEEG,1,'subjselect',{ALLEEG.subject},...
    'variable1','cond','values1',{'flat','low','med','high'},...
    'variable2','group','values2',{});
[STUDY] = std_makedesign(STUDY,ALLEEG,2,'subjselect',{ALLEEG.subject},...
    'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
    'variable2','group','values2',{});
% [STUDY] = std_makedesign(STUDY, ALLEEG, 3,...
%             'subjselect', {ALLEEG.subject},...
%             'variable2','group',...
%             'values2',unique({ALLEEG.group}),...
%             'variable1','',...
%             'values1',{});
% [STUDY] = std_makedesign(STUDY,ALLEEG,3,'subjselect',{ALLEEG.subject},...
%     'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
%     'VARIABLE2','group','values',unique({ALLEEG.group}));
%% SETUP FOR ERSP GENERATION FOR A PARTICULAR CLUSTER
% (07/27/2023) tried 21; need to try 19 (silohette) & 14
% (davies-bouldin)

%## Gather Load-In Parameters & Make Edits
% NOTE:
% (07/18/2023) CL, I don't think it works for trial 'on' condition somehow.
% (07/18/2023) JS, Responding to 07/18/2023 CL comment: there is a switch
% inside the ersp baseline function newtimfbaseln that controls whether
% trialbase works or not. see docs/ERSP_baselining_approaches.docx
%-
% CLUSTER_K_PICKS = [15,19,21];
CLUSTER_K_PICKS = [14,19,21];
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
%%
for k_i = 1:length(CLUSTER_K_PICKS)
    %## Loop Params
    cl_i = CLUSTER_K_PICKS(k_i);
    fprintf('Making data for K=%i...\n',cl_i);
    TMP_STUDY = STUDY;
    spec_data_dir = [save_dir filesep sprintf('k%i',cl_i) filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        mkdir(spec_data_dir)
    end
    %## load chang's algorithmic clustering
    %* load cluster information
    cluster_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster',...
        filesep clustering_method filesep num2str(cl_i)];
    cluster_update = par_load(cluster_dir,sprintf('cluster_update_%i.mat',cl_i));
    TMP_STUDY.cluster = cluster_update;
    %- get inds
    [~,main_cl_inds,~,valid_clusters] = eeglab_get_cluster_comps(TMP_STUDY);
    CLUSTER_PICKS = 1:length(TMP_STUDY.cluster);
    %- gather clustering information
    load_ind_cl = zeros(1,length(CLUSTER_PICKS));
    clust_ind_cl = zeros(1,length(CLUSTER_PICKS));
    for i = 1:length(CLUSTER_PICKS)
        load_ind_cl(i) = i;
        clust_ind_cl(i) = CLUSTER_PICKS(i);
    end
    %- stores
    spec_fpaths = cell(length(TMP_STUDY.design),length(CLUSTER_PICKS));
    spec_ss_fpaths = cell(length(TMP_STUDY.design),length(CLUSTER_PICKS));
    ersp_fpaths = cell(length(TMP_STUDY.design),length(CLUSTER_PICKS));
    ersp_norm_fpaths = cell(length(TMP_STUDY.design),length(CLUSTER_PICKS));
    ersp_normcb_fpaths = cell(length(TMP_STUDY.design),length(CLUSTER_PICKS));
    %- loop designs
    for des_i = 1:length(TMP_STUDY.design)
        store_1 = cell(1,length(CLUSTER_PICKS));
        store_2 = cell(1,length(CLUSTER_PICKS));
        store_3 = cell(1,length(CLUSTER_PICKS));
        store_4 = cell(1,length(CLUSTER_PICKS));
        store_5 = cell(1,length(CLUSTER_PICKS));
        TMP_STUDY.currentdesign = des_i;
        design_char = [];
        for i = 1:length(TMP_STUDY.design(des_i).variable)
            if i==1
                design_char = [TMP_STUDY.design(des_i).variable(1).value{:}];
            else
                design_char = [design_char '_' TMP_STUDY.design(des_i).variable(i).value{:}];
            end
        end
        %## (PARFOR) Compute Specs & ERSPs
        parfor (i = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
            cluster_i = CLUSTER_PICKS(i);
            %- 
            [~,spec_savef,spec_subjcorr_savef] = mim_gen_spec_data(TMP_STUDY,ALLEEG,...
                    averaged_warpto_events,cluster_i,design_char,spec_data_dir);
            store_1{i} = spec_savef;
            store_2{i} = spec_subjcorr_savef;
            %- 
            [~,ersp_savef,ersp_subbase_savef,ersp_subbase_combase_savef] = mim_gen_ersp_data(TMP_STUDY,ALLEEG,averaged_warpto_events,...
                    parameters,cluster_i,design_char,spec_data_dir,...
                    'STAT_PARAMS',ERSP_STAT_PARAMS);
            store_3{i} = ersp_savef;
            store_4{i} = ersp_subbase_savef;
            store_5{i} = ersp_subbase_combase_savef;
        end
        spec_fpaths(des_i,:) = store_1;
        spec_ss_fpaths(des_i,:) = store_2;
        ersp_fpaths(des_i,:)  = store_3;
        ersp_norm_fpaths(des_i,:)  = store_4;
        ersp_normcb_fpaths(des_i,:)  = store_5;
    end
    %## Store Paths in Struct
    gen_data_struct = struct('load_ind_cl',load_ind_cl,...
        'clusst_ind_cl',clust_ind_cl,...
        'ersp_fpaths',ersp_fpaths,...
        'spec_fpaths',spec_fpaths,...
        'spec_ss_fpaths',spec_ss_fpaths,...
        'ersp_norm_fpaths',ersp_norm_fpaths,...
        'ersp_normcb_fpaths',ersp_normcb_fpaths);
    TMP_STUDY.etc.mim_gen_ersp_data = gen_data_struct;
    par_save(gen_data_struct,spec_data_dir,'spec_data_struct.mat');
    [~,~] = parfunc_save_study(TMP_STUDY,ALLEEG,...
                                    TMP_STUDY.filename,spec_data_dir,...
                                    'RESAVE_DATASETS','off');

end
%% (STEP 1) GENERATE ERSP & SPEC DATA FOR-EACH DESIGN & CLUSTER
% delete(pPool);
% % parfor (i = 1:length(CLUSTER_K_PICKS),length(CLUSTER_K_PICKS)*2)
% pPool = parpool(pp, ceil(length(ALLEEG)/2), 'IdleTimeout', 1440);
% for i = 1:length(CLUSTER_K_PICKS)
%     %## Loop Params
%     cl_i = CLUSTER_K_PICKS(i);
%     fprintf('Making data for K=%i...\n',cl_i);
%     TMP_STUDY = STUDY;
%     spec_data_dir = [save_dir filesep sprintf('k%i',cl_i) filesep 'spec_data'];
%     if ~exist(spec_data_dir,'dir')
%         mkdir(spec_data_dir)
%     end 
%     %## load chang's algorithmic clustering
%     %* load cluster information
%     cluster_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster',...
%         filesep clustering_method filesep num2str(cl_i)];
%     cluster_update = par_load(cluster_dir,sprintf('cluster_update_%i.mat',cl_i));
%     TMP_STUDY.cluster = cluster_update;
%     %- get inds
%     [~,main_cl_inds,~,valid_clusters] = eeglab_get_cluster_comps(TMP_STUDY);
%     %- clusters to plot
% %     CLUSTER_PICKS = main_cl_inds; %valid_clusters
%     %## Compute ERSP/Spec data and store in -mat
%     [TMP_STUDY,ersp_data_struct] = mim_gen_ersp_data(TMP_STUDY,ALLEEG,...
%         averaged_warpto_events,parameters,spec_data_dir,...
%         'DO_SPEC_CALC',true,...
%         'DO_ERSP_CALC',true,...
%         'STAT_PARAMS',ERSP_STAT_PARAMS,...
%         'SPEC_PARAMS',SPEC_PARAMS,...
%         'ERSP_PARAMS',ERSP_PARAMS);
%     par_save(ersp_data_struct,spec_data_dir,'spec_data_struct.mat');
%     [~,~] = parfunc_save_study(TMP_STUDY,ALLEEG,...
%                                     TMP_STUDY.filename,spec_data_dir,...
%                                     'RESAVE_DATASETS','off');
% end
%% (STEP 2) PLOT
% delete(pPool);
% pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
%##
clear('ALLEEG');
for i = 1:length(CLUSTER_K_PICKS)
    %## Loop Params
    cl_i = CLUSTER_K_PICKS(i);
    plot_store_dir = [save_dir filesep sprintf('k%i',cl_i) filesep 'plots_out'];
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    spec_data_dir = [save_dir filesep sprintf('k%i',cl_i) filesep 'spec_data'];
    %## Load Study
    % (08/03/2023) JS, this can be optimized in the future by only loding
    % in the file strucutre and maintaining a STUDY file with needed info
    if ~exist([spec_data_dir filesep study_fName_1 '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [TMP_STUDY,~] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',spec_data_dir);
        else
            [TMP_STUDY,~] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',spec_data_dir);
        end
    end
    %- get inds
    [~,main_cl_inds,~,valid_clusters] = eeglab_get_cluster_comps(TMP_STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters
    %## Loop Through Designs
    for des_i = 1:length(TMP_STUDY.design)
        cond_test = TMP_STUDY.design(des_i).variable(1).value;
        fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
        TMP_STUDY.currentdesign = des_i;
        fprintf('Current design: %i\n',TMP_STUDY.currentdesign);
        fprintf('Statistics Parameters:\n');
        disp(STUDY.etc.statistics)
        fprintf('Statistics Fieldtrip Parameters:\n');
        disp(STUDY.etc.statistics.fieldtrip)
        fprintf('Statistics EEGLAB Parameters:\n');
        disp(STUDY.etc.statistics.eeglab)
        fprintf('ERSP Parameters:\n');
        disp(STUDY.etc.erspparams)
        
        parfor (j = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
%         for j = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(j);
            cluster_load_ind = cluster_i;
%             cluster_load_ind = TMP_STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).cluster_n(1,cluster_i);
            %- defaults
            allersp = {};
            alltimes = [];
            allfreqs = [];
            pcond = {};
            pgroup = {};
            pinter = {};
            %- LOAD OPTION 1
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data = par_load(fpath,fname);
    %         allersp = ersp_data.allerspdata;
    %         alltimes = ersp_data.alltimes;
    %         allfreqs = ersp_data.allfreqs;
            %- LOAD OPTION 2
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_norm_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data_norm = par_load(fpath,fname);
    %         allersp = ersp_data_norm.allerspdata;
    %         alltimes = ersp_data_norm.alltimes;
    %         allfreqs = ersp_data_norm.allfreqs;
    %         pcond = ersp_data_norm.pcond;
            %- LOAD OPTION 3
    %         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_normcb_fpaths,'/');
    %         fpath = strjoin(fname(1:end-1),'/');
    %         fname = fname{end};
    %         ersp_data_normcb = par_load(fpath,fname);
    %         allersp = ersp_data_normcb.allerspdata;
    %         alltimes = ersp_data_normcb.alltimes;
    %         allfreqs = ersp_data_normcb.allfreqs;
            %- LOAD OPTION 4
    %         [~, allersp, alltimes, allfreqs, ~, ~] = std_readdata(STUDY,ALLEEG,...
    %             'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
    %             'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
    %             'design',des_i);
    %         %* get stats
    %         [pcond_ersp,pgroup_ersp,pinter_ersp,~,~,~] = std_stat(allersp,...
    %                     'condstats', ERSP_CONDSTATS,...
    %                     'groupstats',ERSP_GROUPSTATS,...
    %                     'method',ERSP_STAT_METHOD,...
    %                     'naccu',ERSP_NACCU,...
    %                     'alpha',ERSP_ALPHA,...
    %                     'mcorrect',ERSP_MCORRECT,'mode',STAT_MODE,'fieldtripalpha',ERSP_ALPHA,...
    %                     'fieldtripmethod',FIELDTRIP_METHOD,'fieldtripmcorrect',ERSP_MCORRECT,'fieldtripnaccu',ERSP_NACCU);
            %## RUN PLOTTING
            fprintf('Plotting Cluster %i for design %i\n',cluster_i,des_i);
            mim_custom_ersp_plots(TMP_STUDY,cond_test,averaged_warpto_events,...
                cluster_i,cluster_load_ind,des_i,plot_store_dir,...
                'CLUSTER_CLIM_MATCH',CLUSTER_CLIM_MATCH,...
                'ALLERSP',allersp,...
                'ALLTIMES',alltimes,...
                'ALLFREQS',allfreqs,...
                'PCOND',pcond,...
                'PGROUP',pgroup,...
                'PINTER',pinter)
        end
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}