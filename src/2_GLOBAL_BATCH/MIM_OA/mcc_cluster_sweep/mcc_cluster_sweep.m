function [error_code] = mcc_cluster_sweep(study_fFull,cluster_load_dir,save_dir,varargin)
%MIM_FEHEADMODEL_DIPFIT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/run_mim_mcc_dipfit_exe.sh
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu

%## TIME
startj = tic;
%## DEFINE DEFAULTS
%- Hard Defines
% POSS_CLUSTER_CHARS = {'R_SenMA','L_SenMA','R_PPA','L_PPA','Cing','R_SupMA','L_SupMA','R_Occ','L_Occ','Caudate'};
% % this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
%* compute measures for spectrum and ersp
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; %false;
%* ERSP PARAMS
STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',10000); 
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','pdf',...
    'logtrials','on',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','on',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
%* loaded data for the custom plots
FORCE_RECALC_ERSP = true;
LOAD_OPTION_PLOTDATA = 'mim_custom_ersp_plots_default';
%- define output
error_code = 0;
%- working directory containing the vol.mat & elec_aligned.mat
errorMsg = 'Value ''study_fFull'' must be PATH. ''study_fFull'' should point to a .study file.'; 
sf_validFcn = @(x) assert(exist(x,'file'),errorMsg);
%- EEG filepath
errorMsg = 'Value ''cluster_load_dir'' must be PATH. ''cluster_load_dir'' should point to the folder containing cluster_update.mat for the resulting kmeans/knearest clustering K.'; 
cld_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%- Output for source.mat filepath
errorMsg = 'Value ''save_dir'' must be PATH. ''save_dir'' should point to a folder that exists. save_dir is the directory where all plots and data will be stored.';
sd_validFcn = @(x) assert(logical(exist(x,'file')),errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'study_fFull',sf_validFcn);
addRequired(p,'cluster_load_dir',cld_validFcn);
addRequired(p,'save_dir',sd_validFcn);
%## OPTIONAL
%## PARAMETER
parse(p,study_fFull,cluster_load_dir,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% fid = fopen([working_dir filesep 'output.txt'],'w');
%- EEGLAB options for opening EEG data
try
    fprintf(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID') '\n']);
    fprintf(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE') '\n']);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([working_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([working_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', 1440);
    disp(pPool);
catch e
    fprintf('Parallel processing failed to start');
    fprintf(['error. identifier: %s\n',...
             'error. %s\n',...
             'error. on working_dir %s\n'],e.identifier,e.message,working_dir);
end
%% LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~ispc
    tmp = strsplit(study_fFull,'/');
else
    tmp = strsplit(study_fFull,'\');
end
study_fPath = strjoin(tmp(1:end-1),filesep);
study_fName = strjoin(tmp(end),filesep);
if ~exist([study_fPath fiesep study_fName],'file')
    msg = 'error. study file does not exist.';
    errID = 'mcc_cluster_sweeps:ImproperInputs';
    baseException = MException(errID,msg);
    throw(baseException);
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',study_fName,'filepath',study_fPath);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',study_fName,'filepath',study_fPath);
    end
end
toc(startj)
%% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
%     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
averaged_warpto_events = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
% averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
ersp_ntimes = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
disp(averaged_warpto_events);
%## ersp plot per cluster per condition
%- params (06/10/2023) JS, removing timerange loading, consolidating single
%trials stats variable
STUDY = pop_statparams(STUDY,'condstats',STAT_PARAMS.condstats,...
        'groupstats',STAT_PARAMS.groupstats,...
        'method',STAT_PARAMS.methd,...
        'singletrials',STAT_PARAMS.singletrials,'mode',STAT_PARAMS.mode,'fieldtripalpha',STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',STAT_PARAMS.fieldtripmethod,'fieldtripmcorrect',STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
tmp = strsplit(ALLEEG(1).filename,'.');
spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
if ~exist(spec_f,'file') || ~exist(topo_f,'file')
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
                    {'specmode',SPEC_PARAMS.specmode,'freqfac',SPEC_PARAMS.freqfac,...
                    'freqrange',SPEC_PARAMS.freqrange,'logtrials',SPEC_PARAMS.logtrials});
    end
end
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
icatimf_f = [ALLEEG(1).filepath FILESEP sprintf('%s.icatimef',ALLEEG(1).subject)];
if ~exist(icatimf_f,'file') || ~FORCE_RECALC_ERSP
    %- load .icatimef load-in parameters
    tmp = load(icatimf_f,'-mat');
    parameters = tmp.parameters;
%     warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
    parameters{find(strcmp(parameters,'baseline'))+1} = [averaged_warpto_events(1),averaged_warpto_events(end)];
    parameters = [parameters, {'trialbase'}, {'off'}];
    cellArray = parameters(1,2:2:length(parameters));
    fields = parameters(1,1:2:length(parameters)-1);
    parameters = cell2struct(cellArray,fields,2);
    % NOTE:
    % (07/18/2023) CL, I don't think it works for trial 'on' condition somehow.
    % (07/18/2023) JS, Responding to 07/18/2023 CL comment: there is a switch
    % inside the ersp baseline function newtimfbaseln that controls whether
    % trialbase works or not. see docs/ERSP_baselining_approaches.docx. so
    % a later step might imploy this freature...
else
    disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
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
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',ersp_ntimes,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                    'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',ersp_ntimes,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
    %- load .icatimef load-in parameters
    tmp = load(icatimf_f,'-mat');
    parameters = tmp.parameters;
%     warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
    parameters{find(strcmp(parameters,'baseline'))+1} = [averaged_warpto_events(1),averaged_warpto_events(end)];
    parameters = [parameters, {'trialbase'}, {'off'}];
    cellArray = parameters(1,2:2:length(parameters));
    fields = parameters(1,1:2:length(parameters)-1);
    parameters = cell2struct(cellArray,fields,2);
end
%% DESIGN CHECK
%## Check to see if STUDY designs have been attached to the STUDY
try
    disp(STUDY.design(1).variable.value(:))
catch e
    msg = 'error. STUDY.design is not filled. set the STUDY designs you''d like to test using std_makedesign.';
    errID = 'mcc_cluster_sweeps:ImproperInputs';
    baseException = MException(errID,msg);
    throw(baseException);
end
%% BIG LOOPS TO CALCULATE EACH CLUSTER'S ERSP
if ~exist([save_dir filesep 'plots'],'dir')
    mkdir([save_dir filesep 'plots']);
end
%## load chang's algorithmic clustering
try
    %* load cluster information
    cluster_update = par_load(cluster_load_dir,sprintf('cluster_update_%i.mat',cl_i));
catch e
    errID = 'mcc_cluster_sweeps:ImproperInputs';
    msg = sprintf(2,'%s',getReport(e));
    baseException = MException(errID,msg);
    throw(baseException);
end
TMP_STUDY.cluster = cluster_update;
ERSP_PARAMS.timerange = [averaged_warpto_events(1),averaged_warpto_events(end)];
%- get inds
[~,main_cl_inds,~,] = eeglab_get_cluster_comps(TMP_STUDY);
%## Skip computation if it's complete, otherwise recompute.
if ~exist([save_dir filesep 'spec_ersp_data'],'dir')
    mkdir([save_dir filesep 'spec_ersp_data'])
end
try
    if max(TMP_STUDY.etc.mim_gen_ersp_data.cluster_n) == length(TMP_STUDY.cluster)
        fprintf('ERSPs & Spectrums already calculated for this clustering K. Skipping Computation..\n')
    else
        [TMP_STUDY] = mim_gen_ersp_data(TMP_STUDY,ALLEEG,averaged_warpto_events,...
            parameters,[save_dir filesep 'spec_ersp_data'],...
            'DO_SPEC_CALC',true,...
            'DO_ERSP_CALC',true,...
            'ERSP_PARAMS',ERSP_PARAMS,...
            'SPEC_PARAMS',SPEC_PARAMS,...
            'STAT_PARAMS',STAT_PARAMS);
    end
catch e
    fprintf(2,'error. %s\n',getReport(e));
    [TMP_STUDY] = mim_gen_ersp_data(TMP_STUDY,ALLEEG,averaged_warpto_events,...
            parameters,[save_dir filesep 'spec_ersp_data']);
end
[TMP_STUDY,~] = parfunc_save_study(TMP_STUDY,ALLEEG,...
                                    TMP_STUDY.filename,TMP_STUDY.filepath,...
                                    'RESAVE_DATASETS','off');
%- clusters to plot
CLUSTER_PICKS = main_cl_inds; %valid_clusters
%##
for des_i = 1:length(TMP_STUDY.design)
    cond_test = TMP_STUDY.design(des_i).variable(1).value;
    fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
    TMP_STUDY.currentdesign = des_i;
    fprintf('Current design: %i\n',TMP_STUDY.currentdesign);
%         parfor (j = 1:length(CLUSTER_PICKS),POOL_SIZE)
    for j = 1:length(CLUSTER_PICKS)
        cluster_i = CLUSTER_PICKS(j);
        cluster_load_ind = TMP_STUDY.etc.mim_gen_ersp_data.cluster_n(des_i,cluster_i);
        switch LOAD_OPTION_PLOTDATA
            case 'mim_custom_ersp_plots_default'
                %- defaults
                allersp = {};
                alltimes = [];
                allfreqs = [];
                pcond = [];
                pgroup = [];
                pinter = [];
            case 'load_opt_1'
                %- LOAD OPTION 1
                fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_fpaths,'/');
                fpath = strjoin(fname(1:end-1),'/');
                fname = fname{end};
                ersp_data = par_load(fpath,fname);
                allersp = ersp_data.allerspdata;
                alltimes = ersp_data.alltimes;
                allfreqs = ersp_data.allfreqs;
                pcond = ersp_data.pcond;
                pgroup = ersp_data.pgroup;
                pinter = ersp_data.pinter;
            case 'load_opt_2'
                %- LOAD OPTION 2
                fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_norm_fpaths,'/');
                fpath = strjoin(fname(1:end-1),'/');
                fname = fname{end};
                ersp_data_norm = par_load(fpath,fname);
                allersp = ersp_data_norm.allerspdata;
                alltimes = ersp_data_norm.alltimes;
                allfreqs = ersp_data_norm.allfreqs;
                pcond = ersp_data_norm.pcond;
                pgroup = ersp_data_norm.pgroup;
                pinter = ersp_data_norm.pgroup;
            case 'load_opt_3'
                %- LOAD OPTION 3
                fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).ersp_normcb_fpaths,'/');
                fpath = strjoin(fname(1:end-1),'/');
                fname = fname{end};
                ersp_data_normcb = par_load(fpath,fname);
                allersp = ersp_data_normcb.allerspdata;
                alltimes = ersp_data_normcb.alltimes;
                allfreqs = ersp_data_normcb.allfreqs;
                pcond = ersp_data_normcb.pcond;
                pgroup = ersp_data_normcb.pgroup;
                pinter = ersp_data_normcb.pinter;
        end
        %## RUN PLOTTING
        fprintf('Plotting Cluster %i for design %i\n',cluster_i,des_i);
        mim_custom_ersp_plots(TMP_STUDY,averaged_warpto_events,cluster_i,cluster_load_ind,des_i,save_dir,...
            'CLUSTER_CLIM_MATCH',CLUSTER_CLIM_MATCH,...
            'ALLERSP',allersp,'ALLTIMES',alltimes,'ALLFREQS',allfreqs,...
            'PCOND',pcond,'PGROUP',pgroup,'PINTER',pinter)
    end
end
end
%% (NOTES)
