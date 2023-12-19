%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_ersps_timewarp.sh

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
%- define the directory to the src folderd
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
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa');
fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams';
dt = 'spca_analysis';
%- study group and saving
SAVE_ALLEEG = false;
%- epoching params
DO_SLIDING_WINDOW = false;
%* sliding window
WINDOW_LENGTH = 6; % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs
%* gait
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-1,4.25];
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
if DO_SLIDING_WINDOW
    SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED';
    TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
else
    SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED';
    TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
end
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
study_fName_1 = 'all_comps_study';
study_fName_2 = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = zeros(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        try
            fNames{cnt} = tmp.name;
            %- Chanlocs fPaths
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            dipfit_norm_fPaths(cnt) = 1;
        catch e
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('%s\n',getReport(e))
        end
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
%         disp(cnt)
%         disp(SUBJ_PICS{group_i}{subj_i})
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
inds = logical(dipfit_norm_fPaths);
fPaths = fPaths(inds);
fNames = fNames(inds);
subjectNames = subjectNames(inds);

%% ===================================================================== %%

%## GENERATE EPOCH MAIN FUNC
tmp = cell(1,length(fPaths));
rmv_subj = zeros(1,length(fPaths));
%## PARFOR LOOP
parfor (subj_i = 1:length(fPaths),floor(length(fPaths)/3))
% for subj_i = 1:length(subjectNames)
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    fprintf('Running subject %s\n',EEG.subject)
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    %## PARSE TRIALS
    epoched_fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED];
    fPath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[TRIAL_TYPES{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    try
        %## EPOCH
        [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,DO_SLIDING_WINDOW,...
            'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,...
            'STD_TIMEWARP',STD_TIMEWARP,...
            'COND_CHARS',TRIAL_TYPES);
        %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
        for i = 1:length(ALLEEG)
            if isfield(ALLEEG(i).event,'trialName')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'trialName');
            end
            if isfield(ALLEEG(i).event,'channel')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'channel');
            end
            if isfield(ALLEEG(i).event,'code')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'code');
            end
            if isfield(ALLEEG(i).event,'bvtime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvtime');
            end
            if isfield(ALLEEG(i).event,'bvmknum')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvmknum');
            end
            if isfield(ALLEEG(i).event,'datetime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'datetime');
            end
        end
        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_ALLEEG
            for i = 1:length(ALLEEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,ALLEEG(i).condition)];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(ALLEEG(i),'savemode','twofiles',...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition);
            end
        end
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## timewarp for across condition
        if ~DO_SLIDING_WINDOW
            timewarp = make_timewarp(ALLEEG,TIMEWARP_EVENTS,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select(ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
        end
        %## STRUCT EDITS
        ALLEEG.urevent = []; % might be needed
        ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
        
        %- checks
        ALLEEG = eeg_checkset(ALLEEG,'eventconsistency');
        ALLEEG = eeg_checkset(ALLEEG);
        ALLEEG = eeg_checkamica(ALLEEG);
        
        %- save
        [ALLEEG] = pop_saveset(ALLEEG,'savemode','twofiles',...
                'filename',fName,...
                'filepath',fPath,...
                'version','6');
        tmp{subj_i} = ALLEEG;
    catch e
        rmv_subj(subj_i) = 1;
        EEG.timewarp = struct([]);
        EEG.urevent = [];
        tmp{subj_i} = []; %EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
    end
end
%% ===================================================================== %%
%## SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s',subjectNames{cellfun(@isempty,tmp)});
tmp = tmp(~cellfun(@isempty,tmp));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i})';
    disp(size(fields(tmp{subj_i})'));
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false); 
    out = [out{:}];
    addFs = fsPrev(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
        end
    end 
%     tmp{subj_i} = EEG;
    tmp{subj_i} = orderfields(EEG);
end
%- CONCATENATE tmp
tmp = cellfun(@(x) [[]; x], tmp);
%##
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName_2,...
                                'filename',study_fName_2,...
                                'filepath',save_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
%% ===================================================================== %%
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_2 '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_2 '.study'],'filepath',save_dir);
% end

%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
FORCE_RECALC_ERSP = false;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false;
DO_SUBJ_PLOTS = true;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
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
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
STUDY.etc.averaged_warpto_events = averaged_warpto_events;
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
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/3))
% for subj_i = 1:length(ALLEEG)
    EEG = ALLEEG(subj_i);
    TMP_STUDY = STUDY;
    EEG = eeg_checkset(EEG,'loaddata');
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
        [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','off',...
                'recompute','on','ersp','on','itc','off',...
                'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
                'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
                'baseline',nan(),'timewarp',timewarp_param,...
                'timewarpms',timewarpms_param}); %ERSP
    end
end
%%
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
