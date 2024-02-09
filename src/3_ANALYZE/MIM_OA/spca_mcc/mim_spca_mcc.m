function [err] = mim_spca_mcc(ica_dir,study_dir,jf_fpath,varargin)
%MIM_SPCA_MCC Summary of this function goes here
%   Detailed explanation goes here
%MIM_FEHEADMODEL_DIPFIT Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%       REQUIRED:
%       ica_dir, CHAR
%           directory where ICA'd EEG data is located for each subject
%           you'd like to run.
%           E.G.,
%               my_ica_folder/
%               ├─ subj_1/
%               │  ├─ clean/
%               ├─ subj_2/
%               │  ├─ clean/
%               ├─ subj_n/
%               │  ├─ clean/
%       study_dir, CHAR
%           directory where the program will store information about the
%           subjects in 'ica_dir'. Please see 'SPCA_PARAMS.force_recalc_ersp' to
%           write-over icatime.f files, or 'EPOCH_PARAMS.do_recalc_epoch' to
%           write-over epoch'd .set files.
%
%       jf_fpath, CHAR
%       JSON Format is as follows:
%             {
%                 "SPCA_PARAMS": {
%                     "iclabel_eye_cutoff": 0.75,
%                     "force_recalc_ersp": true,
%                     "do_timewarp": true,
%                     "do_baseline_corr": false
%                 },
%                 "EPOCH_PARAMS": null,
%                 "ERSP_STAT_PARAMS": null,
%                 "SPEC_PARAMS": null,
%                 "ERSP_PARAMS": null,
%                 "SUBJ_CHARS": null
%             }
% 
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
%% (PARAMETERS) ======================================================== %%
%## hard define
err = 0;
study_fname_rest = 'rest_epoch_study';
study_fname_gait = 'gait_epoch_study';
%## SPCA_PARAMS
SPCA_PARAMS = struct('iclabel_eye_cutoff',0.75,...
    'force_recalc_ersp',true,...
    'do_timewarp',true,...
    'do_baseline_corr',false);
%## EPOCH PARAMS
EPOCH_PARAMS = struct('percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','spca_epoched',...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char','rest',...
    'do_recalc_epoch',true);
%## SPECTRUM CALCULATION PARAMS & STATS
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
%## SUBJECT CHARACTERS LOAD-IN FROM FUNCTION
[SUBJ_CHARS,~,~,~,~,~,~] = mim_dataset_information('yaoa');
SUBJ_CHARS = [SUBJ_CHARS{:}];
%% LOAD JSON & ASSIGN CHECKS
%- connectivity statistics & surrogates params
if ~isempty(jf_fpath)
    fid = fopen(jf_fpath,'r');
    raw = fread(fid,inf); 
    str = char(raw');
    fclose(fid); 
    params = jsondecode(str);
    if ~isempty(params.SPCA_PARAMS)
        tmp_spca = params.SPCA_PARAMS;
        varargin = [varargin, 'SPCA_PARAMS', tmp_spca];
    end
    if ~isempty(params.EPOCH_PARAMS)
        tmp_epoch = params.EPOCH_PARAMS;
        varargin = [varargin, 'EPOCH_PARAMS', tmp_epoch];
    end
    if ~isempty(params.ERSP_STAT_PARAMS)
        tmp_ersp_stat = params.ERSP_STAT_PARAMS;
        varargin = [varargin, 'ERSP_STAT_PARAMS', tmp_ersp_stat];
    end
    if ~isempty(params.SPEC_PARAMS)
        tmp_spec = params.SPEC_PARAMS;
        varargin = [varargin, 'SPEC_PARAMS', tmp_spec];
    end
    if ~isempty(params.ERSP_PARAMS)
        tmp_ersp = params.ERSP_PARAMS;
        varargin = [varargin, 'ERSP_PARAMS', tmp_ersp];
    end
    if ~isempty(params.SUBJ_CHARS)
        tmp_subj_chars = params.SUBJ_CHARS;
        varargin = [varargin, 'SUBJ_CHARS', tmp_subj_chars];
    end
end
%## TIME
%## working directory containing subject ICA .set files
errorMsg = 'Value must be CHAR. Working directory containing all subjects to be epoched & sPCA''d.'; 
id_validFcn = @(x) assert(ischar(x)  && exist(x,'dir'),errorMsg);
errorMsg = 'Value must be CHAR. New directory to store epoched .set files and study information.'; 
sd_validFcn = @(x) assert(ischar(x),errorMsg);
errorMsg = 'Value must be CHAR. File path to a .json file containing parameters to be used.'; 
jf_validFcn = @(x) assert((ischar(x)  && exist(x,'file')) || isempty(x),errorMsg);
%%
fprintf('ica_dir: %s\n',ica_dir);
fprintf('study_dir: %s\n',study_dir);
fprintf('jf_fpath: %s\n',jf_fpath);
if exist(str,'var')
    fprintf('%s',str);
end
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ica_dir',id_validFcn);
addRequired(p,'study_dir',sd_validFcn);
addRequired(p,'jf_fpath',jf_validFcn);
%## OPTIONAL
addParameter(p,'SUBJ_CHARS',SUBJ_CHARS,@iscell);
addParameter(p,'SPCA_PARAMS',SPCA_PARAMS,@(x) validate_struct(x,SPCA_PARAMS));
addParameter(p,'EPOCH_PARAMS',EPOCH_PARAMS,@(x) validate_struct(x,EPOCH_PARAMS));
addParameter(p,'ERSP_STAT_PARAMS',ERSP_STAT_PARAMS,@(x) validate_struct(x,ERSP_STAT_PARAMS));
addParameter(p,'SPEC_PARAMS',SPEC_PARAMS,@(x) validate_struct(x,SPEC_PARAMS));
addParameter(p,'ERSP_PARAMS',ERSP_PARAMS,@(x) validate_struct(x,ERSP_PARAMS));
%## PARAMETERS
parse(p,ica_dir,study_dir,jf_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
SUBJ_CHARS = p.Results.SUBJ_CHARS;
SPCA_PARAMS = p.Results.SPCA_PARAMS;
EPOCH_PARAMS = p.Results.EPOCH_PARAMS;
ERSP_STAT_PARAMS = p.Results.ERSP_STAT_PARAMS;
SPEC_PARAMS = p.Results.SPEC_PARAMS;
ERSP_PARAMS = p.Results.ERSP_PARAMS;
if ~exist(study_dir,'dir')
    mkdir(study_dir);
end
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
%% (GRAB SUBJECT .SET & DIPFIT FILES) ==================================== %%
fnames        = cell(1,length(SUBJ_CHARS));
fpaths        = cell(1,length(SUBJ_CHARS));
chk           = zeros(1,length(SUBJ_CHARS));
for subj_i = 1:length(SUBJ_CHARS)
    %- ICA fPaths
    fpaths{subj_i} = [ica_dir filesep SUBJ_CHARS{subj_i} filesep 'clean'];
    tmp = dir([fpaths{subj_i} filesep '*.set']);
    try
        fnames{subj_i} = tmp.name;
        %- Prints
        fprintf('%s) .set fpath: %s\n',SUBJ_CHARS{subj_i}, fpaths{subj_i});
        chk(subj_i) = 1;
    catch e
        fprintf('Error on subject %s\n',SUBJ_CHARS{subj_i})
        fprintf('%s\n',getReport(e))
    end
end
%- remove subjects without a dipole fit
inds = logical(chk);
fpaths = fpaths(inds);
fnames = fnames(inds);
SUBJ_CHARS = SUBJ_CHARS(inds);

%% ===================================================================== %%
if ~ispc
    chk = exist([study_dir filesep [study_fname_gait '_UNIX.study']],'file') && ...
        exist([study_dir filesep [study_fname_gait '_UNIX.study']],'file');
else
    chk = exist([study_dir filesep [study_fname_gait '.study']],'file') && ...
        exist([study_dir filesep [study_fname_gait '.study']],'file');
end
if EPOCH_PARAMS.do_recalc_epoch || ~chk
    %## GENERATE EPOCH MAIN FUNC
    tmp = cell(1,length(fpaths));
    tmp_rest = cell(1,length(fpaths));
    rmv_subj = zeros(1,length(fpaths));
    %## PARFOR LOOP
    parfor (subj_i = 1:length(fpaths),SLURM_POOL_SIZE)
    % for subj_i = 1:length(subjectNames)
        %## LOAD EEG DATA
    %     EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        [~,EEG,~] = eeglab_loadICA(fnames{subj_i},fpaths{subj_i});
        fprintf('Running subject %s\n',EEG.subject)
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        EEG = iclabel(EEG);
    %         clss = EEG.etc.ic_classification.ICLabel.classes;
        clssts = EEG.etc.ic_classification.ICLabel.classifications;
        bad_eye_ics = find(clssts(:,3) > SPCA_PARAMS.iclabel_eye_cutoff);
        EEG = eeglab_pop_subcomp(EEG,bad_eye_ics,true);
    %     EEG = pop_subcomp(EEG,bad_eye_ics,0,0);
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %## PARSE TRIALS
        epoched_fPath = [study_dir filesep EEG.subject filesep EPOCH_PARAMS.path_ext];
        fPath = [epoched_fPath filesep [EPOCH_PARAMS.gait_trial_chars{:}]];
        fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[EPOCH_PARAMS.gait_trial_chars{:}]);
        if ~exist(fPath,'dir')
            mkdir(fPath)
        end
        %- parse
        try
            %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
            if isfield(EEG.event,'trialName')
                EEG.event = rmfield(EEG.event,'trialName');
            end
            if isfield(EEG.event,'channel')
                EEG.event = rmfield(EEG.event,'channel');
            end
            if isfield(EEG.event,'code')
                EEG.event = rmfield(EEG.event,'code');
            end
            if isfield(EEG.event,'bvtime')
                EEG.event = rmfield(EEG.event,'bvtime');
            end
            if isfield(EEG.event,'bvmknum')
                EEG.event = rmfield(EEG.event,'bvmknum');
            end
            if isfield(EEG.event,'datetime')
                EEG.event = rmfield(EEG.event,'datetime');
            end
            %## EPOCH
            [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,false,...
                'EPOCH_TIME_LIMITS',EPOCH_PARAMS.epoch_time_lims,... 
                'STD_TIMEWARP',EPOCH_PARAMS.tw_stdev,...
                'COND_CHARS',EPOCH_PARAMS.gait_trial_chars);

            %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
            fprintf('\nConcatenating datasets...\n');
            cond_files = struct('fPath',[],'fName',[]);
            ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
            ALLEEG.etc.cond_files = cond_files;
            %## timewarp for across condition
            if ~DO_SLIDING_WINDOW
                timewarp = make_timewarp(ALLEEG,EPOCH_PARAMS.tw_events,'baselineLatency',0, ...
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
                ALLEEG.etc.epoch.epoch_limits = EPOCH_PARAMS.epoch_time_lims;
            end
            %## STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_PARAMS.epoch_time_lims;
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
            %## RESTING STATE
            fPath = [epoched_fPath filesep EPOCH_PARAMS.rest_trial_char];
            fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,'rest');
            if ~exist(fPath,'dir')
                mkdir(fPath)
            end
            tmp_eeg = [];v
            [tmp_eeg,~] = mim_parse_trials(EEG,true,...
                'WINDOW_LENGTH',(EPOCH_PARAMS.epoch_time_lims(2)-EPOCH_PARAMS.epoch_time_lims(1)),...
                'PERCENT_OVERLAP',EPOCH_PARAMS.percent_overlap,...
                'COND_CHARS',{EPOCH_PARAMS.rest_trial_char});
            %- save
            [tmp_eeg] = pop_saveset(tmp_eeg,'savemode','twofiles',...
                    'filename',fName,...
                    'filepath',fPath,...
                    'version','6');
            tmp_rest{subj_i} = tmp_eeg;
        catch e
            rmv_subj(subj_i) = 1;
            EEG.timewarp = struct([]);
            EEG.urevent = [];
            tmp{subj_i} = []; %EEG;
            tmp_rest{subj_i} = [];
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n',...
                     'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
        end
    end
    fprintf('Bugged subjects:\n');
    fprintf('%s\n',SUBJ_CHARS{find(rmv_subj)});
    %% ===================================================================== %%
    %## SAVE BIG STUDY
    fprintf('==== Reformatting Study ====\n');
    %- remove bugged out subjects
    tmp_rest = tmp_rest(~cellfun(@isempty,tmp_rest));
    %## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
    fss = cell(1,length(tmp_rest));
    for subj_i = 1:length(tmp_rest)
        fss{subj_i} = fields(tmp_rest{subj_i})';
        disp(size(fields(tmp_rest{subj_i})'));
    end
    fss = unique([fss{:}]);
    fsPrev = fss;
    for subj_i = 1:length(tmp_rest)
        EEG = tmp_rest{subj_i};
        fs = fields(EEG);
        % delete fields not present in other structs.
        out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false); 
        out = [out{:}];
        addFs = fsPrev(~out);
        if any(~out)
            for j = 1:length(addFs)
                EEG.(addFs{j}) = [];
    %             fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
            end
        end 
    %     tmp_rest{subj_i} = EEG;
        tmp_rest{subj_i} = orderfields(EEG);
    end
    %- CONCATENATE tmp_rest
    tmp_rest = cellfun(@(x) [[]; x], tmp_rest);
    %##
    [STUDY, ALLEEG] = std_editset([],tmp_rest,...
                                    'updatedat','off',...
                                    'savedat','off',...
                                    'name',study_fname_rest,...
                                    'filename',study_fname_rest,...
                                    'filepath',study_dir);
    [STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
    [STUDY_REST,ALLEEG_REST] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY.filename,STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
    %% ===================================================================== %%
    %## SAVE BIG STUDY
    fprintf('==== Reformatting Study ====\n');
    %- remove bugged out subjects
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
    %             fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
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
                                    'name',study_fname_gait,...
                                    'filename',study_fname_gait,...
                                    'filepath',study_dir);
    [STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
    [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY.filename,STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
else
    %% ===================================================================== %%
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '_UNIX.study'],'filepath',study_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname_gait '.study'],'filepath',study_dir);
    end
    if ~ispc
        [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '_UNIX.study'],'filepath',study_dir);
    else
        [STUDY_REST,ALLEEG_REST] = pop_loadstudy('filename',[study_fname_rest '.study'],'filepath',study_dir);
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
SUBJ_CHARS = {ALLEEG.subject};
disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
for subj_i = 1:length(ALLEEG)
    %## 
    tt = tic;
    EEG = ALLEEG(subj_i);
    fprintf('Running Subject %s\n',EEG.subject);
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    EEG = iclabel(EEG);
%         clss = EEG.etc.ic_classification.ICLabel.classes;
    clssts = EEG.etc.ic_classification.ICLabel.classifications;
    bad_eye_ics = find(clssts(:,3) > ICLABEL_EYE_CUTOFF);
    EEG = eeglab_pop_subcomp(EEG,bad_eye_ics,true);
    EEG = pop_subcomp(EEG,bad_eye_ics,0,0);
    icatimf_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
    if ~exist(icatimf_f,'file') || SPCA_PARAMS.force_recalc_ersp % || any(strcmp(EEG.subject,FINISHED_ADULTS))
        TMP_STUDY = STUDY;
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
                    'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                    'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
                    'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
        fprintf('Done calculating timewarped ERSPs for %s\n',EEG.subject);
    else
        fprintf('Timewarped ERSPs already calculated for %s\n',EEG.subject);
    end
    
    %##
    EEG = ALLEEG_REST(subj_i);
    icatimf_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];
    if ~exist(icatimf_f,'file') || SPCA_PARAMS.force_recalc_ersp % || any(strcmp(EEG.subject,FINISHED_ADULTS))
        TMP_STUDY = STUDY_REST;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        EEG = eeglab_pop_subcomp(EEG,bad_eye_ics,true);
        %- overrride datasetinfo to trick std_precomp to run.
        TMP_STUDY.datasetinfo = STUDY_REST.datasetinfo(subj_i);
        TMP_STUDY.datasetinfo(1).index = 1;
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
                    'ntimesout',TIMEWARP_NTIMES,...
                    'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
                    'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES}); %ERSP
        end
        fprintf('Done calculating timewarped ERSPs for %s\n',EEG.subject);
    else
        fprintf('Timewarped ERSPs already calculated for %s\n',EEG.subject);
    end
    
    try
        %## (LOAD REST) ================================================ %%
        epoched_fPath = strsplit(EEG.filepath,filesep);
        epoched_fPath = strjoin(epoched_fPath(1:end-1),filesep);
        fpath = [epoched_fPath filesep 'rest'];
%         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,'rest');
        icatimf_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
        %- load .icatimef load-in parameters
        fprintf('Loading Time-Frequency Data...\n');
        tmp = load(icatimf_f,'-mat');
%         parameters = tmp.parameters;
        %- reshape data [pnts x chans x freq]
        fn = fieldnames(tmp);
        inds = find(contains(fieldnames(tmp),'comp'));
        test = tmp.(fn{inds(1)});
        rest_txf = zeros(size(test,3)*size(test,2),length(inds),size(test,1));
        for i = 1:length(inds)
            rest_txf(:,i,:) = reshape(tmp.(fn{inds(i)}),size(test,3)*size(test,2),1,size(test,1));
        end
        %- average over time, keep magnitude (not power -> would amplify outliers)
        base_txf_mean(1,:,:) = squeeze(mean(abs(rest_txf(10+1:end-10,:,:)),1));
        rest_txf = double.empty;
        par_save(base_txf_mean,fpath,sprintf('rest_avg_txf.mat'));
        %## (LOAD GAIT TIMEWARP) ======================================= %%
        fpath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
%         fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[TRIAL_TYPES{:}]);
        icatimf_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
        %- load .icatimef load-in parameters
        tmp = load(icatimf_f,'-mat');
%         parameters = tmp.parameters;
        %## (APPLY SPCA) =============================================== %%
        fprintf('Running SPCA on All Conditions\n');
        %- sPCA Algorithm
        [ERSP,GPM,gait_avg,output_struct] = apply_spca_cond_timewarp(tmp,base_txf_mean);
        [ERSP_corr, GPM_corr, PSC1, ~,COEFFS] = specPCAdenoising(ERSP);
        
        %## SAVE PCA INFORMATION
        gait_ersp_struct = [];
        gait_ersp_struct.ID         = EEG.subject;
        gait_ersp_struct.Noise_cov  = output_struct.baseline_cov;% noise cov for kernel computation
        gait_ersp_struct.F_Rest     = base_mean;
        gait_ersp_struct.TF         = gait_avg;
        gait_ersp_struct.ERSP_uncor = ERSP;
        gait_ersp_struct.GPM_uncor  = GPM;
        gait_ersp_struct.ERSP       = ERSP_corr;
        gait_ersp_struct.GPM        = GPM_corr;
        gait_ersp_struct.PSC1       = PSC1;
        gait_ersp_struct.chanlocs   = EEG.chanlocs;
        par_save(gait_ersp_struct,fpath,'gait_ersp_spca.mat');
        %## PLOT
        fig = figure(); set(gcf, 'position', [0 0 600 500]);
        plot(tmp.freqs, squeeze(base_mean)', 'k-');
        ylabel('Amplitude (\muV)'), ylabel('Frequency (Hz)');
        grid on; box off
        title('Baseline ERSP (rest)');
        exportgraphics(fig,[fpath filesep 'allcond_baseline_avgs.jpg']);
        %##
        CHAN_INT = randi(size(ERSP,2),1);
        fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],tmp.freqs);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_corr.jpg',CHAN_INT)]);
        %##
        fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],tmp.freqs);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_corr.jpg',CHAN_INT)]);
        %##
        fig = plot_txf(squeeze(ERSP(:,CHAN_INT,:)),[],tmp.freqs);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_ersp_orig.jpg',CHAN_INT)]);
        %##
        fig = plot_txf(squeeze(GPM(:,CHAN_INT,:)),[],tmp.freqs);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_gpm_orig.jpg',CHAN_INT)]);
        %##
        fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],tmp.freqs);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);

        %% DENOISE
        fprintf('\nRunning Denoising with sPCA\n');
        conds = unique({tmp.trialinfo.cond});
        for cond_i = 1:length(conds)
            COND_STR = conds{cond_i};
            [ERSP_corr,GPM_corr,gait_avg,output_struct] = apply_spca_cond_timewarp(tmp,base_txf_mean,...
                'COEFFS',COEFFS,...
                'COND_STR',COND_STR);
            struct_out = struct('ersp_corr',ERSP_corr,...
                'gpm_corr',GPM_corr,...
                'pc1',gait_avg,...
                'coeffs',COEFFS,...
                'apply_spca_cond',output_struct);
            par_save(struct_out,fpath,sprintf('cond%s_spca_ersp.mat',conds{cond_i}));

            %##
            fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) ERSP corrected',SPCA_PARAMS.condition_gait{cond_i}));
            if ~isempty(fpath)
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            end
            %##
            fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) GPM corrected',SPCA_PARAMS.condition_gait{cond_i}));
            if ~isempty(fpath)
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_gpm_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            end
            %##
            fig = plot_txf(squeeze(output_struct.erds_orig(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) ERSP original',SPCA_PARAMS.condition_gait{cond_i}));
            if ~isempty(fpath)
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            end
            %##
            fig = plot_txf(squeeze(output_struct.gpm_orig(:,CHAN_INT,:)),[],tmp.freqs);
            title(sprintf('%s) GPM original',SPCA_PARAMS.condition_gait{cond_i}));
            if ~isempty(fpath)
                exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_ersp_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
            end
            %##
            % fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],tmp.freqs);
            % title(sprintf('%s) PSC1 original',SPCA_PARAMS.condition_gait{cond_i}));
            % exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_psc1_orig.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i})]);
        end
        %## CLEANUP AND FREE UP DRIVE SPACE
        fpath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
        icatimf_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
        delete(icatimf_f);
        fpath = [epoched_fPath filesep 'rest'];
        icatimf_f = [fpath filesep sprintf('%s.icatimef',EEG.subject)];
        delete(icatimf_f);
        fprintf('Subject %s done. time: %0.2f',EEG.subject, toc(tt));
    catch e
        fprintf('\nError occured on subject %s\n%s\n',SUBJ_CHARS{subj_i},getReport(e));
    end
end
err = 1;
end
%% ===================================================================== %%
function [b] = validate_struct(x,DEFAULT_STRUCT)
    b = false;
    struct_name = inputname(2);
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field names
    chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
    if ~all(chk)
        fprintf(2,'\nFields for struct do not match for %s\n',struct_name);
        return
    end
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nValue must be type %s, but is type %s\n',class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end
