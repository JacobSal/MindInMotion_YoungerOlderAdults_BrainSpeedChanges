function [ALLEEG,pathsout] = main_func(EEG,conditions,save_dir,varargin)
%MAIN_FUNC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/06/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%## TIME
tic
%## define DEFAULTS
%-
SAVE_EEG = false;
PATH_EXT = 'dferris';
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}.';
cnd_validFcn = @(x) assert((iscell(x)), errorMsg);
%- MIM specific epoching
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
%- subj_i_epoch
PARSE_TYPE = 'Constant'; %
EPOCH_TIME_LIMITS = [-2,2];
TRIAL_LENGTH = 3*60; % trial length in seconds
PER_OVERLAP = 0.0; % percent overlap between epochs
TRIAL_BEGIN_STR = 'TrialStart';
TRIAL_END_STR = 'TrialEnd';
EVENT_TRIAL_PARSER = 'type';
EVENT_COND_PARSER = 'cond';
%- Connectivity Process
CONN_METHODS = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
STAT_ALPHA = 0.01; % 
% subj_i_genConnMeas
FREQS = (1:100);
CONN_COMPONENTS = EEG.chanlocs.urchan;
CNCTANL_TOOLBOX = 'sift';
DO_BOOTSTRAP = false;
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
ASSIGN_BOOTSTRAP_MEAN = false;
FREQ_BANDS = {FREQS;1:7;7:12;12:28;28:48;48:60};
% subj_i_genConnStats
DO_PHASE_RND = true;
errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
lpr_validFcn = @(x) assert(islogical(x),errorMsg);
%- HARD DEFINES
SUFFIX_PATH_SIFT = 'SIFT';
SUFFIX_PATH_EPOCHED = 'Epoched';
%## define Parser
p = inputParser;
%## define REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'conditions',cnd_validFcn);
addRequired(p,'save_dir',@ischar);
%## define OPTIONAL
addOptional(p,'conn_components',CONN_COMPONENTS,@isnumeric); 
%## define PARAMETER
addParameter(p,'SAVE_EEG',SAVE_EEG,@islogical);
addParameter(p,'PATH_EXT',PATH_EXT,@ischar);
%- MIM specific epoching
addParameter(p,'EVENT_CHAR',EVENT_CHAR,@ischar);
%- subj_i_epoch
addParameter(p,'PARSE_TYPE',PARSE_TYPE,@ischar);
addParameter(p,'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,@isnumeric);
addParameter(p,'TRIAL_LENGTH',TRIAL_LENGTH,@isnumeric);
addParameter(p,'PER_OVERLAP',PER_OVERLAP,(@(x) isnumeric(x) && x <= 1));
addParameter(p,'TRIAL_BEGIN_STR',TRIAL_BEGIN_STR,@ischar);
addParameter(p,'TRIAL_END_STR',TRIAL_END_STR,@ischar);
addParameter(p,'EVENT_TRIAL_PARSER',EVENT_TRIAL_PARSER,@ischar);
addParameter(p,'EVENT_COND_PARSER',EVENT_COND_PARSER,@ischar);
%- subj_i_genConnMeas
addParameter(p,'CONN_METHODS',CONN_METHODS,@iscell);
addParameter(p,'FREQS',FREQS,@isnumeric);
addParameter(p,'CNCTANL_TOOLBOX',CNCTANL_TOOLBOX,@ischar);
addParameter(p,'DO_BOOTSTRAP',DO_BOOTSTRAP,@islogical);
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@isnumeric);
addParameter(p,'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,@isnumeric);
addParameter(p,'ASSIGN_BOOTSTRAP_MEAN',ASSIGN_BOOTSTRAP_MEAN,@islogical);
%- subj_i_genConnStats
addParameter(p,'STAT_ALPHA',STAT_ALPHA,@isnumeric);
addParameter(p,'DO_PHASE_RND',DO_PHASE_RND,lpr_validFcn);
%## SET DEFAULTS
parse(p,EEG,conditions,save_dir,varargin{:});
%- pathing and saving
% OVERRIDE_SOURCE_EEG = false;
SAVE_EEG = p.Results.SAVE_EEG;
PATH_EXT = p.Results.PATH_EXT;
%- MIM specific epoching
EVENT_CHAR = p.Results.EVENT_CHAR;
%- subj_i_epoch
PARSE_TYPE = p.Results.PARSE_TYPE; %
EPOCH_TIME_LIMITS = p.Results.EPOCH_TIME_LIMITS;
TRIAL_LENGTH = p.Results.TRIAL_LENGTH; % trial length in seconds
PER_OVERLAP = p.Results.PER_OVERLAP; % percent overlap between epochs
TRIAL_BEGIN_STR = p.Results.TRIAL_BEGIN_STR;
TRIAL_END_STR = p.Results.TRIAL_END_STR;
EVENT_TRIAL_PARSER = p.Results.EVENT_TRIAL_PARSER;
EVENT_COND_PARSER = p.Results.EVENT_COND_PARSER;
%- Connectivity process
CONN_METHODS = p.Results.CONN_METHODS; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
STAT_ALPHA = p.Results.STAT_ALPHA;
% subj_i_genConnMeas
FREQS = p.Results.FREQS;
conn_components = p.Results.conn_components;
CNCTANL_TOOLBOX = p.Results.CNCTANL_TOOLBOX;
DO_BOOTSTRAP = p.Results.DO_BOOTSTRAP;
WINDOW_LENGTH = p.Results.WINDOW_LENGTH;
WINDOW_STEP_SIZE = p.Results.WINDOW_STEP_SIZE;
ASSIGN_BOOTSTRAP_MEAN = p.Results.ASSIGN_BOOTSTRAP_MEAN;
% subj_i_genConnStats
DO_PHASE_RND = p.Results.DO_PHASE_RND;
%- Get Anatomy Coordinates for different brain areas.
% out = get_Anatomy();
% BRAIN_CHARS = out.BrainAcronym;
% BRAIN_COORDS = out.coords_MNI;
%## PATHSOUT STRUCTURE
pathstruct = struct('filepath',[],...
       'filename',[],...
       'label',[]);
pathsout = struct('condition',[],...
                  'ica',pathstruct,...
                  'source',pathstruct,...
                  'epoch',pathstruct,...
                  'cnctanl',pathstruct,...
                  'META',[]);
for cond_i = 1:length(conditions)
    pathsout(cond_i).condition = conditions{cond_i};
    pathsout(cond_i).ica = pathstruct;
    pathsout(cond_i).source = pathstruct;
    pathsout(cond_i).epoch = pathstruct;
    pathsout(cond_i).cnctanl = pathstruct;
end

%% ===================================================================== %%
%% STEP 1) Reject Components
% (12/7/2022) JS, may want to remove components to improve computation
% speeds

%% STEP 2) EPOCHING
%* empty ALLEEG structure for repopulating
ALLEEG = cell(1,length(conditions)); 
fprintf(1,'\n==== EPOCHING %s  ====\n',EEG.subject);
for cond_i = 1:length(conditions)
    cond_name = conditions{cond_i};
    fprintf(1,'\n==== Processing ''%s'' ====\n',cond_name);    
    %- save path for epoched EEG
    epoch_save_dir = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED filesep PARSE_TYPE filesep conditions{cond_i}];
    if ~strncmp(computer,'PC',2)
        fPath = convertPath2UNIX(epoch_save_dir,PATH_EXT);
    else
        fPath = convertPath2Drive(epoch_save_dir,PATH_EXT);
    end
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %## EPOCH SUBJ DATA
    %- MIM specific function
    if strcmp(cond_name,'rest')
        [TMP_EEG] = MIM_epochData(EEG,PARSE_TYPE,...
                        cond_name,EPOCH_TIME_LIMITS,...
                        'TRIAL_LENGTH',TRIAL_LENGTH,...
                        'PERCENT_OVERLAP',PER_OVERLAP,...
                        'TRIAL_BEG_CHAR',TRIAL_BEGIN_STR,...
                        'TRIAL_END_CHAR',TRIAL_END_STR,...
                        'TRIAL_PARSE_CHAR',EVENT_TRIAL_PARSER,...
                        'CONDITION_PARSE_CHAR',EVENT_COND_PARSER);
        TMP_EEG.timewarp = struct([]);
    else
        %- find trials for desired condition
        tmp = strcmp({EEG.event.(EVENT_COND_PARSER)},cond_name);
        valid_points = [find(tmp,1,'first'),find(tmp,1,'last')];
        %## EPOCH AROUND A SPECIFC EVENT ONLY (NO TIME WARP)
%         TMP_EEG = pop_epoch(EEG,{EVENT_CHAR},EPOCH_TIME_LIMITS,'eventindices',valid_points(1):valid_points(2),'epochinfo','yes');
        
        %##  EPOCH DATA FOR GAIT CYCLES (TIME WARP)
        %seconds to epoch relative to first RHS
        epoch_times = [-1 3];
        events = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
        TMP_EEG = pop_epoch( EEG, {  'RHS'  }, epoch_times, 'eventindices',valid_points(1):valid_points(2), 'newname', sprintf('Merged datasets %s epochs',EEG.subject), 'epochinfo', 'yes');
        %- setup timewarp structure
        timewarp = make_timewarp(TMP_EEG,events,'baselineLatency',0, ...
                'maxSTDForAbsolute',3,...
                'maxSTDForRelative',3);
        %subject specific warpto (later use to help calc grand avg warpto across subjects)
        timewarp.warpto = median(timewarp.latencies);        
        goodepochs  = sort([timewarp.epochs]);
        %probably not needed?
        TMP_EEG = eeg_checkset(TMP_EEG);   
        sedi = setdiff(1:length(TMP_EEG.epoch),goodepochs);
        %reject outlier strides & 
        TMP_EEG = pop_select( TMP_EEG,'notrial',sedi);
        %- store timewarp structure in EEG
        TMP_EEG.timewarp = timewarp;
    end
    %- Update EEG structure
    TMP_EEG.condition = cond_name;
    %- assign values to output structure
    pathsout(cond_i).epoch.label = sprintf('Epoch_%s-%s',PARSE_TYPE,cond_name);
    pathsout(cond_i).epoch.filepath = fPath;
    pathsout(cond_i).epoch.filename = sprintf('%s_%s_EPOCH_TMPEEG.set',TMP_EEG.subject,cond_name);
    pathsout(cond_i).META.parseVariable = PARSE_TYPE;
    %- check to make sure a number isn't the first character
    chk = regexp(cond_name,'\d');
    if any(chk)
        cond_name = sprintf('x%s',cond_name);
    end
    pathsout(cond_i).META.condition = cond_name;
    pathsout(cond_i).META.epoch_times = EPOCH_TIME_LIMITS;
    %- save EEG
    fprintf(1,'Saving Subject %s\n',TMP_EEG.subject); 
    [TMP_EEG] = pop_saveset(TMP_EEG,'savemode','twofiles',...
        'filename',pathsout(cond_i).epoch.filename,...
        'filepath',pathsout(cond_i).epoch.filepath);
    ALLEEG{cond_i} = TMP_EEG;
end
fprintf(1,'\n==== DONE: EPOCHING ====\n');
%- concatenate ALLEEG
val = [];
ALLEEG = cellfun(@(x) [val; x], ALLEEG);
%% STEP 4) GENERATE CONNECTIVITY MEASURES
fprintf(1,'\n==== GENERATING CONNECTIVITY MEASURES FOR SUBJECT DATA ====\n');
%## ASSIGN PATH FOR SIFT
%- create an extra subdirectory for Epoched & SIFT'd data if one does not already exists
% fName = [ALLEEG(1).subject '_cnctanl_ALEEG.set'];
conn_save_dir = [save_dir filesep ALLEEG(1).subject filesep SUFFIX_PATH_SIFT];
%
if ~strncmp(computer,'PC',2)
    fPath = convertPath2UNIX(conn_save_dir,PATH_EXT);
else
    fPath = convertPath2Drive(conn_save_dir,PATH_EXT);
end
if ~exist(fPath,'dir')
    mkdir(fPath)
end
%- assign ALLEEG_cnct to every subject path
for cond_i = 1:length(ALLEEG)
    pathsout(cond_i).cnctanl.filepath = fPath;
    pathsout(cond_i).cnctanl.filename = ALLEEG(cond_i).filename;
end

%## FEVAL Connectivity
fprintf('%s) Processing componets:\n',ALLEEG(1).subject)
fprintf('%i,',conn_components'); fprintf('\n');
%- exit function if there are not enough components
if length(conn_components) < 2
    return;
end
%- Calculate Connectivity
switch CNCTANL_TOOLBOX
    case 'sift'
        [ALLEEG] = cnctanl_connMeas(ALLEEG,conn_components,CONN_METHODS,fPath,...
            'FREQS',FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE);
    case'bsmart'
        slide_step_size = ceil(WINDOW_STEP_SIZE*ALLEEG(1).srate);
        slide_win_len = ceil(WINDOW_LENGTH*ALLEEG(1).srate);
        max_morder = 30;
        ALLEEG = bsmart_cnctanl_connMeas(ALLEEG,conn_components,CONN_METHODS,fPath,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'MAX_MORDER',max_morder,...
            'SLIDE_WIN_LEN',slide_win_len,...
            'SLIDE_STEP_SIZE',slide_step_size);
end

%## (BOOTSTRAPPING) GROUP STATISTICS 
% (09/22/2022), JS, Might want to try and speed up bootstrap by
% adapting stat_surrogateGen.m to use parfor for bootstrapping... If
% possible? doesn't seem built well in the first place, so maybe?
% (10/27/2022), JS, Revist above note again!
% (12/7/2022), JS, need to update this boostrapping to include ALLEEG
if DO_BOOTSTRAP
    fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
    for cond_i=1:length(ALLEEG)
        %- calculate BootStrap distribution
        ALLEEG(cond_i) = cnctanl_groupStats(ALLEEG(cond_i),'BootStrap');    
        %- save BootStrap distribution    
        if SAVE_CONN_BOOTSTRAP
            TMP_EEG = ALLEEG(cond_i).CAT.PConn;
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_BootStrap');
        end
    end
    %## (2) ALTERNATIVE CODE
    %{
    % number of bootstrap samples to draw
    NumSamples = 200;

    % obtain the bootstrap distributions for each condition
    for cond_i=1:length(ALLEEG)
        ALLEEG(cond_i) = pop_stat_surrogateGen(ALLEEG(cond_i),'nogui', ...
            'modelingApproach', ALLEEG(cond_i).CAT.configs.est_fitMVAR, ...
            'connectivityModeling',ALLEEG(cond_i).CAT.configs.est_mvarConnectivity, ...
            'mode',{'Bootstrap' 'nperms' NumSamples 'saveTrialIdx' true}, ...
            'autosave',[], ...
            'verb',1);
        tmp = ALLEEG(cond_i).CAT.PConn;
        par_save(tmp,ALLEEG(cond_i).filepath,fName,[],'_BootStrap');
    end
    %}
    %- assign mean of bootstrap as Conn value
    if ASSIGN_BOOTSTRAP_MEAN
        for cond_i = 1:length(ALLEEG)
            ALLEEG(cond_i).CAT.Conn = stat_getDistribMean(ALLEEG(cond_i).CAT.PConn);
        end
    end
    
    %- clear bootstrap calculation
    EEG.CAT.PConn = [];
    
end

%## ASSIGN OUTPUT VALUES
for cond_i = 1:length(ALLEEG)
    pathsout(cond_i).META.connMethods = CONN_METHODS;
    pathsout(cond_i).META.connComponents = conn_components;
end

%## SAVE
if SAVE_EEG
    for cond_i = 1:length(ALLEEG)
        [~,~] = pop_saveset(ALLEEG(cond_i),...
            'filepath',ALLEEG(cond_i).filepath,'filename',ALLEEG(cond_i).filename,...
            'savemode','twofiles');
    end
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT DATA ====\n');          

%% GENERATE CONNECTIVITY STATISTICS
%## GENERATE PHASE RANDOMIZED DISTRIBUTION    
% see. stat_surrogateGen
% see. stat_surrogateStats
if DO_PHASE_RND
    fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    %## (1) ALTERNATIVE CODE 
    for cond_i=1:length(ALLEEG)
        %- Generate Phase Randomized Distribution
        fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
        [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'PhaseRnd');
        %- Save Phase randomized distribution
        TMP_EEG = ALLEEG(cond_i).CAT.PConn;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_PhaseRnd');
        fprintf('done.\n')
        %{
        [EEG,~] = cnctanl_groupStats(EEG,'PhaseRnd');
        %- Save Phase randomized distribution
        tmp = EEG.CAT.PConn;
        par_save(tmp,EEG.filepath,EEG.filename,'_PhaseRnd');
        fprintf('done.\n')
        %}
        %- Conduct Nonzero Test (Need Phase Randomized Distribution)
        fprintf('\n==== NONZERO STATISTICS ====\n')
        [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'NonZero');
        %- Save Nonzero Stats
        TMP_EEG = ALLEEG(cond_i).CAT.Stats;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_NonZero');
        fprintf('done.\n')
        %{
        [EEG,~] = cnctanl_groupStats(EEG,'NonZero');
        %- Save Phase randomized distribution
        tmp = EEG.CAT.Stats;
        par_save(tmp,EEG.filepath,EEG.filename,'_NonZero');
        fprintf('done.\n')
        %}
    end
    
    %## (2) ALTERNATIVE CODE
    % number of null distribution samples to draw
    % may be unstable on HiPerGator
    %{
    NumSamples = 200;
    for cond_i=1:length(ALLEEG)
        ALLEEG(cond_i) = pop_stat_surrogateGen(ALLEEG(cond_i),'nogui', ...
            'modelingApproach', ALLEEG(cond_i).CAT.configs.est_fitMVAR, ...
            'connectivityModeling',ALLEEG(cond_i).CAT.configs.est_mvarConnectivity, ...
            'mode',{'PhaseRand' 'nperms' NumSamples}, ...
            'autosave',[], ...
            'verb',1);
        %- Save Phase randomized distribution
        tmp = ALLEEG(cond_i).CAT.PConn;
        par_save(tmp,ALLEEG(cond_i).filepath,ALLEEG(cond_i).filename,'_PhaseRnd');
        fprintf('done.\n')
    end
    for cond_i=1:length(ALLEEG)
        ALLEEG(cond_i) = pop_stat_surrogateStats(ALLEEG(cond_i),'nogui',...
                            'statTest',...
                                {'Hnull',...
                                 'testMethod','quantile',...
                                 'tail','right',...
                                 'alpha',STAT_ALPHA,...
                                 'mcorrection','fdr',...
                                 'statcondargs',{}},...
                             'connmethods',CONN_METHODS, ...
                             'verb',1);
        %- Save Nonzero Stats
        tmp = ALLEEG(cond_i).CAT.Stats;
        par_save(tmp,ALLEEG(cond_i).filepath,ALLEEG(cond_i).filename,'_NonZero');
        fprintf('done.\n')
    end
    %}
end
%% CONNECTIVITY MATRICES & VISUALS
disp(ALLEEG);
%## Generate Connectivity Matrix
for conn_i = 1:length(CONN_METHODS)
    pathsout(cond_i).META.(CONN_METHODS{conn_i}).connExtract = struct('freqs',[],...
                                                                        'mats',[],...
                                                                        'sigs',[]);
    for freq_i = 1:length(FREQ_BANDS)
        for cond_i = 1:length(ALLEEG)
            disp(ALLEEG(cond_i).CAT.Conn);
            %- calculate average connnectivity for each component pair across time
            [statMat,extract_sig] = gen_connMatrix(ALLEEG(cond_i),CONN_METHODS{conn_i},FREQ_BANDS{freq_i},STAT_ALPHA);
            pathsout(cond_i).META.(CONN_METHODS{conn_i}).connExtract(freq_i).freqs = FREQ_BANDS{freq_i};
            pathsout(cond_i).META.(CONN_METHODS{conn_i}).connExtract(freq_i).mats = statMat;
            pathsout(cond_i).META.(CONN_METHODS{conn_i}).connExtract(freq_i).sigs = extract_sig;
        end
    end
end

%% WRAP UP 
%- clear STATS & PCONN for space saving.
%## SAVE
if true
    for cond_i = 1:length(ALLEEG)
        if isfield(ALLEEG(cond_i).CAT,'Stats')
            ALLEEG(cond_i).CAT = rmfield(ALLEEG(cond_i).CAT,'Stats');
        end
        if isfield(ALLEEG(cond_i).CAT,'PConn')
            ALLEEG(cond_i).CAT = rmfield(ALLEEG(cond_i).CAT,'PConn');
        end
        ALLEEG(cond_i).etc.CONN_PIPE = pathsout;
        [~,~] = pop_saveset(ALLEEG(cond_i),...
            'filepath',ALLEEG(cond_i).filepath,'filename',ALLEEG(cond_i).filename,...
            'savemode','twofiles');
    end
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');

%% DONE
fprintf('DONE: %s\n',ALLEEG(1).subject);
%## TIME
toc
end

