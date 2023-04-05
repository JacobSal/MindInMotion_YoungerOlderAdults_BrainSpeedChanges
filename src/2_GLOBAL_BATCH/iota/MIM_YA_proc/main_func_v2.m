function [exitcode] = main_func(ALLEEG,conditions,save_dir,varargin)
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
fid = fopen('./main_func_out.txt');
%## define DEFAULTS
%-
SAVE_EEG = false;
PATH_EXT = 'dferris';
errorMsg = 'Value must be of format {CHAR1,CHAR2,...}.';
cnd_validFcn = @(x) assert((iscell(x)), errorMsg);
%- Connectivity Process
CONN_METHODS = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
STAT_ALPHA = 0.01; % 
% subj_i_genConnMeas
FREQS = (1:100);
CONN_COMPONENTS = ALLEEG.chanlocs.urchan;
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
%## define Parser
p = inputParser;
%## define REQUIRED
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'conditions',cnd_validFcn);
addRequired(p,'save_dir',@ischar);
%## define OPTIONAL
addOptional(p,'conn_components',CONN_COMPONENTS,@isnumeric); 
%## define PARAMETER
addParameter(p,'SAVE_EEG',SAVE_EEG,@islogical);
addParameter(p,'PATH_EXT',PATH_EXT,@ischar);
%- MIM specific epoching
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
parse(p,ALLEEG,conditions,save_dir,varargin{:});
%- pathing and saving
% OVERRIDE_SOURCE_EEG = false;
SAVE_EEG = p.Results.SAVE_EEG;
PATH_EXT = p.Results.PATH_EXT;
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
%% STEP 4) GENERATE CONNECTIVITY MEASURES
fprintf(fid,'\n==== GENERATING CONNECTIVITY MEASURES FOR SUBJECT DATA ====\n');
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
fprintf(fid,'%s) Processing componets:\n',ALLEEG(1).subject);
fprintf(fid,'%i,',conn_components'); fprintf(fid,'\n');
%- exit function if there are not enough components
if length(conn_components) < 2
    exitcode = 0;
    exit(exitcode)
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
    fprintf(fid,'\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
    for cond_i=1:length(ALLEEG)
        %- calculate BootStrap distribution
        ALLEEG(cond_i) = cnctanl_groupStats(ALLEEG(cond_i),'BootStrap');    
        %- save BootStrap distribution 
        TMP_EEG = ALLEEG(cond_i).CAT.PConn;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_BootStrap');
    end
    %- assign mean of bootstrap as Conn value
    if ASSIGN_BOOTSTRAP_MEAN
        for cond_i = 1:length(ALLEEG)
            ALLEEG(cond_i).CAT.Conn = stat_getDistribMean(ALLEEG(cond_i).CAT.PConn);
        end
    end
    %- clear bootstrap calculation
    ALLEEG.CAT.PConn = [];
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
fprintf(fid,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT DATA ====\n');          

%% GENERATE CONNECTIVITY STATISTICS
%## GENERATE PHASE RANDOMIZED DISTRIBUTION    
% see. stat_surrogateGen
% see. stat_surrogateStats
if DO_PHASE_RND
    fprintf(fid,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    %## (1) ALTERNATIVE CODE 
    for cond_i=1:length(ALLEEG)
        %- Generate Phase Randomized Distribution
        fprintf(fid,'\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n');
        [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'PhaseRnd');
        %- Save Phase randomized distribution
        TMP_EEG = ALLEEG(cond_i).CAT.PConn;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_PhaseRnd');
        fprintf(fid,'done.\n');
        %- Conduct Nonzero Test (Need Phase Randomized Distribution)
        fprintf(fid,'\n==== NONZERO STATISTICS ====\n');
        [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'NonZero');
        %- Save Nonzero Stats
        TMP_EEG = ALLEEG(cond_i).CAT.Stats;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(cond_i).filepath,fName,'_NonZero');
        fprintf(fid,'done.\n');
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
fprintf(fid,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');

%% DONE
fprintf(fid,'DONE: %s\n',ALLEEG(1).subject);
%## TIME
toc
exitcode = 1;
exit(exitcode)
end

