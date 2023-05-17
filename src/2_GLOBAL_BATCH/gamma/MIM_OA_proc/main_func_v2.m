function [ALLEEG] = main_func_v2(ALLEEG,save_dir,varargin)
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
% errorMsg = 'Value must be of format {CHAR1,CHAR2,...}.';
% cnd_validFcn = @(x) assert((iscell(x)), errorMsg);
%- Connectivity Process
CONN_METHODS = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
STAT_ALPHA = 0.01; % 
% subj_i_genConnMeas
FREQS = (1:100);
CONN_COMPONENTS = ALLEEG(1).chanlocs.urchan;
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
% addRequired(p,'conditions',cnd_validFcn);
addRequired(p,'save_dir',@ischar);
%## define OPTIONAL
addOptional(p,'conn_components',CONN_COMPONENTS,@isnumeric); 
%## define PARAMETER
addParameter(p,'SAVE_EEG',SAVE_EEG,@islogical);
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
parse(p,ALLEEG,save_dir,varargin{:});
%- pathing and saving
% OVERRIDE_SOURCE_EEG = false;
SAVE_EEG = p.Results.SAVE_EEG;
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
%## TABLE VARS
t_fPaths = cell(1,length(ALLEEG));
t_fNames = cell(1,length(ALLEEG));
t_conn_methods = cell(1,length(ALLEEG));
t_conn_components = cell(1,length(ALLEEG));
t_conn_mats = cell(1,length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS));
t_conn_freqs = cell(1,length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS));
t_conn_sigs = cell(1,length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS));
%% ===================================================================== %%
%## STEP 4) GENERATE CONNECTIVITY MEASURES
fprintf(1,'\n==== GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s ====\n',ALLEEG.subject);
%## ASSIGN PATH FOR SIFT
%- make sure path is in right format and make sure it exists
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

%- assign ALLEEG_cnct to every subject path
for trial_i = 1:length(ALLEEG)
    t_fPaths{trial_i} = ALLEEG(trial_i).filepath;
    t_fNames{trial_i} = ALLEEG(trial_i).filename;
    %- assign new fName
%     tmp = strsplit(ALLEEG(trial_i).filename,'.');
%     tmp{1} = [tmp{1} '_sift'];
%     fName = join(tmp,'.');
%     fName = fName{1};
%     t_fNames{trial_i} = fName;
%     ALLEEG(trial_i).filename = fName;
end

%## FEVAL Connectivity
fprintf('%s) Processing componets:\n',ALLEEG(1).subject)
fprintf('%i,',conn_components'); fprintf('\n');
%- exit function if there are not enough components
if length(conn_components) < 2
    return;
end
%- select components from EEG
if length(conn_components) ~= size(ALLEEG.icaweights,1)
    for cond_i = 1:length(ALLEEG)
        ALLEEG(cond_i) = pop_subcomp(ALLEEG(cond_i), sort(conn_components), 0, 1);
    end
end
%- Calculate Connectivity
switch CNCTANL_TOOLBOX
    case 'sift'
        [TMPEEG] = cnctanl_connMeas(ALLEEG,conn_components,CONN_METHODS,save_dir,...
            'FREQS',FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE);
        for trial_i = 1:length(ALLEEG)
            ALLEEG(trial_i).CAT = TMPEEG(trial_i).CAT;
        end
        clear TMPEEG
    case'bsmart'
        slide_step_size = ceil(WINDOW_STEP_SIZE*ALLEEG(1).srate);
        slide_win_len = ceil(WINDOW_LENGTH*ALLEEG(1).srate);
        max_morder = 30;
        ALLEEG = bsmart_cnctanl_connMeas(ALLEEG,conn_components,CONN_METHODS,save_dir,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'MAX_MORDER',max_morder,...
            'SLIDE_WIN_LEN',slide_win_len,...
            'SLIDE_STEP_SIZE',slide_step_size);
end
%% ===================================================================== %%
%## STEP 5.a) (BOOTSTRAPPING) GROUP STATISTICS 
% (09/22/2022), JS, Might want to try and speed up bootstrap by
% adapting stat_surrogateGen.m to use parfor for bootstrapping... If
% possible? doesn't seem built well in the first place, so maybe?
% (10/27/2022), JS, Revist above note again!
% (12/7/2022), JS, need to update this boostrapping to include ALLEEG
if DO_BOOTSTRAP
    fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
    for trial_i=1:length(ALLEEG)
        %- calculate BootStrap distribution
        ALLEEG(trial_i) = cnctanl_groupStats(ALLEEG(trial_i),'BootStrap');    
        %- save BootStrap distribution 
        TMP_EEG = ALLEEG(trial_i).CAT.PConn;
        fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(trial_i).filepath,fName,'_BootStrap');
    end
    %- assign mean of bootstrap as Conn value
    if ASSIGN_BOOTSTRAP_MEAN
        for trial_i = 1:length(ALLEEG)
            ALLEEG(trial_i).CAT.Conn = stat_getDistribMean(ALLEEG(trial_i).CAT.PConn);
        end
    end 
    for trial_i = 1:length(ALLEEG)
        %- clear bootstrap calculation
        ALLEEG(trial_i).CAT.PConn = [];
    end    
end

%## ASSIGN OUTPUT VALUES
for trial_i = 1:length(ALLEEG)
    t_conn_methods{trial_i} = CONN_METHODS;
    t_conn_components{trial_i} = conn_components;
end

%## SAVE
% if SAVE_EEG
%     for trial_i = 1:length(ALLEEG)
%         [~,~] = pop_saveset(ALLEEG(trial_i),...
%             'filepath',ALLEEG(trial_i).filepath,'filename',ALLEEG(trial_i).filename,...
%             'savemode','twofiles');
%     end
% end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT DATA ====\n');          

%% ===================================================================== %%
%## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION    
% see. stat_surrogateGen
% see. stat_surrogateStats
if DO_PHASE_RND
    fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
    %## (1) ALTERNATIVE CODE 
    for trial_i=1:length(ALLEEG)
        %- Generate Phase Randomized Distribution
        fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
        [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'PhaseRnd');
        %- Save Phase randomized distribution
        TMP_EEG = ALLEEG(trial_i).CAT.PConn;
        fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(trial_i).filepath,fName,'_PhaseRnd');
        fprintf('done.\n')
        %- Conduct Nonzero Test (Need Phase Randomized Distribution)
        fprintf('\n==== NONZERO STATISTICS ====\n')
        [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'NonZero');
        %- Save Nonzero Stats
        TMP_EEG = ALLEEG(trial_i).CAT.Stats;
        fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(TMP_EEG,ALLEEG(trial_i).filepath,fName,'_NonZero');
        fprintf('done.\n')
    end
    clear TMP_EEG
end
%% ===================================================================== %%
%## STEP 6) CONNECTIVITY MATRICES & VISUALS
disp(ALLEEG);
%## Generate Connectivity Matrix
cnt = 1;
for conn_i = 1:length(CONN_METHODS)
    for freq_i = 1:length(FREQ_BANDS)
        for trial_i = 1:length(ALLEEG)
            disp(ALLEEG(trial_i).CAT.Conn);
            %- calculate average connnectivity for each component pair across time
            [statMat,extract_sig] = gen_connMatrix(ALLEEG(trial_i),CONN_METHODS{conn_i},FREQ_BANDS{freq_i},STAT_ALPHA);
            t_conn_freqs{cnt} = FREQ_BANDS{freq_i};
            t_conn_mats{cnt} = statMat;
            t_conn_sigs{cnt} = extract_sig;
            cnt = cnt + 1;
        end
    end
end

%% ===================================================================== %%
%## STEP 7) WRAP UP 
%- create table
t_fPaths = repmat(t_fPaths,1,length(t_conn_mats)/length(t_fPaths));
t_fNames = repmat(t_fNames,1,length(t_conn_mats)/length(t_fNames));
t_conn_methods = repmat(t_conn_methods,1,length(t_conn_mats)/length(t_conn_methods));
t_conn_components = repmat(t_conn_components,1,length(t_conn_mats)/length(t_conn_components));
t_out = table(t_fPaths,t_fNames,t_conn_methods,t_conn_components,t_conn_freqs,t_conn_sigs,t_conn_mats);
%## REMOVE FIELDS && SAVE
for trial_i = 1:length(ALLEEG)
    %- clear STATS & PCONN for space saving.
    if isfield(ALLEEG(trial_i).CAT,'Stats')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'Stats');
    end
    if isfield(ALLEEG(trial_i).CAT,'PConn')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'PConn');
    end
    ALLEEG(trial_i).CAT.conn_table = t_out;
%     if false
%         [~] = pop_saveset(ALLEEG(trial_i),...
%             'filepath',ALLEEG(trial_i).filepath,'filename',ALLEEG(trial_i).filename,...
%             'savemode','twofiles');
%    end
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');

%% DONE
fprintf('DONE: %s\n',ALLEEG(1).subject);
%## TIME
toc
end
