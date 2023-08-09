function [ALLEEG,t_out] = main_func(ALLEEG,save_dir,varargin)
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
% errorMsg = 'Value must be of format {CHAR1,CHAR2,...}.';
% cnd_validFcn = @(x) assert((iscell(x)), errorMsg);
%- connectivity methods
CONN_METHODS = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
%- connectivity statistics & surrogates params
%* 
DO_BOOTSTRAP = false;
ASSIGN_BOOTSTRAP_MEAN = false;
%* boolean for phase randomization generation
DO_PHASE_RND = true;
errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
lpr_validFcn = @(x) assert(islogical(x),errorMsg);
%*
STAT_ALPHA = 0.01;
%- connectivity model params
FREQS = (1:100);
CONN_COMPONENTS = (1:size(ALLEEG(1).icaweights,1))';
CNCTANL_TOOLBOX = 'sift';
WINDOW_LENGTH = 0.5; % length of sliding window in time(s)
WINDOW_STEP_SIZE = 0.025; % sliding amount in time(s)
FREQ_BANDS = {FREQS;1:7;7:12;12:28;28:48;48:60};
%## define Parser
p = inputParser;
%## define REQUIRED
addRequired(p,'ALLEEG',@isstruct);
% addRequired(p,'conditions',cnd_validFcn);
addRequired(p,'save_dir',@ischar);
%## define OPTIONAL

%## define PARAMETER
addParameter(p,'CONN_COMPONENTS',CONN_COMPONENTS,@isnumeric);
%- cnctanl_connMeas
addParameter(p,'CONN_METHODS',CONN_METHODS,@iscell);
addParameter(p,'FREQS',FREQS,@isnumeric);
addParameter(p,'CNCTANL_TOOLBOX',CNCTANL_TOOLBOX,@ischar);
addParameter(p,'DO_BOOTSTRAP',DO_BOOTSTRAP,@islogical);
addParameter(p,'WINDOW_LENGTH',WINDOW_LENGTH,@isnumeric);
addParameter(p,'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,@isnumeric);
addParameter(p,'ASSIGN_BOOTSTRAP_MEAN',ASSIGN_BOOTSTRAP_MEAN,@islogical);
%- cnctanl_groupStats
addParameter(p,'STAT_ALPHA',STAT_ALPHA,@isnumeric);
addParameter(p,'DO_PHASE_RND',DO_PHASE_RND,lpr_validFcn);
%## SET DEFAULTS
parse(p,ALLEEG,save_dir,varargin{:});
%- Connectivity process
CONN_METHODS = p.Results.CONN_METHODS; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
STAT_ALPHA = p.Results.STAT_ALPHA;
% subj_i_genConnMeas
FREQS = p.Results.FREQS;
CONN_COMPONENTS = p.Results.CONN_COMPONENTS;
CNCTANL_TOOLBOX = p.Results.CNCTANL_TOOLBOX;
DO_BOOTSTRAP = p.Results.DO_BOOTSTRAP;
WINDOW_LENGTH = p.Results.WINDOW_LENGTH;
WINDOW_STEP_SIZE = p.Results.WINDOW_STEP_SIZE;
ASSIGN_BOOTSTRAP_MEAN = p.Results.ASSIGN_BOOTSTRAP_MEAN;
% subj_i_genConnStats
DO_PHASE_RND = p.Results.DO_PHASE_RND;
%## TABLE VARS
t_fPaths = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
t_fNames = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
% t_conn_methods = cell(length(ALLEEG),1);
t_conn_comps = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
t_conn_mats = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
t_conn_freqs = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
t_conn_sigs = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
t_conn_meas = cell(length(ALLEEG)*length(CONN_METHODS)*length(FREQ_BANDS),1);
%% ===================================================================== %%
%## STEP 4) GENERATE CONNECTIVITY MEASURES
fprintf(1,'\n==== GENERATING CONNECTIVITY MEASURES FOR SUBJECT %s ====\n',ALLEEG(1).subject);
%## ASSIGN PATH FOR SIFT
%- make sure path is in right format and make sure it exists
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

%- assign table values for each trial in ALLEEEG
% for trial_i = 1:length(ALLEEG)
%     t_fPaths{trial_i} = ALLEEG(trial_i).filepath;
%     t_fNames{trial_i} = ALLEEG(trial_i).filename;
% %     t_conn_methods{trial_i} = CONN_METHODS;
%     t_conn_comps{trial_i} = CONN_COMPONENTS;
% end

%## FEVAL Connectivity
fprintf('%s) Processing componets:\n',ALLEEG(1).subject)
fprintf('%i,',CONN_COMPONENTS'); fprintf('\n');
%- exit function if there are not enough components
if length(CONN_COMPONENTS) < 2
    return;
end
%- select components from EEG
if length(CONN_COMPONENTS) ~= size(ALLEEG(1).icaweights,1)
    for cond_i = 1:length(ALLEEG)
        TMP_EEG = ALLEEG(cond_i);
        TMP_EEG = pop_subcomp(TMP_EEG,sort(CONN_COMPONENTS),0,1);
        %- Recalculate ICA Matrices && Book Keeping
        if isempty(TMP_EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',TMP_EEG.subject);
            TMP_EEG.icaact = (TMP_EEG.icaweights*TMP_EEG.icasphere)*TMP_EEG.data(TMP_EEG.icachansind,:);
            TMP_EEG.icaact = reshape(TMP_EEG.icaact,size(TMP_EEG.icaact,1),TMP_EEG.pnts,TMP_EEG.trials);
        end
        ALLEEG(cond_i) = TMP_EEG;
    end
end
%- Calculate Connectivity
switch CNCTANL_TOOLBOX
    case 'sift'
        [TMP_EEG] = cnctanl_conn_meas(ALLEEG,CONN_COMPONENTS,CONN_METHODS,save_dir,...
            'FREQS',FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE);
        for trial_i = 1:length(ALLEEG)
            ALLEEG(trial_i).CAT = TMP_EEG(trial_i).CAT;
        end
    case'bsmart'
        slide_step_size = ceil(WINDOW_STEP_SIZE*ALLEEG(1).srate);
        slide_win_len = ceil(WINDOW_LENGTH*ALLEEG(1).srate);
        max_morder = 30;
        ALLEEG = bsmart_cnctanl_connMeas(ALLEEG,CONN_COMPONENTS,CONN_METHODS,save_dir,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'MAX_MORDER',max_morder,...
            'SLIDE_WIN_LEN',slide_win_len,...
            'SLIDE_STEP_SIZE',slide_step_size);
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY MEASURES FOR SUBJECT %S ====\n',ALLEEG(1).subject);          
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
        %- BootStrap corrected statistics per condition
%         fprintf('\n==== BETWEEN CONDITION STATISTICS ====\n')
%         [ALLEEG(trial_i),~] = cnctanl_groupStats(ALLEEG(trial_i),'BtwnCond');
%         %- Save Nonzero Stats
%         TMP_EEG = ALLEEG(trial_i).CAT.Stats;
%         fName = strsplit(ALLEEG(trial_i).filename,'.'); fName = [fName{1} '.mat'];
%         par_save(TMP_EEG,ALLEEG(trial_i).filepath,fName,'_BtwnCond');
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
    fprintf('\n==== DONE: CALCULATING BOOTSTRAP MEASURES ====\n')
end
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
    fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
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
            %- store outputs into a tabularable ("table able") format.
            t_conn_comps{cnt} = CONN_COMPONENTS;
            t_fPaths{cnt} = ALLEEG(trial_i).filepath;
            t_fNames{cnt} = ALLEEG(trial_i).filename;
            t_conn_freqs{cnt} = FREQ_BANDS{freq_i};
            t_conn_mats{cnt} = statMat;
            t_conn_sigs{cnt} = extract_sig;
            t_conn_meas{cnt} = CONN_METHODS{conn_i};
            cnt = cnt + 1;
        end
    end
end

%% ===================================================================== %%
%## STEP 7) WRAP UP 
%- create table
% t_fPaths = repmat(t_fPaths,length(t_conn_mats)/length(t_fPaths),1);
% t_fNames = repmat(t_fNames,length(t_conn_mats)/length(t_fNames),1);
% t_conn_methods = repmat(t_conn_methods,length(t_conn_mats)/length(t_conn_methods),1);
% t_conn_comps = repmat(t_conn_comps,length(t_conn_mats)/length(t_conn_comps),1);
t_out = table(t_fPaths,t_fNames,t_conn_comps,t_conn_meas,t_conn_freqs,t_conn_sigs,t_conn_mats);
% disp(t_out)
%## REMOVE FIELDS && SAVE
% ALLEEG(trial_i).CAT.conn_table = t_out;
% par_save(t_out,ALLEEG(trial_i).filepath,fName,'_conntable');
for trial_i = 1:length(ALLEEG)
    %- clear STATS & PCONN for space saving.
    if isfield(ALLEEG(trial_i).CAT,'Stats')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'Stats');
    end
    if isfield(ALLEEG(trial_i).CAT,'PConn')
        ALLEEG(trial_i).CAT = rmfield(ALLEEG(trial_i).CAT,'PConn');
    end
end
fprintf(1,'\n==== DONE: GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');

%% DONE
fprintf('DONE: %s\n',ALLEEG(1).subject);
%## TIME
toc
end
