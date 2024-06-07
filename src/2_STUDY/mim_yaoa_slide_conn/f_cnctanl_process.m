%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_slide_conn/run_f_cnctanl_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- connecitivty modeling
CONN_FREQS = (4:50);
CONN_METHODS = {'dDTF08','S'}; %{'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
WINDOW_LENGTH = 10;
WINDOW_STEP_SIZE = 10;
VERBOSITY_LEVEL = 1;
MORDER = 10;
%-
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
RESAMLE_FREQ = 225;
DO_BOOTSTRAP = false;
DO_PHASE_RND = false;
DO_STANDARD_TRIALS = false;
DO_RECOMPUTE = true;
%## PREPARE DATA PARAMS
DEF_PREPDATA = struct('VerbosityLevel',VERBOSITY_LEVEL,...
             'SignalType',{{'Components'}},...
             'VariableNames',[],...
             'Detrend',{{'verb',VERBOSITY_LEVEL,'method',{'linear'},...
                    'plot',false}},...
             'NormalizeData',{{'verb',VERBOSITY_LEVEL,'method',{'time', 'ensemble'}}},...
             'resetConfigs',true,...
             'badsegments',[],...
             'newtrials',[],...
             'equalizetrials',false);
%## ESTIAMTE MODEL ORDER PARAMS
DEF_ESTSELMOD = struct('modelingApproach',{{'Segmentation VAR',...
                        'algorithm',{{'ARfit'}},...
                        'winStartIdx',[],...
                        'winlen',WINDOW_LENGTH,...
                        'winstep',WINDOW_STEP_SIZE,...
                        'taperfcn','rectwin',...
                        'epochTimeLims',[],...
                        'prctWinToSample',100,...
                        'normalize',[],...
                        'detrend',{'method','linear'},...
                        'verb',VERBOSITY_LEVEL}},...
    'morderRange',[1,20],...
    'downdate',true,...
    'RunInParallel',{{'profile','local',...
        'numWorkers',SLURM_POOL_SIZE}},...
    'icselector',{{'aic','hq'}},...
    'winStartIdx',[],...
    'epochTimeLims',[],...
    'prctWinToSample',80,...
    'plot',[],...
    'verb',VERBOSITY_LEVEL);
%- (05/01/2024) changing alg from {'Group Lasso (DAL/SCSA)'} to {'ARfit'}
%could also try {'Vieira-Morf'}
%##
DEF_PLOTORDERCRIT = struct('conditions',{TRIAL_TYPES},    ...
                            'icselector', {DEF_ESTSELMOD.icselector},  ...
                            'minimizer',{{'min'}}, ...
                            'prclim', 90);
%##
DEF_ESTDISPMVAR_CHK = struct('morder',MORDER,...
        'winlen',WINDOW_LENGTH,'winstep',WINDOW_STEP_SIZE,'verb',1);
%##
DEF_ESTFITMVAR = struct('connmethods',{CONN_METHODS}, ...
            'absvalsq',true,           ...
            'spectraldecibels',true,   ...
            'freqs',CONN_FREQS,        ...
            'verb',1);
%## PATHING
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%## soft define
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
conn_fig_dir = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'conn_valid_slide'];
%- load study file
STUDY_FNAME_LOAD = 'all_comps_study';
STUDY_FNAME_SAVE = 'slide_conn_study';
study_fpath = [studies_fpath filesep sprintf('%s',study_dir_name)];
%- load cluster
CLUSTER_K = 11;
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%- create new study directory
if ~exist(conn_fig_dir,'dir')
    mkdir(conn_fig_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
% if ~exist([STUDY_FPATH filesep STUDY_FNAME '.study'],'file')
%     error('ERROR. study file does not exist');
%     exit(); %#ok<UNRCH>
% else
%     %## LOAD STUDY
%     if ~ispc
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',STUDY_FPATH);
%     else
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',STUDY_FPATH);
%     end
%     cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
%     MAIN_STUDY.cluster = cl_struct;
%     [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
% end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[study_fpath filesep sprintf('%s_UNIX.study',STUDY_FNAME_LOAD)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[study_fpath filesep sprintf('%s.study',STUDY_FNAME_LOAD)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,~,~,valid_cls] = eeglab_get_cluster_comps(STUDY);
%% INITIALIZE PARFOR LOOP VARS
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
LOOP_VAR = 1:length(STUDY.datasetinfo);
tmp = cell(1,length(STUDY.datasetinfo));
rmv_subj = zeros(1,length(STUDY.datasetinfo));
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% for subj_i = 1:length(LOOP_VAR)
    %- Parse out components
    components = comps_out(:,subj_i);
    components = sort(components(components ~= 0));
    eeg_savefpath = [fPaths{subj_i} filesep 'conn_slide'];
    if ~exist(eeg_savefpath,'dir')
        mkdir(eeg_savefpath)
    end
    %## LOAD EEG DATA
    if ~exist([eeg_savefpath filesep fNames{subj_i}],'file') || DO_RECOMPUTE
        EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        %- Recalculate ICA Matrices && Book Keeping
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact,size(EEG.icaact,1),EEG.pnts,EEG.trials);
        end
        %- for some reason need to pop this for every loop?
        pop_editoptions('option_computeica', 1);
        %## RESAMPLE
        %- (04/24/2024) JS, the rationale for resampling is majorly based
        %on computation times, however, a argument could be made that by
        %removing higher frequency (i.e., downsampling) data the
        %connectivity model may fit more accurately to lower frequency
        %content. The latter point has not been proven, merely suggested in
        %the literature.
        %(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10044923/)
        % EEG = pop_resample(EEG,RESAMLE_FREQ);
        try
            %## RUN MAIN_FUNC
            [TMP,t_out] = cnctanl_sift_pipe(EEG,components,conn_fig_dir,...
                'DO_PHASE_RND',DO_PHASE_RND,...
                'DO_BOOTSTRAP',DO_BOOTSTRAP,...
                'PREPDATA',DEF_PREPDATA,...
                'ESTSELMOD',DEF_ESTSELMOD,...
                'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
                'ESTFITMVAR',DEF_ESTFITMVAR,...
                'PLOTORDERCRIT',DEF_PLOTORDERCRIT);
            % BIG_CAT = cat(1,TMP(:).CAT);
            EEG.etc.COND_CAT = cat(1,TMP(:).CAT);
            EEG.etc.conn_table = t_out;
            EEG.etc.conn_meta.comps_out = comps_out;
            fName = strsplit(EEG.filename,'.'); fName = [fName{1} '.mat'];
            par_save(t_out,eeg_savefpath,fName,'_conntable');
            [EEG] = pop_saveset(EEG,...
                'filepath',eeg_savefpath,'filename',EEG.filename,...
                'savemode','twofiles');
            tmp{subj_i} = EEG;
            EEG = struct.empty;
            TMP = struct.empty;

        catch e
            rmv_subj(subj_i) = 1;
            EEG.CAT = struct([]);
            tmp{subj_i} = [];
            fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'error. on subject %s\n',...
                     'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
            disp(e);
        end
    else
        tmp{subj_i} = pop_loadset('filepath',eeg_savefpath,'filename',fNames{subj_i});
    end
end
pop_editoptions('option_computeica',0);
%% SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
try
    fprintf('Bugged Subjects:\n');
    fprintf('%s\n',STUDY.datasetinfo(cellfun(@isempty,tmp)).subject);
    bugged_subjs = STUDY.datasetinfo(cellfun(@isempty,tmp)).subject;
    subj_keep = find(~cellfun(@isempty,tmp));
    tmp = tmp(~cellfun(@isempty,tmp));
catch e
    fprintf('\n%s\n',getReport(e))
    bugged_subjs = [];
    subj_keep = find(~cellfun(@isempty,tmp));
    tmp = tmp(~cellfun(@isempty,tmp));
    fprintf('No Bugged Subjects.\n');
end
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding %s %s\n',EEG.subject,addFs{j})
        end
    end 
    tmp{subj_i} = EEG;
end
tmp = cellfun(@(x) [[]; x], tmp);
%##
tmp = eeg_checkset(tmp,'eventconsistency');
[NEW_STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',STUDY_FNAME_SAVE,...
                                'filename',STUDY_FNAME_SAVE,...
                                'filepath',study_fpath);
%## ASSIGN PARAMETERS
NEW_STUDY.etc.a_epoch_process.epoch_chars = TRIAL_TYPES;
NEW_STUDY.etc.d_cnctanl_process.params = struct('DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'PREPDATA',DEF_PREPDATA,...
            'ESTSELMOD',DEF_ESTSELMOD,...
            'ESTDISPMVAR_CHK',DEF_ESTDISPMVAR_CHK,...
            'ESTFITMVAR',DEF_ESTFITMVAR,...
            'PLOTORDERCRIT',DEF_PLOTORDERCRIT,...
            'comps',comps_out);
[NEW_STUDY,ALLEEG] = parfunc_save_study(NEW_STUDY,ALLEEG,...
                                            STUDY_FNAME_SAVE,study_fpath,...
                                            'STUDY_COND',[]);                          
%% Version History
%{
v1.0; (11/11/2022), JS: really need to consider updating bootstrap
    algorithm with parallel computing. Taking ~ 1 day per
    condition for all subjects and the bottle neck is entirely the
    bootstrap.

    Note: validateattributes and assert functions may be helpful
    in more clearly defining function inputs.
        e.g.  DO_PHASE_RND = true;
          errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
          validationFcn = @(x) assert(islogical(x),errorMsg);
v1.0; (12/5/2022) Need to adapt this to include all conditions
    within each SUBJ structure so connectivity can be calculated
    for the ALLEEG structure rather than the EEG structure.
    *** maybe try to ditch the SUBJ strucutre entirely for this
    round?
v1.0.01132023.0 : Initializing versioning for future iterations.
%}

