%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/MIM_YA/run_e_cnctanl_phaserand_surr.sh

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
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = SCRIPT_DIR;
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(SCRIPT_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- connecitivty modeling
CONN_FREQS = (1:64);
% (01/19/2024) JS, trying 4:50 from 1:100(see. Steven Peterson 2019 NeuroImage)
FREQ_BANDS = {CONN_FREQS;4:8;8:13;13:28};
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
CONN_METHODS = {'dDTF08','S'}; %{'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
%## HIGH LEVEL VARS
%- connectivity statistics & surrogates params
%* boolean for model bootstrap distribution
DO_BOOTSTRAP = true;
FORCE_BOOTCALC = false;
ASSIGN_BOOTSTRAP_MEAN = true;
%* boolean for phase randomization generation
DO_PHASE_RND = true;
FORCE_PHASECALC = false;
%*
% conn_components = (1:size(ALLEEG(1).icaweights,1))';
% conn_estimators = {'dDTF08','dDTF','GGC'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
WINDOW_LENGTH = 0.5; % time (s)
WINDOW_STEP_SIZE = 0.025; % time (s)
N_PERMS_PHASE_RND = 200; % number of null distribution samples to draw
N_PERMS_BOOTSRAP = 200;

%- BOOTSTRAP
STAT_BOOTSTRAP_CFG = [];
STAT_BOOTSTRAP_CFG.verb = 1;
STAT_BOOTSTRAP_CFG.mode = {'Bootstrap','nperms',N_PERMS_BOOTSRAP,'saveTrialIdx',false};
%- PHASE RANDOMIZATION
STAT_PHASERND_CFG = [];
STAT_PHASERND_CFG.mode.arg_direct = 1;
STAT_PHASERND_CFG.mode.nperms = N_PERMS_PHASE_RND;
STAT_PHASERND_CFG.mode.arg_selection = 'PhaseRand';
STAT_PHASERND_CFG.modelingApproach = [];%ALLEEG(trial_i).CAT.configs.est_fitMVAR;
STAT_PHASERND_CFG.connectivityModeling = [];%ALLEEG(trial_i).CAT.configs.est_mvarConnectivity;
STAT_PHASERND_CFG.verb = 1;

% (01/29/2024) changing these to 200 to save on computation time, but 2000
% iterations may be more robust.
MORDER = [];
%- datetime override
study_fname_1 = 'epoch_study';
study_dir_fname = '01232023_MIM_YAN32_antsnormalize_iccREMG0p4_powpow0p3_conn';
%- Subject Directory information
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_load_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_fname)];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_fname)];
conn_save_dir = [study_load_dir filesep 'conn_data'];
%- load cluster
CLUSTER_DIR = [STUDIES_DIR filesep sprintf('%s',study_dir_fname) filesep 'cluster'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
CLUSTER_STUDY_DIR = [CLUSTER_DIR filesep 'icrej_5'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
if ~exist(conn_save_dir,'dir')
    mkdir(conn_save_dir);
end
%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
if ~exist([CLUSTER_STUDY_DIR filesep CLUSTER_STUDY_FNAME '.study'],'file')
    error('ERROR. study file does not exist');
    exit(); %#ok<UNRCH>
else
%     if ~ispc
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',study_load_dir);
%     else
%         [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',study_load_dir);
%     end
    %## LOAD STUDY
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '_UNIX.study'],'filepath',CLUSTER_STUDY_DIR);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '.study'],'filepath',CLUSTER_STUDY_DIR);
    end
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
end
%% INITIALIZE PARFOR LOOP VARS
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(MAIN_ALLEEG);
tmp = cell(1,length(MAIN_ALLEEG));
rmv_subj = zeros(1,length(MAIN_ALLEEG));
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
%% CONNECTIVITY MAIN FUNC
fprintf('Computing Connectivity\n');
pop_editoptions('option_computeica', 1);
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% for subj_i = LOOP_VAR
    %## PARAMS
    stat_pr_cfg = STAT_PHASERND_CFG;
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    ALLEEG = cell(length(EEG.etc.cond_files),1);
    %- re-epoch
    for cond_i = 1:length(EEG.etc.cond_files)
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
        ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
        fprintf('done.\n')
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    
    try
        %% ===================================================================== %%
        %## STEP 5.a) (BOOTSTRAPPING) GROUP STATISTICS 
        % (09/22/2022), JS, Might want to try and speed up bootstrap by
        % adapting stat_surrogateGen.m to use parfor for bootstrapping... If
        % possible? doesn't seem built well in the first place, so maybe?
        % (10/27/2022), JS, Revist above note again!
        % (12/7/2022), JS, need to update this boostrapping to include ALLEEG
        if DO_BOOTSTRAP
            fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
            for cond_i=1:length(ALLEEG)
                %- clear PConn
                ALLEEG(cond_i).CAT.PConn  = [];
                fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
                chk = ~exist([ALLEEG(cond_i).filepath filesep fName,'_BootStrap.mat'],'file') || FORCE_BOOTCALC
                if chk
                    [PConn,~] = feval(@stat_surrogateGen,'ALLEEG',ALLEEG(cond_i),STAT_BOOTSTRAP_CFG);
                    ALLEEG(cond_i).CAT.PConn = PConn;
                    %- save BootStrap distribution 
                    bootstrap_dist = ALLEEG(cond_i).CAT.PConn;
                    par_save(bootstrap_dist,ALLEEG(cond_i).filepath,fName,'_BootStrap');
                else
                    ALLEEG(cond_i).CAT.PConn = par_load(ALLEEG(cond_i).filepath,fName,'_BootStrap.mat');
                end
            end
            %- assign mean of bootstrap as Conn value
            if ASSIGN_BOOTSTRAP_MEAN
                for cond_i = 1:length(ALLEEG)
                    ALLEEG(cond_i).CAT.Conn = stat_getDistribMean(ALLEEG(cond_i).CAT.PConn);
                end
            end 
            for cond_i = 1:length(ALLEEG)
                %- clear bootstrap calculation
                ALLEEG(cond_i).CAT.PConn = [];
            end
            fprintf('done.\n');
        end
        %% ===================================================================== %%
        %## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION    
        % see. stat_surrogateGen
        % see. stat_surrogateStats
        if DO_PHASE_RND
            %## (1) ALTERNATIVE CODE 
            for cond_i=1:length(ALLEEG)
                %- Generate Phase Randomized Distribution
                fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
                %- clear PConn
                ALLEEG(cond_i).CAT.PConn  = [];
                fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
                chk = ~exist([ALLEEG(cond_i).filepath filesep fName,'_BootStrap.mat'],'file') || FORCE_PHASECALC
                if chk
                    stat_pr_cfg.modelingApproach = ALLEEG(cond_i).CAT.configs.est_fitMVAR;
                    stat_pr_cfg.connectivityModeling = ALLEEG(cond_i).CAT.configs.est_mvarConnectivity;
                    %- FEVAL
                    [PConn,~] = feval(@stat_surrogateGen,'ALLEEG',ALLEEG(cond_i),stat_pr_cfg);
                    ALLEEG(cond_i).CAT.PConn = PConn;
                    %- Save Phase randomized distribution
                    phasernd_dist = ALLEEG(cond_i).CAT.PConn;
                    
                    par_save(phasernd_dist,ALLEEG(cond_i).filepath,fName,'_PhaseRnd');
                    fprintf('done.\n')
                else
                    ALLEEG(cond_i).CAT.PConn = par_load(ALLEEG(cond_i).filepath,fName,'_PhaseRnd.mat');
                end
            end
        end
    catch e
        rmv_subj(subj_i) = 1;
        EEG.CAT = struct([]);
        tmp{subj_i} = EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
        disp(e);
    end
end
pop_editoptions('option_computeica',0);
%% SAVE BIG STUDY
% fprintf('==== Reformatting Study ====\n');
% %- remove bugged out subjects
% tmp = tmp(~cellfun(@isempty,tmp));
% %## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
% fss = cell(1,length(tmp));
% for subj_i = 1:length(tmp)
%     fss{subj_i} = fields(tmp{subj_i});
% end
% fss = unique([fss{:}]);
% fsPrev = fss;
% for subj_i = 1:length(tmp)
%     EEG = tmp{subj_i};
%     fs = fields(EEG);
%     % delete fields not present in other structs.
%     out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
%     out = [out{:}];
%     addFs = fs(~out);
%     if any(~out)
%         for j = 1:length(addFs)
%             EEG.(addFs{j}) = [];
%             fprintf('%s) Adding %s %s\n',EEG.subject,addFs{j})
%         end
%     end 
%     tmp{subj_i} = EEG;
% end
% tmp = cellfun(@(x) [[]; x], tmp);
% %##
% tmp = eeg_checkset(tmp,'eventconsistency');
% [STUDY, ALLEEG] = std_editset([],tmp,...
%                                 'updatedat','off',...
%                                 'savedat','off',...
%                                 'name',study_fname_1,...
%                                 'filename',study_fname_1,...
%                                 'filepath',study_load_dir);
% %## ASSIGN PARAMETERS
% STUDY.etc.a_epoch_process.epoch_chars = MAIN_STUDY.etc.a_epoch_process;
% STUDY.etc.d_cnctanl_process.epoch_chars = MAIN_STUDY.etc.d_cnctanl_process;
% STUDY.etc.e_cnctanl_surrogates.params = comps_out;
% STUDY.etc.e_cnctanl_surrogates.params = struct();
% [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
%                                             study_fname_1,study_load_dir,...
%                                             'STUDY_COND',[]);
%%
%{
SUBJ_PICS = {MAIN_STUDY.datasetinfo.subject};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
cluster_struct = STUDY.urcluster;
cluster_struct_orig = STUDY.urcluster;
% subj_chars_orig = SUBJ_PICS{1};
% subj_chars_orig = cellfun(@(x) [{} sprintf('Pilot%s',x)],subj_chars_orig);
subj_chars_orig = SUBJ_PICS;
subj_chars = {STUDY.datasetinfo.subject};
subj_keep = zeros(length(subj_chars),1);
for subj_i = 1:length(subj_chars_orig)
    disp(any(strcmp(subj_chars_orig{subj_i},subj_chars)));
    if any(strcmp(subj_chars_orig{subj_i},subj_chars))
        subj_keep(subj_i) = 1;
    end
end
subjs_rmv = find(~subj_keep);
subj_keep = find(subj_keep);
[val,ord] = sort(subj_keep);
for cli = 2:length(cluster_struct)
    si = cluster_struct(cli).sets;
    ci = cluster_struct(cli).comps;
    keep_si = setdiff(si,subjs_rmv);
    tmp = cluster_struct(cli).preclust.preclustdata;
    tmp_preclust = [];
    tmp_si = [];
    tmp_ci = [];
    for i = 1:length(keep_si)
        tmp_si = [tmp_si, repmat(ord(keep_si(i) == val),1,sum(keep_si(i) == si))];
        tmp_ci = [tmp_ci, ci(keep_si(i) == si)];
        
        tmp_preclust = [tmp_preclust; tmp(keep_si(i) == si,:)];
    end
    cluster_struct(cli).sets = tmp_si;
    cluster_struct(cli).comps = tmp_ci;
    cluster_struct(cli).preclust.preclustdata = tmp_preclust;
end
cluster_struct(1).sets = [cluster_struct(2:end).sets];
cluster_struct(1).comps = [cluster_struct(2:end).comps];
%-
STUDY.cluster = cluster_struct;
STUDY.etc.a_epoch_process = MAIN_STUDY.etc.a_epoch_process;
STUDY.etc.d_cnctanl_process.params = comps_out;
STUDY.etc.d_cnctanl_process.params = struct('CONN_METHODS',CONN_METHODS,...
            'MORDER',MORDER,...
            'DO_PHASE_RND',DO_PHASE_RND,...
            'DO_BOOTSTRAP',DO_BOOTSTRAP,...
            'FREQS',CONN_FREQS,...
            'WINDOW_LENGTH',WINDOW_LENGTH,... 
            'WINDOW_STEP_SIZE',WINDOW_STEP_SIZE,...
            'GUI_MODE','nogui',...
            'VERBOSITY_LEVEL',1,...
            'ESTSELMOD_CFG',[],...
            'FREQ_BANDS',FREQ_BANDS);

[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                        STUDY.filename,STUDY.filepath,...
                        'RESAVE_DATASETS','off');
%}
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

