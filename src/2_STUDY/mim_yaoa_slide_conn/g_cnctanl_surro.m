%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_slide_conn/run_g_cnctanl_surro.sh

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
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
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
%## HIGH LEVEL VARS
%- connectivity statistics & surrogates params
%* boolean for model bootstrap distribution
DO_BOOTSTRAP = false;
FORCE_BOOTCALC = false;
ASSIGN_BOOTSTRAP_MEAN = true;
%* boolean for phase randomization generation
DO_PHASE_RND = true;
FORCE_PHASECALC = false;
%## SURROGATE STATISTICS PARAMS
%- BOOTSTRAP
DEF_STAT_BS_CFG = struct('mode',{{'Bootstrap','nperms',200,'saveTrialIdx',false}},...
                           'modelingApproach',[],...
                           'connectivityModeling',[],...
                           'verb',1);
%- PHASE RANDOMIZATION
% DEF_STAT_PR_CFG = struct('mode',struct('arg_direct',1,...
%                                         'nperms',200,...
%                                         'arg_selection','PhaseRand'),...
%                            'modelingApproach',[],...
%                            'connectivityModeling',[],...
%                            'verb',1);
DEF_STAT_PR_CFG = struct('mode',{{'PhaseRand','nperms',200}},...
                           'modelingApproach',[],...
                           'connectivityModeling',[],...
                           'autosave',[], ...
                           'verb',1);
CONN_METHODS = {'dDTF08','S'};
CONN_FREQS = (4:50);
DEF_ESTFITMVAR = struct('connmethods',{CONN_METHODS}, ...
            'absvalsq',true,           ...
            'spectraldecibels',true,   ...
            'freqs',CONN_FREQS,        ...
            'verb',1);
% (01/29/2024) changing these to 200 to save on computation time, but 2000
% iterations may be more robust.
%## DATASET & CLUSTER INFO
DATA_SET = 'MIM_dataset';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%## soft define
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load study file
STUDY_FNAME_LOAD = 'slide_conn_study';
study_fpath = [studies_fpath filesep sprintf('%s',study_dir_name)];
%- load cluster
CLUSTER_K = 11;
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
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
%## PARFOR LOOP
EEG = [];
parfor (subj_i = 1:length(LOOP_VAR),SLURM_POOL_SIZE)
% LOOP_VAR = 1;
% for subj_i = LOOP_VAR
    %## PARAMS
    pop_editoptions('option_computeica', 1); % this is needed or SIFT will bug out
    stat_pr_cfg = DEF_STAT_PR_CFG;
    stat_bs_cfg = DEF_STAT_BS_CFG;
    est_fitmvar_cfg = DEF_ESTFITMVAR;
    %## LOAD EEG DATA
    ALLEEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    %- reassign cat structure
    ALLEEG.CAT = ALLEEG.etc.COND_CAT;
    %% IF NEEDED RECALULATE CONNECTIVITY
    cfg = struct2args(est_fitmvar_cfg);
    for cond_i = 1:length(ALLEEG)
        ALLEEG(cond_i).CAT.configs.est_mvarConnectivity = [];
        % [ALLEEG(cond_i),cfg] = pop_est_mvarConnectivity(ALLEEG(cond_i),'nogui',...
        %     cfg{:});
        [conn,~] = est_mvarConnectivity('ALLEEG',ALLEEG(cond_i),...
            'MODEL',ALLEEG(cond_i).CAT.MODEL,...
            cfg{:});
        ALLEEG(cond_i).CAT.Conn = conn;
        est_fitmvar_cfg.arg_direct = 0;
        ALLEEG(cond_i).CAT.configs.est_mvarConnectivity = est_fitmvar_cfg;
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
        tt = tic();
        for cond_i=1:length(ALLEEG)
            %- clear PConn
            ALLEEG(cond_i).CAT.PConn  = [];
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            chk = ~exist([ALLEEG(cond_i).filepath filesep fName,'_BootStrap.mat'],'file') || FORCE_BOOTCALC;
            if chk

                % [PConn,~] = feval(@stat_surrogateGen,'ALLEEG',ALLEEG(cond_i),DEF_STAT_BS_CFG);
                cfg = struct2args(stat_bs_cfg);
                % [PConn,~] = stat_surrogateGen('ALLEEG',ALLEEG(cond_i),cfg{:});
                [PConn,~] = stat_surrogateGen('EEG',ALLEEG(cond_i),cfg{:});
                
                ALLEEG(cond_i).CAT.PConn = PConn;
                %- save BootStrap distribution 
                bootstrap_dist = ALLEEG(cond_i).CAT.PConn;
                fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
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
        fprintf('%s) Bootstrap Calculation Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
    end
    %% ===================================================================== %%
    %## STEP 5.b) GENERATE PHASE RANDOMIZED DISTRIBUTION    
    % see. stat_surrogateGen
    % see. stat_surrogateStats
    if DO_PHASE_RND
        fprintf(1,'\n==== GENERATING CONNECTIVITY STATISTICS FOR SUBJECT DATA ====\n');
        tt = tic();
        for cond_i=1:length(ALLEEG)
            %- Generate Phase Randomized Distribution
            fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
            %- clear PConn
            ALLEEG(cond_i).CAT.PConn  = [];
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            chk = ~exist([ALLEEG(cond_i).filepath filesep fName,'_PhaseRnd.mat'],'file') || FORCE_PHASECALC;
            if chk
                %- clear PConn
                ALLEEG(cond_i).CAT.PConn  = [];%- FEVAL
                
                %-
                % ALLEEG(cond_i) = pop_stat_surrogateGen(ALLEEG(cond_i),'nogui', ...
                %     'modelingApproach', ALLEEG(cond_i).CAT.configs.est_fitMVAR, ...
                %     'connectivityModeling',ALLEEG(cond_i).CAT.configs.est_mvarConnectivity, ...
                %     'mode',{'PhaseRand','nperms',200}, ...
                %     'autosave',[], ...
                %     'verb',1);
                %-
                % cfg = struct2args(stat_pr_cfg);
                % [ALLEEG(cond_i),~] = pop_stat_surrogateGen(ALLEEG(cond_i),'nogui',cfg{:});
                %-
                stat_pr_cfg.modelingApproach = ALLEEG(cond_i).CAT.configs.est_fitMVAR;
                stat_pr_cfg.connectivityModeling = ALLEEG(cond_i).CAT.configs.est_mvarConnectivity;
                cfg = struct2args(stat_pr_cfg);
                [PConn,~] = stat_surrogateGen('EEG',ALLEEG(cond_i),cfg{:});
                ALLEEG(cond_i).CAT.PConn = PConn;
                %- Save Phase randomized distribution
                phasernd_dist = ALLEEG(cond_i).CAT.PConn;
                fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
                par_save(phasernd_dist,ALLEEG(cond_i).filepath,fName,'_PhaseRnd');
                fprintf('done.\n')
            else
                ALLEEG(cond_i).CAT.PConn = par_load(ALLEEG(cond_i).filepath,fName,'_PhaseRnd.mat');
            end
        end
        fprintf('%s) Phase Randomization Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
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

