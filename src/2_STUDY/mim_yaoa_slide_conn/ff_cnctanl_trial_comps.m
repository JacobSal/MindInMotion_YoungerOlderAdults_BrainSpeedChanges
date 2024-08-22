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
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
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
%## HIGH LEVEL VARS
%-  recalculate params
TRIALS_PROCESS = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
%## SURROGATE STATISTICS PARAMS
CONN_FREQS = (4:50);
CONN_METHODS = {'dDTF08','S'}; %{'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
DEF_ESTFITMVAR = struct('connmethods',{CONN_METHODS}, ...
            'absvalsq',true,           ...
            'spectraldecibels',true,   ...
            'freqs',CONN_FREQS,        ...
            'verb',1);
%## SURROGATE STATISTICS PARAMS

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
% tmp_eeg_store = cell(length(LOOP_VAR));
tmp_fpath_store = cell(length(LOOP_VAR),1);
parfor subj_i = 1:length(LOOP_VAR)
% LOOP_VAR = 2;
% for subj_i = LOOP_VAR
    %## PARAMS
    pop_editoptions('option_computeica', 1); % this is needed or SIFT will bug out
    tmp_est_fitmvar = DEF_ESTFITMVAR;
    %## LOAD EEG DATA
    ALLEEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    %- reassign cat structure
    ALLEEG.CAT = ALLEEG.etc.COND_CAT;
    %% ===================================================================== %%
    %## TRIAL BY TRIAL CONNECTVITY COMPUTATION
    %- NOTE: Next we will compute various dynamical quantities, including connectivity,
    % from the fitted VAR model.
    fprintf('===========\nCONNECTIVITY ESTIMATION\n===============\n');
    tt = tic();
    cfg = struct2args(tmp_est_fitmvar);
    tmp_alleeg = cell(length(TRIALS_PROCESS),1);
    tmp_fpaths = cell(length(TRIALS_PROCESS),1);
    for cond_i = 1:length(TRIALS_PROCESS)
        fprintf('Running condition %s...\n',TRIALS_PROCESS{cond_i})
        EEG = ALLEEG;
        %## Get Trial Indicies & Extract
        inds1 = logical(strcmp({EEG.event.cond}, TRIALS_PROCESS{cond_i}));
        inds2 = logical(strcmp({EEG.event.type}, 'boundary'));
        val_inds = find(inds1 & ~inds2);
        fromt = [EEG.event(val_inds(1)).latency];
        % ind_from = find(EEG.times>=(((fromt-1)/EEG.srate)*1000) & EEG.times<=(((fromt+1)/EEG.srate)*1000));
        tot = [EEG.event(val_inds(end)).latency];
        % ind_to = find(EEG.times>=(((tot-1)/EEG.srate)*1000) & EEG.times<=(((tot+1)/EEG.srate)*1000));
        fprintf('%s) %s length is %0.2fs\n',EEG.subject,TRIALS_PROCESS{cond_i},(tot/EEG.srate)-(fromt/EEG.srate));
        EEG = pop_select(EEG, 'point', [fromt, tot]);
        %- set cat data
        time_crop = EEG.CAT.times >= (fromt/EEG.srate)*1000 & EEG.CAT.times <= (tot/EEG.srate)*1000;
        EEG.CAT.srcdata = EEG.icaact(double(string(EEG.CAT.curComponentNames)),:); %EEG.CAT.srcdata(:,time_crop);
        EEG.CAT.times = EEG.times; %EEG.CAT.times(:,time_crop);
        % EEG.CAT.Conn = []; %?
        %- set winstart times of mvar
        err = 5.5;
        conn_inds = EEG.CAT.MODEL.winStartTimes >= ((fromt)/EEG.srate)-err & EEG.CAT.MODEL.winStartTimes <= ((tot)/EEG.srate)+err;
        EEG.CAT.MODEL.winStartTimes = EEG.CAT.MODEL.winStartTimes(conn_inds);
        EEG.CAT.MODEL.AR = EEG.CAT.MODEL.AR(conn_inds);
        EEG.CAT.MODEL.PE = EEG.CAT.MODEL.PE(conn_inds);
        EEG.CAT.MODEL.RC = EEG.CAT.MODEL.RC(conn_inds);
        EEG.CAT.MODEL.mu = EEG.CAT.MODEL.mu(conn_inds);
        EEG.CAT.MODEL.th = EEG.CAT.MODEL.th(conn_inds);
        %##
        Conn = est_mvarConnectivity('ALLEEG',EEG,'MODEL',EEG.CAT.MODEL,...
                    cfg{:});
        if ~isempty(Conn)
            EEG.CAT.Conn = Conn; 
        end
        % clear any existing visualization GUI config files
        visFields = fieldnames(EEG.CAT.configs);
        visFields = visFields(~cellfun(@isempty,strfind(visFields,'vis_')));
        for k=1:length(visFields)
            EEG.CAT.configs.(visFields{k}) = struct([]);
        end
        
        if ~isempty(cfg)
            % store the configuration structure
            EEG.CAT.configs.('est_mvarConnectivity') = cfg;
        end
        tmp_alleeg{cond_i} = EEG;
        tmp = strsplit(fNames{subj_i},'.');
        tmp{1} = [tmp{1},sprintf('_%s',TRIALS_PROCESS{cond_i})];
        tmp = strjoin(tmp,'.');
        pop_saveset(EEG,'filepath',fPaths{subj_i},'filename',tmp);
        tmp_fpaths{cond_i} = [fPaths{subj_i} filesep tmp];
    end
    tmp_fpath_store{subj_i} = tmp_fpaths;
    % tmp_eeg_store{subj_i} = cat(1,tmp_alleeg{:});
    fprintf('%s) Connectivity Estimation  Done: %0.1f\n\n',ALLEEG(1).subject,toc(tt));
end
par_save(tmp_fpath_store,study_fpath,'conn_slide_condition_fpaths.mat');
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