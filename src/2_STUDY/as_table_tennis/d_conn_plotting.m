%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previousds Version: n/a
%   Summary: s

%- sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/AS/run_d_conn_plotting.sh
 
%{
%## RESTORE MATLAB
% WARNING: restores defdault pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
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
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'AS'];
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
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
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
%% (PARAMETERS) ======================================================== %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'AS_dataset';
% COND_CHARS = {'2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
% EVENT_CHARS = {'Subject_receive'}; %, 'Subject_receive'};
%- rally serve analysis
% COND_CHARS =  {'1Bounce_Human','Serve_Human'};
% EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
%- rally, serve, return analysis
COND_CHARS =  {'2Bounce_Human','1Bounce_Human','Serve_Human'};
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
% dt = '08182023_bounces_1h_2h_2bm_JS';
% dt = '12182023_bounces_1h_2h_2bm_JS_0p25-1';
% dt = '12282023_bounces_1h_2bm_JS_n1-0p5';
% dt = '01182023_subjrec_2bounces_1h_2bm_JS_n5-1p5';
% dt = '01252023_subjrec_2bounces_rally_serve_human_JS_n5-1p5';
% dt = '01292023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75';
dt = '01312023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75';
%## soft define
%- combinations of events and conditions
EVENT_COND_COMBOS = cell(length(COND_CHARS)*length(EVENT_CHARS),1);
cnt = 1;
for cond_i = 1:length(COND_CHARS)
    for event_i = 1:length(EVENT_CHARS)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',COND_CHARS{cond_i},EVENT_CHARS{event_i});
        cnt = cnt + 1;
    end
end
%- path for local data
% study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
study_fName_1 = 'epoch_study';
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'conn'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(EVENT_COND_COMBOS)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
    exit();
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
%     EVENT_COND_COMBOS = MAIN_STUDY.etc.a_epoch_process.epoch_chars;
    %## CUT OUT NON VALID CLUSTERS
    valid_cls = [3,4,5,6,7,8,9,11,12];
    inds = setdiff(1:length(comps_out),valid_cls);
    comps_out(inds,:) = 0;
%     MAIN_ALLEEG(1).etc.conn_meta.comps_out;
end
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
STUDY.etc.a_epoch_process.epoch_params = EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_chars = EVENT_COND_COMBOS;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
%}
%% ===================================================================== %%
% CLUSTER_ITERS = [3,7,5,4];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
CLUSTER_ASSIGNMENTS = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','LTemp','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.

%% 
%- plotting vars
FREQ_BANDS = [];
FREQ_BANDS.theta = (4:8);
FREQ_BANDS.alpha = (8:14);
FREQ_BANDS.beta = (14:30);
FREQ_BANDS.gamma = (30:50);
%-
DO_PLOT = true;
CLIM = [0,0.005];
% CLIM = [-0.005,0.005];
ALPHA = 0.05;
PLOT_CI = true;
FREQSCALE = 'log';
CONN_MEAS = 'dDTF08';
% CONN_MEAS = 'S';
%- initiate vars
fPaths = {MAIN_ALLEEG.filepath};
fNames = {MAIN_ALLEEG.filename};
cl_pairs = [];
for cl_i = main_cl_inds(2:end)
    for cl_j = main_cl_inds(2:end)
        if any(all(ismember(cl_pairs,[cl_i,cl_j]),2)) || cl_i==cl_j
            continue;
        end
        cl_pairs = [cl_pairs; cl_i,cl_j];
    end
end
%-
save_dir = [save_dir filesep CONN_MEAS];
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
save_dir_txf =  [save_dir filesep 'bootstrap_txf'];
if ~exist(save_dir_txf,'dir')
    mkdir(save_dir_txf)
end
%## LOAD TEMP EEG DATA
EEG = MAIN_ALLEEG(1); %pop_loadset('filepath',fPaths{1},'filename',fNames{1});
tmp_cat = EEG.etc.COND_CAT(1);
COND_N = length(EEG.etc.COND_CAT);
T_DIM = size(tmp_cat.Conn.(CONN_MEAS),4);
F_DIM = size(tmp_cat.Conn.(CONN_MEAS),3);
CLUSTER_PICKS = main_cl_inds(2:end);
COND_NAMES = EVENT_COND_COMBOS; %MAIN_STUDY.etc.a_epoch_process.epoch_chars; %unique(EEG.etc.conn_table.t_fNames);
%- bootsrap mean struct
BOOSTRAP_STRUCT = [];
BOOSTRAP_STRUCT.cluster_pairs = cl_pairs;
BOOSTRAP_STRUCT.condition_n = COND_N;
BOOSTRAP_STRUCT.t_dim = T_DIM;
BOOSTRAP_STRUCT.f_dim = F_DIM;
BOOSTRAP_STRUCT.subj_n = length(MAIN_ALLEEG);
BOOSTRAP_STRUCT.cluster_n = length(CLUSTER_PICKS);
BOOSTRAP_STRUCT.bootstrap_mat = [];
%- bootstrap mean mat
BOOTSTRAP_CELL = cell(length(MAIN_ALLEEG));
%%
%{
%## SUBJECT LOOP
parfor (subj_i = 1:length(MAIN_ALLEEG),length(MAIN_ALLEEG))
% for subj_i = 1:length(MAIN_ALLEEG)
    conn_mat_temp = [];
    BOOTSTRAP_MAT = zeros(COND_N,size(comps_out,1),size(comps_out,1),F_DIM,T_DIM);
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
%     clusters = find(comps > 0);
%     clusters = clusters(idxcl(idxcl ~= 0));
    display_names = cell(length(comps),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(~cellfun(@isempty,display_names));
    fprintf('%s) Cluster Names:\n',MAIN_ALLEEG(subj_i).subject);
    fprintf('\t%s\n',display_names{:});
    %## CONDITION LOOP
    for cond_i = 1:COND_N
        %## LOAD EEG DATA
        EEG = MAIN_ALLEEG(subj_i); %pop_loadset('filepath',fPaths{set_i},'filename',fNames{set_i});
        %- load bootsrapped dist
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        chk = strsplit(fName,'.');
        tmp_bs = par_load(fPath,[chk{1}, '_BootStrap.mat']);
        bs_mean = mean(tmp_bs.(CONN_MEAS),5);
        % bs_std = std(tmp_bs.(CONN_MEAS),[],5);
        BOOTSTRAP_MAT(cond_i,idxcl,idxcl,:,:,subj_i) = bs_mean;
        if DO_PLOT
            %- get cond name
            out1 = strsplit(COND_NAMES{cond_i},'.');
            out1 = strsplit(out1{1},'_');
            out1 = strjoin(out1(3:end),' ');
            %- load alleeg
            if ~ispc
                ALLEEG = pop_loadset('filepath',convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath),'filename',EEG.etc.cond_files(cond_i).fName);
            else
                ALLEEG = pop_loadset('filepath',convertPath2Drive(EEG.etc.cond_files(cond_i).fPath),'filename',EEG.etc.cond_files(cond_i).fName);
            end
            ALLEEG.CAT = EEG.etc.COND_CAT(1);
            temp_conn = ALLEEG.CAT.Conn;
            temp_conn.(CONN_MEAS) = bs_mean;
            
            [~,~,~,new_conn] = jsedit_vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',temp_conn,...
                'plotCondDiff',{},...
                'stats',{},...
                'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
                'msubset','all',...
                'MatrixLayout',{'Full','estimator',CONN_MEAS,'clim',CLIM},...
                'thresholding',{'Statistics','plotci',PLOT_CI,'sigthreshmethod','pval','alpha',ALPHA},...
                'transform','linear',...
                'freqscale',FREQSCALE,... 
                'NodeLabels',display_names,...
                'events',{{0,'r',':',2}},...
                'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
                'FrequencyMarkerColor',[0,0,0],...
                'backgroundColor',[1,1,1],...
                'textColor',[0,0,0],...
                'linecolor',[0,0,0],...
                'patchcolor',[0,0,0],...
                'axesFontSize',10,...
                'topoplot','Topoplot'); %,'estimator',CONN_METHODS
            %- plot edits
%             colormap(linspecer)
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
            if ~exist([save_dir_txf filesep sprintf('%s',out1)],'dir')
                mkdir([save_dir_txf filesep sprintf('%s',out1)])
            end
            exportgraphics(fig_i,[save_dir_txf filesep sprintf('%s',out1) filesep sprintf('%s_txf_grid.jpg',EEG.subject)]);
            %## TEST TEST
            %{
            ind_1 = randi(size(tmp_bs.(CONN_MEAS),1));
            ind_2 = randi(size(tmp_bs.(CONN_MEAS),1));
            conn_mat_temp(1,:,:) = squeeze(bs_mean(ind_2,ind_1,:,:));
            conn_mat_temp(2,:,:) = squeeze(bs_mean(ind_1,ind_2,:,:));
            times_bs = tmp_bs.erWinCenterTimes;
            freqs_bs = tmp_bs.freqs;
            vis_TimeFreqCell('ConnMatrix',conn_mat_temp,...
                'alltimes',times_bs,'allfreqs',freqs_bs,...
                'elocs',ALLEEG.chanlocs,...
                'chaninfo',ALLEEG.chaninfo,...
                'topoplot','topoplot',...
                'freqscale',FREQSCALE,...
                'clim',CLIM,...
                'topovec',squeeze(ALLEEG.icawinv(:,ALLEEG.CAT.curComps([ind_1 ind_2])))',...
                'nodelabels',{display_names{ind_1},display_names{ind_2}});
            colormap(linspecer)
            fig_i = get(groot,'CurrentFigure');
            exportgraphics(fig_i,[save_dir_txf filesep sprintf('%s',out1) filesep sprintf('%s_%s_%s_txf_grid.jpg',display_names{ind_1},display_names{ind_2},EEG.subject)]);
            %}
        end
    end
    BOOTSTRAP_CELL{subj_i} = BOOTSTRAP_MAT;
end
BOOTSTRAP_CELL = cat(6,BOOTSTRAP_CELL{:});
BOOSTRAP_STRUCT.bootstrap_mat = BOOTSTRAP_CELL;
par_save(BOOSTRAP_STRUCT,save_dir_txf,sprintf('bootstrap_mean_mat.mat'));
%}
%% ===================================================================== %%
%{
%## BOOTSTRAP ANALYSIS
%- param changes
CLIM = [-0.005,0.005];
%- save directory
save_dir_bootmats =  [save_dir filesep 'bootstrap_mat_files'];
if ~exist(save_dir_bootmats,'dir')
    mkdir(save_dir_bootmats)
end
% for subj_i = 1:length(MAIN_ALLEEG)
parfor (subj_i = 1:length(MAIN_ALLEEG),length(MAIN_ALLEEG))
% for subj_i = 1;
    %## INITATE STRUCT
    freq_inds = [];
    BOOTSTRAP_SUBJ_STRUCT = [];
    BOOTSTRAP_SUBJ_STRUCT.stats = cell(COND_N,COND_N);
    BOOTSTRAP_SUBJ_STRUCT.conditions = cell(COND_N,COND_N);
    BOOTSTRAP_SUBJ_STRUCT.averages = cell(COND_N,COND_N);
    BOOTSTRAP_SUBJ_STRUCT.masked_conn = cell(COND_N,COND_N);
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF = [];
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.FREQ_BANDS = FREQ_BANDS;
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.integrated_cause = cell(COND_N,COND_N);
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.sumed_cause = cell(COND_N,COND_N);
    %## GRAB EEG DATA
    EEG = MAIN_ALLEEG(subj_i);
    ALLEEG = cell(length(EEG.etc.COND_CAT),1);
    bootstrap_pconn = cell(length(EEG.etc.COND_CAT),1);
    bootstrap_maskedconn = cell(length(EEG.etc.COND_CAT),length(EEG.etc.COND_CAT));
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    display_names = cell(length(comps),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(~cellfun(@isempty,display_names));
    fprintf('%s) Cluster Names:\n',MAIN_ALLEEG(subj_i).subject);
    fprintf('\t%s\n',display_names{:});
    % parfor cond_i = 1:length(EEG.etc.cond_files)
    for cond_i = 1:length(EEG.etc.cond_files)
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
        ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
        %- Bootstrap Test Data Handler
        fprintf('\n==== LOADING BOOTSTRAPPED CONNECTIVITY MEASURES ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') 
            error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
        else
            bootstrap_pconn{cond_i} = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
        end
        fprintf('done.\n')
        %- condition override
        out1 = strsplit(COND_NAMES{cond_i},'.');
        out1 = strsplit(out1{1},'_');
        out1 = strjoin(out1(3:end),'_');
        ALLEEG{cond_i}.condition = out1;
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    % Note this function will return a new EEG dataset with the condition
    % differences (Set A - Set B) in the order specified in datasetOrder
    fprintf('\n===================================================\n');
    disp('Between Condition Test')
    fprintf('===================================================\n');
    boot_avgs = cell(COND_N,1);
    for cond_i = 1:length(ALLEEG)
        ALLEEG(cond_i).CAT.Stats = [];
        ALLEEG(cond_i).CAT.PConn = bootstrap_pconn{cond_i};
        bs_mean = mean(bootstrap_pconn{cond_i}.(CONN_MEAS),5);
        boot_avgs{cond_i} = bs_mean;
    end
    boot_stats = cell(COND_N,COND_N);
    cond_tests = cell(COND_N,COND_N);
    masked_conns = cell(COND_N,COND_N);
    cond_pairs = [];
    for cond_i = 1:length(ALLEEG)
        for cond_j = 1:length(ALLEEG)
            if cond_i ~= cond_j
                %## save
                %- get cond name
                out1 = strsplit(COND_NAMES{cond_i},'.');
                out1 = strsplit(out1{1},'_');
                out1 = strjoin(out1(3:end),' ');
                %- get cond name
                out2 = strsplit(COND_NAMES{cond_j},'.');
                out2 = strsplit(out2{1},'_');
                out2 = strjoin(out2(3:end),' ');
                if ~exist([save_dir_bootmats filesep sprintf('%s-%s',out1,out2)],'dir')
                    mkdir([save_dir_bootmats filesep sprintf('%s-%s',out1,out2)])
                end
    %             if any(all(ismember(cond_pairs,[cond_i,cond_j]),2)) %|| cond_i==cond_j
    %                 continue;
    %             end
    %             cond_pairs = [cond_pairs; cond_i,cond_j];
    %             ALLEEG(cond_i).condition = out1;
    %             ALLEEG(cond_j).condition = out2;
                TMP_ALLEEG = [ALLEEG(cond_i);ALLEEG(cond_j)];

                %- Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
                % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
                [tmp_stats,boot_avg,~] = stat_surrogateStats('ALLEEG',TMP_ALLEEG,...
                    'statTest',{'Hab',...
                        'tail','both',... %[left|right|one|both]
                        'testMethod','quantile',...
                        'computeci',true,... 
                        'alpha',ALPHA,...
                        'mcorrection','fdr',... %[fdr|bonferonni|numvars]
                        'statcondargs',{'mode','perm'}},...
                    'connmethods',{CONN_MEAS},...
                    'VerbosityLevel',1);
                % (08/03/2023) JS, changing mcorrection from fdr to bonferoni to be
                % a little less agressive on removing false positives.
                % (08/03/2023) JS, changing back to fdr.
                boot_stats{cond_i,cond_j} = tmp_stats;
                cond_tests{cond_i,cond_j} = {ALLEEG(cond_i).condition,ALLEEG(cond_j).condition};
    %             boot_avgs{cond_i,cond_j} = boot_avg;
                condition_order = {ALLEEG(cond_i).condition,ALLEEG(cond_j).condition};
                plot_boot = boot_avg;
                if strcmp(CONN_MEAS,{'S'})
                    plot_boot.(CONN_MEAS) = real(boot_avg.(CONN_MEAS));
                    for i = 1:length(TMP_ALLEEG)
                        TMP_ALLEEG(i).CAT.Conn.(CONN_MEAS) = real(TMP_ALLEEG(i).CAT.Conn.(CONN_MEAS));
                    end
                end
                [~,~,~,new_conn] = jsedit_vis_TimeFreqGrid('ALLEEG',TMP_ALLEEG,'Conn',plot_boot,...
                    'plotCondDiff',{'condOrder',condition_order},...
                    'stats',tmp_stats,...
                    'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
                    'msubset','all',...
                    'MatrixLayout',{'Full','estimator',CONN_MEAS,'clim',CLIM},...
                    'thresholding',{'Statistics','plotci',PLOT_CI,'sigthreshmethod','pval','alpha',ALPHA},...
                    'transform','linear',...
                    'freqscale',FREQSCALE,... 
                    'NodeLabels',display_names,...
                    'events',{{0,'r',':',2}},...
                    'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
                    'FrequencyMarkerColor',[0,0,0],...
                    'backgroundColor',[1,1,1],...
                    'textColor',[0,0,0],...
                    'linecolor',[0,0,0],...
                    'patchcolor',[0,0,0],...
                    'axesFontSize',10,...
                    'topoplot','Topoplot'); %,'estimator',CONN_METHODS
                %- plot edits
                fig_i = get(groot,'CurrentFigure');
                set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
                masked_conns{cond_i,cond_j} = new_conn;%- save
                exportgraphics(fig_i,[save_dir_bootmats filesep sprintf('%s-%s',out1,out2) filesep sprintf('%s_masked_boot_grid.jpg',EEG.subject)]);
                
                %## PLOT TIME INTEGRATD TRACES FOR CHOICE FREQUENCY BANDS
                noms = fieldnames(FREQ_BANDS);
        %         new_conn = new_conn; %ALLEEG(cond_i).CAT.Conn.(CONN_MEAS);
                freq_temp_1  = zeros(size(new_conn,1),size(new_conn,2),size(new_conn,4),length(noms));
                freq_temp_2 = zeros(size(new_conn,1),size(new_conn,2),size(new_conn,4),length(noms));
                cond_pairs = [];
                for i = 1:size(new_conn,1)
                    for j = 1:size(new_conn,2)
                        if any(all(ismember(cond_pairs,[i,j]),2)) || i==j
                            continue;
                        end
                        cond_pairs = [cond_pairs; i,j];
%                         ff = figure();
%                         title(sprintf('%s) %s-%i to %s-%i',CONN_MEAS,display_names{i},display_names{j},i,j));
%                         hold on;
%                         for nom = 1:length(noms)
%                             freq_inds = FREQ_BANDS.(noms{nom});
%                             y_in = squeeze(new_conn(i,j,freq_inds,:));
%                             x_in = 1:size(y_in,1);
%                             int_causality = trapz(x_in,y_in,1);
%                             freq_temp_1(i,j,:,nom) = int_causality;
%                             plot(ALLEEG(cond_i).CAT.Conn.erWinCenterTimes,int_causality,'DisplayName',sprintf('%s trapz',noms{nom}));
%                             int_causality = sum(y_in,1);
%                             freq_temp_2(i,j,:,nom) = int_causality;
%                             plot(ALLEEG(cond_i).CAT.Conn.erWinCenterTimes,int_causality,'DisplayName',sprintf('%s sum',noms{nom}));
%                         end
%                         xlabel('time (s)');
%                         ylabel('Causality');
%                         legend();
%                         hold off
%                         exportgraphics(ff,[save_dir_bootmats filesep sprintf('%s-%s',out1,out2) filesep,...
%                             sprintf('%s_%s-%s_freq%i-%i_integrated_caus.jpg',EEG.subject,display_names{i},display_names{j},freq_inds(1),freq_inds(end))]);
                    end
                end
                close('all');
                BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.integrated_cause{cond_i,cond_j} = freq_temp_1;
                BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.sumed_cause{cond_i,cond_j} = freq_temp_2;
            end
        end
    end
    BOOTSTRAP_SUBJ_STRUCT.stats = boot_stats;
    BOOTSTRAP_SUBJ_STRUCT.conditions = cond_tests;
    BOOTSTRAP_SUBJ_STRUCT.averages = boot_avgs;
    BOOTSTRAP_SUBJ_STRUCT.masked_conn = masked_conns;
    BOOTSTRAP_SUBJ_STRUCT.subject = EEG.subject;
    BOOTSTRAP_SUBJ_STRUCT.conn_meas = CONN_MEAS;
    par_save(BOOTSTRAP_SUBJ_STRUCT,save_dir_bootmats,sprintf('%s_boot_struct.mat',EEG.subject));
end
%}
%% ===================================================================== %%
%## BASELINE BOOTSTRAP ANALYSIS
%- param changes
CLIM = [-0.005,0.005];
BASELINE_TIME = [-0.5,0];
%- save directory
save_dir_bootmats =  [save_dir filesep 'bootstrap_baseline_mat_files'];
if ~exist(save_dir_bootmats,'dir')
    mkdir(save_dir_bootmats)
end
% for subj_i = 1:length(MAIN_ALLEEG)
parfor (subj_i = 1:length(MAIN_ALLEEG),length(MAIN_ALLEEG))
    %## INITATE STRUCT
    freq_inds = [];
    BOOTSTRAP_SUBJ_STRUCT = [];
    BOOTSTRAP_SUBJ_STRUCT.stats = cell(COND_N,1);
    BOOTSTRAP_SUBJ_STRUCT.conditions = cell(COND_N,1);
    BOOTSTRAP_SUBJ_STRUCT.averages = cell(COND_N,1);
    BOOTSTRAP_SUBJ_STRUCT.masked_conn = cell(COND_N,1);
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF = [];
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.FREQ_BANDS = FREQ_BANDS;
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.integrated_cause = cell(COND_N,1);
    BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.sumed_cause = cell(COND_N,1);
    %## GRAB EEG DATA
    EEG = MAIN_ALLEEG(subj_i);
    ALLEEG = cell(length(EEG.etc.COND_CAT),1);
    bootstrap_pconn = cell(length(EEG.etc.COND_CAT),1);
    bootstrap_maskedconn = cell(length(EEG.etc.COND_CAT),length(EEG.etc.COND_CAT));
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    display_names = cell(length(comps),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(~cellfun(@isempty,display_names));
    fprintf('%s) Cluster Names:\n',MAIN_ALLEEG(subj_i).subject);
    fprintf('\t%s\n',display_names{:});
    % parfor cond_i = 1:length(EEG.etc.cond_files)
    for cond_i = 1:length(EEG.etc.cond_files)
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
        ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
        %- Bootstrap Test Data Handler
        fprintf('\n==== LOADING BOOTSTRAPPED CONNECTIVITY MEASURES ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') 
            error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
        else
            bootstrap_pconn{cond_i} = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
        end
        fprintf('done.\n')
        %- condition override
        out1 = strsplit(COND_NAMES{cond_i},'.');
        out1 = strsplit(out1{1},'_');
        out1 = strjoin(out1(3:end),'_');
        ALLEEG{cond_i}.condition = out1;
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    %-
    fprintf('\n===================================================\n');
    disp('BASELINE Test')
    fprintf('===================================================\n');
    boot_avgs = cell(COND_N,1);
    for cond_i = 1:length(ALLEEG)
        ALLEEG(cond_i).CAT.Stats = [];
        ALLEEG(cond_i).CAT.PConn = bootstrap_pconn{cond_i};
        bs_mean = mean(bootstrap_pconn{cond_i}.(CONN_MEAS),5);
        boot_avgs{cond_i} = bs_mean;
    end
    boot_stats = cell(COND_N,1);
    cond_tests = cell(COND_N,1);
    masked_conns = cell(COND_N,1);
    cond_pairs = [];
    for cond_i = 1:length(ALLEEG)
        %## save
        %- get cond name
        out1 = strsplit(COND_NAMES{cond_i},'.');
        out1 = strsplit(out1{1},'_');
        out1 = strjoin(out1(3:end),' ');
        if ~exist([save_dir_bootmats filesep sprintf('%s',out1)],'dir')
            mkdir([save_dir_bootmats filesep sprintf('%s',out1)])
        end

        %- Statistics for each dynamical measure are now stored in EEG.CAT.Stats.
        % The dimensionality is [num_vars x num_vars x num_freqs x num_times]
        [tmp_stats,~,~] = stat_surrogateStats('ALLEEG',ALLEEG(cond_i),...
                     'statTest',{'Hbase',...
                        'baseline',BASELINE_TIME,...
                        'testMeans',true,...
                        'tail','both',... %[left|right|one|both]
                        'testMethod','quantile',...
                        'computeci',true,...
                        'alpha',ALPHA,...
                        'mcorrection','fdr',...
                        'statcondargs',{'mode','perm'}},...
                    'connmethods',{CONN_MEAS},...
                    'VerbosityLevel',1);
        % (11/08/2023) JS, initiated
        boot_stats{cond_i} = tmp_stats;
        cond_tests{cond_i} = ALLEEG(cond_i).condition;
%         condition_order = {ALLEEG(cond_i).condition,ALLEEG(cond_j).condition};
        if strcmp(CONN_MEAS,{'S'})
            ALLEEG(cond_i).CAT.Conn.(CONN_MEAS) = real(ALLEEG(cond_i).CAT.Conn.(CONN_MEAS));
        end
        [~,~,~,new_conn] = jsedit_vis_TimeFreqGrid('ALLEEG',ALLEEG(cond_i),'Conn',ALLEEG(cond_i).CAT.Conn,...
            'stats',tmp_stats,...
            'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
            'msubset','all',...
            'MatrixLayout',{'Full','estimator',CONN_MEAS,'clim',CLIM},...
            'thresholding',{'Statistics','sigthreshmethod','pval','alpha',ALPHA},...
            'transform','linear',...
            'freqscale',FREQSCALE,... 
            'NodeLabels',display_names,...
            'events',{{0,'r',':',2}},...
            'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
            'FrequencyMarkerColor',[0,0,0],...
            'backgroundColor',[1,1,1],...
            'textColor',[0,0,0],...
            'linecolor',[0,0,0],...
            'patchcolor',[0,0,0],...
            'axesFontSize',10,...
            'topoplot','Topoplot');
        masked_conns{cond_i} = new_conn;
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
        exportgraphics(fig_i,[save_dir_bootmats filesep sprintf('%s',out1) filesep sprintf('%s_%i_masked_boot_grid.jpg',EEG.subject,cond_i)]);
        
        %## PLOT TIME INTEGRATD TRACES FOR CHOICE FREQUENCY BANDS
        noms = fieldnames(FREQ_BANDS);
%         new_conn = new_conn; %ALLEEG(cond_i).CAT.Conn.(CONN_MEAS);
        freq_temp_1  = zeros(size(new_conn,1),size(new_conn,2),size(new_conn,4),length(noms));
        freq_temp_2 = zeros(size(new_conn,1),size(new_conn,2),size(new_conn,4),length(noms));
        cond_pairs = [];
        for i = 1:size(new_conn,1)
            for j = 1:size(new_conn,2)
                if any(all(ismember(cond_pairs,[i,j]),2)) || i==j
                    continue;
                end
                cond_pairs = [cond_pairs; i,j];
%                 ff = figure();
%                 title(sprintf('%s) %s-%i to %s-%i',CONN_MEAS,display_names{i},display_names{j},i,j));
%                 hold on;
%                 for nom = 1:length(noms)
%                     freq_inds = FREQ_BANDS.(noms{nom});
%                     y_in = squeeze(new_conn(i,j,freq_inds,:));
%                     x_in = 1:size(y_in,1);
%                     int_causality = trapz(x_in,y_in,1);
%                     freq_temp_1(i,j,:,nom) = int_causality;
%                     plot(ALLEEG(cond_i).CAT.Conn.erWinCenterTimes,int_causality,'DisplayName',sprintf('%s trapz',noms{nom}));
%                     int_causality = sum(y_in,1);
%                     freq_temp_2(i,j,:,nom) = int_causality;
%                     plot(ALLEEG(cond_i).CAT.Conn.erWinCenterTimes,int_causality,'DisplayName',sprintf('%s sum',noms{nom}));
%                 end
%                 xlabel('time (s)');
%                 ylabel('Causality');
%                 legend();
%                 hold off
%                 exportgraphics(ff,[save_dir_bootmats filesep sprintf('%s',out1) filesep,...
%                     sprintf('%s_%s-%s_freq%i-%i_integrated_caus.jpg',EEG.subject,display_names{i},display_names{j},freq_inds(1),freq_inds(end))]);
            end
        end
        close('all');
        BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.integrated_cause{cond_i,1} = freq_temp_1;
        BOOTSTRAP_SUBJ_STRUCT.FREQ_BAND_INF.sumed_cause{cond_i,1} = freq_temp_2;

    end
    BOOTSTRAP_SUBJ_STRUCT.stats = boot_stats;
    BOOTSTRAP_SUBJ_STRUCT.conditions = cond_tests;
    BOOTSTRAP_SUBJ_STRUCT.averages = boot_avgs;
    BOOTSTRAP_SUBJ_STRUCT.masked_conn = masked_conns;
    BOOTSTRAP_SUBJ_STRUCT.subject = EEG.subject;
    BOOTSTRAP_SUBJ_STRUCT.conn_meas = CONN_MEAS;
    par_save(BOOTSTRAP_SUBJ_STRUCT,save_dir_bootmats,sprintf('%s_boot_struct.mat',EEG.subject));
end
%}
%% ===================================================================== %%
%{
%## PHASERANDOMIZATION ANALYSIS
%- param changes
CLIM = [0,0.005];
%-
save_dir_nonzmats =  [save_dir filesep 'nonzero_mat_files'];
if ~exist(save_dir_nonzmats,'dir')
    mkdir(save_dir_nonzmats)
end
% for subj_i = 1:length(MAIN_ALLEEG)
parfor (subj_i = 1:length(MAIN_ALLEEG),length(MAIN_ALLEEG))
    %## INITATE STRUCT
    freq_inds = [];
    PHASERND_SUBJ_STRUCT = [];
    PHASERND_SUBJ_STRUCT.stats = cell(COND_N,1);
    PHASERND_SUBJ_STRUCT.conditions = cell(COND_N,1);
    PHASERND_SUBJ_STRUCT.averages = cell(COND_N,1);
    PHASERND_SUBJ_STRUCT.masked_conn = cell(COND_N,1);
    PHASERND_SUBJ_STRUCT.FREQ_BAND_INF = [];
    PHASERND_SUBJ_STRUCT.FREQ_BAND_INF.FREQ_BANDS = FREQ_BANDS;
    PHASERND_SUBJ_STRUCT.FREQ_BAND_INF.integrated_cause = cell(COND_N,1);
    PHASERND_SUBJ_STRUCT.FREQ_BAND_INF.sumed_cause = cell(COND_N,1);
    %## GRAB EEG DATA
    EEG = MAIN_ALLEEG(subj_i);
    ALLEEG = cell(length(EEG.etc.COND_CAT),1);
    phasernd_pconn = cell(length(EEG.etc.COND_CAT),1);
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    display_names = cell(length(comps),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(~cellfun(@isempty,display_names));
    fprintf('%s) Cluster Names:\n',MAIN_ALLEEG(subj_i).subject);
    fprintf('\t%s\n',display_names{:});
    % parfor cond_i = 1:length(EEG.etc.cond_files)
    for cond_i = 1:length(EEG.etc.cond_files)
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        ALLEEG{cond_i} = pop_loadset('filepath',fPath,'filename',fName);
        ALLEEG{cond_i}.CAT = EEG.etc.COND_CAT(cond_i);
        %- PhaseRnd test
        fprintf('\n==== LOADING PHASERAND STATISTICS ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file')
            error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
        else
            phasernd_pconn{cond_i} = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
        end
        fprintf('done.\n')
        %- condition override
%         ALLEEG{cond_i}.condition = sprintf('cond_%i',cond_i);
    end
    ALLEEG = cellfun(@(x) [[],x],ALLEEG);
    %## 3) Test for non-zero connectivity
    %     We are testing with respect to a phase-randomized null
    %     distribution. A p-value for rejection of the null hypothesis
    %     can be obtained by computing the probability that the
    %     observed connectivity is a random sample from the null distribution
    nonzero_stats = cell(COND_N,1);
    cond_tests = cell(COND_N,1);
    masked_conns = cell(COND_N,1);
    for cond_i = 1:length(ALLEEG)
        fprintf('\n===================================================\n');
        disp('NonZero Test')
        fprintf('===================================================\n');
        %##
        %- get cond name
        out1 = strsplit(COND_NAMES{cond_i},'.');
        out1 = strsplit(out1{1},'_');
        out1 = strjoin(out1(3:end),' ');
        if ~exist([save_dir_nonzmats filesep sprintf('%s',out1)],'dir')
            mkdir([save_dir_nonzmats filesep sprintf('%s',out1)])
        end
        %- set params
        ALLEEG(cond_i).CAT.PConn = phasernd_pconn{cond_i};
        ALLEEG(cond_i).CAT.Stats = [];
        %- compute stats
        [tmp_stats,~,~] = stat_surrogateStats('ALLEEG',ALLEEG(cond_i),...
                         'statTest',{'Hnull',...
                            'tail','one',... %[left|right|one|both]
                            'testMethod','quantile',... %'quantile','condstat'
                            'computeci',false,... %[true|false}
                            'alpha',0.95,... %0-1
                            'mcorrection','none',... % ['none','fdr','bonferroni','numvars']
                            'statcondargs',{'mode','perm'}},...
                        'connmethods',{CONN_MEAS},...
                        'VerbosityLevel',1);
        % (11/08/2023) JS, initiated
        %## PLOT
        if strcmp(CONN_MEAS,{'S'})
            ALLEEG(cond_i).CAT.Conn.(CONN_MEAS) = real(ALLEEG(cond_i).CAT.Conn.(CONN_MEAS));
        end
        [~,~,~,new_conn] = jsedit_vis_TimeFreqGrid('ALLEEG',ALLEEG(cond_i),'Conn',ALLEEG(cond_i).CAT.Conn,...
            'stats',tmp_stats,...
            'vismode','TimeXFrequency',... %'TimeXFrequency','TimeXCausality','FrequencyXCausality');
            'msubset','all',...
            'MatrixLayout',{'Full','estimator',CONN_MEAS,'clim',CLIM},...
            'thresholding',{'Statistics','sigthreshmethod','pval','alpha',ALPHA},...
            'transform','linear',...
            'freqscale',FREQSCALE,... 
            'NodeLabels',display_names,...
            'events',{{0,'r',':',2}},...
            'FrequencyMarkers',[0,1.3863,1.9459,2.8904,3.3673,3.8501,4.3307],...
            'FrequencyMarkerColor',[0,0,0],...
            'backgroundColor',[1,1,1],...
            'textColor',[0,0,0],...
            'linecolor',[0,0,0],...
            'patchcolor',[0,0,0],...
            'axesFontSize',10,...
            'topoplot','Topoplot');
        
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[0.05,0.3,0.7,0.7]);
        exportgraphics(fig_i,[save_dir_nonzmats filesep sprintf('%s',out1) filesep sprintf('%s_%i_masked_boot_grid.jpg',EEG.subject,cond_i)]);
        %## PLOT TIME INTEGRATD TRACES FOR CHOICE FREQUENCY BANDS
        noms = fieldnames(FREQ_BANDS);
        %## CLEANUP
        masked_conns{cond_i} = new_conn;
        nonzero_stats{cond_i} = tmp_stats;
        cond_tests{cond_i} = ALLEEG(cond_i).condition;
        close('all');
    end
    PHASERND_SUBJ_STRUCT.stats = nonzero_stats;
    PHASERND_SUBJ_STRUCT.conditions = cond_tests;
    PHASERND_SUBJ_STRUCT.masked_conn = masked_conns;
    PHASERND_SUBJ_STRUCT.subject = EEG.subject;
    PHASERND_SUBJ_STRUCT.conn_meas = CONN_MEAS;
    par_save(PHASERND_SUBJ_STRUCT,save_dir_nonzmats,sprintf('%s_phasernd_struct.mat',EEG.subject));
end
%}
%%
for i = 1:length(MAIN_ALLEEG)
    disp(sum(strcmp('Serve_Human',{MAIN_ALLEEG(i).event.bounces})))
end
%%
% for cond_i = COND_N
%     
%     %## BOOTSTRAP TEST PLOT
%     cl_pairs = [];
%     cnt = 1;
%     save_dir_cond = [save_dir_txf filesep sprintf('%i',cond_i)];
%     if ~exist(save_dir_cond,'dir')
%         mkdir(save_dir_cond)
%     end
%     for cl_i = CLUSTER_PICKS
%         for cl_j = CLUSTER_PICKS
%             if any(all(ismember(cl_pairs,[cl_i,cl_j]),2)) || cl_i==cl_j
%                 continue;
%             end
%             fprintf('Plotting sets for connection %i to %i\n',cl_i,cl_j);
%             comp_i = comps_out(cl_i,:);
%             comp_j = comps_out(cl_j,:);
%             sets_i = find(comp_i ~= 0);
%             sets_j = find(comp_j ~= 0);
%             set_vals = intersect(sets_i,sets_j);
%             bootstrapped_out.conn_set_num = [bootstrapped_out.conn_set_num; length(set_vals)];
%             connmat_out = cell(length(set_vals),1);
%             conn_mat_temp = [];
%             for i = 1:lenght(set_vals)
%                 set_i = set_vals(i);
%                 cl_1 = comp_i(set_i);
%                 cl_2 = comp_j(set_i); 
%                 tmp = comps_out(:,set_i);
%                 tmp = sort(tmp(tmp ~= 0));
%                 comp_1 = find(tmp == cl_1);
%                 comp_2 = find(tmp == cl_2);
%                 fprintf('loading subject %s\n',MAIN_ALLEEG(set_i).subject)
%                 %##
%                  EEG = MAIN_ALLEEG(subj_i); %pop_loadset('filepath',fPaths{set_i},'filename',fNames{set_i});
%                 %- load bootsrapped dist
%                 if ispc
%                     fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
%                 else
%                     fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
%                 end
%                 if ~ispc
%                     ALLEEG = pop_loadset('filepath',convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath),'filename',EEG.etc.cond_files(cond_i).fName);
%                 else
%                     ALLEEG = pop_loadset('filepath',convertPath2Drive(EEG.etc.cond_files(cond_i).fPath),'filename',EEG.etc.cond_files(cond_i).fName);
%                 end
%                 ALLEEG.CAT = EEG.etc.COND_CAT(1);
%                 fName = EEG.etc.cond_files(cond_i).fName;
%                 chk = strsplit(fName,'.');
%                 tmp_bs = par_load(fPath,[chk{1}, '_BootStrap.mat']);
%                 %##
%                 times_bs = tmp_bs.erWinCenterTimes;
%                 freqs_bs = tmp_bs.freqs;
%     %                     rnd_i = randi(size(tmp_bs.(CONN_MEAS),5));
%                 bs_mean = mean(tmp_bs.(CONN_MEAS),5);
%     %                     bs_std = std(tmp_bs.(CONN_MEAS),[],5);
%                 conn_in = bs_mean;
%                 conn_mat_temp(1,:,:) = squeeze(conn_in(comp_1,comp_2,:,:));
%                 conn_mat_temp(2,:,:) = squeeze(conn_in(comp_2,comp_1,:,:));
%                 vis_TimeFreqCell('ConnMatrix',ConnMatrix,...
%                     'alltimes',times_bs,'allfreqs',freqs_bs,...
%                     'elocs',ALLEEG.chanlocs,...
%                     'chaninfo',ALLEEG.chaninfo,...
%                     'topoplot','topoplot',...
%                     'colormap',linspecer,...
%                     'freqscale','log',...
%                     'clim',[0,0.01],...
%                     'topovec',squeeze(ALLEEG.icawinv(:,ALLEEG.CAT.curComps([comp_1 comp_2])))',...
%                     'nodelabels',{num2str(cl_1),num2str(cl_2)});
%                 fig = get(groot,'CurrentFigure');
%                 exportgraphics(fig,[save_dir_cond filesep sprintf('%i_%i_%s_bootstrapped_mean_txf.jpg',cl_i,cl_j,EEG.subject)])
%             end
%             cl_pairs = [cl_pairs; cl_i,cl_j];
%             cnt = cnt + 1;
%         end
%     end
% end
%{
%%
if ~ispec
    addpath(convertPath2UNIX('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
else
    addpath(convertPath2Drive('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
end
%%
input1 = squeeze(bootstrapped_out.matrix{1}(1,:,:,:));
input2 = squeeze(bootstrapped_out.matrix{1}(2,:,:,:));
t_test_type = 1;
pval_preselc = 100; %0.05;
n_iters = 2000;
[mask,tscore,pValue] = clusterLevelPermutationTest(input1,input2,t_test_type,pval_preselc,n_iters);
%%
mask_in = pValue<0.05;
mean_mat = mean(input1,3);
contour_in = mean_mat.*mask_in;
%##
figure();
colormap(linspecer)
contourf(contour_in)
colorbar();
%}
