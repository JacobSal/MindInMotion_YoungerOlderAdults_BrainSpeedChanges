%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: s

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/AS/run_d_conn_plotting.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
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
% COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
% EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
COND_CHARS =  {'2Bounce_Human','2Bounce_BM'};
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
% dt = '08182023_bounces_1h_2h_2bm_JS';
% dt = '12182023_bounces_1h_2h_2bm_JS_0p25-1';
dt = '12282023_bounces_1h_2bm_JS_n1-0p5';
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
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
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
else
    if ~ispc
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(MAIN_STUDY);
end
%%
% CLUSTER_ITERS = [3,7,5,4];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12];
CLUSTER_ASSIGNMENTS = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
%% 
%- plotting vars
FREQ_BANDS = [];
FREQ_BANDS.theta = (4:8);
FREQ_BANDS.alpha = (8:14);
FREQ_BANDS.beta = (14:30);
FREQ_BANDS.gamma = (30:50);
FREQ_BANDS.broad_4to100 = (4:100);
% fThetaInds=find(EEG.CAT.Conn.freqs>=4 & EEG.CAT.Conn.freqs<=8);
% fAlphaInds=find(EEG.CAT.Conn.freqs>=8 & EEG.CAT.Conn.freqs<=13);
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
FREQS_VEC = tmp_cat.Conn.freqs;
CLUSTER_PICKS = main_cl_inds(2:end);
tmp = unique(EEG.etc.conn_table.t_fNames);
COND_NAMES = cell(length(tmp),1);
for i = 1:length(tmp)
    out = strsplit(tmp{i},'.');
    out = strsplit(out{1},'_');
    out = strjoin(out(3:end),' ');
    COND_NAMES{i} = out;
end
%% MODEL ORDER & VALIDATION DATA
[tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(save_dir,COND_CHARS);
fprintf('HQ minimum information crit model order: %0.2f\n',mean(tbl_summary_out.min_modorder_info_crit_hq_line));
fprintf('AIC minimum information crit model order: %0.2f\n',mean(tbl_summary_out.min_modorder_info_crit_aic_line));
fprintf('HQ mean histogram model order: %0.2f\n',mean(tbl_summary_out.mean_hist_hq_amnts));
fprintf('AIC mean histogram model order: %0.2f\n',mean(tbl_summary_out.mean_hist_aic_amnts));

%% ===================================================================== %%
%## LOAD BOOSTRAP & PHASE RANDOMIZED INFORMATION
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\_save\01042024_indv_morder';
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\conn\01082024_dDTF08';
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\conn\S';
save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\conn\dDTF08';

save_dir_nonzmats =  [save_dir filesep 'nonzero_mat_files'];
% save_dir_bootmats =  [save_dir filesep 'bootstrap_mat_files'];
save_dir_bootmats =  [save_dir filesep 'bootstrap_baseline_mat_files'];
BOOTSTRAP_STRUCT = cell(length(MAIN_ALLEEG),1);
PHASERND_STRUCT = cell(length(MAIN_ALLEEG),1);
for subj_i = 1:length(MAIN_ALLEEG)
    EEG = MAIN_ALLEEG(subj_i);
    BOOTSTRAP_STRUCT{subj_i} = par_load(save_dir_bootmats,sprintf('%s_boot_struct.mat',EEG.subject));
    PHASERND_STRUCT{subj_i} = par_load(save_dir_nonzmats,sprintf('%s_phasernd_struct.mat',EEG.subject));
end
%-
BOOTSTRAP_STRUCT = cat(1,BOOTSTRAP_STRUCT{:});
PHASERND_STRUCT = cat(1,PHASERND_STRUCT{:});
CONN_MEAS = BOOTSTRAP_STRUCT(1).conn_meas;
FREQ_BANDS = BOOTSTRAP_STRUCT(1).FREQ_BAND_INF.FREQ_BANDS;
%% ===================================================================== %%
%## BOOTSTRAP EXTRACTION
BOOT_CONN_AVG_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
BOOT_CONN_FREQ_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),T_DIM,length(fieldnames(FREQ_BANDS)),length(BOOTSTRAP_STRUCT));
%% ===================================================================== %%
%## NONZERO EXTRACTION
NON_CONN_AVG_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(PHASERND_STRUCT));
NON_CONN_FREQ_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),T_DIM,length(fieldnames(FREQ_BANDS)),length(PHASERND_STRUCT));
%% CREATE CONNECTION SUMMARY STATISTICS
% CONNECT_TABLE = table;
% cluster_ints = [3,4,5,6,7,8,9,10,11,12];
% cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
cluster_ints = [3,7,5,4];
cluster_names = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
ALPHA = 0.1;
COND_PAIR_ORDER = [1,2];
CONN_MEAS = PHASERND_STRUCT(1).conn_meas;
FREQ_NAMES = fieldnames(FREQ_BANDS);
allfreqs = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.freqs;
alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
rej_subj = zeros(size(comps_out,2),1); 
%## REJECT SUBJECTS
%- method 1 (not so good)
% bools = zeros(length(cluster_ints),length(cluster_ints),size(NON_CONN_MASK_STORE,5));
% rej_subj = zeros(size(NON_CONN_MASK_STORE,5),1);
% for subj_i = 1:size(NON_CONN_MASK_STORE,5)
%     tmp = squeeze(NON_CONN_MASK_STORE(cluster_ints,cluster_ints,:,:,subj_i,:));
%     tmp(tmp>0) = 1;
% %             tmp = ~tmp;
%     bools(:,:,subj_i) = any(tmp,[3,4,5]);
% %     disp(subj_i);
%     disp(bools(:,:,subj_i));
%     if ~all(bools(:,:,subj_i),[1,2])
%         fprintf('Reject: %i\n',subj_i);
%         rej_subj(subj_i) = 1;
%     end
% end
%- method 2
% rej_subj = zeros(size(comps_out,2),1); 
% for subj_i = 1:size(comps_out,2)
%     tmp = comps_out(cluster_ints,subj_i)>0;
%     if ~all(tmp)
%         fprintf('Reject: %i\n',subj_i);
%         rej_subj(subj_i) = 1;
%     end
% end

%## NONZERO STATS
subj_chars = {PHASERND_STRUCT(~rej_subj).subject};
NON_CONN_MASK_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(PHASERND_STRUCT),length(COND_NAMES));
for subj_i = 1:length(PHASERND_STRUCT)
    for cond_i = 1:length(COND_NAMES)
        comps = squeeze(comps_out(:,subj_i));
        [tmpcl,idxcl] = sort(comps);
        idxcl = idxcl(tmpcl~=0);
        nonz_in = (PHASERND_STRUCT(subj_i).stats{cond_i}.(CONN_MEAS).pval<ALPHA);
        boot_in = BOOTSTRAP_STRUCT(subj_i).averages{cond_i,1};
        %- assign
        for i = 1:length(idxcl)
            for j = 1:length(idxcl)
                %- method 1
                %- method 2
        %         tmp = PHASERND_STRUCT(subj_i).masked_conn{cond_i};
                NON_CONN_MASK_STORE(idxcl(i),idxcl(j),:,:,subj_i,cond_i) = squeeze(nonz_in(i,j,:,:)).*squeeze(boot_in(i,j,:,:));
            end
        end
    end
end
% NON_CONN_MASK_STORE = NON_CONN_MASK_STORE(cluster_ints,cluster_ints,:,:,:,:);
NON_CONN_MASK_STORE = NON_CONN_MASK_STORE(cluster_ints,cluster_ints,:,:,~rej_subj,:);
NON_CONN_MASK_STORE(isnan(NON_CONN_MASK_STORE)) = 0;

%## BOOTSTRAP STATS (VERY CONSERVATIVE)
BOOT_CONN_MASK_STORE = zeros(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(BOOTSTRAP_STRUCT),length(COND_NAMES));
for subj_i = 1:length(BOOTSTRAP_STRUCT)
    for cond_i = 1:length(COND_NAMES)
        comps = squeeze(comps_out(:,subj_i));
        [tmpcl,idxcl] = sort(comps);
        idxcl = idxcl(tmpcl~=0);
        %- for Hab (between condition) bootstraps
%         tmp = (PHASERND_STRUCT(subj_i).stats{cond_i}.(CONN_MEAS).pval<ALPHA).*(BOOTSTRAP_STRUCT(subj_i).stats{COND_PAIR_ORDER(1),COND_PAIR_ORDER(2)}.(CONN_MEAS).pval<ALPHA).*BOOTSTRAP_STRUCT(subj_i).averages{cond_i,1};
%         tmp = BOOTSTRAP_STRUCT(subj_i).masked_conn{COND_PAIR_ORDER(1),COND_PAIR_ORDER(2)};
        %- for baseline bootstraps
        tmp = BOOTSTRAP_STRUCT(subj_i).masked_conn{cond_i};
        for i = 1:length(idxcl)
            for j = 1:length(idxcl)
                BOOT_CONN_MASK_STORE(idxcl(i),idxcl(j),:,:,subj_i,cond_i) = tmp(i,j,:,:);
            end
        end
    end
end
% BOOT_CONN_MASK_STORE = BOOT_CONN_MASK_STORE(cluster_ints,cluster_ints,:,:,:,:);
BOOT_CONN_MASK_STORE = BOOT_CONN_MASK_STORE(cluster_ints,cluster_ints,:,:,~rej_subj,:);
BOOT_CONN_MASK_STORE(isnan(BOOT_CONN_MASK_STORE)) = 0;
%##

%## average across subjects before averaging across time and summing
%frequency***
NON_CONN_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),F_DIM,T_DIM,length(COND_NAMES));
BOOT_CONN_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),F_DIM,T_DIM,length(COND_NAMES));
NON_FREQBAND_TAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),size(NON_CONN_MASK_STORE,5),length(COND_NAMES));
BOOT_FREQBAND_TAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),size(NON_CONN_MASK_STORE,5),length(COND_NAMES));
for cond_i = 1:size(NON_CONN_MASK_STORE,6)
    %- median
    NON_CONN_SUBJAVG(:,:,:,:,cond_i) = squeeze(median(squeeze(NON_CONN_MASK_STORE(:,:,:,:,:,cond_i)),5));
    BOOT_CONN_SUBJAVG(:,:,:,:,cond_i) = squeeze(median(squeeze(BOOT_CONN_MASK_STORE(:,:,:,:,:,cond_i)),5));
    %- mean
%     NON_CONN_MASK_STORE(:,:,:,:,:,cond_i) = mean(NON_CONN_MASK_STORE(:,:,:,:,:,cond_i),5);
%     BOOT_CONN_MASK_STORE(:,:,:,:,:,cond_i) = mean(BOOT_CONN_MASK_STORE(:,:,:,:,:,cond_i),5);
    for freq_i = 1:length(FREQ_NAMES)
        %- time-before subj method 1
        NON_FREQBAND_TAVG(:,:,freq_i,:,cond_i) = squeeze(mean(squeeze(trapz(NON_CONN_MASK_STORE(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,:,cond_i),3)),3));
        BOOT_FREQBAND_TAVG(:,:,freq_i,:,cond_i) = squeeze(mean(squeeze(trapz(BOOT_CONN_MASK_STORE(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,:,cond_i),3)),3));
        
    end
end
%## FREQ-TIME AVERAGES
NON_FREQBAND_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),T_DIM,length(COND_NAMES));
BOOT_FREQBAND_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),T_DIM,length(COND_NAMES));
NON_FREQBAND_TAVG_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),length(COND_NAMES));
BOOT_FREQBAND_TAVG_SUBJAVG = zeros(length(cluster_ints),length(cluster_ints),length(FREQ_NAMES),length(COND_NAMES));

for freq_i = 1:length(FREQ_NAMES)
    for cond_i = 1:size(NON_CONN_SUBJAVG,5)
        for i = 1:size(NON_CONN_SUBJAVG,1)
            for j = 1:size(NON_CONN_SUBJAVG,2)
                %- freq method 1
                NON_FREQBAND_SUBJAVG(i,j,freq_i,:,cond_i) = squeeze(trapz(NON_CONN_SUBJAVG(i,j,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3));
                BOOT_FREQBAND_SUBJAVG(i,j,freq_i,:,cond_i) = squeeze(trapz(BOOT_CONN_SUBJAVG(i,j,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3));
                %- freq method 2
        %         NON_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = sum(NON_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3);
        %         BOOT_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = sum(BOOT_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3);
                %- time method 1
                NON_FREQBAND_TAVG_SUBJAVG(i,j,freq_i,cond_i) = squeeze(mean(NON_FREQBAND_SUBJAVG(i,j,freq_i,:,cond_i),4));
                BOOT_FREQBAND_TAVG_SUBJAVG(i,j,freq_i,cond_i) = squeeze(mean(BOOT_FREQBAND_SUBJAVG(i,j,freq_i,:,cond_i),4));
                %- time method 2
        %         NON_FREQBAND_TAVG_SUBJAVG(:,:,freq_i,cond_i) = squeeze(median(NON_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i),4));
        %         BOOT_FREQBAND_TAVG_SUBJAVG(:,:,freq_i,cond_i) = squeeze(median(BOOT_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i),4));
                %- freq-time method 3
%                 tmp1 = squeeze(NON_CONN_SUBJAVG(i,j,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i));
%                 tmp2 = squeeze(BOOT_CONN_SUBJAVG(i,j,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i));
%                 NON_FREQBAND_TAVG_SUBJAVG(i,j,freq_i,cond_i) = mean(tmp1(:));
%                 BOOT_FREQBAND_TAVG_SUBJAVG(i,j,freq_i,cond_i) = mean(tmp2(:));
            end
        end
        %- freq method 1
%         NON_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = squeeze(trapz(NON_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3));
%         BOOT_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = squeeze(trapz(BOOT_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3));
%         %- freq method 2
% %         NON_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = sum(NON_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3);
% %         BOOT_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i) = sum(BOOT_CONN_SUBJAVG(:,:,FREQ_BANDS.(FREQ_NAMES{freq_i}),:,cond_i),3);
%         %- time method 1
%         NON_FREQBAND_TAVG_SUBJAVG(:,:,freq_i,cond_i) = squeeze(mean(NON_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i),4));
%         BOOT_FREQBAND_TAVG_SUBJAVG(:,:,freq_i,cond_i) = squeeze(mean(BOOT_FREQBAND_SUBJAVG(:,:,freq_i,:,cond_i),4));
    end
end
cl_pairs = [];
for cl_i = 1:length(cluster_ints)
    for cl_j = 1:length(cluster_ints)
        if any(all(ismember(cl_pairs,[cl_i,cl_j]),2)) || cl_i==cl_j
            continue;
        end
        cl_pairs = [cl_pairs; cl_i,cl_j];
    end
end
VALS_OUT = struct('cluster_inds',cluster_ints,...
    'cluster_chars',{cluster_names},...
    'subjs_rejected',rej_subj,...
    'freq_bands',FREQ_BANDS,...
    'freq_names',{fieldnames(FREQ_BANDS)},...
    'conn_freqs',allfreqs,...
    'conn_times',alltimes,...
    'subjects',{subj_chars},...
    'cond_names',{COND_NAMES},...
    'cluster_pairs',cl_pairs,...
    'nonz_conn_cli_clj_f_t_s_c',NON_CONN_MASK_STORE,...
    'boot_conn_cli_clj_f_t_s_c',BOOT_CONN_MASK_STORE,...
    'nonz_conn_cli_clj_fi_s_c',NON_FREQBAND_TAVG,...
    'boot_conn_cli_clj_fi_s_c',BOOT_FREQBAND_TAVG,...
    'nonz_conn_cli_clj_f_t_c',NON_CONN_SUBJAVG,...
    'boot_conn_cli_clj_f_t_c',BOOT_CONN_SUBJAVG,...
    'nonz_conn_cli_clj_fi_t_c',NON_FREQBAND_SUBJAVG,...
    'boot_conn_cli_clj_fi_t_c',BOOT_FREQBAND_SUBJAVG,...
    'nonz_conn_cli_clj_fi_c',NON_FREQBAND_TAVG_SUBJAVG,...
    'boot_conn_cli_clj_fi_c',BOOT_FREQBAND_TAVG_SUBJAVG);
% clear BOOT_FREQBAND_SUBJAVG NON_FREQBAND_TAVG_SUBJAVG BOOT_FREQBAND_TAVG_SUBJAVG NON_CONN_MASK_STORE BOOT_CONN_MASK_STORE NON_CONN_SUBJAVG BOOT_CONN_SUBJAVG NON_FREQBAND_SUBJAVG

%##
%{
par_save(VALS_OUT.nonz_conn_cli_clj_fi_s_c,destination_folder,'nonz_cli_clj_fi_s_c.mat');
par_save(VALS_OUT.boot_conn_cli_clj_fi_s_c,destination_folder,'boot_cli_clj_fi_s_c.mat');
par_save(VALS_OUT.nonz_conn_cli_clj_fi_c,destination_folder,'nonz_cli_clj_fi_c.mat');
par_save(VALS_OUT.boot_conn_cli_clj_fi_c,destination_folder,'boot_cli_clj_fi_c.mat');
%}
%% ===================================================================== %%
%## NONZERO CONNECTIVITY MATRICIES
COLOR_LIMITS=[0,20e-4];
net_vals_freq_chars = {'theta','alpha','beta','all'};
for cond_i = 1:length(COND_NAMES)
    for freq_i = 1:length(VALS_OUT.freq_names)
        FREQS_SUBPATH = sprintf('%s',VALS_OUT.freq_names{freq_i});
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cl_chars = cell(length(VALS_OUT.cluster_chars),1);
        for i = 1:length(VALS_OUT.cluster_chars)
            cl_chars{i} = sprintf('(N=%i) %s',length(VALS_OUT.subjects),VALS_OUT.cluster_chars{i});
        end
        %## PLOT
%         tmp = squeeze(VALS_OUT.nonz_conn_cli_clj_fi_c(:,:,freq_i,cond_i));
        tmp = zeros(length(cluster_ints),length(cluster_ints));
        for i = 1:length(cluster_ints)
            for j = 1:length(cluster_ints)
                tt = squeeze(NON_CONN_MASK_STORE(i,j,:,:,:,cond_i));
                tt = mean(tt,3);
                tt = tt(VALS_OUT.freq_bands.(VALS_OUT.freq_names{freq_i}),:);
                tt = mean(tt(:));
                disp(tt)
                tmp(i,j) = tt;
            end
        end
        disp(tmp)
%         tmp = net_vals_ave.(net_vals_freq_chars{freq_i})*10^4;
        %- nan out diagnol & zeros
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
        %- plot
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        %- create title
        title(sprintf('(%s) ''%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,VALS_OUT.cond_names{cond_i}));
        hnd.YDisplayLabels = cl_chars;
        hnd.XDisplayLabels = cl_chars;
        hnd.ColorLimits = COLOR_LIMITS;
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        hnd.NodeChildren(3).Title.Interpreter = 'none';
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[10,100,720,620])
        exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s_cond_nonz_mat.jpg',strjoin(strsplit(VALS_OUT.cond_names{cond_i},' '),'_'))],'Resolution',300);
    end
end
%% ===================================================================== %%
%## DIFFERENCE NONZERO CONNECTIVITY MATRICIES
cond_1 = 1;
cond_2 = 2;
COLOR_LIMITS=[0,20e-4];
for freq_i = 1:length(VALS_OUT.freq_names)
    FREQS_SUBPATH = sprintf('%s',VALS_OUT.freq_names{freq_i});
    save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
    if ~exist(save_sub_figs,'dir')
        mkdir(save_sub_figs);
    end
    %- plot
    cl_chars = cell(length(VALS_OUT.cluster_chars),1);
    for i = 1:length(VALS_OUT.cluster_chars)
        cl_chars{i} = sprintf('(N=%i) %s',length(VALS_OUT.subjects),VALS_OUT.cluster_chars{i});
    end
    %## PLOT
    tmp = squeeze(VALS_OUT.nonz_conn_cli_clj_fi_c(:,:,freq_i,cond_1))-squeeze(VALS_OUT.nonz_conn_cli_clj_fi_c(:,:,freq_i,cond_2));
    %- nan out diagnol & zeros
    I = eye(size(tmp));
    I = (I == 0);
    tmp = tmp.*I;
    tmp(tmp == 0) = nan();
    %- plot
    figure;
    hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
    %- create title
    title(sprintf('(%s) ''%s-%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,VALS_OUT.cond_names{cond_1},VALS_OUT.cond_names{cond_2}));
    hnd.YDisplayLabels = cl_chars;
    hnd.XDisplayLabels = cl_chars;
    hnd.ColorLimits = COLOR_LIMITS;
    hnd.GridVisible = 'off';
    hnd.CellLabelFormat = '%0.1g';
    hnd.NodeChildren(3).Title.Interpreter = 'none';
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Position',[10,100,720,620])
    exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s-%s_nonz_diff_mat.jpg',strjoin(strsplit(VALS_OUT.cond_names{cond_1},' '),'_'),strjoin(strsplit(VALS_OUT.cond_names{cond_2},' '),'_'))],'Resolution',300);
end
%% ===================================================================== %%
%## BOOTSTAT CONSERVATIVE CONNECTIVITY MATRICIES
COLOR_LIMITS=[0,20e-4];
for cond_i = 1:length(COND_NAMES)
    for freq_i = 1:length(VALS_OUT.freq_names)
        FREQS_SUBPATH = sprintf('%s',VALS_OUT.freq_names{freq_i});
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cl_chars = cell(length(VALS_OUT.cluster_chars),1);
        for i = 1:length(VALS_OUT.cluster_chars)
            cl_chars{i} = sprintf('(N=%i) %s',length(VALS_OUT.subjects),VALS_OUT.cluster_chars{i});
        end
        %## PLOT
        tmp = squeeze(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_i));
        %- nan out diagnol & zeros
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
        %- plot
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        %- create title
        title(sprintf('(%s) ''%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,VALS_OUT.cond_names{cond_i}));
        hnd.YDisplayLabels = cl_chars;
        hnd.XDisplayLabels = cl_chars;
        hnd.ColorLimits = COLOR_LIMITS;
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        hnd.NodeChildren(3).Title.Interpreter = 'none';
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[10,100,720,620])
        exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s_bootconserv_bootavg_mat.jpg',strjoin(strsplit(VALS_OUT.cond_names{cond_i},' '),'_'))],'Resolution',300);
    end
end
%% ===================================================================== %%
%## DIFFERENCE BOOTSTAT CONSERVATIVE CONNECTIVITY MATRICIES
COLOR_LIMITS=[-20e-4,20e-4];
cond_1 = 1;
cond_2 = 2;
for freq_i = 1:length(VALS_OUT.freq_names)
    FREQS_SUBPATH = sprintf('%s',VALS_OUT.freq_names{freq_i});
    save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
    if ~exist(save_sub_figs,'dir')
        mkdir(save_sub_figs);
    end
    %- plot
    cl_chars = cell(length(VALS_OUT.cluster_chars),1);
    for i = 1:length(VALS_OUT.cluster_chars)
        cl_chars{i} = sprintf('(N=%i) %s',length(VALS_OUT.subjects),VALS_OUT.cluster_chars{i});
    end
    %## PLOT
    tmp = squeeze(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_1))-squeeze(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_2));
    %- nan out diagnol & zeros
    I = eye(size(tmp));
    I = (I == 0);
    tmp = tmp.*I;
    tmp(tmp == 0) = nan();
    %- plot
    figure;
    hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
    %- create title
    title(sprintf('(%s) ''%s-%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,VALS_OUT.cond_names{cond_1},VALS_OUT.cond_names{cond_2}));
    hnd.YDisplayLabels = cl_chars;
    hnd.XDisplayLabels = cl_chars;
    hnd.ColorLimits = COLOR_LIMITS;
    hnd.GridVisible = 'off';
    hnd.CellLabelFormat = '%0.1g';
    hnd.NodeChildren(3).Title.Interpreter = 'none';
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Position',[10,100,720,620])
    exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s-%s_diff_bootconserv_mat.jpg',strjoin(strsplit(VALS_OUT.cond_names{cond_1},' '),'_'),strjoin(strsplit(VALS_OUT.cond_names{cond_2},' '),'_'))],'Resolution',300);
end
%% ===================================================================== %%
%## (LOWER TRIANGLE) BOOTSTAT CONSERVATIVE CONNECTIVITY MATRICIES
COLOR_LIMITS=[-20e-4,20e-4];
cond_1 = 1;
cond_2 = 2;
for freq_i = 1:length(VALS_OUT.freq_names)
    FREQS_SUBPATH = sprintf('%s',VALS_OUT.freq_names{freq_i});
    save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
    if ~exist(save_sub_figs,'dir')
        mkdir(save_sub_figs);
    end
    %- plot
    cl_chars = cell(length(VALS_OUT.cluster_chars),1);
    for i = 1:length(VALS_OUT.cluster_chars)
        cl_chars{i} = sprintf('(N=%i) %s',length(VALS_OUT.subjects),VALS_OUT.cluster_chars{i});
    end
    %## PLOT
    tmp1 = zeros(size(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_1)));
    tmp2 = zeros(size(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_1)));
    for i = 1:size(VALS_OUT.cluster_pairs)
        pair = VALS_OUT.cluster_pairs(i,:);
        tmp11 = VALS_OUT.boot_conn_cli_clj_fi_c(pair(1),pair(2),freq_i,cond_1);
        tmp12 = VALS_OUT.boot_conn_cli_clj_fi_c(pair(2),pair(1),freq_i,cond_1);
        tmp1(pair(2),pair(1)) = tmp11-tmp12;
        tmp21 = VALS_OUT.boot_conn_cli_clj_fi_c(pair(1),pair(2),freq_i,cond_2);
        tmp22 = VALS_OUT.boot_conn_cli_clj_fi_c(pair(2),pair(1),freq_i,cond_2);
        tmp2(pair(2),pair(1)) = tmp21-tmp22;
    end
%     tmp = squeeze(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_1))-squeeze(VALS_OUT.boot_conn_cli_clj_fi_c(:,:,freq_i,cond_2));
    tmp = tmp1 - tmp2;
    %- nan out diagnol & zeros
    I = eye(size(tmp));
    I = (I == 0);
    tmp = tmp.*I;
    tmp(tmp == 0) = nan();
    %- plot
    figure;
    hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
    %- create title
    title(sprintf('(%s) ''%s-%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,VALS_OUT.cond_names{cond_1},VALS_OUT.cond_names{cond_2}));
    hnd.YDisplayLabels = cl_chars;
    hnd.XDisplayLabels = cl_chars;
    hnd.ColorLimits = COLOR_LIMITS;
    hnd.GridVisible = 'off';
    hnd.CellLabelFormat = '%0.1g';
    hnd.NodeChildren(3).Title.Interpreter = 'none';
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Position',[10,100,720,620])
    exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s-%s_lowertri_bootavg_mat.jpg',strjoin(strsplit(VALS_OUT.cond_names{cond_1},' '),'_'),strjoin(strsplit(VALS_OUT.cond_names{cond_2},' '),'_'))],'Resolution',300);
end
%% ===================================================================== %%
%## TXF COMPARISONS
if ~ispc
    addpath(convertPath2UNIX('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
else
    addpath(convertPath2Drive('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
end

%-
save_sub_figs = [save_dir filesep 'cluster_level_permutation'];
if ~exist(save_sub_figs,'dir')
    mkdir(save_sub_figs);
end
%##
%-
cond_1 = 1;
cond_2 = 2;
REPEATED_MEAS = 0.05;
ALPHA = 0.05;
NUM_ITERS = 200;
%-
allfreqs = VALS_OUT.conn_freqs; %MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.freqs;
alltimes = VALS_OUT.conn_times; %MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
PLOT_STRUCT = struct('figure_position',[100,100,350,350],...
    'xtick_label','Time (s)',...
    'ytick_label','Frequency (Hz)',...
    'clim',[-5e-3,5e-3],...
    'font_size',12,...
    'freq_lims',[4,60],...
    'time_lims',[alltimes(1),alltimes(end)],...
    'subplot_width',0.20,...
    'subplot_height',0.7,...
    'subplot_shift',0.24,...
    'colorbar_shift',0.1,...
    'plot_inchs',[3,3,10,5],...
    'colormap',linspecer,...
    'event_times',[0],...
    'event_chars',{{'subject hit'}},...
    'title','',...
    'subplot_titles',{{'HUMAN','BM','STATS'}},...
    'cbar_intv',1e-3,...
    'cbar_label','Connectivity');
for i = 1:length(VALS_OUT.cluster_inds)
    for j = 1:length(VALS_OUT.cluster_inds)
        %-
        sig1_in = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,:,:,:,cond_1));
        sig1_in(isnan(sig1_in)) = 0;
        %-
        sig2_in = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,:,:,:,cond_2));
        sig2_in(isnan(sig2_in)) = 0;
        combase = cat(3,sig1_in,sig2_in);
        combase = median(combase,3);
%         combase = mean(combase,3);
        tmp1 = [];
        tmp2 = [];
        for subj_i = 1:size(sig1_in,3)
            tmp1 = cat(3,tmp1,sig1_in(:,:,subj_i)-combase);
            tmp2 = cat(3,tmp2,sig2_in(:,:,subj_i)-combase);
        end
        %-
        [mask,tscore,pval] = clusterLevelPermutationTest(sig1_in,sig2_in,...
            REPEATED_MEAS,ALPHA,NUM_ITERS);
%         [mask,tscore,pval] = clusterLevelPermutationTest(tmp1,tmp2,...
%             REPEATED_MEAS,ALPHA,NUM_ITERS);
%         allersp = {median(tmp1,3),median(tmp2,3),mask};
        allersp = {median(sig1_in,3),median(sig2_in,3),mask};
        fig = plot_multipane_txf(allersp,alltimes,allfreqs,...
            'PLOT_STRUCT',PLOT_STRUCT);
        exportgraphics(fig,[save_sub_figs filesep sprintf('%sto%s_txf_comboplot.jpg',VALS_OUT.cluster_chars{i},VALS_OUT.cluster_chars{j})]);
    end
end

%% ===================================================================== %%
%## SIGNAL COMPARISONS
% freq_1 = 1;
cond_1 = 1;
cond_2 = 2;
% COLOR_LIMITS = [-20e-4,20e-4];
% COLOR_LIMITS = [0,1];
COLOR_LIMITS = [0,5e-3];
REPEATED_MEAS = 0.05;
ALPHA = 0.05;
NUM_ITERS = 200;
alltimes = VALS_OUT.conn_times; %MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
%##
% for freq_i = 1:length(VALS_OUT.freq_names)
    for i = 1:length(VALS_OUT.cluster_inds)
        for j = 1:length(VALS_OUT.cluster_inds)
%             dat_in1 = VALS_OUT.nonz_conn_cli_clj_fi_t_c(i,j,freq_i,:,cond_1);
%             dat_in2 = VALS_OUT.nonz_conn_cli_clj_fi_t_c(i,j,freq_i,:,cond_2);
%             dat_in1 = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,VALS_OUT.freq_bands.(VALS_OUT.freq_names{freq_i}),:,:,cond_1));
%             dat_in2 = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,VALS_OUT.freq_bands.(VALS_OUT.freq_names{freq_i}),:,:,cond_2));
            dat_in1 = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,:,:,:,cond_1));
            dat_in2 = squeeze(VALS_OUT.nonz_conn_cli_clj_f_t_s_c(i,j,:,:,:,cond_2));
            %- method 1
%             dat_in1 = squeeze(sum(dat_in1,1));
%             dat_in2 = squeeze(sum(dat_in1,1));
            %- method 2
            dat_in1 = squeeze(trapz(dat_in1,1));
            dat_in2 = squeeze(trapz(dat_in2,1));
%             y1_mean = mean(dat_in1,2); %rand(1,10); % your mean vector;
            y2_mean = median(dat_in2,2);
            y1_mean = median(dat_in1,2);
            %- stats method
%             dat_in1 = (trapz(dat_in1,1));
%             dat_in2 = (trapz(dat_in2,1));
%             dat_in1 = reshape(dat_in1,size(dat_in1,2),size(dat_in1,1),size(dat_in1,3));
%             dat_in2 = reshape(dat_in2,size(dat_in2,2),size(dat_in2,1),size(dat_in2,3));
%             y1_mean = squeeze(mean(dat_in1,3)); %rand(1,10); % your mean vector;
%             y2_mean = squeeze(mean(dat_in2,3)); %rand(1,10); % your mean vector;
%             [mask,tscore,pval] = clusterLevelPermutationTest(dat_in1,dat_in2,...
%                 REPEATED_MEAS,ALPHA,NUM_ITERS);
            
            fig = figure();
            hold on;
%             curve1 = quantile(tmp,0.95,2);
%             curve2 = quantile(tmp,0.25,2);
            %-
%             plot(alltimes,curve1,'Color',[COLORS_MEAN(freq_i,:),0.5],'DisplayName','95 percentile');
%             plot(alltimes,y_med,'Color',[COLORS_MEAN(freq_i,:),0.5],'DisplayName','median');
            plot(alltimes,y1_mean,'LineWidth',2,'Color','r',...
                'DisplayName',sprintf('%s',VALS_OUT.cond_names{cond_1}));
            plot(alltimes,y2_mean,'LineWidth',2,'Color','g',...
                'DisplayName',sprintf('%s',VALS_OUT.cond_names{cond_2}));
            if ~isempty(mask)
                plot(alltimes,mask,'LineWidth',2,'Color','k',...
                    'DisplayName',sprintf('%s',VALS_OUT.cond_names{cond_2}));
            end
            legend('Location','southeast');
            ylim(COLOR_LIMITS);
            title(sprintf('%s-%s',VALS_OUT.cluster_chars{i},VALS_OUT.cluster_chars{j}));
            hold off;

        end
    end
% end
%% ===================================================================== %%
%## time-frequency statistics?
cond_1 = 1;
cond_2 = 2;
ALPHA = 0.05;
% cl_1 = 3;
% cl_2 = 1;
save_sub_figs = [save_dir filesep 'txf_plots'];
if ~exist(save_sub_figs,'dir')
    mkdir(save_sub_figs);
end
for cl_i = 1:size(cl_pairs,1)
% for cl_i = 1
    cl_1 = cl_pairs(cl_i,1)==CLUSTER_ITERS;
    cl_2 = cl_pairs(cl_i,2)==CLUSTER_ITERS;
    nonz_conn_1 = NON_CONN_AVG_STORE;
    nonz_conn_2 = NON_CONN_AVG_STORE;
    nonz_stat_1 = NON_CONN_AVG_STORE;
    nonz_stat_2 = NON_CONN_AVG_STORE;
    for subj_i = 1:length(PHASERND_STRUCT)
        comps = squeeze(comps_out(:,subj_i));
        [tmpcl,idxcl] = sort(comps);
        idxcl = idxcl(tmpcl~=0);
        tmpcl = tmpcl(tmpcl~=0);
        nonz_conn_1(idxcl,idxcl,:,:,subj_i) = PHASERND_STRUCT(subj_i).masked_conn{cond_1}; %.(CONN_MEAS); %human
        nonz_conn_2(idxcl,idxcl,:,:,subj_i) = PHASERND_STRUCT(subj_i).masked_conn{cond_2}; %.(CONN_MEAS %ball machine
        nonz_stat_1(idxcl,idxcl,:,:,subj_i) = PHASERND_STRUCT(subj_i).stats{cond_1}.(CONN_MEAS).pval; %.(CONN_MEAS); %human
        nonz_stat_2(idxcl,idxcl,:,:,subj_i) = PHASERND_STRUCT(subj_i).stats{cond_2}.(CONN_MEAS).pval; %.(CONN_MEAS %ball machine
    end
    tmp1_o = [];
    tmp2_o = [];
    for subj_i = 1:length(PHASERND_STRUCT)
        tmp1 = squeeze(nonz_conn_1(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,subj_i));
        tmp2 = squeeze(nonz_conn_2(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,subj_i));
        
        if ~all(isnan(tmp1),[1,2])
            tmp1_o = cat(3,tmp1_o,tmp1);
        end
        if ~all(isnan(tmp2),[1,2])
            tmp2_o = cat(3,tmp2_o,tmp2);
        end
    end
    tmp1 = tmp1_o;
    tmp2 = tmp2_o;
    disp_nonz_stat1 = {nanmean(squeeze(nonz_stat_1(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:))<ALPHA,3)};
    disp_nonz_stat2 = {nanmean(squeeze(nonz_stat_2(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:))<ALPHA,3)};

    % tmp11 = zeros(size(tmp1,2),size(tmp1,1),size(tmp1,3));
    % tmp22 = zeros(size(tmp1,2),size(tmp1,1),size(tmp1,3));
    % for subj_i = 1:size(tmp1,3)
    %     for t_i = 1:size(tmp1,2)
    %         tmp11(t_i,:,subj_i) = tmp1(:,t_i,subj_i);
    %         tmp22(t_i,:,subj_i) = tmp2(:,t_i,subj_i);
    %     end
    %     
    % end
    % allersp = {tmp11;tmp22};
    allersp = {tmp1;tmp2};
    %## time-frequency statistics?
    boot_conn_1 = BOOT_CONN_AVG_STORE;
    boot_conn_2 = BOOT_CONN_AVG_STORE;
    boot_stat_1 = BOOT_CONN_AVG_STORE;
    boot_stat_2 = BOOT_CONN_AVG_STORE;
    boot_ci_up = BOOT_CONN_AVG_STORE;
    boot_ci_low = BOOT_CONN_AVG_STORE;
    
    for subj_i = 1:length(BOOTSTRAP_STRUCT)
        comps = squeeze(comps_out(:,subj_i));
        [tmpcl,idxcl] = sort(comps);
        idxcl = idxcl(tmpcl~=0);
        tmpcl = tmpcl(tmpcl~=0);
    %     boot_conn_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).averages{cond_1,1}; %.(CONN_MEAS);
    %     boot_conn_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).averages{cond_2,1}; %.(CONN_MEAS
        boot_conn_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).masked_conn{cond_1,cond_2}; %.(CONN_MEAS);
%             boot_conn_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).masked_conn{cond_2,cond_1}; %.(CONN_MEAS
        boot_stat_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).stats{cond_1,cond_2}.(CONN_MEAS).pval; %.(CONN_MEAS);
%             boot_stat_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).stats{cond_2,cond_1}.(CONN_MEAS).pval; %.(CONN_MEAS
        boot_ci_up(idxcl,idxcl,:,:,subj_i) = squeeze(BOOTSTRAP_STRUCT(subj_i).stats{cond_1,cond_2}.(CONN_MEAS).ci(1,:,:,:,:));
        boot_ci_low(idxcl,idxcl,:,:,subj_i) = squeeze(BOOTSTRAP_STRUCT(subj_i).stats{cond_1,cond_2}.(CONN_MEAS).ci(2,:,:,:,:));
    end
    %##
    
    %##
    % allersp_boot = {squeeze(boot_conn_1(cl_1,cl_2,:,:,:));squeeze(boot_conn_2(cl_1,cl_2,:,:,:))};
    allersp_boot = {squeeze(boot_conn_1(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:))};
    tmp = squeeze(boot_stat_1(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:))<ALPHA;
    disp_boot_stat = {nanmean(tmp,3)};
    boot_ci_up = nanmean(squeeze(boot_ci_up(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:)),3);
    boot_ci_low = nanmean(squeeze(boot_ci_low(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,:,:)),3);
    %##
    %## ERSP PARAMS
    ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
        'groupstats','off',... %['on'|'off']
        'method','perm',... % ['param'|'perm'|'bootstrap']
        'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
        'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
        'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
        'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
        'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
        'fieldtripnaccu',2000);
    ERSP_PARAMS = struct('subbaseline','off',...
        'timerange',[],...
        'ersplim',[-2,2],...
        'freqfac',4,...
        'cycles',[3,0.8],...
        'freqrange',[1,100]);
    %## POP PARAMS
    MAIN_STUDY = pop_erspparams(MAIN_STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
          'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
    MAIN_STUDY = pop_statparams(MAIN_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    %- get stats parameters
    MAIN_STUDY = std_makedesign(MAIN_STUDY,MAIN_ALLEEG,1,'variable1','bounces','values1',{'2Bounce_Human','2Bounce_BM'});
    stats = MAIN_STUDY.etc.statistics;
    stats.paired = {'off'}; 
    [pcond, pgroup, pinter] = std_stat(allersp, stats);

    %##
    allfreqs = FREQS_VEC;
    alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes; %BOOTSTRAP_STRUCT(1).averages{1,1}.erWinCenterTimes;
    end_time = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.winCenterTimes(end); %BOOTSTRAP_STRUCT(1).averages{1,1}.erWinCenterTimes;
   
    allpcond = pcond{1};
    allersp_mean = cellfun(@(x) nanmean(x,3),allersp,'UniformOutput',false);
    allersp_boot = cellfun(@(x) nanmean(x,3),allersp_boot,'UniformOutput',false);
    allersp_diff = {allersp_mean{1}-allersp_mean{2}};
%     allersp_mean = [allersp_mean; disp_nonz_stat1; disp_nonz_stat2; allersp_diff; allersp_boot; disp_boot_stat ];
    allersp_mean = [allersp_mean; boot_ci_up; boot_ci_low; allersp_diff; allersp_boot; disp_boot_stat ];
    %%
    FIGURE_POSITION = [10,10,1020,480];
    clim_max = [-0.0005,0.0005]; %[-0.001,0.001];
    colormap_ersp = linspecer;
    SUB_FREQ_LIMS = [3,100];
    TIME_0_CHAR = 'subject hit';
%     alltitles = {'HUMAN','BALL MACHINE','NZ STAT HUMAN','NZ STAT BM','HUMAN - BM','BOOTSTRAP MASK','BOOTSTRAP STATS',};
    alltitles = {'HUMAN','BALL MACHINE','BOOT CI UPPER','BOOT CI LOWER','HUMAN - BM','BOOTSTRAP MASK','BOOTSTRAP STATS',};
    
    XTICK_LABEL = 'time (s)';
    YTICK_LABEL = 'Frequency (Hz)';
    FONT_SIZE = 12;
    SUBPLOT_WIDTH = 0.08;
    SUBPLOT_HEIGHT = 0.7;
    SHIFT_AMNT = 0.10;
    STATS_TITLE = 'CLUSTER STATS';
    if length(clim_max) == 2
        clim_ersp = clim_max;
    else
        clim_ersp = [-clim_max,clim_max];d
    end
    %##
    fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    set(fig,'Units','inches','Position',[3 3 18 5])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    horiz_shift = 0;
    hold on;
    for j = 1:length(allersp_mean)
        subplot(1,length(allersp_mean)+1,j); %,'position',[0.01+horiz_shift,0.1,0.5,0.5])
        ax = gca;
        tftopo(allersp_mean{j},alltimes,allfreqs,'limits',... 
            [nan nan nan nan clim_ersp],...
            'logfreq','native');
        hold on;
        colormap(colormap_ersp);
        %- adjust subplot position and height
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
    %         disp(get(ax,'Position'));
        %- add vertical line
        xline(ax,0,'k--');
        %- set ylims
        ylim(log(SUB_FREQ_LIMS))
        if SUB_FREQ_LIMS(2) <= 50
            set(ax,'YTick',log([4.01,8,13,30,50])); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end  
        %- set color lims
        set(ax,'clim',clim_ersp);
        %- set x-axis & y-axis labels
        if j == 1
            ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
            xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        else
            xlabel('','FontSize',FONT_SIZE);
            ylabel('','fontsize',FONT_SIZE,'fontweight','bold');
        end
        %- set x-axis ticks
        set(ax,'XTick',0,'XTickLabel',TIME_0_CHAR);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- title
        title(alltitles{j});  
        horiz_shift = horiz_shift + SHIFT_AMNT;
    end
    %## Add Stats To Plot
    if ~isempty(allpcond)
        subplot(1,length(allersp_mean)+1,length(allersp_mean)+1) % add one subplot for stats
        tftopo(double(allpcond),alltimes,allfreqs,'limits',... 
            [nan nan nan nan  clim_ersp],...
            'logfreq','native')
        colormap(colormap_ersp);
        ax = gca;
        %-
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
        disp(get(ax,'Position'));
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.04;
        c.Limits = clim_ersp;
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',FONT_SIZE);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        set(hL,'Rotation',0);
        %- add vertical line
        xline(ax,0,'k--');
        %- set ylims
        ylim(log(SUB_FREQ_LIMS))
        if SUB_FREQ_LIMS(2) <= 50
            set(ax,'YTick',log([4.01,8,13,30,50])); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end  
        %- set color lims
        set(ax,'clim',clim_ersp);
        %- set y-axis labels
        xlabel('','FontSize',FONT_SIZE);
        ylabel(sprintf(''),'fontsize',FONT_SIZE,'fontweight','bold');
        %- set x-axis ticks
        set(ax,'XTick',0,'XTickLabel',TIME_0_CHAR);
        set(ax,'XTick',end_time(end))
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- title
        title(STATS_TITLE)
    else
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.05;
        c.Limits = clim_ersp;
    end
    hold off;
    fig = get(groot,'CurrentFigure');
    exportgraphics(fig,[save_sub_figs filesep sprintf('%sto%s_txf_comboplot.jpg',CLUSTER_ASSIGNMENTS{cl_1},CLUSTER_ASSIGNMENTS{cl_2})]);
end
