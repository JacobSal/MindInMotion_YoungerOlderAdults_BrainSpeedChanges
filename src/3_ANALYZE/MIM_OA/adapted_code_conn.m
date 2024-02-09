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
COND_CHARS =  {'2Bounce_Human','2Bounce_BM'};
% COND_CHARS =  {'1Bounce_Human','Serve_Human'};
% EVENT_CHARS = {'Subject_hit'}; 
EVENT_CHARS = {'Subject_receive'}; 
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
% dt = '06152023_bounces_1h2h2bm_JS';
% dt = '07272023_bounces_1h_2h_2bm_JS';
% dt = '08182023_bounces_1h_2h_2bm_JS';
% dt = '12182023_bounces_1h_2h_2bm_JS_0p25-1';
% dt = '12282023_bounces_1h_2bm_JS_n1-0p5';
dt = '01182023_subjrec_2bounces_1h_2bm_JS_n5-1p5';
% dt = '01252023_subjrec_2bounces_rally_serve_human_JS_n5-1p5';
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
    [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
    EVENT_COND_COMBOS = MAIN_STUDY.etc.a_epoch_process.epoch_chars;
    %## CUT OUT NON VALID CLUSTERS
    inds = setdiff(1:length(comps_out),valid_cls);
    comps_out(inds,:) = 0;
end
% mim_gen_cluster_figs(MAIN_STUDY,MAIN_ALLEEG,save_dir,...
%     'CLUSTERS_TO_PLOT',valid_cls);
%%
% CLUSTER_ITERS = [3,7,5,4];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
% CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
% CLUSTER_ITERS = [3,4,5,6,7,8,9,11,12];
% CLUSTER_ASSIGNMENTS = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','Cing','LSuppMotor'}; 
% (01/24/2024) JS, updated to only include clusters considered in AS
% neuroimage paper.
%%
save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\conn';
destination_folder = save_dir;  %change as appropriate
if ~exist(destination_folder,'dir')
    mkdir(destination_folder);
end
%##
% CONN_MEASURES = {'S'}; %{'dDTF08','S'};
CONN_MEASURES = {'dDTF08'}; %{'dDTF08','S'};
FREQ_CROP = [4:50];
nonz_folder = 'nonzero_mat_files';
% boot_folder = 'bootstrap_baseline_mat_files';
boot_folder = 'bootstrap_mat_files';
HAB_BOOT_INDS = [2,1]; %[1,2];
% cluster_ints = [3,4,5,6,7,8,9,10,11,12];
% cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
cluster_ints = [3,4,5,6,7,8,9,10,11,12];
sub_ints = [3,4,5,7];
cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','LTemp','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.

usefdr = 0;
ALPHA = 0.1;
BASELINE_TIME = [-0.5,0];
TIME_LIMS=[-0.5 0.5];
TIME_LIMS_AVE=[0 0.5];
% cluster_ints = [3,7,5,4];
% cluster_names = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
% COND_PAIR_ORDER = [1,2];
save_dir_bootmats =  [save_dir filesep CONN_MEASURES{1} filesep 'bootstrap_baseline_mat_files'];
BOOTSTRAP_STRUCT = par_load(save_dir_bootmats,sprintf('%s_boot_struct.mat',MAIN_ALLEEG(1).subject));
CONN_MEAS = BOOTSTRAP_STRUCT.conn_meas;
FREQ_BANDS = BOOTSTRAP_STRUCT.FREQ_BAND_INF.FREQ_BANDS;
FREQ_NAMES = fieldnames(FREQ_BANDS);
allfreqs = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.freqs;
% alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
cluster_struct = MAIN_STUDY.cluster;
subj_inds = 1:length(MAIN_ALLEEG);
fAllInds=find(allfreqs>=FREQ_CROP(1) & allfreqs<=FREQ_CROP(end));
tInds=find(MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes>=TIME_LIMS(1) & MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes<=TIME_LIMS(2));
alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
alltimes = alltimes(tInds);
%## LOAD TEMP EEG DATA
COND_NAMES = MAIN_STUDY.etc.a_epoch_process.epoch_chars; %unique(EEG.etc.conn_table.t_fNames);
% tmp = unique(MAIN_ALLEEG(1).etc.conn_table.t_fNames);
% COND_NAMES = cell(length(tmp),1);
% for i = 1:length(tmp)
%     out = strsplit(tmp{i},'.');
%     out = strsplit(out{1},'_');
%     out = strjoin(out(3:end),' ');
%     COND_NAMES{i} = out;
% end
%% MODEL ORDER & VALIDATION DATA
%{
[tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(save_dir,COND_CHARS);
%-
fprintf('HQ minimum information crit model order: %0.2f\n',median(tbl_summary_out.min_modorder_info_crit_hq_line));
fprintf('HQ minimum information crit model order: %0.2f\n',iqr(tbl_summary_out.min_modorder_info_crit_hq_line));
%-
fprintf('AIC minimum information crit model order: %0.2f\n',mean(tbl_summary_out.min_modorder_info_crit_aic_line));
fprintf('HQ mean histogram model order: %0.2f\n',mean(tbl_summary_out.mean_hist_hq_amnts));
fprintf('AIC mean histogram model order: %0.2f\n',mean(tbl_summary_out.mean_hist_aic_amnts));
%-
fprintf('Consistency: %0.2f\n',mean(tbl_summary_out.mean_perc_cons))
fprintf('ACF Whiteness: %0.2f\n',mean(tbl_summary_out.mean_white_sig_acf))
fprintf('LJB Whiteness: %0.2f\n',mean(tbl_summary_out.mean_white_sig_ljb))
fprintf('BOXP Whiteness: %0.2f\n',mean(tbl_summary_out.mean_white_sig_boxp))
fprintf('LIMCL Whiteness: %0.2f\n',mean(tbl_summary_out.mean_white_sig_limcl))
fprintf('Stability: %0.2f\n',mean(tbl_summary_out.mean_stability_ind))
writetable(tbl_summary_out,[save_dir filesep 'model_crit_summary.xlsx']);
%}
%%
%## REJECT SUBJECTS
rej_subj = zeros(size(comps_out,2),1); 
for subj_i = 1:size(comps_out,2)
    tmp = comps_out(sub_ints,subj_i)>0;
    if ~all(tmp)
        fprintf('Reject: %i\n',subj_i);
        rej_subj(subj_i) = 1;
    end
end
subj_inds = find(~rej_subj);
% subj_chars = {PHASERND_STRUCT.subject};
%## 
for conn_i = 1:length(CONN_MEASURES)
    for cond_i=1:length(COND_NAMES)
        fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data' filesep COND_NAMES{cond_i}];
        if ~exist(fpath,'dir')
            mkdir(fpath);
        end
        connStruct_boot=cell(length(cluster_struct),length(cluster_struct));
        connStruct_nonz_pconn = cell(length(cluster_struct),length(cluster_struct));
        connStruct_boot_pconn = cell(length(cluster_struct),length(cluster_struct));
        connStruct_boot_pconn_ci = cell(length(cluster_struct),length(cluster_struct));
        subj_cl_ics=zeros(length(cluster_struct),length(cluster_struct));
        %##
        for subj_i=1:length(subj_inds)
            EEG = MAIN_ALLEEG(subj_inds(subj_i));
%             CAT = EEG.etc.COND_CAT(cond_i);
            %- laod statistics
            save_dir_bootmats = [save_dir filesep CONN_MEASURES{conn_i} filesep boot_folder];
            save_dir_nonzmats =  [save_dir filesep CONN_MEASURES{conn_i} filesep nonz_folder];
            BOOTSTRAP_STRUCT = par_load(save_dir_bootmats,sprintf('%s_boot_struct.mat',EEG.subject));
            PHASERND_STRUCT = par_load(save_dir_nonzmats,sprintf('%s_phasernd_struct.mat',EEG.subject));
%             boot_in = BOOTSTRAP_STRUCT.masked_conn{cond_i};
            %- baseline boot test
%             boot_in = BOOTSTRAP_STRUCT.stats{cond_i}.(CONN_MEASURES{conn_i}).pval;
%             boot_in = BOOTSTRAP_STRUCT.masked_conn{cond_i};
            boot_in = BOOTSTRAP_STRUCT.averages{cond_i};
            %- Hab boot test
%             boot_in = BOOTSTRAP_STRUCT.stats{HAB_BOOT_INDS(1,1),HAB_BOOT_INDS(1,2)}.(CONN_MEASURES{conn_i}).pval;
%             boot_in = BOOTSTRAP_STRUCT.masked_conn{HAB_BOOT_INDS(1,1),HAB_BOOT_INDS(1,2)};
%             boot_in = BOOTSTRAP_STRUCT.averages{HAB_BOOT_INDS(1,1),HAB_BOOT_INDS(1,2)};
            %- nonzero test
            nonz_in = PHASERND_STRUCT.masked_conn{cond_i};
%             nonz_in = PHASERND_STRUCT.stats{cond_i}.(CONN_MEASURES{conn_i}).pval;
            %- determine IC to Cluster assignments
            icaDat2Add=[]; icaEst4_crossVal=[]; ic_nums=[];
            for j=1:length(cluster_ints)
                if ~any(cluster_struct(cluster_ints(j)).sets==subj_inds(subj_i))
                else
                    ic_nums=[ic_nums j];
                end
            end
%             disp(ic_nums)
            %- assign data to cells and save
            subj_cl_ics(ic_nums,ic_nums)=subj_cl_ics(ic_nums,ic_nums)+1;
            for j=1:length(ic_nums)
                for k=1:length(ic_nums)
                    fprintf('%s) Assiging edge %i->%i\n',EEG.subject,ic_nums(j),ic_nums(k))
                    %- option 1
%                     connStruct_boot{ic_nums(j),ic_nums(k)}(subj_cl_ics(ic_nums(j),ic_nums(k)),:,:)=squeeze(nonz_in(j,k,:,:));
                    %- option 2
%                     connStruct_boot{ic_nums(j),ic_nums(k)}(numConns_boot(ic_nums(j),ic_nums(k)),:,:)=squeeze(CAT.Conn.(CONN_MEASURES{conn_i})(j,k,:,:));
                    %- option 3
%                     connStruct_boot{ic_nums(j),ic_nums(k)}(numConns_boot(ic_nums(j),ic_nums(k)),:,:)=squeeze(boot_in(j,k,:,:));
                    %- option 4
                    connStruct_boot{ic_nums(j),ic_nums(k)}(subj_cl_ics(ic_nums(j),ic_nums(k)),:,:)=real(squeeze(boot_in(j,k,fAllInds,tInds)).*squeeze(PHASERND_STRUCT.stats{cond_i}.(CONN_MEASURES{conn_i}).pval(j,k,fAllInds,tInds)));
                    %- mask saves
%                     connStruct_boot_pconn{ic_nums(j),ic_nums(k)}(subj_cl_ics(ic_nums(j),ic_nums(k)),:,:)=squeeze(BOOTSTRAP_STRUCT.stats{HAB_BOOT_INDS(1,1),HAB_BOOT_INDS(1,2)}.(CONN_MEASURES{conn_i}).pval(j,k,:,:));
%                     connStruct_boot_pconn{ic_nums(j),ic_nums(k)}(numConns_boot(ic_nums(j),ic_nums(k)),:,:)=squeeze(BOOTSTRAP_STRUCT.stats{cond_i}.(CONN_MEASURES{conn_i}).pval(j,k,:,:));
%                     connStruct_nonz_pconn{ic_nums(j),ic_nums(k)}(subj_cl_ics(ic_nums(j),ic_nums(k)),:,:)=squeeze(PHASERND_STRUCT.stats{cond_i}.(CONN_MEASURES{conn_i}).pval(j,k,:,:));
    %                 disp('Hi!');
                end
            end
        end
        par_save(connStruct_boot,fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
%         par_save(connStruct_boot_pconn,fpath,sprintf('connStruct%s_bootpconn.mat',CONN_MEASURES{conn_i}));
%         par_save(connStruct_nonz_pconn,fpath,sprintf('connStruct%s_nonzpconn.mat',CONN_MEASURES{conn_i}));
        %###
        % Set basetime to NaN if you don't want to significant mask
        % Otherwise set basetime to the limits of your cycle, ex. stride
        
        baselines=[];
        connStruct = zeros(length(cluster_struct),length(cluster_struct),length(fAllInds),length(alltimes));
        for j=1:length(cluster_ints)
            for k=1:length(cluster_ints)
                %## BASELINE
                %- Calculate and subtract baseline
                baseidx=find(alltimes>=BASELINE_TIME(1) & alltimes<=BASELINE_TIME(2));
                baseVals = median(connStruct_boot{j,k}(:,:,baseidx),3);
%                 baseVals = mean(connStruct_boot{j,k}(:,:,baseidx),3);
                curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(alltimes)]);
                curr_ersp = permute(curr_ersp,[2 3 1]);
                %- Use bootstat & bootstrap and significance mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',ALPHA,'dimaccu',2);
%                 pboot = bootstat(curr_ersp,'mean(arg1,3);','boottype','shuffle',...
%                     'label','ERSP','bootside','both','naccu',200,...
%                     'basevect',baseidx,'alpha',ALPHA,'dimaccu',2);
%                 [p_fdr,p_mask] = fdr(pboot,ALPHA,'parametric');
                fprintf('\n');
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
                %- store
                curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
                connStruct(j,k,:,:)=squeeze(curr_maskedersp);
                %## NO BASELINE (SIFT STATS)
%                 curr_ersp = permute(connStruct_boot{j,k},[2 3 1]);
%                 %- use SIFT output (either boot, nonz, or both), must be stats!
% %                 mask_in = permute(connStruct_nonz_pconn{j,k}<ALPHA,[2 3 1]);
%                 mask_in = permute(connStruct_boot_pconn{j,k}<ALPHA,[2 3 1]);
%                 fprintf('Significant TxF points: %i\n',sum(mask_in,[1,2,3]));
%                 curr_maskedersp = curr_ersp;
%                 for subj_i = 1:size(curr_ersp,3)
%                     curr_maskedersp(:,:,subj_i) = curr_maskedersp(:,:,subj_i).*mask_in(:,:,subj_i);
%                 end
%                 fprintf('TxF masked sum: %i\n',sum(curr_maskedersp,[1,2,3]));
%                 %- average subjects (median)
%                 curr_maskedersp = median(curr_maskedersp,3);
%                 %- store
%                 curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
%                 connStruct(j,k,:,:)=squeeze(curr_maskedersp);
                
            end
        end
        par_save(connStruct,fpath,sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
        
    end
end
%% NO BASELINE (Condition Difference Bootstat)
%{
for conn_i = 1:length(CONN_MEASURES)
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data'];
    if ~exist([fpath filesep 'btwn_cond'],'dir')
        mkdir([fpath filesep 'btwn_cond']);
    end
    connStruct = zeros(length(cluster_struct),length(cluster_struct),length(allfreqs),length(alltimes));
    diff_curr_ersp = cell(length(cluster_struct),length(cluster_struct));
    connStruct_boot_1 = par_load([fpath filesep COND_NAMES{1}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
    connStruct_boot_2 = par_load([fpath filesep COND_NAMES{2}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
    for j=1:length(cluster_ints)
        for k=1:length(cluster_ints)
            %- Calculate and subtract baseline
            diff_curr_ersp{j,k} = connStruct_boot_1{j,k}-connStruct_boot_2{j,k};
            curr_ersp = permute(diff_curr_ersp{j,k},[2 3 1]);
            %- Use bootstat & bootstrap and significance mask
            pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                'label','ERSP','bootside','both','naccu',200,...
                'alpha',ALPHA,'dimaccu',2);
            fprintf('\n');
            curr_ersp = median(curr_ersp,3);
            curr_maskedersp = curr_ersp;
            curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
            %- store
            curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
            connStruct(j,k,:,:)=squeeze(curr_maskedersp);
        end
    end
    par_save(diff_curr_ersp,[fpath filesep 'btwn_cond'],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
    par_save(connStruct,[fpath filesep 'btwn_cond'],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
end
%}
%%
%Create theta and alpha band networks (4-8, 8-13 Hz), averaging 1st second
%of activity; then do statistical testing using bootstrap distributions
%with same significance mask as average, this will determine which edges
%are significantly different from the other conditions
%-
EEG = MAIN_ALLEEG(1);
CONN_DIMS = length(cluster_ints);
NET_VALS = struct('theta',{cell(CONN_DIMS,CONN_DIMS)},...
    'alpha',{cell(CONN_DIMS,CONN_DIMS)},...
    'beta',{cell(CONN_DIMS,CONN_DIMS)},...
    'all',{cell(CONN_DIMS,CONN_DIMS)});
NETVALS_AVE = struct('theta',{zeros(CONN_DIMS,CONN_DIMS)},...
    'alpha',{zeros(CONN_DIMS,CONN_DIMS)},...
    'beta',{zeros(CONN_DIMS,CONN_DIMS)},...
    'all',{zeros(CONN_DIMS,CONN_DIMS)});
CAT = EEG.etc.COND_CAT(1);
tInds=find(alltimes>=TIME_LIMS_AVE(1) & alltimes<=TIME_LIMS_AVE(2));
fThetaInds=find(CAT.Conn.freqs>=FREQ_BANDS.theta(1) & CAT.Conn.freqs<=FREQ_BANDS.theta(end));
fAlphaInds=find(CAT.Conn.freqs>=FREQ_BANDS.alpha(1) & CAT.Conn.freqs<=FREQ_BANDS.alpha(end));
fBetaInds=find(CAT.Conn.freqs>=FREQ_BANDS.beta(1) & CAT.Conn.freqs<=FREQ_BANDS.beta(end));
% fAllInds=find(CAT.Conn.freqs>=FREQ_CROP(1) & CAT.Conn.freqs<=FREQ_CROP(end));
%##
for conn_i = 1:length(CONN_MEASURES)
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data'];
    if ~exist(fpath,'dir')
        mkdir(fpath);
    end
    for cond_i=1:length(COND_NAMES)
        net_vals = NET_VALS;
        net_vals_ave = NETVALS_AVE;
        connStruct_boot = par_load([fpath filesep COND_NAMES{cond_i}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
        connStruct = par_load([fpath filesep COND_NAMES{cond_i}],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
%         connStruct_boot = par_load([fpath filesep 'btwn_cond'],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
%         connStruct = par_load([fpath filesep 'btwn_cond'],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
        
        connStruct = real(connStruct);
        for j=1:length(cluster_ints)
            for k=1:length(cluster_ints)
                %## BASELINE
                %- baseline using mean of period
                baseidx=find(alltimes>=BASELINE_TIME(1) & alltimes<=BASELINE_TIME(2));
                baseVals=mean(connStruct_boot{j,k}(:,:,baseidx),3);
                curr_ersp = real(connStruct_boot{j,k}-repmat(real(baseVals), [1, 1, length(alltimes)]));
                %- no baseline
%                 curr_ersp = real(connStruct_boot{j,k});
                %## EXTRACT FREQ BANDS & TIMES
                bootDat_theta=curr_ersp(:,fThetaInds,tInds);
                bootDat_alpha=curr_ersp(:,fAlphaInds,tInds);
                bootDat_beta=curr_ersp(:,fBetaInds,tInds);
                bootDat_all=curr_ersp(:,:,tInds);
                for subj_i=1:size(bootDat_theta,1)
                    A=bootDat_theta(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fThetaInds,tInds))==0)=0; %mask using average mask
                    net_vals.theta{j,k}=[net_vals.theta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fThetaInds,tInds);
                    net_vals_ave.theta(j,k)= mean(tmp_dat(:));

                    A=bootDat_alpha(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fAlphaInds,tInds))==0)=0; %mask using average mask
                    net_vals.alpha{j,k}=[net_vals.alpha{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fAlphaInds,tInds);
                    net_vals_ave.alpha(j,k)= mean(tmp_dat(:));

                    A=bootDat_beta(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fBetaInds,tInds))==0)=0; %mask using average mask
                    net_vals.beta{j,k}=[net_vals.beta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fBetaInds,tInds);
                    net_vals_ave.beta(j,k)= mean(tmp_dat(:));

                    A=bootDat_all(subj_i,:,:);
                    A(squeeze(connStruct(j,k,:,tInds))==0)=0; %mask using average mask
                    net_vals.all{j,k}=[net_vals.all{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,:,tInds);
                    net_vals_ave.all(j,k)= mean(tmp_dat(:));
                end
            end
        end
        par_save(net_vals,fpath,sprintf('netVals_%s_aveTime_sbjs.mat',strjoin(strsplit(COND_NAMES{cond_i},' '),'_')));
        par_save(net_vals_ave,fpath,sprintf('netVals_%s_aveTime.mat',strjoin(strsplit(COND_NAMES{cond_i},' '),'_')));
    end
end
%% PERMUTATE SUBJECTS
%Permuting subset of subjects to test the min, average, and max samples per connection.
%This helps estimate how increases in subject sample size affect the sample size at each connection.

% subjectInds=[1:7 9:13 15:17 19:20 22:33];
% numConns_boot=cell(16,16);
% cond='allSVZ'; group = 'all';
% RootData='/usr/local/VR_connectivity/Data/';
goodClusts=3:12;
sub_cl_inds=3:12;
% load('/usr/local/VR_connectivity/Data/STUDYtopo.mat');
% load('/usr/local/VR_connectivity/Data/cluster_1222018.mat'); %_pruned.mat');
cl_struct = MAIN_STUDY.cluster;
allSubjs = 1:length(MAIN_ALLEEG);
subj_cl_ics = cell(length(cl_struct),length(cl_struct));
for subj_i=1:length(MAIN_ALLEEG)
%     EEG = MAIN_ALLEEG(subj_i);
    %Now figure out which clusters this subject is missing
    icaDat2Add=[]; icaEst4_crossVal=[]; icNums=[];
    for j=1:length(goodClusts)
        if ~any(allSubjs(cl_struct(goodClusts(j)).sets)==subj_i)
        else
            icNums=[icNums j];
        end
    end
    
    icNums=[icNums 9:16]; %add in muscles
    for j=icNums
        for k=icNums
            subj_cl_ics{j,k}=[subj_cl_ics{j,k} subj_i];
        end
    end
end

%% From a range of 1-29 subjects, do 200 subject permutations and average the result
numPerms=200; 
finalPermMean=zeros(1,length(MAIN_ALLEEG)); finalPermStd=finalPermMean;
finalPermMean_min=finalPermMean; finalPermStd_min=finalPermMean;
finalPermMean_max=finalPermMean; finalPermStd_max=finalPermMean;

for numSubjs=1:length(MAIN_ALLEEG)
    permVals=zeros(1,numPerms);
    permVals_min=zeros(1,numPerms);
    permVals_max=zeros(1,numPerms);
    disp(numSubjs)
    for i=1:numPerms
        tmpInds=randperm(length(allSubjs),numSubjs);
        newSubjInds=allSubjs(tmpInds);
        
        %Check the number of connections for each off-diagonal and average
        %results
        numConns=[];
        for j=1:length(sub_cl_inds)
            for k=1:length(sub_cl_inds)
                if j~=k
                    tmpSum=0;
                    for p=newSubjInds
                        tmpSum=tmpSum+sum(p==subj_cl_ics{j,k});
                    end
                    numConns=[numConns tmpSum];
                end
            end
        end
        permVals(1,i)=mean(numConns);
        permVals_min(1,i)=min(numConns);
        permVals_max(1,i)=max(numConns);
    end
    finalPermMean(1,numSubjs)=mean(permVals); finalPermStd(1,numSubjs)=std(permVals);
    finalPermMean_min(1,numSubjs)=mean(permVals_min); finalPermStd_min(1,numSubjs)=std(permVals_min);
    finalPermMean_max(1,numSubjs)=mean(permVals_max); finalPermStd_max(1,numSubjs)=std(permVals_max);
end

%%
figure; hold on;
errorbar(1:length(MAIN_ALLEEG),finalPermMean,finalPermStd,'g');
errorbar(1:length(MAIN_ALLEEG),finalPermMean_min,finalPermStd_min,'r');
errorbar(1:length(MAIN_ALLEEG),finalPermMean_max,finalPermStd_max,'b');
hold off;
xlim([0 30]); ylim([-1 20]);

Fit_mean = polyfit(1:length(MAIN_ALLEEG),finalPermMean,2);
Fit_min = polyfit(1:length(MAIN_ALLEEG),finalPermMean_min,2);
Fit_max = polyfit(1:length(MAIN_ALLEEG),finalPermMean_max,2);
finalPerm.meanVal=finalPermMean;
finalPerm.maxVal=finalPermMean_max;
finalPerm.minVal=finalPermMean_min;
finalPerm.meanValstd=finalPermStd;
finalPerm.maxValstd=finalPermStd_max;
finalPerm.minValstd=finalPermStd_min;
%% Simulation of IC number vs. number of connections (should show quadratic behavior)
numConns=zeros(1,length(MAIN_ALLEEG));
for i=1:length(MAIN_ALLEEG)
    for j=1:i
        for k=1:i
            if j~=k
                numConns(i)=numConns(i)+1;
            end
        end
    end
end
figure; plot(1:length(MAIN_ALLEEG),numConns);
Fit_ICnum = polyfit(1:length(MAIN_ALLEEG),numConns,2);
finalPerm.simFitICnum=numConns;
finalPerm.x=1:length(MAIN_ALLEEG);

% save('/usr/local/VR_connectivity/Data/regularGroupConn/plots/numSubjsNumConns_fit.mat','finalPerm');