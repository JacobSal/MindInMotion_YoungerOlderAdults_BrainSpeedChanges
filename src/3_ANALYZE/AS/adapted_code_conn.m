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
%- ball machine vs human rally
% COND_CHARS =  {'2Bounce_Human','2Bounce_BM'};
% EVENT_CHARS = {'Subject_receive'}; 
%- server vs rally human
% COND_CHARS =  {'1Bounce_Human','Serve_Human'};
% EVENT_CHARS = {'Subject_hit'}; 
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
% dt = '01312023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75';
dt = '02012023_subjrec_2bounces_rally_serve_human_epochfix_JS_n1p5-1p5';
% dt = '02012023_subjrec_2bounces_rally_serve_human_epochfixfix_JS_n1p5-1p5';
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
EVENT_COND_COMBOS{cnt} = 'StandingBaseline';
%- path for local data
% study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
study_fName_1 = 'epoch_study';
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'conn'];
% save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'conn_data'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
if ~ispc
    addpath(convertPath2UNIX('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
else
    addpath(convertPath2Drive('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
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
%     EVENT_COND_COMBOS = MAIN_STUDY.etc.a_epoch_process.epoch_chars;
    %## CUT OUT NON VALID CLUSTERS
    valid_cls = [3,4,5,6,7,8,9,11,12];
    inds = setdiff(1:length(comps_out),valid_cls);
    comps_out(inds,:) = 0;
end
%%
% SUBJ_PICS = {{'02','03','04','05','09','11','15','16','18','19','21','22',...
%             '23','24','25','27','28','29','30','31','32','33','35','36','38'}};
% cluster_struct = MAIN_STUDY.urcluster;
% cluster_struct_orig = MAIN_STUDY.urcluster;
% subj_chars_orig = SUBJ_PICS{1};
% subj_chars_orig = cellfun(@(x) [{} sprintf('Pilot%s',x)],subj_chars_orig);
% subj_chars = {MAIN_STUDY.datasetinfo.subject};
% subj_keep = zeros(length(subj_chars),1);
% for subj_i = 1:length(subj_chars_orig)
%     disp(any(strcmp(subj_chars_orig{subj_i},subj_chars)));
%     if any(strcmp(subj_chars_orig{subj_i},subj_chars))
%         subj_keep(subj_i) = 1;
%     end
% end
% subjs_rmv = find(~subj_keep);
% subj_keep = find(subj_keep);
% [val,ord] = sort(subj_keep);
% for cli = 2:length(cluster_struct)
%     si = cluster_struct(cli).sets;
%     ci = cluster_struct(cli).comps;
%     keep_si = setdiff(si,subjs_rmv);
%     tmp = cluster_struct(cli).preclust.preclustdata;
%     tmp_preclust = [];
%     tmp_si = [];
%     tmp_ci = [];
%     for i = 1:length(keep_si)
%         tmp_si = [tmp_si, repmat(ord(keep_si(i) == val),1,sum(keep_si(i) == si))];
%         tmp_ci = [tmp_ci, ci(keep_si(i) == si)];
%         tmp_preclust = [tmp_preclust; tmp(keep_si(i) == si,:)];
%     end
%     cluster_struct(cli).sets = tmp_si;
%     cluster_struct(cli).comps = tmp_ci;
%     cluster_struct(cli).preclust.preclustdata = tmp_preclust;
% end
% cluster_struct(1).sets = [cluster_struct(2:end).sets];
% cluster_struct(1).comps = [cluster_struct(2:end).comps];
% MAIN_STUDY.cluster = cluster_struct;
% [MAIN_STUDY,MAIN_ALLEEG] = parfunc_save_study(MAIN_STUDY,MAIN_ALLEEG,...
%                         MAIN_STUDY.filename,MAIN_STUDY.filepath,...
%                         'RESAVE_DATASETS','off');
%%
%{
sub_ints = [3,6,7,8,9,12];
% mim_gen_cluster_figs(MAIN_STUDY,MAIN_ALLEEG,save_dir,...
%     'CLUSTERS_TO_PLOT',sub_ints);
DIP_SING_POS=[16 582 420 360];
% different clusters have different colors
MAIN_STUDY.etc.dipparams.centrline = 'off';
std_dipplot_CL(MAIN_STUDY,MAIN_ALLEEG,'clusters',sub_ints,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',DIP_SING_POS,'color','w')
camzoom(1.2^2);
% camzoom(1.2^2);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.eps')],'ContentType','vector','Resolution',300);
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
view([45,0,0])
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.eps')],'ContentType','vector','Resolution',300);
view([0,-45,0])
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.eps')],'ContentType','vector','Resolution',300);
%}
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
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\AS_dataset\_studies\12282023_bounces_1h_2bm_JS_n1-0p5\_figs\conn';
% 
destination_folder = save_dir;  %change as appropriate
if ~exist(destination_folder,'dir')
    mkdir(destination_folder);
end
%##
% CONN_MEASURES = {'S'}; %{'dDTF08','S'};
CONN_MEASURES = {'dDTF08'}; %{'dDTF08','S'};
FREQ_CROP = [4:50];
FREQ_BANDS = struct('theta',(4:8),...
    'alpha',(8:13),...
    'beta',(14:30),...
    'all',FREQ_CROP);
nonz_folder = 'nonzero_mat_files';
% boot_folder = 'bootstrap_baseline_mat_files';
boot_folder = 'bootstrap_mat_files';
HAB_BOOT_INDS = [2,1]; %[1,2];
% cluster_ints = [3,4,5,6,7,8,9,10,11,12];
% cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RSuppMotor','LPPa-Oc','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
% cluster_ints = [3,4,5,6,7,8,9,10,11,12];
% sub_ints = [3,4,5,7];
% cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','LTemp','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.
cluster_ints = [3,4,5,6,7,8,9,11,12];
% sub_ints = [3,4,5,7];
sub_ints = [3,6,7,8,9,12];
cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.

DO_FDR_BOOTSTAT = false;
REPEATED_MEAS = 1;
ALPHA = 0.05;
NUM_ITERS = 200;
SIFT_ALPHA = 0.95;
% BASELINE_TIME = [-0.75,-0.15];
% TIME_LIMS=[-0.75 0.75];
% TIME_LIMS_AVE=[0.15 0.75];
BASELINE_TIME = [-0.5,0];
TIME_LIMS=[-0.5 0.5];
TIME_LIMS_AVE=[0 0.5];
% do baseline 0 0.64
% cluster_ints = [3,7,5,4];
% cluster_names = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.
% COND_PAIR_ORDER = [1,2];
% save_dir_bootmats =  [save_dir filesep CONN_MEASURES{1} filesep 'bootstrap_baseline_mat_files'];
% BOOTSTRAP_STRUCT = par_load(save_dir_bootmats,sprintf('%s_boot_struct.mat',MAIN_ALLEEG(1).subject));
% CONN_MEAS = BOOTSTRAP_STRUCT.conn_meas;
% FREQ_BANDS = BOOTSTRAP_STRUCT.FREQ_BAND_INF.FREQ_BANDS;
% FREQ_NAMES = fieldnames(FREQ_BANDS);
allfreqs = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.freqs;
% alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
cluster_struct = MAIN_STUDY.cluster;
fAllInds=find(allfreqs>=FREQ_CROP(1) & allfreqs<=FREQ_CROP(end));
tInds=find(MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes>=TIME_LIMS(1) & MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes<=TIME_LIMS(2));
alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes;
alltimes = alltimes(tInds);
tInds_ave = find(alltimes>=TIME_LIMS_AVE(1) & alltimes<=TIME_LIMS_AVE(2));

%## LOAD TEMP EEG DATA
COND_NAMES = EVENT_COND_COMBOS; %MAIN_STUDY.etc.a_epoch_process.epoch_chars; %unique(EEG.etc.conn_table.t_fNames);
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
% [tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(save_dir,COND_CHARS);
[tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table([STUDIES_DIR filesep sprintf('%s',dt) filesep 'conn_data'],COND_CHARS);

%-
fprintf('HQ median information crit model order: %0.2f\n',median(tbl_summary_out.min_modorder_info_crit_hq_line));
fprintf('HQ iqr information crit model order: %0.2f\n',iqr(tbl_summary_out.min_modorder_info_crit_hq_line));
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
subj_inds = 1:length(MAIN_ALLEEG);
%-
% rej_subj = zeros(size(comps_out,2),1); 
% for subj_i = 1:size(comps_out,2)
%     tmp = comps_out(sub_ints,subj_i)>0;
%     if ~all(tmp)
%         fprintf('Reject: %i\n',subj_i);
%         rej_subj(subj_i) = 1;
%     end
% end
% subj_inds = find(~rej_subj);
conn_vals = cell(length(COND_NAMES),1);
%## 
conn_i = 1;
trial_nums = zeros(length(subj_inds),length(COND_NAMES));
for cond_i=1:length(COND_NAMES)
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data' filesep COND_NAMES{cond_i}];
    if ~exist(fpath,'dir')
        mkdir(fpath);
    end
    connStruct_boot=cell(length(cluster_struct),length(cluster_struct));
%     maskStruct_boot=cell(length(cluster_struct),length(cluster_struct));
%     connStruct_nonz_pconn = cell(length(cluster_struct),length(cluster_struct));
%     connStruct_boot_pconn = cell(length(cluster_struct),length(cluster_struct));
%     connStruct_boot_pconn_ci = cell(length(cluster_struct),length(cluster_struct));
    subj_cl_ics=zeros(length(cluster_struct),length(cluster_struct));
    %##
    for subj_i=1:length(subj_inds)
        EEG = MAIN_ALLEEG(subj_inds(subj_i));
        CAT = EEG.etc.COND_CAT(cond_i);
        fprintf('%s) Number of trials: %i',EEG.subject,CAT.trials);
        trial_nums(subj_i,cond_i) = CAT.trials;
        %- determine IC to Cluster assignments
        icaDat2Add=[]; icaEst4_crossVal=[]; ic_nums=[]; cl_num=[];
        for j=1:length(cluster_ints)
            if ~any(cluster_struct(cluster_ints(j)).sets==subj_inds(subj_i))
            else
                ic_nums=[ic_nums j];
                cl_num=[cl_num cluster_ints(j)];
            end
        end
%             disp(ic_nums)
        %- assign data to cells and save
        subj_cl_ics(ic_nums,ic_nums)=subj_cl_ics(ic_nums,ic_nums)+1;
        for j=1:length(ic_nums)
            for k=1:length(ic_nums)
                fprintf('\tAssiging edge %i->%i\n',ic_nums(j),ic_nums(k))
                connStruct_boot{ic_nums(j),ic_nums(k)}(subj_cl_ics(ic_nums(j),ic_nums(k)),:,:)=real(squeeze(CAT.Conn.(CONN_MEASURES{conn_i})(j,k,fAllInds,tInds)));

            end
        end
    end
    
    par_save(connStruct_boot,fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
end
%%
for cond_i = 1:length(COND_NAMES)
    fpath = [destination_folder filesep CONN_MEASURES{1} filesep 'R_data' filesep COND_NAMES{4}];
    base_conn = par_load(fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{1}));
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data' filesep COND_NAMES{cond_i}];
    connStruct_boot = par_load(fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
    for j=1:length(cluster_ints)
        for k=1:length(cluster_ints)
            basestand = median(base_conn{j,k}(:,:,:),3);
            connStruct_boot{j,k} = connStruct_boot{j,k}-repmat(basestand, [1, 1, length(alltimes)]);
        end
    end
    par_save(connStruct_boot,fpath,sprintf('connStruct%s_basecorr.mat',CONN_MEASURES{conn_i}));
end
%%
% fpath = [destination_folder filesep CONN_MEASURES{1} filesep 'R_data' filesep COND_NAMES{4}];
% base_conn = par_load(fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{1}));
for cond_i = 1:length(COND_NAMES)
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data' filesep COND_NAMES{cond_i}];
%     connStruct_boot = par_load(fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
    connStruct_boot = par_load(fpath,sprintf('connStruct%s_basecorr.mat',CONN_MEASURES{conn_i}));
    baselines=[];
    connStruct = zeros(length(cluster_struct),length(cluster_struct),length(fAllInds),length(alltimes));
    for j=1:length(cluster_ints)
        for k=1:length(cluster_ints)
            %## BASELINE
            %- Calculate and subtract baseline
            baseidx=find(alltimes>=BASELINE_TIME(1) & alltimes<=BASELINE_TIME(2));
            baseVals = median(connStruct_boot{j,k}(:,:,baseidx),3);
%             baseVals = median(base_conn{j,k}(:,:,:),3);
%                 baseVals = mean(connStruct_boot{j,k}(:,:,baseidx),3);
            if DO_FDR_BOOTSTAT
                %## GROUPSIFT BASELINE TEST (SUBTRACTION PERM)
%                 curr_ersp_compare = connStruct_boot{j,k};
                curr_ersp_compare = repmat(baseVals, [1, 1, length(alltimes)]);
                curr_ersp_compare = permute(curr_ersp_compare,[3 2 1]);
                curr_ersp = connStruct_boot{j,k}; %-repmat(baseVals, [1, 1, length(alltimes)]);
                curr_ersp = permute(curr_ersp,[3 2 1]);
%                 [mask,tscore,pval] = clusterLevelPermutationTest(curr_ersp(tInds_ave,:,:),curr_ersp_orig(tInds_ave,:,:),...
%                                     1,ALPHA,NUM_ITERS); 
                [mask,tscore,pval] = clusterLevelPermutationTest(curr_ersp,curr_ersp_compare,...
                                    REPEATED_MEAS,ALPHA,NUM_ITERS); 
                pval_fdr = fdr(pval);
                fprintf('\nPercent un-masked: %0.1f%%\n',100*(sum(pval_fdr<ALPHA,[1,2])/(size(pval_fdr,1)*size(pval_fdr,2))));
                curr_ersp = permute(curr_ersp,[2,1,3]);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(pval_fdr>=ALPHA) = 0;
%                     curr_maskedersp(squeeze(sum(maskStruct_boot{j,k},1))<2) = 0;
            else
                curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(alltimes)]);
                curr_ersp = permute(curr_ersp,[2 3 1]);
                %- Use bootstat & bootstrap and significance mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',2000,...
                    'basevect',baseidx,'alpha',ALPHA,'dimaccu',2);
                fprintf('\n');
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
            end
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
% base_conn_cond = zeros(length(MAIN_ALLEEG),length(COND_NAMES),length(cluster_ints),length(cluster_ints));
% base_corr_cond = zeros(length(MAIN_ALLEEG),length(COND_NAMES),length(cluster_ints),length(cluster_ints));
% base_conn_cond = cell(length(COND_NAMES),length(cluster_ints),length(cluster_ints));
% base_corr_cond = cell(length(COND_NAMES),length(cluster_ints),length(cluster_ints));
conn_vec = {cell(length(COND_NAMES),length(cluster_ints),length(cluster_ints))};
subj_vec = {cell(length(COND_NAMES),length(cluster_ints),length(cluster_ints))};
cond_vec = {cell(length(COND_NAMES),length(cluster_ints),length(cluster_ints))};
store_s = struct('subjects',subj_vec,...
    'condition',cond_vec,...
    'base_conn_alpha',conn_vec,...
    'base_corr_alpha',conn_vec,...
    'post_conn_alpha',conn_vec,...
    'base_conn_theta',conn_vec,...
    'base_corr_theta',conn_vec,...
    'post_conn_theta',conn_vec,...
    'base_conn_beta',conn_vec,...
    'base_corr_beta',conn_vec,...
    'post_conn_beta',conn_vec);

subjt = [];
condt = [];
valt_a1 = [];
valt_a2 = [];
valt_a3 = [];
valt_t1 = [];
valt_t2 = [];
valt_t3 = [];
valt_b1 = [];
valt_b2 = [];
valt_b3 = [];
jt = [];
kt = [];
%##
% fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data'];
% baseconn_boot = par_load([fpath filesep COND_NAMES{4}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{1}));
% baseconn = par_load([fpath filesep COND_NAMES{4}],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{1}));
for conn_i = 1:length(CONN_MEASURES)
    fpath = [destination_folder filesep CONN_MEASURES{conn_i} filesep 'R_data'];
    if ~exist(fpath,'dir')
        mkdir(fpath);
    end
    for cond_i=1:length(COND_NAMES)
        net_vals = NET_VALS;
        net_vals_ave = NETVALS_AVE;
        connStruct_boot = par_load(fpath,sprintf('connStruct%s_basecorr.mat',CONN_MEASURES{conn_i}));
%         connStruct_boot = par_load([fpath filesep COND_NAMES{cond_i}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
        connStruct = par_load([fpath filesep COND_NAMES{cond_i}],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
%         connStruct_boot = par_load([fpath filesep 'btwn_cond'],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
%         connStruct = par_load([fpath filesep 'btwn_cond'],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
        
        connStruct = real(connStruct);
        for j=1:length(cluster_ints)
            for k=1:length(cluster_ints)
                %## BASELINE
                %- baseline using mean of period
                baseidx=find(alltimes>=BASELINE_TIME(1) & alltimes<=BASELINE_TIME(2));
%                 baseVals=mean(connStruct_boot{j,k}(:,:,baseidx),3);
                baseVals = median(connStruct_boot{j,k}(:,:,baseidx),3);
%                 baseVals = median(baseconn_boot{j,k}(:,:,:),3);
%                 tmp = squeeze(median(baseVals,1));
%                 baseVals = median(base_conn{j,k}(:,:,:),3);
                curr_ersp = real(connStruct_boot{j,k}-repmat(real(baseVals), [1, 1, length(alltimes)]));
%                 baseVals = median(curr_ersp(:,:,baseidx),3);
%                 curr_ersp = real(curr_ersp-repmat(real(baseVals), [1, 1, length(alltimes)]));
                non_base_ersp = connStruct_boot{j,k};
                %- no baseline
%                 curr_ersp = real(connStruct_boot{j,k});
                %## EXTRACT FREQ BANDS & TIMES
                bootDat_theta=curr_ersp(:,fThetaInds,tInds);
                bootDat_alpha=curr_ersp(:,fAlphaInds,tInds);
                bootDat_beta=curr_ersp(:,fBetaInds,tInds);
                bootDat_all=curr_ersp(:,:,tInds);
                
                for subj_i=1:size(bootDat_theta,1)
                    %-
                    tt = squeeze(baseVals(subj_i,:));
                    ttt = connStruct_boot{j,k};
                    ttt = squeeze(ttt(subj_i,:,tInds));
                    
%                     base_conn_cond(subj_i,cond_i,j,k) = mean(tmptmp(:));
%                     base_corr_cond(subj_i,cond_i,j,k) = mean(bootDat_all(:));
%                     base_conn_cond{cond_i,j,k} = [base_conn_cond{cond_i,j,k}, mean(tmptmp(:))];
%                     base_corr_cond{cond_i,j,k} = [base_corr_cond{cond_i,j,k}, mean(bootDat_all(:))];
%                     base_conn_cond{cond_i,j,k} = [base_conn_cond{cond_i,j,k}, mean(tt(fBetaInds))];
%                     base_corr_cond{cond_i,j,k} = [base_corr_cond{cond_i,j,k}, mean(bootDat_beta(:))];
%                     store_s.subjects{cond_i,j,k} = [store_s.subjects{cond_i,j,k}, subj_i];
%                     store_s.condition{cond_i,j,k} = [store_s.condition{cond_i,j,k}, cond_i];
%                     store_s.base_conn_alpha{cond_i,j,k} = [store_s.base_conn_alpha{cond_i,j,k}, mean(tt(fAlphaInds))];
%                     store_s.base_corr_alpha{cond_i,j,k} = [store_s.base_corr_alpha{cond_i,j,k}, mean(bootDat_alpha(:))];
%                     store_s.post_conn_alpha{cond_i,j,k} = [store_s.post_conn_alpha{cond_i,j,k}, mean(ttt(fAlphaInds,:),'all')];
%                     store_s.base_conn_theta{cond_i,j,k} = [store_s.base_conn_theta{cond_i,j,k}, mean(tt(fThetaInds))];
%                     store_s.base_corr_theta{cond_i,j,k} = [store_s.base_corr_theta{cond_i,j,k}, mean(bootDat_theta(:))];
%                     store_s.post_conn_theta{cond_i,j,k} = [store_s.post_conn_theta{cond_i,j,k}, mean(ttt(fThetaInds,:),'all')];
%                     store_s.base_conn_beta{cond_i,j,k} = [store_s.base_conn_beta{cond_i,j,k}, mean(tt(fBetaInds))];
%                     store_s.base_corr_beta{cond_i,j,k} = [store_s.base_corr_beta{cond_i,j,k}, mean(bootDat_beta(:))];
%                     store_s.post_conn_beta{cond_i,j,k} = [store_s.post_conn_beta{cond_i,j,k}, mean(ttt(fBetaInds,:),'all')];
%                     %-
%                     subjt = [subjt; subj_i];
%                     condt = [condt; cond_i];
%                     jt = [jt; j];
%                     kt = [kt; k];
%                     valt_a1 = [valt_a1; mean(tt(fAlphaInds))];
%                     valt_a2 = [valt_a2; mean(bootDat_alpha(:))];
%                     valt_a3 = [valt_a3; mean(ttt(fAlphaInds,:),'all')];
%                     valt_t1 = [valt_t1; mean(tt(fThetaInds))];
%                     valt_t2 = [valt_t2; mean(bootDat_theta(:))];
%                     valt_t3 = [valt_t3; mean(ttt(fThetaInds,:),'all')];
%                     valt_b1 = [valt_b1; mean(tt(fBetaInds))];
%                     valt_b2 = [valt_b2; mean(bootDat_beta(:))];
%                     valt_b3 = [valt_b3; mean(ttt(fBetaInds,:),'all')];
                    %-
                    A=bootDat_theta(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fThetaInds,tInds))==0)=0; %mask using average mask
                    net_vals.theta{j,k}=[net_vals.theta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fThetaInds,tInds);
                    net_vals_ave.theta(j,k)= mean(tmp_dat(:));
                    
                    t1 = non_base_ersp(subj_i,fThetaInds,tInds);
                    t2 = non_base_ersp(subj_i,fThetaInds,baseidx);
                    store_s.base_conn_theta{cond_i,j,k} = [store_s.base_conn_theta{cond_i,j,k}, mean(t2(:))];
                    store_s.base_corr_theta{cond_i,j,k} = [store_s.base_corr_theta{cond_i,j,k}, mean(squeeze(mean(squeeze(A),1)))];
                    store_s.post_conn_theta{cond_i,j,k} = [store_s.post_conn_theta{cond_i,j,k}, mean(t1(:))];
                    valt_t1 = [valt_t1; mean(t2(:))];
                    valt_t2 = [valt_t2; mean(squeeze(mean(squeeze(A),1)))];
                    valt_t3 = [valt_t3; mean(t1(:))];
                    

                    A=bootDat_alpha(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fAlphaInds,tInds))==0)=0; %mask using average mask
                    net_vals.alpha{j,k}=[net_vals.alpha{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fAlphaInds,tInds);
                    net_vals_ave.alpha(j,k)= mean(tmp_dat(:));
                    
                    t1 = non_base_ersp(subj_i,fAlphaInds,tInds);
                    t2 = non_base_ersp(subj_i,fAlphaInds,baseidx);
                    store_s.base_conn_alpha{cond_i,j,k} = [store_s.base_conn_alpha{cond_i,j,k}, mean(t2(:))];
                    store_s.base_corr_alpha{cond_i,j,k} = [store_s.base_corr_alpha{cond_i,j,k}, mean(squeeze(mean(squeeze(A),1)))];
                    store_s.post_conn_alpha{cond_i,j,k} = [store_s.post_conn_alpha{cond_i,j,k}, mean(t1(:))];
                    valt_a1 = [valt_a1; mean(t2(:))];
                    valt_a2 = [valt_a2; mean(squeeze(mean(squeeze(A),1)))];
                    valt_a3 = [valt_a3; mean(t1(:))];

                    A=bootDat_beta(subj_i,:,:);
                    A(squeeze(connStruct(j,k,fBetaInds,tInds))==0)=0; %mask using average mask
                    net_vals.beta{j,k}=[net_vals.beta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,fBetaInds,tInds);
                    net_vals_ave.beta(j,k)= mean(tmp_dat(:));
                    
                    t1 = non_base_ersp(subj_i,fBetaInds,tInds);
                    t2 = non_base_ersp(subj_i,fBetaInds,baseidx);
                    store_s.base_conn_beta{cond_i,j,k} = [store_s.base_conn_beta{cond_i,j,k}, mean(t2(:))];
                    store_s.base_corr_beta{cond_i,j,k} = [store_s.base_corr_beta{cond_i,j,k}, mean(squeeze(mean(squeeze(A),1)))];
                    store_s.post_conn_beta{cond_i,j,k} = [store_s.post_conn_beta{cond_i,j,k}, mean(t1(:))];
                    valt_b1 = [valt_b1; mean(t2(:))];
                    valt_b2 = [valt_b2; mean(squeeze(mean(squeeze(A),1)))];
                    valt_b3 = [valt_b3; mean(t1(:))];

                    A=bootDat_all(subj_i,:,:);
                    A(squeeze(connStruct(j,k,:,tInds))==0)=0; %mask using average mask
                    net_vals.all{j,k}=[net_vals.all{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                    tmp_dat = connStruct(j,k,:,tInds);
                    net_vals_ave.all(j,k)= mean(tmp_dat(:));
                    
                    store_s.subjects{cond_i,j,k} = [store_s.subjects{cond_i,j,k}, subj_i];
                    store_s.condition{cond_i,j,k} = [store_s.condition{cond_i,j,k}, cond_i];
                    
                    
                    
                    %-
                    subjt = [subjt; subj_i];
                    condt = [condt; cond_i];
                    jt = [jt; j];
                    kt = [kt; k];
                    
                    
                    
                end
            end
        end
        par_save(net_vals,fpath,sprintf('netVals_%s_aveTime_sbjs.mat',strjoin(strsplit(COND_NAMES{cond_i},' '),'_')));
        par_save(net_vals_ave,fpath,sprintf('netVals_%s_aveTime.mat',strjoin(strsplit(COND_NAMES{cond_i},' '),'_')));
    end
end

%% corrected conn validation
sub_ints = [7,1,5,6,9,4];
% sub_ints = [1,2,3,4,5,6,7,8,9];
cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.
COLOR_LIMITS = [-0.00001,0.00001];
for cond_i = 1:length(COND_NAMES)
    datain = par_load(fpath,sprintf('netVals_%s_aveTime.mat',strjoin(strsplit(COND_NAMES{cond_i},' '),'_')));
%     tmp = squeeze(nanmedian(base_conn_cond(:,cond_i,sub_ints,sub_ints),1));
%         tmp = squeeze(base_corr_cond(cond_i,sub_ints,sub_ints));
%     tmp = squeeze(median([base_conn_cond(cond_i,sub_ints,sub_ints)]));
    fn = fieldnames(datain);
    for i = 1:length(fn)
        tmp = datain.(fn{i});
        tmp = tmp(sub_ints,sub_ints);
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
        %- plot
    %     figure;
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        %- create title
        hnd.YDisplayLabels = cluster_names(sub_ints);
        hnd.XDisplayLabels = cluster_names(sub_ints);
        hnd.ColorLimits = COLOR_LIMITS;
        hnd.GridVisible = 'off';
        hnd.FontName = 'Times New Roman';
        hnd.FontSize = 16;
        hnd.CellLabelFormat = '%0.1g';
        hnd.NodeChildren(3).Title.Interpreter = 'none';
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[10,100,720,620])
        exportgraphics(fig_i,[destination_folder filesep sprintf('%s_condheatmap_%s.jpg',COND_NAMES{cond_i},fn{i})],'Resolution',300);
        close(fig_i);
    end
end
%% corrected conn validation
sub_ints = [7,1,5,6,9,4];
% sub_ints = [1,2,3,4,5,6,7,8,9];
cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.
% COLOR_LIMITS = [-0.0001,0.0001];
COLOR_LIMITS = [0,1];
% base_conn_cond(base_conn_cond==0) = nan();
fn = fieldnames(store_s);
fn = fn(3:end);
% fn = {'base_corr_alpha','base_corr_theta','base_corr_beta'};
for data_i = 1:length(fn)
    for cond_i = 1:4
    %     tmp = squeeze(nanmedian(base_conn_cond(:,cond_i,sub_ints,sub_ints),1));
%         tmp = squeeze(base_corr_cond(cond_i,sub_ints,sub_ints));
        tmp = squeeze(store_s.(fn{data_i})(cond_i,sub_ints,sub_ints));
        tmp = cellfun(@median,tmp)*10^4;
    %     tmp = squeeze(median([base_conn_cond(cond_i,sub_ints,sub_ints)]));
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
        %- plot
    %     figure;
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        %- create title
        hnd.YDisplayLabels = cluster_names(sub_ints);
        hnd.XDisplayLabels = cluster_names(sub_ints);
        hnd.ColorLimits = COLOR_LIMITS;
        hnd.GridVisible = 'off';
        hnd.FontName = 'Times New Roman';
        hnd.FontSize = 16;
        hnd.CellLabelFormat = '%0.2g';
        hnd.NodeChildren(3).Title.Interpreter = 'none';
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'Position',[10,100,720,620])
        exportgraphics(fig_i,[destination_folder filesep sprintf('%s_condheatmap_%s.jpg',COND_NAMES{cond_i},fn{data_i})],'Resolution',300);
        close(fig_i);
        pause(1);
    end
end
%% GLM TEST
% subjt = categorical(subjt);
subjt = categorical(subjt);
condt = categorical(condt);
jt = categorical(jt);
kt = categorical(kt);
valt_a1 = double(valt_a1);
valt_a2 = double(valt_a2);
valt_a3 = double(valt_a3);
valt_t1 = double(valt_t1);
valt_t2 = double(valt_t2);
valt_t3 = double(valt_t3);
valt_b1 = double(valt_b1);
valt_b2 = double(valt_b2);
valt_b3 = double(valt_b3);
cond_table = table(subjt,condt,jt,kt,(valt_a1),(valt_a2),(valt_a3),...
    (valt_t1),(valt_t2),(valt_t3),(valt_b1),(valt_b2),(valt_b3));
table_vars = cond_table.Properties.VariableNames;
table_vars = table_vars(5:end);
for i = 1:length(table_vars)
    %## (MAIN TEST) CONDITION MIXED EFFECTS
%     modelspec = sprintf('%s~1+condt+subjt',table_vars{i});
%     mdl_one = fitlm(cond_table,modelspec);
%     modelspec = sprintf('%s~1+condt',table_vars{i});
%     mdl_one = fitlm(cond_table,modelspec);
%     modelspec = sprintf('%s~1+condt+(condt|subjt)',table_vars{i});
    modelspec = sprintf('%s~1+condt+jt*kt',table_vars{i});
    mdl_one = fitlme(cond_table,modelspec);
    modelspec = sprintf('%s~1',table_vars{i});
    mdl_comp = fitlme(cond_table,modelspec);
    [stats] = anova(mdl_one);
    terrain_mixc_f = stats{2,5};
    disp(mdl_one)
    disp(anova(mdl_one));
    fid = fopen([save_dir filesep sprintf('%s_mdl_fitlme.txt',table_vars{i})],'wt');
    %- converst Summary anova to char array
    txt = evalc('anova(mdl_one)');
    fprintf(fid,'ANOVA EVAL\n');
    fprintf(fid,'%s',txt);
    fprintf(fid,'\n');
    %- Convert summary to char array
    fprintf(fid,'FITGLME EVAL\n');
    txt = evalc('mdl_one');
    fprintf(fid,'%s',txt);
%     for j = 1:length(vals_tn)
%         fprintf(fid,'%i: %s\n',j,unique(cond_table.cond_t(j));
%     end
    %- multiple comparisons
    for y = 1:numel(mdl_one.CoefficientNames)
        fprintf(fid,'Contrast position %i: %s\n', y, char(mdl_one.CoefficientNames{y}));
    end
    p_mcomp = [];
%     [pVal] = coefTest(mdl_one, [0,1,0,0]); fprintf(fid,'multicomp intercept:cat_1_2, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
%     [pVal] = coefTest(mdl_one, [0,0,1,0]); fprintf(fid,'multicomp intercept:cat_1_3, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
%     [pVal] = coefTest(mdl_one, [0,0,0,1]); fprintf(fid,'multicomp intercept:cat_1_4, F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
%     [pVal] = coefTest(mdl_one, [0,1,2,0]); fprintf(fid,'multicomp cat_1_2:cat_1_3,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
%     [pVal] = coefTest(mdl_one, [0,1,0,2]); fprintf(fid,'multicomp cat_1_2:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
%     [pVal] = coefTest(mdl_one, [0,0,1,2]); fprintf(fid,'multicomp cat_1_3:cat_1_4,   F(%i, %i)=%0.2f,\t p=%0.7f\n', df1, df2, F, pVal); p_mcomp = [p_mcomp pVal];
    
    [h,crit_p,adj_ci_cvrg,adj_p_speed] = fdr_bh(p_mcomp,0.05,'pdep','no');
    fprintf(fid,'False Discovery Corrected:\n');
    fprintf(fid,'%0.4f\n',adj_p_speed); fprintf(fid,'%i\n',h); fprintf(fid,'critical fdr_p: %0.4f\n',crit_p);
    %- cohens f2 test
    R21 = mdl_comp.SSR/mdl_comp.SST;
    R22 = mdl_one.SSR/mdl_one.SST; 
    cohens_f2 = (R22-R21)/(1-R22);
    fprintf(fid,'cohens f2: %0.4f\n',cohens_f2);
    fprintf('Cohens f2: %0.4f\n',cohens_f2);
    fclose(fid);
end
%% PERMUTATE SUBJECTS
%Permuting subset of subjects to test the min, average, and max samples per connection.
%This helps estimate how increases in subject sample size affect the sample size at each connection.
goodClusts=valid_cls;
sub_cl_inds=sub_ints;
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