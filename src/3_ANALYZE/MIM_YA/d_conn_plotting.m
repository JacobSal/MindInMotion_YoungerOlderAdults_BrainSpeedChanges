%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM/run_d_conn_plotting.sh

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
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_YA'];
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
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75', '1p0','flat','low','med','high'};
% TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
%- datetime override
dt = '07222023_MIM_YAN33_subset_prep_verified_gait_conn'; %'16022023';
%- connectiviy specific
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % AS (06/22/2023)
CONN_MEAS_ANLYZ = 'dDTF08';
ALPHA = 0.05;
CLUSTER_INF_FPATH = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_YAN33_subset_prep_verified_gait_conn\cluster\dipole_1_scalp_0_ersp_0_spec_0\14\cluster_update_14.mat';
%## soft define
%- path for local data
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs' filesep 'conn'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_INF_FPATH);
    else
        cluster_dir = convertPath2Drive(CLUSTER_INF_FPATH);
    end
%     fprintf('Using cluster information from...\n%s\n',CLUSTER_INF_FPATH);
    cluster_update = par_load(CLUSTER_INF_FPATH,[]);
    STUDY.cluster = cluster_update;
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
%%
%## CONN MAT PARAMS
CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12,13,14,15];
CLUSTER_ASSIGNMENTS = {'Cingulate','Cuneus','L Occipital','Postcentral','Frontal_Mid','Paracentral','R Frontal Sup','R Parietal Sup','L Insula','R Occipital','L Parietal Sup','Post Cingulum','L Postcentral','N/A'}; % (06/27/2023) JS, unsure on these as of yet.
% SUBJS_TEST = 1:length(ALLEEG);
%## VARIABLE INITIALIZATION
tmp = ALLEEG(1).etc.conn_table;
% CONN_METHODS = unique(tmp.t_conn_meas);
conn_conds = unique(tmp.t_fNames);
conn_comps = [tmp.t_conn_comps{1}];
idx = find(strcmp(CONN_MEAS_ANLYZ,tmp.t_conn_meas));
%- get unique frequencies tested
conn_freqs = tmp.t_conn_freqs(idx);
tmpsz = cellfun(@size,conn_freqs,'UniformOutput',false); tmpsz = cellfun(@max,conn_freqs); tmpsz = max(tmpsz);
store_freqs = zeros(length(conn_freqs),tmpsz);
for i = 1:length(conn_freqs)
    store_freqs(i,1:length(conn_freqs{i})) = conn_freqs{i};
end
uniq_freqs = unique(store_freqs,'rows');
%- empty storage for values
mat_out_nan = nan(length(STUDY.cluster),length(STUDY.cluster),length(ALLEEG),length(conn_conds),size(uniq_freqs,1));
%## LOOP
for subj_i = 1:length(ALLEEG)
    tmp = ALLEEG(subj_i).etc.conn_table;
%     conn_meas = unique(tmp.t_conn_meas);
    conn_conds = unique(tmp.t_fNames);
    conn_comps = [tmp.t_conn_comps{1}];
    idx = find(strcmp(CONN_MEAS_ANLYZ,tmp.t_conn_meas));
    mats = tmp.t_conn_mats(idx);
    cnt = 1;
    mat_nan = nan(length(STUDY.cluster),length(STUDY.cluster));
    chk = (comps_out(1:end,subj_i)>0);
    if any(chk)
        for freq_i = 1:size(uniq_freqs,1)
            for cond_i = 1:length(conn_conds)
                fprintf('%s) Frequencies: %sHz\n',ALLEEG(subj_i).subject,...
                    num2str(uniq_freqs(freq_i,uniq_freqs(freq_i,:)>0)));
                meanMat = squeeze(mats{cnt}(1,:,:));  %stdvMat = squeeze(mats(2,:,:)); %medMat = squeeze(mats(3,:,:));
                clust_idx = zeros(1,length(conn_comps));
                comp_idx = zeros(1,length(conn_comps));
                fprintf('%s) Number of absent connections: %i\n',ALLEEG(subj_i).subject,...
                    sum(isnan(meanMat(:))));
                for i = 1:length(conn_comps)
                    chk = find(conn_comps(i) == comps_out(1:end,subj_i));
                    if ~isempty(chk)
                        clust_idx(i) = chk;
                        comp_idx(i) = i;
                    end
                end
                clust_idx = clust_idx(clust_idx ~= 0);
                comp_idx = comp_idx(comp_idx ~= 0);
                if ~isempty(comp_idx)
                    val_in = meanMat(comp_idx,comp_idx);
                    val_in(isnan(val_in)) = 0;
                    mat_nan(clust_idx,clust_idx) = val_in;
                else
                    cnt = cnt + 1;
                    continue;
                end 
                %- store
                mat_out_nan(:,:,subj_i,cond_i,freq_i) = mat_nan;
%                 disp(cnt)
                cnt = cnt + 1;
            end
        end
    else
        cnt = cnt + 1;
        continue;
    end
end

%% CONNECTIVITY MATRICIES BASELINED
BASELINE_INT = 2;
for freq_i = 1:size(uniq_freqs,1)
%     baseline = nanmean(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    baseline = nanmedian(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    test_conds = setdiff(1:length(conn_conds),BASELINE_INT);
    for cond_i = test_conds
        tmp = uniq_freqs(freq_i,uniq_freqs(freq_i,:)>0);
        FREQS_SUBPATH = sprintf('%i-%i',...
                    tmp(1),tmp(end));
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cnt = 1;
        for i = 1:length(STUDY.cluster)
            N = length(STUDY.cluster(i).sets);
    %         clusterNames{i} = sprintf('(N=%i) Cluster %i',N,i);
            if any(i == CLUSTER_ITERS)
                %- assign name if available
                idx = find(i == CLUSTER_ITERS);
                clusterNames{idx} = sprintf('(N=%i) %s',N,CLUSTER_ASSIGNMENTS{(i == CLUSTER_ITERS)});
                cnt = cnt + 1;            
            end
        end
        %## PLOT
        %- assign matrix
%         tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3); %nanmean(mat_nan,3);
%         tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
        tmp = nanmedian(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
    %     tmp = log(tmp);
        %- delte unused clusters
        tmp = tmp(CLUSTER_ITERS,:);
        tmp = tmp(:,CLUSTER_ITERS);
        tmp_mat_nan = mat_nan;
        tmp_mat_nan = squeeze(tmp_mat_nan(CLUSTER_ITERS,:,:));
        tmp_mat_nan = squeeze(tmp_mat_nan(:,CLUSTER_ITERS,:));
        %- plot
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        title(sprintf('''%s''-''%s'' Connectivity mean across clusters',conn_conds{cond_i},conn_conds{BASELINE_INT}));
        hnd.YDisplayLabels = clusterNames;
        hnd.XDisplayLabels = clusterNames;
        hnd.ColorLimits = [-0.01,0.01];
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        fig_i = get(groot,'CurrentFigure');
%         saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMatbased_%s.fig','clusters',conn_conds{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMatbased_%s.jpg','clusters',conn_conds{cond_i})]);
%         close(fig_i)
    end
end
%% CONNECTIVITY MATRICIES
BASELINE_INT = 1;
for freq_i = 1:size(uniq_freqs,1)
    baseline = nanmean(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    for cond_i = 1:length(conn_conds)
%         sub_cond_i = cond_i + find(strcmp(load_trials{1},TRIAL_TYPES))-1; % offset for weirdness
        tmp = uniq_freqs(freq_i,uniq_freqs(freq_i,:)>0);
        FREQS_SUBPATH = sprintf('%i-%i',...
                    tmp(1),tmp(end));
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cnt = 1;
        for i = 1:length(STUDY.cluster)
            N = length(STUDY.cluster(i).sets);
    %         clusterNames{i} = sprintf('(N=%i) Cluster %i',N,i);
            if any(i == CLUSTER_ITERS)
                %- assign name if available
                idx = find(i == CLUSTER_ITERS);
                clusterNames{idx} = sprintf('(N=%i) %s',N,CLUSTER_ASSIGNMENTS{(i == CLUSTER_ITERS)});
                cnt = cnt + 1;            
            end
        end
        %## PLOT
        %- assign matrix
        tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3); %nanmean(mat_nan,3);
%         tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),1) - baseline;
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
    %     tmp = log(tmp);
        %- delte unused clusters
        tmp = tmp(CLUSTER_ITERS,:);
        tmp = tmp(:,CLUSTER_ITERS);
        tmp_mat_nan = mat_nan;
        tmp_mat_nan = squeeze(tmp_mat_nan(CLUSTER_ITERS,:,:));
        tmp_mat_nan = squeeze(tmp_mat_nan(:,CLUSTER_ITERS,:));
        %- plot
        figure;
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        title(sprintf('''%s'' Connectivity mean across clusters',conn_conds{cond_i}));
        hnd.YDisplayLabels = clusterNames;
        hnd.XDisplayLabels = clusterNames;
        hnd.ColorLimits = [0,0.05];
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        fig_i = get(groot,'CurrentFigure');
%         saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.fig','clusters',conn_conds{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.jpg','clusters',conn_conds{cond_i})]);
%         close(fig_i)
    end
end
%% TEST TIMEFREQ GRID
%{
STAT_CHAR = {'bootstrap','nonzeros'};
if ~exist([save_dir filesep 'conn_mats'],'dir')
    mkdir([save_dir filesep 'conn_mats'])
end
parfor subj_i = 1:length(ALLEEG)
% for subj_i = 1:length(ALLEEG)
    [conn_subj_out]=as_cnct_stat_valid(ALLEEG(subj_i),STUDY,...
        CLUSTER_ITERS,CLUSTER_ASSIGNMENTS,squeeze(comps_out(:,subj_i)),save_dir,...
        'ALPHA',ALPHA,...
        'CONN_METHODS',{CONN_MEAS_ANLYZ});
    par_save(conn_subj_out,[save_dir filesep 'conn_mats'],sprintf('%s_connmat.mat',ALLEEG(subj_i).subject));            
end
%}
%% NONZERO STATISTICS MASK GENERATION
if ~exist([save_dir filesep 'nz_bs_test'],'dir')
    mkdir([save_dir filesep 'nz_bs_test'])
end
if ~exist([save_dir filesep 'pr_conn_mats'],'dir')
    mkdir([save_dir filesep 'pr_conn_mats'])
end
parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/2))
    %- generate display names based on CLUSTER_ASSIGNMENTS
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    display_names = cell(length(tmpcl),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{idxcl(i)} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(idxcl);
%     cluster_inds = idxcl;
    %-
    [~,phasernd_conn] = cnctanl_nz_test(ALLEEG(subj_i),...
        'ALPHA',ALPHA,...
        'CONN_METHODS',{CONN_MEAS_ANLYZ},...
        'SAVE_DIR',[save_dir filesep 'nz_bs_test'],...
        'DISPLAYNAMES',display_names);
    fprintf('Saving %s...',ALLEEG(subj_i).subject);
    par_save(phasernd_conn,[save_dir filesep 'pr_conn_mats'],sprintf('%s_connmat.mat',ALLEEG(subj_i).subject));
end
%% LOAD THEN AVERAGE NONZEROED CONNECTIVITY MATRICIES
%## LOAD
meth_i = 1; % method iter
freq_dim = length(ALLEEG(1).etc.COND_CAT(1).Conn.freqs);
FREQ_BANDS = {1:freq_dim;1:7;7:12;12:28;28:48;48:60};
conn_store = nan(length(STUDY.cluster),length(STUDY.cluster),length(FREQ_BANDS),length(ALLEEG),length(ALLEEG(1).etc.COND_CAT));
for subj_i = 1:length(ALLEEG)
    conn_subj_out = par_load([save_dir filesep 'pr_conn_mats'],sprintf('%s_connmat.mat',ALLEEG(subj_i).subject));
    for freq_i = 1:length(FREQ_BANDS)
        for cond_i = 1:size(conn_subj_out,1)
            %- generate display names based on CLUSTER_ASSIGNMENTS
            comps = squeeze(comps_out(:,subj_i));
            [tmpcl,idxcl] = sort(comps);
            idxcl = idxcl(tmpcl~=0);
            display_names = cell(length(tmpcl),1);
            for i = 1:length(idxcl)
                if any(idxcl(i) == CLUSTER_ITERS)
                    display_names{idxcl(i)} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
                end
            end
            display_names = display_names(idxcl);
            cluster_inds = idxcl;
            %- nan mask
            tmp = conn_subj_out{cond_i,meth_i};
            tmp(tmp == 0) = nan();
            %- color limits handle
            %* sum across frequencies (recreate connectivity trace
            % previously decomposed using fourier transform)
            tmp = squeeze(tmp(:,:,FREQ_BANDS{freq_i},:));
            tmp = squeeze(nansum(tmp,3));
            %* average across time
            tmp = squeeze(nanmean(tmp,3));
%             tmp = squeeze(nanmedian(tmp,3));
            %- store
            conn_store(cluster_inds,cluster_inds,freq_i,subj_i,cond_i) = tmp;
        end
    end
end
%% 3D Box Plot
%{
%## ANOVAS AGGREGATE ALL
freq_i = 1;
unravel_out = squeeze(conn_store(:,:,:,:,freq_i));
unravel_out = unravel_out(:);
% unravel_out(isnan(unravel_out)) = 0;
% group_labs_org = cell(length(unravel_out),1);
% group_labs_subj = cell(length(unravel_out),1);
% group_labs_clust_ij = cell(length(unravel_out),1);
group_labs_org = zeros(length(unravel_out),1);
group_cli = zeros(length(unravel_out),1);
group_clj = zeros(length(unravel_out),1);
cnt = 1;
for cond_i = 1:size(conn_store,4)
    for subj_i = 1:size(conn_store,3)
        for clust_i = 1:size(conn_store,2)
            for clust_j = 1:size(conn_store,1)
                group_cli(cnt) = clust_i;
                group_clj(cnt) = clust_j;
                if clust_i == clust_j
                    unravel_out(cnt) = nan();
                end
                %- group by condition and subject
        %         group_labs{i} = sprintf('s%i_c%i',subj_i,cond_i);
                %- group by condition
%                 group_labs_org{cnt} = cond_i; %sprintf('c%i',cond_i);
%                 group_labs_subj{cnt} = sprintf('s%i',subj_i);
                %- group by component connections and condition
%                 group_labs_clust_ij{cnt} = str2double(sprintf('%i%i',clust_j,clust_i)); %sprintf('ci%i_cj%i',clust_j,clust_i);
                cnt = cnt+1;
            end
        end
    end
end
cond_i = 1;
boxPlot3D(unravel_out,group_cli,group_clj,[ 0.25 0.5 0.75])
%}
%% 3D BOXPLOTS
BASELINE_INT = 2;
for freq_i = 1:length(FREQ_BANDS)
    for cond_i = 1:size(conn_store,5)
        baseline = squeeze(conn_store(:,:,freq_i,:,BASELINE_INT));
        FREQS_SUBPATH = sprintf('%i-%i',...
                    FREQ_BANDS{freq_i}(1),FREQ_BANDS{freq_i}(end));
        save_sub_figs = [save_dir filesep 'box3d_plots' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        tmp_in = squeeze(conn_store(:,:,freq_i,:,cond_i));
        %- (1) baseline to average
        tmp_in = tmp_in - nanmean(baseline,3);
        %- (2) baseline per subject
%         tmp_in = tmp_in - baseline;
%         for subj_i = 1:size(tmp_in,3)
%             tmp_in(:,:,subj_i) = tmp_in(:,:,subj_i) - baseline(:,:,subj_i);
%         end
%         tmp_in = reshape(tmp_in,[size(tmp_in,3),size(tmp_in,1),size(tmp_in,2)]);
        
        %- delete unused clusters
        tmp_in = tmp_in(:,CLUSTER_ITERS,:);
        tmp_in = tmp_in(:,:,CLUSTER_ITERS);
        %- (PLOT)
        custom_boxPlot3D(tmp_in)
        %- plot edits
        fig_i = get(groot,'CurrentFigure');
        tmp = strsplit(conn_conds{cond_i},'_');
        tmp = strsplit(strjoin(tmp(2:end),' '),'.');
        title(sprintf('(%s) condition %s',FREQS_SUBPATH,tmp{1}));
        fig_i.Children(2).YTick = [1:length(CLUSTER_ASSIGNMENTS)];
        fig_i.Children(2).YTickLabel = CLUSTER_ASSIGNMENTS;
        fig_i.Children(2).YTickLabelRotation = 45;
        fig_i.Children(2).XTick = [1:length(CLUSTER_ASSIGNMENTS)];
        fig_i.Children(2).XTickLabel = CLUSTER_ASSIGNMENTS;
        fig_i.Children(2).XTickLabelRotation = 45;
%         saveas(fig_i,[save_sub_figs filesep sprintf('nonzero_3dbox_%i.fig',cond_i)]);
        saveas(fig_i,[save_sub_figs filesep sprintf('nonzero_3dbox_%i.jpg',cond_i)]);
    end
end
%% CONNECTIVITY MATRICIES
%## TIME
tic
%## CALC BOOTSTRAPPED MEAN
for freq_i = 1:length(FREQ_BANDS)
    for cond_i = 1:size(conn_store,5)
        %## SAVE_PATH
        FREQS_SUBPATH = sprintf('%i-%i',...
                    FREQ_BANDS{freq_i}(1),FREQ_BANDS{freq_i}(end));
        save_sub_figs = [save_dir filesep 'nz_cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cnt = 1;
        for i = 1:length(STUDY.cluster)
            N = length(STUDY.cluster(i).sets);
            if any(i == CLUSTER_ITERS)
                %- assign name if available
                idx = find(i == CLUSTER_ITERS);
                clusterNames{idx} = sprintf('(N=%i) %s',N,CLUSTER_ASSIGNMENTS{(i == CLUSTER_ITERS)});
                cnt = cnt + 1;            
            end
        end
        %## Extract
        tmp = squeeze(conn_store(:,:,freq_i,:,cond_i));
        %## average across subjects
%         tmp = squeeze(nanmean(tmp,3));
        tmp = squeeze(nanmedian(tmp,3));
        %## PLOT
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
    %     tmp = log(tmp);
        %- delte unused clusters
        tmp = tmp(CLUSTER_ITERS,:);
        tmp = tmp(:,CLUSTER_ITERS);
        %- plot
        figure('Color','w','Name',sprintf('Mean Nonzero Masked'));
        hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
        tmp = strsplit(conn_conds{cond_i},'_');
        tmp = strsplit(tmp{2},'.');
        title(sprintf('(%s) Nonzero Masked %s',FREQS_SUBPATH,tmp{1}));
        hnd.YDisplayLabels = clusterNames;
        hnd.XDisplayLabels = clusterNames;
        hnd.ColorLimits = [0,0.01];
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_sub_figs filesep sprintf('nonzero_TimeFreqChart_%i.fig',cond_i)]);
        saveas(fig_i,[save_sub_figs filesep sprintf('nonzero_TimeFreqChart_%i.jpg',cond_i)]);
%         close(fig_i)
    end
end
close all
%%
if ~exist([save_dir filesep 'bs_conn_mats'],'dir')
    mkdir([save_dir filesep 'bs_conn_mats'])
end
parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/2))
    %- generate display names based on CLUSTER_ASSIGNMENTS
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    display_names = cell(length(tmpcl),1);
    for i = 1:length(idxcl)
        if any(idxcl(i) == CLUSTER_ITERS)
            display_names{idxcl(i)} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
        end
    end
    display_names = display_names(idxcl);
%     cluster_inds = idxcl;
    %-
    [~,bootstrap_conn,~] = cnctanl_bs_test(ALLEEG(subj_i),...
        'ALPHA',ALPHA,...
        'CONN_METHODS',{CONN_MEAS_ANLYZ},...
        'SAVE_DIR',[save_dir filesep 'nz_bs_test'],...
        'DISPLAYNAMES',display_names);
    fprintf('Saving %s...',ALLEEG(subj_i).subject);
    par_save(bootstrap_conn,[save_dir filesep 'bs_conn_mats'],sprintf('%s_connmat.mat',ALLEEG(subj_i).subject));
end

%% AVERAGING BOOTSTRAPPED & NONZEROED CONNECTIVITY MATRICIES
STAT_CHARS = {'bootstrap'};
%## LOAD
stat_i = 1;
freq_dim = length(ALLEEG(1).etc.COND_CAT(1).Conn.freqs);
FREQ_BANDS = {1:freq_dim;1:7;7:12;12:28;28:48;48:60};
conn_store = nan(length(STUDY.cluster),length(STUDY.cluster),length(FREQ_BANDS),length(ALLEEG),length(ALLEEG(1).etc.COND_CAT),length(ALLEEG(1).etc.COND_CAT));
for subj_i = 1:length(ALLEEG)
    conn_subj_out = par_load([save_dir filesep 'bs_conn_mats'],sprintf('%s_connmat.mat',ALLEEG(subj_i).subject));
    for eeg_i = 1:2
        for freq_i = 1:length(FREQ_BANDS)
            done = [];
            for cond_i = 1:size(conn_subj_out,1)
                for cond_j = 1:size(conn_subj_out,2)
                    if any((cond_j==done)) || cond_i==cond_j
                        continue;
                    end
                    %- generate display names based on CLUSTER_ASSIGNMENTS
                    comps = squeeze(comps_out(:,subj_i));
                    [tmpcl,idxcl] = sort(comps);
                    idxcl = idxcl(tmpcl~=0);
                    display_names = cell(length(tmpcl),1);
                    for i = 1:length(idxcl)
                        if any(idxcl(i) == CLUSTER_ITERS)
                            display_names{idxcl(i)} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{(idxcl(i) == CLUSTER_ITERS)},comps(idxcl(i)));
                        end
                    end
                    display_names = display_names(idxcl);
                    cluster_inds = idxcl;
                    %- nan mask
                    tmp = conn_subj_out{cond_i,cond_j};
                    tmp(tmp == 0) = nan();
                    %- color limits handle
                    %* sum across frequencies (recreate connectivity trace
                    % previously decomposed using fourier transform)
                    tmp = squeeze(tmp(:,:,FREQ_BANDS{freq_i},:,eeg_i));
                    tmp = squeeze(nansum(tmp,3));
                    %* average across time
                    tmp = squeeze(nanmean(tmp,3));
%                     tmp = squeeze(nanmedian(tmp,3));
                    %- store
                    conn_store(cluster_inds,cluster_inds,freq_i,subj_i,cond_i,cond_j,eeg_i) = tmp;
                    done = [done cond_i];
                end
            end
        end
    end
end
% bs_struct = struct('dims',{'from_clusters','to_clusters','frequencies','condition_1','ocndition_2','eeg_masked'},...
%     'conditions',{conn_conds},...
%     'cluster_nums',CLUSTER_ITERS,...
%     'cluster_names',{CLUSTER_ASSIGNMENTS},...
%     'frequency_bands',uniq_freqs,...
%     'bootstrap_masked',mat_out_nan);
% par_save(bs_struct,save_dir,'all_subj_connectivity.mat');
%% CONNECTIVITY MATRICIES
%## TIME
tic
%## CALC BOOTSTRAPPED MEAN
done = [];
for eeg_i = 1:2
    for freq_i = 1:length(FREQ_BANDS)
        for cond_i = 1:size(conn_store,6)
            for cond_j = 1:size(conn_store,7)
                if any((cond_j == done))
                    continue;
                end
                %## SAVE_PATH
                FREQS_SUBPATH = sprintf('%i-%i',...
                            FREQ_BANDS{freq_i}(1),FREQ_BANDS{freq_i}(end));
                save_sub_figs = [save_dir filesep 'bs_nz_cluster_mats' filesep FREQS_SUBPATH];
                if ~exist(save_sub_figs,'dir')
                    mkdir(save_sub_figs);
                end
                %- plot
                cnt = 1;
                for i = 1:length(STUDY.cluster)
                    N = length(STUDY.cluster(i).sets);
                    if any(i == CLUSTER_ITERS)
                        %- assign name if available
                        idx = find(i == CLUSTER_ITERS);
                        clusterNames{idx} = sprintf('(N=%i) %s',N,CLUSTER_ASSIGNMENTS{(i == CLUSTER_ITERS)});
                        cnt = cnt + 1;            
                    end
                end
                %## Extract
                tmp = squeeze(conn_store(:,:,:,freq_i,cond_i,cond_j,eeg_i));
                %## average across subjects
%                 tmp = squeeze(nanmean(tmp,3));
                tmp = squeeze(nanmedian(tmp,3));
                %## PLOT
                I = eye(size(tmp));
                I = (I == 0);
                tmp = tmp.*I;
                tmp(tmp == 0) = nan();
            %     tmp = log(tmp);
                %- delte unused clusters
                tmp = tmp(CLUSTER_ITERS,:);
                tmp = tmp(:,CLUSTER_ITERS);
                %- plot
                figure;
                hnd = heatmap(tmp,'Colormap',jet); %,'CellLabelColor', 'None');
                title(sprintf('%s: ',STAT_CHARS{stat_i},conn_conds{cond_i}));
                hnd.YDisplayLabels = clusterNames;
                hnd.XDisplayLabels = clusterNames;
                hnd.ColorLimits = [0,0.05];
                hnd.GridVisible = 'off';
                hnd.CellLabelFormat = '%0.1g';
                fig_i = get(groot,'CurrentFigure');
                saveas(fig_i,[save_sub_figs filesep sprintf('%s_eeg%i_TimeFreqChart_%i-%i.jpg',STAT_CHARS{stat_i},eeg_i,cond_i,cond_j)]);
                saveas(fig_i,[save_sub_figs filesep sprintf('%s_eeg%i_TimeFreqChart_%i-%i.jpg',STAT_CHARS{stat_i},eeg_i,cond_i,cond_j)]);
                close(fig_i)
            end
        end
    end
end

%## TIME
toc
%% ANOVAS AGGREGATE ALL
freq_i = 1;
unravel_out = squeeze(mat_out_nan(:,:,:,:,freq_i));
unravel_out = unravel_out(:);
% unravel_out(isnan(unravel_out)) = 0;
group_labs_org = cell(length(unravel_out),1);
group_labs_subj = cell(length(unravel_out),1);
group_labs_clust_ij = cell(length(unravel_out),1);
cnt = 1;
for cond_i = 1:size(mat_out_nan,4)
    for subj_i = 1:size(mat_out_nan,3)
        for clust_i = 1:size(mat_out_nan,2)
            for clust_j = 1:size(mat_out_nan,1)
                %- group by condition and subject
        %         group_labs{i} = sprintf('s%i_c%i',subj_i,cond_i);
                %- group by condition
                group_labs_org{cnt} = sprintf('c%i',cond_i);
                group_labs_subj{cnt} = sprintf('s%i',subj_i);
                %- group by component connections and condition
                group_labs_clust_ij{cnt} = sprintf('ci%i_cj%i',clust_j,clust_i);
                cnt = cnt+1;
            end
        end
    end
end

%%
%- subselect component
% CLUST_I = 2;
% CLUST_J = 3;
for i = 1:length(CLUSTER_ASSIGNMENTS)
    clust_i = CLUSTER_ITERS(i);
    for j = 1:length(CLUSTER_ASSIGNMENTS)
        clust_j = CLUSTER_ITERS(j);
        %## LOOP PATHS
        figs_save_dir = [save_dir filesep 'anovas' filesep [conn_conds{:}] filesep sprintf('%s_%s',CLUSTER_ASSIGNMENTS{i},CLUSTER_ASSIGNMENTS{j})];
        if ~exist(figs_save_dir,'dir')
            mkdir(figs_save_dir)
        end
        %## LOOP MEAT
        comp_pick = sprintf('ci%i_cj%i',clust_i,clust_j);
        idx = strcmp(comp_pick,group_labs_clust_ij);
        unravel_in = unravel_out(idx);
        %- ANOVA-N analysis
        % group_labs = {group_labs_org,group_labs_subj,group_labs_clust_ij};
        group_labs = {group_labs_org(idx),group_labs_subj(idx)};
        % group_labs = {group_labs_org};
        [P,T,stats,terms] = anovan(unravel_in,group_labs,'display','on','alpha',0.05,'sstype',3);
        % [P,T,stats,terms] = anovan(unravel_out,group_labs,'display','on','alpha',0.05,'sstype',3);
%         writecell(T,[figs_save_dir filesep 'anovaresults.txt'])
        % Write string to file
        tblStr = cell(size(T,1),1);
        for col_i = 1:size(T,1)
            if col_i == 1
                tblStr{col_i} = sprintf('%-6s   %-9s   %-9s   %-9s   %-9s   %-9s   %-9s\n',T{col_i,:});
            else
                tblStr{col_i} = sprintf('%-6s   %-9.3g   %-9.0f   %-9.0f   %-9.3g   %-9.3g   %-9.3g\n',T{col_i,:});
            end
        end
        fid = fopen([figs_save_dir filesep 'anovaTable.txt'], 'wt');
        fileCleanup = onCleanup(@()fclose(fid));
        formatSpec = '%s\n';
        cellfun(@(x) fprintf(fid, formatSpec, x), tblStr)
%         fprintf(fid, formatSpec, tblStr);
        clear('fileCleanup')
        %- multiple comparissons test
%         figure();
%         hold on;
%         [c,m,h,gnames] = multcompare(stats);
%         hold off;
        %- Violin Plot;
        % violin_data = cell(1,length(load_trials));
        violin_data = nan(size(mat_out_nan,3),length(conn_conds));
        violin_group = cell(1,length(conn_conds));
        for cond_i = 1:length(conn_conds)
        %     violin_data{cond_i} = squeeze(mat_out_nan(CLUST_I,CLUST_J,:,cond_i,freq_i));
            violin_data(:,cond_i) = squeeze(mat_out_nan(clust_i,clust_j,:,cond_i,freq_i))';
            violin_group{cond_i} = conn_conds{cond_i};
        end
        figure;
        hold on
        violinplot(violin_data,violin_group)
        % violinplot(unravel_in,group_labs);
        % violinplot(unravel_out,group_labs);
        hold off;
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[figs_save_dir filesep sprintf('Violin_avg.fig')]);
        saveas(fig_i,[figs_save_dir filesep sprintf('Violin_avg.jpg')]);
    end
end
