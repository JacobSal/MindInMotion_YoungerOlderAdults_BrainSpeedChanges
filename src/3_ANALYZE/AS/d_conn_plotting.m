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
dt = '08182023_bounces_1h_2h_2bm_JS';
%- connectiviy specific
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % AS (06/22/2023)
CONN_MEAS_ANLYZ = 'dDTF08';
ALPHA = 0.05;
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
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
%%
%## CONN MAT PARAMS
% CLUSTER_ITERS = [2,3,4,5,6,7,8,9];
% CLUSTER_ASSIGNMENTS = {'RPPa','LPPa','PFC','LSM','ROc','ACC','RSM','Cerr'};
% CLUSTER_ITERS = [2,3,4,5,6,7,8,9];
% CLUSTER_ASSIGNMENTS = {'Cune','RSM','Cerr','ACC','RPPa','LPPa','LSM','LOc'};
CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12];
% CLUSTER_ITERS = [2,3,4,5,6,7,8,9,10,11];
CLUSTER_ASSIGNMENTS = {'RPPa','Cuneus','Precuneus','RSuppMotor','LPPa','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
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

%% MAKE CONDITIONS SUBJ STACKS
% if ~exist([save_dir filesep 'groupsift'],'dir')
%     mkdir([save_dir filesep 'groupsift'])
% end
% parfor subj_i = 1:length(ALLEEG)
%     EEG = ALLEEG(subj_i);
%     [EEG,stats_bootstrap,cond_bootstrap,stats_nonzero] = as_cnct_stat_test(EEG);
%     for cond_i = 1:length(EEG)
%         EEG{cond_i}.CAT.PConn = [];
%         pop_saveset(EEG(cond_i),'filename',EEG(cond_i).filename,'filepath',[save_dir filesep 'groupsift']);
%     end
% end
% tmp_struct = struct('dims',{'from_clusters','to_clusters','conditions','frequencies'},...
%     'conditions',{conn_conds},...
%     'cluster_nums',CLUSTER_ITERS,...
%     'cluster_names',{CLUSTER_ASSIGNMENTS},...
%     'frequency_bands',uniq_freqs,...
%     'connectivity_mat',mat_out_nan);
% par_save(tmp_struct,save_dir,'all_subj_connectivity.mat');
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
%% (NONZERO) STATISTICS MASK GENERATION
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
%% (NONZERO)LOAD & AVERAGE NONZEROED CONNECTIVITY MATRICIES
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
%% (NONZERO) CONNECTIVITY MATRICIES
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
%% 3D BOXPLOTS
BASELINE_INT = 1;
for freq_i = 1:length(FREQ_BANDS)
    for cond_i = 1:size(conn_store,5)
        baseline = squeeze(conn_store(:,:,freq_i,:,BASELINE_INT));
        FREQS_SUBPATH = sprintf('%i-%i',...
                    FREQ_BANDS{freq_i}(1),FREQ_BANDS{freq_i}(end));
        save_sub_figs = [save_dir filesep 'nz_cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        tmp_in = squeeze(conn_store(:,:,freq_i,:,cond_i));
        tmp_in = tmp_in - baseline;
        tmp_in = reshape(tmp_in,[size(tmp_in,3),size(tmp_in,1),size(tmp_in,2)]);
        
        %- delete unused clusters
        tmp_in = tmp_in(:,CLUSTER_ITERS,:);
        tmp_in = tmp_in(:,:,CLUSTER_ITERS);
        %- (PLOT)
        custom_boxPlot3D(tmp_in)
        %- plot edits
        fig_i = get(groot,'CurrentFigure');
        tmp = strsplit(conn_conds{cond_i},'_');
        tmp = strsplit(tmp{2},'.');
        title(sprintf('(%s) condition %s',FREQS_SUBPATH,tmp{1}));
        fig_i.Children(2).YTick = 1:length(CLUSTER_ASSIGNMENTS);
        fig_i.Children(2).YTickLabel = CLUSTER_ASSIGNMENTS;
        fig_i.Children(2).YTickLabelRotation = 45;
        fig_i.Children(2).XTick = 1:length(CLUSTER_ASSIGNMENTS);
        fig_i.Children(2).XTickLabel = CLUSTER_ASSIGNMENTS;
        fig_i.Children(2).XTickLabelRotation = 45;
        saveas(fig_i,[save_sub_figs filesep sprintf('3d_boxplot_cond%i.fig',cond_i)]);
        saveas(fig_i,[save_sub_figs filesep sprintf('3d_boxplot_cond%i.jpg',cond_i)]);
    end
end
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