%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: s

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
COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
EVENT_CHARS = {'Subject_hit'}; %, 'Subject_receive'};
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
dt = '06152023_bounces_1h2h2bm_JS';
%- connectiviy specific
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % AS (06/22/2023)
conn_meas = 'dDTF08';
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
% CLUSTER_ITERS = [3,4,5,6,7,8,9,10,11,12];
CLUSTER_ITERS = [2,3,4,5,6,7,8,9,10,11];
CLUSTER_ASSIGNMENTS = {'RPPa','Cuneus','Precuneus','RSuppMotor','LPPa','LSM','RSM','LTemp','Cing','LSuppMotor'}; % (06/27/2023) JS, unsure on these as of yet.
% SUBJS_TEST = 1:length(ALLEEG);
CONN_MEAS_ANLYZ = 'dDTF';
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
mat_out_nan = nan(length(main_cl_inds),length(main_cl_inds),length(ALLEEG),length(conn_conds),size(uniq_freqs,1));
%## LOOP
for subj_i = 1:length(ALLEEG)
    tmp = ALLEEG(subj_i).etc.conn_table;
%     conn_meas = unique(tmp.t_conn_meas);
    conn_conds = unique(tmp.t_fNames);
    conn_comps = [tmp.t_conn_comps{1}];
    idx = find(strcmp(CONN_MEAS_ANLYZ,tmp.t_conn_meas));
    mats = tmp.t_conn_mats(idx);
%     sigs = tmp.t_conn_sigs;
    cnt = 1;
    mat_nan = nan(length(main_cl_inds),length(main_cl_inds));
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
%{
% hlpr = fields(ALLEEG(subj_i).etc.js_processing(cond_i).META);
for freq_i = FREQ_INTS
    fprintf('==== Processing Frequency Iter %i ====',freq_i);
    %- loop through each condition   
    for cond_i = 1:length(load_trials)
        %## LOOP PARAMS
        sub_cond_i = cond_i + find(strcmp(load_trials{1},TRIAL_TYPES))-1; % offset for weirdness
        clusterNames = cell(1,length(CLUSTER_ASSIGNMENTS));
        name_trial = load_trials{cond_i};
        mat_nan = nan(length(STUDY.cluster),length(STUDY.cluster),length(ALLEEG));
        ALLEEG = ALLEEGS{cond_i};
        %## LOOP PATHING
        
        %## LOOP MEAT
        %- loop through subjects in each condition
        for subj_i = 1:length(ALLEEG)
            chk = (comps_out(subj_i,1:end)>0);
            if any(chk)
                fprintf('%s) Frequencies: %sHz\n',ALLEEG(subj_i).subject,...
                    FREQS_SUBPATH);
                statMat = ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).mats;
                comps = ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.connComponents';
                meanMat = squeeze(statMat(1,:,:));
    %             stdvMat = squeeze(statMat(2,:,:));
    %             medMat = squeeze(statMat(3,:,:));
                clust_idx = zeros(1,length(comps));
                comp_idx = zeros(1,length(comps));
                fprintf('%s) Number of absent connections: %i\n',ALLEEG(subj_i).subject,...
                    sum(isnan(meanMat(:))));
                for i = 1:length(comps)
                    chk = find(comps(i) == comps_out(subj_i,1:end));
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
                    mat_nan(clust_idx,clust_idx,subj_i) = val_in;
                else
                    continue;
                end                    
            else
                continue;
            end
        end
        %- store
        mat_out_nan(:,:,:,cond_i,freq_i) = mat_nan;
    end
end
%}
%% CONNECITIVITY CHORD PLOTS
BASELINE_INT = 1;
for freq_i = size(uniq_freqs,1)
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
        for i = main_cl_inds
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
%         tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
        I = eye(size(tmp));
        I = (I == 0);
        tmp = tmp.*I;
        tmp(tmp == 0) = nan();
        %- delte unused clusters
        tmp = tmp(CLUSTER_ITERS,:);
        tmp = tmp(:,CLUSTER_ITERS);
        group_1 = zeros(length(tmp(:)),1);
        group_2 = zeros(length(tmp(:)),1);
        cnt = 1;
        for i = 1:size(tmp,1)
            for j = 1:size(tmp,2)
                group_1(cnt) = j;
                group_2(cnt) = i;
                cnt = cnt + 1;
            end
        end
        data_in = tmp(:);
        sign_in = sign(tmp(:));
%         data_in(isnan(data_in)) = 0;
        data_in = abs(data_in);
        sign_in(isnan(sign_in)) = 1;
        chord_plot_in = table(group_1,group_2,data_in,sign_in);
        h = chordPlot(clusterNames,chord_plot_in);
       
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.fig','clusters',conn_conds{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.jpg','clusters',conn_conds{cond_i})]);
        close(fig_i)
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
%         tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
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
        hnd = heatmap(tmp,'Colormap',jet); %,'CellLabelColor', 'None');
        title(sprintf('''%s'' Connectivity mean across clusters',conn_conds{cond_i}));
        hnd.YDisplayLabels = clusterNames;
        hnd.XDisplayLabels = clusterNames;
        hnd.ColorLimits = [0,0.05];
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.fig','clusters',conn_conds{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.jpg','clusters',conn_conds{cond_i})]);
%         close(fig_i)
    end
end
%% Condition Comparisons
%- if using comp selection
for freq_i = 1:size(uniq_freqs,1)
    store_val_raw = zeros(2,size(mat_out_nan,3),length(conn_conds));
    done = [];
    %## PATHING
    tmp = uniq_freqs(freq_i,uniq_freqs(freq_i,:)>0);
        FREQS_SUBPATH = sprintf('%i-%i',...
                    tmp(1),tmp(end));
    save_sub_figs = [save_dir filesep 'avg_across_conditions' filesep FREQS_SUBPATH];
    if ~exist(save_sub_figs,'dir')
        mkdir(save_sub_figs);
    end
    %## LOOP MEAT
    baseline =  nanmean(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    for i = 1:length(CLUSTER_ITERS)
        clust_i = CLUSTER_ITERS(i);
        for j = 1:length(CLUSTER_ITERS)
            clust_j = CLUSTER_ITERS(j);
            if clust_i == clust_j || any((j == done))
                continue;
            end
            for cond_i = 1:length(conn_conds)
%                 store_val_raw(1,:,cond_i) = squeeze(mat_out_nan(clust_i,clust_j,:,cond_i,freq_i))-baseline(clust_i,clust_j);
%                 store_val_raw(2,:,cond_i) = squeeze(mat_out_nan(clust_j,clust_i,:,cond_i,freq_i))-baseline(clust_i,clust_j);
                store_val_raw(1,:,cond_i) = squeeze(mat_out_nan(clust_i,clust_j,:,cond_i,freq_i));
                store_val_raw(2,:,cond_i) = squeeze(mat_out_nan(clust_j,clust_i,:,cond_i,freq_i));
            end
            two_clusts = {CLUSTER_ASSIGNMENTS{i},CLUSTER_ASSIGNMENTS{j}};
            subj_clusts_comps = [comps_out(clust_i,:)',comps_out(clust_j,:)'];
            %- plot function
            cnctanl_plot_cond(store_val_raw,conn_conds,...
                        two_clusts,subj_clusts_comps,save_sub_figs);
            close all
        end
        done = [done i];
    end
end
%% VISUALIZE TIME FREQ GRID (components to component)
%## TIMEFREQ AVG PLOT
for cond_i = 1:length(conn_conds)
    %- loop through subjects in each condition
    for subj_i = 1:length(ALLEEG)
        %## LOAD SUBJECT CAT.Stats
%         EEG = cnctanl_loadCAT(ALLEEG(subj_i),'NonzeroTest');
        EEG = ALLEEG(subj_i);
        EEG.CAT = EEG.etc.COND_CAT(cond_i);
        %- Phase Randomization Permutation Test Data Handler
        fprintf('\n==== LOADING PHASE RANDOMIZED CONNECTIVITY MEASURES ====\n')
        if ispc
            fPath = convertPath2Drive(EEG.etc.cond_files(cond_i).fPath);
        else
            fPath = convertPath2UNIX(EEG.etc.cond_files(cond_i).fPath);
        end
        fName = EEG.etc.cond_files(cond_i).fName;
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file') 
            error('%s does not exist.\nRun GLOBAL_BATCH to generate phase randomized permutation test values',[fPath filesep [chk{1}, '_PhaseRnd.mat']]);
        else
            EEG.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
        end
        fprintf('done.\n')
        %- Nonzero Statistics Data Handler
        fprintf('\n==== LOADING NONZERO STATISTICS ====\n')
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file')
            error('%s does not exist.\nRun GLOBAL_BATCH to generate nonzero test values',[fPath filesep [chk{1}, '_NonZero.mat']]);
        else
            EEG.CAT.Stats = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
        end
        fprintf('done.\n')
        %## LOOP PATHS
        figs_save_dir = [save_dir filesep 'indvidual_TxF_conn' filesep conn_conds{cond_i} filesep sprintf('%s',EEG.subject)];
        if ~exist(figs_save_dir,'dir')
            mkdir(figs_save_dir)
        end
        %## LOOP MEAT
        %- generate display names based on CLUSTER_ASSIGNMENTS
        orig_cl = squeeze(comps_out(:,subj_i));
        orig_cl = orig_cl(2:end);
        [tmpcl,idxcl] = sort(orig_cl);
        display_names = cell(1,length(tmpcl));
%         for i = idxcl
%             if i == 1
%                 continue;
%             else
%                 display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{i-1},orig_cl(i));
%             end    
%         end
        for i = 1:length(idxcl)
            display_names{idxcl(i)} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{idxcl(i)},orig_cl(idxcl(i)));
        end
        idxcl = idxcl(tmpcl~=0);
        idx = ~cellfun(@isempty, display_names);
        display_names = display_names(idx);
        display_names = display_names(idxcl);
        %- assign component indicies from EEG structure
        comp_indxs = EEG.CAT.curComps;
        %- plot
        cnctanl_vis_timefreq(EEG,comp_indxs,...
            'CONN_MEASURES',{conn_meas},...
            'SAVE_DIR',figs_save_dir,...
            'DISPLAY_NAMES',display_names,...
            'COLOR_LIM',[0,0.005])
        close all
    end
end
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
%% Helper Code
%{
delete STUDY ALLEEG
cond_i = 1;
STUDY = STUDIES{cond_i};
ALLEEG = ALLEEGS{cond_i};
%-
path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
mniMRI = fullfile(path2BEM, 'standard_mri.mat');
mniVol = fullfile(path2BEM, 'standard_vol.mat');

for subj_i = 1:length(ALLEEG)
     ALLEEG(subj_i).dipfit.mrifile = mniMRI;
     ALLEEG(subj_i).dipfit.hdmfile = mniVol;
     ALLEEG(subj_i).dipfit.coordformat = 'MNI';
end

%% CONDITION STATISTICS FOR ICA EPOCHS
ALLEEG = eeg_checkset(ALLEEG,'loaddata');
for subj_i = 1:length(ALLEEG)
    if isempty(ALLEEG(subj_i).icaact)
        fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
        ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
    end
end
%% REORGANIZE DIPFITS
fprintf('==== Reorganizing dipfit indices ====\n');
for subj_i = 1:length(ALLEEG)
    dipfit_i = ALLEEGS{cond_i}(subj_i).dipfit.model;
    vals = [];
    %- extract cluster number and associated component number
    for clust_i = 2:length(STUDIES{cond_i}.cluster)
        idx = (STUDIES{cond_i}.cluster(clust_i).sets == subj_i);
        if any(idx)
            vals = [vals; clust_i, find(idx), STUDIES{cond_i}.cluster(clust_i).comps(idx)];
        end
    end
    %- 
    [~,idx] = sort(vals);
    for val_i = 1:size(vals,1)
        STUDIES{cond_i}.cluster(vals(idx(val_i,3),1)).comps(vals(idx(val_i,3),2)) = val_i;
    end
    STUDIES{cond_i}.datasetinfo(subj_i).comps = idx(:,3);
end
all_sets = [];
all_comps = [];
for clust_i = 2:length(STUDIES{cond_i}.cluster)
    all_sets = [all_sets, STUDIES{cond_i}.cluster(clust_i).sets];
    all_comps = [all_sets, STUDIES{cond_i}.cluster(clust_i).comps];
end
STUDIES{cond_i}.cluster(1).comps = all_comps;
STUDIES{cond_i}.cluster(1).sets = all_sets;

%% COMPARISON PLOTS
std_dipplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster),'mode','multicolor');
% view([45,0,0])
% view([0,-45,0])
% view([0,0,45])
%- Spec plot
specMin = 15;
specMax = 40;
std_specplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Topo plot
std_topoplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster));
%}