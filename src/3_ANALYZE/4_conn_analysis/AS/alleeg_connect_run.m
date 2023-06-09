%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_test/3_paper_MIM_HOA/run_alleeg_connect_run.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    DO_UNIX = false;
    PATH_EXT = 'M';
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '_test' filesep '3_paper_MIM_HOA'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
setWorkspace
%% PARPOOL SETUP
if DO_UNIX
%     eeg_options;
    pop_editoptions('option_parallel',1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number o f workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i',pp.NumWorkers);
    fprintf('Number of threads: %i',pp.NumThreads);
    %- make meta data directory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', 1440);
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
%% ===================================================================== %%
%## PROCESSING PARAMS
%- subjstruct and saving
LOAD_STUDY = true;
%## POST-PROCESSING PARAMS
%- datetime override
dt = '16022023';
%- hard define
% CONN_METHODS = {'dDTF','GGC','ffDTF'}; % MIM_YA (12/26/2022)
CONN_METHODS = {'dDTF','GGC','dDTF08'}; % MIM_YA (12/26/2022)
conn_meas = 'dDTF08';
% load_trials = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
% load_trials = {'rest','0p25','0p5','0p75','1p0'}; 
load_trials = {'flat','low','med','high'};
% load_trials = {'flat'};
%- soft define
subjinfDir = [SUBJINF_DIR filesep sprintf('%s',dt)];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%- create new subject directory
if ~exist(subjinfDir,'dir')
    mkdir(subjinfDir);
end
%% ==== Load Study ==== %%
STUDIES = cell(1,length(load_trials));
ALLEEGS = cell(1,length(load_trials));
%- Create STUDY & ALLEEG structs
for cond_i = 1:length(load_trials)
    study_fName = sprintf('%s_MIM_study',load_trials{cond_i});
    fprintf('==== Loading STUDY %s ====\n',study_fName);
    if ~exist([load_dir filesep study_fName '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if DO_UNIX
            [STUDIES{cond_i},ALLEEGS{cond_i}] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',load_dir);
        else
            [STUDIES{cond_i},ALLEEGS{cond_i}] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);
        end
    end
    %## NOTE: DO NOT REORGANIZE DIPFIT INDICIES
end
%## Extract component store
MAIN_STUDY = STUDIES{1};
MAIN_ALLEEG = ALLEEGS{1};
%- extract component array
comps_store = zeros(length(MAIN_ALLEEG),length(MAIN_STUDY.cluster));
for clus_i = 2:length(MAIN_STUDY.cluster)
    sets_i = MAIN_STUDY.cluster(clus_i).sets;
    for j = 1:length(sets_i)
        comps_store(sets_i(j),clus_i) = MAIN_STUDY.cluster(clus_i).comps(j);
    end
end
%%
%## CONN MAT PARAMS
% FREQ_INTS = [1,2,3,4,5];
FREQ_INTS = (1:length(MAIN_ALLEEG(1).etc.js_processing(1).META.(conn_meas).connExtract));
% CLUSTER_ITERS = [2,3,4,5,6,7,8,9];
% CLUSTER_ASSIGNMENTS = {'RPPa','LPPa','PFC','LSM','ROc','ACC','RSM','Cerr'};
CLUSTER_ITERS = [2,3,4,5,6,7,8,9];
CLUSTER_ASSIGNMENTS = {'Cune','RSM','Cerr','ACC','RPPa','LPPa','LSM','LOc'};
mat_out_nan = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),length(MAIN_ALLEEG),length(load_trials),length(FREQ_INTS));
SUBJS_TEST = 1:length(MAIN_ALLEEG);
%## LOOP
% hlpr = fields(MAIN_ALLEEG(subj_i).etc.js_processing(cond_i).META);
for freq_i = FREQ_INTS
    fprintf('==== Processing Frequency Iter %i ====',freq_i);
    %- loop through each condition   
    for cond_i = 1:length(ALLEEG(1).CAT)
        %## LOOP PARAMS
        sub_cond_i = cond_i + find(strcmp(load_trials{1},TRIAL_TYPES))-1; % offset for weirdness
        clusterNames = cell(1,length(CLUSTER_ASSIGNMENTS));
        name_trial = load_trials{cond_i};
        mat_nan = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),length(MAIN_ALLEEG));
        MAIN_ALLEEG = ALLEEGS{cond_i};
        %## LOOP PATHING
        
        %## LOOP MEAT
        %- loop through subjects in each condition
        for subj_i = 1:length(MAIN_ALLEEG)
            chk = (comps_store(subj_i,1:end)>0);
            if any(chk)
                fprintf('%s) Frequencies: %sHz\n',MAIN_ALLEEG(subj_i).subject,...
                    FREQS_SUBPATH);
                statMat =  MAIN_ALLEEG(subj_i).CAT(i).Conn.dDTF; MAIN_ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).mats;
                comps = MAIN_ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.connComponents';
                meanMat = squeeze(statMat(1,:,:));
    %             stdvMat = squeeze(statMat(2,:,:)); 
    %             medMat = squeeze(statMat(3,:,:));
                clust_idx = zeros(1,length(comps));
                comp_idx = zeros(1,length(comps));
                fprintf('%s) Number of absent connections: %i\n',MAIN_ALLEEG(subj_i).subject,...
                    sum(isnan(meanMat(:))));
                for i = 1:length(comps)
                    chk = find(comps(i) == comps_store(subj_i,1:end));
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
%% CONNECITIVITY CHORD PLOTS
BASELINE_INT = 1;
for freq_i = FREQ_INTS
    baseline = nanmean(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    for cond_i = 1:length(load_trials)
        sub_cond_i = cond_i + find(strcmp(load_trials{1},TRIAL_TYPES))-1; % offset for weirdness
        
        FREQS_SUBPATH = sprintf('%i-%i',...
                    MAIN_ALLEEG(1).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(1),...
                    MAIN_ALLEEG(1).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(end));
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cnt = 1;
        for i = 1:length(MAIN_STUDY.cluster)
            N = length(MAIN_STUDY.cluster(i).sets);
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
        tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
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
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.fig','clusters',load_trials{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.jpg','clusters',load_trials{cond_i})]);
        close(fig_i)
    end
end
%% CONNECTIVITY MATRICIES
BASELINE_INT = 1;
for freq_i = FREQ_INTS
    baseline = nanmean(mat_out_nan(:,:,:,BASELINE_INT,freq_i),3);
    for cond_i = 1:length(load_trials)
        sub_cond_i = cond_i + find(strcmp(load_trials{1},TRIAL_TYPES))-1; % offset for weirdness
        
        FREQS_SUBPATH = sprintf('%i-%i',...
                    MAIN_ALLEEG(1).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(1),...
                    MAIN_ALLEEG(1).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(end));
        save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
        if ~exist(save_sub_figs,'dir')
            mkdir(save_sub_figs);
        end
        %- plot
        cnt = 1;
        for i = 1:length(MAIN_STUDY.cluster)
            N = length(MAIN_STUDY.cluster(i).sets);
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
        tmp = nanmean(mat_out_nan(:,:,:,cond_i,freq_i),3) - baseline;
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
        title(sprintf('''%s'' Connectivity mean across clusters',load_trials{cond_i}));
        hnd.YDisplayLabels = clusterNames;
        hnd.XDisplayLabels = clusterNames;
        hnd.ColorLimits = [0,0.05];
        hnd.GridVisible = 'off';
        hnd.CellLabelFormat = '%0.1g';
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.fig','clusters',load_trials{cond_i})]);
        saveas(fig_i,[save_sub_figs filesep sprintf('%s_meanMat_%s.jpg','clusters',load_trials{cond_i})]);
        close(fig_i)
    end
end
%% Condition Comparisons
%- if using comp selection
for freq_i = FREQ_INTS
    store_val_raw = zeros(2,size(mat_out_nan,3),length(load_trials));
    done = [];
    %## PATHING
    FREQS_SUBPATH = sprintf('%i-%i',...
                MAIN_ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(1),...
                MAIN_ALLEEG(subj_i).etc.js_processing(sub_cond_i).META.(conn_meas).connExtract(freq_i).freqs(end));
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
            for cond_i = 1:length(load_trials)
                store_val_raw(1,:,cond_i) = squeeze(mat_out_nan(clust_i,clust_j,:,cond_i,freq_i))-baseline(clust_i,clust_j);
                store_val_raw(2,:,cond_i) = squeeze(mat_out_nan(clust_j,clust_i,:,cond_i,freq_i))-baseline(clust_i,clust_j);
            end
            two_clusts = {CLUSTER_ASSIGNMENTS{i},CLUSTER_ASSIGNMENTS{j}};
            subj_clusts_comps = [comps_store(:,clust_i),comps_store(:,clust_j)];
            %- plot function
            cnctanl_plot_cond(store_val_raw,load_trials,...
                        two_clusts,subj_clusts_comps,save_sub_figs);
            close all
        end
        done = [done i];
    end
end
%% VISUALIZE TIME FREQ GRID (components to component)
%## TIMEFREQ AVG PLOT
for cond_i = 1:length(load_trials)
    sub_cond_i = cond_i;
%     sub_cond_i = cond_i + 5; % offset for weirdness
    MAIN_ALLEEG = ALLEEGS{sub_cond_i};
    %- loop through subjects in each condition
    for subj_i = 1:length(MAIN_ALLEEG)
        %## LOAD SUBJECT CAT.Stats
        EEG = cnctanl_loadCAT(MAIN_ALLEEG(subj_i),'NonzeroTest');
        %## LOOP PATHS
        figs_save_dir = [save_dir filesep 'indvidual_TxF_conn' filesep load_trials{sub_cond_i} filesep sprintf('%s',EEG.subject)];
        if ~exist(figs_save_dir,'dir')
            mkdir(figs_save_dir)
        end
        %## LOOP MEAT
        %- generate display names based on CLUSTER_ASSIGNMENTS
        orig_cl = squeeze(comps_store(subj_i,:));
        [tmpcl,idxcl] = sort(orig_cl);
        display_names = cell(1,length(tmpcl));
        for i = idxcl
            if i == 1
                continue;
            else
                display_names{i} = sprintf('%s_ic%i',CLUSTER_ASSIGNMENTS{i-1},orig_cl(i));
            end    
        end
        idxcl = idxcl(tmpcl~=0);
        idx = ~cellfun(@isempty, display_names);
        display_names = display_names(idx);
        display_names = display_names(idxcl-1);
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
        figs_save_dir = [save_dir filesep 'anovas' filesep [load_trials{:}] filesep sprintf('%s_%s',CLUSTER_ASSIGNMENTS{i},CLUSTER_ASSIGNMENTS{j})];
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
        violin_data = nan(size(mat_out_nan,3),length(load_trials));
        violin_group = cell(1,length(load_trials));
        for cond_i = 1:length(load_trials)
        %     violin_data{cond_i} = squeeze(mat_out_nan(CLUST_I,CLUST_J,:,cond_i,freq_i));
            violin_data(:,cond_i) = squeeze(mat_out_nan(clust_i,clust_j,:,cond_i,freq_i))';
            violin_group{cond_i} = load_trials{cond_i};
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
delete MAIN_STUDY MAIN_ALLEEG
cond_i = 1;
MAIN_STUDY = STUDIES{cond_i};
MAIN_ALLEEG = ALLEEGS{cond_i};
%-
path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
mniMRI = fullfile(path2BEM, 'standard_mri.mat');
mniVol = fullfile(path2BEM, 'standard_vol.mat');

for subj_i = 1:length(MAIN_ALLEEG)
     MAIN_ALLEEG(subj_i).dipfit.mrifile = mniMRI;
     MAIN_ALLEEG(subj_i).dipfit.hdmfile = mniVol;
     MAIN_ALLEEG(subj_i).dipfit.coordformat = 'MNI';
end

%% CONDITION STATISTICS FOR ICA EPOCHS
MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
for subj_i = 1:length(MAIN_ALLEEG)
    if isempty(MAIN_ALLEEG(subj_i).icaact)
        fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
        MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
    end
end
%% REORGANIZE DIPFITS
fprintf('==== Reorganizing dipfit indices ====\n');
for subj_i = 1:length(MAIN_ALLEEG)
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
std_dipplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),'mode','multicolor');
% view([45,0,0])
% view([0,-45,0])
% view([0,0,45])
%- Spec plot
specMin = 15;
specMax = 40;
std_specplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Topo plot
std_topoplot(MAIN_STUDY,MAIN_ALLEEG,'clusters',2:length(MAIN_STUDY.cluster));
%}