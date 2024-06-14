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
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS SCRIPT_DIR STUDY_DIR%#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = SCRIPT_DIR;
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(SCRIPT_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (PARAMETERS) ======================================================== %%
%## PATHS
%-
ATLAS_PATH = [PATHS.src_dir filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
%- ball machine vs human rally
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
CONN_METHODS = {'dDTF08','S'};
%- datetime override
% study_dir_name = '01232023_MIM_YAN32_antsnormalize_iccREMG0p4_powpow0p3_conn';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%## soft define
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
conn_fig_dir = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'conn_valid_slide'];
%- load study file
STUDY_FNAME_LOAD = 'slide_conn_study';
study_fpath = [studies_fpath filesep sprintf('%s',study_dir_name)];
%- load cluster
CLUSTER_K = 11;
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%- 
save_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K) filesep 'conn_figs'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
if ~ispc
    addpath(convertPath2UNIX('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
else
    addpath(convertPath2Drive('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
end

%% LOAD EPOCH STUDY
%- Create STUDY & ALLEEG structs
%## LOAD STUDY
% if ~ispc
%     [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_LOAD '_UNIX.study'],'filepath',study_fpath);
% else
%     [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME_LOAD '.study'],'filepath',study_fpath);
% end
% cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
% MAIN_STUDY.cluster = cl_struct;
% [comps_out,main_cl_inds,outlier_cl_inds,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);

%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[study_fpath filesep sprintf('%s_UNIX.study',STUDY_FNAME_LOAD)]);
    MAIN_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[study_fpath filesep sprintf('%s.study',STUDY_FNAME_LOAD)]);
    MAIN_STUDY = tmp.STUDY;
end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
MAIN_STUDY.cluster = cl_struct;
[comps_out,~,~,valid_cls] = eeglab_get_cluster_comps(MAIN_STUDY);
%## CUT OUT NON VALID CLUSTERS
inds = setdiff(1:length(comps_out),valid_cls);
comps_out(inds,:) = 0;
clusters = valid_cls;
%-
fPaths = {MAIN_STUDY.datasetinfo.filepath};
fNames = {MAIN_STUDY.datasetinfo.filename};
% condition_gait = unique({MAIN_STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
% subject_chars = {MAIN_STUDY.datasetinfo.subject};
%% MODEL ORDER & VALIDATION DATA
%{
% [tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(save_dir,COND_CHARS);
[tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table([STUDIES_DIR filesep sprintf('%s',study_dir_fname) filesep 'conn_data'],conditions);

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
%% ===================================================================== %%
[MAIN_STUDY,centroid] = std_centroid(MAIN_STUDY,MAIN_ALLEEG,double(string(clusters)),'dipole');
txt_store = cell(length(clusters),1);
atlas_name_store = cell(length(clusters),1);
for k_i = 1:length(clusters)
    k = double(string(clusters(k_i)));
    %## ANATOMY
    dip1 = MAIN_STUDY.cluster(k).all_diplocs;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'single';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    cfg.verbose = 0;
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if ~isempty(sub_indx)
                atlas_name = names{sub_indx};
            end
        else
            atlas_name = names{indx};
        end
    end
    %## ANATOMY
    
    dip1 = MAIN_STUDY.cluster(k).centroid.dipole.posxyz;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name_ct = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'multiple';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    cfg.verbose = 0;
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if ~isempty(sub_indx)
                atlas_name_ct = names{sub_indx};
            end
        else
            atlas_name_ct = names{indx};
        end
    end
    txt_store{k} = [sprintf('CL%i: N=%i\n',k,length(MAIN_STUDY.cluster(k).sets)),...
    sprintf('CL%i: %s\n',k,atlas_name),...
    sprintf('Dip Center: [%0.1f,%0.1f,%0.1f]\n',MAIN_STUDY.cluster(k).dipole.posxyz),...
    sprintf('CENTROID: CL%i: %s\n',k,atlas_name_ct),...
    sprintf('CENTROID: Dip %0.1f,%0.1f,%0.1f]\n\n',MAIN_STUDY.cluster(k).centroid.dipole.posxyz)];
    atlas_name_store{k_i} = sprintf('CL%i: %s\n',k,atlas_name);
end
cellfun(@(x) disp(x),txt_store);
%% ===================================================================== %%
tmp = load([MAIN_STUDY.datasetinfo(1).filepath filesep MAIN_STUDY.datasetinfo(1).filename],'-mat');
FREQ_BOUND = [4,60];
TIME_BOUNDS = [0,3000]; % in seconds
FREQ_INDS = tmp.etc.COND_CAT.Conn.freqs >= FREQ_BOUND(1) & tmp.etc.COND_CAT.Conn.freqs <= FREQ_BOUND(2);
conn_freqs = tmp.etc.COND_CAT.Conn.freqs(FREQ_INDS);
TIME_INDS = tmp.etc.COND_CAT.Conn.winCenterTimes >= TIME_BOUNDS(1) & tmp.etc.COND_CAT.Conn.winCenterTimes <= TIME_BOUNDS(2);
BOOT_ACCU_N = 2000;
%## SEPERATE REST AND GAIT SECTIONS
condition_base = 'rest';
CONN_MEAS = 'dDTF08';
subject_chars = {MAIN_STUDY.datasetinfo.subject};
condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
all_conds = {condition_base,condition_gait{:}};
%-
subj_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
cond_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
comps_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
conn_cond_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
% gait_conn_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
gait_conn_all_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
comps_char_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
win_times_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
group_c = cell(length(MAIN_STUDY.datasetinfo)*length(condition_gait),1);
cnt = 1;
for subj_i = 1:length(MAIN_STUDY.datasetinfo)
    %-
    % EEG_full = MAIN_ALLEEG(subj_i);
    % chk = exist([MAIN_STUDY.datasetinfo(subj_i).filepath filesep sprintf('rest_conn_%s',CONN_MEAS)],'file') &&...
    %     exist([MAIN_STUDY.datasetinfo(subj_i).filepath filesep sprintf('gait_conn_%s',CONN_MEAS)],'file');
    %-
    EEG_full = load([MAIN_STUDY.datasetinfo(subj_i).filepath filesep MAIN_STUDY.datasetinfo(subj_i).filename],'-mat');
    %- get rest indicies
    inds1 = logical(strcmp({EEG_full.event.cond}, condition_base));
    inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
    val_inds = find(inds1 & ~inds2);
    FROM = [EEG_full.event(val_inds(1)).latency];
    TO = [EEG_full.event(val_inds(end)).latency];
    fprintf('%s) Rest length is %0.2fs\n',subject_chars{subj_i},(TO-FROM)/1000);
    %- generate EEG set
    % EEG_BASE = pop_select(EEG_full, 'point', [FROM; TO]');
    % tmp = strsplit(EEG_full.filename,'.');
    % tmp{1} = [tmp{1},'_rest'];
    % tmp = strjoin(tmp,'.');
    % pop_saveset(EEG_BASE,'filepath',EEG_full.filepath,'filename',tmp);
    %- extract connectivity
    winstarts = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
    inds = winstarts > (FROM/1000) & winstarts < (TO/1000);
    tmp_conn = EEG_full.etc.COND_CAT.Conn.(CONN_MEAS);
    rest_conn = tmp_conn(:,:,:,inds);
    % par_save(rest_conn,MAIN_STUDY.datasetinfo(subj_i).filepath,sprintf('rest_conn_%s',CONN_MEAS));
    %## get gait EEG
    EEG_GAIT = cell(length(condition_gait),1);
    TO = zeros(length(condition_gait),1);
    FROM = zeros(length(condition_gait),1);
    tmp_gait_conn = cell(length(condition_gait),1);
    for cond_i = 1:length(condition_gait)
        inds1 = logical(strcmp({EEG_full.event.cond}, condition_gait{cond_i}));
        inds2 = logical(strcmp({EEG_full.event.type}, 'boundary'));
        val_inds = find(inds1 & ~inds2);
        FROM(cond_i) = [EEG_full.event(val_inds(1)).latency];
        TO(cond_i) = [EEG_full.event(val_inds(end)).latency];
        % EEG_GAIT{cond_i} = pop_select(EEG_full, 'point', [FROM; TO]');
        % print
        winstarts = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
        inds = winstarts > (FROM(cond_i)/1000) & winstarts < (TO(cond_i)/1000);
        tmp_conn = EEG_full.etc.COND_CAT.Conn.(CONN_MEAS);
        gait_conn = tmp_conn(:,:,:,inds);
        tmp_gait_conn{cond_i} = gait_conn;
        %-
        fprintf('\n%s) Condition %s''s length is %0.2fs\n',subject_chars{subj_i},...
            condition_gait{cond_i},(TO(cond_i)-FROM(cond_i))/1000);
    end
    % EEG_GAIT =  cellfun(@(x) [[]; x], EEG_GAIT);
    % EEG_GAIT = pop_mergeset(EEG_GAIT,1:length(EEG_GAIT),1);
    % tmp = strsplit(EEG_full.filename,'.');
    % tmp{1} = [tmp{1},'_gait'];
    % tmp = strjoin(tmp,'.');
    % pop_saveset(EEG_GAIT,'filepath',EEG_full.filepath,'filename',tmp);
    FROM = min(FROM);
    TO = max(TO);
    winstarts = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
    inds = winstarts > (FROM/1000) & winstarts < (TO/1000);
    tmp_conn = EEG_full.etc.COND_CAT.Conn.(CONN_MEAS);
    gait_conn = tmp_conn(:,:,:,inds);
    
    % par_save(gait_conn,MAIN_STUDY.datasetinfo(subj_i).filepath,sprintf('gait_conn_%s',CONN_MEAS));
    %## STORE
    cond_c{cnt} = 'rest';
    subj_c{cnt} = MAIN_STUDY.datasetinfo(subj_i).subject;
    win_times_c{cnt} = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
    comps_c{cnt} = EEG_full.etc.COND_CAT.curComps;
    comps_char_c{cnt} = EEG_full.etc.COND_CAT.curComponentNames;
    conn_cond_c{cnt} = rest_conn;
    group_c{cnt} = EEG_full.group;
    cnt = cnt+1;
    %-
    cond_c{cnt} = 'all_gait_median';
    subj_c{cnt} = MAIN_STUDY.datasetinfo(subj_i).subject;
    win_times_c{cnt} = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
    comps_c{cnt} = EEG_full.etc.COND_CAT.curComps;
    comps_char_c{cnt} = EEG_full.etc.COND_CAT.curComponentNames;
    conn_cond_c{cnt} = median(gait_conn,4);
    group_c{cnt} = EEG_full.group;
    cnt = cnt+1;
    %-
    cond_c{cnt} = 'all_gait_mean';
    subj_c{cnt} = MAIN_STUDY.datasetinfo(subj_i).subject;
    win_times_c{cnt} = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
    comps_c{cnt} = EEG_full.etc.COND_CAT.curComps;
    comps_char_c{cnt} = EEG_full.etc.COND_CAT.curComponentNames;
    conn_cond_c{cnt} = mean(gait_conn,4);
    group_c{cnt} = EEG_full.group;
    cnt = cnt+1;
    %-
    for cond_i = 1:length(condition_gait)
        subj_c{cnt} = MAIN_STUDY.datasetinfo(subj_i).subject;
        cond_c{cnt} = condition_gait{cond_i};
        conn_cond_c{cnt} = tmp_gait_conn{cond_i};
        win_times_c{cnt} = EEG_full.etc.COND_CAT.Conn.winCenterTimes;
        comps_c{cnt} = EEG_full.etc.COND_CAT.curComps;
        comps_char_c{cnt} = EEG_full.etc.COND_CAT.curComponentNames;
        group_c{cnt} = EEG_full.group;
        cnt = cnt+1;
    end
end
%-
conn_table = table(subj_c,comps_c,conn_cond_c,...
    cond_c,comps_char_c,win_times_c);
par_save(conn_table,save_dir,'conn_table.mat');
% clear MAIN_ALLEEG
%% ===================================================================== %%
%## TERRAIN vs SPEED vs REST
conds_test = {'rest',{'0p25','0p5','0p75','1p0'},{'flat','low','med','high'}};
cluster_struct = MAIN_STUDY.cluster;
conn_clust_mat = cell(length(cluster_struct),length(cluster_struct),length(conds_test));
subj_clust_mat = cell(length(cluster_struct),length(cluster_struct),length(conds_test));
% subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
for cond_i = 1:length(conds_test)
    subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
    for subj_i = 1:length(subject_chars)
        if length(conds_test{cond_i}) > 1 && iscell(conds_test{cond_i})
            out = cellfun(@(x) strcmp(conn_table.cond_c,x) & strcmp(subject_chars{subj_i},conn_table.subj_c),conds_test{cond_i},'UniformOutput',false);
            subj_ind_t = logical(sum(cat(2,out{:}),2));
        else
            subj_ind_t = strcmp(conn_table.cond_c,conds_test{cond_i}) & strcmp(subject_chars{subj_i},conn_table.subj_c);
            ic_nums = str2double(conn_table.comps_char_c{subj_ind_t});
        end
        tmp = find(subj_ind_t);
        ic_nums = str2double(conn_table.comps_char_c{tmp(1)});
        cl_nums = zeros(length(ic_nums),1);
        for ic_i = 1:length(ic_nums)
            for cc = 2:length(cluster_struct)  % exclude parent cluster
                ind = cluster_struct(cc).sets == subj_i & cluster_struct(cc).comps == ic_nums(ic_i);
                if any(ind)
                    cl_nums(ic_i) = cc;
                end
            end
        end
        %- assign data to cells and save
        subj_cl_ics(cl_nums,cl_nums)=subj_cl_ics(cl_nums,cl_nums)+1;
        for j=1:length(cl_nums)
            for k=1:length(cl_nums)
                fprintf('\tAssiging edge %i->%i\n',cl_nums(j),cl_nums(k))
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=real(squeeze(CAT.Conn.(CONN_MEASURES{conn_i})(j,k,fAllInds,tInds)));
                tmp = conn_table.conn_cond_c(subj_ind_t);
                if length(tmp)>1
                    tmp = {cat(4,tmp{:})};
                end
                subj_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:) = {subject_chars{subj_i}};
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(mean(tmp{subj_ind_t}(j,k,FREQ_INDS,:),4));
                conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(median(tmp{1}(j,k,FREQ_INDS,:),4));
            end
        end
    end
end
%% ===================================================================== %%
%## ALL GAIT vs REST
conds_test = {'all_gait_median','rest'};
cluster_struct = MAIN_STUDY.cluster;
conn_clust_mat = cell(length(cluster_struct),length(cluster_struct),length(conds_test));
% subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
for cond_i = 1:length(conds_test)
    subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
    for subj_i = 1:length(subject_chars)
        subj_ind_t = strcmp(conn_table.cond_c,conds_test{cond_i}) & strcmp(subject_chars{subj_i},conn_table.subj_c);
        ic_nums = str2double(conn_table.comps_char_c{subj_ind_t});
        cl_nums = zeros(length(ic_nums),1);
        for ic_i = 1:length(ic_nums)
            for cc = 2:length(cluster_struct)  % exclude parent cluster
                ind = cluster_struct(cc).sets == subj_i & cluster_struct(cc).comps == ic_nums(ic_i);
                if any(ind)
                    cl_nums(ic_i) = cc;
                end
            end
        end
        %- assign data to cells and save
        subj_cl_ics(cl_nums,cl_nums)=subj_cl_ics(cl_nums,cl_nums)+1;
        for j=1:length(cl_nums)
            for k=1:length(cl_nums)
                fprintf('\tAssiging edge %i->%i\n',cl_nums(j),cl_nums(k))
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=real(squeeze(CAT.Conn.(CONN_MEASURES{conn_i})(j,k,fAllInds,tInds)));
                tmp = conn_table.conn_cond_c(subj_ind_t);
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(mean(tmp{subj_ind_t}(j,k,FREQ_INDS,:),4));
                conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(median(tmp{1}(j,k,FREQ_INDS,:),4));
            end
        end
    end
end
%% ===================================================================== %%
%{
%## EACH GAIT CONDITION SEPERATED
conds_test = all_conds;
cluster_struct = MAIN_STUDY.cluster;
conn_clust_mat = cell(length(cluster_struct),length(cluster_struct),length(all_conds));
% subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
for cond_i = 1:length(all_conds)
    subj_cl_ics = zeros(length(cluster_struct),length(cluster_struct));
    for subj_i = 1:length(subject_chars)
        subj_ind_t = strcmp(conn_table.cond_c,all_conds{cond_i}) & strcmp(subject_chars{subj_i},conn_table.subj_c);
        ic_nums = str2double(conn_table.comps_char_c{subj_ind_t});
        cl_nums = zeros(length(ic_nums),1);
        for ic_i = 1:length(ic_nums)
            for cc = 2:length(cluster_struct)  % exclude parent cluster
                ind = cluster_struct(cc).sets == subj_i & cluster_struct(cc).comps == ic_nums(ic_i);
                if any(ind)
                    cl_nums(ic_i) = cc;
                end
            end
        end
        %- assign data to cells and save
        subj_cl_ics(cl_nums,cl_nums)=subj_cl_ics(cl_nums,cl_nums)+1;
        for j=1:length(cl_nums)
            for k=1:length(cl_nums)
                fprintf('\tAssiging edge %i->%i\n',cl_nums(j),cl_nums(k))
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=real(squeeze(CAT.Conn.(CONN_MEASURES{conn_i})(j,k,fAllInds,tInds)));
                tmp = conn_table.conn_cond_c(subj_ind_t);
                % conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(mean(tmp{subj_ind_t}(j,k,FREQ_INDS,:),4));
                conn_clust_mat{cl_nums(j),cl_nums(k),cond_i}(subj_cl_ics(cl_nums(j),cl_nums(k)),:,:)=squeeze(median(tmp{1}(j,k,FREQ_INDS,:),4));
            end
        end
    end
end
%}
%%
%## (OPTION 1) could potentially sub sample each cluster connection,
%calculate a median then do that a bunch and compare our data to that
%surrogate to determine if certain frequencies and edges are random?
%{
for cond_i = 1:length(all_conds)
    for cl_i = 1:length(cluster_struct)
        fprintf('Performing Stats for Condition %i & Cluster %i\n',cond_i,cl_i);
        tmp = allersp_sb{cond_i,1};
        tmp_mean = mean(tmp,3);
        boot_freq = 1:size(tmp,1);
        boot_subj = 1:size(tmp,3);
        boot_surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
        surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
        %- scramble time samples and calculate the average across
        %all times and all frequencies and store that value.
        for n = 1:BOOT_NITERS
            boot_time = randi(size(tmp,2),[size(tmp,2),1]); % random time samples
            tmpSurro = mean(tmp(boot_freq,boot_time,boot_subj),3);
            surro(:,:,n) = tmpSurro; % save 2000 iterations of surrogates 
        end
        %- Pull length(subject) surrogate averages from distribution then calc mean across
        %surrogates 
        for n = 1:BOOT_NITERS
            bootIdx  = randi(BOOT_NITERS,[size(tmp,3),1]);
            tmpSurro = mean(surro(:,:,bootIdx),3);
            boot_surro(:,:,n) = tmpSurro;
        end
        pvalMap = stat_surrogate_pvals(boot_surro,tmp_mean,'both');
        pvalMap(pvalMap>1)=1; 
        [p_masked, ~, ~, ~] = fdr_bh(pvalMap,BOOT_ALPHA,'pdep',1);
        % debri removal
        [labelMap,~] = bwlabeln(p_masked);
        tmpDisp = sort(labelMap(:),'descend');
        %             [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
        [occurrence,idx,~] = histcounts(tmpDisp,unique(tmpDisp));
        kMask = ismember(labelMap,idx((occurrence<BOOT_CLUST_THRESH)));
        finalMask = p_masked-kMask;
        clust_ersp{cond_i} = tmp_mean; 
        tmp = clust_ersp{cond_i}; 
        tmp(~finalMask) = 0;
        clust_maskedersp{cond_i} = tmp;
    end
end
%}
%%
ttest_pairs = zeros(length(cluster_struct)^2,2);
cnt = 1;
for i = 1:length(cluster_struct)
    for j = 1:length(cluster_struct)
        if i ~= j 
            % chk = intersect([i,j],ttest_pairs,'rows');
            chk1 = ttest_pairs(:,1) == i;
            chk2 = ttest_pairs(:,2) == j;
            chk_to = chk1 & chk2;
            % chk1 = ttest_pairs(:,1) == j;
            % chk2 = ttest_pairs(:,2) == i;
            % chk_from = chk1 & chk2;
            if ~any(chk_to)
                ttest_pairs(cnt,:) = [i,j];
                cnt = cnt + 1;
            end
        end
    end
end
ttest_pairs = ttest_pairs(~all(ttest_pairs == 0,2),:);
%%
% ttest_pairs = [1,2;3,4;5,6];
% cond_ttest = [1,2;1,3;1,4;2,3;2,4;3,4];
% cond_ttest = [1,2];
cond_ttest = [1,2;1,3;2,3];
REPEATED_MEAS = true;
ALPHA = 0.05;
NUM_ITERS = 500;
pval_clust_mat = cell(size(conn_clust_mat,1),size(conn_clust_mat,2),length(cond_ttest));
for i = 1:length(ttest_pairs)
    for c_i = 1:size(cond_ttest,1)
        % in_1 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),1};
        % in_2 = cat(3,conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),2:end});
        in_1 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,1)};
        in_2 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,2)};
        if size(in_1,1) >1 & size(in_2,1) > 1 && ~isempty(in_1) && ~isempty(in_2)
            in_1 = permute(in_1,[2,1]);
            in_1 = reshape(in_1,[1,size(in_1,1),size(in_1,2)]);
            in_2 = permute(in_2,[2,1]);
            in_2 = reshape(in_2,[1,size(in_2,1 ),size(in_2,2)]);
            [mask,tscore,pval] = clusterLevelPermutationTest(in_1,in_2,...
                                REPEATED_MEAS,ALPHA,NUM_ITERS);
            pval_fdr = fdr(pval);
            pval_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),c_i} = pval_fdr;
            fprintf('\nPercent un-masked: %0.1f%%\n',100*(sum(pval_fdr<ALPHA,[1,2])/(size(pval_fdr,1)*size(pval_fdr,2))));
            % curr_ersp = permute(curr_ersp,[2,1,3]);
            % curr_ersp = median(curr_ersp,3);
            % curr_maskedersp = curr_ersp;
            % curr_maskedersp(pval_fdr>=ALPHA) = 0;
        end
    end
end
%% AVERAGE BAND POWER
conds_test = {'rest','speed','terrain'};
freq_bands.alpha = [9,13];
freq_bands.theta = [4,8];
freq_bands.beta = [14,30];
f_fields =fields(freq_bands);
table_sz = size(conn_clust_mat,1)*size(conn_clust_mat,2)*length(conds_test)*length(fields(freq_bands));
% freq_clust_mat = cell();
freq_band = categorical(repmat({''},table_sz,1));
freq_band_vals = cell(table_sz,1);
freq_band_mean = nan(table_sz,1);
cluster_i = nan(table_sz,1);
cluster_j = nan(table_sz,1);
cond_c = categorical(repmat({''},table_sz,1));
subj_c = categorical(repmat({''},table_sz,1));
cnt = 1;
for i = 1:length(ttest_pairs)
    for freq_i = 1:length(f_fields)
        freq_inds = conn_freqs >= freq_bands.(f_fields{freq_i})(1) & conn_freqs <= freq_bands.(f_fields{freq_i})(2);
        for c_i = 1:size(conds_test,2)
            in_1 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),c_i};
            subj_1 = subj_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),c_i};
            if size(in_1,1) >1 && ~isempty(in_1)
                for subj_i = 1:length(subj_1)
                    % freq_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),c_i,freq_i} = mean(median(in_1(:,freq_inds),1));
                    freq_band(cnt) = categorical(f_fields(freq_i));
                    freq_band_vals{cnt} = in_1(subj_i,freq_inds);
                    freq_band_mean(cnt) = mean(in_1(subj_i,freq_inds));
                    % freq_band_subj_med{cnt} = median(in_1(:,freq_inds),1);
                    % freq_band_mean(cnt) = mean(median(in_1(:,freq_inds),1));
                    cluster_i(cnt) = ttest_pairs(i,1);
                    cluster_j(cnt) = ttest_pairs(i,2);
                    cond_c(cnt) = categorical(conds_test(c_i));
                    subj_c(cnt) = categorical(subj_1(subj_i));
                    cnt = cnt + 1;
                end
            end
        end
    end
end
%##
freq_band_table = table(freq_band,freq_band_vals,freq_band_mean,cluster_i,cluster_j,cond_c,subj_c);
freq_band_table = rmmissing(freq_band_table);
%% PLOT VIOLINS & STATS ================================================ %%
%## LEFT OFF HERE, work on making nice summary figure for these stats. Good
%direction to go with the analysis?
% freq_band_table = freq_band_table(all(isempty(freq_band_table{:}),2))
clusters_i = unique(freq_band_table.cluster_i);
clusters_j = unique(freq_band_table.cluster_j);
freqs = unique(freq_band_table.freq_band);
for i = 1:length(clusters_i)
    for j = 1:length(clusters_j)
        inds = freq_band_table.cluster_i == clusters_i(i) & freq_band_table.cluster_j == clusters_j(j);
        table_in = freq_band_table(inds,:);
        if ~isempty(table_in)
            %- test lienar model
            STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{[]},...
                          'regress_xvals',0);
            for freq_i = 1:length(freqs)
                inds = table_in.freq_band == freqs(freq_i);
                sub_table_in = table_in(inds,:);
                %-
                % mod = 'freq_band_mean ~ 1 + cond_c + freq_band + cond_c:freq_band';
                mod = 'freq_band_mean ~ 1 + cond_c';
                lme_out = fitlme(sub_table_in,mod,'FitMethod','ML');
                anova_out = anova(lme_out);
                %- test normality
                [theta_h,theta_p] = lillietest(lme_out.residuals);
                %-
                STATS_STRUCT.anova{freq_i} = anova_out.pValue(2);
                STATS_STRUCT.pvals{freq_i} = [1,lme_out.Coefficients.pValue(2),lme_out.Coefficients.pValue(3)];
                STATS_STRUCT.pvals_pairs{freq_i} = {[1,1],[1,2],[1,3]};
            end
            %-
            %## PLOT
            ax = group_violin(table_in,'freq_band_mean','cond_c','freq_band',...
                'STATS_STRUCT',STATS_STRUCT);
            % drawnow;
            % exportgraphics(ax,[save_dir filesep sprintf('cl%i-cl%i_cond%i-cond%i_spedterrainvsrest_conn_clusterttests.tiff',ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,1),cond_ttest(c_i,2))],'Resolution',300)
            % close(ax);
        end
    end
end
%% MULTI-CLUSTER PLOT OF ALL SUBJECTS ================================== %%
%-
color_dark = linspecer(length(cond_ttest));
color_light = color_dark*0.5;
IM_RESIZE = 0.5;
MAX_VT_DIM = 4;
VERTICAL_SHIFT =  0.2;
HORIZONTAL_SHIFT = 0.3;
HORIZ_START = 0.08;
VERTICAL_START = 0.75;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%##
vert_shift = 0;
vt = 1;
for i = 1:length(ttest_pairs)
    for c_i = 1:size(cond_ttest,1)
        horiz_shift = HORIZ_START;
        %-
        in_1 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,1)};
        in_2 = conn_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,2)};
        pvals = pval_clust_mat{ttest_pairs(i,1),ttest_pairs(i,2),c_i};
        clim = [min([in_1;in_2],[],'all')-std([in_1;in_2],[],'all'),max([in_1;in_2],[],'all')+std([in_1;in_2],[],'all')];
        if size(in_1,1) >1 & size(in_2,1) > 1 && ~isempty(in_1) && ~isempty(in_2)
            if vt > MAX_VT_DIM || c_i == 1
                fig = figure('color','white','renderer','Painters');
                sgtitle(sprintf('ttest: %i <=> %i',ttest_pairs(i,1),ttest_pairs(i,2)),'FontName','Arial','FontSize',14,'FontWeight','bold','Interpreter','none');
                set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
                set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
                hold on;
                set(gca,AXES_DEFAULT_PROPS{:})
                %-
                vert_shift = 0;
                vt = 1;
            end
            %-
            axes();
            hold on;
            %## IN_1
            in_1_mean = mean(in_1);
            in_1_median = median(in_1);
            subjs = plot(conn_freqs,in_1,'color',[0,0,0,0.15],'linestyle','-','linewidth',2,'displayname','subject conn');
            mean_plot = plot(conn_freqs,in_1_mean,'color',color_dark(cond_ttest(c_i,1),:),'linestyle','-','linewidth',4,'displayname','mean conn');
            median_plot = plot(conn_freqs,in_1_median,'color',color_dark(cond_ttest(c_i,1),:),'linestyle','--','linewidth',4,'displayname','median conn');
            %-
            ax = gca;
            plot([0 40],[0 0],'--','color','black');
            xlim([3 60]);
            ylim([clim(1) clim(2)]);
            xlabel('Frequency(Hz)');
            ylabel('connectivity');
            set(ax,'FontName','Arial',...
                'FontSize',12,...
                'FontWeight','bold');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName','Arial','FontSize',10,...
                'FontWeight','bold');
            % title(sprintf('Cond %s',conds_test{cond_ttest(c_i,1)}))
            title(sprintf('Cond %i',cond_ttest(c_i,1)))
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[horiz_shift,VERTICAL_START-vert_shift,0.4*IM_RESIZE,0.25*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + HORIZONTAL_SHIFT;

            %## IN_2
            axes();
            hold on;
            in_2_mean = mean(in_2);
            in_2_median = median(in_2);
            subjs = plot(conn_freqs,in_2,'color',[0,0,0,0.15],'linestyle','-','linewidth',2,'displayname','subject conn');
            mean_plot = plot(conn_freqs,in_2_mean,'color',color_dark(cond_ttest(c_i,2),:),'linestyle','-','linewidth',4,'displayname','mean conn');
            median_plot = plot(conn_freqs,in_2_median,'color',color_dark(cond_ttest(c_i,2),:),'linestyle','-','linewidth',4,'displayname','median conn');
            %-
            ax = gca;
            plot([0 40],[0 0],'--','color','black');
            xlim([3 60]);
            ylim([clim(1) clim(2)]);
            xlabel('Frequency(Hz)');
            ylabel('');
            set(ax,'FontName','Arial',...
                'FontSize',12,...
                'FontWeight','bold');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName','Arial','FontSize',10,...
                'FontWeight','bold');
            % title(sprintf('Cond %s',conds_test{cond_ttest(c_i,2)}))
            title(sprintf('Cond %i',cond_ttest(c_i,2)))
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[horiz_shift,VERTICAL_START-vert_shift,0.4*IM_RESIZE,0.25*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + HORIZONTAL_SHIFT;

            %## PVALS
            axes();
            hold on;
            subjs = plot(conn_freqs,pvals<0.05,'color',[0,0,0],'linestyle','-','linewidth',2,'displayname','subject conn');
            %-
            ax = gca;
            plot([0 40],[0 0],'--','color','black');
            xlim([3 60]);
            ylim([0 1.5]);
            xlabel('Frequency(Hz)');
            ylabel('');
            set(ax,'FontName','Arial',...
                'FontSize',12,...
                'FontWeight','bold');
            xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
            set(ax,'FontName','Arial','FontSize',10,...
                'FontWeight','bold');
            title(sprintf('Cluster Ttest'));
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[horiz_shift,VERTICAL_START-vert_shift,0.4*IM_RESIZE,0.25*IM_RESIZE]);  %[left bottom width height]
            horiz_shift = horiz_shift + HORIZONTAL_SHIFT;
            vert_shift = vert_shift + VERTICAL_SHIFT;
            vt = vt + 1;
            if vt > MAX_VT_DIM || c_i == size(cond_ttest,1)
                hold off;
                exportgraphics(fig,[save_dir filesep sprintf('cl%i-cl%i_cond%i-cond%i_spedterrainvsrest_conn_clusterttests.tiff',ttest_pairs(i,1),ttest_pairs(i,2),cond_ttest(c_i,1),cond_ttest(c_i,2))],'Resolution',300)
                close(fig);
            end
        end
    end
end


%%
%%
% fpath = [destination_folder filesep CONN_MEASURES{1} filesep 'R_data' filesep COND_NAMES{4}];
% base_conn = par_load(fpath,sprintf('connStruct%s_boot.mat',CONN_MEASURES{1}));
for cond_i = 1:length(conditions)
    fpath = [conn_fig_dir filesep CONN_MEASURES{conn_i} filesep 'R_data' filesep conditions{cond_i}];
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
conn_vec = {cell(length(conditions),length(cluster_ints),length(cluster_ints))};
subj_vec = {cell(length(conditions),length(cluster_ints),length(cluster_ints))};
cond_vec = {cell(length(conditions),length(cluster_ints),length(cluster_ints))};
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
    fpath = [conn_fig_dir filesep CONN_MEASURES{conn_i} filesep 'R_data'];
    if ~exist(fpath,'dir')
        mkdir(fpath);
    end
    for cond_i=1:length(conditions)
        net_vals = NET_VALS;
        net_vals_ave = NETVALS_AVE;
        connStruct_boot = par_load(fpath,sprintf('connStruct%s_basecorr.mat',CONN_MEASURES{conn_i}));
%         connStruct_boot = par_load([fpath filesep COND_NAMES{cond_i}],sprintf('connStruct%s_boot.mat',CONN_MEASURES{conn_i}));
        connStruct = par_load([fpath filesep conditions{cond_i}],sprintf('connStruct%s_baseSub.mat',CONN_MEASURES{conn_i}));
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
        par_save(net_vals,fpath,sprintf('netVals_%s_aveTime_sbjs.mat',strjoin(strsplit(conditions{cond_i},' '),'_')));
        par_save(net_vals_ave,fpath,sprintf('netVals_%s_aveTime.mat',strjoin(strsplit(conditions{cond_i},' '),'_')));
    end
end

%% corrected conn validation
sub_ints = [7,1,5,6,9,4];
% sub_ints = [1,2,3,4,5,6,7,8,9];
cluster_names = {'RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','SuppMotor','LFrontal'}; % (06/27/2023) JS, unsure on these as of yet.
COLOR_LIMITS = [-0.00001,0.00001];
for cond_i = 1:length(conditions)
    datain = par_load(fpath,sprintf('netVals_%s_aveTime.mat',strjoin(strsplit(conditions{cond_i},' '),'_')));
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
        exportgraphics(fig_i,[conn_fig_dir filesep sprintf('%s_condheatmap_%s.jpg',conditions{cond_i},fn{i})],'Resolution',300);
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
        exportgraphics(fig_i,[conn_fig_dir filesep sprintf('%s_condheatmap_%s.jpg',conditions{cond_i},fn{data_i})],'Resolution',300);
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