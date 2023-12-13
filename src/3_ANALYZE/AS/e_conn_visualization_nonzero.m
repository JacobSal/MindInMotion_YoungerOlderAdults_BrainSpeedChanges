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
%-
DO_PLOT = true;
CLIM = [0,0.005];
% CLIM = [-0.005,0.005];
ALPHA = 0.05;
PLOT_CI = true;
FREQSCALE = 'log';
CONN_MEAS = 'dDTF08';
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
COND_NAMES = unique(EEG.etc.conn_table.t_fNames);
%% ===================================================================== %%
%## LOAD BOOSTRAP & PHASE RANDOMIZED INFORMATION
save_dir_nonzmats =  [save_dir filesep 'nonzero_mat_files'];
save_dir_bootmats =  [save_dir filesep 'bootstrap_mat_files'];
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
%% ===================================================================== %%
%## BOOTSTRAP EXTRACTION
CONN_MEAS = BOOTSTRAP_STRUCT(1).conn_meas;
FREQ_BANDS = BOOTSTRAP_STRUCT(1).FREQ_BAND_INF.FREQ_BANDS;
cond_1 = 1;
cond_2 = 2;
extract_mat_1 = zeros(F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
extract_mat_2 = zeros(F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
BOOT_CONN_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
BOOT_CONN_AVG_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
BOOT_CONN_MASK_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),F_DIM,T_DIM,length(BOOTSTRAP_STRUCT));
BOOT_CONN_FREQ_STORE = nan(length(MAIN_STUDY.cluster),length(MAIN_STUDY.cluster),T_DIM,length(fieldnames(FREQ_BANDS)),length(BOOTSTRAP_STRUCT));
for subj_i = 1:length(BOOTSTRAP_STRUCT)
    fprintf('assigning subject data %s...\n',BOOTSTRAP_STRUCT.subject);
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    %- time x frequency information
    BOOT_CONN_MASK_STORE(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).masked_conn{cond_1,cond_2};
    BOOT_CONN_STORE(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).masked_conn{cond_1,cond_2};
end
%% ===================================================================== %%
%## CONNECTIVITY MATRICIES
% freq_1 = 3;
% cl_1 = 1;
% cl_2 = 2;
cond_1 = 1;
cond_2 = 2;
%-
freq_conn_1 = BOOT_CONN_FREQ_STORE;
freq_conn_2 = BOOT_CONN_FREQ_STORE;
for subj_i = 1:length(BOOTSTRAP_STRUCT)
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    %-
    freq_conn_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).FREQ_BAND_INF.integrated_cause{cond_1,cond_2};
    freq_conn_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).FREQ_BAND_INF.integrated_cause{cond_2,cond_1};
end
%##
freq_names = fieldnames(FREQ_BANDS);
freq_bands = FREQ_BANDS;
COLOR_LIMITS=[-0.0001,0.0001];
boot_mat = freq_conn_1;
cond_i = 1;
% boot_mat = freq_conn_2;
% cond_i = 2;
cluster_ints = [3,7,5,4];
cluster_names = {'RPPa-Oc','LPPa-Oc','Precuneus','Cuneus'}; % (06/27/2023) JS, unsure on these as of yet.

for freq_i = 1:size(freq_conn_1,4)
    FREQS_SUBPATH = sprintf('%s',freq_names{freq_i});
    save_sub_figs = [save_dir filesep 'cluster_mats' filesep FREQS_SUBPATH];
    if ~exist(save_sub_figs,'dir')
        mkdir(save_sub_figs);
    end
    %- plot
    cnt = 1;
    clusterNames = {};
    for i = 1:length(MAIN_STUDY.cluster)
        N = length(MAIN_STUDY.cluster(i).sets);
        if any(i == cluster_ints)
            %- assign name if available
            idx = find(i == cluster_ints);
            clusterNames{idx} = sprintf('(N=%i) %s',N,cluster_names{(i == cluster_ints)});
            cnt = cnt + 1;            
        end
    end
    %## PLOT
    %- average acroos time
    tmp = squeeze(nanmean(squeeze(boot_mat(:,:,:,freq_i,:)),3));
    %- average across subject
    tmp = nanmean(tmp,3);
    I = eye(size(tmp));
    I = (I == 0);
    tmp = tmp.*I;
    tmp(tmp == 0) = nan();
%     tmp = log(tmp);
    %- delte unused clusters
    tmp = tmp(cluster_ints,:);
    tmp = tmp(:,cluster_ints);
    %- plot
    figure;
    hnd = heatmap(tmp,'Colormap',linspecer); %,'CellLabelColor', 'None');
    %- create title
    out1 = strsplit(COND_NAMES{cond_i},'.');
    out1 = strsplit(out1{1},'_');
    out1 = strjoin(out1(3:end),' ');
    title(sprintf('(%s) ''%s'' Connectivity Mean Across Clusters',FREQS_SUBPATH,out1));
    hnd.YDisplayLabels = clusterNames;
    hnd.XDisplayLabels = clusterNames;
    hnd.ColorLimits = COLOR_LIMITS;
    hnd.GridVisible = 'off';
    hnd.CellLabelFormat = '%0.1g';
    hnd.NodeChildren(3).Title.Interpreter = 'none';
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'Position',[10,100,720,620])
    exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s_mean_mat_nobase_%s.jpg','clusters',out1)],'Resolution',300);
%         exportgraphics(fig_i,[save_sub_figs filesep sprintf('%s_meanMatbased_%s.pdf','clusters',conn_conds{cond_i})],'Resolution',300);


end
%% ===================================================================== %%
%## SIGNAL COMPARISONS
% freq_1 = 1;
cl_1 = 7;
cl_2 = 3;
cond_1 = 1;
cond_2 = 2;
%-
freq_conn_1 = BOOT_CONN_FREQ_STORE;
freq_conn_2 = BOOT_CONN_FREQ_STORE;
for subj_i = 1:length(BOOTSTRAP_STRUCT)
    %- get components and cluster assignments
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    %-
    freq_conn_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).FREQ_BAND_INF.integrated_cause{cond_1,cond_2};
    freq_conn_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).FREQ_BAND_INF.integrated_cause{cond_2,cond_1};
end
COLORS_MEAN = linspecer(size(freq_conn_1,4))-0.15;
COLORS_STD = COLORS_MEAN+0.15;
%- mean of frequency band across time
%{
temp_mean_avg = squeeze(BOOT_CONN_FREQ_STORE(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,freq_1,:));
temp_mean_avg = nanmean(temp_mean_avg,2);
figure;
plot(temp_mean_avg);
fig_i = figure();
hold on;
%}
%- create title
in1 =CLUSTER_ASSIGNMENTS{cl_1}; %CLUSTER_ITERS(i)};
in2 =CLUSTER_ASSIGNMENTS{cl_2}; %CLUSTER_ITERS(j)};
freq_names = fieldnames(FREQ_BANDS);
figure();
hold on;
title(sprintf('Connectivity Across Time between %s & %s',in1,in2));
for freq_i = 1:size(freq_conn_1,4)
%     fprintf('here\n');
    tmp = squeeze(freq_conn_1(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,freq_i,:));
%     tmp = squeeze(freq_conn_2(CLUSTER_ITERS(cl_1),CLUSTER_ITERS(cl_2),:,freq_i,:));
%     tmp = cat(1,tmp{:});
    y_mean = nanmean(tmp,2); %rand(1,10); % your mean vector;
    y_med = nanmedian(tmp);
    x = (1:size(tmp,1)); %1:numel(y);
    xx = (1:size(tmp,1)); %1:numel(y);
    std_dev = nanstd(tmp,[],2);
    curve1 = y_mean + std_dev;
%     x(curve1==0) = [];
%     curve1(curve1==0) = [];
%     interp1(x,curve1,xx,'linear')
    curve2 = y_mean - std_dev;
    x2 = [x; fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
%     pp = fill(x2', inBetween, COLORS_STD(freq_i,:));
%     pp.FaceAlpha = 0.25;
%     pp.EdgeAlpha = 0.25;
%     plot(x,curve1);
%     plot(x,curve2);
    pp = plot(x,y_mean,'LineWidth',2,'Color',[COLORS_MEAN(freq_i,:),0.5],...
        'DisplayName',sprintf('%s',freq_names{freq_i}));
end
legend();
hold off;
% legend(tmp_leg{:})
% exportgraphics(fig_i,[save_sub_figs filesep sprintf('across_time_plot_%s_%s.jpg',in1,in2)],'Resolution',300);

%% ===================================================================== %%
%## time-frequency statistics?
cond_1 = 1;
cond_2 = 2;
cl_1 = 3;
cl_2 = 4;
STUDY = MAIN_STUDY;
boot_conn_1 = BOOT_CONN_AVG_STORE;
boot_conn_2 = BOOT_CONN_AVG_STORE;
for subj_i = 1:length(BOOTSTRAP_STRUCT)
    comps = squeeze(comps_out(:,subj_i));
    [tmpcl,idxcl] = sort(comps);
    idxcl = idxcl(tmpcl~=0);
    tmpcl = tmpcl(tmpcl~=0);
    boot_conn_1(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).averages{cond_1,1}; %.(CONN_MEAS);
    boot_conn_2(idxcl,idxcl,:,:,subj_i) = BOOTSTRAP_STRUCT(subj_i).averages{cond_2,1}; %.(CONN_MEAS
end
allersp = {squeeze(boot_conn_1(cl_1,cl_2,:,:,:));squeeze(boot_conn_2(cl_1,cl_2,:,:,:))};
% allfreqs = BOOTSTRAP_STRUCT(1).averages{1,1}.freqs;
allfreqs = FREQS_VEC;
alltimes = MAIN_ALLEEG(1).etc.COND_CAT(1).Conn.erWinCenterTimes; %BOOTSTRAP_STRUCT(1).averages{1,1}.erWinCenterTimes;
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
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
    'groupstats',ERSP_STAT_PARAMS.groupstats,...
    'method',ERSP_STAT_PARAMS.method,...
    'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
    'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
    'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
    'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
%- get stats parameters
stats = STUDY.etc.statistics;
stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
% if isempty(STUDY.design(STUDY.currentdesign).variable)
%     stats.paired = { };
% else
%     stats.paired = { STUDY.design(STUDY.currentdesign).variable(:).pairing };
% end
stats.paired = {'on'}; %{'BM','Human'};
%- get ersp params
params = STUDY.etc.erspparams;
params.plottf = [];
% select specific time and freq
% -----------------------------
if ~isempty(params.plottf)
    if length(params.plottf) < 3
        params.plottf(3:4) = params.plottf(2);
        params.plottf(2)   = params.plottf(1);
    end
    [~, fi1] = min(abs(allfreqs-params.plottf(1)));
    [~, fi2] = min(abs(allfreqs-params.plottf(2)));
    [~, ti1] = min(abs(alltimes-params.plottf(3)));
    [~, ti2] = min(abs(alltimes-params.plottf(4)));
    for index = 1:length(allersp(:))
        allersp{index} = mean(mean(allersp{index}(fi1:fi2,ti1:ti2,:,:),1),2);
        allersp{index} = reshape(allersp{index}, [1 size(allersp{index},3) size(allersp{index},4) ]);
    end

    % prepare channel neighbor matrix for Fieldtrip
    statstruct = std_prepare_neighbors(STUDY, ALLEEG);
    stats.fieldtrip.channelneighbor = statstruct.etc.statistics.fieldtrip.channelneighbor;

    params.plottf = { params.plottf(1:2) params.plottf(3:4) };
    [pcond, pgroup, pinter] = std_stat(allersp, stats);
    if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end % single subject STUDY
else
    [pcond, pgroup, pinter] = std_stat(allersp, stats);
    if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
            (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1))
        pcond = {}; pgroup = {}; pinter = {};
        disp('No statistics possible for single subject STUDY');
    end % single subject STUDY
end
%% ===================================================================== %%
% cl_1 = 3;
% cl_2 = 4;
%- mean of bootstrap comparison 
% temp_mean_tval = squeeze(BOOT_CONN_STORE(cl_1,cl_2,:,:,:));
% temp_mean_tval(temp_mean_tval==0) = nan();
% temp_mean_tval = nanmean(temp_mean_tval,3);
% temp_mean_tval(isnan(temp_mean_tval)) = 0;
%- mean of boostraps mean
% temp_mean_avg = BOOT_CONN_AVG_STORE(cl_1,cl_2,:,:,:);
% temp_mean_avg(temp_mean_avg==0) = nan();
% temp_mean_avg = nanmean(temp_mean_avg,3);
% temp_mean_avg(isnan(temp_mean_avg)) = 0;
%- masking consolidation?
% temp_mean_avg = BOOT_CONN_AVG_STORE(cl_1,cl_2,:,:,:);
% temp_mean_avg(temp_mean_avg==0) = nan();
% temp_mean_avg = nanmean(temp_mean_avg,3);
% temp_mean_avg(isnan(temp_mean_avg)) = 0;
%- pvalue from paired cluster based stats
allpcond = pcond{1};
allersp_mean = cellfun(@(x) nanmean(x,3),allersp,'UniformOutput',false);
%-
FIGURE_POSITION = [10,100,620,480];
clim_max = [0,0.0005]; %[-0.001,0.001];
SUB_FREQ_LIMS = [3,60];
alltitles = {'HUMAN','BALL MACHINE'};
XTICK_LABEL = 'time (s)';
YTICK_LABEL = 'Frequency (Hz)';
FONT_SIZE = 12;
SUBPLOT_WIDTH = 0.15;
SUBPLOT_HEIGHT = 0.7;
SHIFT_AMNT = 0.175;
STATS_TITLE = 'CUSTOM STATS';
if length(clim_max) == 2
    clim_ersp = clim_max;
else
    clim_ersp = [-clim_max,clim_max];
end
%##
fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
set(fig,'Units','inches','Position',[3 3 14 5])
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
%     for i = 1:length(warping_times)
%         xline(ax,warping_times(i),'k--');
%     end
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
%     xrng = get(ax,'XLim');
%     if warping_times(1) < xrng(1)
%         warping_times(1) = xrng(1);
%     end
%     if warping_times(end) > xrng(end)
%         warping_times(end) = xrng(end);
%     end
%     set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
%     xtickangle(45)
%     ax.XAxis.FontSize = FONT_SIZE;
%         ylim()
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
%     for i = 1:length(warping_times)
%         xline(ax,warping_times(i),'k--');
%     end
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
    %- set x-axis labels
%     xrng = get(ax,'XLim');
%     if warping_times(1) < xrng(1)
%         warping_times(1) = xrng(1);
%     end
%     if warping_times(end) > xrng(end)
%         warping_times(end) = xrng(end);
%     end
%     set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
%     xtickangle(45)
%     ax.XAxis.FontSize = FONT_SIZE;
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
%% ===================================================================== %%
if ~ispec
    addpath(convertPath2UNIX('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
else
    addpath(convertPath2Drive('M:\jsalminen\GitHub\par_EEGProcessing\submodules\groupSIFT'));
end
%% ===================================================================== %%
% input1 = squeeze(bootstrapped_out.matrix{1}(1,:,:,:));
% input2 = squeeze(bootstrapped_out.matrix{1}(2,:,:,:));
ind_1 = 1;
ind_2 = 3;
ind_3 = 4;
ind_4 = 5;
cl_1 = CLUSTER_ITERS(ind_1);
cl_2 = CLUSTER_ITERS(ind_2);
cl_3 = CLUSTER_ITERS(ind_3);
cl_4 = CLUSTER_ITERS(ind_4);
input_1 = BOOT_CONN_STORE(cl_1,cl_2,:,:,:);
input_2 = BOOT_CONN_STORE(cl_3,cl_4,:,:,:);
t_test_type = 1;
pval_preselc = 100; %0.05;
n_iters = 2000;
[mask,tscore,pValue] = clusterLevelPermutationTest(input1,input2,t_test_type,pval_preselc,n_iters);

mask_in = pValue<0.05;
mean_mat = mean(input1,3);
contour_in = mean_mat.*mask_in;
%##
figure();
colormap(linspecer)
contourf(contour_in)
colorbar();
%% CLUSTER_LEVEL_PERMUTATION =========================================== %%
cl_pairs = [];
for cl_i = main_cl_inds(2:end)
    for cl_j = main_cl_inds(2:end)
        if any(all(ismember(cl_pairs,[cl_i,cl_j]),2)) || cl_i==cl_j
            continue;
        end
        cl_pairs = [cl_pairs; cl_i,cl_j];
    end
end
PERMUTE_STATS_OUT = zeros(length(STUDY.cluster),length(STUDY.cluster),length(FREQ_BANDS));
for ind_1 = 1:length(CLUSTER_ITERS)
    for ind_2 = 1:length(CLUSTER_ITERS)
        for ind_3 = 1:length(CLUSTER_ITERS)
            for ind_4 = 1:length(CLUSTER_ITERS)
                cl_1 = CLUSTER_ITERS(ind_1);
                cl_2 = CLUSTER_ITERS(ind_2);
                cl_3 = CLUSTER_ITERS(ind_3);
                cl_4 = CLUSTER_ITERS(ind_4);
                input_1 = BOOT_CONN_STORE(cl_1,cl_2);
                input_2 = BOOT_CONN_STORE(cl_3,cl_4);
                t_test_type = 1;
                pval_preselc = 100; %0.05;
                n_iters = 2000;
                [mask,tscore,pValue] = clusterLevelPermutationTest(input1,input2,t_test_type,pval_preselc,n_iters);

                mask_in = pValue<0.05;
                mean_mat = mean(input1,3);
                contour_in = mean_mat.*mask_in;
            end
        end
    end
end
%##
figure();
colormap(linspecer)
contourf(contour_in)
colorbar();