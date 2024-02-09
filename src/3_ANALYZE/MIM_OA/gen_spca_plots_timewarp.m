%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_plots_timewarp.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
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
%- define the directory to the src folderd
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_OA'];
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
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('oa');
% fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
% fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
% cluster_study_dir = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
% cluster_study_dir = '01232023_MIM_YAN32_antsnormalize_iccREMG0p4_powpow0p3_conn';
cluster_study_dir = '01232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3';
% cluster_study_dir = '12012023_OAYA104_icc0p65-0p4_changparams';
study_fName_1 = 'epoch_study';
spca_study_dir = '01122024_spca_analysis';
study_fName_2 = 'epoch_study';
%- study group and saving
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir) filesep 'spca'];
load_dir_1 = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_2 = [STUDIES_DIR filesep sprintf('%s',spca_study_dir)];
% OUTSIDE_DATA_DIR = [STUDIES_DIR filesep ica_orig_dir]; % JACOB,SAL(02/23/2023)
%- load cluster
CLUSTER_DIR = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir) filesep 'cluster'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
CLUSTER_STUDY_DIR = [CLUSTER_DIR filesep 'icrej_5'];
CLUSTER_K = 12;
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[CLUSTER_STUDY_DIR filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[CLUSTER_STUDY_DIR filesep sprintf('%s.study',CLUSTER_STUDY_FNAME)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
% [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(TMP_STUDY);
% fPaths = {TMP_STUDY.datasetinfo.filepath};
% fNames = {TMP_STUDY.datasetinfo.filename};
% condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {STUDY.datasetinfo.subject};
%-
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
group_id = zeros(length(subject_chars),1);
for subj_i = 1:length(STUDY.datasetinfo)
    tmp = regexp(subject_chars{subj_i},'\d','match');
    group_id(subj_i) = str2num(tmp{1});
end
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
%% ===================================================================== %%

%* ERSP PARAMS
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
    'freqrange',[1,200]);
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
%## ersp plot per cluster per condition
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%%
%##
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%##
condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% condition_pairs = {{'flat','low','med','high'},...
%     {'0p25','0p5','0p75','1p0'}};
%## LOAD ICATIMEF OPTS
%- option 1
% icatimef_f = [TMP_STUDY.datasetinfo(1).filepath filesep sprintf('%s.icatimef',TMP_STUDY.datasetinfo(1).subject)];
% fprintf('Loading Resting Time-Frequency Data...\n');
% tmp = load(icatimef_f,'-mat');
%- option 2 (doesn't work if you delete .icatimef after spca creation)
% icatimef_f = [load_dir_2 filesep TMP_STUDY.datasetinfo(1).subject filesep 'GAIT_EPOCHED' filesep  sprintf('%s.icatimef',TMP_STUDY.datasetinfo(1).subject)];
% fprintf('Loading Resting Time-Frequency Data...\n');
% tmp = load(icatimef_f,'-mat');
%- option 3 (hardcode)
% timef_params = {'cycles',[3,0.800],'padratio',1,'alpha',NaN,...
%     'freqscale','log','nfreqs',200,'ntimesout',159,'timewarp',[],...
%     'timewarpms',[0,259,703,977,1408],'freqs',[3,250],...
%     'plotersp','off','plotitc','off','plotphase','off','baseline',NaN};
%- option 4
icatimef_f = [STUDY.datasetinfo(1).filepath filesep sprintf('%s.icatimef',STUDY.datasetinfo(1).subject)];
[~, timef_params, hardcode_times, hardcode_freqs, ~ ] = std_readfile( icatimef_f,'components',1);
%##
% allfreqs = 1:size(spca_table.tf_erspcorr_c{1},2); %4:100; %1:size(allersp_com{1},2);
% alltimes = 1:size(spca_table.tf_erspcorr_c{1},1);
COLORS_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
% hardcode_freqs = [3.00000000000000,3.06742258053202,3.13636042918590,3.20684760039063,3.27891891392104,3.35260997209831,3.42795717737705,3.50499775032772,3.58376974802305,3.66431208283782,3.74666454167101,3.83086780560009,3.91696346997695,4.00499406497544,4.09500307660079,4.18703496817112,4.28113520228174,4.37735026326317,4.47572768014374,4.57631605012836,4.67916506260494,4.78432552369029,4.89184938132775,5.00178975094877,5.11420094171129,5.22913848332777,5.34665915349618,5.46682100594745,5.58968339912332,5.71530702549861,5.84375394156257,5.97508759847400,6.10937287340531,6.24667610159107,6.38706510909672,6.53060924632382,6.67737942226828,6.82744813954852,6.98088953022080,7.13777939239961,7.29819522770088,7.46221627952689,7.62992357221147,7.80139995104498,7.97673012319891,8.15600069957009,8.33930023756540,8.52671928484803,8.71835042406688,8.91428831859121,9.11462975927315,9.31947371226118,9.52892136788816,9.74307619065805,9.96204397035612,10.1859328743077,10.4148535008116,10.6489189337742,10.8882447985713,11.1329493191659,11.3831533765093,11.6389805682547,11.9005572698126,12.1680126967792,12.4414789687669,12.7210911746699,13.0069874393964,13.2993089921002,13.5982002359469,13.9038088194464,14.2162857093901,14.5357852654259,14.8624653163106,15.1964872378750,15.5380160327415,15.8872204118333,16.2442728777155,16.6093498098094,16.9826315515215,17.3643024993308,17.7545511938786,18.1535704131050,18.5615572674787,18.9787132973675,19.4052445725961,19.8413617942425,20.2872803987215,20.7432206642077,21.2094078194496,21.6860721550307,22.1734491371293,22.6717795238361,23.1813094840861,23.7022907192622,24.2349805875331,24.7796422309847,25.3365447056091,25.9059631142147,26.4881787423239,27.0834791971242,27.6921585495426,28.3145174795132,28.9508634245091,29.6015107314125,30.2667808117985,30.9470023007079,31.6425112189893,32.3536511392884,33.0807733557695,33.8242370576498,34.5844095066342,35.3616662183387,36.1563911477894,36.9689768790923,37.7998248193646,38.6493453970245,39.5179582645380,40.4060925057219,41.3141868477056,42.2426898776570,43.1920602643787,44.1627669848850,45.1552895560700,46.1701182715835,47.2077544440297,48.2687106526091,49.3535109963264,50.4626913528889,51.5967996434230,52.7563961031407,53.9420535580883,55.1543577081158,56.3939074162048,57.6613150042995,58.9572065557859,60.2822222247693,61.6370165523021,63.0222587897190,64.4386332292388,65.8868395419959,67.3675931236693,68.8816254478788,70.4296844275240,72.0125347842437,73.6309584261788,75.2857548342250,76.9777414569664,78.7077541144847,80.4766474112439,82.2852951582543,84.1345908047237,86.0254478794103,87.9588004418943,89.9356035439920,91.9568337015388,94.0234893767758,96.1365914715780,98.2971838317667,100.506333762756,102.765132556788,105.074696032019,107.436165083718,109.850706247854,112.319512277352,114.843802731298,117.424824577382,120.063852807891,122.762191069532,125.521172307423,128.342159423546,131.226545950008,134.175756737426,137.191248658784,140.274511329112,143.427067841337,146.650475518672,149.946326683911,153.316249446019,156.761908504399,160.285005971229,163.887282212286,167.570516706663,171.336528925812,175.187179232337,179.124369798993,183.150045548333,187.266195113475,191.474851820462,195.778094692702,200.178049477977,204.676889698534,209.276837724781,213.980165873109,218.789197528387,223.706308291685,228.733927153790,233.874537695100,239.130679312479,244.504948473686,250.000000000000];
% hardcode_times = [58,82,106,130,156,180,204,230,254,278,302,328,352,376,402,426,450,474,500,524,548,574,598,622,646,672,696,720,746,770,794,820,844,868,892,918,942,966,992,1016,1040,1064,1090,1114,1138,1164,1188,1212,1236,1262,1286,1310,1336,1360,1384,1410,1434,1458,1482,1508,1532,1556,1582,1606,1630,1654,1680,1704,1728,1754,1778,1802,1826,1852,1876,1900,1926,1950,1974,2000,2024,2048,2072,2098,2122,2146,2172,2196,2220,2244,2270,2294,2318,2344,2368,2392,2416,2442,2466,2490,2516,2540,2564,2588,2614,2638,2662,2688,2712,2736,2762,2786,2810,2834,2860,2884,2908,2934,2958,2982,3006,3032,3056,3080,3106,3130,3154,3178,3204,3228,3252,3278,3302,3326,3352,3376,3400,3424,3450,3474,3498,3524,3548,3572,3596,3622,3646,3670,3696,3720,3744,3768,3794,3818,3842,3868,3892,3916,3942];
COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2];
COLOR_LIM_ERR = 0.05;
FREQ_BOUND = [4,60];
TIME_BOUND = [timef_params.timewarpms(1),timef_params.timewarpms(end)];
freq_crop = find(hardcode_freqs>=FREQ_BOUND(1) & hardcode_freqs<=FREQ_BOUND(2));
time_crop = find(hardcode_times>=TIME_BOUND(1) & hardcode_times<=TIME_BOUND(2));
% allfreqs = hardcode_freqs(freq_crop);
% alltimes = hardcode_times(time_crop);
psd_freqs = hardcode_freqs(freq_crop);
timewarp_times = [0,249,723,985,1449];
timewarp_chars = {'RHS','LTO','LHS','RTO','RHS'};
cond_iters = {1:4,5:8};
COLOR_PRCTILE= [15,95];
TERRAIN_REF_CHAR = 'flat';
SPEED_REF_CHAR = '1p0';
SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
ff = @chk_cell;
PLOT_STRUCT = struct('figure_position_inch',[3,3,14,6],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',timewarp_times,...
    'xticklabel_chars',{timewarp_chars},...
    'clim',[],...
    'font_size',12,...
    'font_name','Arial',...
    'freq_lims',FREQ_BOUND,...
    'time_lims',TIME_BOUND,...
    'subplot_width',0.15,...
    'subplot_height',0.65,...
    'shift_amnt',0.175,...
    'stats_title','F Statistic Mask',...
    'figure_title','');
SAVE_STATS = false;
%%
spca_table = par_load(CLUSTER_STUDY_DIR,'spca_cluster_table.mat');
%%
paramsersp = timef_params; 
alltimes = hardcode_times(time_crop);
allfreqs = hardcode_freqs(freq_crop);
CLUSTER_PICKS = main_cl_inds(2:end);
parfor (ii = 1:length(CLUSTER_PICKS),SLURM_POOL_SIZE)
    %## PARAMS
    cl_i = CLUSTER_PICKS(ii)
    alltimes = hardcode_times(time_crop);
    allfreqs = hardcode_freqs(freq_crop);
    PLOT_STRUCT_PAR = PLOT_STRUCT;
    TMP_STUDY = STUDY;
    %## ANATOMY
%     dip1 = mean(TMP_STUDY.cluster(cl_i).all_diplocs);
%     dip1 = TMP_STUDY.cluster(cl_i).dipole.posxyz;
    dip1 = TMP_STUDY.cluster(cl_i).all_diplocs;
%     dip1 = ([TMP_STUDY.cluster(cl_i).all_diplocs;TMP_STUDY.cluster(cl_i+CLUSTER_K+1).all_diplocs]);
%     dip1 = mean([TMP_STUDY.cluster(cl_i).all_diplocs;TMP_STUDY.cluster(cl_i+CLUSTER_K+1).all_diplocs]);
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    cfg              = [];
    cfg.roi        = dip1;
    cfg.output     = 'single';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if isempty(sub_indx)
                atlas_name = names{indx};
            end
        else
            atlas_name = names{indx};
        end
    end
    fprintf('Cluster %i) Anatomy Label: %s\n',cl_i,atlas_name);
    %%
    for des_i = 1:length(cond_iters)
        %% REPOP STATS
        TMP_STUDY = pop_statparams(TMP_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
            'groupstats',ERSP_STAT_PARAMS.groupstats,...
            'method',ERSP_STAT_PARAMS.method,...
            'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
            'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
            'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
            'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
        %%
        con_con = cond_iters{des_i};
        allersp = cell(length(con_con),1); 
        allgpm = cell(length(con_con),1); 
        subjs = zeros(length(con_con),1);
        trial_order = cell(length(con_con),1);
        allersp_orig = cell(length(con_con),1);
        allgpm_orig = cell(length(con_con),1);
        spec_ersp = cell(length(con_con),1);
        spec_gpm = cell(length(con_con),1);
        spec_ersp_orig = cell(length(con_con),1);
        spec_gpm_orig = cell(length(con_con),1);
        p1 = cell(length(con_con),1);
        p2 = cell(length(con_con),1);
        p3 = cell(length(con_con),1);
        p4 = cell(length(con_con),1);
        cnt = 1;
        for cond_i = 1:length(con_con)
            %- get cluster & condition indices
            inds_cl = cellfun(@(x) ff(x,cl_i),spca_table.cluster_c);
            inds_cond = strcmp(spca_table.cond_c,condition_gait{con_con(cond_i)});
            trial_order{cond_i} = condition_gait{con_con(cond_i)};
            inds = inds_cl & inds_cond;
            fprintf('True subjects in cluster: %i, Alg subjects in cluster: %i\n',length(TMP_STUDY.cluster(cl_i).sets),sum(inds))
            %- extract ersp
            tmp = cat(3,spca_table.tf_erspcorr_c{inds});
            tmp = permute(tmp,[2,1,3]);
            allersp{cond_i} = tmp;
%             allersp{cond_i} = tmp(time_crop,freq_crop,:);
            spec_ersp{cond_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
            %- extract gpm
            tmp = cat(3,spca_table.tf_gpmcorr_c{inds});
            tmp = permute(tmp,[2,1,3]);
            allgpm{cond_i} = tmp;
%             allgpm{cond_i} = tmp(time_crop,freq_crop,:);
            spec_gpm{cond_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
            %- extract original ersp
            tmp = cat(3,spca_table.tf_ersporig_c{inds});
            tmp = permute(tmp,[2,1,3]);
            allersp_orig{cond_i} = tmp;
%             allersp_orig{cond_i} = tmp(time_crop,freq_crop,:);
            spec_ersp_orig{cond_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
            %- extract original gpm
            tmp = cat(3,spca_table.tf_gpmorig_c{inds});
            tmp = permute(tmp,[2,1,3]);
            allgpm_orig{cond_i} = tmp;
%             allgpm_orig{cond_i} = tmp(time_crop,freq_crop,:);
            spec_gpm_orig{cond_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
        end
        %##
%         alltitles = std_figtitle('condnames',trial_order, ...
%             'clustname', sprintf('%s',atlas_name));
        if any(strcmp(trial_order,SPEED_REF_CHAR))
            condnames = SPEED_OVERRIDE_CHARS;
        else
            condnames = trial_order;
        end
%         alltitles = std_figtitle('condnames',condnames, ...
%                 'clustname', sprintf('CL%i',cl_i));
        alltitles = std_figtitle('condnames',condnames);
        %% SUBJECT-SPECIFIC WITHIN CONDITION BASELINE
        [allersp_sb_f,allersp_sb,~,~] = eeglab_baseln(allersp,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        [allgpm_sb_f,allgpm_sb,~,~] = eeglab_baseln(allgpm,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        [allerspo_sb_f,allerspo_sb,~,~] = eeglab_baseln(allersp_orig,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        [allgpmo_sb_f,allgpmo_sb,~,~] = eeglab_baseln(allgpm_orig,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
        %## COMMON BASELINE
        [allersp_com_f,allersp_com,~,~] = eeglab_baseln(allersp,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        [allgpm_com_f,allgpm_com,~,~] = eeglab_baseln(allgpm,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        [allerspo_com_f,allerspo_com,~,~] = eeglab_baseln(allersp_orig,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        [allgpmo_com_f,allgpmo_com,~,~] = eeglab_baseln(allgpm_orig,hardcode_times,hardcode_freqs,...
            TIME_BOUND,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',false);
        %## CROP NON-BASELINED FOR PLOTTINGS
        allersp_f = allersp;
        allgpm_f = allgpm;
        allersp_orig_f = allersp_orig;
        allgpm_orig_f = allgpm_orig;
        allersp = cellfun(@(x) x(freq_crop,time_crop,:),allersp,'uniformoutput',false);
        allgpm = cellfun(@(x) x(freq_crop,time_crop,:),allgpm,'uniformoutput',false);
        allersp_orig = cellfun(@(x) x(freq_crop,time_crop,:),allersp_orig,'uniformoutput',false);
        allgpm_orig = cellfun(@(x) x(freq_crop,time_crop,:),allgpm_orig,'uniformoutput',false);
        %% NO BASELINE
%         pcond_ersp_crop = {{}};
%         pcond_gpm_crop = {{}};
%         pcond_ersporig_crop = {{}};
%         pcond_gpmorig_crop = {{}};
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm,allfreqs,alltimes);
        [pcond_ersporig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp_orig,allfreqs,alltimes);
        [pcond_gpmorig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm_orig,allfreqs,alltimes);
        %- calculate per condition means
        for cond_i = 1:length(con_con)
            p1{cond_i} = squeeze(mean(allersp{cond_i},3));
            p2{cond_i} = squeeze(mean(allgpm{cond_i},3));
            p3{cond_i} = squeeze(mean(allersp_orig{cond_i},3));
            p4{cond_i} = squeeze(mean(allgpm_orig{cond_i},3));
        end
        %##
        PLOT_STRUCT_PAR.alltitles = alltitles;
        PLOT_STRUCT_PAR.figure_title = 'ERSP corrected';
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_ersp.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM corrected';
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpm.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'ERSP original';
        fig = plot_txf_conds_tftopo(p3,alltimes,allfreqs,pcond_ersporig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_ersporig.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM original';
        fig = plot_txf_conds_tftopo(p4,alltimes,allfreqs,pcond_gpmorig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpmorig.tiff',cl_i,des_i)],'Resolution',900);
        %% WITHIN CONDITION BASELINED ERSPs
%         pcond_ersp_crop = {{}};
%         pcond_gpm_crop = {{}};
%         pcond_ersporig_crop = {{}};
%         pcond_gpmorig_crop = {{}};
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp_sb,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm_sb,allfreqs,alltimes);
        [pcond_ersporig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allerspo_sb,allfreqs,alltimes);
        [pcond_gpmorig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpmo_sb,allfreqs,alltimes);
        %-
        for cond_i = 1:length(con_con)
            p1{cond_i} = squeeze(mean(allersp_sb{cond_i},3));
            p2{cond_i} = squeeze(mean(allgpm_sb{cond_i},3));
            p3{cond_i} = squeeze(mean(allerspo_sb{cond_i},3));
            p4{cond_i} = squeeze(mean(allgpmo_sb{cond_i},3));
        end
        %##
        PLOT_STRUCT_PAR.alltitles = alltitles;
        PLOT_STRUCT_PAR.figure_title = 'ERSP corrected sub-base';
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_ersp_sb.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM corrected sub-base';
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpm_sb.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'ERSP original corrected sub-base';
        fig = plot_txf_conds_tftopo(p3,alltimes,allfreqs,pcond_ersporig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_erspo_sb.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM original corrected sub-base';
        fig = plot_txf_conds_tftopo(p4,alltimes,allfreqs,pcond_gpmorig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpmo_sb.tiff',cl_i,des_i)],'Resolution',900);
        %% ACROSS CONDITIONS BASELINED ERSPS
%         pcond_ersp_crop = {{}};
%         pcond_gpm_crop = {{}};
%         pcond_ersporig_crop = {{}};
%         pcond_gpmorig_crop = {{}};
        [pcond_ersp_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allersp_com,allfreqs,alltimes);
        [pcond_gpm_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpm_com,allfreqs,alltimes);
        [pcond_ersporig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allerspo_com,allfreqs,alltimes);
        [pcond_gpmorig_crop, ~, ~] = ersp_stats_conds(TMP_STUDY,allgpmo_com,allfreqs,alltimes);
        %- subject specific plots
        for cond_i = 1:length(con_con)
            p1{cond_i} = squeeze(mean(allersp_com{cond_i},3));
            p2{cond_i} = squeeze(mean(allgpm_com{cond_i},3));
            p3{cond_i} = squeeze(mean(allerspo_com{cond_i},3));
            p4{cond_i} = squeeze(mean(allgpmo_com{cond_i},3));
        end
        %##
        PLOT_STRUCT_PAR.alltitles = alltitles;
        PLOT_STRUCT_PAR.figure_title = 'ERSP corrected common-base';
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_ersp_com.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM corrected common-base';
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpm_com.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'ERSP original corrected common-base';
        fig = plot_txf_conds_tftopo(p3,alltimes,allfreqs,pcond_ersporig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_erspo_sb.tiff',cl_i,des_i)],'Resolution',900);
        %##
        PLOT_STRUCT_PAR.figure_title = 'GPM original corrected common-base';
        fig = plot_txf_conds_tftopo(p4,alltimes,allfreqs,pcond_gpmorig_crop{1},...
            'PLOT_STRUCT',PLOT_STRUCT_PAR);
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_gpmo_sb.tiff',cl_i,des_i)],'Resolution',900);
        %% ============================================================= %%
        data_to_proc = {allersp_f,allgpm_f,allersp_orig_f,allgpm_orig_f,...
            allersp_sb_f,allgpm_sb_f,allerspo_sb_f,allgpmo_sb_f,...
            allersp_com_f,allgpm_com_f,allerspo_com_f,allgpmo_com_f};
        data_chars = {'allersp','allgpm','allersp_orig','allgpm_orig',...
            'allersp_sb','allgpm_sb','allerspo_sb','allgpmo_sb',...
            'allersp_com','allgpm_com','allerspo_com','allgpmo_com'};
%         if any(strcmp(trial_order,SPEED_REF_CHAR))
%             condnames = SPEED_OVERRIDE_CHARS;
%         else
%             condnames = trial_order;
%         end
%         alltitles = std_figtitle('condnames',condnames, ...
%                 'clustname', sprintf('CL%i',cl_i));
        for data_i = 1:length(data_to_proc)
            allersp_in = data_to_proc{data_i};
            chk_1 = strcmp(TERRAIN_REF_CHAR,trial_order);
            chk_2 = strcmp(SPEED_REF_CHAR,trial_order);
            cond_chars = trial_order;
            if any(chk_1)
                refErspCond = TERRAIN_REF_CHAR;
                refErspCond_fext = TERRAIN_REF_CHAR;
                refErspCond_ind = find(chk_1);
            elseif any(chk_2)
                refErspCond_ind = find(chk_2);
                refErspCond = SPEED_OVERRIDE_CHARS{refErspCond_ind};
                refErspCond_fext = SPEED_REF_CHAR;
            else
                error('Condition for reference ersp not found in STUDY design: %i',des_i)
            end
            inds_to_comp = setdiff(1:length(alltitles),refErspCond_ind);
            if ~isempty(refErspCond)
                % mask differenec ersps- check that it's sig. different from zer
                erspDiff = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
                erspDiff_wind = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
                %- calculate pairwise statistics between conditions of interest
                for c = inds_to_comp
                    %-
                    fprintf('Computing Pair Stat for %s - %s...\n',refErspCond,cond_chars{c})
                    curr_ersp = allersp_in{c,1};
                    ref_ersp = allersp_in{refErspCond_ind,1};
                    [tmp, ~, ~] = ersp_stats_conds(TMP_STUDY,{curr_ersp;ref_ersp},hardcode_freqs,hardcode_times);
                    erspDiff(c).raw = mean(curr_ersp-ref_ersp,3);
                    erspDiff(c).masked = erspDiff(c).raw.*tmp{1,1};
                    erspDiff(c).pcond = tmp{1,1};
                    %-
                    curr_ersp_wind = allersp_in{c,1}(freq_crop,time_crop,:);
                    ref_ersp_wind = allersp_in{refErspCond_ind,1}(freq_crop,time_crop,:);
                    [tmp, ~, ~] = ersp_stats_conds(TMP_STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
                    
                    erspDiff_wind(c).raw = mean(curr_ersp_wind-ref_ersp_wind,3);
                    erspDiff_wind(c).masked = erspDiff_wind(c).raw.*tmp{1,1};
                    erspDiff_wind(c).pcond = tmp{1,1};
                end
            end
            if SAVE_STATS
                if ~exist([save_dir filesep 'stats_out'])
                    mkdir([save_dir filesep 'stats_out'])
                end
                par_save(erspDiff_wind,[save_dir filesep 'stats_out'],sprintf('%i_%s_erspdiff_wind.mat',cl_i,data_chars{data_i}));
                par_save(erspDiff,[save_dir filesep 'stats_out'],sprintf('%i_%s_erspdiff.mat',cl_i,data_chars{data_i}));
            end
            %-
%             ersp_raw = {erspDiff.raw};
%             ersp_pcond = {erspDiff.pcond};
%             ersp_masked = {erspDiff.masked};
%             alltimes = hardcode_times;
%             allfreqs = hardcode_freqs;
            %-
            ersp_raw = {erspDiff_wind.raw};
            ersp_pcond = {erspDiff_wind.pcond};
            ersp_masked = {erspDiff_wind.masked};
            alltimes = hardcode_times(time_crop);
            allfreqs = hardcode_freqs(freq_crop);
            %-
            ersp_raw = ersp_raw(~cellfun(@isempty,ersp_raw));
            ersp_pcond = ersp_pcond(~cellfun(@isempty,ersp_pcond));
            ersp_masked = ersp_masked(~cellfun(@isempty,ersp_masked));
            plot_alltitles = cell(size(inds_to_comp));
            PLOT_STRUCT_PAR.figure_title = sprintf('%s based',data_chars{data_i});
            for j = 1:length(inds_to_comp)
                cc = inds_to_comp(j);
                plot_alltitles{j} = sprintf('%s - %s',alltitles{cc},refErspCond);
            end
            PLOT_STRUCT_PAR.alltitles = plot_alltitles;
            [fig] = plot_txf_mask_contourf(ersp_raw,alltimes,allfreqs,ersp_masked,ersp_pcond,...
                'PLOT_STRUCT',PLOT_STRUCT_PAR);
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spcadiff_%s.tiff',cl_i,des_i,data_chars{data_i})],'Resolution',900);
        end
        %% ============================================================= %%
        %## PARAMS
        %{
        try
            TMP_STUDY.etc = rmfield(TMP_STUDY.etc,'statistics');
        catch e
            fprintf('error. trouble removing field ''statistics'' from STUDY\n');
        end
        TMP_STUDY = pop_statparams(TMP_STUDY, 'groupstats','off','condstats', 'on',...
                    'method','perm',...
                    'singletrials','off','mode','fieldtrip','fieldtripalpha',0.05,...
                    'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',2000);
        stats = TMP_STUDY.etc.statistics;
        stats.paired{1} = 'on'; % Condition stats
        stats.paired{2} = 'off'; % Group stats
        %- 
        switch des_i
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.5','0.75','1.0'};
        end
        %## PSD STATS & FIGURE
        data_in = {spec_ersp,spec_gpm,spec_gpm_orig,spec_ersp_orig};
        data_in_char = {'ersp','gpm','gpmorig','ersporig'};
        for d_i = 1:length(data_in)
            data = data_in{d_i};
            SPEC_XLIM = [4,60];
            SPEC_YLIM = [prctile([data{:}],15,'all')-0.2,prctile([data{:}],95,'all')+0.2];
            [spec_pcond,spec_group,~,~,~,~] = std_stat(data, stats);
            %- figure
            fig = figure;
            hold on;
            for i = 1:length(data)
                y = data{i}';
                JackKnife_sung(psd_freqs,mean(y),[mean(y)-std(y)/sqrt(size(y,1))],[mean(y)+std(y)/sqrt(size(y,1))],...
                    color_dark(i,:),color_light(i,:));
                plot(psd_freqs,mean(y),'color',color_dark(i,:),'linewidth',2);
            end  

            ax = gca;
            xlim(SPEC_XLIM);ylim(SPEC_YLIM);
            axsignif = highlight_CL(ax, psd_freqs, spec_pcond{1}, 'background', 'Frequency(Hz)');
            plot(SPEC_XLIM,[0 0],'--','color','black');
            xlabel('Frequency(Hz)');ylabel('10*log_{10}(Power)');
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','bold')
            xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
            title(sprintf('PSD %s',data_in_char{d_i}));
            hold on;
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_psd_%s.jpg',cl_i,des_i,data_in_char{d_i})]);
        end
        %}
    end
end
%%

% hardcode_freqs = [3.00000000000000,3.06742258053202,3.13636042918590,3.20684760039063,3.27891891392104,3.35260997209831,3.42795717737705,3.50499775032772,3.58376974802305,3.66431208283782,3.74666454167101,3.83086780560009,3.91696346997695,4.00499406497544,4.09500307660079,4.18703496817112,4.28113520228174,4.37735026326317,4.47572768014374,4.57631605012836,4.67916506260494,4.78432552369029,4.89184938132775,5.00178975094877,5.11420094171129,5.22913848332777,5.34665915349618,5.46682100594745,5.58968339912332,5.71530702549861,5.84375394156257,5.97508759847400,6.10937287340531,6.24667610159107,6.38706510909672,6.53060924632382,6.67737942226828,6.82744813954852,6.98088953022080,7.13777939239961,7.29819522770088,7.46221627952689,7.62992357221147,7.80139995104498,7.97673012319891,8.15600069957009,8.33930023756540,8.52671928484803,8.71835042406688,8.91428831859121,9.11462975927315,9.31947371226118,9.52892136788816,9.74307619065805,9.96204397035612,10.1859328743077,10.4148535008116,10.6489189337742,10.8882447985713,11.1329493191659,11.3831533765093,11.6389805682547,11.9005572698126,12.1680126967792,12.4414789687669,12.7210911746699,13.0069874393964,13.2993089921002,13.5982002359469,13.9038088194464,14.2162857093901,14.5357852654259,14.8624653163106,15.1964872378750,15.5380160327415,15.8872204118333,16.2442728777155,16.6093498098094,16.9826315515215,17.3643024993308,17.7545511938786,18.1535704131050,18.5615572674787,18.9787132973675,19.4052445725961,19.8413617942425,20.2872803987215,20.7432206642077,21.2094078194496,21.6860721550307,22.1734491371293,22.6717795238361,23.1813094840861,23.7022907192622,24.2349805875331,24.7796422309847,25.3365447056091,25.9059631142147,26.4881787423239,27.0834791971242,27.6921585495426,28.3145174795132,28.9508634245091,29.6015107314125,30.2667808117985,30.9470023007079,31.6425112189893,32.3536511392884,33.0807733557695,33.8242370576498,34.5844095066342,35.3616662183387,36.1563911477894,36.9689768790923,37.7998248193646,38.6493453970245,39.5179582645380,40.4060925057219,41.3141868477056,42.2426898776570,43.1920602643787,44.1627669848850,45.1552895560700,46.1701182715835,47.2077544440297,48.2687106526091,49.3535109963264,50.4626913528889,51.5967996434230,52.7563961031407,53.9420535580883,55.1543577081158,56.3939074162048,57.6613150042995,58.9572065557859,60.2822222247693,61.6370165523021,63.0222587897190,64.4386332292388,65.8868395419959,67.3675931236693,68.8816254478788,70.4296844275240,72.0125347842437,73.6309584261788,75.2857548342250,76.9777414569664,78.7077541144847,80.4766474112439,82.2852951582543,84.1345908047237,86.0254478794103,87.9588004418943,89.9356035439920,91.9568337015388,94.0234893767758,96.1365914715780,98.2971838317667,100.506333762756,102.765132556788,105.074696032019,107.436165083718,109.850706247854,112.319512277352,114.843802731298,117.424824577382,120.063852807891,122.762191069532,125.521172307423,128.342159423546,131.226545950008,134.175756737426,137.191248658784,140.274511329112,143.427067841337,146.650475518672,149.946326683911,153.316249446019,156.761908504399,160.285005971229,163.887282212286,167.570516706663,171.336528925812,175.187179232337,179.124369798993,183.150045548333,187.266195113475,191.474851820462,195.778094692702,200.178049477977,204.676889698534,209.276837724781,213.980165873109,218.789197528387,223.706308291685,228.733927153790,233.874537695100,239.130679312479,244.504948473686,250.000000000000];
% hardcode_times = [58,82,106,130,156,180,204,230,254,278,302,328,352,376,402,426,450,474,500,524,548,574,598,622,646,672,696,720,746,770,794,820,844,868,892,918,942,966,992,1016,1040,1064,1090,1114,1138,1164,1188,1212,1236,1262,1286,1310,1336,1360,1384,1410,1434,1458,1482,1508,1532,1556,1582,1606,1630,1654,1680,1704,1728,1754,1778,1802,1826,1852,1876,1900,1926,1950,1974,2000,2024,2048,2072,2098,2122,2146,2172,2196,2220,2244,2270,2294,2318,2344,2368,2392,2416,2442,2466,2490,2516,2540,2564,2588,2614,2638,2662,2688,2712,2736,2762,2786,2810,2834,2860,2884,2908,2934,2958,2982,3006,3032,3056,3080,3106,3130,3154,3178,3204,3228,3252,3278,3302,3326,3352,3376,3400,3424,3450,3474,3498,3524,3548,3572,3596,3622,3646,3670,3696,3720,3744,3768,3794,3818,3842,3868,3892,3916,3942];
%##
% fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
% exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
%##
