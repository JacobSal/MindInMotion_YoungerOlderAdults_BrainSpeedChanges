%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/run_spca_dd_tw_plots_gclim.sh

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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
ica_orig_dir = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
spca_study_dir = '03232024_spca_analysis_OA';
% ICLABEL_EYE_CUTOFF = 0.75;
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
% save_dir = [studies_fpath filesep sprintf('%s',study_dir_name)];
load_dir = [studies_fpath filesep sprintf('%s',spca_study_dir)];
OUTSIDE_DATA_DIR = [studies_fpath filesep ica_orig_dir]; % JACOB,SAL(02/23/2023)
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt'];
% cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%- create new study directory
save_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K) filesep 'spca_ersp_plots'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
%- extract study chars

% condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {STUDY.datasetinfo.subject};
group_id = zeros(length(subject_chars),1);
for subj_i = 1:length(STUDY.datasetinfo)
    tmp = regexp(subject_chars{subj_i},'\d','match');
    group_id(subj_i) = str2num(tmp{1});
end
%-
% fPaths = {STUDY.datasetinfo.filepath};
% fNames = {STUDY.datasetinfo.filename};
%% ===================================================================== %%
%-
ERSP_STAT_PARAMS_COND = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
%-
ERSP_STAT_PARAMS_GROUP = ERSP_STAT_PARAMS_COND;
ERSP_STAT_PARAMS_GROUP.groupstats = 'on';
ERSP_STAT_PARAMS_GROUP.condstats = 'off';
%- 
ERSP_STAT_PARAMS_GC = ERSP_STAT_PARAMS_COND;
ERSP_STAT_PARAMS_GC.groupstats = 'on';
%-
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable1','cond','values1',{'flat','low','med','high'},...
%             'variable2','group','values2',{}},...
%             {'subjselect',{},...
%             'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
%             'variable2','group','values2',{}}};
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'}}};
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},'pairing','on',...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},'pairing','on',...
%             'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},'pairing','on',...
            'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},'pairing','on',...
            'variable2','group','values2',{'H1000''s','H2000''s','H3000''s'}}};
% condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); 
% % (09/18/2024) JS, reorders it weirdly. Manually override
%## ersp plot per cluster per condition
args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
STUDY = pop_statparams(STUDY,args{:});
args = eeglab_struct2args(ERSP_PARAMS);
STUDY = pop_erspparams(STUDY,args{:});
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,[],des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%%
%##
ATLAS_PATH = [PATHS.submods_dir filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
%##
CLUSTER_PICKS = main_cl_inds(2:end);
% ERSP_CLIM_OVERRIDE = {[-.5,.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],...
%     [-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],...
%     [-0.5,0.5],[-0.5,0.5]};
% cluster_titles = {'Precuneus','Right Posterior Parietal',...
%     'Left Occipital','Left Supplementary Motor','Left Sensorimotor','Left Posterior Parietal',...
%     'Eye','Left Temporal','Mid/Posterior Cingulate','Right Sensorimotor','Right Temporal'};
%-
cluster_titles = {'Right Occipital','Left Occipital','Mid Cingulate',...
    'Right Sensorimotor','Left Supplementary','Precuneus','Left Temporal','Left Sensorimotor',...
    'Right Posterior Parietal','Left Posterior Parietal','Right Temporal'};
ERSP_CLIM_OVERRIDE = {[-1,1],[-1,1],[-1,1],...
    [-1,1],[-1,1],[-1,1],[-1,1],[-1,1],...
    [-1,1],[-1,1],[-1,1]};
GPM_CLIM_OVERRIDE = {[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],...
    [-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],...
    [-0.5,0.5],[-0.5,0.5],[-0.5,0.5]};
%##
groups_ind = [1,2,3];
groups = unique({STUDY.datasetinfo.group});
% group_chars = {'Young Adult',{'Older Adult','High Mobility'},{'Older Adult','Low Mobility'}};
group_chars = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
% g_chars_topo = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
%## WEIRD WAY TO BE CONSISTENT WITH CONDITION NAME CONVENTION
tmp_load_chars = cell(1,length(unique({STUDY.datasetinfo(1).trialinfo.cond})));
cnt = 1;
for des_i = length(STUDY.design):-1:1
    for i = 1:length(STUDY.design(des_i).variable)
        if strcmp(STUDY.design(des_i).variable(i).label,'cond')
            tmp_load_chars(1,cnt:length(STUDY.design(des_i).variable(i).value)+cnt-1) = STUDY.design(des_i).variable(i).value;
            cnt = cnt + length(STUDY.design(des_i).variable(i).value);
        end
    end
end
%## LOAD ICATIMEF OPTS
%- option 5 (update?)
tmpf = par_load([load_dir filesep sprintf('%s',STUDY.datasetinfo(1).subject) filesep 'GAIT_EPOCHED' filesep [tmp_load_chars{:}]],'gait_ersp_spca.mat');
timef_params = tmpf.icatimefopts;
timef_params.timewarpms = tmpf.warptimes;
tmpf = par_load([load_dir filesep sprintf('%s',STUDY.datasetinfo(1).subject) filesep 'GAIT_EPOCHED' filesep [tmp_load_chars{:}]],'condmed_spca_ersp.mat');
hardcode_times = tmpf.times;
hardcode_freqs = tmpf.freqs;
%## hard define
BOOT_NITERS = 2000;
BOOT_ALPHA = 0.05;
BOOT_CLUST_THRESH = 150;
%## ERSP PLOT PARAMS
COLORS_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
% hardcode_freqs = [3.00000000000000,3.06742258053202,3.13636042918590,3.20684760039063,3.27891891392104,3.35260997209831,3.42795717737705,3.50499775032772,3.58376974802305,3.66431208283782,3.74666454167101,3.83086780560009,3.91696346997695,4.00499406497544,4.09500307660079,4.18703496817112,4.28113520228174,4.37735026326317,4.47572768014374,4.57631605012836,4.67916506260494,4.78432552369029,4.89184938132775,5.00178975094877,5.11420094171129,5.22913848332777,5.34665915349618,5.46682100594745,5.58968339912332,5.71530702549861,5.84375394156257,5.97508759847400,6.10937287340531,6.24667610159107,6.38706510909672,6.53060924632382,6.67737942226828,6.82744813954852,6.98088953022080,7.13777939239961,7.29819522770088,7.46221627952689,7.62992357221147,7.80139995104498,7.97673012319891,8.15600069957009,8.33930023756540,8.52671928484803,8.71835042406688,8.91428831859121,9.11462975927315,9.31947371226118,9.52892136788816,9.74307619065805,9.96204397035612,10.1859328743077,10.4148535008116,10.6489189337742,10.8882447985713,11.1329493191659,11.3831533765093,11.6389805682547,11.9005572698126,12.1680126967792,12.4414789687669,12.7210911746699,13.0069874393964,13.2993089921002,13.5982002359469,13.9038088194464,14.2162857093901,14.5357852654259,14.8624653163106,15.1964872378750,15.5380160327415,15.8872204118333,16.2442728777155,16.6093498098094,16.9826315515215,17.3643024993308,17.7545511938786,18.1535704131050,18.5615572674787,18.9787132973675,19.4052445725961,19.8413617942425,20.2872803987215,20.7432206642077,21.2094078194496,21.6860721550307,22.1734491371293,22.6717795238361,23.1813094840861,23.7022907192622,24.2349805875331,24.7796422309847,25.3365447056091,25.9059631142147,26.4881787423239,27.0834791971242,27.6921585495426,28.3145174795132,28.9508634245091,29.6015107314125,30.2667808117985,30.9470023007079,31.6425112189893,32.3536511392884,33.0807733557695,33.8242370576498,34.5844095066342,35.3616662183387,36.1563911477894,36.9689768790923,37.7998248193646,38.6493453970245,39.5179582645380,40.4060925057219,41.3141868477056,42.2426898776570,43.1920602643787,44.1627669848850,45.1552895560700,46.1701182715835,47.2077544440297,48.2687106526091,49.3535109963264,50.4626913528889,51.5967996434230,52.7563961031407,53.9420535580883,55.1543577081158,56.3939074162048,57.6613150042995,58.9572065557859,60.2822222247693,61.6370165523021,63.0222587897190,64.4386332292388,65.8868395419959,67.3675931236693,68.8816254478788,70.4296844275240,72.0125347842437,73.6309584261788,75.2857548342250,76.9777414569664,78.7077541144847,80.4766474112439,82.2852951582543,84.1345908047237,86.0254478794103,87.9588004418943,89.9356035439920,91.9568337015388,94.0234893767758,96.1365914715780,98.2971838317667,100.506333762756,102.765132556788,105.074696032019,107.436165083718,109.850706247854,112.319512277352,114.843802731298,117.424824577382,120.063852807891,122.762191069532,125.521172307423,128.342159423546,131.226545950008,134.175756737426,137.191248658784,140.274511329112,143.427067841337,146.650475518672,149.946326683911,153.316249446019,156.761908504399,160.285005971229,163.887282212286,167.570516706663,171.336528925812,175.187179232337,179.124369798993,183.150045548333,187.266195113475,191.474851820462,195.778094692702,200.178049477977,204.676889698534,209.276837724781,213.980165873109,218.789197528387,223.706308291685,228.733927153790,233.874537695100,239.130679312479,244.504948473686,250.000000000000];
% hardcode_times = [58,82,106,130,156,180,204,230,254,278,302,328,352,376,402,426,450,474,500,524,548,574,598,622,646,672,696,720,746,770,794,820,844,868,892,918,942,966,992,1016,1040,1064,1090,1114,1138,1164,1188,1212,1236,1262,1286,1310,1336,1360,1384,1410,1434,1458,1482,1508,1532,1556,1582,1606,1630,1654,1680,1704,1728,1754,1778,1802,1826,1852,1876,1900,1926,1950,1974,2000,2024,2048,2072,2098,2122,2146,2172,2196,2220,2244,2270,2294,2318,2344,2368,2392,2416,2442,2466,2490,2516,2540,2564,2588,2614,2638,2662,2688,2712,2736,2762,2786,2810,2834,2860,2884,2908,2934,2958,2982,3006,3032,3056,3080,3106,3130,3154,3178,3204,3228,3252,3278,3302,3326,3352,3376,3400,3424,3450,3474,3498,3524,3548,3572,3596,3622,3646,3670,3696,3720,3744,3768,3794,3818,3842,3868,3892,3916,3942];
% COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2];
% COLOR_LIM_ERR = 0.05;
FREQ_BOUND = [4,60];
% FREQ_BOUND_PSD = [4,100];
TIME_BOUND = [timef_params.timewarpms(1),timef_params.timewarpms(end)];
freq_crop = find(hardcode_freqs>=FREQ_BOUND(1) & hardcode_freqs<=FREQ_BOUND(2));
% freq_crop_psd = find(hardcode_freqs>=FREQ_BOUND_PSD(1) & hardcode_freqs<=FREQ_BOUND_PSD(2));
time_crop = find(hardcode_times>=TIME_BOUND(1) & hardcode_times<=TIME_BOUND(2));
% psd_freqs = hardcode_freqs(freq_crop_psd);
timewarp_times = [0,249,723,985,1449];
timewarp_chars = {'RHS','LTO','LHS','RTO','RHS'};
cond_iters = {1:4,5:8};
% COLOR_PRCTILE= [15,95];
SAVE_RESOLUTION = 1000;
TERRAIN_REF_CHAR = 'flat';
% SPEED_REF_CHAR = '1p0';
SPEED_REF_CHAR = '1p0';
SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
ff = @chk_cell;
PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Events',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',timewarp_times,...
    'xticklabel_chars',{timewarp_chars},...
    'xticklabel_angle',45,...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[],...
    'time_lims',[],...
    'stats_title','F Stat (p<0.05)',...
    'figure_title','',...
    'contourf_grain',ceil((500/pi())),...
    'group_titles',{{}},...
    'group_titles_shift_x',0.0,...
    'group_titles_shift_y',0.65,...
    'subplot_width',0.13,...
    'subplot_height',0.16,... %(02/17/2024) was 0.2
    'subplot_shift_x',0.035,...
    'subplot_shift_y',0.05,...
    'subplot_init_y',0.7,...
    'subplot_init_x',0.06,...
    'colorbar_shift_x',0.145,...
    'colorbar_shift_y',0,...
    'colorbar_label_shift_x',0,...
    'colorbar_label_shift_y',0.005,...
    'colorbar_label','\Delta Power (dB)',...
    'colorbar_fontsize',8,...
    'colorbar_fontweight','bold',...
    'display_bandmarks',true,...
    'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
    'bandmarks_shift_y',[-0.42,-0.2,0.05,0.325],...
    'bandmarks_shift_x',-0.03,...
    'bandmarks_fontsize',8,...
    'bandmarks_fontweight','bold');
SAVE_STATS = false;
%%
% spca_table = par_load(cluster_study_fpath,'spca_cluster_table.mat');
spca_table = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],'spca_cluster_table.mat');
spca_table = spca_table(~all(cellfun(@isempty,spca_table{:,:}),2),:);
% paramsersp = timef_params; 
alltimes = hardcode_times(time_crop);
allfreqs = hardcode_freqs(freq_crop);

% (09/18/2024) JS, reorders it weirdly. Manually override
%##
% condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};

%%
parfor (ii = 1:length(CLUSTER_PICKS),SLURM_POOL_SIZE)
% for ii = 1:length(CLUSTER_PICKS)
    %## PARAMS
    cl_i = CLUSTER_PICKS(ii);
    alltimes = hardcode_times(time_crop);
    allfreqs = hardcode_freqs(freq_crop);
    tmp_plot_struct = PLOT_STRUCT;
    tmp_study = STUDY;
    c_chars = {};
    g_chars = {};
    g_inds = [];
    allerspo_sb_f = {};
    refErspCond = {};
    erspDiff = {};
    erspDiff_wind = {};
    ersp_raw = {};
    ersp_pcond ={};
    %% ANATOMY
%     dip1 = TMP_STUDY.cluster(cl_i).all_diplocs;
%     atlas_name = 'error';
%     cfg              = [];
%     cfg.roi        = dip1;
%     cfg.output     = 'single';
%     cfg.atlas      = atlas;
%     cfg.inputcoord = 'mni';
%     %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
%     cfg.sphere = 3;
%     label_i = ft_volumelookup(cfg, atlas);
%     if ~isempty(label_i)
%         counts = sum([label_i.count],2);
%         [val, indx] = max(counts);
%         names = label_i(1).name;
%         if strcmp(names(indx),'no_label_found')
%             sub_indx = find(counts ~= 0 & counts < val);
%             if isempty(sub_indx)
%                 atlas_name = names{indx};
%             end
%         else
%             atlas_name = names{indx};
%         end
%     end
%     fprintf('Cluster %i) Anatomy Label: %s\n',cl_i,atlas_name);
    %%
    STORAGE = cell(length(cond_iters),20);
    for des_i = 1:length(cond_iters)
        %## REPOP STATS
        %-
        args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
        % args = eeglab_struct2args(ERSP_STAT_PARAMS_GROUP);
        tmp_study = pop_statparams(tmp_study,args{:});
        %-
%         args = eeglab_struct2args(ERSP_STAT_PARAMS_GC);
%         TMP_STUDY = pop_statparams(TMP_STUDY,args{:});
        %##
        % c_chars = cond_iters{des_i};
        s_chars = {STUDY.datasetinfo(STUDY.cluster(cl_i).sets).subject};
        for i = 1:length(STUDY.design(des_i).variable)
            if strcmp(STUDY.design(des_i).variable(i).label,'cond')
                c_chars = STUDY.design(des_i).variable(i).value;
            elseif strcmp(STUDY.design(des_i).variable(i).label,'group')
                g_chars = STUDY.design(des_i).variable(i).value;
            end
        end
        allersp = cell(length(c_chars),length(groups_ind)); 
        allgpm = cell(length(c_chars),length(groups_ind)); 
        subjs = zeros(length(c_chars),length(groups_ind));
        trial_order = cell(length(c_chars),length(groups_ind));
        allersp_orig = cell(length(c_chars),length(groups_ind));
        allgpm_orig = cell(length(c_chars),length(groups_ind));
        spec_ersp = cell(length(c_chars),length(groups_ind));
        spec_gpm = cell(length(c_chars),length(groups_ind));
        spec_ersp_orig = cell(length(c_chars),length(groups_ind));
        spec_gpm_orig = cell(length(c_chars),length(groups_ind));
        cnt = 1;
        for group_i = 1:length(groups_ind)
            try
                if isnan(g_chars)
                    g_inds = (1:length(STUDY.cluster(cl_i).sets));
                else
                    fprintf('g_chars has some weird value in it\n')
                end
            catch
                g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});                   
            end
            %-
            cl_chars = s_chars(g_inds);
            for cond_i = 1:length(c_chars)
                %- get cluster & condition indices
                inds_cl = cellfun(@(x) ff(x,cl_i),spca_table.cluster_c);
                inds_cond = strcmp(spca_table.cond_c,c_chars{cond_i});
                inds_grp = cellfun(@(x) x == groups_ind(group_i),spca_table.group_c);
                inds = inds_cl & inds_cond & inds_grp;
                %-
                g_inds = cellfun(@(x) strcmp(x,g_chars{groups_ind(group_i)}),{tmp_study.datasetinfo(tmp_study.cluster(cl_i).sets).group});           
                % fprintf('CL%i) True subjects of group %s in cluster: %i, Alg subjects in cluster: %i\n',cl_i,groups{groups_ind(group_i)},sum(g_inds),sum(inds))
                fprintf('CL%i) Condition: %s, Group: %s\n\t Cluster N=%i, Subjects Added=%i\n',...
                    cl_i,c_chars{cond_i},g_chars{group_i},sum(g_inds),sum(inds));
                % fprintf('CL%i) STUDY Subjects (N=%i): %s\n',cl_i,length(cl_chars),sprintf('%s,',cl_chars{:}));
                % fprintf('CL%i) SPCA Subjects (N=%i): %s\n',cl_i,sum(inds),sprintf('%s,',spca_table.subj_c{inds}));
                subj_miss = [setdiff(cl_chars,spca_table.subj_c(inds)'),setdiff(spca_table.subj_c(inds)',cl_chars)];
                fprintf('CL%i) Subjects not in agreement between SPCA & STUDY (N=%i):\n %s\n',cl_i,length(subj_miss),sprintf('%s,',subj_miss{:}));
                trial_order(cond_i,group_i) = c_chars(cond_i); %unique({spca_table.cond_c{inds}}); %condition_gait{c_chars(cond_i)};
                %- extract ersp
                tmp = cat(3,spca_table.tf_erspcorr_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allersp{cond_i,group_i} = tmp;
                spec_ersp{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
                %- extract gpm
                tmp = cat(3,spca_table.tf_gpmcorr_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allgpm{cond_i,group_i} = tmp;
                spec_gpm{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
                %- extract original ersp
                tmp = cat(3,spca_table.tf_ersporig_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allersp_orig{cond_i,group_i} = tmp;
                spec_ersp_orig{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
                %- extract original gpm
                tmp = cat(3,spca_table.tf_gpmorig_c{inds});
                tmp = permute(tmp,[2,1,3]);
                allgpm_orig{cond_i,group_i} = tmp;
                spec_gpm_orig{cond_i,group_i} = squeeze(mean(tmp(freq_crop,time_crop,:),1));
            end
        end
        %## SUBJECT-SPECIFIC WITHIN CONDITION BASELINE
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
        %## STORE
        STORAGE{des_i,1} = c_chars;
        STORAGE{des_i,2} = subjs;
        STORAGE{des_i,3} = trial_order;
        %-
        STORAGE{des_i,4} = allersp;
        STORAGE{des_i,5} = allgpm;
        STORAGE{des_i,6} = allersp_orig;
        STORAGE{des_i,7} = allgpm_orig;
        %-
        STORAGE{des_i,8} = allersp_f;
        STORAGE{des_i,9} = allgpm_f;
        STORAGE{des_i,10} = allersp_orig_f;
        STORAGE{des_i,11} = allgpm_orig_f;
        %-
        STORAGE{des_i,12} = allersp_com;
        STORAGE{des_i,13} = allgpm_com;
        STORAGE{des_i,14} = allersp_com_f;
        STORAGE{des_i,15} = allgpm_com_f;
        %-
        STORAGE{des_i,16} = allersp_sb;
        STORAGE{des_i,17} = allgpm_sb;
        STORAGE{des_i,18} = allersp_sb_f;
        STORAGE{des_i,19} = allgpm_sb_f;
    end
    %% CLIM
    % HIGH = 99;
    % LOW = 100-HIGH;
    %## GPM PLOTS
    % INT=16;
    % data = [];
    % for des_i = 1:length(cond_iters)
    %     tmp = STORAGE{des_i,INT};
    %     tmp = cellfun(@(x) mean(x,3),tmp,'UniformOutput',false);
    %     data = cat(3,data,tmp{:});
    % end
    % bound = max([abs(prctile(data,LOW,'all')), abs(prctile(data,HIGH,'all'))]);
    % GPM_CLIM = sort([-round(bound,1),round(bound,1)]);
    GPM_CLIM = GPM_CLIM_OVERRIDE{ii};
    %## ERSP PLOTS
    % INT=12;
    % data = [];
    % for des_i = 1:length(cond_iters)
    %     tmp = STORAGE{des_i,INT};
    %     tmp = cellfun(@(x) mean(x,3),tmp,'UniformOutput',false);
    %     data = cat(3,data,tmp{:});
    % end
    % bound = max([abs(prctile([data],LOW,'all')),abs(prctile([data],HIGH,'all'))]);
    % ERSP_CLIM = sort([-round(bound,1),round(bound,1)]);
    ERSP_CLIM = ERSP_CLIM_OVERRIDE{ii};
    %%
    for des_i = 1:length(cond_iters)
        %## VARS
        c_chars = STORAGE{des_i,1};
        subjs = STORAGE{des_i,2};
        trial_order = STORAGE{des_i,3};
        %-
        allersp = STORAGE{des_i,4};
        allgpm = STORAGE{des_i,5};
        allersp_orig = STORAGE{des_i,6};
        allgpm_orig = STORAGE{des_i,7};
        %-
        allersp_f  = STORAGE{des_i,8};
        allgpm_f = STORAGE{des_i,9};
        allersp_orig_f = STORAGE{des_i,10};
        allgpm_orig_f = STORAGE{des_i,11};
        %-
        allersp_com = STORAGE{des_i,12};
        allgpm_com = STORAGE{des_i,13};
        allersp_com_f = STORAGE{des_i,14};
        allgpm_com_f = STORAGE{des_i,15};
        %-
        allersp_sb = STORAGE{des_i,16};
        allgpm_sb = STORAGE{des_i,17};
        allersp_sb_f = STORAGE{des_i,18};
        allgpm_sb_f = STORAGE{des_i,19};
        %-
        p1 = cell(length(c_chars),1);
        p2 = cell(length(c_chars),1);
        % p3 = cell(length(c_chars),1);
        % p4 = cell(length(c_chars),1);
        %##
        if any(strcmp(trial_order,SPEED_REF_CHAR))
            condnames = SPEED_OVERRIDE_CHARS;
        else
            condnames = trial_order;
        end
        alltitles = std_figtitle('condnames',condnames);
        %% BOOTSTRAPING ALLERSP W/ SUBBASELINE (GOOD TO LOOK AT)
        clust_ersp = cell(size(allersp_sb,1),size(allersp_sb,2));
        clust_maskedersp = cell(size(allersp_sb,1),size(allersp_sb,2));
        for group_i = 1:size(allersp_sb,2)
            for cond_i = 1:size(allersp_sb,1)
                fprintf('Performing Stats for Condition %i & Cluster %i\n',cond_i,cl_i);
                tmp = allersp_sb{cond_i,group_i};
                tmp_mean = mean(tmp,3);
                boot_freq = 1:size(tmp,1);
                boot_subj = 1:size(tmp,3);
                boot_surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
                surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
                %- scramble time samples and calculate the average across
                %all times and all frequencies and store that value.
                for n = 1:BOOT_NITERS
                    boot_time = randi(size(tmp,2),[size(tmp,2),1]); % random time samples
                    tmpSurro = mean(tmp(boot_freq,boot_time,boot_subj),3); % calculate mean across subjects for randomized time samples
                    surro(:,:,n) = tmpSurro; % save 2000 iterations of surrogates that represent completely random ERSP data 
                end
                %## VISUALIZE DISTRIBUTION?

                %- Pull length(subject) surrogate averages from distribution then calc mean across
                %surrogates 
                for n = 1:BOOT_NITERS
                    bootIdx  = randi(BOOT_NITERS,[size(tmp,3),1]); % randomly pick subjects
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
                clust_ersp{cond_i,group_i} = tmp_mean; 
                tmp = clust_ersp{cond_i,group_i}; 
                tmp(~finalMask) = 0;
                clust_maskedersp{cond_i,group_i} = tmp;
            end
        end
        % tmp_plot_struct.subplot_width = 0.13;
        % tmp_plot_struct.subplot_height = 0.16;
        % tmp_plot_struct.horiz_shift_amnt = 0.17;
        % tmp_plot_struct.subplot_shift_x = 0.22;
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.figure_title = cluster_titles{ii};
        tmp_plot_struct.clim = GPM_CLIM;
        tmp_plot_struct.group_titles = group_chars(groups_ind); %{groups{groups_ind},'Group Stats'};
        [fig] = plot_txf_mask_contourf(clust_ersp,alltimes,allfreqs,clust_maskedersp,clust_maskedersp,{},...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_group_bootstraps_ersp_sb.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        close(fig);
        tmp_plot_struct = PLOT_STRUCT;
        %%
        % p1 = cellfun(@(x) std(x,[],3),allersp_sb,'UniformOutput',false);
        % HIGH = 99.52;
        % LOW = 100-HIGH;
        % %## GPM PLOTS
        % INT=5;
        % data = cat(3,p1{:});
        % tmp_mask = cell(size(p1));
        % tmp_mask = cellfun(@(x) ones(size(p1{1})),tmp_mask,'UniformOutput',false)
        % bound = max([abs(prctile([data],LOW,'all')),abs(prctile([data],HIGH,'all'))]);
        % std_bound = std(data,[],'all');
        % % tmp_clim = sort([-round(bound,1)-0.5*std_bound,round(bound,1)+0.5*std_bound]);
        % tmp_clim = sort([prctile([data],LOW,'all')-0.5*std_bound,prctile([data],HIGH,'all')+0.5*std_bound]);
        % 
        % PLOT_STRUCT_PAR.subplot_width = 0.13;
        % PLOT_STRUCT_PAR.subplot_height = 0.16;
        % PLOT_STRUCT_PAR.horiz_shift_amnt = 0.17;
        % PLOT_STRUCT_PAR.vert_shift_amnt = 0.22;
        % PLOT_STRUCT_PAR.alltitles = alltitles;
        % PLOT_STRUCT_PAR.clim = tmp_clim;
        % PLOT_STRUCT_PAR.group_titles = {groups{groups_ind},'Group Stats'};
        % [fig] = plot_txf_mask_contourf(p1,alltimes,allfreqs,tmp_mask,tmp_mask,{},...
        %     'PLOT_STRUCT',PLOT_STRUCT_PAR);
        % drawnow;
        % exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_group_stdev_ersp_sb.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        % close(fig);
        % PLOT_STRUCT_PAR = PLOT_STRUCT;
        %%
        %{
        pcond_p1 = cell(size(allersp_sb,2),1);
        pgroup_p1 = cell(size(allersp_sb,1),1);
        for group_i = 1:size(allersp_sb,2)
            data = allersp_sb(:,group_i);
            [stat,cfg] = eeglab_ft_stats(data);
            pcond_p1{:,group_i} = squeeze(stat.prob);
            % pgroup_p1 = [];
        end
        %}
        %% ACROSS CONDITIONS BASELINED ERSPS
        [pcond_ersp_crop, pgroup_ersp_crop, ~] = ersp_stats_conds(tmp_study,allersp_com,allfreqs,alltimes);
        [pcond_gpm_crop, pgroup_gpm_crop, ~] = ersp_stats_conds(tmp_study,allersp_sb,allfreqs,alltimes);
        %- subject specific plots
        %- calculate per condition means
        p1 = cellfun(@(x) mean(x,3),allersp_com,'UniformOutput',false);
        p2 = cellfun(@(x) mean(x,3),allersp_sb,'UniformOutput',false);
        %- debug data
        %{
        pcond_ersp_crop = cell(1,3);
        pcond_ersp_crop = cellfun(@(x) zeros(size(p1{1})),pcond_ersp_crop,'UniformOutput',false);  
        %}
        %## CONTOURF (ERSP COM)
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = ERSP_CLIM;
        tmp_plot_struct.figure_title = cluster_titles{ii};
        %-
        ersp_masked = cellfun(@(x) x > 0.05, pcond_ersp_crop,'UniformOutput',false);
        ersp_masked = repmat(ersp_masked,[4,1]);
        pcond_ersp_crop = repmat(pcond_ersp_crop,[4,1]);
        [fig] = plot_txf_mask_contourf(p1,alltimes,allfreqs,ersp_masked,pcond_ersp_crop,{},...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_mask_groupersp_com.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        close(fig);
        tmp_plot_struct = PLOT_STRUCT;
        %## TFTOPO (ERSP COM)
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = ERSP_CLIM;
        tmp_plot_struct.figure_title = cluster_titles{ii};
        %-
        fig = plot_txf_conds_tftopo(p1,alltimes,allfreqs,pcond_ersp_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_fstat_groupersp_com.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        close(fig);
        tmp_plot_struct = PLOT_STRUCT;
        %## CONTOURF (ERSP SB)
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = GPM_CLIM;
        tmp_plot_struct.figure_title = cluster_titles{ii};
        %-
        ersp_masked = cellfun(@(x) x > 0.05, pcond_gpm_crop,'UniformOutput',false);
        ersp_masked = repmat(ersp_masked,[4,1]);
        pcond_gpm_crop = repmat(pcond_gpm_crop,[4,1]);
        [fig] = plot_txf_mask_contourf(p2,alltimes,allfreqs,ersp_masked,pcond_gpm_crop,{},...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_mask_groupersp_sb.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        close(fig);
        tmp_plot_struct = PLOT_STRUCT;
        %##
        tmp_plot_struct = PLOT_STRUCT;
        tmp_plot_struct.alltitles = alltitles;
        tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind},'Group Stats'};
        tmp_plot_struct.clim = GPM_CLIM;
        tmp_plot_struct.figure_title = cluster_titles{ii};
        %-
        fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_gpm_crop,...
            'PLOT_STRUCT',tmp_plot_struct);
        drawnow;
        exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spca_fstat_groupersp_sb.tiff',cl_i,des_i)],'Resolution',SAVE_RESOLUTION);
        close(fig);
        tmp_plot_struct = PLOT_STRUCT;
        %% ============================================================= %%
%         data_to_proc = {allersp_f,allgpm_f,allersp_orig_f,allgpm_orig_f,...
%             allersp_sb_f,allgpm_sb_f,allerspo_sb_f,allgpmo_sb_f,...
%             allersp_com_f,allgpm_com_f,allerspo_com_f,allgpmo_com_f};
%         data_chars = {'allersp','allgpm','allersp_orig','allgpm_orig',...
%             'allersp_sb','allgpm_sb','allerspo_sb','allgpmo_sb',...
%             'allersp_com','allgpm_com','allerspo_com','allgpmo_com'};
        data_to_proc = {allersp_com_f,allerspo_sb_f};
        data_chars = {'allersp_com','allersp_sb'};
        data_clim = {GPM_CLIM,GPM_CLIM};
        args = eeglab_struct2args(ERSP_STAT_PARAMS_COND);
%         args = eeglab_struct2args(ERSP_STAT_PARAMS_GC);
        tmp_study = pop_statparams(tmp_study,args{:});
        for data_i = 1:length(data_to_proc)
            allersp_in = data_to_proc{data_i};
            allersp_out = cell(size(allersp_in));
            allersp_pcond = cell(size(allersp_in));
            allersp_mask = cell(size(allersp_in));
            plot_alltitles = cell(size(allersp_in));
            for group_i = 1:size(allersp_in,2)
                chk_1 = any(strcmp(TERRAIN_REF_CHAR,trial_order),2);
                chk_2 = any(strcmp(SPEED_REF_CHAR,trial_order),2);
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
                    erspDiff = struct('raw',cell(length(inds_to_comp),1),...
                        'masked',cell(length(inds_to_comp),1),....
                        'pcond',cell(length(inds_to_comp),1),...
                        'pgroup',cell(length(inds_to_comp),1));
                    erspDiff_wind = struct('raw',cell(length(inds_to_comp),1),...
                        'masked',cell(length(inds_to_comp),1),...
                        'pcond',cell(length(inds_to_comp),1),...
                        'pgroup',cell(length(inds_to_comp),1));
                    %- calculate pairwise statistics between conditions of interest
                    for c = inds_to_comp
                        %-
                        fprintf('Computing Pair Stat for %s - %s...\n',refErspCond,cond_chars{c})
                        curr_ersp = allersp_in{c,group_i};
                        ref_ersp = allersp_in{refErspCond_ind,group_i};
                        [tmp,tmpg, ~] = ersp_stats_conds(tmp_study,{curr_ersp;ref_ersp},hardcode_freqs,hardcode_times);
                        erspDiff(c).raw = mean(curr_ersp-ref_ersp,3);
                        erspDiff(c).masked = erspDiff(c).raw.*tmp{1,1};
                        erspDiff(c).pcond = tmp{1,1};
%                         erspDiff(c).pgroup = tmpg{1,1};
                        %-
                        curr_ersp_wind = allersp_in{c,group_i}(freq_crop,time_crop,:);
                        ref_ersp_wind = allersp_in{refErspCond_ind,group_i}(freq_crop,time_crop,:);
                        [tmp, tmpg, ~] = ersp_stats_conds(tmp_study,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
                        erspDiff_wind(c).raw = mean(curr_ersp_wind-ref_ersp_wind,3);
                        erspDiff_wind(c).masked = erspDiff_wind(c).raw.*tmp{1,1};
                        erspDiff_wind(c).pcond = tmp{1,1};
%                         erspDiff_wind(c).pgroup = tmpg{1,1};
                    end
                end
                if SAVE_STATS
                    if ~exist([save_dir filesep 'stats_out'])
                        mkdir([save_dir filesep 'stats_out'])
                    end
                    par_save(erspDiff_wind,[save_dir filesep 'stats_out'],sprintf('%i_%s_grouperspdiff_wind.mat',cl_i,data_chars{data_i}));
                    par_save(erspDiff,[save_dir filesep 'stats_out'],sprintf('%i_%s_grouperspdiff.mat',cl_i,data_chars{data_i}));
                end
                %-
    %             ersp_raw = {erspDiff.raw};
    %             ersp_pcond = {erspDiff.pcond};
%                 ersp_pgroup = {erspDiff.pgroup};
    %             ersp_masked = {erspDiff.masked};
    %             alltimes = hardcode_times;
    %             allfreqs = hardcode_freqs;
                %-
                ersp_raw = {erspDiff_wind.raw}';
                ersp_pcond = {erspDiff_wind.pcond}';
                ersp_pgroup = {erspDiff_wind.pgroup};
                ersp_masked = {erspDiff_wind.masked}';
                alltimes = hardcode_times(time_crop);
                allfreqs = hardcode_freqs(freq_crop);
                %-
                ersp_raw = ersp_raw(~cellfun(@isempty,ersp_raw));
                ersp_pcond = ersp_pcond(~cellfun(@isempty,ersp_pcond));
                ersp_masked = ersp_masked(~cellfun(@isempty,ersp_masked));
                %-
                allersp_out(1:length(ersp_raw),group_i) = ersp_raw;
                allersp_pcond(1:length(ersp_raw),group_i) = ersp_pcond;
                allersp_mask(1:length(ersp_raw),group_i) = ersp_masked;
                for j = 1:length(inds_to_comp)
                    cc = inds_to_comp(j);
                    plot_alltitles{j,group_i} = sprintf('%s - %s',alltitles{cc},refErspCond);
                end
            end
            %##
            tmp1 = [];
            for i = 1:size(allersp_out,1)
                tmp1 = [tmp1; allersp_out(i,~cellfun(@isempty,allersp_out(i,:)))];
            end
            tmp2 = [];
            for i = 1:size(allersp_pcond,1)
                tmp2 = [tmp2; allersp_pcond(i,~cellfun(@isempty,allersp_pcond(i,:)))];
            end
            tmp3 = [];
            for i = 1:size(allersp_mask,1)
                tmp3 = [tmp3; allersp_mask(i,~cellfun(@isempty,allersp_mask(i,:)))];
            end
            tmp4 = [];
            for i = 1:size(plot_alltitles,1)
                tmp4 = [tmp4; plot_alltitles(i,~cellfun(@isempty,plot_alltitles(i,:)))];
            end
            % allersp_out = allersp_out(tmp);
            % allersp_pcond = allersp_pcond(~cellfun(@isempty,allersp_pcond));
            % allersp_mask = allersp_mask(~cellfun(@isempty,allersp_mask));
            % plot_alltitles = plot_alltitles(~cellfun(@isempty,plot_alltitles));
            % tmp_plot_struct.subplot_width = 0.13;
            % tmp_plot_struct.subplot_height = 0.16;
            % tmp_plot_struct.horiz_shift_amnt = 0.185;
            % tmp_plot_struct.subplot_shift_x = 0.22;
            % PLOT_STRUCT_PAR.alltitles = alltitles;
            % PLOT_STRUCT_PAR.clim = GPM_CLIM;
            % PLOT_STRUCT_PAR.group_titles = {groups{groups_ind},'Group Stats'};
            % tmp_plot_struct.figure_title = cluster_titles{ii};
            % tmp_plot_struct.alltitles = tmp4;
            % tmp_plot_struct.clim = data_clim{data_i};
            % tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind(group_i)},'Group Stats'};
            % tmp_plot_struct.colorbar_label_xshift = 0.14;
            % tmp_plot_struct.subplot_init_x  = 0.2;
            % tmp_plot_struct = PLOT_STRUCT;
            tmp_plot_struct.subplot_init_x = 0.25;
            tmp_plot_struct.subplot_shift_x = 0.05;
            tmp_plot_struct.alltitles = tmp4;
            tmp_plot_struct.group_titles = group_chars(groups_ind); %{group_chars{groups_ind},'Group Stats'};
            tmp_plot_struct.clim = data_clim{data_i};
            tmp_plot_struct.figure_title = cluster_titles{ii};
            [fig] = plot_txf_mask_contourf(tmp1,alltimes,allfreqs,tmp3,tmp2,{},...
                'PLOT_STRUCT',tmp_plot_struct);
            drawnow;
            % exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_group%s_spcadiff_%s.tiff',cl_i,des_i,groups{group_i},data_chars{data_i})],'Resolution',1000);
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_group_spcadiff_ref%s_%s.tiff',cl_i,des_i,refErspCond_fext,data_chars{data_i})],'Resolution',SAVE_RESOLUTION);
%             exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_spcadiff_%s.jpg',cl_i,des_i,data_chars{data_i})],'Resolution',300);
            close(fig);
            tmp_plot_struct = PLOT_STRUCT;
        end
    end
end
%%

% hardcode_freqs = [3.00000000000000,3.06742258053202,3.13636042918590,3.20684760039063,3.27891891392104,3.35260997209831,3.42795717737705,3.50499775032772,3.58376974802305,3.66431208283782,3.74666454167101,3.83086780560009,3.91696346997695,4.00499406497544,4.09500307660079,4.18703496817112,4.28113520228174,4.37735026326317,4.47572768014374,4.57631605012836,4.67916506260494,4.78432552369029,4.89184938132775,5.00178975094877,5.11420094171129,5.22913848332777,5.34665915349618,5.46682100594745,5.58968339912332,5.71530702549861,5.84375394156257,5.97508759847400,6.10937287340531,6.24667610159107,6.38706510909672,6.53060924632382,6.67737942226828,6.82744813954852,6.98088953022080,7.13777939239961,7.29819522770088,7.46221627952689,7.62992357221147,7.80139995104498,7.97673012319891,8.15600069957009,8.33930023756540,8.52671928484803,8.71835042406688,8.91428831859121,9.11462975927315,9.31947371226118,9.52892136788816,9.74307619065805,9.96204397035612,10.1859328743077,10.4148535008116,10.6489189337742,10.8882447985713,11.1329493191659,11.3831533765093,11.6389805682547,11.9005572698126,12.1680126967792,12.4414789687669,12.7210911746699,13.0069874393964,13.2993089921002,13.5982002359469,13.9038088194464,14.2162857093901,14.5357852654259,14.8624653163106,15.1964872378750,15.5380160327415,15.8872204118333,16.2442728777155,16.6093498098094,16.9826315515215,17.3643024993308,17.7545511938786,18.1535704131050,18.5615572674787,18.9787132973675,19.4052445725961,19.8413617942425,20.2872803987215,20.7432206642077,21.2094078194496,21.6860721550307,22.1734491371293,22.6717795238361,23.1813094840861,23.7022907192622,24.2349805875331,24.7796422309847,25.3365447056091,25.9059631142147,26.4881787423239,27.0834791971242,27.6921585495426,28.3145174795132,28.9508634245091,29.6015107314125,30.2667808117985,30.9470023007079,31.6425112189893,32.3536511392884,33.0807733557695,33.8242370576498,34.5844095066342,35.3616662183387,36.1563911477894,36.9689768790923,37.7998248193646,38.6493453970245,39.5179582645380,40.4060925057219,41.3141868477056,42.2426898776570,43.1920602643787,44.1627669848850,45.1552895560700,46.1701182715835,47.2077544440297,48.2687106526091,49.3535109963264,50.4626913528889,51.5967996434230,52.7563961031407,53.9420535580883,55.1543577081158,56.3939074162048,57.6613150042995,58.9572065557859,60.2822222247693,61.6370165523021,63.0222587897190,64.4386332292388,65.8868395419959,67.3675931236693,68.8816254478788,70.4296844275240,72.0125347842437,73.6309584261788,75.2857548342250,76.9777414569664,78.7077541144847,80.4766474112439,82.2852951582543,84.1345908047237,86.0254478794103,87.9588004418943,89.9356035439920,91.9568337015388,94.0234893767758,96.1365914715780,98.2971838317667,100.506333762756,102.765132556788,105.074696032019,107.436165083718,109.850706247854,112.319512277352,114.843802731298,117.424824577382,120.063852807891,122.762191069532,125.521172307423,128.342159423546,131.226545950008,134.175756737426,137.191248658784,140.274511329112,143.427067841337,146.650475518672,149.946326683911,153.316249446019,156.761908504399,160.285005971229,163.887282212286,167.570516706663,171.336528925812,175.187179232337,179.124369798993,183.150045548333,187.266195113475,191.474851820462,195.778094692702,200.178049477977,204.676889698534,209.276837724781,213.980165873109,218.789197528387,223.706308291685,228.733927153790,233.874537695100,239.130679312479,244.504948473686,250.000000000000];
% hardcode_times = [58,82,106,130,156,180,204,230,254,278,302,328,352,376,402,426,450,474,500,524,548,574,598,622,646,672,696,720,746,770,794,820,844,868,892,918,942,966,992,1016,1040,1064,1090,1114,1138,1164,1188,1212,1236,1262,1286,1310,1336,1360,1384,1410,1434,1458,1482,1508,1532,1556,1582,1606,1630,1654,1680,1704,1728,1754,1778,1802,1826,1852,1876,1900,1926,1950,1974,2000,2024,2048,2072,2098,2122,2146,2172,2196,2220,2244,2270,2294,2318,2344,2368,2392,2416,2442,2466,2490,2516,2540,2564,2588,2614,2638,2662,2688,2712,2736,2762,2786,2810,2834,2860,2884,2908,2934,2958,2982,3006,3032,3056,3080,3106,3130,3154,3178,3204,3228,3252,3278,3302,3326,3352,3376,3400,3424,3450,3474,3498,3524,3548,3572,3596,3622,3646,3670,3696,3720,3744,3768,3794,3818,3842,3868,3892,3916,3942];
%##
% fig = plot_txf(squeeze(PSC1(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
% exportgraphics(fig,[fpath filesep sprintf('chan%i_allcond_psc1_orig.jpg',CHAN_INT)]);
%##
