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
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        STUDY_DIR = getenv('STUDY_DIR');
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
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
% [tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(save_dir,COND_CHARS);
conditions = {'walking_rest_all'};
[tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table([study_fpath filesep 'conn_valid_slide'],conditions);
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
%% ===================================================================== %%
connectivityMeasure = 'dDTF08'; % Alternatively, 'RPDC'
freqVec   = EEG.CAT.Conn.freqs;
deltaFreq = mean(diff(freqVec));
percentileThreshold = 95; 
EEG(end) = pop_vis_TimeFreqGrid(EEG(end),GUI_MODE, ...
                        'plotCondDiff',false,   ...
                        'vismode','TimeXFrequency', ...
                        'MatrixLayout', ...
                            {'Partial' ...
                             'triu' connectivityMeasure 'ut_clim' 99 ... % Was 100 (04/18/2021 Makoto)
                             'tril' connectivityMeasure 'lt_clim' 99 ... % Was 100 (04/18/2021 Makoto)
                             'diag' 'none' 'd_clim' 100   ...
                             'clim' 99},  ...
                         'clim',100,        ...
                         'timeRange',[],    ...
                         'freqValues', freqVec(1):deltaFreq:freqVec(end),   ... % Was [1 40] (04/18/2021 Makoto)
                         'windows',[],      ...
                         'pcontour',[],     ...
                         'thresholding',    ...
                            {'Simple'       ...
                            'prcthresh' [percentileThreshold 3]  ...
                            'absthresh' []},    ...
                        'baseline',[] , ... % Was [-1 -0.25] (04/18/2021 Makoto)
                        'fighandles',[],        ...
                        'smooth',false,         ...
                        'xord',[],'yord',[],    ...
                        'plotorder',[],         ...
                        'topoplot','dipole',    ...
                        'topoplot_opts',{},     ...
                        'customTopoMatrix',[],  ...
                        'dipplot',  ...
                            {'mri' '' 'coordformat' 'mni' 'dipplotopt' {}}, ...
                        'nodelabels',ComponentNames,        ...
                        'foilines',[3 7 15],    ...
                        'foilinecolor',[0.7 0.7 0.7] ,  ...
                        'events',{{0 'r' ':' 2}},       ...
                        'freqscale','linear',           ...
                        'transform','linear',           ...
                        'yTickLoc','right',             ...
                        'titleString','',               ...
                        'titleFontSize',12,             ...
                        'axesFontSize',11,              ...
                        'textColor',[1 1 1] ,           ...
                        'linecolor',[1 1 1] ,           ...
                        'patchcolor',[1 1 1] ,          ...
                        'colormap',jet(64),             ...
                        'backgroundColor',[0 0 0]);