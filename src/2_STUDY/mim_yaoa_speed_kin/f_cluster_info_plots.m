%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_f_cluster_info_plots.sh

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
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = true;
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
DO_SUBJ_PLOTS = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result?
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% (08/03/2023) JS, changing fieldtripnaccu to 10,000; changing
% fieldtripmcorrect to cluster; method to perm; (these are the parameters
% set inside CL's PlotAndSaveERSP_CL_V3.m...
% pipeline although this doesn't align with her YA manuscript methods?
% (08/06/2023) JS, changing fieldtripnaccu to 2000 again and mcorrect to fdr...
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
%## Soft Define
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- datset name
DATA_SET = 'MIM_dataset';
%- cluster directory for study
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- study info
SUB_GROUP_FNAME = 'group_spec';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 12;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## LOAD STUDY
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
end
cl_struct = par_load([cluster_study_fpath filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {STUDY.datasetinfo.subject};
%-
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
CLUSTER_PICKS = main_cl_inds(2:end);
%## CREATE PATHS
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
if ~exist(spec_data_dir,'dir')
    mkdir(spec_data_dir);
    % error('spec_data dir does not exist');
end
%## LOAD STUDY
% cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
% if ~isempty(SUB_GROUP_FNAME)
%     spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
%     plot_store_dir = [cluster_dir filesep 'plots_out' filesep SUB_GROUP_FNAME];
% else
%     spec_data_dir = [cluster_dir filesep 'spec_data'];
%     plot_store_dir = [cluster_dir filesep 'plots_out'];
% end
% if ~exist(spec_data_dir,'dir')
%     error('spec_data dir does not exist');
% end
% if ~exist(plot_store_dir,'dir')
%     mkdir(plot_store_dir);
% end
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',spec_data_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',spec_data_dir);
% end
% cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
% STUDY.cluster = cl_struct;
% [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
% CLUSTER_PICKS = main_cl_inds(2:end);
%% NEW DIPOLE IMPLEMENTATION
HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
if ~ispc
    HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
else
    HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
end
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
hires_mesh = [fpath filesep fname '_dipplotvol.mat'];
hires_mri = [fpath filesep fname '_dipplotmri.mat'];
mri = load(hires_mri);
mri = mri.mri;
vol = hires_mesh;
%## default mri & vol
%{
tmp = strsplit(path,';');
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
mri = MNI_MRI;
mri = load(mri);
mri = mri.mri;
vol = MNI_VOL;
%}
%-
%     transform = [1,0,0,-99;...
%                      0,1,0,-135;...
%                      0,0,1,-73;...
%                      0,0,0,1];
save_dir = cluster_dir;
DIPPLOT_STRUCT = struct('rvrange',[0,30],... % this is a value from 0 to 100 (e.g., rv = 0.15 is 15)
        'summary','off',...
        'mri',mri,...
        'coordformat','MNI',...
        'transform',[],...
        'image','mri',...
        'plot','on',...
        'color',{{[0,0,1]}},...
        'view',[1,1,1],...
        'mesh','off',...
        'meshdata',vol,...
        'axistight','on',... % display the closest MRI slice to distribution
        'gui','off',...
        'num','off',...
        'cornermri','on',...
        'drawedges','off',...
        'projimg','off',...
        'projlines','off',...
        'projwidth',1,...
        'projcol',{{[0,0,1]}},...
        'dipolesize',30,...
        'dipolelength',0,...
        'pointout','off',...
        'sphere',1,...
        'spheres','off',...
        'normlen','off',...
        'dipnames',{{}},...
        'holdon','on',...
        'camera','auto',...
        'density','off');
%## ALL DIPOLES FOR CLUSTERS
[fig] = eeglab_dipplot(STUDY,ALLEEG,CLUSTER_PICKS,...
    'PLOT_TYPE','all_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
pause(2);
 % camzoom(1.2^2);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.tiff')],'Resolution',1000);
saveas(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
view([45,0,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_coronal.tiff')],'Resolution',1000);
view([0,-45,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.tiff')],'Resolution',1000);
drawnow;
close(fig);
%## AVERAGE DIPOLE FOR CLUSTERS
[fig] = eeglab_dipplot(STUDY,ALLEEG,CLUSTER_PICKS,...
    'PLOT_TYPE','average_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
%     zoom(1.05)
pause(2);
% camzoom(1.2^2);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.tiff')],'Resolution',1000);
saveas(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.fig')]);
view([45,0,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.tiff')],'Resolution',1000);
view([0,-45,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.tiff')],'Resolution',1000);
pause(2);
close(fig);
%##
for i = 1:length(CLUSTER_PICKS)
%         cluster_i = CLUSTER_PICKS(i);
%         [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
%             'PLOT_TYPE','all_nogroup',...
%             'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
%         pause(2);
% %         camzoom(1.2^2);
%         exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cluster_i)],'Resolution',1000);
%         view([45,0,0])
%         exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cluster_i)],'Resolution',1000);
%         view([0,-45,0])
%         exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cluster_i)],'Resolution',1000);
%         pause(2);
%         close(fig);
    cluster_i = CLUSTER_PICKS(i);
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_group',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    % GROUP_MARKERS = {'.','diamond','^','pentagram','hexagram'};
    % for ii = 1:length(fig.Children(end).Children)
    %     fig.Children(end).Children(ii+3).Marker = GROUP_MARKERS{gg};     
    % end
    pause(2);
%         camzoom(1.2^2);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cluster_i)],'Resolution',1000);
    view([45,0,0])
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cluster_i)],'Resolution',1000);
    view([0,-45,0])
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cluster_i)],'Resolution',1000);
    pause(2);
    close(fig);
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}