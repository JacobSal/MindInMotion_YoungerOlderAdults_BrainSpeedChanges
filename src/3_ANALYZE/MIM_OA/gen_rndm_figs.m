%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_c_ersp_gen_expanded_v20.sh

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
%     SLURM_POOL_SIZE = 2;
%     pp = parcluster('local');
%     pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
end
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
DO_SUBJ_PLOTS = true;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
do_multivariate_data = 1;
evaluate_method = 'min_rv';
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
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
% SPEC_PARAMS = struct('freqrange',[1,200],...
%     'subject','',...
%     'specmode','psd',...
%     'freqfac',4,...
%     'logtrials','on',...
%     'comps','all',...
%     'plot_freqrange',[4,60]);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
%## Soft Define
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% ===================================================================== %%
%## ADMIN SET
%## GET EEGLAB PATH
tmp = strsplit(path,';');
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
%## ATLAS
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'brainnetome' filesep 'BNA_MPM_thr25_1.25mm.nii'],...
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat'],...
    [ATLAS_PATH filesep 'vtpm' filesep 'vtpm.mat'],...
    [ATLAS_PATH filesep 'yeo' filesep 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'],...
    [ATLAS_PATH filesep 'brainweb' filesep 'brainweb_discrete.mat']}; % also a discrete version of this
SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\dipole_1_scalp_0_ersp_0_spec_0';
%- convert SUB_DIR
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true,true};
CLUSTER_K_PICKS = [14,14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'subjrejs_minics6' filesep '14'],...
    [SUB_DIR filesep 'subjrejs_minics5' filesep '14']};
CLUSTER_FILES = {'cluster_update_14.mat','cluster_update_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'subjrejs_minics6'],...
    [SUB_DIR filesep 'subjrejs_minics5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
CLUSTER_CLIM_MATCH = [];
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};   
%%
%## assign a subject path
SUBJECT_ID = 'NH3128';
subject_fpath = [DATA_DIR filesep DATA_SET filesep SUBJECT_ID filesep 'MRI'];
vol = load([subject_fpath filesep 'vol.mat']);
vol = vol.vol;
%% (PLOT 3)
mesh = load([subject_fpath filesep 'mesh.mat']);
mesh = mesh.mesh;
TISSUES_TO_PLOT = [3,5,6]; %[3,5,6]; %[1,2,5]; %3,1,2]; %4; %[1,3,4]; %1:length(mesh.tissuelabel);
TISSUE_ALPHA = [0.2,0.3,0.3];
MATERIAL = [0.5,0.1,0.1,1,0.01];
VIDEO_ANGLE = 30;
mesh_cols = linspecer(length(mesh.tissuelabel)); %linspecer(length(TISSUES_TO_PLOT));
% 1: air
% 2: csf
% 3: gray
% 4: scalp
% 5: skull
% 6: white
figure;
tmp_mesh = mesh;
alpha_shift = 0.9;
SHIFT = 0.2;
hold on;
for i = 1:length(TISSUES_TO_PLOT)
    tmp_mesh.hex = mesh.hex(mesh.tissue==TISSUES_TO_PLOT(i),:);
    % tmp_mesh.tri = tri;
    % tmp_mesh = rmfield(tmp_mesh,'hex');
    ft_plot_mesh(tmp_mesh, 'surfaceonly','yes',...
        'facecolor', mesh_cols(TISSUES_TO_PLOT(i),:),...
        'edgecolor', 'none',...
        'vertexcolor', 'none',...
        'facealpha', TISSUE_ALPHA(i),...
        'edgealpha', TISSUE_ALPHA(i),...
        'material',MATERIAL,... %material(MATERIAL),...
        'maskstyle','opacity');
end
hold off;
%- save specific view
% view([90,90,90]);
% camzoom(1.2);
fig = get(groot,'CurrentFigure');
set(fig,'Color','w')
set(fig,'Units','inches','Position',[3 3 6 5])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
% exportgraphics(fig,[subject_fpath filesep sprintf('mesh_colormaped.jpg')],'Resolution',300);
%- make video
view([0,0,90]);
ViewZ = zeros(360*2,2);
ViewZ(:,1) = 0.5:0.5:360;
ViewZ(:,2) = repmat(VIDEO_ANGLE,size(ViewZ,1),1);
OptionZ = [];
OptionZ.FrameRate = 30;
OptionZ.Duration = 8;
OptionZ.Periodic = false;
CaptureFigVid(ViewZ,[subject_fpath filesep sprintf('mesh_colormaped.mp4')],OptionZ)

%% (PLOT)
% SUBJECT_ID = 'NH3128';
% subject_fpath = [DATA_DIR filesep DATA_SET filesep SUBJECT_ID filesep 'MRI'];
source_in = ft_read_mri([subject_fpath filesep  sprintf('%s_masks_contr.nii',SUBJECT_ID)]);
source_in.coordsys = 'acpc';
segmented = source_in;
segmented.white = source_in.anatomy == 1;
segmented.gray = source_in.anatomy == 2;
segmented.csf = (source_in.anatomy == 3 | source_in.anatomy == 8);% csf + ventricles
segmented.skull = source_in.anatomy == 4;
segmented.scalp = (source_in.anatomy == 5 | source_in.anatomy == 7);%skin and eye
segmented.air = source_in.anatomy == 6; 
segmented = rmfield(segmented,'anatomy');
seg_i_headreco = ft_datatype_segmentation(segmented,'segmentationstyle','indexed');
%-
% cfg = [];
% cfg.method = 'slice';
% cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
% cfg.anaparameter = 'anatomy';
% cfg.funcolormap  = linspecer(7);
% cfg.position = [];
% cfg.downsample = 1;
% cfg.colorbartext = 'slice';
% cfg.colorbar = 'yes';
% cfg.funcolormap = linspecer;
% cfg.renderer  = 'painter';
%-
% seg_cmap = [[1,1,1];mesh_cols(3,:);...
%     mesh_cols(6,:);mesh_cols(6,:);...
%     mesh_cols(6,:);mesh_cols(6,:);mesh_cols(3,:);];
% 1: air
% 2: csf
% 3: gray red
% 4: scalp 
% 5: skull blue
% 6: white green
cfg              = [];
cfg.method = 'ortho';
cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
cfg.anaparameter = 'anatomy';
cfg.funcolormap  = [[1,1,1];mesh_cols]; %linspecer(6); % distinct color per tissue
cfg.location     = 'center';
% cfg.atlas        = ft_read_atlas(ATLAS_FPATHS{2}); %seg_i_headreco;    % the segmentation can also be used as atlas
cfg.renderer     = 'painter';
cfg.crosshair = 'no';
hold on;
ft_sourceplot(cfg,seg_i_headreco);
hold off;
fig = get(groot,'CurrentFigure');
set(fig,'Color','w')
set(fig,'Units','inches','Position',[3 3 6 5])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
exportgraphics(fig,[subject_fpath filesep sprintf('segmentation_colormaps.jpg')],'Resolution',300);
    
% fig.Children.FontSize = 13;
% ax = gca; %fig_i.CurrentAxes;
% set(ax,'LineWidth',1)
% set(ax,'FontName','Arial','FontSize',12,'FontWeight','bold')
% set(ax,'OuterPosition',[0 0 1 1]);
% set(ax,'Position',[0.15 0.2 0.8 0.7]);  %[left bottom width height] This I what I added, You need to play with this
% %     set(ax,'XTick',sort(xticks));
% set(ax,'XTickLabel',trial_names_terrain);
%% ===================================================================== %%
%## MAKE EDITS TO DIPPLOTS FIGURE
DIP_SING_POS=[16 582 420 360];
tmp = load('-mat','M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\_save\10022023_cluster_slowwalkers_included\icrej_5\14\spec_data\temp_study_rejics5.study');
STUDY = tmp.STUDY;
[STUDY,ALLEEG] = pop_loadstudy('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\_save\10022023_cluster_slowwalkers_included\icrej_5\14\spec_data\temp_study_rejics5.study');
save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\_save\10022023_cluster_slowwalkers_included\icrej_5\14';
%- get inds
[~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
%- clusters to plot
CLUSTERS_IN_FIG = main_cl_inds(2:end);
CLUSTERS_TO_KEEP = valid_clusters;
%%
STUDY.etc.dipparams.centrline = 'off';
std_dipplot_CL(STUDY,ALLEEG,'clusters',CLUSTERS_TO_KEEP,'figure','off','mode','together_averaged_only','spheres','off','projlines','off');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',DIP_SING_POS,'color','w')
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_avgdipspc_top.pdf')],'ContentType','vector','Resolution',300);
saveas(fig_i,[save_dir filesep sprintf('new_dipplot_avgdipspc_top.fig')]);
view([45,0,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_avgdipspc_coronal.pdf')],'ContentType','vector','Resolution',300);
view([0,-45,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_avgdipspc_sagittal.pdf')],'ContentType','vector','Resolution',300);
%-
STUDY.etc.dipparams.centrline = 'off';
std_dipplot_CL(STUDY,ALLEEG,'clusters',CLUSTERS_TO_KEEP,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',DIP_SING_POS,'color','w')
camzoom(1);
% camzoom(1.2^2);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_top.pdf')],'ContentType','vector','Resolution',300);
saveas(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_top.fig')]);
view([45,0,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_coronal.pdf')],'ContentType','vector','Resolution',300);
view([0,-45,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_sagittal.pdf')],'ContentType','vector','Resolution',300);
%%
COLOR_OPTS = linspecer(length(CLUSTERS_TO_KEEP)+1);
colors = cell(1,size(COLOR_OPTS,1));
for i = 1:size(COLOR_OPTS,1)
    colors{i} = COLOR_OPTS(i,:);
end
for i = 1:length(CLUSTERS_TO_KEEP)
    cluster_i = CLUSTERS_TO_KEEP(i);
    STUDY.etc.dipparams.centrline = 'off';
    std_dipplot_CL(STUDY,ALLEEG,'clusters',cluster_i,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
    fig_i = get(groot,'CurrentFigure');
    set(fig_i,'position',DIP_SING_POS,'color','w')
    for j = 1:length(fig_i.Children(2).Children)
        try
            fig_i.Children(2).Children(j).Color = colors{i};
        catch
            fprintf('Can''t change the color of this child\n')
        end
    end
%     camzoom(1);
    % exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.jpg')],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('%inew_dipplot_alldipspc_top.pdf',cluster_i)],'ContentType','vector','Resolution',300);
%     saveas(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_top.fig')]);
    view([45,0,0])
    % exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.jpg')],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('%inew_dipplot_alldipspc_coronal.pdf',cluster_i)],'ContentType','vector','Resolution',300);
    view([0,-45,0])
    % exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.jpg')],'Resolution',300);
    exportgraphics(fig_i,[save_dir filesep sprintf('%inew_dipplot_alldipspc_sagittal.pdf',cluster_i)],'ContentType','vector','Resolution',300);
end

%% (TOPO)
if ~isfield(STUDY.cluster,'topo') 
    STUDY.cluster(1).topo = [];
end
for i = 1:length(CLUSTERS_TO_KEEP) % For each cluster requested
    clust_i = CLUSTERS_TO_KEEP(i);
    disp(clust_i)
    if isempty(STUDY.cluster(clust_i).topo)
        % Using this custom modified code to allow taking average within participant for each cluster
        STUDY = std_readtopoclust_CL(STUDY,ALLEEG,clust_i);
    end
end
TOPO_ALL_POS=[16 100 1240 920];
figure;
std_topoplot_CL(STUDY,CLUSTERS_TO_KEEP,'together');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',TOPO_ALL_POS,'color','w')
for c = 2:length(fig_i.Children)
    fig_i.Children(c).Title.Interpreter = 'none';
    fig_i.Children(c).FontSize = 13;
    fig_i.Children(c).FontName = 'Arial';
end
% saveas(fig_i,fullfile(save_dir,'Cluster_topo_avg.jpg'));
exportgraphics(fig_i,[save_dir filesep sprintf('new_cluster_topo_avg.pdf')],'ContentType','vector','Resolution',300);
saveas(fig_i,fullfile(save_dir,'new_cluster_topo_avg.fig'));
%{
% OFFSET = 4;
% del_rng_min = 1;
% del_rng_max = 1;
% inds_del = [];
% uiopen('M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\_save\10022023_cluster_slowwalkers_included\icrej_5\14\dipplot_alldipspc_top.fig');
% fig_i = gcf;
% for cluster_i = CLUSTERS_IN_FIG
%     del_rng_max = del_rng_max+length(STUDY.cluster(cluster_i).sets)+OFFSET;
%     if any(cluster_i == CLUSTERS_TO_KEEP)
%         fprintf('keeping %i\n',cluster_i)
%     else
%         fprintf('min: %i\n',del_rng_min);
%         fprintf('max: %i\n',del_rng_max);
%         fprintf('cluster: %i\n',cluster_i)
%         inds_del = [inds_del, (del_rng_min:del_rng_max)];
%     end
%     del_rng_min = del_rng_max+1;
% end
% delete(fig_i.Children.Children(inds_del))
% view([0,0,90]);
% view([90,0,0]);
% DELETE_INDS = 
%}