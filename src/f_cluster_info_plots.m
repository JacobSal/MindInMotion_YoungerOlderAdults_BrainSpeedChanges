%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_f_cluster_info_plots.sh

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
ADD_ALL_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = STUDY_DIR;
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
    STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PATHS) ============================================================= %%
%## DATASET
DATA_SET = 'MIM_dataset';
%## STUDY INFO
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
% STUDY_DNAME = '01192025_mim_yaoa_nopowpow_crit_speed';
% STUDY_FNAME = 'epoch_study';
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
%## CLUSTER
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_dir filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- save dir
save_dir = [cluster_k_dir filesep 'topo_dip_inf'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ================================================================== %%
%## LOAD STUDY
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
% end
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
[comps_out,main_cl_inds,outlier_cl_inds,valid_clusters] = eeglab_get_cluster_comps(STUDY);
cluster_inds = main_cl_inds(1:end);
% cluster_inds = [3,4,6,7,8,9,10,13]; % 01192025_mim_yaoa_nopowpow_crit_speed, rb3
%- save dir
save_dir = [save_dir filesep 'figure_gen']; %sprintf('cl%i-cl%i',min(cluster_inds),max(cluster_inds))];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD TOPO & DIP INFO ================================================ %%
% dipfit_structs = par_load(STUDY.filepath,'dipfit_structs.mat');
% topo_cells = par_load(STUDY.filepath,'topo_cells.mat');
%--
% CALC_STRUCT = struct('cluster_inds',(2:length(STUDY.cluster)),...
%     'save_inf',true,...
%     'recalculate',true);
% [STUDY,dipfit_structs,topo_cells] = eeglab_get_topodip(STUDY,...
%     'CALC_STRUCT',CALC_STRUCT,...
%     'ALLEEG',ALLEEG);
%% (ANATOMY) =========================================================== %%
% addpath([PATHS.submods_dir filesep 'AAL3']);
ANATOMY_STRUCT = struct('atlas_fpath',{{[PATHS.submods_dir filesep 'AAL3' filesep 'AAL3v1.nii'],...
    [PATHS.submods_dir filesep 'AAL3' filesep 'ROI_MNI_V7.nii']}},...
    'group_chars',{unique({STUDY.datasetinfo.group})},...
    'cluster_inds',(cluster_inds),...
    'anatomy_calcs',{{'all aggregate','all centroid'}},... % ('all calcs','group centroid','all centroid','group aggregate','all aggregate')
    'save_dir',cluster_k_dir,...
    'save_inf',true);
[STUDY,anat_struct,~,~,txt_out] = eeglab_get_anatomy(STUDY,...
    'ANATOMY_STRUCT',ANATOMY_STRUCT,...
    'ALLEEG',ALLEEG);
%% TOPO ================================================================ %%
groups = unique({STUDY.datasetinfo.group});
AX_HORIZ_SHIFT = 0.4;
AX_VERT_SHIFT = 0.05;
DIP_IM_DPI = 1200;
AX_INIT_HORIZ_TOPO = 0.085;
AX_INIT_VERT_TOPO = 0.75;
FIGURE_POSITION =[1,1,6.5,9];
PG_SIZE = [6.5,9];
FONT_NAME = 'Arial';
TOPO_FONTSIZE = 7;
HZ_DIM = 3;
g_chars = unique({STUDY.datasetinfo.group});
g_chars_topo = {'Young Adults','Older High Mobility Adults','Older Low Mobility Adults'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
% AX_H = 0.25;
% AX_W = 0.225;
% IM_RESIZE = 0.7;
IM_RESIZE = 0.225;
horiz_shift = 0;
vert_shift = 0;
hz = 1;
%## TOPO PLOT
fig = figure('color','white');
set(fig,'Units','inches','Position',FIGURE_POSITION)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
set(gca,AXES_DEFAULT_PROPS{:});
hold on;
for k_i = 1:length(cluster_inds)
    cl_i = double(string(cluster_inds(k_i)));
    %## TOPO PLOT
    [~,h] = std_topoplot_CL(STUDY,cl_i,'together');
    set(h,'color','w')
    set(h,'Units','normalized');    
    ax = get(h,'CurrentAxes');
    colormap(h,linspecer);   
    pos = get(ax,'Position');
    opos = get(ax,'OuterPosition');
    ac = get(ax,'Children');
    ax1 = axes('Parent',fig,...
        'DataAspectRatio',get(ax,'DataAspectRatio'));
    copyobj(ac,ax1);
    set(ax1,AXES_DEFAULT_PROPS{:});
    colormap(ax1,linspecer)
    g_counts = cell(length(groups),1);
    for g_i = 1:length(groups)
        g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
        if length(g_chars_topo{g_i}) == 1 || ischar(g_chars_topo{g_i})
            g_counts{g_i} =sprintf('%s N=%i',g_chars_topo{g_i},sum(g_inds));
        else
            g_counts{g_i} =sprintf('%s\n%s N=%i',g_chars_topo{g_i}{1},g_chars_topo{g_i}{2},sum(g_inds));
        end
    end
    ax1.Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    ax1.Title.Interpreter = 'none';
    ax1.FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    ax1.FontName = FONT_NAME;
    ax1.FontWeight = 'bold';
    ax1.OuterPosition = [0,0,1,1];
    ax1.Units = 'Normalized';
    pp = get(ax1,'Position');    
    ax1.Position = [AX_INIT_HORIZ_TOPO+horiz_shift,AX_INIT_VERT_TOPO+vert_shift,pp(3)*IM_RESIZE,pp(4)*IM_RESIZE];  %[left bottom width height]
    close(h);
    % std_topoplot_CL(STUDY,cl_i,'together');
    % colormap(linspecer); %colormap_ersp)
    % fig_i = get(groot,'CurrentFigure');
    % set(fig_i,'color','w')
    % g_counts = cell(length(g_chars),1);
    % for group_i = 1:length(g_chars)
    %     g_inds = cellfun(@(x) strcmp(x,g_chars{group_i}),{STUDY.datasetinfo(STUDY.cluster(cl_i).sets).group});
    %     if length(g_chars_topo{group_i}) == 1 || ischar(g_chars_topo{group_i})
    %         g_counts{group_i} =sprintf('%s N=%i',g_chars_topo{group_i},sum(g_inds));
    %     else
    %         g_counts{group_i} =sprintf('%s\n%s N=%i',g_chars_topo{group_i}{1},g_chars_topo{group_i}{2},sum(g_inds));
    %     end
    % end
    % obj_i = findobj(fig_i,'type','Axes');
    % obj_i(1).Title.String = [{sprintf('CL%i',cl_i)};g_counts]; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
    % obj_i(1).Title.Interpreter = 'none';
    % obj_i(1).FontSize = TOPO_FONTSIZE; %PLOT_STRUCT.font_size;
    % obj_i(1).FontName = FONT_NAME;
    % obj_i(1).FontWeight = 'bold';
    % obj_i(1).OuterPosition = [0,0,1,1];
    % obj_i(1).Units = 'Normalized';
    % obj_i(1).Position = [AX_INIT_HORIZ_TOPO+horiz_shift,AX_INIT_VERT_TOPO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE];  %[left bottom width height]
    %##
    if hz < HZ_DIM
        horiz_shift = horiz_shift + pp(3)*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
    else
        vert_shift = vert_shift - (pp(4)*IM_RESIZE + AX_VERT_SHIFT*IM_RESIZE);
        horiz_shift = 0;
        hz = 0;
    end
    hz = hz + 1;
end
hold off;
drawnow;
exportgraphics(fig,[save_dir filesep 'topo_plots.pdf'],'ContentType','vector')
exportgraphics(fig,[save_dir filesep 'topo_plots.tif'],'ContentType','image','Resolution',600)
savefig(fig,[save_dir filesep 'topo_plots.fig'])
% close(fig);
%% NEW DIPOLE IMPLEMENTATION
% HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
% if ~ispc
%     HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
% else
%     HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
% end
HIRES_TEMPLATE = [PATHS.data_dir filesep '_resources' filesep 'mni_icbm152_nlin_sym_09a' filesep 'mni_icbm152_t1_tal_nlin_sym_09a.nii'];
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
        'axistight','off',... % display the closest MRI slice to distribution
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
DIPPLOT_STRUCT.axistight = 'on';
%## ALL DIPOLES FOR CLUSTERS
[fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_inds,...
    'PLOT_TYPE','all_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
pause(2);
% camzoom(1.2^2);
% exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
saveas(fig,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
view([45,0,0])
% exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_coronal.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_coronal.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
view([0,-45,0])
% exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
drawnow;
close(fig);
%## AVERAGE DIPOLE FOR CLUSTERS
[fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_inds,...
    'PLOT_TYPE','average_nogroup',...
    'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
%     zoom(1.05)
pause(2);
% camzoom(1.2^2);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
saveas(fig,[save_dir filesep sprintf('dipplot_avgdipspc_top.fig')]);
view([45,0,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
view([0,-45,0])
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.tiff')],'Resolution',DIP_IM_DPI);
exportgraphics(fig,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.pdf')],'ContentType','vector','Resolution',DIP_IM_DPI);
pause(2);
close(fig);
%##

for i = 1:length(cluster_inds)
    cluster_i = cluster_inds(i);
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_nogroup',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    pause(2);
%         camzoom(1.2^2);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    saveas(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_top.fig',cluster_i)]);
    view([45,0,0])
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_coronal.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_coronal.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    view([0,-45,0])
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_sagittal.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_nogroup_sagittal.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    pause(2);
    close(fig);
    cluster_i = cluster_inds(i);
    [fig] = eeglab_dipplot(STUDY,ALLEEG,cluster_i,...
        'PLOT_TYPE','all_group',...
        'DIPPLOT_STRUCT',DIPPLOT_STRUCT);
    % GROUP_MARKERS = {'.','diamond','^','pentagram','hexagram'};
    % for ii = 1:length(fig.Children(end).Children)
    %     fig.Children(end).Children(ii+3).Marker = GROUP_MARKERS{gg};     
    % end
    pause(2);
%         camzoom(1.1^2);
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    saveas(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.fig',cluster_i)]);
    view([45,0,0])
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    view([0,-45,0])
    % exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cluster_i)],'Resolution',DIP_IM_DPI);
    exportgraphics(fig,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.pdf',cluster_i)],'ContentType','vector','Resolution',DIP_IM_DPI);
    pause(2);
    close(fig);

    %##
    
end