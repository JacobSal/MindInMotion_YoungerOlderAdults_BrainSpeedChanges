%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_terrain_clin/run_spca_d_tw_plots_clim.sh

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
    STUDY_DIR = fileparts(SCRIPT_DIR);
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
BOOT_NITERS = 2000;
BOOT_ALPHA = 0.05;
BOOT_CLUST_THRESH = 1000;
%- datset name
DATA_SET = 'MIM_dataset';
% study_dir_name = '04162024_MIM_OAN57_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
spca_dir_name = '03232024_spca_analysis_OA';
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
groups = {'YA','HOA','FOA'};
SUB_GROUP_FNAME = 'all_spec';
SUB_GROUP_GFNAME = 'group_spec';
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
save_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- design parameters
% STUDY_DESI_PARAMS = {{'subjselect',{},...
%             'variable2','cond','values2',{'flat','low','med','high'},...
%             'variable1','group','values1',{'H2000''s','H3000''s'}},...
%             {'subjselect',{},...
%             'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
%             'variable1','group','values1',{'H2000''s','H3000''s'}}};
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'}}};
%- spca dir
spca_dir_fpath = [studies_fpath filesep spca_dir_name];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% ===================================================================== %%
%## LOAD STUDY
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
gspec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_GFNAME];
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
%* assign studies
STUDY.design = STUDY_DESI_PARAMS;
%## CLUSTER_INFO
[comps_out,main_cl_inds,outlier_cl_inds,valid_cluster,~,nonzero_cluster] = eeglab_get_cluster_comps(STUDY);
main_cl_inds = main_cl_inds(2:end);
CLUSTER_PICKS = nonzero_cluster; %valid_cluster; %main_cl_inds(2:end);
fprintf('Clusters with more than 50%% of subjects:'); fprintf('%i,',valid_cluster(1:end-1)); fprintf('%i',valid_cluster(end)); fprintf('\n');
fprintf('Main cluster numbers:'); fprintf('%i,',main_cl_inds(1:end-1)); fprintf('%i',main_cl_inds(end)); fprintf('\n');
%% ===================================================================== %%
%## PARAMS
ATLAS_PATH = [PATHS.submods_dir,...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
atlas_i = 1;
% [STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
atlas_name_store = cell(length(main_cl_inds),1);
titles_store = cell(length(main_cl_inds),1);
for k_i = 1:length(main_cl_inds)
    cl_i = main_cl_inds(k_i);
    %## ANATOMY
    STUDY.cluster(cl_i).centroid.dipole.posxyz = mean(STUDY.cluster(cl_i).all_diplocs);
    dip1 = STUDY.cluster(cl_i).all_diplocs;
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
    atlas_name_store{cl_i} = atlas_name;
    titles_store{cl_i} = sprintf('CL%i: %s\n[%0.1f,%0.1f,%0.1f]\n',cl_i,atlas_name,STUDY.cluster(cl_i).centroid.dipole.posxyz);
end
%% ===================================================================== %%
% caption1 = ['\\textbf{Figure %i) %s area''s topographies and power spectral density plots.} (\\textbf{A})',...
%     ' Presented are the average scalp topography and the dipole locations for each component in',...
%     ' the cluster. Dipole locations are overlayed on the Montreal Neurological Institute template.',...
%     ' Average PSDs for the cluster across (\\textbf{B}) speeds and (\\textbf{C}) terrain difficulties.',...
%     ' On the left is the original PSD, and on the right is the flattened PSD. Shaded gray areas are,',...
%     ' frequencies where there is a significant differences between conditions.',...
%     ' Shaded color areas (i.e., green, yellow, orange, and red) are the standard error of the PSDs',...
%     ' across components in the cluster. Vertical dashed lines mark the frequency bands of interest:',...
%     ' theta (4-8Hz), alpha (8-13Hz), and beta (13-30Hz). Average power within each frequency band of',...
%     ' interest for each (\\textbf{D}) speed and (\\textbf{E}) terrain condition. (\\textbf{D}) Slope and',...
%     ' correlation-coefficient (\\textbf{R^2}) of the general linear model is presented if it is significant',...
%     ' across speeds. (textbf{E}) Black horizontal bars represent that significance level of the Post-Hoc test',...
%     ' comparing each terrain to the flat condition (*p<0.05, **p<0.01, ***p<0.001).'];
% caption2 = ['\\textbf{Figure %i) %s gait-phase modulation event-related spectral perturbations.}',...
%     ' x-axes: time of gait cycle with event markers for right heel strike (RHS),',...
%     ' left toe off (LTO), left heel strike (LHS), and right toe off (RTO). y-axes:',...
%     ' frequency in Hertz where bands of interest are denoted as theta (4-8Hz), alpha (8-13Hz),',...
%     ' beta (13-30Hz), and low-gamma (30-50Hz). Colors represent significant increases',...
%     ' (red, synchronization) and decreases (blue, desynchronization) in spectral power (dB)',...
%     ' from the mean spectrum for all gait cycles. Cluster-based bootstrapping statistics are included',...
%     ' in the far-right panel (p<0.05). Additionally, each condition is bootstrapped across time and alpha masked (i.e., more transparent)',...
%     ' if outside the 95% confidence interval of the surrogate distribution. Stats are corrected for multiple comparisons using the false discovery ratio.'];
% caption3 = ['\\textbf{Figure %i) %s common baselined event-related spectral perturbations.} (\\textbf{A}),',...
%     ' and their relative changes from baseline conditions, 1.0m/s for speed and flat for terrain (\\textbf{B}).',...
%     ' x-axes: time of gait cycle with event markers for right heel strike (RHS), left toe off (LTO),',...
%     ' left heel strike (LHS), and right toe off (RTO). y-axes: frequency in Hertz. Shown are frequency',...
%     ' bands of theta (4-8Hz), alpha (8-13Hz), beta (13-30Hz), and low-gamma (30-50Hz).',...
%     ' (\\textbf{A}) Colors represent significant increases (reds, synchronization) and decreases',...
%     ' (blues, desynchronization) in spectral power (dB) from the mean spectrum for all conditions.',...
%     ' (\\textbf{B}) Each condition.s ERSP is significance masked (p<0.05) using nonparametric bootstrapping',...
%     ' corrected for multiple comparisons using the false discovery ratio.'];
caption1 = ['Figure %i) %s area''s topographies and power spectral density plots. (A)',...
    ' Presented are the average scalp topography and the dipole locations for each component in',...
    ' the cluster. Dipole locations are overlayed on the Montreal Neurological Institute template.',...
    ' Average PSDs for the cluster across (B) speeds and (C) terrain difficulties.',...
    ' On the left is the original PSD, and on the right is the flattened PSD. Shaded gray areas are,',...
    ' frequencies where there is a significant differences between conditions.',...
    ' Shaded color areas (i.e., green, yellow, orange, and red) are the standard error of the PSDs',...
    ' across components in the cluster. Vertical dashed lines mark the frequency bands of interest:',...
    ' theta (4-8Hz), alpha (8-13Hz), and beta (13-30Hz). Average power within each frequency band of',...
    ' interest for each (D) speed and (E) terrain condition. (D) Slope and',...
    ' correlation-coefficient (R^2) of the general linear model is presented if it is significant',...
    ' across speeds. (E) Black horizontal bars represent that significance level of the Post-Hoc test',...
    ' comparing each terrain to the flat condition (*p<0.05, **p<0.01, ***p<0.001).'];
caption2 = ['Figure %i) %s gait-phase modulation event-related spectral perturbations.',...
    ' x-axes: time of gait cycle with event markers for right heel strike (RHS),',...
    ' left toe off (LTO), left heel strike (LHS), and right toe off (RTO). y-axes:',...
    ' frequency in Hertz where bands of interest are denoted as theta (4-8Hz), alpha (8-13Hz),',...
    ' beta (13-30Hz), and low-gamma (30-50Hz). Colors represent significant increases',...
    ' (red, synchronization) and decreases (blue, desynchronization) in spectral power (dB)',...
    ' from the mean spectrum for all gait cycles. Cluster-based bootstrapping statistics are included',...
    ' in the far-right panel (p<0.05). Additionally, each condition is bootstrapped across time and alpha masked (i.e., more transparent)',...
    ' if outside the 95% confidence interval of the surrogate distribution. Stats are corrected for multiple comparisons using the false discovery ratio.'];
caption3 = ['Figure %i) %s common baselined event-related spectral perturbations. (A),',...
    ' and their relative changes from baseline conditions, 1.0m/s for speed and flat for terrain (B).',...
    ' x-axes: time of gait cycle with event markers for right heel strike (RHS), left toe off (LTO),',...
    ' left heel strike (LHS), and right toe off (RTO). y-axes: frequency in Hertz. Shown are frequency',...
    ' bands of theta (4-8Hz), alpha (8-13Hz), beta (13-30Hz), and low-gamma (30-50Hz).',...
    ' (A) Colors represent significant increases (reds, synchronization) and decreases',...
    ' (blues, desynchronization) in spectral power (dB) from the mean spectrum for all conditions.',...
    ' (B) Each condition.s ERSP is significance masked (p<0.05) using nonparametric bootstrapping',...
    ' corrected for multiple comparisons using the false discovery ratio.'];
%%
%## Create PowerPoint
%- note: this essentially uses VBA functionality and many of the properites
%and methods overlap. (see.
%https://learn.microsoft.com/en-us/office/vba/api/);
%- start powerpoint
ppt = actxserver('powerpoint.application');
ppt.Visible = 1;
ppt.Presentation.invoke;
p = ppt.Presentation.Add;
layout = ppt.ActivePresentation.SlideMaster.CustomLayouts.Item(7);
% p.ApplyTheme('C:\Users\jsalminen\Desktop\figure_gen_teplate.potx')
%- note
% SlideHeight = 540 (7.5in)
% SlideWidth = 960 (13.3333in)
FONT_SIZE = 10;
FONT_NAME = 'Arial';
SLIDE_W = 6.5;
SLIDE_H = 9;
COVNERT_PPI_DPI = (540/7.5); % 72dpi
% scale_h = (SLIDE_H/COVNERT_PPI_DPI);
% scale_w = (SLIDE_W/COVNERT_PPI_DPI);
p.PageSetup.SlideWidth = COVNERT_PPI_DPI*SLIDE_W;
p.PageSetup.SlideHeight = COVNERT_PPI_DPI*SLIDE_H;
p.Slides.AddSlide(1,layout)
slides = p.Slides;
figure_cnt = 3;
for k_i = flip(main_cl_inds,2) %length(cl_names):-1:1 %1:length(cl_names)
    %%
    atlas_name = titles_store{k_i};
    spca_fpath = [cluster_dir filesep 'spca'];
    spec_fpath = [spec_data_dir filesep 'psd_calcs'];
    %## (SLIDE 8) Cluster Level GROUP GPMs (not bootstrapped?)
    % IM_DPI = 300;
    % im_scale = COVNERT_PPI_DPI/IM_DPI;
    % TOP_DIST = 50;
    % LEFT_DIST = 50;
    % des_i = 2;
    % %-
    % newSlide = slides.AddSlide(1,layout);
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    % Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    % %-
    % vertical_move = 0;
    % for g_i = 1:length(groups)
    %     try
    %         im_f = [spca_fpath filesep sprintf('cl%i_des%i_group%s_spcadiff_allersp_com.tiff',k_i,des_i,groups{g_i})];
    %         tmp = imread(im_f);
    %         newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    %         vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    %     catch e
    %         fprintf('%s\n',getReport(e));
    %     end
    % end
    % % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+vertical_move+im_scale*125,SLIDE_W*COVNERT_PPI_DPI,70);
    % % Title1.TextFrame.TextRange.Text = sprintf(caption2,figure_cnt-1,atlas_name_store{k_i});
    % % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    % 
    % %## (SLIDE 7) Cluster Level GROUP GPMs (not bootstrapped?)
    % IM_DPI = 300;
    % im_scale = COVNERT_PPI_DPI/IM_DPI;
    % TOP_DIST = 50;
    % LEFT_DIST = 50;
    % des_i = 1;
    % %-
    % newSlide = slides.AddSlide(1,layout);
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    % Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    % %-
    % vertical_move = 0;
    % for g_i = 1:length(groups)
    %     try
    %         im_f = [spca_fpath filesep sprintf('cl%i_des%i_group%s_spcadiff_allersp_com.tiff',k_i,des_i,groups{g_i})];
    %         tmp = imread(im_f);
    %         newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    %         vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    %     catch e
    %         fprintf('%s\n',getReport(e));
    %     end
    % end
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+vertical_move+im_scale*125,SLIDE_W*COVNERT_PPI_DPI,70);
    % Title1.TextFrame.TextRange.Text = sprintf(caption2,figure_cnt-1,atlas_name_store{k_i});
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;

    %## (SLIDE 7) Cluster Level GROUP GPMs (bootstrapped?)
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 50;
    LEFT_DIST = 50;
    des_i = 1;
    %-
    newSlide = slides.AddSlide(1,layout);
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %-
    vertical_move = 0;
    try
        im_f = [spca_fpath filesep sprintf('cl%i_des%i_group_bootstraps_ersp_sb.tiff',k_i,des_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
        vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    catch e
        fprintf('%s\n',getReport(e));
    end

    %## (SLIDE 7) Cluster Level GROUP GPMs (bootstrapped?)
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 50;
    LEFT_DIST = 50;
    des_i = 2;
    %-
    newSlide = slides.AddSlide(1,layout);
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %-
    vertical_move = 0;
    try
        im_f = [spca_fpath filesep sprintf('cl%i_des%i_group_bootstraps_ersp_sb.tiff',k_i,des_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
        vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    catch e
        fprintf('%s\n',getReport(e));
    end

    %## (SLIDE 6) Cluster Level GROUP ERSP COM
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 50;
    LEFT_DIST = 5;
    des_i = 2;
    %-
    newSlide = slides.AddSlide(1,layout);
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %-
    vertical_move = 0;
    try
        im_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_groupersp_com.tiff',k_i,des_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
        vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    catch e
        fprintf('%s\n',getReport(e));
    end
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+vertical_move+im_scale*125,SLIDE_W*COVNERT_PPI_DPI,70);
    % Title1.TextFrame.TextRange.Text = sprintf(caption2,figure_cnt-1,atlas_name_store{k_i});
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;

    %## (SLIDE 5) Cluster Level GROUP ERSP COM
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 50;
    LEFT_DIST = 5;
    des_i = 1;
    %-
    newSlide = slides.AddSlide(1,layout);
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %-
    vertical_move = 0;
    try
        im_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_groupersp_com.tiff',k_i,des_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
        vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125; 
    catch e
        fprintf('%s\n',getReport(e));
    end
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+vertical_move+im_scale*125,SLIDE_W*COVNERT_PPI_DPI,70);
    % Title1.TextFrame.TextRange.Text = sprintf(caption2,figure_cnt-1,atlas_name_store{k_i});
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %## (SLIDE 4) Common Baselined ERSPs & Differences
    % IM_DPI = 300;
    % im_scale = COVNERT_PPI_DPI/IM_DPI;
    % TOP_DIST = 50;
    % LEFT_DIST = 0;
    %## ERSP COMMON BASELINED
    % newSlide = slides.AddSlide(1,layout);
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    % % Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    % % unicode =double(sprintf(caption3,figure_cnt,atlas_name_store{k_i}));
    % % for i = 1:length(unicode)
    % %     Title1.TextFrame.TextRange.InsertSymbol('Arial',)
    % % end
    % caption_height = 0;
    % vertical_move = 0;
    % for des_i = 1:length(STUDY.design)
    %     try
    %         im_f = [spca_fpath filesep sprintf('cl%i_des%i_spca_ersp_com.jpg',k_i,des_i)];
    %         tmp = imread(im_f);
    %         newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    %         vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125;
    %     catch e
    %         fprintf('%s\n',getReport(e));
    %     end
    % end
    % caption_height = vertical_move;
    %## DIFFERENCE ERSPS
    % IM_DPI = 300;
    % im_scale = COVNERT_PPI_DPI/IM_DPI;
    % vertical_move = vertical_move + im_scale*125;
    % for des_i = 1:length(STUDY.design)
    %     try
    %         im_f = [spca_fpath filesep sprintf('cl%i_des%i_spcadiff_allersp_com.jpg',k_i,des_i)];
    %         tmp = imread(im_f);
    %         newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    %         vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125;
    %     catch e
    %         fprintf('%s\n',getReport(e));
    %     end
    % end
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',size(tmp,2)*im_scale,TOP_DIST+caption_height+im_scale*125,COVNERT_PPI_DPI*6.5-size(tmp,2)*im_scale,70);
    % Title1.TextFrame.TextRange.Text = sprintf(caption3,figure_cnt,atlas_name_store{k_i});
    % % newSlide.TextRange.MathZones
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    %## (SLIDE 3) Cluster Level GPMs (bootstrapped?)
    % IM_DPI = 300;
    % im_scale = COVNERT_PPI_DPI/IM_DPI;
    % TOP_DIST = 50;
    % LEFT_DIST = 50;
    % %-
    % newSlide = slides.AddSlide(1,layout);
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,400,70);
    % Title1.TextFrame.TextRange.Text = sprintf('(N=%i) %s',length(STUDY.cluster(k_i).sets),atlas_name);
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    % %-
    % vertical_move = 0;
    % for des_i = 1:length(STUDY.design)
    %     try
    %         im_f = [spca_fpath filesep sprintf('cl%i_des%i_bootstraps_ersp_sb.jpg',k_i,des_i)];
    %         tmp = imread(im_f);
    %         newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST+vertical_move,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    %         vertical_move = vertical_move + size(tmp,1)*im_scale-im_scale*125;  
    %     catch e
    %         fprintf('%s\n',getReport(e));
    %     end
    % end
    % Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+vertical_move+im_scale*125,SLIDE_W*COVNERT_PPI_DPI,70);
    % Title1.TextFrame.TextRange.Text = sprintf(caption2,figure_cnt-1,atlas_name_store{k_i});
    % Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    % Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    
    %## (SLIDE 2) GROUP VIOLIN STATS
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 30;
    LEFT_DIST = 0;
    newSlide = slides.AddSlide(1,layout);
    
    try
        im_f = [gspec_data_dir filesep 'psd_calcs' filesep sprintf('Group_Violins_cl%i.tiff',k_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    catch e
        fprintf('%s\n',getReport(e));
    end
    figure_cnt = figure_cnt + 3;
    %## (SLIDE 1) GROUP PSD STATS
    IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 30;
    LEFT_DIST = 0;
    newSlide = slides.AddSlide(1,layout);
    %- title slide
    try
        im_f = [gspec_data_dir filesep 'psd_calcs' filesep sprintf('Group_PSDs_cl%i.tiff',k_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    catch e
        fprintf('%s\n',getReport(e));
    end
    figure_cnt = figure_cnt + 3;
    %## (SLIDE TITLE) TOPO, DIPOLE, AND PSD STATS
    IM_DPI = 300;
    CATCH_IM_DPI = 300;
    im_scale = COVNERT_PPI_DPI/IM_DPI;
    TOP_DIST = 30;
    LEFT_DIST = 0;
    newSlide = slides.AddSlide(1,layout);
    %- title slide
    subj_inds = STUDY.cluster(k_i).sets;
    subj_char = {STUDY.datasetinfo(subj_inds).subject};
    subj_char = strjoin(subj_char,', ');
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,0,475,70);
    Title1.TextFrame.TextRange.Text = sprintf('%s\nSubjects: %s',atlas_name,subj_char);
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    try
        % im_f = [spec_fpath filesep sprintf('cl%i_topo-dips-psd-violins.tiff',k_i)];
        im_f = [gspec_data_dir filesep 'psd_calcs' filesep sprintf('Group_speed_violin_psds_cl%i.tiff',k_i)];
        tmp = imread(im_f);
        newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
    catch e
        fprintf('%s\n',getReport(e));
        try
            im_scale = COVNERT_PPI_DPI/CATCH_IM_DPI;
            im_f = [gspec_data_dir filesep 'psd_calcs' filesep sprintf('TOPO_DIP_cl%i.tiff',k_i)];
            tmp = imread(im_f);
            newSlide.Shapes.AddPicture(im_f, 'msoFalse', 'msoTrue',LEFT_DIST,TOP_DIST,size(tmp,2)*im_scale,size(tmp,1)*im_scale); %Left, top, width, height
        catch e
            fprintf('%s\n',getReport(e));
        end
    end
    Title1 = newSlide.Shapes.AddTextbox('msoTextOrientationHorizontal',0,TOP_DIST+size(tmp,1)*im_scale,SLIDE_W*COVNERT_PPI_DPI,70);
    Title1.TextFrame.TextRange.Text = sprintf(caption1,figure_cnt-2,atlas_name_store{k_i});
    Title1.TextFrame.TextRange.Font.Name = FONT_NAME;
    Title1.TextFrame.TextRange.Font.Size = FONT_SIZE;
    figure_cnt = figure_cnt + 3;
end
ppt.ActivePresentation.SaveAs([cluster_dir filesep sprintf('spcasimple_report.pptx')]);
% Close Powerpoint and delete the object
ppt.ActivePresentation.Close;
% ppt.Quit;
% ppt.delete;
fprintf('\nResults powerpoint saved here:\n%s\n',cluster_dir);
