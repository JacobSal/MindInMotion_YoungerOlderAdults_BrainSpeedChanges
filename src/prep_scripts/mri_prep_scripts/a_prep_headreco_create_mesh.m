% Use the segmentation outcome from headreco
% Chang Liu - 2022-05-04
% Editted 2022-5-16 to create loop
% Editted by Yiru to create automatic transfer to M drive 
% 2022-11-1 create mesh for all MiM participants
% 2023-02-09 create mesh for all MiM participants - including older adults
%% Initialization
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
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
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% DEFINE SOURCE DIRECTORY & CD ======================================== %%
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '1_BATCH_PREP' filesep 'MIM_OA'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = true;
setWorkspace
%-
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('oa');
SUBJ_PICS = {{'H2012_FU','H2018_FU'},{'H3120','NH3129'}};
% GROUP_NAMES = {'H2000''s','H3000''s'};
% SUBJ_ITERS = {1:length(SUBJ_PICS{1}),1:length(SUBJ_PICS{2})};
%- 
MIM_RDRIVE = 'R:\Ferris-Lab\share\MindInMotion\Data\';
MIM_MDRIVE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset';
all_subjStr = [SUBJ_PICS{:}];
%-
COND_VALS = [1.65,0.33,0.33,0.01,0.126,2.5*10^(-14)];
%% Loop through all participants
for i = 1:length(all_subjStr)
    subjStr =  all_subjStr{i};
    source = fullfile(MIM_RDRIVE,subjStr,'MRI','Headmodel','skull_0.01'); %change to where the files are generated
    mkdir(source)
    Main_folder_headreco = fullfile(MIM_RDRIVE,subjStr,'\MRI\Segmentation\headreco\');    
    if ~exist(fullfile(source,'vol.mat'),'file') && exist(Main_folder_headreco, 'dir')
        cd(fullfile(Main_folder_headreco,['m2m_',subjStr]));
        gunzip([subjStr,'_masks_contr.nii.gz'])
        allMask = ft_read_mri([subjStr,'_masks_contr.nii']);
        allMask.coordsys = 'acpc';

        segmented = allMask;
        segmented.white = allMask.anatomy == 1;
        segmented.gray = allMask.anatomy == 2;
        segmented.csf = (allMask.anatomy == 3 | allMask.anatomy == 8);% csf + ventricles
        segmented.skull = allMask.anatomy == 4;
        segmented.scalp = (allMask.anatomy == 5 | allMask.anatomy == 7);%skin and eye
        segmented.air = allMask.anatomy == 6; 
        % segmented.eye = allMask.anatomy == 7;

        segmented = rmfield(segmented,'anatomy');
        seg_i_headreco = ft_datatype_segmentation(segmented,'segmentationstyle','indexed');

        cfg              = [];
        cfg.funparameter = 'tissue'; % They did an update on May.2021 in source code but not the tutorial
        cfg.anaparameter = 'anatomy';
        cfg.funcolormap  = turbo(7); % distinct color per tissue
        cfg.location     = 'center';
        cfg.atlas        = seg_i_headreco;    % the segmentation can also be used as atlas

        ft_sourceplot(cfg, seg_i_headreco);

        %% Test if fieldtrip can make headmodel from headreco segmentation --------------------------------------
        % Create mesh
        cfg        = [];
        cfg.shift  = 0.3;
        cfg.method = 'hexahedral';
        mesh = ft_prepare_mesh(cfg,segmented);

        % Create head model
        cfg        = [];
        cfg.method = 'simbio';
        cfg.conductivity = zeros(1,5);
        scale = 1;

        % shell conducantces are changed to match BEM model conduactances 
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'csf'))) = COND_VALS(1)*scale;
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'gray'))) = COND_VALS(2)*scale;
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'scalp'))) = COND_VALS(3)*scale;
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'skull'))) = COND_VALS(4)*scale;   
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'white'))) = COND_VALS(5)*scale;
        cfg.conductivity(find(strcmp(mesh.tissuelabel,'air'))) = COND_VALS(6)*scale;

        vol        = ft_prepare_headmodel(cfg, mesh);
        save(fullfile(source,'vol.mat'),'vol','-v7.3') 


        %% Align electrodes
        clear fid
        load(fullfile(MIM_RDRIVE,subjStr,'\MRI\Processed_fiducials','mri_acpc.mat'))
        load(fullfile(MIM_RDRIVE,subjStr,'\MRI\Processed_fiducials','ctf_fiducials.mat')); 
        fiducial_locs = ctf_fiducials;
        chanloc_scan_folder = fullfile(MIM_RDRIVE,subjStr,'\HeadScan\CustomElectrodeLocations.txt');
        chanlocs = readtable(chanloc_scan_folder);% Same output text file from getchalocs.
        chanlocs.Properties.VariableNames = {'labels','X','Y','Z'};
        elec.chanpos(:,1) = [chanlocs.X];
        elec.chanpos(:,2) = [chanlocs.Y];
        elec.chanpos(:,3) = [chanlocs.Z];
        elec.elecpos      = elec.chanpos;
        elec.label(:,1)   = [chanlocs.labels]';

        % Convert the fiducial position from voxel into CTF 
        nas = fiducial_locs.nas;
        lpa = fiducial_locs.lpa;
        rpa = fiducial_locs.rpa;

        vox2head = mri_acpc.transform;

        nas = ft_warp_apply(vox2head, nas, 'homogenous');
        lpa = ft_warp_apply(vox2head, lpa, 'homogenous');
        rpa = ft_warp_apply(vox2head, rpa, 'homogenous');

        % create a structure similar to a template set of electrodes
        fid.chanpos       = [nas; lpa; rpa];       % CTF head coordinates of fiducials
        fid.label         = {'nas','lhj','rhj'};    % use the same labels as those in elec
        fid.unit          = 'mm';                  % use the same units as those in mri
        fid.elecpos       = fid.chanpos;           % otherwise the electroderealign cannot find elecpos

        % alignment
        cfg               = [];
        cfg.viewmode      = 'surface';
        cfg.method        = 'fiducial';
        cfg.headshape     = vol;
        cfg.elec          = elec;                  % the electrodes we want to align
        cfg.elecstyle     = {'facecolor','red'};
        cfg.headmodelstyle = {'facecolor','b','edgecolor', 'none', 'facealpha', 0.4};
        cfg.template      = fid;                   % the template we want to align to
        cfg.fiducial      = {'nas', 'lhj', 'rhj'};  % labels of fiducials in fid and in elec
        elec_aligned_init = ft_electroderealign(cfg);

        figure;
        hold on;
        ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
        camlight
        ft_plot_sens(elec_aligned_init,'style','.r');
        ft_plot_sens(fid,'style','xb');%plot fiducial points
        savefig(fullfile( source,'elec_aligned_init.fig'))
        save(fullfile(source,'elec_aligned_init.mat'),'elec_aligned_init')
        % a dialog window to check inital alignment is good
    %     keyboard 

        cfg               = [];
        cfg.method        = 'project';
        cfg.elec          = elec_aligned_init;
        cfg.headshape     = vol;
        elec_aligned      = ft_electroderealign(cfg);

        figure;
        hold on;
        ft_plot_mesh(mesh,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.5)
        camlight
        ft_plot_sens(elec_aligned ,'style','.r');
        ft_plot_sens(fid,'style','xb');%plot fiducial points
        savefig(fullfile(source,'elec_aligned.fig'))
        save(fullfile(source,'elec_aligned.mat'), 'elec_aligned')
    end
end

