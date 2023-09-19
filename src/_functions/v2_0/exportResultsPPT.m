% Create study for cluster ICs. This code only works for cluster without
% using ERSP. Precompute ERSP needed to be done on Hipergator
% Chang Liu - 2021-11-23 - V1

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%   NJacobsen notes
%   When timewarping data, save values as EEG.timewarp = timewarp;
%   EEG.timewarp.medianlatency = median(timewarp.latencies(:,:));%Warping to the median latency of my 5 events
%   By default, std_ersp will use the median of all subject's
%   timewarp.latencies(:,:) as 'timewarpms' unless individual subject 
%   warpto is indiciated using 'timewarpms', 'subject tw matrix'
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230417.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_b_ic_clustering_refine.sh

%{
%## RESTORE MATLABs
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
    SLURM_POOL_SIZE = 1;
end
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
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
SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
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
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
% NOTE: see. statcondfieldtrip.m or std_stat.m
COND_EVENT_CHAR = 'cond';
%- clustering parameters
MIN_ICS_SUBJ = [2,3,4,5,6,7,8]; % iterative clustering
DISTS_FROM_MEDIAN = [10,20,30,40,50,60,70]; % distance in units (mm)
K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
cluster_ks = K_RANGE(1):K_RANGE(2);
% (08/21/2023) JS, this currenlty doesn't do anything but take up more
% memory.
REPEATED_CLUSTERING_STD = 3;
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',[],...
    'filepath',[],...
    'outliers',inf(),...
    'maxiter',200);
%- 
%- custom params
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
STD_PRECLUST_COMMAND = {'dipoles','weight',clustering_weights.dipoles};
%- iterative clustering parameters
n_iterations = 50;
outlier_sigma = 3;
%- datetime override
% dt = '05182023_MIM_OA_subset_N85_oldpipe';
% dt = '05192023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06122023_MIM_OAN79_subset_prep_verified_gait';
% dt = '06282023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07112023_MIM_OAN79_subset_prep_verified_gait';
% dt = '07152023_MIM_OAN79_subset_prep_verified_gait';
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
for i = 1:length(cluster_ks)
    clust_i = cluster_ks(i);
    for j = 1:length(DISTS_FROM_MEDIAN)
        dist_i = DISTS_FROM_MEDIAN(j);
        cluster_dir = [save_dir filesep clustering_method filesep num2str(clust_i) filesep sprintf('dist%i',dist_i)];
        input_folder = 
    end
end

inputFolder = 'E:\Splitbelt Processing\STUDY-CustomHDM-EEG\Figures\PSD vs SLA\LME\fits\jpg'; %change directories to folder with images
plotPSD =0; % [0|1]
plotERSP =0; % [0|1]
plotLME = 1;
clear files
project_dir = 'R:\Ferris-Lab\jacobsen.noelle\Split Belt Pilot Study\STUDY-CustomHDM-EEG\Figures'; % Figure folder
dip_dir = [project_dir,'\Dipoles'];
topo_dir = [project_dir,'\Topo'];
outputFolder = project_dir ;

if plotPSD ==1
    files = dir( fullfile(inputFolder, '*.jpg')); % List all current folder contents ending with .png. Resulting names will
elseif plotERSP ==1
files = dir( fullfile(inputFolder, '*.png')); %ersp or Psd files
elseif plotLME==1
   files = dir( fullfile(inputFolder, '*.jpg'));
   CEfolder = 'E:\Splitbelt Processing\STUDY-CustomHDM-EEG\Figures\PSD vs SLA\LME\Coeff\jpg';
   CEfiles =  dir( fullfile('E:\Splitbelt Processing\STUDY-CustomHDM-EEG\Figures\PSD vs SLA\LME\Coeff\jpg','*.jpg'));
    
end



dip_files= dir( fullfile(dip_dir , '*.jpg'));
topo_files = dir( fullfile(topo_dir , '*.png'));
files = [files; topo_files;dip_files];
T = struct2table(files); % convert the struct array to a table
sortedT = sortrows(T, 'name'); % sort the table by 'name'
files = table2struct(sortedT); % change it back to struct array if necessary

% appear in the order returned by the operating system.
% Prompt the user for the PowerPoint file to amend
[fn, pn] = uigetfile('*.pptx', 'Select PowerPoint File To Amend');
filename = fullfile(pn, fn);

% getplotParams;
% p = 21; %set study design

%% Test 1
% Create Common Object Model (COM) server so MATLAB can export data to
% PowerPoint
% add study results to powerpoint

g = actxserver('powerpoint.application');
g.Visible = 1;% Open PowerPoint and make it visible
Presentation = g.Presentation;
Presentation = invoke(Presentation, 'open', filename);
slide_count = get(Presentation.Slides, 'Count');% Get current number of slides
% Export all PNGs in the current directory to the PowerPoint file specified
% above. The following slides will be added to the END of the PowerPoint
% file. All slides will have a common title.
% EMG_labels = [{'sterno_L'},{'splen_cap_L'},{'trap_lat_L'},{'trap_med_L'},{'trap_med_R'},{'trap_lat_R'},{'splen_cap_R'},{'sterno_L'}];
i = 0;
%for p = 2
   for CL = [3:10 12:length(STUDY.cluster)]
   
%     for C=1:length(channels_to_plot)
%         chan = channels_to_plot{1,C};
%         chan = cellstr(chan);
%         chan = chan{1,1};
        %       %image index from file list
        D1_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_side']));
        D2_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_top']));
        D3_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_back']));
        topo_ind = find(contains({files.name},['topo_Cls',num2str(CL)]));
      
        %load images
        slide_count = int32(double(slide_count)+1);
        slide = invoke(Presentation.Slides, 'Add', slide_count(end), 11); %slide title
        try %check if there's an anatomical label
            label = cellstr(STUDY.cluster(CL).label);
            label = label{1,1};
            set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL),'-',label));
        catch %no label
            set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL)));
%             label = cellstr(EMG_labels{1,C});
%             label = label{1,1};
%             set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('Chan ',num2str(C),'-',label));
        end
        D1 = fullfile(dip_dir, files(D1_ind).name);
        D2 = fullfile(dip_dir, files(D2_ind).name);
        D3 = fullfile(dip_dir, files(D3_ind).name);
        topo = fullfile(topo_dir, files(topo_ind).name);
        
        
        %add images to slide - PSD, topo, dipoles
        if plotPSD ==1
%       PSD_ind = find(contains({files.name},[plotParams(p).figname,chan,'_100Hz']));
        PSD_ind = find(contains({files.name},[plotParams(p).figname,num2str(CL),'_30Hz']));
        PSD = fullfile(inputFolder, files(PSD_ind).name);
        Image{i+1} = slide.Shapes.AddPicture(PSD, 'msoFalse', 'msoTrue',118.75,166,712.9, 374.4); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',400,36,108, 108); %Left, top, width, height
        end
        
        %add images to next slide - ERSPs
        %       slide_count = int32(double(slide_count)+1);
        %       slide = invoke(Presentation.Slides, 'Add', slide_count(end), 11);
        %       set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL),'-',STUDY.cluster(CL).label{:}));
              %Image{i+1} = slide.Shapes.AddPicture(ERSP, 'msoFalse', 'msoTrue',0,122.4,950.6, 374.4); %Left, top, width, height
        
              if plotERSP ==1
%          %add ERSPs with dipoles and topo
%          ERSP_ind = find(contains({files.name},['ERSP_B3_A1_perturb_Cls',num2str(CL),'_condVsBaseline_perm0.05_cluster_fieldtrip_contour']));
        ERSP_ind = find(contains({files.name},['Cls' num2str(CL)]) & contains({files.name},'ERSP') & contains({files.name}, 'masked'));
        ERSP = fullfile(inputFolder, files(ERSP_ind).name);
        Image{i+1} = slide.Shapes.AddPicture(ERSP, 'msoFalse', 'msoTrue',0,152.4,950.6, 374.4); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
        Image{i+1} = slide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',400,36,108, 108); %Left, top, width, height
              end

              if plotLME ==1
                  fig_ind = find(contains({files.name},['CL' num2str(CL)]) & contains({files.name},'PSDvsDSR'));
                  fig = fullfile(inputFolder, files(fig_ind).name);
                  Image{i+1} = slide.Shapes.AddPicture(fig, 'msoFalse', 'msoTrue',0,152.4,575, 374.4); %Left, top, width, height
                  Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
                  Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
                  Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
                  Image{i+1} = slide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',400,36,108, 108); %Left, top, width, height
                  fig_ind = find(contains({CEfiles.name},['CL' num2str(CL)]) & contains({CEfiles.name},'PSDvsDSR'));
                  CEFig = fullfile(CEfolder, CEfiles(fig_ind).name);
                  Image{i+1} = slide.Shapes.AddPicture(CEFig, 'msoFalse', 'msoTrue',518,152.4,425, 374.4); %Left, top, width, height
              
              end

                  %
%           %add images to next slide - ERSPs
%         slide_count = int32(double(slide_count)+1);
%         slide = invoke(Presentation.Slides, 'Add', slide_count(end), 11); %slide title
%          try %check if there's an anatomical label
%             label = cellstr(STUDY.cluster(CL).label);
%             label = label{1,1};
%             set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL),'-',label));
%         catch %no label
%             set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL)));
%          end
%         ERSP_ind = find(contains({files.name},['ERSP_B3vA1_Cls',num2str(CL)]));
%         ERSP = fullfile(inputFolder, files(ERSP_ind).name);
%         Image{i+1} = slide.Shapes.AddPicture(ERSP, 'msoFalse', 'msoTrue',0,152.4,950.6, 374.4); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
%         Image{i+1} = slide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',400,36,108, 108); %Left, top, width, height
        
        % %% lots of PSD, comparison
        %        %image index from file list
        % %       D1_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_side']));
        % %       D2_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_top']));
        % %       D3_ind = find(contains({files.name},['CLS_',num2str(CL),'_comps_cents_back']));
        %       %topo_ind = find(contains({pngfiles.name},['topo_Cls',num2str(CL)]));
        % %      PSD_ind =  find(contains({files.name},[plotParams(p).figname,num2str(CL),'_100Hz']));
        % %       PSD1_ind = find(contains({files.name},['adapt1_B3_PSD_Cls',num2str(CL),'_50Hz']));
        % %       PSD2_ind = find(contains({files.name},['adapt2_PSD_Cls',num2str(CL),'_50Hz']));
        %       %ERSP_ind = find(contains({pngfiles.name},['ERSP_Cls',num2str(CL)]));
        %
        %       %load images
        %       slide_count = int32(double(slide_count)+1);
        %       slide = invoke(Presentation.Slides, 'Add', slide_count(end), 11); %slide title
        %       %set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL)));
        %       try %check if there's an anatomical label
        %       label = cellstr(STUDY.cluster(CL).label);
        %       label = label{1,1};
        %       set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL),'-',label));
        %       catch %no label
        %         set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL)));
        %       end
        % %       D1 = fullfile(project_dir, files(D1_ind).name);
        % %       D2 = fullfile(project_dir, files(D2_ind).name);
        % %       D3 = fullfile(project_dir, files(D3_ind).name);
        % %       PSD = fullfile(project_dir, files(PSD_ind).name);
        % %       PSD1 = fullfile(project_dir, files(PSD1_ind).name);
        % %       PSD2 = fullfile(project_dir, files(PSD2_ind).name);
        % %       topo = fullfile(project_dir, pngfiles(topo_ind).name);
        % %       ERSP = fullfile(project_dir, pngfiles(ERSP_ind).name);
        %
        %       %add images to slide
        %       Image{i+1} = slide.Shapes.AddPicture(PSD, 'msoFalse', 'msoTrue',0,166,950.4, 374.4); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',518, 36,144, 108); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',669.2, 36,108, 108); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',792,36,144, 108); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(topo, 'msoFalse', 'msoTrue',36,126,126, 126); %Left, top, width, height
        %
        %
        % %
        % %       slide_count = int32(double(slide_count)+1);
        % %       slide = invoke(Presentation.Slides, 'Add', slide_count(end), 11);
        % %       set(slide.Shapes.Title.Textframe.Textrange, 'Text', strcat('CL ',num2str(CL),'-',STUDY.cluster(CL).label{:}));
        %       %Image{i+1} = slide.Shapes.AddPicture(ERSP, 'msoFalse', 'msoTrue',0,122.4,950.6, 374.4); %Left, top, width, height
        %
        % %% PSD comparison
        % %       Image{i+1} = slide.Shapes.AddPicture(PSD1, 'msoFalse', 'msoTrue',178,100.64,522,224.64); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(PSD2, 'msoFalse', 'msoTrue',178,311,522,224.64); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D1, 'msoFalse', 'msoTrue',593.8, 36,126, 88); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D2, 'msoFalse', 'msoTrue',719.8, 36,98, 88); %Left, top, width, height
        % %       Image{i+1} = slide.Shapes.AddPicture(D3, 'msoFalse', 'msoTrue',817.8,36,126, 88); %Left, top, width, height
        % %
        
    end
%end
outfile = fullfile(outputFolder, ['Results-',date,'.pptx']);% Save the amended PowerPoint presentation to the current directory
Presentation.SaveAs(outfile);
g.Quit;% Close PowerPoint as a COM automation server
g.delete;
fprintf('\nResults powerpoint saved here:%s',outputFolder);