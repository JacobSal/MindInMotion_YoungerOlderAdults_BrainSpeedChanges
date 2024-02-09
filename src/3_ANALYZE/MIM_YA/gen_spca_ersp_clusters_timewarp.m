%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_ersp_clusters.sh

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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('ya');
subject_chars = [SUBJ_PICS{:}];
% fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
% fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC)]));
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
cluster_study_dir = '12082023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
study_fName_1 = 'epoch_study';
spca_study_dir = '01122024_spca_analysis';
study_fName_2 = 'epoch_study';
%- study group and saving
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];

save_dir = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_1 = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_2 = [STUDIES_DIR filesep sprintf('%s',spca_study_dir)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% ===================================================================== %%
% SUB_DIR = [load_dir_1 filesep 'cluster'];
% if ~ispc
%     SUB_DIR = convertPath2UNIX(SUB_DIR);
% else
%     SUB_DIR = convertPath2Drive(SUB_DIR);
% end
tmp_dir = 'R:\Ferris-Lab\jsalminen\Experiments_Funding\Experiment_8_MIM_OA_YA_CRUNCH\data_saves\10022023_MIM_OAYA_N112_CRUNCH_gait';
SUB_DIR = [tmp_dir filesep 'cluster'];
%## USER SET
LOAD_DIFFERENT_STUDY = {true};
CLUSTER_K_PICKS = [12];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_5' filesep '12']};
CLUSTER_FILES = {'cl_inf_12.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
SUB_GROUP_FNAME = []; %'H3000'; %[]; %'H2000';
SUB_GROUP_FNAME_REGEX = []; %'H3000''s'; %[]; %'H2000''s';
CLUSTER_CLIM_MATCH = [];
k_i = 1;
%- convert cluster directory
if ~ispc
    cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
else
    cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
end
if ~isempty(SUB_GROUP_FNAME_REGEX)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
    plot_store_dir = [cluster_dir filesep 'plots_out' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    plot_store_dir = [cluster_dir filesep 'plots_out'];
end
if ~exist(spec_data_dir,'dir')
    error('spec_data dir does not exist');
end
if ~exist(plot_store_dir,'dir')
    mkdir(plot_store_dir);
end
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAMES{k_i})]);
    TMP_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[spec_data_dir filesep sprintf('%s.study',CLUSTER_STUDY_FNAMES{k_i})]);
    TMP_STUDY = tmp.STUDY;
end
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(TMP_STUDY);
% fPaths = {TMP_STUDY.datasetinfo.filepath};
% fNames = {TMP_STUDY.datasetinfo.filename};
condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
subj_chars = {TMP_STUDY.datasetinfo.subject};
fPaths = cell(length(subj_chars),1);
fNames = cell(length(subj_chars),1);
group_id = zeros(length(subj_chars),1);
for subj_i = 1:length(TMP_STUDY.datasetinfo)
    fPaths{subj_i} = [tmp_dir filesep subj_chars{subj_i} filesep 'GAIT_EPOCHED' filesep condition_gait{:}];
    fNames{subj_i} = TMP_STUDY.datasetinfo(subj_i).filename;
    tmp = regexp(subj_chars{subj_i},'\d','match');
    group_id(subj_i) = str2num(tmp{1});
end
%%
condition_gait = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%-
subj_i = 2;
cond_i = 1;
spca_fpath = [load_dir_2 filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep [condition_gait{:}]];
spca_ersp = par_load(spca_fpath,sprintf('cond%s_spca_ersp.mat',condition_gait{cond_i}));
%-
% allersp_clusters = cell(length(TMP_STUDY.cluster),length(condition_gait),length(GROUP_NAMES));
% allgpm_clusters = cell(length(TMP_STUDY.cluster),length(condition_gait),length(GROUP_NAMES));
allersp_clusters = zeros(length(TMP_STUDY.cluster),length(condition_gait),length(GROUP_NAMES),length(subject_chars),size(spca_ersp.ersp_corr,1),size(spca_ersp.ersp_corr,3));
allgpm_clusters = zeros(length(TMP_STUDY.cluster),length(condition_gait),length(GROUP_NAMES),length(subject_chars),size(spca_ersp.ersp_corr,1),size(spca_ersp.ersp_corr,3));
for subj_i = 1:length(subject_chars)
    %## LOAD EEG DATA
    try
        spca_fpath = [load_dir_2 filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep [condition_gait{:}]];
        EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        fprintf('Running subject %s\n',EEG.subject)
        g_i = group_id(subj_i);
        for cond_i = 1:length(condition_gait)
            ic_keep = EEG.etc.urreject.ic_keep;
            ic_rej = EEG.etc.urreject.ic_rej;
            spca_ersp = par_load(spca_fpath,sprintf('cond%s_spca_ersp.mat',condition_gait{cond_i}));
            for cl_i = 1:length(TMP_STUDY.cluster)
                set_i = (TMP_STUDY.cluster(cl_i).sets == subj_i);
                comp_i = TMP_STUDY.cluster(cl_i).comps(set_i);
                if ~isempty(comp_i) && length(comp_i) < 2
                    fprintf('%s) assigning IC %i(%i) to cluster %i in condition %s...\n',subject_chars{subj_i},ic_keep(comp_i),comp_i,cl_i,condition_gait{cond_i});
                    allersp_clusters(cl_i,cond_i,g_i,subj_i,:,:) = squeeze(spca_ersp.ersp_corr(:,ic_keep(comp_i),:));
                    allgpm_clusters(cl_i,cond_i,g_i,subj_i,:,:) = squeeze(spca_ersp.gpm_corr(:,ic_keep(comp_i),:));
%                     allgpm_clusters{cl_i,cond_i,g_i} = [allgpm_clusters{cl_i}, spca_ersp.gpm_corr(:,ic_keep(chk),:)];
                end
            end
        end
    catch e
        fprintf('\nError occured on subject %s\n%s\n',subject_chars{subj_i},getReport(e));
    end
end
par_save(allersp_clusters,spec_data_dir,'spca_ersp_clusters.mat');
par_save(allgpm_clusters,spec_data_dir,'spca_gpm_clusters.mat');