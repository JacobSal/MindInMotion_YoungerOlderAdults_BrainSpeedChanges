%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_BATCH_PREP/MIM_OA/run_ants_norm_dips.sh

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
run_dir = [source_dir filesep '1_BATCH_PREP' filesep 'MIM_OA'];
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
%% ===================================================================== %%
SUBJ_1YA = {'H1002','H1004','H1007','H1009',...
    'H1010','H1011','H1012','H1013','H1017',...
    'H1018','H1019','H1020','H1022','H1024',...
    'H1025','H1026','H1027','H1029','H1030','H1031',...
    'H1032','H1033','H1034','H1035','H1036',...
    'H1037','H1038','H1039','H1041','H1042',...
    'H1044','H1045','H1046','H1047','H1048'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2002','H2007','H2008',...
    'H2013','H2015','H2017','H2020','H2021',...
    'H2022','H2023','H2025','H2026','H2027',...
    'H2033','H2034','H2037','H2038','H2039',...
    'H2042','H2052','H2059','H2062','H2082',...
    'H2090','H2095','H2111','H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3053',...
    'H3063','H3072','H3077','H3103',...
    'H3107',...
    'NH3006','NH3007','NH3008','NH3010','NH3021',...
    'NH3026','NH3030','NH3036','NH3040',...
    'NH3041','NH3043','NH3054',...
    'NH3055','NH3058','NH3059','NH3066',...
    'NH3068','NH3069','NH3070','NH3074',...
    'NH3076','NH3086','NH3090','H3092','NH3102',...
    'NH3104','NH3105','NH3106','NH3108','NH3110',...
    'NH3112','NH3113','NH3114','NH3123','NH3128',...
    };
SUBJ_SLOW_WALKERS = {'H3042','H3046','H3047','H3073',...
    'H3092','NH3025','NH3051','NH3056','NH3071','NH3082'};
SUBJ_NO_MRI = {'H2010','H2012','H2018','H2036','H2041',...
    'H2072','H3018','H3120','NH3002','NH3009','NH3027','NH3129'};
SUBJ_MISSING_COND = {'H3024','NH3028'};
SUBJ_DONT_INC = {'NH3004','NH3023'};
SUBJ_PICS = {SUBJ_1YA,SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_1YA),1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
subj_names = {};
for i = 1:length(SUBJ_PICS)
    subj_names = [subj_names, SUBJ_PICS{i}];
end
%%
%## hard define
%- datset name
THRESHOLD_RV_BRAIN = 0.15;
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
OA_PREP_FPATH = '08202023_OAN82_iccRX0p60_iccREMG0p3_newparams'; 
%## soft define
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
parfor (subj_i = 1:length(subj_names),floor(length(subj_names)/2))
% for subj_i = 1:length(subj_names)
    subj_name = subj_names{subj_i};
    mri_path = [DATA_DIR filesep DATA_SET filesep subj_name filesep 'MRI'];
    in_fpath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'head_model'];
    out_fpath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'clean'];
    try
        %-
        norm_pos = readtable([in_fpath filesep 'dip_pos_outf.csv']);
        norm_pos = [norm_pos{:,:}];
        norm_chan = par_load(in_fpath, 'dip_pos.mat');
        %-
        fname = dir([out_fpath filesep '*.set']);
        EEG = pop_loadset('filepath',out_fpath,'filename',fname(1).name);
        %-
    %     trans_mat = load([mri_path filesep 'ants0GenericAffine.mat']);
    %     R_mat = reshape(trans_mat.AffineTransform_double_3_3,3,4);
    %     R_mat = [R_mat; 0 0 0 1];
    %     F_mat = eye(3,3);
    %     F_mat = [F_mat, trans_mat.fixed]; F_mat = [F_mat; 0 0 0 1];
        %-
        ants_mri = ft_read_mri([mri_path filesep 'antsWarped.nii.gz']);
        voxinds = round(ft_warp_apply(pinv(ants_mri.transform), norm_pos ));
        %-
        tmp = load([in_fpath filesep 'dipfit_struct']);
        dipfit_fem = tmp.SAVEVAR;
        %## Reformat Dipfit Structure
        empty_dip_struct = struct('posxyz',[nan(),nan(),nan()],'momxyz',[nan(),nan(),nan()],'rv',nan(),'diffmap',nan(),'sourcepot',nan(),'datapot',nan());
        EEG.dipfit_fem = [];
        EEG.dipfit_fem.model = empty_dip_struct;
        dipfit_fem_pos = zeros(length([dipfit_fem.component]),3);
        for i=1:length([dipfit_fem.component])
            %- 
            if ~isempty(dipfit_fem(i).dip)
                EEG.dipfit_fem.model(i).posxyz = dipfit_fem(i).dip.pos;
                EEG.dipfit_fem.model(i).momxyz = reshape(dipfit_fem(i).dip.mom, 3, length(dipfit_fem(i).dip.mom)/3)';
                if ~isempty(dipfit_fem(i).dip.rv)
                    EEG.dipfit_fem.model(i).rv     = dipfit_fem(i).dip.rv;
                else
                    EEG.dipfit_fem.model(i).rv     = nan();
                end
                %- 
                EEG.dipfit_fem.model(i).diffmap = dipfit_fem(i).Vmodel - dipfit_fem(i).Vdata;
                EEG.dipfit_fem.model(i).sourcepot = dipfit_fem(i).Vmodel;
                EEG.dipfit_fem.model(i).datapot   = dipfit_fem(i).Vdata;
                dipfit_fem_pos(i,:) = dipfit_fem(i).dip.pos;
            else
                EEG.dipfit_fem.model(i) = empty_dip_struct;
            end
        end
        %##
        for i=1:size(norm_chan,1)
            EEG.dipfit_fem.model(norm_chan.chan(i)).mnipos = norm_pos(i,:);
            EEG.dipfit_fem.model(norm_chan.chan(i)).mni_voxinds = voxinds(i,:);
            EEG.dipfit_fem.model(norm_chan.chan(i)).pos_old = [norm_chan{i,2:end}];
            EEG.dipfit_fem.model(norm_chan.chan(i)).posxyz = norm_pos(i,:);
%             disp(EEG.dipfit_fem.model(norm_chan.chan(i)));
        end
        %## SAVE
        dipfit_fem_norm = EEG.dipfit_fem;
        EEG.dipfit = EEG.dipfit_fem;
        %## Validitiy check
        IC_RV = vertcat(EEG.dipfit.model.rv);
        IC_POSXYZ = vertcat(EEG.dipfit.model.posxyz);
        ICs_RVthreshold_keep = find(IC_RV <= THRESHOLD_RV_BRAIN & all(~isnan(IC_POSXYZ),2));

        EEG = eeg_checkset(EEG);
        pop_saveset(EEG,'filepath',EEG.filepath,'filename',EEG.filename);
        par_save(dipfit_fem_norm,out_fpath,'dipfit_fem_norm_ants.mat');
    catch e
        fprintf('%s\n',getReport(e));
    end
   
end

