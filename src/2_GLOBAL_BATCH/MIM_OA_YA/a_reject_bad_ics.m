%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/MIM_OA/run_a_reject_bad_ics.sh

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
study_savedir = datetime;
study_savedir.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'MIM_OA'];
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
%## (MIND IN MOTION) DATASET SPECIFIC PARAMS (05/24/2023)
SUBJ_1YA = {'H1002','H1004','H1007','H1009',...
    'H1010','H1011','H1012','H1013','H1017',...
    'H1018','H1019','H1020','H1022','H1024',...
    'H1025','H1026','H1027','H1029','H1030','H1031',...
    'H1032','H1033','H1034','H1034','H1036',...
    'H1037','H1038','H1039','H1041','H1042',...
    'H1044','H1045','H1046','H1047','H1048'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2007','H2008',...
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
    'NH3076','NH3086','NH3090','NH3102',...
    'NH3104','NH3105','NH3106','NH3108','NH3110',...
    'NH3112','NH3113','NH3114','NH3123','NH3128',...
    };
SUBJ_ERROR = {'H2002'}; % (11/22/2023) 08202023_OAN82_iccRX0p65_iccREMG0p4_changparams
SUBJ_SLOW_WALKERS = {'H3042','H3046','H3047','H3073',...
    'H3092','NH3025','NH3051','NH3056','NH3071','NH3082'};
SUBJ_NO_MRI = {'H2010','H2012','H2018','H2036','H2041',...
    'H2072','H3018','H3120','NH3002','NH3009','NH3027','NH3129'};
SUBJ_MISSING_COND = {'H3024','NH3028'};
% SUBJ_UNKNOWN_ERR = {'NH3108','NH3030','NH3040','NH3025'};
% (08/21/2023) JS, 
% NH3108 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3030 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3040 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% NH3025 seems to bug out do to an indexing error during
% cropping (endTime = EEG.times( round(ExactCropLatencies))/1000;)
% (08/22/2023) JS, NH3108 bug seems to be related to entry errors in the
% Trial_Cropping_V2_test.xlsx sheet used to remove bad time ranges
% identified during collection. (fixed)
% NH3030 bug was due to how the CropTrialCheckFunc_checkLoadsol.m
% interpreted subject characters. It would consider 'NH3030_FU' as
% 'NH3030'. Changed from 'contains' to 'strcmp' func. (fixed)
% NH3040 bug was due to an entry error in Trial_Cropping_V2_test.xlsx (fixed)
SUBJ_DONT_INC = {'NH3004','NH3023'};
% (08/20/2023) JS, NH3004 has no headscan; NH3023 has no headscan clicks;
%- (OA) Subject Picks
SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%-
fprintf('Total subjects processing: %i\n',sum(cellfun(@(x) length({x{:}}),SUBJ_PICS)));
fprintf('Total subjects unable to be processed: %i\n',sum([length(SUBJ_NO_MRI),length(SUBJ_DONT_INC),length(SUBJ_SLOW_WALKERS)]));
%% (PARAMETERS) ======================================================== %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study group and saving
SESSION_NUMBER = '1';
DO_SAVE_ICA = false;
DO_SAVE_REDUCED_EEG = true;
%- Study & EEG Directory Naming
study_savedir = '10302023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p1';
study_fName = 'all_comps_study';
%- Subjects FULL ICA folder path
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
%% DEFINE PATHS
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_savedir)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
conditions      = cell(1,length([SUBJ_ITERS{:}]));
groups          = cell(1,length([SUBJ_ITERS{:}]));
sessions        = cell(1,length([SUBJ_ITERS{:}]));
subj_chars    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
chanlocs_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
% dipfit_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = cell(1,length([SUBJ_ITERS{:}]));
% vol_fPaths = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        try
            fNames{cnt} = tmp.name;
            %- Chanlocs fPaths
    %         chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
            chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
    %         dipfit_fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'head_model' filesep 'dipfit_struct.mat'];
%             dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm.mat'];
            dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm_ants.mat'];
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
    %         fprintf('DIPFIT Exists: %i\n',exist(dipfit_fPaths{cnt},'file'));
            fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{cnt},'file'));
        catch e
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('%s\n',getReport(e))
            dipfit_norm_fPaths{cnt} = [];
        end
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subj_chars{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = GROUP_NAMES{group_i};
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
inds = logical(cellfun(@(x) exist(x,'file'),dipfit_norm_fPaths));
chanlocs_fPaths = chanlocs_fPaths(inds);
dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
fPaths = fPaths(inds);
fNames = fNames(inds);
sessions = sessions(inds);
groups = groups(inds);
conditions = conditions(inds);
subj_chars = subj_chars(inds);
%%
ALLEEG = cell(length(fNames),1);
parfor(subj_i = 1:length(fNames),SLURM_POOL_SIZE)
    try
        EEG = format_eeg_ica(fNames{subj_i},fPaths{subj_i},subj_chars{subj_i},save_dir,...
                            'CONDITION',conditions{subj_i},...
                            'GROUP',groups{subj_i},...
                            'SESSION',sessions{subj_i},...
                            'DO_SAVE_ICA',DO_SAVE_ICA,...
                            'CHANLOCS_FPATHS',[],...
                            'HEADMODEL_METHOD','fem',...
                            'DIPFIT_MAT_F',dipfit_norm_fPaths{subj_i});    
        ALLEEG{subj_i} = EEG;
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'report. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
        exit();
    end
end
%- format cell output
[ALLEEG] = format_alleeg_cell(ALLEEG);
%%
%## vars
THRESH_BRAIN_SCORE = 8;
tmp_rmv_subjs = zeros(1,length(ALLEEG));
%## DIPOLE REJECTION
brain_ics_count = zeros(length(ALLEEG),1);
brain_ics_scores = cell(length(ALLEEG),1);
powpow_rej_count = zeros(length(ALLEEG),1);
all_rej_count = zeros(length(ALLEEG),1);
powpow_rej_nums = cell(length(ALLEEG),1);
brain_ic_nums = cell(length(ALLEEG),1);
all_rej_nums = cell(length(ALLEEG),1);
subj_str = cell(length(ALLEEG),1);
good_diplocs_orig = cell(length(ALLEEG),1);
bad_diplocs_orig = cell(length(ALLEEG),1);
bad_diplocs_mni = cell(length(ALLEEG),1);
good_diplocs_mni = cell(length(ALLEEG),1);
parfor(subj_i = 1:length(ALLEEG),SLURM_POOL_SIZE)
    try
        fprintf('Rejecting IC''s for subject %s...\n',ALLEEG(subj_i).subject);
        %- use rejection criteria to determine bad ic's
        reject_struct = mim_reject_ics(ALLEEG(subj_i),ALLEEG(subj_i).filepath);
        %- log good & bad components
        chk = (reject_struct.IC_all_brain >= THRESH_BRAIN_SCORE & reject_struct.IC_all_brain ~= 9);
        chk_w_powpow = unique(cat(1,find(chk),reject_struct.IC_powpow_rej));
        tmp_bad = setdiff(find((1:size(ALLEEG(subj_i).icaweights,1))),chk_w_powpow);
        tmp_good = chk_w_powpow';
        ALLEEG(subj_i).etc.urreject = [];
        ALLEEG(subj_i).etc.urreject.crit = [];
        ALLEEG(subj_i).etc.urreject.ic_keep = [];
        ALLEEG(subj_i).etc.urreject.ic_rej = [];
        ALLEEG(subj_i).etc.urreject.dipfit = [];
        if isempty(tmp_good)
            tmp_good = 0;
            fprintf('** Subject %s has 0 brain components\n',ALLEEG(subj_i).subject);
        else
            ALLEEG(subj_i).etc.urreject.crit = reject_struct;
            ALLEEG(subj_i).etc.urreject.ic_keep = tmp_good;
            ALLEEG(subj_i).etc.urreject.ic_rej = tmp_bad;
            ALLEEG(subj_i).etc.urreject.dipfit = ALLEEG(subj_i).dipfit;
            fprintf('** Subject %s has %i brain components\n',ALLEEG(subj_i).subject, length(tmp_good));
        end
        %- table printouts
        brain_ics_count(subj_i) = length(tmp_good);
        powpow_rej_count(subj_i) = length(reject_struct.IC_powpow_rej);
        all_rej_count(subj_i) = length(tmp_bad);
        tmp = reject_struct.IC_powpow_rej;
        if ~isempty(tmp)
            powpow_rej_nums{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %reject_struct.IC_powpow_rej;
        else
            powpow_rej_nums{subj_i} = '';
        end
        if ~isempty(tmp_good)
            brain_ic_nums{subj_i} = ['[' sprintf('%i,',tmp_good(1:end-1)) sprintf('%i',tmp_good(end)) ']']; %sprintf('%i,',tmp_good); %tmp_good;
        else
            brain_ic_nums{subj_i} = '';
        end
        if ~isempty(tmp_bad)
            all_rej_nums{subj_i} = ['[' sprintf('%i,',tmp_bad(1:end-1)) sprintf('%i',tmp_bad(end)) ']']; %sprintf('%i,',tmp_bad); %tmp_bad
        else
            all_rej_nums{subj_i} = '';
        end
        tmp = reject_struct.IC_all_brain;
        tmp = tmp(chk);
        if ~isempty(tmp)
            brain_ics_scores{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %sprintf('%i,',tmp(chk)); %tmp(chk);
        else
            brain_ics_scores{subj_i} = '';
        end
        subj_str{subj_i} = ALLEEG(subj_i).subject;
        %- remove IC's from struct
        if length(ALLEEG(subj_i).etc.urreject.ic_keep) < 2 
            fprintf('** Subject %s rejected.\n',ALLEEG(subj_i).subject);
            tmp_rmv_subjs(subj_i) = 1;
        else
            [ALLEEG(subj_i)] = rmv_ics(ALLEEG(subj_i))
            %- dipfit mods
            fprintf('making dipfit modifications\n');
            tmp = {ALLEEG(subj_i).dipfit.model(tmp_good).pos_old};
            if ~isempty(tmp)
                print_tmp = [];
                for i = 1:length(tmp)
                    print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
                end
                good_diplocs_orig{subj_i} = print_tmp;
            else
                good_diplocs_orig{subj_i} = '';
            end
            %- 
            tmp = {ALLEEG(subj_i).dipfit.model(tmp_bad).pos_old};
            if ~isempty(tmp)
                print_tmp = [];
                for i = 1:length(tmp)
                    print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
                end
                bad_diplocs_orig{subj_i} = print_tmp;
            else
                bad_diplocs_orig{subj_i} = '';
            end
            %-
            tmp = {ALLEEG(subj_i).dipfit.model(tmp_bad).mnipos};
            if ~isempty(tmp)
                print_tmp = [];
                for i = 1:length(tmp)
                    print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
                end
                bad_diplocs_mni{subj_i} = print_tmp; 
            else
                bad_diplocs_mni{subj_i} = '';
            end
            %-
            tmp = {ALLEEG(subj_i).dipfit.model(tmp_good).mnipos};
            if ~isempty(tmp)
                print_tmp = [];
                for i = 1:length(tmp)
                    print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
                end
                good_diplocs_mni{subj_i} = print_tmp; 
            else
                good_diplocs_mni{subj_i} = '';
            end
            %-
            ALLEEG(subj_i).dipfit.model = ALLEEG(subj_i).dipfit.model(goodinds);
        end
        %## SAVE EEG
        if DO_SAVE_REDUCED_EEG
            fprintf(1,'Saving Subject %s\n',ALLEEG(subj_i).subject);
            [ALLEEG(subj_i)] = pop_saveset(ALLEEG(subj_i),'savemode','twofiles',...
                'filename',ALLEEG(subj_i).filename,...
                'filepath',ALLEEG(subj_i).filepath);
        end
        fprintf('DONE. rejecting IC''s for subject %s.\n',ALLEEG(subj_i).subject);
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'report. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
        exit();
    end
end
ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
%% WRITE TO XLSX
fprintf('Writing subject rejection criteria table...\n');
tmp_table = table(brain_ics_count,powpow_rej_count,all_rej_count,powpow_rej_nums,brain_ic_nums,...
    all_rej_nums,brain_ics_scores,good_diplocs_orig,bad_diplocs_orig,good_diplocs_mni,bad_diplocs_mni,'RowNames',subj_str);
writetable(tmp_table,[save_dir filesep 'rejection_crit.xlsx'],'WriteRowNames',true,'WriteVariableNames',true) 
fprintf('DONE. Writing subject rejection criteria table.\n');
%%
fprintf('\n==== Making Study Modifications ====\n');
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName,...
                                'filename',study_fName,...
                                'filepath',save_dir);
% make sure all .mat files have a .fdt file associated with it.
% (03/08/23) why were these turned off? <-- to save memory!! (03/16/2023),
% use eeglab_options to set memory options so it doesn't conflict.
% (08/28/22) updatedat turnned off 
% (08/28/22) savedat turned off
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                    study_fName,save_dir,...
                                    'RESAVE_DATASETS','on');