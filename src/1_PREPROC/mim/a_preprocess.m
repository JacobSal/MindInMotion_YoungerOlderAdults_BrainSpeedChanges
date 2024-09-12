%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 

%- run script
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROCE/mim/run_a_preprocess.sh

%- run amica
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROCE/mim/b_run_singlenode_amica.sh

%- mim dipfit
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROCE/mim/c_run_mim_mcc_dipfit.sh

%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = fileparts(fileparts(STUDY_DIR));
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
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
SUBJ_PICS = {{'H2012_FU','H2018_FU'},{'H3120','NH3129'}};
GROUP_NAMES = {'H2000''s','H3000''s'};
SUBJ_ITERS = {1:length(SUBJ_PICS{1}),1:length(SUBJ_PICS{2})};
%% (PROCESSING PARAMS) ================================================= %%
%## hard define
%- dataset name
DATA_SET = 'MIM_dataset';
%- datetime override
% OA_PREP_FNAME = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FNAME = '07122023_OAN79_iccRX0p9_iccREMG0p3'; % JACOB,SAL(07/12/2023)
% OA_PREP_FNAME = '07142023_OAN79_iccRX0p55_iccREMG0p3_changparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(10/26/2023)
% OA_PREP_FNAME = '10302023_OAN82_iccRX0p60_iccREMG0p4_newparams'; % JACOB,SAL(10/30/2023)
% OA_PREP_FNAME = 'EMG_ANALYSIS'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '08202023_OAN82_iccRX0p60_iccREMG0p3_newparams'; % JACOB,SAL(07/14/2023)
% OA_PREP_FNAME = '11262023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(07/14/2023)
OA_PREP_FNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
%## soft define
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET]; % JACOB,SAL(02/23/2023)
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',OA_PREP_FNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (FIND FILES) ======================================================== %%
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
cnt = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        if ~isempty(tmp)
            fNames{cnt} = tmp.name;
        else
            fNames{cnt} = [];
        end
        %- Prints
        fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
        fprintf('EEG Exists: %i\n',exist([fPaths{cnt} filesep fNames{cnt}],'file')>1)
        stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
        cnt = cnt + 1;
    end
end
tmp = ~cellfun(@isempty,fNames);
fNames = fNames(tmp);
tmp = ~cellfun(@isempty,fPaths);
fPaths = fPaths(tmp);
tmp = ~cellfun(@isempty,subjectNames);
subjectNames = subjectNames(tmp);
%% SET POOLSIZE
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(fPaths)]);
else
    POOL_SIZE = 1;
end
%% (PARFOR) GENERATE CONNECTIVITY METRICS 
LOOP_VAR = 1:length(fPaths);
amica_cmd = cell(length(fPaths),1);
params = cell(length(fPaths),1);
subj_pick = 'H3092';
% parfor (subj_i = LOOP_VAR,POOL_SIZE)
for subj_i = find(strcmp(subjectNames,'NH3113'))
    fprintf('Running subject %s...\n',subjectNames{subj_i})
    %## PREP for MAIN_FUNC
    if ~exist([save_dir filesep subjectNames{subj_i}],'dir')
        mkdir([save_dir filesep subjectNames{subj_i}]);
    end
    %## RUN MAIN_FUNC
    try
        [EEG,amica_cmd{subj_i},params{subj_i}] = main_func(subjectNames{subj_i},fPaths{subj_i},...
                [save_dir filesep subjectNames{subj_i}],STUDIES_DIR);
        fprintf('%s\n',amica_cmd{subj_i}{2})
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,subjectNames{subj_i},getReport(e));
    end
end
