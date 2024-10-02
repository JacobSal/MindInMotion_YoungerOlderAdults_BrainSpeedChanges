function [err] = mcc_epoch_process(ica_dir,study_dir,jf_fpath,varargin)
%MIM_SPCA_MCC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%       REQUIRED:
%       ica_dir, CHAR
%           directory where ICA'd EEG data is located for each subject
%           you'd like to run.
%           E.G.,
%               my_ica_folder/
%               ├─ subj_1/
%               │  ├─ clean/
%               ├─ subj_2/
%               │  ├─ clean/
%               ├─ subj_n/
%               │  ├─ clean/
%       study_dir, CHAR
%           directory where the program will store information about the
%           subjects in 'ica_dir'.
%
%       jf_fpath, CHAR
%           JSON Format is as follows:
%               {
% 	                 "EPOCH_PARAMS": {
% 		                "percent_overlap": 0,
% 		                "epoch_event_char": "RHS",
% 		                "epoch_time_lims": [-0.5,4.5],
% 		                "tw_stdev": 3,
% 		                "tw_events": ["RHS","LTO","LHS","RTO","RHS"],
% 		                "path_ext": "gait_epoched",
% 		                "gait_trial_chars": ["0p25","0p5","0p75","1p0","flat","low","med","high"],
% 		                "rest_trial_char": [],
% 		                "do_recalc_epoch": true
% 	                },
% 	                "SUBJ_CHARS": null
%               }
% 
%   OUT: 
%   IMPORTANT: 
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/run_mim_mcc_dipfit_exe.sh
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
%% (PARAMETERS) ======================================================== %%
%## hard define
cat_logo();
str = '';
err = 0;
ADD_CLEANING_SUBMODS = false;
ADD_DIPFIT_COMPILE_SUBMODS = false;
SAVE_SEPERATE_EEG = false;
DO_SLIDING_WINDOW = false;
study_fname_icared = 'cont_icreduced_study';
study_fname_gait = 'epoch_study';
%## EPOCH PARAMS
EPOCH_PARAMS = struct('percent_overlap',0,...
    'epoch_event_char','RHS',...
    'epoch_time_lims',[-0.5,4.5],...
    'tw_stdev',3,...
    'tw_events',{{'RHS','LTO','LHS','RTO','RHS'}},...
    'path_ext','gait_epoched',...
    'gait_trial_chars',{{'0p25','0p5','0p75','1p0','flat','low','med','high'}},...
    'rest_trial_char',{{}},...
    'do_recalc_epoch',true);
%## SUBJECT CHARACTERS LOAD-IN FROM FUNCTION
[SUBJ_CHARS,GROUP_NAMES,~,~,~,~,~] = mim_dataset_information('yaoa_spca');
SUBJ_CHARS = [SUBJ_CHARS{:}];
%% LOAD JSON & ASSIGN CHECKS
%- connectivity statistics & surrogates params
if ~isempty(jf_fpath)
    fid = fopen(jf_fpath,'r');
    raw = fread(fid,inf); 
    str = char(raw');
    fclose(fid); 
    params = jsondecode(str);
    if ~isempty(params.EPOCH_PARAMS)
        tmp_epoch = params.EPOCH_PARAMS;
        varargin = [varargin, 'EPOCH_PARAMS', tmp_epoch];
    end
    if ~isempty(params.SUBJ_CHARS)
        tmp_subj_chars = params.SUBJ_CHARS;
        varargin = [varargin, 'SUBJ_CHARS', tmp_subj_chars];
    end
end
%## working directory containing subject ICA .set files
errorMsg = 'Value must be CHAR. Working directory containing all subjects to be epoched & sPCA''d.'; 
id_validFcn = @(x) assert(ischar(x)  && exist(x,'dir'),errorMsg);
errorMsg = 'Value must be CHAR. New directory to store epoched .set files and study information.'; 
sd_validFcn = @(x) assert(ischar(x),errorMsg);
errorMsg = 'Value must be CHAR. File path to a .json file containing parameters to be used.'; 
jf_validFcn = @(x) assert((ischar(x) && exist(x,'file')) || isempty(x),errorMsg);
%% (PARSER) ============================================================ %%
fprintf('ica_dir: %s\n',ica_dir);
fprintf('study_dir: %s\n',study_dir);
fprintf('jf_fpath: %s\n',jf_fpath);
fprintf('%s\n',str);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ica_dir',id_validFcn);
addRequired(p,'study_dir',sd_validFcn);
addRequired(p,'jf_fpath',jf_validFcn);
%## OPTIONAL
addParameter(p,'SUBJ_CHARS',SUBJ_CHARS,@iscell);
addParameter(p,'EPOCH_PARAMS',EPOCH_PARAMS,@(x) validate_struct(x,EPOCH_PARAMS));
%## PARAMETERS
parse(p,ica_dir,study_dir,jf_fpath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
SUBJ_CHARS = p.Results.SUBJ_CHARS;
EPOCH_PARAMS = p.Results.EPOCH_PARAMS;
if ~exist(study_dir,'dir')
    mkdir(study_dir);
end
%% PARPOOL SETUP ======================================================= %%
[PATHS,SLURM_POOL_SIZE,ALLEEG,STUDY,CURRENTSTUDY,CURRENTSET] = set_workspace();
%% (GRAB SUBJECT .SET & DIPFIT FILES) ==================================== %%
fnames        = cell(1,length(SUBJ_CHARS));
fpaths        = cell(1,length(SUBJ_CHARS));
conditions    = cell(1,length(SUBJ_CHARS));
groups        = cell(1,length(SUBJ_CHARS));
sessions      = cell(1,length(SUBJ_CHARS));
chk           = zeros(1,length(SUBJ_CHARS));
for subj_i = 1:length(SUBJ_CHARS)
    %- ICA fPaths
    fpaths{subj_i} = [ica_dir filesep SUBJ_CHARS{subj_i} filesep 'clean'];
    tmp = dir([fpaths{subj_i} filesep '*.set']);
    conditions(subj_i) = join(EPOCH_PARAMS.gait_trial_chars,'-');
    val = regexp(SUBJ_CHARS{subj_i},'\d*');
    gi = contains(GROUP_NAMES,SUBJ_CHARS{subj_i}(val));
    groups{subj_i} = GROUP_NAMES{gi};
    sessions{subj_i} = '1'; %- (04/12/2024) default session number 
    try
        fnames{subj_i} = tmp.name;
        %- Prints
        fprintf('%s) .set fpath: %s\n',SUBJ_CHARS{subj_i}, fpaths{subj_i});
        chk(subj_i) = 1;
    catch e
        fprintf('Error on subject %s\n',SUBJ_CHARS{subj_i})
        fprintf('%s\n',getReport(e))
    end
end
%- remove subjects without a dipole fit
inds = logical(chk);
fpaths = fpaths(inds);
fnames = fnames(inds);
SUBJ_CHARS = SUBJ_CHARS(inds);

%% Create STUDY & ALLEEG structs ======================================= %%
try
    fprintf('Creating ALLEEG...\n');
    [MAIN_ALLEEG] = mim_create_alleeg(fnames,fpaths,SUBJ_CHARS,study_dir,...
                        conditions,groups,sessions);
catch e
    fprintf('\n%s\n',getReport(e));
    exit();
end
try
    fprintf('Generating STUDY...\n');
    [MAIN_STUDY,MAIN_ALLEEG] = mim_create_study(MAIN_ALLEEG,study_fname_icared,study_dir);
catch e
    fprintf('\n%s\n',getReport(e));
    exit();
end
%## Update
fpaths = {MAIN_ALLEEG.filepath};
fnames = {MAIN_ALLEEG.filename};

%% ===================================================================== %%
%## GENERATE EPOCH MAIN FUNC
alleeg_fpaths = cell(1,length(fpaths));
tmp = cell(1,length(fpaths));
rmv_subj = zeros(1,length(fpaths));
%## PARFOR LOOP
parfor (subj_i = 1:length(fpaths),SLURM_POOL_SIZE)
% for subj_i = 1:length(subjectNames)
    %## LOAD EEG DATA
    EEG = pop_loadset('filepath',fpaths{subj_i},'filename',fnames{subj_i});
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    %## PARSE TRIALS
    epoched_fPath = [study_dir filesep EEG.subject filesep EPOCH_PARAMS.path_ext];
    fPath = [epoched_fPath filesep [EPOCH_PARAMS.gait_trial_chars{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG.set',EEG.subject,[EPOCH_PARAMS.gait_trial_chars{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    try
        %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
        if isfield(EEG.event,'trialName')
            EEG.event = rmfield(EEG.event,'trialName');
        end
        if isfield(EEG.event,'channel')
            EEG.event = rmfield(EEG.event,'channel');
        end
        if isfield(EEG.event,'code')
            EEG.event = rmfield(EEG.event,'code');
        end
        if isfield(EEG.event,'bvtime')
            EEG.event = rmfield(EEG.event,'bvtime');
        end
        if isfield(EEG.event,'bvmknum')
            EEG.event = rmfield(EEG.event,'bvmknum');
        end
        if isfield(EEG.event,'datetime')
            EEG.event = rmfield(EEG.event,'datetime');
        end
        %## EPOCH
        [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,false,...
            'EPOCH_TIME_LIMITS',EPOCH_PARAMS.epoch_time_lims,... 
            'STD_TIMEWARP',EPOCH_PARAMS.tw_stdev,...
            'COND_CHARS',EPOCH_PARAMS.gait_trial_chars);

        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
        fprintf('\nConcatenating datasets...\n');
        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_SEPERATE_EEG
            for i = 1:length(ALLEEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,ALLEEG(i).condition)];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(ALLEEG(i),'savemode','twofiles',...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition);
            end
            alleeg_fpaths{subj_i} = cond_files;
        end
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## timewarp for across condition
        if ~DO_SLIDING_WINDOW
            timewarp = make_timewarp(ALLEEG,EPOCH_PARAMS.tw_events,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select(ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_PARAMS.epoch_time_lims;
        end
        %## STRUCT EDITS
        ALLEEG.urevent = []; % might be needed
        ALLEEG.etc.epoch.epoch_limits = EPOCH_PARAMS.epoch_time_lims;
        %- checks
        ALLEEG = eeg_checkset(ALLEEG,'eventconsistency');
        ALLEEG = eeg_checkset(ALLEEG);
        ALLEEG = eeg_checkamica(ALLEEG);
        %- save
        [ALLEEG] = pop_saveset(ALLEEG,'savemode','twofiles',...
                'filename',fName,...
                'filepath',fPath,...
                'version','6');
        tmp{subj_i} = ALLEEG;
    catch e
        rmv_subj(subj_i) = 1;
        EEG.timewarp = struct([]);
        EEG.urevent = [];
        tmp{subj_i} = []; %EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
    end
end
%% ===================================================================== %%
%## SAVE BIG STUDY
fprintf('==== Reformatting Study ====\n');
%- remove bugged out subjects
fprintf('Bugged Subjects: %s',MAIN_ALLEEG(cellfun(@isempty,tmp)).subject);
tmp = tmp(~cellfun(@isempty,tmp));
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(tmp));
for subj_i = 1:length(tmp)
    fss{subj_i} = fields(tmp{subj_i})';
    disp(size(fields(tmp{subj_i})'));
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(tmp)
    EEG = tmp{subj_i};
    fs = fields(EEG);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false); 
    out = [out{:}];
    addFs = fsPrev(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
%             fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j})
        end
    end 
    tmp{subj_i} = orderfields(EEG);
end
%- CONCATENATE tmp
tmp = cellfun(@(x) [[]; x], tmp);
%##
[STUDY, ALLEEG] = std_editset([],tmp,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fname_gait,...
                                'filename',study_fname_gait,...
                                'filepath',study_dir);
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
STUDY.etc.a_epoch_process = EPOCH_PARAMS;
STUDY.etc.a_epoch_process.epoch_alleeg_fpaths = alleeg_fpaths;
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
end
%% (SUBFUNCTIONS) ====================================================== %%
function [PATHS,SLURM_POOL_SIZE,ALLEEG,STUDY,CURRENTSTUDY,CURRENTSET] = set_workspace()
    %## TIME
    TT = tic;
    %## ================================================================= %%
    TMP_PWD = dir(['.' filesep]);
    TMP_PWD = TMP_PWD(1).folder;
    fprintf(1,'Current folder: %s\n',TMP_PWD);
    %## datetime
    dt         = datetime;
    dt.Format  = 'ddMMyyyy';
    %## PARAMS ========================================================== %%
    tmp = strsplit(TMP_PWD,filesep);
    src_ind = find(strcmp(tmp,'src'));
    if ~ispc
        %- Add source directory where setup functions are 
        src_dir = [filesep strjoin(tmp(1:src_ind),filesep)];
    else
        %- Add source directory where setup functions are
        src_dir = strjoin(tmp(1:src_ind),filesep);
    end
    if ~ispc
        %- Add submodules directory where packages are 
        submodules_dir = [filesep strjoin(tmp(1:src_ind-1),filesep) filesep 'submodules'];
    else
        %- Add submodules directory where packages are 
        submodules_dir = [strjoin(tmp(1:src_ind-1),filesep) filesep 'submodules'];
    end
    %## FUNCTIONS FOLDER
    FUNC_FPATH = [src_dir filesep '_functions' filesep 'v2_0'];
    %##
    fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n-FUNCTIONS: %s\n',src_dir,submodules_dir,FUNC_FPATH);
    %## HARDCODE PATHS STRUCT =========================================== %%
    PATHS = [];
    %- submods path
    PATHS.submods_dir = submodules_dir;
    %- src folder
    PATHS.src_dir = src_dir;
    %- _data folder
    PATHS.data_dir = fullfile(src_dir,'_data');
    %- EEGLAB folder
    PATHS.eeglab_dir = [submodules_dir filesep 'eeglab'];
    %- _functions folder
    PATHS.functions_dir = FUNC_FPATH;
    %## ADDPATH for FIELDTRIP =========================================== %%
    ft_defaults;
    %- always start eeglab last.
    ALLEEG=[]; STUDY=[]; CURRENTSET=0; CURRENTSTUDY=0;
    eeglab;
    %## PARPOOL SETUP =================================================== %%
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
        mkdir([TMP_PWD filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')])
        pp.JobStorageLocation = [TMP_PWD filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')];
        %- create your p-pool (NOTE: gross!!)
        pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
    else
        pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
        'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
        'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
        SLURM_POOL_SIZE = 1;
    end
    %## TIME
    fprintf('Done with workplace setup: %0.2f',toc(TT));
end

%% (SUBFUNCTIONS) ====================================================== %%
function [p] = unix_genpath(d)
    %GENPATH Generate recursive toolbox path.
    %   P = GENPATH returns a character vector containing a path name 
    %   that includes all the folders and subfolders below MATLABROOT/toolbox, 
    %   including empty subfolders.
    %
    %   P = GENPATH(FOLDERNAME) returns a character vector containing a path 
    %   name that includes FOLDERNAME and all subfolders of FOLDERNAME, 
    %   including empty subfolders.
    %   
    %   NOTE 1: GENPATH will not exactly recreate the original MATLAB path.
    %
    %   NOTE 2: GENPATH only includes subfolders allowed on the MATLAB
    %   path.
    %
    %   See also PATH, ADDPATH, RMPATH, SAVEPATH.

    %   Copyright 1984-2018 The MathWorks, Inc.
    %------------------------------------------------------------------------------

    % String Adoption
    if nargin > 0
        d = convertStringsToChars(d);
    end

    if nargin==0
      p = genpath(fullfile(matlabroot,'toolbox'));
      if length(p) > 1, p(end) = []; end % Remove trailing pathsep
      return
    end

    % initialise variables
    classsep = '@';  % qualifier for overloaded class directories
    packagesep = '+';  % qualifier for overloaded package directories
    p = '';           % path to be returned

    % Generate path based on given root directory
    files = dir(d);
    if isempty(files)
      return
    end

    % Add d to the path even if it is empty.
    p = [p d pathsep];

    % set logical vector for subdirectory entries in d
    isdir = logical(cat(1,files.isdir));
    %
    % Recursively descend through directories which are neither
    % private nor "class" directories.
    %
    dirs = files(isdir); % select only directory entries from the current listing

    for i=1:length(dirs)
       dirname = dirs(i).name;
       if    ~strcmp( dirname,'.')          && ...
             ~strcmp( dirname,'..')         && ...
             ~strncmp( dirname,classsep,1) && ...
             ~strncmp( dirname,packagesep,1) && ...
             ~strcmp( dirname,'private')    && ...
             ~strcmp( dirname,'resources') && ...
             ~strcmp( dirname,'__archive')
          p = [p genpath([d filesep dirname])]; % recursive calling of this function.
       end
    end
end
%% ===================================================================== %%
function [args] = eeglab_struct2args(struct)
    %EEGLAB_STRUCT2ARGS Summary of this function goes here
    %   Detailed explanation goes here
    %## Define Parser
    p = inputParser;
    %## REQUIRED
    addRequired(p,'struct',@isstruct);
    parse(p,struct);
    %##
    fn = fieldnames(struct);
    sc = struct2cell(struct);
    args = cell(length(fn)+length(sc),1);
    cnt = 1;
    for i = 1:length(fn)
        args{cnt} = fn{i};
        args{cnt+1} = sc{i};
        cnt = cnt + 2;
    end
end

