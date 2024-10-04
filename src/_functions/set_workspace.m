%   Project Title: SET WORKSPACE FOR SCRIPTS
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220810.0
%   Previous Version: n/a
%   Summary: this script is an initializer and workspace variable setup for
%   all scripts in this repository

%% FLEXIBLE HANDLING OF SRC FOLDER
%## TIME
TT = tic;
%## folders
if ~exist('ADD_CLEANING_SUBMODS','var')
    ADD_CLEANING_SUBMODS = false;
end
if ~exist('ADD_DIPFIT_COMPILE_SUBMODS','var')
    ADD_DIPFIT_COMPILE_SUBMODS = false;
end
TMP_PWD = dir(['.' filesep]);
TMP_PWD = TMP_PWD(1).folder;
fprintf(1,'Current folder: %s\n',TMP_PWD);
%## datetime
dt         = datetime;
dt.Format  = 'ddMMyyyy';
% ----------------------------------------------------------------------- %
%% PARAMS
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
FUNC_FPATH = [src_dir filesep '_functions'];
%##
path(path,src_dir)
path(path,FUNC_FPATH);
fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n-FUNCTIONS: %s\n',src_dir,submodules_dir,FUNC_FPATH);
% ----------------------------------------------------------------------- %
%% HARDCODE PATHS STRUCT
PATHS = [];
%- Cleaning SUBMODS
%#ok<*UNRCH>
% if ADD_CLEANING_SUBMODS
%     SUBMODULES = {'eeglab','sift','fieldtrip','spm12','postamicautility',...
%                     'bemobil_pipeline','bidstool',...
%                     'cleanline','firfilt','iclabel','limo_tools','powpowcat','bvatool',...
%                     'eeglab_specparam','viewprops','bcilab','clean_rawdata','dipfit',...
%                     'trimoutlier','icanclean'}; 
%     SUBMODULES_GENPATH = {'cleanline'};
%     SUBMODULES_ITERS = (1:length(SUBMODULES));
% elseif ADD_DIPFIT_COMPILE_SUBMODS
%     SUBMODULES = {'fieldtrip','eeglab','postamicautility'};
%     SUBMODULES_GENPATH = {};
%     SUBMODULES_ITERS = (1:length(SUBMODULES));
% else
%     %- Conn SUBMODS
%     SUBMODULES = {'fieldtrip','eeglab','sift','postamicautility',...
%         'iclabel','viewprops','powpowcat'};
%     SUBMODULES_GENPATH = {};
%     SUBMODULES_ITERS = (1:length(SUBMODULES));
% end
if ADD_CLEANING_SUBMODS
    SUBMODULES = {'eeglab','SIFT','fieldtrip','spm12','postAmicaUtility',...
                    'Granger_Geweke_Causality','MindInMotion','bemobil-pipeline-master','bids-matlab-tools5.3.1',...
                    'Cleanline2.00','firfilt','ICLabel','LIMO3.2','PowPowCAT3.0','bva-io1.7',...
                    'EEGLAB-specparam-master','iCanClean','Viewprops1.5.4',...
                    'trimOutlier-master','Gait Tracking With x-IMU'};
    SUBMODULES_GENPATH = {'Cleanline2.00','Gait Tracking With x-IMU'};
    SUBMODULES_ITERS = (1:length(SUBMODULES));
elseif ADD_DIPFIT_COMPILE_SUBMODS
    SUBMODULES = {'fieldtrip','eeglab','postAmicaUtility'};
    SUBMODULES_GENPATH = {};
    SUBMODULES_ITERS = (1:length(SUBMODULES));
else
    %- Conn SUBMODS
    SUBMODULES = {'fieldtrip','eeglab','SIFT','postAmicaUtility',...
        'Granger_Geweke_Causality','AAL3'...
        'ICLabel','Viewprops1.5.4','PowPowCAT3.0'};
    SUBMODULES_GENPATH = {};
    SUBMODULES_ITERS = (1:length(SUBMODULES));
end
%## add submodules
if ispc
    DELIM = ';';
else
    DELIM = ':';
end
PATHS.PATHS = cell(length(SUBMODULES),1);
for ss = SUBMODULES_ITERS
    if any(strcmp(SUBMODULES{ss},SUBMODULES_GENPATH))
        fprintf('Adding submodule using genpath(): %s...\n',[submodules_dir filesep SUBMODULES{ss}]);
        a_ftmp = unix_genpath([submodules_dir filesep SUBMODULES{ss}]);
        a_ftmp = split(a_ftmp,DELIM); a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
        cellfun(@(x) path(path,x),a_ftmp);
        cellfun(@(x) fprintf('Adding functions in: %s...\n',x),a_ftmp);
    else
        fprintf('Adding submodule: %s...\n',[submodules_dir filesep SUBMODULES{ss}]);
        path(path,[submodules_dir filesep SUBMODULES{ss}]);
    end
    %- spm12 exception
    if strcmp(SUBMODULES{ss},'spm12')
        path(path,[submodules_dir filesep SUBMODULES{ss} filesep SUBMODULES{ss}]);
    end
    PATHS.PATHS{ss} = [submodules_dir filesep SUBMODULES{ss}];
end
%## special paths
% %- current dir
% PATHS.curr_dir = TMP_PWD;
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
a_ftmp = unix_genpath(PATHS.functions_dir);
a_ftmp = split(a_ftmp,DELIM); a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
cellfun(@(x) path(path,x),a_ftmp);
cellfun(@(x) fprintf('Adding functions in: %s...\n',x),a_ftmp);
a_ftmp = char.empty;
% ----------------------------------------------------------------------- %
%% ADDPATH for FIELDTRIP Toolboxbemobil
if contains('fieldtrip',SUBMODULES)
    ft_defaults;
end
%% INITIALIZE MIM & EEGLAB
%start EEGLAB if necessary
if contains('SIFT',SUBMODULES)
    StartSIFT;
end
%- always start eeglab last.
ALLEEG=[]; STUDY=[]; CURRENTSET=0; CURRENTSTUDY=0;
eeglab;

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
    if exist(STUDY_DIR,'var')
        mkdir([STUDY_DIR filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')])
        pp.JobStorageLocation = [STUDY_DIR filesep '_slurm_scratch' filesep getenv('SLURM_JOB_ID')];
    else
        fprintf('Slurm scratch location is not defined...\n\n');
    end
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', Inf);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0,'option_saveversion6',1, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%## TIME
fprintf('Done with workplace setup: %0.2f\n\n',toc(TT));
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
             ~strcmp( dirname,'__archive') && ...
             ~strcmp( dirname,'_compiles')
          p = [p genpath([d filesep dirname])]; % recursive calling of this function.
       end
    end
end