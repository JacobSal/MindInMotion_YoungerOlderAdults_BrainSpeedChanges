%   Project Title: SET WORKSPACE FOR SCRIPTS
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220810.0
%   Previous Version: n/a
%   Summary: this script is an initializer and workspace variable setup for
%   all scripts in this repository

%## TIME
tic

%% FLEXIBLE HANDLING OF SRC FOLDER
% restoredefaultpath;
%## OVERRIDE USING "global"
% global FUNCTIONS_VERSION ADD_CLEANING_SUBMODS ADD_DIPFIT_COMPILE_SUBMODS
if ~exist('FUNCTIONS_VERSION','var')
    FUNCTIONS_VERSION = 'v2_0';
end
if ~exist('ADD_CLEANING_SUBMODS','var')
    ADD_CLEANING_SUBMODS = false;
end
if ~exist('ADD_DIPFIT_COMPILE_SUBMODS','var')
    ADD_DIPFIT_COMPILE_SUBMODS = false;
end
% if ~exist('ADD_SIFT','var')
%     ADD_SIFT = false;
% end
% Current Directory
tmp                 = dir(['.' filesep]);
fprintf(1,'Current folder: %s\n',tmp(1).folder);
%## datetime
dt                  = datetime;
dt.Format = 'ddMMyyyy';

% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% tmp = dir(cfpath);
% if ~exist([tmp(1).folder filesep '_functions'],'dir')
%     error('''_functions'' folder does not exist. Please make one');
% end
% if ~exist([tmp(1).folder filesep '_functions' filesep FUNCTIONS_VERSION],'dir')
%     error('A version CHAR is not assigned, or the wrong version CHAR is assigned. see. setWorkspace.m');
% end
% ----------------------------------------------------------------------- %
%% PARAMS
cfname_path    = mfilename('fullpath');
cfpath = strsplit(cfname_path,filesep);
if ~ispc
    %- Add source directory where setup functions are 
    path4src = [filesep strjoin(cfpath(1:end-1),filesep)];
else
    %- Add source directory where setup functions are
    path4src = strjoin(cfpath(1:end-1),filesep);
end
if ~ispc
    %- Add submodules directory where packages are 
    submodules_dir = [filesep strjoin(cfpath(1:end-2),filesep) filesep 'submodules'];
else
    %- Add submodules directory where packages are 
    submodules_dir = [strjoin(cfpath(1:end-2),filesep) filesep 'submodules'];
end
fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n',path4src,submodules_dir);
% ----------------------------------------------------------------------- %
%% HARDCODE PATHS STRUCT
PATHS = [];
%## ALL SUBMODS
%{
SUBMODULES = {'groupSIFT','SIFT','fieldtrip','spm12','eeglab','postAmicaUtility',...
                'bsmart','Granger_Geweke_Causality','MindInMotion','bemobil-pipeline-master','bids-matlab-tools5.3.1',...
                'Cleanline2.00','firfilt','ICLabel','LIMO3.2','PowPowCAT3.0','bva-io1.7',...
                'EEGLAB-specparam-master','iCanClean','Viewprops1.5.4',...
                'trimOutlier-master'};
%}
%- Cleaning SUBMODS
if ADD_CLEANING_SUBMODS
    SUBMODULES = {'eeglab','SIFT','fieldtrip','spm12','postAmicaUtility',...
                    'Granger_Geweke_Causality','MindInMotion','bemobil-pipeline-master','bids-matlab-tools5.3.1',...
                    'Cleanline2.00','firfilt','ICLabel','LIMO3.2','PowPowCAT3.0','bva-io1.7',...
                    'EEGLAB-specparam-master','iCanClean','Viewprops1.5.4',...
                    'trimOutlier-master'};
    SUBMODULES_GENPATH = {'Cleanline2.00'};
    SUBMODULES_ITERS = (1:length(SUBMODULES));
elseif ADD_DIPFIT_COMPILE_SUBMODS
    SUBMODULES = {'fieldtrip','eeglab','postAmicaUtility'};
    SUBMODULES_GENPATH = {};
    SUBMODULES_ITERS = (1:length(SUBMODULES));
else
    %- Conn SUBMODS
    SUBMODULES = {'fieldtrip','eeglab','SIFT','postAmicaUtility',...
        'Granger_Geweke_Causality',...
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
    PATHS.PATHS{ss} = [submodules_dir filesep SUBMODULES{ss}];
end
%## special paths
%- src folder
PATHS.path4src = path4src;
%- _data folder
PATHS.localStorage = fullfile(path4src,'_data');
%- EEGLAB folder
PATHS.path4EEGlab = [submodules_dir filesep 'eeglab'];
%- _functions folder
PATHS.path4functions = [PATHS.path4src filesep '_functions' filesep FUNCTIONS_VERSION];
a_ftmp = unix_genpath(PATHS.path4functions);
a_ftmp = split(a_ftmp,DELIM); a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
cellfun(@(x) path(path,x),a_ftmp);
cellfun(@(x) fprintf('Adding functions in: %s...\n',x),a_ftmp);
clear a_ftmp
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
%% Colors
color.lightGreen  = [186,228,179]/255;color.Green = [44 162 95]/255;
color.lightRed = [252 174 145]/255;color.Red = [222,45,38]/255;
color.lightBlue = [189 215 231]/255;color.Blue = [49 130 189]/255;
color.lightGray = [0.7 0.7 0.7];color.darkGray = [0.3 0.3 0.3];
color.RedGradient = [255,245,240; 254,224,210; 252,187,161 ; 252,146,114; 251,106,74;...
    239,101,72;252,141,89;253,187,132;239,59,44; 220,40,30;203,24,29;165,15,21 ;103,0,13]/255;
color.BlueGradient = [247,251,255;222,235,247;198,219,239;158,202,225;158,188,218;107,174,214;66,146,198;...
    33,113,181;8,81,156;8,48,107]/255;
color.GreenGradient = [247,252,245; 229,245,224; 217,240,163; 199,233,192; 173,221,142; 161,217,155; 116,196,118; 
    65,171,93; 35,132,67;
    35,139,69; 0,109,44; 0,68,27]/255;
color.White = [1 1 1];color.Black = [0 0 0];
color.subject = [141,211,199;255,255,179;190,186,218;251,128,114;128,177,211;253,180,98;179,222,105;252,205,229;217,217,217];
color.ptb = [66,146,198;158,188,218;198,219,239;199,233,192; 116,196,118; 35,139,69;170 170 170]/255;
color.ptb = [27, 3, 163;27, 133, 189;135, 206, 251;178, 203, 32;88, 170, 61;0, 107, 62;170 170 170]/255;
color.subject = [color.subject;175, 191, 192;132, 126, 137;214, 162, 173;160, 175, 132;74, 102, 112;];
color.Yellow = [254,196,79]/255;
color.Orange = [255,127,0]/255;
color.terrain = [color.darkGray;color.Green;color.Yellow;color.Orange;color.Red];
color.terrain_shade = [217,217,217;161,217,155;255,255,204;254,232,200;252,187,161]./255;
color.seg = [0 0 0;200 200 200;213,62,79;175,141,195;230,245,152;252,141,89;153,213,148;50,136,189]/255;
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);
color.speed = [29, 120, 116;90,180,172;	213, 163, 32;	144, 109, 36]/255;
color.speed_shade = [132, 162, 161;205,232,230;238, 212, 144;	232, 210, 164]/255;
%## CLEAR MEM

%## TIME
toc

%% SUBFUNCTIONS
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
%%
%     ftmp = dir(directory);
%     ftmp = ftmp(([ftmp.isdir] == 1));
%     a_ftmp = cell(length(ftmp),1);
%     for i = 1:length(ftmp)
%         if ~(strcmp(ftmp(i).name,'.') || strcmp(ftmp(i).name,'..'))
%             a_ftmp{i} = [ftmp(i).folder filesep ftmp(i).name]; 
%         end
%     end
%     a_ftmp = a_ftmp(~cellfun(@isempty,a_ftmp));
%     a_ftmp{end+1} = directory;
%     a_ftmp = a_ftmp(~contains(a_ftmp,'__archive'));
%     a_ftmp = a_ftmp(~contains(a_ftmp,'private'));