%   Project Title: SET WORKSPACE FOR SCRIPTS
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230110.0
%   Previous Version: n/a
%   Summary: this script is an initializer and workspace variable setup for
%   all scripts in this repository

%## TIME
tic

%% FLEXIBLE HANDLING OF SRC FOLDER
global SUBMODULES_DIR
% import mlreportgen.ppt.*
FUNCTIONS_VERSION   = 'v1_0';
ADD_SIFT            = true;
ADD_BSMART          = false;
ADD_GGC             = true;
ADD_GT              = true;
% Current Directory
tmp                 = dir(['.' filesep]);
fprintf(1,'Current folder: %s\n',tmp(1).folder);
%## datetime
dt                  = datetime;
dt.Format = 'ddMMyyyy';
% ----------------------------------------------------------------------- %
%% PARAMS
try
    if strncmp(computer,'PC',2)
        DO_UNIX = false;
    else
        DO_UNIX = true;
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end
cfname_path = mfilename('fullpath');
cfpath = strsplit(cfname_path,filesep);
if DO_UNIX
    %- Add source directory where setup functions are 
    path4src = [filesep strjoin(cfpath(1:end-1),filesep)];
else
    %- Add source directory where setup functions are
    path4src = strjoin(cfpath(1:end-1),filesep);
end
if ~exist('SUBMODULES_DIR','var')
    if DO_UNIX
        %- Add submodules directory where packages are 
        submodules_dir = [filesep strjoin(cfpath(1:end-2),filesep) filesep 'submodules'];
    else
        %- Add submodules directory where packages are 
        submodules_dir = [strjoin(cfpath(1:end-2),filesep) filesep 'submodules'];
    end
else
    submodules_dir = SUBMODULES_DIR;
end
fprintf(1,'Using pathing:\n-WORKSPACE: %s\n-SUBMODULES: %s\n',path4src,submodules_dir);
% ----------------------------------------------------------------------- %
%% HARDCODE PATHS STRUCT
PATHS = [];
PATHS.path4src = path4src;
PATHS.localStorage = fullfile(path4src,'_data');
% set groupSIFT toolbox path
PATHS.path4groupSIFT = [submodules_dir filesep 'groupSIFT'];
% set SIFT toolbox path
PATHS.path4SIFT = [submodules_dir filesep 'SIFT'];
% set FieldTrip toolbox path
PATHS.path4ft = [submodules_dir filesep 'fieldtrip'];
% set spm toolbox path
PATHS.path4spm = [submodules_dir filesep 'spm12'];
% set EEGlab path
PATHS.path4EEGlab = [submodules_dir filesep 'eeglab'];
% set postAmicaUtility path
PATHS.path4postAmicaUtility = [submodules_dir filesep 'postAmicaUtility'];
% set viewprops path
PATHS.path4viewprops = [submodules_dir filesep 'viewprops'];
% set bsmart path
PATHS.path4bsmart = [submodules_dir filesep 'bsmart'];
% set bsmart path
PATHS.path4GGC = [submodules_dir filesep 'Granger_Geweke_Causality'];
% set Gait-Tracking-With-x-IMU-master path
PATHS.path4GaitTracking = [submodules_dir filesep 'Gait Tracking With x-IMU'];
% ----------------------------------------------------------------------- %
%% ADDPATH for FIELDTRIP Toolbox
fprintf(1,'adding path: %s\n',PATHS.path4ft)
path(path,PATHS.path4ft)
ft_defaults;
%% ADDPATH for spm12
fprintf(1,'adding path: %s\n',PATHS.path4spm)
path(path,PATHS.path4spm)
%% INITIALIZE MIM & EEGLAB
%start EEGLAB if necessary\
path(path,PATHS.path4EEGlab)
% add pop_prop_extended.m
path(path,PATHS.path4viewprops);
% add postAmicaUtility
path(path,PATHS.path4postAmicaUtility);
% add SIFT Toolbox
if ADD_SIFT
    fprintf('adding path: %s\n',PATHS.path4SIFT)
    path(path,PATHS.path4SIFT)
    StartSIFT;
end

% add BSMART Toolbox
if ADD_BSMART
    fprintf('adding path: %s\n',PATHS.path4bsmart)
    path(path,PATHS.path4bsmart)
%     bsmart;
end

% add Granger_Geweke_Causality (Hualou Liang)
if ADD_GGC
    fprintf('adding path: %s\n',PATHS.path4bsmart)
    path(path,PATHS.path4GGC);
end

% add Gait-Tracking-With-x-IMU-master path
if ADD_GT
    fprintf('adding path: %s\n',PATHS.path4GaitTracking);
    path(path,genpath(PATHS.path4GaitTracking));
end
% ----------------------------------------------------------------------- %
%- Add _functions folder
PATHS.path4functions = [PATHS.path4src filesep '_functions' filesep FUNCTIONS_VERSION];
fprintf('adding path: %s\n',PATHS.path4functions);
path(path,PATHS.path4functions);
% ----------------------------------------------------------------------- %
%- always start eeglab last.
eeglab;
%% COLORS
color.lightGreen  = [186,228,179]/255;
color.Green = [44 162 95]/255;

color.lightRed = [252 174 145]/255;
color.Red = [222,45,38]/255;

color.lightBlue = [189 215 231]/255;
color.Blue = [49 130 189]/255;

color.lightGray = [0.7 0.7 0.7];
color.darkGray = [0.3 0.3 0.3];

color.RedGradient = [255,245,240; 254,224,210; 252,187,161 ; 252,146,114; 251,106,74;...
    239,101,72;252,141,89;253,187,132;239,59,44; 220,40,30;203,24,29;165,15,21 ;103,0,13]/255;

color.BlueGradient = [247,251,255;222,235,247;198,219,239;158,202,225;158,188,218;107,174,214;66,146,198;...
    33,113,181;8,81,156;8,48,107]/255;
color.GreenGradient = [247,252,245; 229,245,224; 217,240,163; 199,233,192; 173,221,142; 161,217,155; 116,196,118; 
    65,171,93; 35,132,67;
    35,139,69; 0,109,44; 0,68,27]/255;

color.White = [1 1 1];
color.Black = [0 0 0];

color.subject = [141,211,199;255,255,179;190,186,218;251,128,114;128,177,211;253,180,98;179,222,105;252,205,229;217,217,217];
color.ptb = [66,146,198;158,188,218;198,219,239;199,233,192; 116,196,118; 35,139,69;170 170 170]/255;
color.ptb = [27, 3, 163;27, 133, 189;135, 206, 251;178, 203, 32;88, 170, 61;0, 107, 62;170 170 170]/255;
color.subject = [color.subject;175, 191, 192;132, 126, 137;214, 162, 173;160, 175, 132;74, 102, 112;];
color.Yellow = [254,196,79]/255;
color.Orange = [255,127,0]/255;
color.terrain = [color.darkGray;color.Green;color.Yellow;color.Orange;color.Red];
color.terrain_shade = [217,217,217;229,245,224;255,255,204;254,232,200;252,187,161]./255;
color.seg = [0 0 0;200 200 200;213,62,79;175,141,195;230,245,152;252,141,89;153,213,148;50,136,189]/255;

%## TIME
toc