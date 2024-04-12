%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 
%- compile exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/_compiled/a_epoch_process/run_compile_mcc.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% REQUIRED SETUP 4 ALL SCRIPTS ======================================== %%
STUDY_DIR = getenv('STUDY_DIR');
SCRIPT_DIR = getenv('SCRIPT_DIR');
addpath(STUDY_DIR);
cd(SCRIPT_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%% SET WORKSPACE ======================================================= %%'
%## 
tic
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%- set EEGLAB path, _functions path, and others...
set_workspace
%##
com.mathworks.mwswing.MJUtilities.initJIDE;  % Initialize JIDE's usage within Matlab
pop_editoptions('option_storedisk',1,'option_savetwofiles',1, ...
    'option_single',1,'option_memmapdata',0,'option_computeica',0,...
    'option_saveversion6',1,'option_scaleicarms',1,'option_rememberfolder',1);
%% (HELPER CODE) DEFINE DEPENDENCIES
fprintf('Adding necessary paths to MATLAB path\n');
dep = {SCRIPT_DIR,...
    [PATHS.submods_dir filesep 'fieldtrip/contrib/spike'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/fileexchange'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/signal'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/stats'],...
    [PATHS.submods_dir filesep 'fieldtrip/fileio'],...
    [PATHS.submods_dir filesep 'fieldtrip/forward'],...
    [PATHS.submods_dir filesep 'fieldtrip/inverse'],...
    [PATHS.submods_dir filesep 'fieldtrip/plotting'],...
    [PATHS.submods_dir filesep 'fieldtrip/preproc'],...
    [PATHS.submods_dir filesep 'fieldtrip/utilities'],...
    [PATHS.submods_dir filesep 'postAmicaUtility'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/freesurfer'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/simbio']};
cellfun(@(x) path(path,x),dep);
%% COMPILE
files_to_compile={[SCRIPT_DIR filesep 'mcc_epoch_process.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_create_alleeg.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_create_study.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_parse_trials.m'],...
    };
data_to_include={'./eeg_options.txt','./eeg_optionsbackup.txt'};
fprintf('Compiling...\n');
%- cd to source directory
mkdir([SCRIPT_DIR filesep '_out']);
% mcc('-m',files_to_compile{:},'-d','./_out','-v','-R','-singleCompThread')
% mcc('-m',files_to_compile{:},'-d',[SCRIPT_DIR filesep '_out'],'-v','-R','-singleCompThread','-a',data_to_include{:})
mcc('-m',files_to_compile{:},'-d',[SCRIPT_DIR filesep '_out'],'-v','-R','-singleCompThread');
%## TOC
toc
exit();