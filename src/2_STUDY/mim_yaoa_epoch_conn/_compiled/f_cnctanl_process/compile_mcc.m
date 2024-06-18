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
%- run exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/_compiled/a_epoch_process/run_mcc_epoch_process.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    STUDY_DIR = getenv('STUDY_DIR');
    SCRIPT_DIR = getenv('SCRIPT_DIR');
    SRC_DIR = getenv('SRC_DIR');
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = filesep(filesep(STUDY_DIR));
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
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
    [PATHS.submods_dir filesep 'fieldtrip/external/simbio',...
    unix_genpath([PATHS.submods_dir filesep 'eeglab/plugins']),...
    };
cellfun(@(x) path(path,x),dep);
%% COMPILE
%- these are files that you should include if the compiler has trouble
%finding dependencies simply based on pathing. 
files_to_compile={[SCRIPT_DIR filesep 'mcc_epoch_process.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_create_alleeg.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_create_study.m'],...
    [PATHS.functions_dir filesep 'MIM/mim_parse_trials.m'],...
    [PATHS.submods_dir filesep 'eeglab/eeglab.m'],...
    [PATHS.submods_dir filesep 'eeglab/functions/popfunc/pop_reref.m']};
%- files that need to be included for later recall (e.g., .mat and .txt files)
data_to_include={[SCRIPT_DIR filesep 'eeg_options.txt'],...
    [SCRIPT_DIR filesep 'eeg_optionsbackup.txt'],...
    [PATHS.submods_dir filesep 'eeglab/plugins/ICLabel/netICL_lite.mat'],...
    [PATHS.submods_dir filesep 'eeglab/plugins/ICLabel/netICL.mat']};
fprintf('Compiling...\n');
%- cd to source directory
mkdir([SCRIPT_DIR filesep '_out']);
eval(['mcc -m ' strjoin(files_to_compile,' ')...
    ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' strjoin(data_to_include,' -a ')...
    ' -R -singleCompThread -v'])
%## TOC
toc