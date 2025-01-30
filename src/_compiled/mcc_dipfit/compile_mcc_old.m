%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)

%- compile exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/_compiled/mcc_dipfit/run_compile_mcc_dipfit.sh

%- run exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/1_PREPROC/mim/mcc_dipfit/c_run_mim_mcc_dipfit.sh

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
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
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
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (HELPER CODE) DEFINE DEPENDENCIES
fprintf('Adding necessary paths to MATLAB path\n');
dep = {SCRIPT_DIR,...
    [PATHS.submods_dir filesep 'fieldtrip/contrib/spike'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/fileexchange'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/signal'],...
    [PATHS.functions_dir filesep '_override/simbio'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/stats'],...
    [PATHS.submods_dir filesep 'fieldtrip/fileio'],...
    [PATHS.submods_dir filesep 'fieldtrip/forward'],...
    [PATHS.submods_dir filesep 'fieldtrip/inverse'],...
    [PATHS.submods_dir filesep 'fieldtrip/plotting'],...
    [PATHS.submods_dir filesep 'fieldtrip/preproc'],...
    [PATHS.submods_dir filesep 'fieldtrip/utilities'],...
    [PATHS.submods_dir filesep 'postAmicaUtility'],...
    [PATHS.submods_dir filesep 'fieldtrip/external/freesurfer'],...
    };
%(12/05/2024) JS, removed "[PATHS.submods_dir filesep
%'fieldtrip/external/simbio']" as I think its messing with parallel
%processing
cellfun(@(x) path(path,x),dep);
%% COMPILE
%- these are files that you should include if the compiler has trouble
%finding dependencies simply based on pathing. 
files_to_compile={[SCRIPT_DIR filesep 'mcc_dipfit.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_sourcemodel.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_leadfield.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_dipolefitting.m'],...
    [PATHS.submods_dir filesep 'eeglab/plugins/dipfit/eeglab2fieldtrip.m'],...
    [PATHS.submods_dir filesep 'postAmicaUtility/pop_loadmodout.m'],...
    [PATHS.submods_dir filesep 'fieldtrip/ft_prepare_headmodel.m']};
%- files that need to be included for later recall (e.g., .mat and .txt files)
data_to_include={[PATHS.submods_dir filesep 'fieldtrip/external/simbio/calc_stiff_matrix_val.mexa64']};
fprintf('Compiling...\n');
%- cd to source directory
mkdir([SCRIPT_DIR filesep '_out']);
if length(data_to_include) > 1
    eval(['mcc -m ' strjoin(files_to_compile,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' strjoin(data_to_include,' -a ')...
        ' -R -singleCompThread -v'])
else
    eval(['mcc -m ' strjoin(files_to_compile,' ')...
        ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' data_to_include{1},...
        ' -R -singleCompThread -v'])
end
%## TOC
toc