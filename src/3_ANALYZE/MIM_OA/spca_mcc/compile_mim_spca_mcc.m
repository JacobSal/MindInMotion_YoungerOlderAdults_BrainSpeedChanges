%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 
%- compile exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc/compile_mim_spca_mcc.sh
%- run test exe
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit/_out/run_mim_mcc_dipfit_exe.sh
%
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
%% REQUIRED SETUP 4 ALL SCRIPTS
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s/n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_OA' filesep 'spca_mcc'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_DIPFIT_COMPILE_SUBMODS
ADD_DIPFIT_COMPILE_SUBMODS = false;
setWorkspace
com.mathworks.mwswing.MJUtilities.initJIDE;  % Initialize JIDE's usage within Matlab
pop_editoptions('option_storedisk',1,'option_savetwofiles',1, ...
    'option_single',1,'option_memmapdata',0,'option_computeica',0,...
    'option_saveversion6',1,'option_scaleicarms',1,'option_rememberfolder',1);
%% (HELPER CODE) DEFINE DEPENDENCIES
% dep = matlab.codetools.requiredFilesAndProducts('calc_stiff_matrix_val.m');
% dep = cellfun(@fileparts,dep);
% dep = unique(dep);
%     '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/adminfunc',...
%     '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/guifunc',...
%     '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/popfunc',...
%     '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/sigprocfunc',...
%     '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/studyfunc',...
fprintf('Adding necessary paths to MATLAB path\n');
dep = {'/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/i_HEADMODEL/2_dipole_fit/MIM/mcc_dipfit',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/contrib/spike',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/_override/simbio',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/fileexchange',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/signal',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/stats',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/fileio',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/forward',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/inverse',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/plotting',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/preproc',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/utilities',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/postAmicaUtility',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/freesurfer',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/simbio'};
cellfun(@(x) path(path,x),dep);
%% COMPILE
files_to_compile={'/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/spca_mcc/mim_spca_mcc.m',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/sPCA/plot_txf.m',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/sPCA/apply_spca_cond_timewarp.m',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_functions/v2_0/eeglab/eeglab_pop_subcomp.m',...
    '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/eeglab/functions/studyfunc/pop_precomp.m'};
fprintf('Compiling...\n');
%- cd to source directory
cd(run_dir)
% mcc('-m','./mim_mcc_dipfit.m','-a',dep{:},'-d','./_out','-v')
% mcc('-m','./mim_mcc_dipfit.m','-a','/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip/external/simbio/*','./_out','-v')
mcc('-m',files_to_compile{:},'-d','./_out','-v')
%## TOC
toc
exit();
%%
%{
ica_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
study_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\test_spca'
jf_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\3_ANALYZE\MIM_OA\spca_mcc\spca_mcc_params.json';
addpath('M:\jsalminen\GitHub\par_EEGProcessing\src\3_ANALYZE\MIM_OA\spca_mcc');
[err] = mim_spca_mcc(ica_dir,study_dir,jf_fpath,'SUBJ_CHARS',{'H1004'})
%}