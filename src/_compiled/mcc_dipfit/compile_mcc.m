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
%%
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
    unix_genpath([PATHS.submods_dir filesep 'eeglab/plugins']),...
    };
%(12/05/2024) JS, removed "[PATHS.submods_dir filesep
%'fieldtrip/external/simbio']" as I think its messing with parallel
%processing
cellfun(@(x) path(path,x),dep);
%% COMPILE
%- these are files that you should include if the compiler has trouble
%finding dependencies simply based on pathing. 
%- method 1
fout = matlab.codetools.requiredFilesAndProducts([SCRIPT_DIR filesep 'mcc_dipfit.m']);
par_save(fout,SCRIPT_DIR,'dependencies.mat');
files_to_compile = par_load(SCRIPT_DIR,'dependencies.mat');
keepinds = zeros(length(files_to_compile),1);
paths_to_add = {};
data_to_include = {};
% rejections = {'horzcat2.m'};
extra = {};
for i = 1:length(files_to_compile)
    if ~ispc
        files_to_compile{i} = convertPath2UNIX(files_to_compile{i});
    else
        files_to_compile{i} = convertPath2Drive(files_to_compile{i});
    end
    tmp = strsplit(files_to_compile{i},filesep);
    if any(contains(tmp,'@')) || any(contains(tmp,'+')) || any(contains(tmp,'private'))
        fprintf('Skipping: %s\n',strjoin(tmp,filesep));
    else
        paths_to_add = [paths_to_add, strjoin(tmp(1:end-1),filesep)];
        regchk = regexp(files_to_compile{i},'(?<=\.)[^.]*$','match');
        if strcmp(regchk,'m')
            keepinds(i) = 1;
        elseif strcmp(regchk,'mat')
            data_to_include = [data_to_include, files_to_compile(i)];
        elseif strcmp(regchk,'mexw64')
            if ~ispc
                tmp = strsplit(files_to_compile{i},'.');
                tmp{end}= 'mexa64';
                tmp = strjoin(tmp,'.');
                files_to_compile{i} = tmp;
            end
            disp(i);
            keepinds(i) = 1;
        end
    end
end
%## FINAL CHECK FOR DUPLICATES
files_to_compile = files_to_compile(logical(keepinds));
out = cellfun(@(x) strsplit(x,filesep),files_to_compile,'UniformOutput',false);
out_n = cellfun(@(x) regexp(x{end},'[^.]*(?=\.)','match'),out);
out_e = cellfun(@(x) regexp(x{end},'(?<=\.)[^.]*$','match'),out);
[~,ind,ic] = unique(out_n);
ind_keep = [];
for o = 1:length(ind)
    tmp = out{ind(o)};
    t_n = out_n{ind(o)};
    t_e = out_e{ind(o)};
    ii = find(strcmp(t_n,out_n));
    if length(ii) > 1
        t_e = out_e(ii);
        i1 = strcmp(t_e,'m');
        if ~ispc
            i2 = strcmp(t_e,'mexa64');
        else
            i2 = strcmp(t_e,'mexw64');
        end
        if any(i2)
            ind_keep = [ind_keep, ii(i2)];
        else
            fprintf('Unable to reconcile multiple file matches for: %s\n',t_n);
            fprintf('Values: '); fprintf('%s, ',t_e{:}); fprintf('\n');
        end
    else
        ind_keep = [ind_keep, ind(o)];
    end
end
filese_to_compile = files_to_compile(ind_keep);
% files_to_compile = {[SCRIPT_DIR filesep 'mcc_epoch_process.m'],...
%     files_to_compile{:},...
%     [PATHS.submods_dir filesep 'eeglab/plugins/ICLabel/matconvnet/matlab/+dagnn_bc/Conv.m'],...
%     };
paths_to_add = [unique(paths_to_add)];
cellfun(@(x) path(path,x),paths_to_add);
% disp(files_to_compile);
% disp(data_to_include);
%- method 2
% iclabel_pckg = strsplit(mcc_gen_filepaths([PATHS.submods_dir filesep 'eeglab/plugins/ICLabel']),pathsep);
% files_to_compile={[SCRIPT_DIR filesep 'mcc_epoch_process.m'],...
%     [PATHS.functions_dir filesep 'MIM/mim_create_alleeg.m'],...
%     [PATHS.functions_dir filesep 'MIM/mim_create_study.m'],...
%     [PATHS.functions_dir filesep 'MIM/mim_parse_trials.m'],...
%     [PATHS.submods_dir filesep 'eeglab/eeglab.m'],...
%     [PATHS.submods_dir filesep 'eeglab/functions/popfunc/pop_reref.m'],...
%     iclabel_pckg{:}};
%     % [PATHS.submods_dir filesep 'eeglab/plugins/ICLabel/matconvnet/matlab/+dagnn_bc/Conv.m']};
%- files that need to be included for later recall (e.g., .mat and .txt files)
data_to_include={data_to_include{:}...
    [SCRIPT_DIR filesep 'eeg_options.txt'],...
    [SCRIPT_DIR filesep 'eeg_optionsbackup.txt']};
%% ===================================================================== %%
fprintf('Compiling...\n');
%- cd to source directory
mkdir([SCRIPT_DIR filesep '_out']);
eval(['mcc -m ' strjoin(files_to_compile,' ')...
    ' -d ' [SCRIPT_DIR filesep '_out'] ' -a ' strjoin(data_to_include,' -a ')...
    ' -R -singleCompThread -v'])
%## TOC
toc