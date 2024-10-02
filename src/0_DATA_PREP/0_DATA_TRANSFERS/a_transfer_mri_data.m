%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh

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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
%% ===================================================================== %%
%## PARAMS
%- datetime override 
% dt = '03232023_AS_Bishoy';
dt = 'test_ls_alg';
study_fName = sprintf('copy_study');
%- soft define
% save_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset' filesep '_studies' filesep sprintf('%s',dt)];
save_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset', filesep '_studies' filesep 'arkaprava_data'];
load_dir = [PATHS.src_dir filesep '_data' filesep 'MIM_dataset'];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% (FIND FILES) ======================================================== %%
%## (HELPER SCRIPT) TRANSFER ELEC_ALIGNED.MAT & VOL.MAT & VALIDATION FIGURES FROM M:\ TO R:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        folder_to = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01'];
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'vol.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_mesh.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_sens_1.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_plot_sens_2.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        file_from = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'ft_sourceplot.fig'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
    end
end
%% (HELPER SCRIPT) TRANSFER HEADMODEL & ELECTRODE_ALIGNED DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        fprintf('%s) Transfering files for headmodel for subject...\n',SUBJ_PICS{group_i}{subj_i})
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI'];
        if exist(folder_to,'file')
            delete(folder_to)
        end
        if ~exist(folder_to,'dir')
            mkdir(folder_to)
        end
        %## BARE NECESSARY REQUIRED
        %- CustomElectrodeLocations.txt
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\HeadScan
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.txt'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        %- *_masks_contr.nii.gz
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Segmentation\headreco
        file_from = dir([R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep  'm2m_*' filesep '*_masks_contr.nii.gz']);
        if ~isempty(file_from)
            try
                copyfile([file_from.folder filesep file_from.name],folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: *_masks_contr.nii.gz\n',SUBJ_PICS{group_i}{subj_i});
        end
%         file_from = dir([R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'm2m_*' filesep '*_MRI_acpc_rs.nii']);
%         if ~isempty(file_from)
%             copyfile([file_from.folder filesep file_from.name],folder_to)
%         else
%             fprintf('%s) File does not exist: *_MRI_acpc_rs.nii\n',SUBJ_PICS{group_i}{subj_i});
%         end
        %- mri_acpc.mat
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Processed_fiducials
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        %- mri_acpc.mat
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Processed_fiducials
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc_rs.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        %- ctf_fiducials.mat
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'ctf_fiducials.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
        end
        %## REGENERATED VOL & ELEC ALIGNED
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned_init.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
%         %## custom electrode locations and headmodels for SKULL_0.0042
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned_init.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'vol.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',SUBJ_PICS{group_i}{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},file_from);
%         end
        fprintf('Done.\n')
    end
end
%% (HELPER SCRIPT) TRANSFER MRI & Fiducial DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'MRI'];
%         if exist(folder_to,'dir')
%             rmdir(folder_to,'s')
%         end
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' ];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'HeadModel'];
%         delete(folder_to)
%         if exist(folder_to,'dir')
%             rmdir(folder_to)
%         end
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        %- custom electrode locations for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',SUBJ_PICS{group_i}{subj_i},file_from)
        end
        %- custom headmodel for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',SUBJ_PICS{group_i}{subj_i},file_from)
        end
%         file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Raw' filesep ''];
%         copyfile(file_from,folder_to);
        %- custom electrode locations from digital 3D scan
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.txt'];    
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc rostral-caudal reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc_rs.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to ctf reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'ctf_fiducials.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri segmentation using SimNIBS aligned to acpc/s coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep sprintf('%s_MRI_acpc_rs.nii',SUBJ_PICS{group_i}{subj_i})];
%         copyfile(tmp,folder_to);
        %-----------------------------------------------------------------%
        %{
        %- electrode coordinates from HEADSCAN aligned to MRI mesh
        tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'elec_aligned.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
                continue;
            end
        else
            fprintf('copying file %s.\n',tmp);
            copyfile(tmp,folder_to);
        end
        
        %- subject specific headmodel volume 
        tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'vol.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',SUBJ_PICS{group_i}{subj_i},tmp);
                continue;
            end
        else
            fprintf('copying file %s.\n',tmp);
            copyfile(tmp,folder_to);
        end
        %}
    end
end
%% (HELPER SCRIPT) TRANSFER CHANLOCS DATA FROM R:\ TO M:\
BAD_SUBJS = {}; %{'NH3004','NH3009'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        file_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan'];
        delete(folder_to)
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        copyfile(file_from,folder_to);
    end
end
%% (HELPER SCRIPT) TRANSFER FILES FROM TEH AMICA FOLDER TO THE CLEAN FOLDER FOR ALLS
%SUBJECTS
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];
subjs = [SUBJ_PICS{:}];
% dt = '07042023_OA_prep_verified';
% dt = '04182023_YA_N37_prep_verified';
dt = '05192023_YAN33_OAN79_prep_verified';
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
for subj_i = 1:length(subjs)
    folder_from = [save_dir filesep subjs{subj_i} filesep 'amica'];
    folder_to = [save_dir filesep subjs{subj_i} filesep 'clean'];
    if ~exist([folder_to filesep 'W'],'file')
        fprintf('Transfering folder for subject %s...\n',subjs{subj_i});
        transfer_folder(folder_from,folder_to,'.');
    end
end
%% (HELPER SCRIPT) TRANSFER MERGED EEG DATA FROM R:\ TO M:\
BAD_SUBJS = {'NH3004','NH3009'};
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data\';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        % R:\Ferris-Lab\share\MindInMotion\Data\H3039\EEG\Trials
        folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Trials'];
%         folder_from = [R_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Merged'];
%         folder_to = [M_MIND_IN_MOTION_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'Merged'];
%         transfer_folder(folder_from,folder_to,'*EEGandLSandIMU.set');
%         transfer_folder(folder_from,folder_to,'*EEGandLSandIMU.fdt');
        tmp = dir([folder_to filesep '*.set']);
        if (isempty(tmp) && exist(folder_from,'dir')) && ~any(strcmp(SUBJ_PICS{group_i}{subj_i},BAD_SUBJS))
            transfer_folder(folder_from,folder_to,'*.set');
            transfer_folder(folder_from,folder_to,'*.fdt');
        end
    end
end
%% HELPER SCRIPT
EMAIL_CHAR = 'jsalminen@ufl.edu';
LOOP_VAR = 1:length(fPaths);
avg_ref_pca_reduction = 1; % length({'EEG'})
for subj_i = LOOP_VAR
    %## RUN MAIN_FUNC
    amica_out_fPath = [save_dir filesep subjectNames{subj_i} filesep 'amica'];
    float_fPath = [save_dir filesep subjectNames{subj_i} filesep 'clean'];
    set_fName = sprintf('%s_cleanEEG.set',subjectNames{subj_i});
    if ~exist(amica_out_fPath,'dir')
        mkdir(amica_out_fPath);
    end
    if exist(float_fPath,'dir')
        %- find .set file
        tmp = dir([float_fPath filesep '*.set']);
        %- load EEG
        EEG = pop_loadset('filepath',float_fPath,'filename',tmp.name);
        %- set the .fdt (float) file
        tmp = split(tmp.name,'.');
        float_fName = [tmp{1} '.fdt'];
        %- create bash and param files
        [EEG,cmd_out] = mim_prep_hpg_amica(EEG,[float_fPath filesep float_fName],float_fPath,EMAIL_CHAR,avg_ref_pca_reduction);
        fprintf('%s\n',cmd_out{2});
    end
    if ~ispc
        system([cmd_out,''])
    else
        fprintf('run this on unix\n');
    end
end
%% (HELPER SCRIPT) 
%{
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];
subjs = [SUBJ_PICS{:}];
% dt = '07042023_OA_prep_verified';
% dt = '04182023_YA_N37_prep_verified';
dt = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
load_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\spca_analysis';
for subj_i = 1:length(subjs)
    folder_from = [load_dir filesep subjs{subj_i} filesep 'clean'];
    folder_to = [save_dir filesep subjs{subj_i} filesep 'spca_subsample'];
    fprintf('Transfering folder for subject %s...\n',subjs{subj_i});
    try
        % if dirOut doesn't exist, create it.
        if ~exist(folder_to,'dir')
            fprintf('Making directory: %s\n',folder_to);
            mkdir(folder_to);
        end
        % dont copy if float/set files already exist
        dirFrom = dir([folder_from filesep 'cond*.mat']);
        for fi = 1:length(dirFrom)
            copyfile([dirFrom(1).folder filesep dirFrom(fi).name],folder_to);
            delete([dirFrom(1).folder filesep dirFrom(fi).name])
        end
        dirFrom = dir([folder_from filesep 'chan*.jpg']);
        for fi = 1:length(dirFrom)
            copyfile([dirFrom(1).folder filesep dirFrom(fi).name],folder_to);
            delete([dirFrom(1).folder filesep dirFrom(fi).name])
        end
        copyfile([folder_from filesep 'gait_ersp_spca.mat'],folder_to);
        delete([folder_from filesep 'gait_ersp_spca.mat'])
    catch e
        fprintf('Error. on Subject %s:\n\n%s\n',subjs{subj_i},getReport(e));
    end
end
%}
%% (HELPER SCRIPT) TRANSFER SPCA RESULTS
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET];
subjs = [SUBJ_PICS{:}];
% dt = '07042023_OA_prep_verified';
% dt = '04182023_YA_N37_prep_verified';
% dt = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% load_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% save_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\spca_analysis';

load_dir = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\spca_analysis';
save_dir = 'R:\Ferris-Lab\jsalminen\Experiments_Funding\Experiment_6_MIM_OA\data_set_saves\01122024_spca_analysis';
for subj_i = 1:length(subjs)
    folder_from = [load_dir filesep subjs{subj_i} filesep 'spca_subsample'];
    folder_to = [save_dir filesep subjs{subj_i} filesep 'spca_subsample'];
    fprintf('Transfering folder for subject %s...\n',subjs{subj_i});
    try
        % if dirOut doesn't exist, create it.
        if ~exist(folder_to,'dir')
            fprintf('Making directory: %s\n',folder_to);
            mkdir(folder_to);
        end
        % dont copy if float/set files already exist
        dirFrom = dir([folder_from filesep 'cond*.mat']);
        for fi = 1:length(dirFrom)
            copyfile([dirFrom(1).folder filesep dirFrom(fi).name],folder_to);
            delete([dirFrom(1).folder filesep dirFrom(fi).name])
        end
        dirFrom = dir([folder_from filesep 'chan*.jpg']);
        for fi = 1:length(dirFrom)
            copyfile([dirFrom(1).folder filesep dirFrom(fi).name],folder_to);
            delete([dirFrom(1).folder filesep dirFrom(fi).name])
        end
        copyfile([folder_from filesep 'gait_ersp_spca.mat'],folder_to);
        delete([folder_from filesep 'gait_ersp_spca.mat'])
    catch e
        fprintf('Error. on Subject %s:\n\n%s\n',subjs{subj_i},getReport(e));
    end
end