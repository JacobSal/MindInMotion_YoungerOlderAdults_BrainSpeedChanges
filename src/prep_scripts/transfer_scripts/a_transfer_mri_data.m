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
ADD_ALL_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR; % change this if in sub folder
        SRC_DIR = STUDY_DIR;
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
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
SUBJ_PICS = {{'H3046','H3047','H3073','H3077','H3092', ...
    'NH3023','NH3025','NH3027',' NH3028', ...
    'NH3051','NH3056','NH3071','NH3082','NH3123'}};
%% ===================================================================== %%
%## PARAMS
%- datetime override
study_fName = sprintf('copy_study');
%- soft define
save_dir_bluedrive = [PATHS.data_dir filesep 'MIM_dataset', filesep '_studies' filesep study_fName];
load_dir_bluedrive = [PATHS.data_dir filesep 'MIM_dataset'];
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = load_dir_bluedrive;
%- create new study directory
if ~exist(save_dir_bluedrive,'dir')
    mkdir(save_dir_bluedrive);
end
%% (TRANSFER MRI R->BLUE) ============================================== %%
%## MRI .NII ACPC_RS
subj_chars = [SUBJ_PICS{:}];
for subj_i = 1:length(subj_chars)
    file_from = dir([R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep, ...
        'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'm2m_*' filesep '*_MRI_acpc_rs.nii']);
    folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI'];
    %--
    % if exist(folder_to,'file')
    %     delete(folder_to)
    % end
    if ~exist(folder_to,'dir')
        mkdir(folder_to)
    end
    %-- transfer if exist
    if ~isempty(file_from)
        copyfile([file_from.folder filesep file_from.name],folder_to)
    else
        fprintf('%s) File does not exist: *_MRI_acpc_rs.nii\n',subj_chars{subj_i} );
    end
end
%## CUSTOM ELECTRODES .TXT
for subj_i = 1:length(subj_chars)
    %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\HeadScan
    file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep, ...
        'HeadScan' filesep 'CustomElectrodeLocations.txt'];
    folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI'];
    %--
    % if exist(folder_to,'file')
    %     delete(folder_to)
    % end
    if ~exist(folder_to,'dir')
        mkdir(folder_to)
    end
    %-- transfer if exist
    if exist(file_from,'file')
        try
            copyfile(file_from,folder_to)
        catch
            fprintf('file exists.\n')
        end
    else
        fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
    end
end
%% (HELPER SCRIPT) TRANSFER HEADMODEL & ELECTRODE_ALIGNED DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
        fprintf('%s) Transfering files for headmodel for subject...\n',subj_chars{subj_i})
        folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI'];
        if exist(folder_to,'file')
            delete(folder_to)
        end
        if ~exist(folder_to,'dir')
            mkdir(folder_to)
        end
        %## BARE NECESSARY REQUIRED
        %- CustomElectrodeLocations.txt
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\HeadScan
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.txt'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
        end
        %- *_masks_contr.nii.gz
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Segmentation\headreco
        file_from = dir([R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep  'm2m_*' filesep '*_masks_contr.nii.gz']);
        if ~isempty(file_from)
            try
                copyfile([file_from.folder filesep file_from.name],folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: *_masks_contr.nii.gz\n',subj_chars{subj_i});
        end
        % file_from = dir([R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'm2m_*' filesep '*_MRI_acpc_rs.nii']);
        % if ~isempty(file_from)
        %     copyfile([file_from.folder filesep file_from.name],folder_to)
        % else
        %     fprintf('%s) File does not exist: *_MRI_acpc_rs.nii\n',subj_chars{subj_i});
        % end
        %- mri_acpc.mat
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Processed_fiducials
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
        end
        %- mri_acpc.mat
        %R:\Ferris-Lab\share\MindInMotion\Data\NH3040\MRI\Processed_fiducials
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc_rs.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
        end
        %- ctf_fiducials.mat
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'ctf_fiducials.mat'];
        if exist(file_from,'file')
            try
                copyfile(file_from,folder_to)
            catch
                fprintf('file exists.\n')
            end
        else
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
        end
        %## REGENERATED VOL & ELEC ALIGNED
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned_init.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         %## custom electrode locations and headmodels for SKULL_0.0042
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned_init.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'elec_aligned.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.0042' filesep 'vol.mat'];
%         if exist(file_from,'file')
% %             fprintf('%s) copying file %s',subj_chars{subj_i},file_from);
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
        fprintf('Done.\n')
    end
end
% %## (HELPER SCRIPT) TRANSFER ELEC_ALIGNED.MAT & VOL.MAT & VALIDATION FIGURES FROM M:\ TO R:\
% R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
% M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
% %- Loop through directory
% for group_i = 1:size(SUBJ_PICS,2)
%     for subj_i = 1:length(SUBJ_PICS{group_i})
%         folder_to = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01'];
%         if ~exist(folder_to,'dir')
%             mkdir(folder_to);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'elec_aligned.mat'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'vol.mat'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'ft_plot_mesh.fig'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'ft_plot_sens_1.fig'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'ft_plot_sens_2.fig'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%         file_from = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'ft_sourceplot.fig'];
%         if exist(file_from,'file')
%             copyfile(file_from,folder_to)
%         else
%             fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},file_from);
%         end
%     end
% end
%% (HELPER SCRIPT) TRANSFER MRI & Fiducial DATA FROM R:\ TO M:\
R_MIND_IN_MOTION_DIR = 'R:\Ferris-Lab\share\MindInMotion\Data';
M_MIND_IN_MOTION_DIR = [DATA_DIR filesep DATA_SET]; %'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset'
%- Loop through directory
for group_i = 1:size(SUBJ_PICS,2)
    for subj_i = 1:length(SUBJ_PICS{group_i})
%         folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'EEG' filesep 'MRI'];
%         if exist(folder_to,'dir')
%             rmdir(folder_to,'s')
%         end
%         folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' ];
        folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'HeadModel'];
%         delete(folder_to)
%         if exist(folder_to,'dir')
%             rmdir(folder_to)
%         end
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        %- custom electrode locations for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',subj_chars{subj_i},file_from)
        end
        %- custom headmodel for SKULL_0.01
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
        if exist(file_from,'file')
            copyfile(file_from,folder_to)
        else
            fprintf('Subject %s does not have file: %s',subj_chars{subj_i},file_from)
        end
%         file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Raw' filesep ''];
%         copyfile(file_from,folder_to);
        %- custom electrode locations from digital 3D scan
%         tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.txt'];    
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc rostral-caudal reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc_rs.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to acpc reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'mri_acpc.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri aligned to ctf reference coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Processed_fiducials' filesep 'ctf_fiducials.mat' ];
%         if exist(tmp,'file')
%             fprintf('copying file %s.\n',tmp);
%             copyfile(tmp,folder_to);
%         end
        %- mri segmentation using SimNIBS aligned to acpc/s coordinate system
%         tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep sprintf('%s_MRI_acpc_rs.nii',subj_chars{subj_i})];
%         copyfile(tmp,folder_to);
        %-----------------------------------------------------------------%
        %{
        %- electrode coordinates from HEADSCAN aligned to MRI mesh
        tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'elec_aligned.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'elec_aligned.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},tmp);
                continue;
            end
        else
            fprintf('copying file %s.\n',tmp);
            copyfile(tmp,folder_to);
        end
        
        %- subject specific headmodel volume 
        tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Segmentation' filesep 'headreco' filesep 'vol.mat'];
        if ~exist(tmp,'file')
            fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},tmp);
            try
                tmp = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'MRI' filesep 'Headmodel' filesep 'skull_0.01' filesep 'vol.mat'];
                fprintf('copying file %s.\n',tmp);
                copyfile(tmp,folder_to);
            catch
                fprintf('%s) File does not exist: %s\n',subj_chars{subj_i},tmp);
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
        file_from = [R_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
        folder_to = [M_MIND_IN_MOTION_DIR filesep subj_chars{subj_i} filesep 'EEG' filesep 'HeadScan'];
        delete(folder_to)
        if ~exist(folder_to,'dir')
            mkdir(folder_to);
        end
        copyfile(file_from,folder_to);
    end
end