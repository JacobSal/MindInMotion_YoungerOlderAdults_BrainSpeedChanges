%   Project Title: MIM OA & YA SPEED & KINETICS ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_sts_b_epoch_eeg_kin.sh

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
ADD_ALL_SUBMODS = true;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
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
    STUDY_DIR = fileparts(SCRIPT_DIR); % change this if in sub folder
    SRC_DIR = STUDY_DIR;
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
addpath(STUDY_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%% (DATASET INFORMATION) =============================================== %%
% [SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('test');
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
%## UNIT TEST
% SUBJ_PICS = {{'H1004','H1007'},{'H2013','H2018_FU'},{'H3120','NH3129'}};
% GROUP_NAMES = {'H1000''s','H2000''s','H3000''s'};
% SUBJ_ITERS = {[1,2],[1,2],[1,2]};
%% Store fNames and fPaths
subj_chars          = [SUBJ_PICS{:}];
%% (PROC 1)
% tmp_alleeg = cell(length(fPaths),1);
%## PATHING UPDATES
% path(unix_genpath([PATHS.submods_dir filesep 'Gait Tracking With x-IMU']),path);
path(unix_genpath([PATHS.submods_dir filesep 'gait_tracking_w_imu']),path);
% data_fpath = '/blue/ark007/ark007/Dan lab IMU data/imu_data_full_dataset';
data_fpath = 'M:\jsalminen\GitHub\MIND_IN_MOTION_PRJ\_data\MIM_dataset';
%##
% parfor subj_i = 1:length(subj_chars)
for subj_i = 39
    trial_fpath = [data_fpath filesep subj_chars{subj_i} filesep 'IMU' filesep 'Raw'];
    %## Find All Trials of Interest
    files_imu = dir([trial_fpath filesep '*.csv']);
    tmp_savedir = [save_dir filesep subj_chars{subj_i} filesep 'alignd_data'];
    if ~exist(tmp_savedir,'dir')
        mkdir(tmp_savedir);
    end
    tmp_biom = cell(size(files_imu,1),1);
    [imu_table,rec_start_struct] =  convert_imu_to_table([files_imu(1).folder]);
    for trial_i = 1:size(files_imu,1)
        %- Load trial
        nn = strsplit(files_imu(trial_i).name,'.');
        nn = strsplit(nn{1},'-');
        nn = nn{end};
        [tbl_out,out_struct] = imu_get_pos_coords_tbl(imu_table{trial_i},tmp_savedir,nn);
        % nn = strsplit(fileList(trial_i).name,'.')
        % fh = imu_valid_plots(EEG,tmp_savedir,sprintf('%s_%s_',subj_chars{subj_i},nn{1}));
        % close(fh);
        % tmp_biom.nbchan = length(tmp_biom.chanlocs);
        tmp_biom{trial_i} = tbl_out;
        % pop_eegplot( EEG, 1, 1, 1);
        % fprintf('percent data points nan: %0.2g\n',sum(logical(isnan(tmp_biom{trial_i}.data)),'all')/(size(tmp_biom{trial_i}.data,1)*size(tmp_biom{trial_i}.data,2)));
    end
    % tmpt = cat(1,imu_table{:});    
    % [tbl_out,out_struct] = imu_get_pos_coords_tbl(tmpt,tmp_savedir,'all');
    tmp_biom = cat(1,tmp_biom{:});
    par_save(tmp_biom,[tmp_savedir filesep 'imu_table_orient.mat']);
    writetable(tmp_biom,[tmp_savedir filesep 'imu_table_orient.xlsx']);

end