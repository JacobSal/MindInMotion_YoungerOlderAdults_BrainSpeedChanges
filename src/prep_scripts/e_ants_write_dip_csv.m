%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 

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
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace

%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca_speed');
subj_names = [SUBJ_PICS{:}];
% SUBJ_PICS = {{'H3046','H3047','H3073','H3077','H3092', ...
%     'NH3023','NH3025','NH3027',' NH3028', ...
%     'NH3051','NH3056','NH3071','NH3082','NH3123'}};
% subj_chars = [SUBJ_PICS{:}];
%%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams'; 
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
ica_data_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
for subj_i = 1:length(subj_chars)
% for subj_i = find(strcmp(subj_names,'H3113'))
    subj_name = subj_chars{subj_i};
    dipfit_fPath = [ica_data_dir filesep subj_name filesep 'head_model'];
    try
        dip_struct = par_load(dipfit_fPath,'dipfit_struct.mat');
        coords = zeros(length({dip_struct.dip}),3);
        chan = zeros(length({dip_struct.dip}),1);
        for i = 1:length({dip_struct.dip})
            if ~isempty(dip_struct(i).dip)
                coords(i,:) = dip_struct(i).dip.pos;
                chan(i,:) = i;
            end
        end
        coords(~all(coords,2),:) = [];
        chan(~all(chan,2),:) = [];
        %-
        x_coord = -coords(:,1); 
        y_coord = -coords(:,2);
        z_coord = coords(:,3);
        % (05/09/2024) JS, flipping x,y signs to convert from LPS space to
        % RAS spcae. ANTs uses LPS but fieldtrip uses RAS (acpc as origin
        % still).
        %-
        t_in = table(x_coord,y_coord,z_coord);
    %     csvwrite([dipfit_fPath filesep 'dip_pos.csv'],coords);
        writetable(t_in,[dipfit_fPath filesep 'dip_pos.csv']);
        t_in = table(chan,x_coord,y_coord,z_coord);
        par_save(t_in,dipfit_fPath,'dip_pos.mat');
    catch e
        fprintf('\nSubject %s had an error occur.\n',subj_name);
        fprintf('\n%s\n',getReport(e));
    end
        
end