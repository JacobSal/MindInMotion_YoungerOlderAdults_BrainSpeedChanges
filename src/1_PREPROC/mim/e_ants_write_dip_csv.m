%   Project Title: MIM PREPROCESSING SCRIPTS
%
%   Code Designer: Jacob salminen
%   Summary: 


%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
global ADD_CLEANING_SUBMODS STUDY_DIR SCRIPT_DIR %#ok<GVMIS>
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
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa');
subj_names = [SUBJ_PICS{:}];
%%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
%- eeglab_cluster.m spectral params
% OA_PREP_FPATH = '05192023_YAN33_OAN79_prep_verified'; % JACOB,SAL(04/10/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p4_changparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p65_iccREMG0p3_newparams'; % JACOB,SAL(09/26/2023)
% OA_PREP_FPATH = '08202023_OAN82_iccRX0p60_iccREMG0p3_newparams'; 
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams'; 
% OA_PREP_FPATH = '01132024_antsnorm_iccREEG0p65_iccREMG0p4_skull0p0042';
%## soft define
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
for subj_i = 1:length(subj_names)
% for subj_i = find(strcmp(subj_names,'H3113'))
    subj_name = subj_names{subj_i};
    dipfit_fPath = [OUTSIDE_DATA_DIR filesep subj_name filesep 'head_model'];
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
        x_coord = coords(:,1); 
        y_coord = coords(:,2);
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