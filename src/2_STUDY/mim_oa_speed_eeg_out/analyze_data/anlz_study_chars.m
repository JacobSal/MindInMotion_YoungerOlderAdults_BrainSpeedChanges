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
    STUDY_DIR = fileparts(SCRIPT_DIR);
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%% Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
%- statistics & conditions
SPEED_VALS = {'0.25','0.5','0.75','1.0';
              '0p25','0p5','0p75','1p0'};
TERRAIN_VALS = {'flat','low','med','high'};
COLORS_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%##
TRIAL_TYPES = [SPEED_VALS(2,:),TERRAIN_VALS];
%- datset name
DATA_SET = 'MIM_dataset';
%- cluster directory for study
% study_dir_name = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%- study info
SUB_GROUP_FNAME = [];
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',study_dir_name) filesep 'cluster'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
%% ================================================================== %%
%## LOAD STUDY
cluster_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
%## LOAD STUDY
% if ~ispc
%     tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
%     STUDY = tmp.STUDY;
% else
%     tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
%     STUDY = tmp.STUDY;
% end
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
end
cl_struct = par_load(cluster_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
save_dir = [spec_data_dir filesep 'psd_calcs'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
POSS_CLUSTER_CHARS = {};
% this is a matrix of integers matching the cluster number for clustering K=i to the index in the POSS_CLUSTER_CHARS
% SUB_GROUP_FNAME = 'H3000';  %'H2000';
% SUB_GROUP_FNAME_REGEX = 'H3000''s';  %'H2000''s';
SUB_GROUP_FNAME = []; 
SUB_GROUP_FNAME_REGEX = [];
CLUSTER_CLIM_MATCH = [];
%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
speed_alleeg = cell(length(ALLEEG),2);
for i = 1:length(ALLEEG)
    ss = ALLEEG(i).subject;
    ind = speed_table.Var1==ss;
    chk1 = strcmp(ALLEEG(i).group,SUB_GROUP_FNAME_REGEX) || isempty(SUB_GROUP_FNAME_REGEX);
    if any(ind) && chk1
        speed_alleeg{i,1} = speed_table.Var1(ind);
        speed_alleeg{i,2} = double(speed_table.Var2(ind));
    end
end
speed_alleeg = speed_alleeg(~cellfun(@isempty,speed_alleeg(:,1)),:);
%% STUDY STATS
subj_inds = 1:length(ALLEEG);
subj_chars = {STUDY.datasetinfo.subject}';
trial_counts = cell(length(subj_chars),length(TRIAL_TYPES));
% icc_noise_comps_rmv = cell(length(subj_chars),1);
channel_rejs = cell(length(subj_chars),1);
sample_rej_perc = cell(length(subj_chars),1);
for subj_i=1:length(subj_inds)
    EEG = ALLEEG(subj_i);
    for cond_i = 1:length(TRIAL_TYPES)
        tmp = sum(strcmp({STUDY.datasetinfo(subj_i).trialinfo.cond},TRIAL_TYPES{cond_i}));
        fprintf('%s) Number of %s trials: %i\n',EEG.subject,TRIAL_TYPES{cond_i},tmp);
        trial_counts{subj_i,cond_i} = tmp;
    end
%     icc_noise_comps_rmv{subj_i,1} = EEG.etc.iCanClean.numNoiseCompsRemovedOnAvg;
    channel_rejs{subj_i,1} = sum(~EEG.etc.channel_mask(strcmp({EEG.urchanlocs.type},'EEG')));
    sample_rej_perc{subj_i,1} = (length(EEG.etc.clean_sample_mask)-sum(EEG.etc.clean_sample_mask))/length(EEG.etc.clean_sample_mask);
end
tbl_out = table(subj_chars,trial_counts(:,1),trial_counts(:,2),trial_counts(:,3),...
    trial_counts(:,4),trial_counts(:,5),trial_counts(:,6),trial_counts(:,7),trial_counts(:,8),channel_rejs,sample_rej_perc,'VariableNames',{'subj_chars','0p25','0p5','0p75','1p0','flat','low','med','high','channel_rejs','sample_rej_perc'});
writetable(tbl_out,[cluster_study_fpath filesep 'study_characterisics.xlsx'])
%%
CLUSTER_STUDY_NAME = 
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
end
TRIALS_PROCESS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
trial_times = zeros(length(ALLEEG),length(TRIALS_PROCESS));
for subj_i = 1:length(ALLEEG)
    for cond_i = 1:length(TRIALS_PROCESS)
        fprintf('Running condition %s...\n',TRIALS_PROCESS{cond_i})
        EEG = ALLEEG(subj_i);
        %## Get Trial Indicies & Extract
        inds1 = logical(strcmp({EEG.event.cond}, TRIALS_PROCESS{cond_i}));
        inds2 = logical(strcmp({EEG.event.type}, 'boundary'));
        val_inds = find(inds1 & ~inds2);
        %-
        tmp_epochn = unique([tmp_event.epoch]);
        tmp_event = EEG.event(val_inds);
        lats = [tmp_event.latency];
        %-
        % tmp_epochn = unique([tmp_event.epoch]);
        % lats = zeros(length(tmp_epochn),1);
        % for epoch_i = 1:length(tmp_epochn)
        %     inds = [tmp_event.epoch] == tmp_epochn(epoch_i);
        %     tmp = tmp_event(inds);
        %     % lat = (max([tmp.latency])-min([tmp.latency]))/EEG.srate;
        %     lats(epoch_i) = (max([tmp.latency])-min([tmp.latency]))/EEG.srate;
        % end
        % trial_times(subj_i,cond_i) = mean(lats)*length(tmp_epochn);
        %-
        % fromt = [EEG.event(val_inds(1)).latency];
        % ind_from = find(EEG.times>=(((fromt-1)/EEG.srate)*1000) & EEG.times<=(((fromt+1)/EEG.srate)*1000));
        % tot = [EEG.event(val_inds(end)).latency];
        % ind_to = find(EEG.times>=(((tot-1)/EEG.srate)*1000) & EEG.times<=(((tot+1)/EEG.srate)*1000));
        % fprintf('%s) %s length is %0.2fs\n',EEG.subject,TRIALS_PROCESS{cond_i},(tot/EEG.srate)-(fromt/EEG.srate));
        % trial_times(subj_i,cond_i) = (tot/EEG.srate)-(fromt/EEG.srate);

    end
end
