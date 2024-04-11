%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_a_epoch_process.sh

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
%% REQUIRED SETUP 4 ALL SCRIPTS ======================================== %%
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current user: %s\n',USER_NAME);
PWD_DIR = mfilename('fullpath');
if contains(PWD_DIR,'LiveEditorEvaluationHelper')
    PWD_DIR = matlab.desktop.editor.getActiveFilename;
    PWD_DIR = fileparts(PWD_DIR);
else
    try
        PWD_DIR = matlab.desktop.editor.getActiveFilename;
        PWD_DIR = fileparts(PWD_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        PWD_DIR = dir(['.' filesep]);
        PWD_DIR = PWD_DIR(1).folder;
    end
end
addpath(PWD_DIR);
cd(PWD_DIR)
fprintf(1,'Current folder: %s\n',PWD_DIR);
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
ADD_CLEANING_SUBMODS = false;
%- set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
ica_orig_dir = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
% cluster_study_dir = '12012023_OAYA104_icc0p65-0p2_changparams';
cluster_study_dir = '01232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3';
% cluster_study_dir = '01232023_MIM_YAN32_antsnormalize_iccREMG0p4_powpow0p3_conn';

study_fName_1 = 'epoch_study';
spca_study_dir = '01122024_spca_analysis';
study_fName_2 = 'epoch_study';
ICLABEL_EYE_CUTOFF = 0.75;
%- study group and saving
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_1 = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_2 = [STUDIES_DIR filesep sprintf('%s',spca_study_dir)];
OUTSIDE_DATA_DIR = [STUDIES_DIR filesep ica_orig_dir]; % JACOB,SAL(02/23/2023)
%- load cluster
CLUSTER_DIR = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir) filesep 'cluster'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
CLUSTER_STUDY_DIR = [CLUSTER_DIR filesep 'icrej_5'];
CLUSTER_K = 12;
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    tmp = load('-mat',[CLUSTER_STUDY_DIR filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_FNAME)]);
    TMP_STUDY = tmp.STUDY;
else
    tmp = load('-mat',[CLUSTER_STUDY_DIR filesep sprintf('%s.study',CLUSTER_STUDY_FNAME)]);
    TMP_STUDY = tmp.STUDY;
end
cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
TMP_STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(TMP_STUDY);
condition_gait = unique({TMP_STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {TMP_STUDY.datasetinfo.subject};
%-
fPaths = {TMP_STUDY.datasetinfo.filepath};
fNames = {TMP_STUDY.datasetinfo.filename};
%%
CLUSTER_PICKS = main_cl_inds(2:end);
subj_i = 2;
cond_i = 1;
tmpd = dir([load_dir_2 filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep '*' filesep sprintf('cond%s_spca_ersp.mat',condition_gait{cond_i})]);
gait_epoch_subf = strsplit(tmpd.folder,filesep);
gait_epoch_subf = gait_epoch_subf{end};
spca_fpath = [load_dir_2 filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep gait_epoch_subf];
spca_ersp = par_load(spca_fpath,sprintf('cond%s_spca_ersp.mat',condition_gait{cond_i}));
T_DIM = size(spca_ersp.ersp_corr,1);
F_DIM = size(spca_ersp.ersp_corr,3);
CL_DIM = length(TMP_STUDY.cluster);
%%
loop_store = cell(length(subject_chars),1);
parfor (subj_i = 1:length(subject_chars),SLURM_POOL_SIZE)
% for subj_i = 1:length(subject_chars)
    %## INITIATION
    %-
    cond_c          = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    subj_c          = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    group_c         = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    cluster_c       = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    comp_c          = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_erspcorr_c   = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_gpmcorr_c    = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_ersporig_c   = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_gpmorig_c    = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_pc1_c        = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    tf_coeff_c      = cell(length(condition_gait)*length(CLUSTER_PICKS),1);
    cnt = 1;
    %## LOAD EEG DATA
    try
        %- get eye ic numbers
        fpath = [OUTSIDE_DATA_DIR filesep subject_chars{subj_i} filesep 'clean'];
        tmp = dir([fpath filesep '*.set']);
        [~,EEG_full,~] = eeglab_loadICA(tmp.name,tmp.folder);
        fprintf('Running subject %s\n',EEG_full.subject)
        EEG_full = iclabel(EEG_full);
        clssts = EEG_full.etc.ic_classification.ICLabel.classifications;
        bad_eye_ics = find(clssts(:,3) > ICLABEL_EYE_CUTOFF);
        %- figure out unmixing for cluster assignment
        ics_orig = 1:size(EEG_full.icaweights,2);
        tmp_cut = ics_orig;
        tmp_cut(bad_eye_ics) = [];
        [valc,ordc] = sort(tmp_cut);
        unmix = [valc; ordc];
        %- load spca data
        spca_fpath = [load_dir_2 filesep subject_chars{subj_i} filesep 'GAIT_EPOCHED' filesep gait_epoch_subf];
        EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
        fprintf('Running subject %s\n',EEG.subject)
        ic_keep = EEG.etc.urreject.ic_keep;
        ic_rej = EEG.etc.urreject.ic_rej;
        %- perhaps a new way of determining clusters
        tmp1 = cat(1,EEG.dipfit.model.posxyz);
        tmp2 = cat(1,EEG_full.dipfit.model.posxyz);
        dips = zeros(size(tmp2,1),1);
        for i = 1:size(tmp1,1)
            ind = find(all(tmp1(i,:)==tmp2,2));
            if ~isempty(ind)
                dips(i) = ind;
            end
        end
        dips = dips(dips~=0);
        %- sanity checks
        subj_ics = [];
        subj_cls = [];
        subj_unmix = [];
        for k = 2:length(TMP_STUDY.cluster)
%             cl_i = CLUSTER_PICKS(i)
            set_i = (TMP_STUDY.cluster(k).sets == subj_i);
            comp_i = TMP_STUDY.cluster(k).comps(set_i);
            if ~isempty(comp_i)
                subj_ics = [subj_ics, comp_i];
                subj_cls = [subj_cls, repmat(k,1,length(comp_i))];
                for j = 1:length(comp_i)
                    subj_unmix = [subj_unmix, unmix(2,ic_keep(comp_i(j)) == unmix(1,:))];
                end
            end
        end
        spca_ersp = par_load(spca_fpath,sprintf('cond%s_spca_ersp.mat',condition_gait{1}));
        chk = size(spca_ersp.ersp_corr,2) == (length([ic_keep,ic_rej])-length(bad_eye_ics));
        chk2 = (length(ic_keep)==length(TMP_STUDY.datasetinfo(subj_i).comps))
        fprintf('IC numbers check: %i\n',chk2&&chk);
        fprintf('IC''s in clusters:\n'); fprintf('ICs: '); fprintf('%i,',subj_ics); fprintf('\n');
        fprintf('Clusters:\t'); fprintf('%i,',subj_cls); fprintf('\n');
        fprintf('New ICs:\t'); fprintf('%i,',dips(subj_ics)); fprintf('\n');
        fprintf('Unmix ICs:\t'); fprintf('%i,',subj_unmix); fprintf('\n');
%         g_i = group_id(subj_i);
        for cond_i = 1:length(condition_gait)
            spca_ersp = par_load(spca_fpath,sprintf('cond%s_spca_ersp.mat',condition_gait{cond_i}));
            for i = 1:length(CLUSTER_PICKS)
                cl_i = CLUSTER_PICKS(i);
                set_i = (TMP_STUDY.cluster(cl_i).sets == subj_i);
                comp_i = TMP_STUDY.cluster(cl_i).comps(set_i);
                if ~isempty(comp_i) && length(comp_i) < 2
                    tmp_c = unmix(2,ic_keep(comp_i) == unmix(1,:));
                    fprintf('%s) assigning IC(spca index,study index) %i(%i,%i) to cluster %i in condition %s...\n',subject_chars{subj_i},ic_keep(comp_i),comp_i,tmp_c,cl_i,condition_gait{cond_i});
                    tf_erspcorr_c{cnt} = squeeze(spca_ersp.ersp_corr(:,tmp_c,:));
                    tf_gpmcorr_c{cnt} = squeeze(spca_ersp.gpm_corr(:,tmp_c,:));
                    tf_ersporig_c{cnt} = squeeze(spca_ersp.apply_spca_cond.erds_orig(:,tmp_c,:));
                    tf_gpmorig_c{cnt} = squeeze(spca_ersp.apply_spca_cond.gpm_orig(:,tmp_c,:));
                    tf_pc1_c{cnt} = squeeze(spca_ersp.apply_spca_cond.psc1(:,tmp_c,:));
                    tf_coeff_c{cnt} = spca_ersp.coeffs;
%                     tf_freqs{cnt} = [];
%                     tf_times{cnt} = [];
                    cond_c{cnt} = condition_gait{cond_i};
                    subj_c{cnt} = subject_chars{subj_i};
                    tmp = regexp(subject_chars{subj_i},'\d','match');
                    group_c{cnt} = str2num(tmp{1});
                    cluster_c{cnt} = cl_i;
                    comp_c{cnt} = sprintf('study_ic, %i; unmix_ic, %i; keep_ic, %i; cluster_ic, %i',comp_i,tmp_c,ic_keep(comp_i),comp_i);
                                %struct('eeglab_study',comp_i,...
%                                     'unmix',tmp_c,...
%                                     'ic_keep',ic_keep(comp_i));
                    cnt = cnt +1;
                end
            end
        end
        loop_store{subj_i} = table(subj_c,group_c,cluster_c,cond_c,comp_c,tf_erspcorr_c,tf_gpmcorr_c,tf_ersporig_c,tf_gpmorig_c,tf_pc1_c,tf_coeff_c);
    catch e
        fprintf('\nError occured on subject %s\n%s\n',subject_chars{subj_i},getReport(e));
    end
end
% tt = table(subj_c,group_c,cluster_c,cond_c,tf_erspcorr_c,tf_gpmcorr_c,tf_erporig_c,tf_gpmorig_c,tf_pc1_c,tf_coeff_c);
tt = cat(1,loop_store{:});
tt = tt(~all(cellfun(@isempty,tt{:,:}),2),:);
par_save(tt,CLUSTER_STUDY_DIR,'spca_cluster_table.mat');
