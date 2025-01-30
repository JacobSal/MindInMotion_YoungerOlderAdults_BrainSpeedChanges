%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/run_spca_f_psds_cluster.sh

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
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = SCRIPT_DIR;
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
    STUDY_DIR = SCRIPT_DIR;
    SRC_DIR = fileparts(fileparts(STUDY_DIR));
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SCRIPT_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
%- file regexps
ICA_REGEXP = '%s_cleanEEG_EMG_HP3std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10.set';
%- datset name
DATA_SET = 'MIM_dataset';
%- studies paths
ICA_DNAME = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams';
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
SPCA_STUDY_DNAME = '01102025_mim_yaoa_spca_calcs';
STUDY_FNAME_GAIT = 'spca_gait_epoch_study';
STUDY_FNAME_REST = 'spca_rest_slide_study';
ICLABEL_EYE_CUTOFF = 0.75;
%- study group and saving
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
spca_dir = [studies_fpath filesep sprintf('%s',SPCA_STUDY_DNAME)];
ica_data_dir = [studies_fpath filesep ICA_DNAME]; % JACOB,SAL(02/23/2023)
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb5']; % rb10 & rb3 available
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%% ===================================================================== %%
%## LOAD CLUSTER STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
%% ===================================================================== %%
%## LOAD SPCA STUDIES
%- gait
if ~ispc
    tmp = load('-mat',[spca_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
else
    tmp = load('-mat',[spca_dir filesep sprintf('%s.study',STUDY_FNAME_GAIT)]);
    STUDY_GAIT = tmp.STUDY;
end
% %- rest
% if ~ispc
%     tmp = load('-mat',[spca_dir filesep sprintf('%s_UNIX.study',STUDY_FNAME_REST)]);
%     STUDY_REST = tmp.STUDY;
% else
%     tmp = load('-mat',[spca_dir filesep sprintf('%s.study',STUDY_FNAME_REST)]);
%     STUDY_REST = tmp.STUDY;
% end
%% TEST LOAD
%-
gait_conds = unique({STUDY.datasetinfo(1).trialinfo.cond});
subj_chars = {STUDY.datasetinfo.subject};
eeg_fpaths = {STUDY.datasetinfo.filepath};
eeg_fnames = {STUDY.datasetinfo.filename};
spca_eeg_fpaths = {STUDY_GAIT.datasetinfo.filepath};
spca_eeg_fnames = {STUDY_GAIT.datasetinfo.filename};
%-
% subj_i = 35;
% cond_i = 1;
% spca_fpath = STUDY_GAIT.datasetinfo(subj_i).filepath;
% spca_psd = par_load(spca_fpath,sprintf('cond%s_spca_psd.mat',gait_conds{cond_i}));
% T_DIM = size(spca_psd.psd_corr_based,1);
% F_DIM = size(spca_psd.psd_corr_based,3);
% CL_DIM = length(STUDY.cluster);
%-
def_spca_struct = struct('subj_c',{''}, ...
    'group_n',[], ...
    'cluster_n',[], ...
    'cond_c',{''}, ...
    'comp_c',{''}, ...
    'psd_orig_avg_c',{[]}, ...
    'psd_orig_baselined_c',{[]}, ...
    'psd_rest_c',{[]}, ...
    'psd_corr_based_c',{[]}, ...
    'psd_corr_unbase_c',{[]}, ...
    'psd_corr_psc1',{[]}, ...
    'tf_coeff_c',{[]});
%% TABLE CONSTRUCTION ================================================== %%
struct_store = cell(length(subj_chars),1);
parfor subj_i = 1:length(subj_chars)
% for subj_i = 1:3
    %## INITIATION
    tmp_study = STUDY;
    spca_struct = repmat(def_spca_struct,[1,length(main_cl_inds)*length(gait_conds)]);
    cnt = 1;
    %## LOAD EEG DATA
    try
        %- load spca eeg
        % EEG = pop_loadset('filepath',spca_eeg_fpaths{subj_i},'filename',spca_eeg_fnames{subj_i});
        EEG = load([spca_eeg_fpaths{subj_i} filesep spca_eeg_fnames{subj_i}],'-mat');
        unmix_mat = EEG.etc.spca.unmix_mat;
        bad_eye_ics = EEG.etc.spca.eye_ic_rej;
        %- load ica eeg
        % EEG = pop_loadset('filepath',eeg_fpaths{subj_i},'filename',eeg_fnames{subj_i});
        EEG_ICA = load([spca_eeg_fpaths{subj_i} filesep spca_eeg_fnames{subj_i}],'-mat');
        fprintf('Running subject %s\n',EEG.subject)        
        ic_keep = EEG.etc.urreject.ic_keep;
        ic_rej = EEG.etc.urreject.ic_rej;
        %- perhaps a new way of determining clusters 
        tmp1 = cat(1,EEG.dipfit.model.posxyz);
        tmp2 = cat(1,EEG_ICA.dipfit.model.posxyz);
        dips = zeros(size(tmp2,1),1);
        for i = 1:size(tmp1,1)
            ind = find(all(tmp1(i,:)==tmp2,2));
            if ~isempty(ind)
                dips(i) = ind;
            end
        end
        dips = dips(dips~=0);
        
        %## SANITY CHECKS & IC TRANSFORMS
        %- determine component transformation after eye ic removal
        subj_ics = zeros(1,length(tmp_study.cluster));
        subj_cls = zeros(1,length(tmp_study.cluster));
        subj_unmix = zeros(1,length(tmp_study.cluster));
        kk = 1;
        jj = 1;
        for k = 2:length(tmp_study.cluster)
            set_i = (tmp_study.cluster(k).sets == subj_i);
            comp_i = tmp_study.cluster(k).comps(set_i);
            if ~isempty(comp_i)
                subj_ics(kk:kk+length(comp_i)-1) = comp_i;
                subj_cls(kk:kk+length(comp_i)-1) = repmat(k,1,length(comp_i));
                kk = kk + length(comp_i);
                for j = 1:length(comp_i)
                    subj_unmix(jj) = unmix_mat(2,ic_keep(comp_i(j)) == unmix_mat(1,:));
                    jj = jj + 1;
                end
            end
        end
        subj_ics = subj_ics(subj_ics~=0);
        subj_cls = subj_cls(subj_cls~=0);
        subj_unmix = subj_unmix(subj_unmix~=0);
        %- load spca data
        spca_psd = par_load(spca_eeg_fpaths{subj_i},sprintf('cond%s_spca_psd.mat',gait_conds{1}));
        %- chk ic transforms
        chk = size(spca_psd.psd_corr_based,2) == (length([ic_keep,ic_rej])-length(bad_eye_ics));
        chk2 = (length(ic_keep)==length(tmp_study.datasetinfo(subj_i).comps));
        fprintf('IC numbers check: %i\n',chk2&&chk);
        fprintf('IC''s in clusters: %s\n',strjoin(string(subj_ics),','));
        fprintf('Clusters: %s\t',strjoin(string(subj_cls),',')); 
        fprintf('New ICs: %s\t',strjoin(string(dips(subj_ics)),',')); 
        fprintf('Unmix ICs: %s\t',strjoin(string(subj_unmix),','));

        %## EXTRACT SPCA PSD DATA
        for cond_i = 1:length(gait_conds)
            spca_psd = par_load(spca_eeg_fpaths{subj_i},sprintf('cond%s_spca_psd.mat',gait_conds{cond_i}));
            for i = 1:length(main_cl_inds)
                cl_i = main_cl_inds(i);
                set_i = (tmp_study.cluster(cl_i).sets == subj_i);
                comp_i = tmp_study.cluster(cl_i).comps(set_i);
                if ~isempty(comp_i) && length(comp_i) < 2
                    tmp_c = unmix_mat(2,ic_keep(comp_i) == unmix_mat(1,:));
                    fprintf('%s) assigning IC(spca index,study index) %i(%i,%i) to cluster %i in condition %s...\n',subj_chars{subj_i},ic_keep(comp_i),comp_i,tmp_c,cl_i,gait_conds{cond_i});
                    spca_struct(cnt).psd_corr_based_c = squeeze(spca_psd.psd_corr_based(:,tmp_c,:));
                    spca_struct(cnt).psd_orig_avg_c = squeeze(spca_psd.psd_orig_avg(:,tmp_c,:));
                    spca_struct(cnt).psd_orig_baselined_c = squeeze(spca_psd.psd_orig_baselined(:,tmp_c,:));
                    spca_struct(cnt).psd_rest_c = squeeze(spca_psd.baseline_psd(:,tmp_c,:));
                    spca_struct(cnt).psd_corr_unbase_c = squeeze(spca_psd.psd_corr_based(:,tmp_c,:))+squeeze(spca_psd.baseline_psd(:,tmp_c,:));
                    spca_struct(cnt).psd_corr_psc1 = squeeze(spca_psd.psd_corr_psc1(:,tmp_c,:));
                    % tf_coeff_c{cnt} = spca_psd.coeffs;
                    %- validation
                    % figure;
                    % hold on;
                    % plot(freqs,squeeze(spca_psd.apply_spca_cond.erds_orig(:,tmp_c,:)));
                    % plot(freqs,squeeze(rest_psd(:,tmp_c,:)));
                    % plot(freqs,squeeze(spca_psd.apply_spca_cond.erds_orig(:,tmp_c,:))+squeeze(rest_psd(:,tmp_c,:)));
                    % plot(freqs,squeeze(spca_psd.apply_spca_cond.psc1(:,tmp_c,:)));
                    % hold off;
                    %-
                    spca_struct(cnt).cond_c = gait_conds{cond_i};
                    spca_struct(cnt).subj_c = subj_chars{subj_i};
                    tmp = regexp(subj_chars{subj_i},'\d','match');
                    spca_struct(cnt).group_n = double(string(tmp{1}));
                    spca_struct(cnt).cluster_n = cl_i;
                    spca_struct(cnt).comp_c = sprintf('study_ic, %i; unmix_ic, %i; keep_ic, %i; cluster_ic, %i',comp_i,tmp_c,ic_keep(comp_i),comp_i);
                    cnt = cnt + 1;
                end
            end
        end
        %- remove empty indices
        chk = ~cellfun(@isempty,{spca_struct.subj_c});
        spca_struct = spca_struct(chk);
        struct_store{subj_i} = spca_struct;    
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,subj_chars{subj_i},getReport(e));
    end
end
struct_store = struct_store(~cellfun(@isempty,struct_store));
tt = cat(2,struct_store{:});
tt = struct2table(tt);
par_save(tt,cluster_k_dir,'spca_cluster_table_psd.mat');

