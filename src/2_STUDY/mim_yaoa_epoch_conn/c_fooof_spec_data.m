%   Project Title: MIM YOUNGER AND OLDER ADULTS CONNECTIVITY SWING VS.
%   STANCE
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
%##
FORCE_RECALC_SPEC = false;
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
% cluster_study_dir = '12012023_OAYA104_icc0p65-0p2_changparams';
% cluster_study_dir = '01232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3';
cluster_study_dir = '03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% cluster_study_dir = '01232023_MIM_YAN32_antsnormalize_iccREMG0p4_powpow0p3_conn';
ICLABEL_EYE_CUTOFF = 0.75;
%- study group and saving
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
load_dir_1 = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir)];
%- load cluster
CLUSTER_DIR = [STUDIES_DIR filesep sprintf('%s',cluster_study_dir) filesep 'cluster'];
CLUSTER_STUDY_FNAME = 'temp_study_rejics5';
CLUSTER_STUDY_DIR = [CLUSTER_DIR filesep 'icrej_5'];
CLUSTER_K = 12;
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%##
CLUSTER_CLIM_MATCH = [];
SUB_GROUP_FNAME = ['spec_of_oh'];
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H2000''s','H3000''s'}}};
%% ===================================================================== %%
%## LOAD STUDY
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '_UNIX.study'],'filepath',CLUSTER_STUDY_DIR);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAME '.study'],'filepath',CLUSTER_STUDY_DIR);
end
cl_struct = par_load([CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)],sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
condition_gait = unique({STUDY.datasetinfo(1).trialinfo.cond}); %{'0p25','0p5','0p75','1p0','flat','low','med','high'};
subject_chars = {STUDY.datasetinfo.subject};
%-
fPaths = {STUDY.datasetinfo.filepath};
fNames = {STUDY.datasetinfo.filename};
CLUSTER_PICKS = main_cl_inds(2:end);
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
%## ersp plot per cluster per condition
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');

%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
tmp = strsplit(ALLEEG(1).filename,'.');
spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
    parfor (subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/2))
        EEG = ALLEEG(subj_i);
        TMP_STUDY = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %- overrride datasetinfo to trick std_precomp to run.
        TMP_STUDY.datasetinfo = TMP_STUDY.datasetinfo(subj_i);
        TMP_STUDY.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(TMP_STUDY, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_PARAMS.specmode,'freqfac',SPEC_PARAMS.freqfac,...
                    'freqrange',SPEC_PARAMS.freqrange,'logtrials',SPEC_PARAMS.logtrials});
    end
end
 %% MAKE DESIGNS
%## NOTE (07/18/2023) JS, scripts/functions adapted from CL, AS, NJ scripts:
% PlotAndSaveERSP_CL_V3.m, Stats_Plotting_ERSPs_local_allCondBaseline_V3.m,
% Figures_Plotting_ERSPs_local_allCondBaseline_V3.m, plotERSP_parfor.m
%-
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%## Get Cluster Information
cluster_dir = [CLUSTER_STUDY_DIR filesep sprintf('%i',CLUSTER_K)];
if ~isempty(SUB_GROUP_FNAME)
    spec_data_dir = [cluster_dir filesep 'spec_data' filesep SUB_GROUP_FNAME];
else
    spec_data_dir = [cluster_dir filesep 'spec_data'];
end
if ~exist(spec_data_dir,'dir')
    mkdir(spec_data_dir)
end
%%
%- stores
cnt = 1;
cnt2 = 1;
spec_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
spec_ss_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
ersp_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
ersp_norm_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
ersp_normcb_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
ersp_singtrial_fpaths = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
load_ind_cl = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
clust_ind_cl = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
des_ind = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
cnt_ind = cell(length(STUDY.design)*length(CLUSTER_PICKS),1);
%- loop designs
for des_i = 1:length(STUDY.design)
    store_1 = cell(length(CLUSTER_PICKS),1);
    store_2 = cell(length(CLUSTER_PICKS),1);
    store_3 = cell(length(CLUSTER_PICKS),1);
    store_4 = cell(length(CLUSTER_PICKS),1);
    store_5 = cell(length(CLUSTER_PICKS),1);
    store_6 = cell(length(CLUSTER_PICKS),1);
    store_7 = cell(length(CLUSTER_PICKS),1);
    store_8 = cell(length(CLUSTER_PICKS),1);
    store_9 = cell(length(CLUSTER_PICKS),1);
    store_10 = cell(length(CLUSTER_PICKS),1);
    STUDY.currentdesign = des_i;
    design_char = [];
    for i = 1:length(STUDY.design(des_i).variable)
        if i == 1
            design_char = [STUDY.design(des_i).variable(1).value{:}];
        else
            design_char = [design_char '_' STUDY.design(des_i).variable(i).value{:}];
        end
    end
    %## (PARFOR) Compute Specs & ERSPs
    parfor (i = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
%     for i = 1:length(CLUSTER_PICKS)
        cluster_i = CLUSTER_PICKS(i);
        %- generate power spectrums for each cluster and condition
        [~,spec_savef,spec_subjcorr_savef] = mim_gen_spec_data(STUDY,ALLEEG,...
                des_i,cluster_i,design_char,spec_data_dir);
        store_1{i} = spec_savef;
        store_2{i} = spec_subjcorr_savef;
        store_7{i} = i;
        store_8{i} = cluster_i;
        store_9{i} = des_i;
    end
    spec_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_1;
    spec_ss_fpaths(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_2;
    load_ind_cl(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_7;
    clust_ind_cl(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_8;
    des_ind(cnt:cnt+length(CLUSTER_PICKS)-1,1) = store_9;
    for k = cnt:cnt+length(CLUSTER_PICKS)-1
        cnt_ind{k,1} = k;
    end
    cnt = cnt + length(CLUSTER_PICKS);
end
%## Store Paths in Struct
gen_data_struct = struct('spec_fpaths',spec_fpaths,...
    'spec_ss_fpaths',spec_ss_fpaths,...
    'load_ind_cl',load_ind_cl,...
    'clust_ind_cl',clust_ind_cl,...
    'des_ind',des_ind,...
    'cnt_ind',cnt_ind);
STUDY.etc.mim_gen_ersp_data = gen_data_struct;
par_save(gen_data_struct,spec_data_dir,'spec_data_struct.mat');
[~,~] = parfunc_save_study(STUDY,ALLEEG,...
                                STUDY.filename,spec_data_dir,...
                                'RESAVE_DATASETS','off');
%%