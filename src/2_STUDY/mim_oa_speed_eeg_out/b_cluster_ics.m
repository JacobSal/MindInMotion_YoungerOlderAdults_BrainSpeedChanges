%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/run_b_cluster_ics.sh

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
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
%## statistics & conditions
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
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
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
%- clustering parameters
MIN_ICS_SUBJ = 5; %[2,3,4,5,6,7,8]; % iterative clustering
DO_K_ICPRUNE = true;
CLUSTER_STRUCT = struct('algorithm','kmeans',...
    'clust_k_num',[10,11,12],...
    'clust_k_evals',(5:25),...
    'clust_k_maxiter',10000,...
    'clust_k_replicates',1000,...
    'clust_k_repeat_iters',50,...
    'clust_k_repeat_std',3,...
    'clust_k_robust_maxiter',5);
%- custom params
STD_PRECLUST_COMMAND = {'dipoles','weight',1};
%- datetime override
% dt = '03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% cluster_study_dir = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04232024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
study_dir_name = '04232024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej';
%## soft define
DATA_DIR = [PATHS.src_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
study_fName_1 = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name) filesep 'iclabel_cluster'];
% save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
end
CLUSTER_STRUCT.filename = STUDY.filename;
CLUSTER_STRUCT.filepath = STUDY.filepath;
%% (SET PARAMS)
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
tmp_group_orig = cell(length(ALLEEG),1);
tmp_group_unif = cell(length(ALLEEG),1);
for subj_i = 1:length(ALLEEG)
    tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
    tmp_group_unif{subj_i} = 'Older Adults';
end
%-
%{
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = tmp_group_orig{subj_i};
    STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
end
%}
%- NOTE: partly adapt from bemobil_repeated_clustering
numIC = zeros(length(STUDY.datasetinfo),1);
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
tmp = strsplit(ALLEEG(1).filename,'.');
spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
% pPool = parpool(pp, SLURM_POOL_SIZE);
if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
    fprintf('Calculating Spectograms...\n');
    %- override variables for the stats
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = tmp_group_orig{subj_i};
        STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
        for in_i = 1:length(STUDY.datasetinfo(subj_i).trialinfo)
            STUDY.datasetinfo(subj_i).trialinfo(in_i).group = tmp_group_orig{subj_i};
        end
    end
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
        TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
        fprintf('SUBJECT: %s\n',TMP_STUDY.datasetinfo.subject);
        fprintf('GROUP: %s\n',TMP_STUDY.datasetinfo.group);
        disp(STUDY.datasetinfo(subj_i));
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
%% (STEP 4) CALCULATE CLUSTER SOLUTIONS AFTER PRUNING SUBJECTS WITH LOW NUMBER OF ICs
%## GET EEGLAB PATH
% if ~ispc
%     tmp = strsplit(path,':');
% else
%     tmp = strsplit(path,';');
% end
% % tmp = strsplit(path,';');
% b1 = regexp(tmp,'eeglab','end');
% b2 = tmp(~cellfun(@isempty,b1));
% PATH_EEGLAB = b2{1}; %b2{1}(1:b1{1});
% fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
% %- set default paths for boundary element head model
% PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
% MNI_MRI = [PATH_EEGLAB_BEM filesep 'standard_mri.mat'];
% MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
%##
HIRES_TEMPLATE = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
if ~ispc
    HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
else
    HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
end
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
MNI_VOL = [fpath filesep fname '_dipplotvol.mat'];
MNI_MRI = [fpath filesep fname '_dipplotmri.mat'];
% MNI_MESH = [fpath filesep fname '_mesh.mat'];
mri = load(MNI_MRI);
mri = mri.mri;
vol = MNI_VOL;
%% (KMEANS K SEARCH CRIT) ============================================== %%
%{
eva1 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'silhouette'); % this is the method to find optimal number of clusters (confirmed by EEGlab mailing list)
eva2 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'CalinskiHarabasz');
eva3 = evalclusters(STUDY.etc.preclust.preclustdata, cluster_idx, 'DaviesBouldin');

figure();
subplot(1,3,1)
plot(eva1.InspectedK,eva1.CriterionValues,'-o');hold on;
plot(eva1.OptimalK,eva1.CriterionValues(eva1.InspectedK == eva1.OptimalK),'o','color','r');
xlabel('Clusters');ylabel('Silhouette');
subplot(1,3,2)
plot(eva2.InspectedK,eva2.CriterionValues,'-o');hold on;
plot(eva2.OptimalK,eva2.CriterionValues(eva2.InspectedK == eva2.OptimalK),'o','color','r');
xlabel('Clusters');ylabel('CalinskiHarabasz');
subplot(1,3,3)
plot(eva3.InspectedK,eva3.CriterionValues,'-o');hold on;
plot(eva3.OptimalK,eva3.CriterionValues(eva3.InspectedK == eva3.OptimalK),'o','color','r');
xlabel('Clusters');ylabel('DaviesBouldin');

saveas(gcf, [save_dir filesep 'cluster_kmeans_eval.fig'])
saveas(gcf, [save_dir filesep 'cluster_kmeans_eval.jpg'])
%}
%%
if DO_K_ICPRUNE 
%     for i = 1:length(MIN_ICS_SUBJ)
    parfor (i = 1:length(MIN_ICS_SUBJ),length(MIN_ICS_SUBJ))
        tmp_dir = [save_dir filesep sprintf('icrej_%i',MIN_ICS_SUBJ(i))];
        if ~exist(tmp_dir,'dir')
            mkdir(tmp_dir)
        end
        ics_subjs = cellfun(@(x) size(x,2),{ALLEEG.icawinv});
        idx = ics_subjs > MIN_ICS_SUBJ(i);
        TMP_ALLEEG = ALLEEG(idx);
        %##
        [TMP_STUDY, TMP_ALLEEG] = std_editset([],TMP_ALLEEG,...
                                        'updatedat','off',...
                                        'savedat','off',...
                                        'name',sprintf('temp_study_rejics%i',MIN_ICS_SUBJ(i)),...
                                        'filename',sprintf('temp_study_rejics%i.study',MIN_ICS_SUBJ(i)),...
                                        'filepath',tmp_dir);
        [TMP_STUDY,TMP_ALLEEG] = std_checkset(TMP_STUDY,TMP_ALLEEG);
        [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                            TMP_STUDY.filename,TMP_STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
        %##
        [TMP_STUDY,TMP_ALLEEG] = std_preclust(TMP_STUDY,TMP_ALLEEG,1,STD_PRECLUST_COMMAND);         
        %- store essential info in STUDY struct for later reading
        TMP_STUDY.etc.clustering.preclustparams.clustering_weights = STD_PRECLUST_COMMAND;
        TMP_STUDY.etc.clustering.preclustparams.freqrange = SPEC_PARAMS.freqrange;
        %## (Step 2) CALCULATE REPEATED CLUSTERED SOLUTIONS
        all_solutions = mim_cluster_process(TMP_STUDY,TMP_ALLEEG,...
            'CLUSTER_STRUCT',CLUSTER_STRUCT);
        %## EVALUATE CLUSTER SOLUTIONS
        tmp = TMP_STUDY;
        CLUST_MIN = length(tmp.datasetinfo)/2;
        GROUP_MIN = 15;
        g_chars = unique({tmp.datasetinfo.group});
        cl_chks = zeros(length(CLUSTER_STRUCT.clust_k_num),length(CLUSTER_STRUCT.clust_k_repeate_iters))
        for j = 1:length(CLUSTER_STRUCT.clust_k_num)
            for k = 1:length(CLUSTER_STRUCT.clust_k_repeate_iters)
                %## Calculate dipole positions
                cluster_update = cluster_comp_dipole(TMP_ALLEEG, all_solutions{j}.solutions{k});
                tmp.cluster = cluster_update;
                %## REMOVE BASED ON ICLABEL
                [tmptmp,~,~] = cluster_iclabel_reduce(tmp,TMP_ALLEEG);
                %##
                cl_bools = zeros(length(inds),length(g_chars)+1)
                inds = find(cellfun(@(x) contains(x,'Outlier','IgnoreCase',true) ||...
                    contains(x,'Parentcluster'),{tmp.cluster.name}));
                for l = 1:length(inds)
                    cl_i = inds(l);
                    if length(tmptmp.cluster(cl_i)) < CLUST_MIN
                        cl_bools(l,1) = true;
                    end
                    for g_i = 1:length(g_chars)
                        g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{tmptmp.datasetinfo(tmptmp.cluster(cl_i).sets).group});
                        if sum(g_inds) < GROUP_MIN
                            cl_bools(l,g_i+1) = true;
                        end
                    end
                end
                t_bools = sum(all(cl_bools,2),1)*4;
                tt_bools = sum(any(cl_bools,2),1);
                cl_chks(j,k) = t_bools+tt_bools;
                cl_chks(j,k) = sum(cl_bools,[2,1]);
            end
        end
        TMP_STUDY.etc.cluster_vars = [];
        for j = 1:length(CLUSTER_STRUCT.clust_k_num)
            %-
            clust_i = CLUSTER_STRUCT.clust_k_num(j);
            cluster_dir = [tmp_dir filesep num2str(clust_i)];
            if ~exist(cluster_dir,'dir')
                mkdir(cluster_dir)
            end
            %## Calculate dipole positions
            cluster_update = cluster_comp_dipole(TMP_ALLEEG, all_solutions{j}.solutions{1});
            TMP_STUDY.cluster = cluster_update;
            %## REMOVE BASED ON RV
            % [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
            % [TMP_STUDY,~,~] = cluster_rv_reduce(TMP_STUDY,TMP_ALLEEG);
            %## REMOVE BASED ON ICLABEL
            [TMP_STUDY,~,~] = cluster_iclabel_reduce(TMP_STUDY,TMP_ALLEEG);
            %## REMOVE BASED ON PCA ALG
            %{
            for subj_i = 1:length(TMP_ALLEEG)
                EEG = eeg_checkset(TMP_ALLEEG(subj_i),'loaddata');
                if isempty(EEG.icaact)
                    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
                end
                TMP_ALLEEG(subj_i) = EEG;
            end
            try
                [TMP_STUDY,TMP_ALLEEG,~] = cluster_pca_reduce(TMP_STUDY,TMP_ALLEEG);
            catch e
                fprintf('%s',getReport(e));
                error('cluster_pca_reduce error.');
            end
            for subj_i = 1:length(TMP_ALLEEG)
                TMP_ALLEEG(subj_i).filename = sprintf('%s_pcareduced_comps.set',TMP_ALLEEG(subj_i).subject);
                TMP_ALLEEG(subj_i) = pop_saveset(TMP_ALLEEG(subj_i),...
                    'filename',TMP_ALLEEG(subj_i).filename,...
                    'filepath',TMP_ALLEEG(subj_i).filepath);

            end
            %}
            cluster_update = TMP_STUDY.cluster;
            %- get cluster centroid and residual variance
            cluster_update = cluster_comp_dipole(TMP_ALLEEG, cluster_update);
            TMP_STUDY.cluster = cluster_update;
            %## Look up cluster centroid Brodmann area
            [~,atlas_names,~] = add_anatomical_labels(TMP_STUDY,TMP_ALLEEG);
            for k = 1:length(cluster_update)
                cluster_update(k).analabel = atlas_names{k,2};
            end
            par_save(cluster_update,cluster_dir,sprintf('cl_inf_%i.mat',clust_i));
            TMP_STUDY.etc.cluster_vars(j).fpath = [cluster_dir filesep sprintf('cl_inf_%i.mat',clust_i)];
            TMP_STUDY.etc.cluster_vars(j).inf = {'type','kmeans','k',clust_i,'params',CLUSTER_STRUCT};
        end
        
        %- save
        [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                            TMP_STUDY.filename,TMP_STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
    end
end