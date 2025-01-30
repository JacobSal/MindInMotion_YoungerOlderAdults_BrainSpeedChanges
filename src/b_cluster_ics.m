%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/MIND_IN_MOTION_PRJ/MindInMotion_YoungerOlderAdult_KinEEGCorrs/src/_bash_sh_files/run_b_cluster_ics.sh

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
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
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
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
%- clustering parameters
MIN_ICS_SUBJ = 5; %[2,3,4,5,6,7,8]; % iterative clustering
% CLUSTER_STRUCT = struct('algorithm','kmeans',...
%     'clust_k_num',[7,10,11,12],... %(10-12 seems ideal, see. Makoto's PCA test)
%     'clust_k_evals',[7,10,11,12],...
%     'clust_k_maxiter',10000,...
%     'clust_k_replicates',1000,...
%     'clust_k_repeat_iters',1,...
%     'clust_k_repeat_std',3,... % (inf | INT), INT uses robust_kmeans_CL
%     'clust_k_robust_maxiter',3,... % how many times to redo kmean clustering with removal of unmatching IC's
%     'clust_k_empty_action','drop',...
%     'do_eval_clusters',true);
CLUSTER_STRUCT = struct('algorithm','kmeans',...
    'clust_k_num',11,... %(10-12 seems ideal, see. Makoto's PCA test)
    'clust_k_evals',11,...
    'clust_k_maxiter',10000,...
    'clust_k_replicates',1000,...
    'clust_k_repeat_iters',1,...
    'clust_k_repeat_std',3,... % (inf | INT), INT uses robust_kmeans_CL
    'clust_k_robust_maxiter',10,... % how many times to redo kmean clustering with removal of unmatching IC's
    'clust_k_empty_action','drop',...
    'do_eval_clusters',true);
%{
%- unit test
CLUSTER_STRUCT = struct('algorithm','kmeans',...
    'clust_k_num',[11],... %(10-12 seems ideal, see. Makoto's PCA test)
    'clust_k_evals',(10:12),...
    'clust_k_maxiter',10000,...
    'clust_k_replicates',1000,...
    'clust_k_repeat_iters',10,...
    'clust_k_repeat_std',3,... % (inf | INT), INT uses robust_kmeans_CL
    'clust_k_robust_maxiter',5,...
    'clust_k_empty_action','drop',...
    'do_eval_clusters',true);
%}
%-
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
%- custom params
STD_PRECLUST_COMMAND = {'dipoles','weight',1};
%% (PATHS) ============================================================= %%
%- datetime override
% STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
STUDY_DNAME = '01192025_mim_yaoa_nopowpow_crit_speed';
STUDY_FNAME = 'epoch_study';
%## soft define
studies_dir = [PATHS.data_dir filesep DATA_SET filesep '_studies'];
save_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME) filesep '__iclabel_cluster_kmeansalt_rb3'];
load_dir = [studies_dir filesep sprintf('%s',STUDY_DNAME)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep STUDY_FNAME '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[STUDY_FNAME '.study'],'filepath',load_dir);
    end
end
%% (SET PARAMS) ======================================================== %%
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
%-
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%- NOTE: partly adapt from bemobil_repeated_clustering
numIC = zeros(length(STUDY.datasetinfo),1);
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS ============================= %%
tmp = strsplit(ALLEEG(1).filename,'.');
spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
% pPool = parpool(pp, SLURM_POOL_SIZE);
if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
    fprintf('Calculating Spectograms...\n');
    parfor subj_i = 1:length(ALLEEG)
        tmp_spec_params = SPEC_PARAMS;
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
                    {'specmode',tmp_spec_params.specmode, ...
                    'freqfac',tmp_spec_params.freqfac,...
                    'freqrange',tmp_spec_params.freqrange, ...
                    'logtrials',tmp_spec_params.logtrials});
    end
end
%% (IC PRUNING) ======================================================== %%
%## Remove subjects if they have too few ICs
%*
if ~ispc
    tp = strsplit(path,':');
else
    tp = strsplit(path,';');
end
b1 = contains(tp,'AAL3','IgnoreCase',true);
b2 = tp(b1);
try
    ind = regexp(b2(1),'AAL3','end');
    path_aal3 = b2{1}(1:ind{1});
    fprintf('ALL3 path: %s\n',path_aal3);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('AAL3 path not found.\n');
    end
end
addpath(path_aal3);
for i = 1:length(MIN_ICS_SUBJ)
% parfor (i = 1:length(MIN_ICS_SUBJ),length(MIN_ICS_SUBJ))
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
    %## REDUCE CLUSTER & GET INFORMATION
    TMP_STUDY.etc.cluster_vars = [];
    for j = 1:length(CLUSTER_STRUCT.clust_k_num)
    % parfor j = 1:length(CLUSTER_STRUCT.clust_k_num)
        fprintf('Running cluster struct generation for k=%i\n',CLUSTER_STRUCT.clust_k_num(j));
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
            % [TMP_STUDY,TMP_ALLEEG,~] = cluster_pca_reduce_simple(TMP_STUDY,TMP_ALLEEG);
            [TMP_STUDY,TMP_ALLEEG,~] = cluster_pca_reduce_simple(TMP_STUDY,TMP_ALLEEG);
        catch e
            fprintf('%s',getReport(e));
            error('cluster_pca_reduce error.');
        end
        % for subj_i = 1:length(TMP_ALLEEG)
        %     TMP_ALLEEG(subj_i).filename = sprintf('%s_pcareduced_comps.set',TMP_ALLEEG(subj_i).subject);
        %     TMP_ALLEEG(subj_i) = pop_saveset(TMP_ALLEEG(subj_i),...
        %         'filename',TMP_ALLEEG(subj_i).filename,...
        %         'filepath',TMP_ALLEEG(subj_i).filepath);
        % 
        % end
        %}
        %## ANATOMY CALCULATION
        CALC_STRUCT = struct('cluster_inds',(2:length(TMP_STUDY.cluster)),...
            'save_inf',false,...
            'recalculate',true);
        [TMP_STUDY,dipfit_structs,topo_cells] = eeglab_get_topodip(TMP_STUDY,...
            'CALC_STRUCT',CALC_STRUCT,...
            'ALLEEG',TMP_ALLEEG);
        ANATOMY_STRUCT = struct('atlas_fpath',{{[path_aal3 filesep 'AAL3v1.nii']}},...
            'group_chars',{unique({TMP_STUDY.datasetinfo.group})},...
            'cluster_inds',3:length(TMP_STUDY.cluster),...
            'anatomy_calcs',{{'all aggregate','all centroid'}},... % ('all calcs','group centroid','all centroid','group aggregate','all aggregate')
            'save_dir',cluster_dir,...
            'save_inf',true);
        %(01/16/2025) JS, 'all centroid' option has a bug where indexing
        %doesn't seem to want to work on line 535 in poi_box
        [TMP_STUDY,anat_struct,~,~,txt_out] = eeglab_get_anatomy(TMP_STUDY,...
            'ALLEEG',TMP_ALLEEG,...
            'ANATOMY_STRUCT',ANATOMY_STRUCT);
        %-
        cl_tmp = TMP_STUDY.cluster;
        atlas_char = 'AAL3v1.nii';
        inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'aggregate label for all');
        at_out = {anat_struct(inds).anatomy_label};
        at_ind = [anat_struct(inds).cluster];
        [cl_tmp(at_ind).agg_anat_all] = at_out{:};
        %-
        % inds = strcmp({anat_struct.atlas_label},atlas_char) & strcmp({anat_struct.calculation},'centroid label for all');
        % at_out = {anat_struct(inds).anatomy_label};
        % at_ind = [anat_struct(inds).cluster];
        % [cl_tmp(at_ind).ctr_anat_all] = at_out{:};
        %## OLD ANATOMY CALC (USES MNI ATLAS)
        [~,atlas_names,~] = add_anatomical_labels(TMP_STUDY,TMP_ALLEEG);
        for k = 1:length(cl_tmp)
            cl_tmp(k).analabel = atlas_names{k,2};
        end
        %## UPDATE
        TMP_STUDY.cluster = cl_tmp;
        
        par_save(cl_tmp,cluster_dir,sprintf('cl_inf_%i.mat',clust_i));
        TMP_STUDY.etc.cluster_vars(j).fpath = [cluster_dir filesep sprintf('cl_inf_%i.mat',clust_i)];
        TMP_STUDY.etc.cluster_vars(j).inf = {'type','kmeans','k',clust_i,'params',CLUSTER_STRUCT};
        fprintf('done.\n');
    end
    
    %- save
    [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                        TMP_STUDY.filename,TMP_STUDY.filepath,...
                                        'RESAVE_DATASETS','off');
end