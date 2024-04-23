%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_terrain_clin/run_b_cluster_ics.sh

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
global ADD_CLEANING_SUBMODS %#ok<GVMIS>
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
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = false;
%## statistics & conditions
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result? 
SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
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
MIN_ICS_SUBJ = [5]; %[2,3,4,5,6,7,8]; % iterative clustering
% K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
CLUSTER_SWEEP_VALS = [12]; %[10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
% DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = true;
% DO_K_SWEEPING = false;
% (08/21/2023) JS, this currenlty doesn't do anything but take up more
% memory.
REPEATED_CLUSTERING_STD = 3;
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',[],...
    'filepath',[],...
    'outliers',inf(),...
    'maxiter',200);
%- custom params
STD_PRECLUST_COMMAND = {'dipoles','weight',1};
%- datetime override
% dt = '03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% cluster_study_dir = '03232023_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
% study_dir_name = '04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
study_dir_name = '04162024_MIM_OAN57_antsnormalize_iccREMG0p4_powpow0p3_skull0p01';
%## soft define
STUDIES_DIR = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
study_fname = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
save_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',study_dir_name)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fname '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fname '.study'],'filepath',load_dir);
    end
end
CLUSTER_PARAMS.filename = STUDY.filename;
CLUSTER_PARAMS.filepath = STUDY.filepath;
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

%-
STUDY_COND_DESI = {{'subjselect',{},...
            'variable1','cond','values1',{'flat','low','med','high'},...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',{'0p25','0p5','0p75','1p0'},...
            'variable2','group','values2',{}}};
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
        %- remove components outside the brain
        % [TMP_STUDY, TMP_ALLEEG] = std_editset(TMP_STUDY,TMP_ALLEEG,...
        %                                 'commands',{'inbrain','on'});                  
        %- store essential info in STUDY struct for later reading
        % freqrange = [4 60];
        TMP_STUDY.etc.clustering.preclustparams.clustering_weights = STD_PRECLUST_COMMAND;
        TMP_STUDY.etc.clustering.preclustparams.freqrange = SPEC_PARAMS.freqrange;
        %## (Step 2) CALCULATE REPEATED CLUSTERED SOLUTIONS
        all_solutions = mim_cluster_process(TMP_STUDY,TMP_ALLEEG,tmp_dir,...
            'CLUSTER_PARAMS',CLUSTER_PARAMS,...
            'REPEATED_CLUSTERING_STD',REPEATED_CLUSTERING_STD,...
            'MAX_REPEATED_ITERATIONS',MAX_REPEATED_ITERATIONS,...
            'KS_TO_TEST',CLUSTER_SWEEP_VALS);
        
        %## REMOVE DIPOLES OUTSIDE THE BRAIN AREA
        TMP_STUDY.etc.cluster_vars = [];
        for j = 1:length(CLUSTER_SWEEP_VALS)
            %-
            clust_i = CLUSTER_SWEEP_VALS(j);
            cluster_dir = [tmp_dir filesep num2str(clust_i)];
            if ~exist(cluster_dir,'dir')
                mkdir(cluster_dir)
            end
            %## Calculate dipole positions
            cluster_update = cluster_comp_dipole(TMP_ALLEEG, all_solutions{j}.solutions{1});
            TMP_STUDY.cluster = cluster_update;
            %## (PLOT) Looking for dipoles outside of brain.
            vol = load(MNI_VOL);
            try
                vol = vol.vol;
            catch
                vol = vol.mesh;
            end
            rmv_dip = [];
            fig = figure('color','w');
            ft_plot_mesh(vol.bnd(3));
            hold on;
            for c = 2:length(TMP_STUDY.cluster)
                for d = 1:size(TMP_STUDY.cluster(c).all_diplocs,1)
                    depth = ft_sourcedepth(TMP_STUDY.cluster(c).all_diplocs(d,:), vol);
                    if depth > 0
                        rmv_dip = [rmv_dip;[c,d]];
                        plot3(TMP_STUDY.cluster(c).all_diplocs(d,1),TMP_STUDY.cluster(c).all_diplocs(d,2),TMP_STUDY.cluster(c).all_diplocs(d,3),'*-');
                    end
                end
            end
            hold off;
            drawnow;
            saveas(fig,[cluster_dir filesep sprintf('ics_out_of_brain.fig')]);
            if ~isempty(rmv_dip)
                sets_ob = zeros(size(rmv_dip,1),1);
                comps_ob = zeros(size(rmv_dip,1),1);
                clusts = unique(rmv_dip(:,1));
                cnt = 1;
                for c_i = 1:length(clusts)
                    inds = clusts(c_i)==rmv_dip(:,1);
                    d_i = rmv_dip(inds,2);
                    sets_ob(cnt:cnt+length(d_i)-1) = TMP_STUDY.cluster(clusts(c_i)).sets(d_i);
                    comps_ob(cnt:cnt+length(d_i)-1) = TMP_STUDY.cluster(clusts(c_i)).comps(d_i);
                    TMP_STUDY.cluster(clusts(c_i)).sets(d_i) = [];
                    TMP_STUDY.cluster(clusts(c_i)).comps(d_i) = [];
                    cnt = cnt + length(d_i);
                end
                TMP_STUDY.cluster(end+1).sets = sets_ob';
                TMP_STUDY.cluster(end).comps = comps_ob';
                TMP_STUDY.cluster(end).name = 'Outlier cluster_outside-brain';
                TMP_STUDY.cluster(end).parent = TMP_STUDY.cluster(2).parent;
                TMP_STUDY.cluster(end).algorithm = 'ft_sourcedepth < 0';
            end
            %## REMOVE BASED ON RV
            % [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
            [TMP_STUDY,~,~] = cluster_rv_reduce(TMP_STUDY,TMP_ALLEEG);
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
            TMP_STUDY.etc.cluster_vars(j).inf = {'type','kmeans','k',clust_i,'params',CLUSTER_PARAMS};
        end
        
        %- save
        [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
                                            TMP_STUDY.filename,TMP_STUDY.filepath,...
                                            'RESAVE_DATASETS','off');
    end
end