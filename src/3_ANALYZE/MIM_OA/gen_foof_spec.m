%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%- run .sh
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_ersp_plots.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '3_ANALYZE' filesep 'MIM_OA'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
%     SLURM_POOL_SIZE = 2;
%     pp = parcluster('local');
%     pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
end
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
DATA_SET = 'MIM_dataset';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = true;
FORCE_RECALC_ERSP = true;
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false; 
DO_SUBJ_PLOTS = false;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%##
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
do_multivariate_data = 1;
evaluate_method = 'min_rv';
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
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
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% (08/03/2023) JS, changing fieldtripnaccu to 10,000; changing
% fieldtripmcorrect to cluster; method to perm; (these are the parameters
% set inside CL's PlotAndSaveERSP_CL_V3.m...
% pipeline although this doesn't align with her YA manuscript methods?
% (08/06/2023) JS, changing fieldtripnaccu to 2000 again and mcorrect to fdr...
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200]);
% (08/03/2023) JS, turning subbaseline to off to align with methods set
% inside CL's PlotAndSaveERSP_CL_V3.m...
%- datetime override
dt = '07222023_MIM_OAN79_subset_prep_verified_gait_conn';
%## Soft Define
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## ADMIN SET
ATLAS_PATH = [PATH_ROOT filesep 'par_EEGProcessing' filesep 'submodules',...
    filesep 'fieldtrip' filesep 'template' filesep 'atlas'];
% see. https://www.fieldtriptoolbox.org/template/atlas/
ATLAS_FPATHS = {[ATLAS_PATH filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [ATLAS_PATH filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [ATLAS_PATH filesep 'brainnetome' filesep 'BNA_MPM_thr25_1.25mm.nii'],...
    [ATLAS_PATH filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat'],...
    [ATLAS_PATH filesep 'vtpm' filesep 'vtpm.mat'],...
    [ATLAS_PATH filesep 'yeo' filesep 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'],...
    [ATLAS_PATH filesep 'brainweb' filesep 'brainweb_discrete.mat']}; % also a discrete version of this
%- convert SUB_DIR
SUB_DIR = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster';
if ~ispc
    SUB_DIR = convertPath2UNIX(SUB_DIR);
else
    SUB_DIR = convertPath2Drive(SUB_DIR);
end
%## USER SET
LOAD_DIFFERENT_STUDY = {true,true};
CLUSTER_K_PICKS = [14,14];
CLUSTER_STUDY_FNAMES = {'temp_study_rejics6','temp_study_rejics5'};
CLUSTER_DIRS = {[SUB_DIR filesep 'icrej_6' filesep '14'],...
    [SUB_DIR filesep 'icrej_5' filesep '14']};
CLUSTER_FILES = {'cl_inf_14.mat','cl_inf_14.mat'};
CLUSTER_STUDY_DIRS = {[SUB_DIR filesep 'icrej_6'],...
    [SUB_DIR filesep 'icrej_5']};
POSS_CLUSTER_CHARS = {};
%% (STEP 2) PLOT
%##
for k_i = 1:length(CLUSTER_K_PICKS)
% parfor (k_i = 1:length(CLUSTER_K_PICKS),length(CLUSTER_K_PICKS))
    fprintf('Loading Cluster K=%i',CLUSTER_K_PICKS(k_i));
    %## Loop Params
    clust_i = CLUSTER_K_PICKS(k_i);
    if ~ispc
        cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
    else
        cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
    end
    plot_store_dir = [cluster_dir filesep 'plots_out'];
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    spec_data_dir = [cluster_dir filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        error('error. path %s doesn''t exist',spec_data_dir);
    end
    if ~exist([spec_data_dir filesep CLUSTER_STUDY_FNAMES{k_i} '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '_UNIX.study'],'filepath',spec_data_dir);
        else
            [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_FNAMES{k_i} '.study'],'filepath',spec_data_dir);
        end
    end
    
    %## CALCULATE GRANDAVERAGE WARPTOs
    for subj_i = 1:length(ALLEEG)
        %- assign percondition timewarping
        ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    %     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
    end
    allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    % allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    for subj_i = 1:length(ALLEEG)
        allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
    end
    % grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
    averaged_warpto_events = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
    %## (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
    TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
    ERSP_TIMERANGE=[averaged_warpto_events(1), averaged_warpto_events(end)];
    STUDY.etc.averaged_warpto_events = averaged_warpto_events;
    fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
    disp(averaged_warpto_events);
    %## RE-POP PARAMS
    STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
    STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
          'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_TIMERANGE);
    %## Cluster Update
    cluster_update = par_load(cluster_dir,CLUSTER_FILES{k_i});
    STUDY.cluster = cluster_update;
    %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps(STUDY);
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); %valid_clusters; %main_cl_inds(2:end); %valid_clusters
    %## PLOT cluster based information
%     mim_gen_cluster_figs(STUDY,TMP_ALLEEG,CLUSTER_DIRS{k_i},...
%         'CLUSTERS_TO_PLOT',main_cl_inds);
    %% Loop Through Designs
    for des_i = 1:length(STUDY.design)
        cond_test = STUDY.design(des_i).variable(1).value;
        fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
        STUDY.currentdesign = des_i;
        fprintf('Current design: %i\n',STUDY.currentdesign);
        fprintf('Statistics Parameters:\n');
        disp(STUDY.etc.statistics)
        fprintf('Statistics Fieldtrip Parameters:\n');
        disp(STUDY.etc.statistics.fieldtrip)
        fprintf('Statistics EEGLAB Parameters:\n');
        disp(STUDY.etc.statistics.eeglab)
        fprintf('ERSP Parameters:\n');
        disp(STUDY.etc.erspparams)
        cl_inds = [STUDY.etc.mim_gen_ersp_data.clust_ind_cl];
        des_inds = [STUDY.etc.mim_gen_ersp_data.des_ind];
%         des_cls = TMP_STUDY.etc.mim_gen_ersp_data.clust_ind_cl([TMP_STUDY.etc.mim_gen_ersp_data.des_ind] == des_i);
        parfor (j = 1:length(CLUSTER_PICKS),length(CLUSTER_PICKS))
%         for j = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(j);
            cluster_load_ind = find(logical(cl_inds == cluster_i) & logical(des_inds == des_i));
%             cluster_load_ind = cluster_i;
%             cluster_load_ind = TMP_STUDY.etc.mim_gen_ersp_data(des_i,cluster_i).cluster_n(1,cluster_i);
            %- defaults
            allersp = {};
            alltimes = [];
            allfreqs = [];
            pcond = {};
            pgroup = {};
            pinter = {};
            %- LOAD OPTION 1
            
        end
    end
end

% STUDY_walking = importdata(fullfile(MainDirectory,studyname ));
%% FOOOF SETUP & PYTHON
%{
% Check which python is being used
pyversion

% The print out from above should tell you which Python you are calling
%  It should show that you are using Python version 3.X
%  If you are using anaconda, it should show your Python is in the anaconda folder
%  If either of these things are not right, reset which Python you are using, as below

% Set python version to use
%  Note: you must do this first thing after opening Matlab (relaunch if you need to)
%  You should only ever have to run this at most, once.
%  You might need to change the path to where your python or anaconda install is
%    For example, your anaconda folder might be `anaconda3` instead of `anaconda`
%    or your anaconda path might be somewhere else, for example, '/opt/anaconda3/bin/python'
%## MANUALLY SET TO A VERSION LOWER THAT 3.10
pyversion('C:\Users\jsalminen\AppData\Local\Programs\Python\Python37\python.EXE');
% pyversion('C:\Users\jsalminen\AppData\Local\Microsoft\WindowsApps\python37.EXE');
% pyversion('/cygdrive/c/Users/jsalminen/AppData/Local/Microsoft/WindowsApps/python3.8')
pe = pyenv;
%}
%% Setup Fooof
SAVE_DATA = true;
settings = struct();  % Use defaults
settings.peak_width_limits = [1 8];%default [1 6] % Amanda used [1 8] in her paper
settings.min_peak_height = 0.05;
% settings.peak_threshold = 2;
settings.max_n_peaks = 3; % originally set to be 2 - 2023-06-07
% the settings are consitent with fooof on github
f_range = [3, 40];
theta_band = [4, 8];
alpha_band = [8 12];
beta_band  = [12 30];
cond_terrains = {'flat','low','med','high'};
cond_speeds = {'0p25','0p5','0p75','1p0'};
COND_CHARS = {cond_terrains,cond_speeds};
% spec_data_cl2_0p250p50p751p0_subbaselined_commonbase
% spec_data_cl2_flatlowmedhigh_subtractmean
keywords = {'Terrain','Speed'};
save_path = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\icrej_6\14';
data_fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\07222023_MIM_OAN79_subset_prep_verified_gait_conn\cluster\icrej_6\14\spec_data';
mkdir(fullfile(save_path, 'Fooof_Plots'))
%%
% if ~isfield(STUDY.cluster,'topo'), STUDY.cluster(1).topo = []; end
% for clus = 3:length(STUDY.cluster) % For each cluster requested
%     if isempty(STUDY.cluster(clus).topo)
%         STUDY = std_readtopoclust_CL(STUDY,ALLEEG, clus);% Using this custom modified code to allow taking average within participant for each cluster
%     end
% end
for g = 1:2
    for k = main_cl_inds(2:end)      
        file_mat = [data_fpath filesep sprintf('spec_data_cl%i_%s_subtractmean.mat',k,[COND_CHARS{g}{:}])]; %['readSPEC_',num2str(g),'_',keywords{g},'_subSubjectMean.mat'];
%         file_mat = [data_fpath filesep sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',k,[COND_CHARS{g}{:}])]; %['readSPEC_',num2str(g),'_',keywords{g},'_subSubjectMean.mat'];
        tmp = par_load(file_mat,[]);
        % Note: spec_subj_mean_stats separate young and older adults
        specdata = {tmp.specdata}';
        specfreqs = {tmp.specfreqs}';
        specfreqs = specfreqs{1};
        %% Run fooof
        % Input should be in linear spacing   
        i_ind = 0;
        for group = 1:size(specdata,2) % in case there is young and old adult group
            for cond = 1:size(specdata,1) % different level of terrains
                specdata_nolog = 10.^(specdata{cond,group}/10);
                % Run FOOOF
                return_model = true;
                for i = 1:size(specdata{cond,group},2)
                    fooof_results{g}{cond,group}{i} = fooof(specfreqs, specdata_nolog(:,i), f_range, settings, return_model);
                    fooof_group_results_org{g}{k}(i_ind + i).subID = cluster_update(k).sets(i);
                    fooof_group_results_org{g}{k}(i_ind + i).compID = cluster_update(k).comps(i);
                    fooof_group_results_org{g}{k}(i_ind + i).study = g;%1 = terrain, 2 = speed
                    fooof_group_results_org{g}{k}(i_ind + i).cond = cond;
                    fooof_group_results_org{g}{k}(i_ind + i).group = group;
                    fooof_group_results_org{g}{k}(i_ind + i).cluster = k;
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_exp = fooof_results{g}{cond,group}{i}.aperiodic_params(2);
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_offset = fooof_results{g}{cond,group}{i}.aperiodic_params(1);
                    fooof_group_results_org{g}{k}(i_ind + i).central_freq = fooof_results{g}{cond,group}{i}.peak_params(:,1);
                    fooof_group_results_org{g}{k}(i_ind + i).power = fooof_results{g}{cond,group}{i}.peak_params(:,2);
                    fooof_group_results_org{g}{k}(i_ind + i).r_squared = fooof_results{g}{cond,group}{i}.r_squared;

                    %------------ Compute average power after flatten curve
%                     fooof_diff = fooof_results{g}{cond,group}{i}.power_spectrum - fooof_results{g}{cond,group}{i}.ap_fit;
                    % Super important, the output is already logged, the
                    % only difference is the magnitude by 10
                    fooof_diff = 10*(fooof_results{g}{cond,group}{i}.power_spectrum) - 10*(fooof_results{g}{cond,group}{i}.ap_fit);
                    fooof_freq = fooof_results{g}{cond,group}{i}.freqs;
                    fooof_group_results_org{g}{k}(i_ind + i).theta_avg_power = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_avg_power = mean(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).beta_avg_power = mean(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));

                    % data structure needs to be freq x subject
                    fooof_diff_store{g}{k}{cond}(:,i) = fooof_diff';
                    fooof_apfit_store{g}{k}{cond}(:,i) = 10*(fooof_results{g}{cond,group}{i}.ap_fit);
                    
                    % - store original spec data
                    spec_data_original{g}{k}{cond} = specdata{cond,group}(specfreqs >= f_range(1) & specfreqs <= f_range(2),:);
                end
                i_ind = i_ind + size(specdata{cond,group},2);

            end
        end

    %         writematrix(specdata_nolog,'D:\GoogleDriveUFL\Postdoc_Code\FOOOF\specdata_nolog.csv')
        %% Sanity check - 
        plot_all_comp = 1;
        if plot_all_comp
            log_freq = 0;
            single_fig = 1;
            cond = 1; group = 1;
            figure();set(gcf,'color','white','position',[100 300 800 800]);
            for i = 1:length(fooof_results{g}{cond,group})

                if i <= 15
                    subplot(5,3,i)
                    fooof_plot_CL(fooof_results{g}{cond,group}{i},log_freq,single_fig,0,'outline_dot');
                    set(gcf,'color','white');
                    xlabel('Frequency(Hz)');
                    ylabel('log10(Power)');
                    title(['SubID:',num2str(cluster_update(k).sets(i)),' CompID:',num2str(cluster_update(k).comps(i))]);
                    saveas(gcf, fullfile(save_path, 'Fooof_Plots',['cluster_',keywords{g},'_',num2str(k),'_1.fig']))
                elseif i <= 30
                    if i == 16; figure();set(gcf,'color','white');end;
                    subplot(5,3,i-15)
                    fooof_plot_CL(fooof_results{g}{cond,group}{i},log_freq,single_fig,0,'outline_dot');
                    set(gcf,'color','white');
                    xlabel('Frequency(Hz)');
                    ylabel('log10(Power)');
                    title(['SubID:',num2str(cluster_update(k).sets(i)),' CompID:',num2str(cluster_update(k).comps(i))]);
                    saveas(gcf, fullfile(save_path, 'Fooof_Plots',['cluster_',keywords{g},'_',num2str(k),'_2.fig']))
                end
            end
        end

        plot_flattened_curve = 1
        if plot_flattened_curve
            log_freq = 0;
            single_fig = 1;
            cond = 1; group = 1;
            figure();set(gcf,'color','white','position',[100 300 800 800]);
            for i = 1:length(fooof_results{g}{cond,group})
                if i <= 15
                    subplot(5,3,i)
                    fooof_plot_CL(fooof_results{g}{cond,group}{i},log_freq,single_fig,0,'flatten_curve');
                    set(gcf,'color','white');
                    xlabel('Frequency(Hz)');
                    ylabel('log10(Power)');
                    title(['SubID:',num2str(cluster_update(k).sets(i)),' CompID:',num2str(cluster_update(k).comps(i))]);
                    saveas(gcf, fullfile(save_path, 'Fooof_Plots',['FLATTEN cluster_',keywords{g},'_',num2str(k),'_1.fig']))
%                     saveas(gcf, fullfile('M:\liu.chang1\STUDY-MiM-Imagined-2022_03_14_HY_NoBaseCorr\ERSP_Plots\FOOOF_plots',['FLATTEN cluster_',keyword,'_',num2str(k),'_1.fig']))
                elseif i <= 30
                    if i == 16; figure();set(gcf,'color','white','position',[100 300 800 800]);end;
                    subplot(5,3,i-15)
                    fooof_plot_CL(fooof_results{g}{cond,group}{i},log_freq,single_fig,0,'flatten_curve');
                    set(gcf,'color','white');
                    xlabel('Frequency(Hz)');
                    ylabel('log10(Power)');
                    title(['SubID:',num2str(cluster_update(k).sets(i)),' CompID:',num2str(cluster_update(k).comps(i))]);
                    saveas(gcf, fullfile(save_path, 'Fooof_Plots',['FLATTEN cluster_',keywords{g},'_',num2str(k),'_2.fig']))
%                     saveas(gcf, fullfile('M:\liu.chang1\STUDY-MiM-Imagined-2022_03_14_HY_NoBaseCorr\ERSP_Plots\FOOOF_plots',['FLATTEN cluster_',keyword,'_',num2str(k),'_2.fig']))
                end
            end
        end

    end
    
end
if SAVE_DATA
    save(fullfile(save_path,...
        ['fooof_results_summary,.mat']),'fooof_group_results_org','fooof_diff_store');
end
keyboard
%% Sanity check - time series plots from aperiodic subtraction
k = 3;
figure('color','white');
data_min = [min(fooof_diff_store{g}{k}{1});min(fooof_diff_store{g}{k}{2});min(fooof_diff_store{g}{k}{3});...
min(fooof_diff_store{g}{k}{4});];
data_max = [max(fooof_diff_store{g}{k}{1});min(fooof_diff_store{g}{k}{2});max(fooof_diff_store{g}{k}{3});...
max(fooof_diff_store{g}{k}{4});];

for i = 1:4
    subplot(4,1,i)
    data = fooof_diff_store{g}{k}{i};
    plot(fooof_freq,data,'color',color.terrain(i,:));
    ylabel('log10(Power)')
    ylim([min(data_min,[],'all') max(data_max,[],'all')]);
end
xlabel('Frequency(Hz)');
title(['Cluster ',num2str(k)]);

%% Determine peak occurs at different band, don't think I am using this
% Looks like 
% Define frequency bands of interest
for g = 1:2
for k = 3:length(fooof_group_results_org{g})
    for i = 1:length(fooof_group_results_org{g}{k})
        if ~isempty(fooof_group_results_org{g}{k}(i).central_freq)
            fooof_group_results_org{g}{k}(i).theta = [];
            fooof_group_results_org{g}{k}(i).alpha = [];
            fooof_group_results_org{g}{k}(i).beta = [];
            for j = 1:length(fooof_group_results_org{g}{k}(i).central_freq)
                cf = fooof_group_results_org{g}{k}(i).central_freq(j);
                if cf > theta_band(1) & cf <= theta_band(2)
                    fooof_group_results_org{g}{k}(i).theta = [fooof_group_results_org{g}{k}(i).theta; cf fooof_group_results_org{g}{k}(i).power(j)];
                elseif cf > alpha_band(1) & cf <= alpha_band(2)
                    fooof_group_results_org{g}{k}(i).alpha = [fooof_group_results_org{g}{k}(i).alpha; cf fooof_group_results_org{g}{k}(i).power(j)];
                elseif cf > beta_band(1) & cf <= beta_band(2)
                    fooof_group_results_org{g}{k}(i).beta = [fooof_group_results_org{g}{k}(i).beta; cf fooof_group_results_org{g}{k}(i).power(j)];
                end
            end
            if length(fooof_group_results_org{g}{k}(i).theta) > 1
                [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).theta(:,1)-6));
                temp_power = fooof_group_results_org{g}{k}(i).theta(:,2);
                fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).theta(indx,1) temp_power(indx)];
            end
            if length(fooof_group_results_org{g}{k}(i).alpha) > 1
                [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).alpha(:,1)-10));
                temp_power = fooof_group_results_org{g}{k}(i).alpha(:,2);
                fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).alpha(indx,1) temp_power(indx)];
            end
            if length(fooof_group_results_org{g}{k}(i).beta) > 1
                [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).beta(:,1)-20));
                temp_power = fooof_group_results_org{g}{k}(i).beta(:,2);
                fooof_group_results_org{g}{k}(i).beta_center = [fooof_group_results_org{g}{k}(i).beta(indx,1) temp_power(indx)];
            end
        end
    end
end
end
% save('M:\liu.chang1\STUDY-MiM-Imagined-2022_03_14_HY_NoBaseCorr\ERSP_Plots\FOOOF_plots\fooof_results_all.mat','fooof_group_results_org');

%% Create table from group results, take mean across participants ICs    
C1_table = table;
for g = 1:2
    for k = 3:length(fooof_group_results_org{g})
        temp_table_C1 = fooof_group_results_org{g}{k};
        for i = 1:length(temp_table_C1)
            if isempty(temp_table_C1(i).alpha)
                temp_table_C1(i).alpha = [NaN, NaN];
            end
            if isempty(temp_table_C1(i).beta)
                temp_table_C1(i).beta = [NaN, NaN];
            end
            if isempty(temp_table_C1(i).alpha_center)
                temp_table_C1(i).alpha_center = [NaN, NaN];
            end

            if isempty(temp_table_C1(i).beta_center)
                temp_table_C1(i).beta_center = [NaN, NaN];
            end

            [~,idx_a] = max(temp_table_C1(i).alpha(:,2));
            [~,idx_b] = max(temp_table_C1(i).beta(:,2));

            C1_table = vertcat(C1_table,table(temp_table_C1(i).subID,temp_table_C1(i).compID,temp_table_C1(i).study,...
                temp_table_C1(i).cond,temp_table_C1(i).group,temp_table_C1(i).cluster,temp_table_C1(i).aperiodic_exp,temp_table_C1(i).aperiodic_offset,...
                temp_table_C1(i).r_squared,temp_table_C1(i).alpha(idx_a,1),temp_table_C1(i).alpha(idx_a,2),temp_table_C1(i).beta(idx_b,1),temp_table_C1(i).beta(idx_b,2),...
                temp_table_C1(i).alpha_center(1,1),temp_table_C1(i).alpha_center(1,2),temp_table_C1(i).beta_center(1,1),temp_table_C1(i).beta_center(1,2),...
                temp_table_C1(i).theta_avg_power(1,1),temp_table_C1(i).alpha_avg_power(1,1),temp_table_C1(i).beta_avg_power(1,1),...
                'VariableNames',...
                {'subID','compID','study','cond','group','cluster','aperiodic_exp','aperiodic_offset','r_squared','alpha_cf','alpha_p',...
                'beta_cf','beta_p','alpha_center','alpha_centerP','beta_center','beta_centerP',...
                'theta_avg_power','alpha_avg_power','beta_avg_power'}));
        end
    end
end
C1_table.subID = categorical(C1_table.subID);
C1_table.cond = categorical(C1_table.cond);
C1_table.cluster = categorical(C1_table.cluster);
C1_table.study = categorical(C1_table.study);

% C1_table.study = categorical(C1_table.study);

grp_C1_table = grpstats(C1_table,["subID","study","cond","cluster"],'nanmedian','DataVars',["r_squared","alpha_cf","alpha_p",...
            "beta_cf","beta_p","alpha_center","alpha_centerP","beta_center","beta_centerP","theta_avg_power","alpha_avg_power","beta_avg_power"]);
grp_C1_table.Properties.VariableNames = {'subID','study','cond','cluster','GroupCount','med_r_squared','med_alpha_cf','med_alpha_p',...
            'med_beta_cf','med_beta_p','med_alpha_center','med_alpha_centerP','med_beta_center','med_beta_centerP',...
            'med_theta_avg_power','med_alpha_avg_power','med_beta_avg_power'};

%% Sanity check: Plot distribution of aperiodic params (exp), central frequency, and goodness of fit
g = 1;
for k = 3%3:length(STUDY_imagined.cluster)
    temp_table = C1_table(C1_table.cluster == num2str(k) & C1_table.study == num2str(g),: );
    temp_table.cond = categorical(temp_table.cond);
    figure();set(gcf,'color','white');
    subplot(2,2,1)
    boxchart(temp_table.cond,temp_table.aperiodic_exp);
    ylabel('Aperodic exponent');

%     subplot(2,2,2)
%     boxchart(temp_table.cond,temp_table.center_frequency);
%     xlabel('Central frequency');ylabel('# Peaks')

    subplot(2,2,3)
    boxchart(temp_table.cond,temp_table.r_squared);
    ylabel('R squared');

end
%% --------------------
%% STATS Section: 
%% Preliminary stats on average power (substracted background)
% 
% Create stats table
T_stats = table;
for g = 1:2
    for k = 3:length(fooof_diff_store{g})
   
        t1 = C1_table(C1_table.cluster == num2str(k) & C1_table.study == num2str(g) ,:);
        t1.cond = categorical(t1.cond);
        t1.log_alpha_avg_power = log(t1.alpha_avg_power+5);
        t1.log_theta_avg_power = log(t1.theta_avg_power+5);
        % t2 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,'2')& grp_C1_table.cond ~= 1,:);
        % t2.cond = categorical(t2.cond);

        lme_theta_avg_power = fitlme(t1,'theta_avg_power ~ cond + (1|subID)');
        Th = anova(lme_theta_avg_power);[h,p] = lillietest(lme_theta_avg_power.residuals)
        lme_alpha_avg_power= fitlme(t1,'alpha_avg_power ~ cond + (1|subID)');
        A = anova(lme_alpha_avg_power);[h,p] = lillietest(lme_alpha_avg_power.residuals)
        lme_beta_avg_power = fitlme(t1,'beta_avg_power ~ cond + (1|subID)');
        B = anova(lme_beta_avg_power);[h,p] = lillietest(lme_beta_avg_power.residuals)
%         lme_log_theta_avg_power = fitlme(t1,'log_theta_avg_power ~ cond + (1|subID)');
%         log_Th = anova(lme_log_theta_avg_power);[h,p] = lillietest(lme_log_theta_avg_power.residuals)
%         lme_log_alpha_avg_power= fitlme(t1,'log_alpha_avg_power ~ cond + (1|subID)');
%         log_A = anova(lme_log_alpha_avg_power);[h,p] = lillietest(lme_log_alpha_avg_power.residuals)
        
        temp_stats = table;
        temp_stats.study = g;
        temp_stats.cluster = k;
        temp_stats.theta_anova = Th.pValue(2);
        temp_stats.theta_F = Th.FStat(2);
        temp_stats.theta_F_DF2 = Th.DF2(2);
        temp_stats.alpha_anova = A.pValue(2);
        temp_stats.alpha_F = A.FStat(2);
        temp_stats.alpha_F_DF2 = A.DF2(2);
        temp_stats.beta_anova = B.pValue(2);
        temp_stats.beta_F = B.FStat(2);
        temp_stats.beta_F_DF2 = B.DF2(2);
        temp_stats.theta_cond2 = lme_theta_avg_power.Coefficients.pValue(2);
        temp_stats.theta_cond3 = lme_theta_avg_power.Coefficients.pValue(3);
        temp_stats.theta_cond4 = lme_theta_avg_power.Coefficients.pValue(4);

        temp_stats.alpha_cond2 = lme_alpha_avg_power.Coefficients.pValue(2);
        temp_stats.alpha_cond3 = lme_alpha_avg_power.Coefficients.pValue(3);
        temp_stats.alpha_cond4 = lme_alpha_avg_power.Coefficients.pValue(4);

        temp_stats.beta_cond2 = lme_beta_avg_power.Coefficients.pValue(2);
        temp_stats.beta_cond3 = lme_beta_avg_power.Coefficients.pValue(3);
        temp_stats.beta_cond4 = lme_beta_avg_power.Coefficients.pValue(4);

        % - add stats for aperiod fit/offset
        lme_ap_exp = fitlme(t1,'aperiodic_exp ~ cond + (1|subID)');
        ap_exp = anova(lme_ap_exp);
        temp_stats.ap_exp_anova = ap_exp.pValue(2);
        temp_stats.ap_exp2 = lme_ap_exp.Coefficients.pValue(2);
        temp_stats.ap_exp3 = lme_ap_exp.Coefficients.pValue(3);
        temp_stats.ap_exp4 = lme_ap_exp.Coefficients.pValue(4);
        
        lme_ap_offset = fitlme(t1,'aperiodic_offset ~ cond + (1|subID)');
        ap_offset = anova(lme_ap_offset);
        temp_stats.ap_offset_anova = ap_offset.pValue(2);
        temp_stats.ap_offset2 = lme_ap_offset.Coefficients.pValue(2);
        temp_stats.ap_offset3 = lme_ap_offset.Coefficients.pValue(3);
        temp_stats.ap_offset4 = lme_ap_offset.Coefficients.pValue(4);
        
       
        if g == 2
            % use continuous variable
            t1.cond = double(t1.cond)*0.25;
            lme_theta_avg_power_num = fitlme(t1,'theta_avg_power ~ cond + (1|subID)');
            Th_num = anova(lme_theta_avg_power_num);[h,p] = lillietest(lme_theta_avg_power_num.residuals)
            lme_alpha_avg_power_num = fitlme(t1,'alpha_avg_power ~ cond  + (1|subID)');
            A_num  = anova(lme_alpha_avg_power_num);[h,p] = lillietest(lme_alpha_avg_power_num.residuals)
            lme_beta_avg_power_num = fitlme(t1,'beta_avg_power ~ cond + (1|subID)');
            B_num  = anova(lme_beta_avg_power_num);[h,p] = lillietest(lme_beta_avg_power_num.residuals)
            
            temp_stats.Th_num = Th_num.pValue(2);
            temp_stats.A_num = A_num.pValue(2);
            temp_stats.B_num = B_num.pValue(2);
            temp_stats.theta_F_num = Th_num.FStat(2);
            temp_stats.theta_F_DF2_num = Th_num.DF2(2);
            temp_stats.alpha_F_num = A_num.FStat(2);
            temp_stats.alpha_F_DF2_num = A_num.DF2(2);
            temp_stats.beta_F_num = B_num.FStat(2);
            temp_stats.beta_F_DF2_num = B_num.DF2(2);
            
            % store the linear fit for all these data
            temp_stats.Th_slope = lme_theta_avg_power_num.Coefficients.Estimate(2);
            temp_stats.Th_intercept = lme_theta_avg_power_num.Coefficients.Estimate(1);
            temp_stats.A_slope = lme_alpha_avg_power_num.Coefficients.Estimate(2);
            temp_stats.A_intercept = lme_alpha_avg_power_num.Coefficients.Estimate(1);
            temp_stats.B_slope = lme_beta_avg_power_num.Coefficients.Estimate(2);
            temp_stats.B_intercept = lme_beta_avg_power_num.Coefficients.Estimate(1);
            temp_stats.Th_num_R2 = lme_theta_avg_power_num.Rsquared.Adjusted;
            temp_stats.A_num_R2 = lme_alpha_avg_power_num.Rsquared.Adjusted;
            temp_stats.B_num_R2 = lme_beta_avg_power_num.Rsquared.Adjusted;
        else 
            temp_stats.Th_num = NaN;
            temp_stats.A_num = NaN;
            temp_stats.B_num = NaN;
            temp_stats.theta_F_num = NaN;
            temp_stats.theta_F_DF2_num = NaN;
            temp_stats.alpha_F_num = NaN;
            temp_stats.alpha_F_DF2_num = NaN;
            temp_stats.beta_F_num = NaN;
            temp_stats.beta_F_DF2_num = NaN;
            
            temp_stats.Th_slope = NaN;
            temp_stats.Th_intercept = NaN;
            temp_stats.A_slope = NaN;
            temp_stats.A_intercept = NaN;
            temp_stats.B_slope = NaN;
            temp_stats.B_intercept = NaN;
            
            temp_stats.Th_num_R2 = NaN;
            temp_stats.A_num_R2 = NaN;
            temp_stats.B_num_R2 = NaN;
        end
        
        T_stats = vertcat(T_stats,temp_stats);
    end
end
T_stats.study = categorical(T_stats.study);
T_stats.cluster = categorical(T_stats.cluster);

%% Perform time series stats on the flattened curve
iter = 200; % in eeglab, the fdr stats will automatically *20
try
    STUDY.etc = rmfield(STUDY.etc,'statistics')
end
% STUDY = pop_statparams(STUDY,'groupstats','on','condstats','on','statistics','perm',...
%     'singletrials','off','mode','eeglab','effect','main','alpha',NaN,'mcorrect','fdr','naccu',iter);% If not using mcorrect, use none, Not sure why, if using fdr correction, none of these are significant
% 
STUDY = pop_statparams(STUDY, 'groupstats','off','condstats', 'on',...
            'method','perm',...
            'singletrials','off','mode','fieldtrip','fieldtripalpha',NaN,...
            'fieldtripmethod','montecarlo','fieldtripmcorrect','fdr','fieldtripnaccu',iter*20);

stats = STUDY.etc.statistics;
stats.paired{1} = 'on'; % Condition stats
stats.paired{2} = 'off'; % Group stats

% -----------------
% fooof_diff_store needs to be freq x subject, and condition by row
for g = 1:2
    for k = 3:length(fooof_diff_store{g})
        k
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(fooof_diff_store{g}{k}', stats);
        pcond{g}{k} = temp_pcond;
        pgroup{g}{k} = temp_pcond;
        pinter{g}{k} = temp_pinter;
        statcond{g}{k} = temp_statcond;
        statgroup{g}{k} = temp_statgroup;
        statinter{g}{k} = temp_statinter;
        for k0 = 1:length(pcond{g}{k})
            pcond{g}{k}{k0}(:,2) = pcond{g}{k}{k0}(:,1)<0.05;            
        end
        for k0 = 1:length(pgroup{g}{k})
            if ~isempty(pgroup{g}{k}{k0})
                pgroup{g}{k}{k0}(:,2) = pgroup{g}{k}{k0}(:,1)<0.05;  
            end
        end
        for k0 = 1:length(pinter{g}{k})
            if ~isempty(pinter{g}{k}{k0})
                pinter{g}{k}{k0}(:,2) = pinter{g}{k}{k0}(:,1)<0.05;
            end
        end            
    end
end

for g = 1:2
    for k = 3:length(spec_data_original{g})
        k
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(spec_data_original{g}{k}', stats);
        pcond_org{g}{k} = temp_pcond;
        pgroup_org{g}{k} = temp_pcond;
        pinter_org{g}{k} = temp_pinter;
        statcond_org{g}{k} = temp_statcond;
        statgroup_org{g}{k} = temp_statgroup;
        statinter_org{g}{k} = temp_statinter;
        for k0 = 1:length(pcond_org{g}{k})
            pcond_org{g}{k}{k0}(:,2) = pcond_org{g}{k}{k0}(:,1)<0.05;            
        end
        for k0 = 1:length(pgroup_org{g}{k})
            if ~isempty(pgroup_org{g}{k}{k0})
                pgroup_org{g}{k}{k0}(:,2) = pgroup_org{g}{k}{k0}(:,1)<0.05;  
            end
        end
        for k0 = 1:length(pinter_org{g}{k})
            if ~isempty(pinter_org{g}{k}{k0})
                pinter_org{g}{k}{k0}(:,2) = pinter_org{g}{k}{k0}(:,1)<0.05;
            end
        end            
    end
end

%% Paper Figures 
%% Figure after flattened the curve ------------------------------------------------
g = 1;
switch g
    case 1
        color_dark = color.terrain(2:end,:);
        color_light = color.terrain_shade(2:end,:);
    case 2
        color_dark = color.speed;
        color_light = color.speed_shade;
end
% %- colors
% COLORS_MAPS = linspecer(size(speed_chars,2));
% custom_yellow = [254,223,0]/255;
% COLORS_MAPS = [COLORS_MAPS(3,:);custom_yellow;COLORS_MAPS(4,:);COLORS_MAPS(2,:)];
SUBPLOT_HEIGHT = 4;
SUBPLOT_WIDTH = 5;
figure('color','white','position',[100 300 800 800]);
for k = 3:length(fooof_group_results_org{g})
    subplot(SUBPLOT_HEIGHT,SUBPLOT_WIDTH,k)
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
            color_dark(i,:),color_light(i,:));% input need to be a row vector
    end
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2);
    end  
   
    ax = gca;       
    axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
    xlim([4 40]);
    plot([0 40],[0 0],'--','color','black');
    xlabel('Frequency(Hz)');ylabel('10*log10(Power)');
    set(gca,'fontsize',10);
    title(['Cluster ',num2str(k)]);
    xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
end

%% Figure Plot original spec_data and the flattened data side by side ------------------------------------------------
% Updated 2023-06-08 as the official figures
g = 2;
switch g
    case 1
        color_dark = color.terrain(2:end,:);
        color_light = color.terrain_shade(2:end,:);
    case 2
        color_dark = color.speed;
        color_light = color.speed_shade;
end

% figure('color','white','position',[100 300 600 1000],'Renderer','Painters');
j = 1;
for k = 3:length(fooof_group_results_org{g})
    % hardcode to preset the axis limit
    switch k
        case {3,12} % sensorimotor area
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-2 6.5];
        case {7,9} % posterior area
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-2 8];
        case 14 % cingulate
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-1 4];
        case {6,13} % supplementary motor
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-1 4];
        case {11,10} % occipital
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-2 8];
        case 4 % caudate
            ylim_value_origPSD = [-30 -10];
            ylim_value_flatPSD = [-1 4];
    end

    if mod(j,2) == 1 
        figure('color','white','position',[100 300 450 210],'Renderer','Painters');
        j = 1;
    end

    % -------------- original PSD -----------
    subplot(1,2,j)
    title('Original PSD')
    for i = 1:length(spec_data_original{g}{k})
        data = spec_data_original{g}{k}{i}';
        JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
            color_dark(i,:),color_light(i,:));% input need to be a row vector
    end
    for i = 1:length(spec_data_original{g}{k})
        data = spec_data_original{g}{k}{i}';
        plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
    end  
    % plot the aperiodic line
    for i = 1:length(spec_data_original{g}{k})
        aperiodic_fit = fooof_apfit_store{g}{k}{i}';
        plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',0.5);
    end
    ylim(ylim_value_origPSD);
    ax = gca;       
    axsignif = highlight_CL(ax, fooof_freq, pcond_org{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
    xlim([4 40]);
    
    plot([0 40],[0 0],'--','color','black');
    xlabel('Frequency(Hz)');
    ylabel('10*log(PSD)');
    set(gca,'fontsize',10);
%     title(['Cluster ',num2str(k)]);
    xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    
    % -------------- Flattened PSD -----------
    subplot(1,2,j+1)
    title('Flattened PSD');
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
            color_dark(i,:),color_light(i,:));% input need to be a row vector
    end
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
    end  
    ylim(ylim_value_flatPSD);
    ax = gca;       
    axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
    xlim([4 40]);
    plot([0 40],[0 0],'-','color',[0.5 0.5 0.5]);
    ylabel('');
    set(gca,'fontsize',10);
    
    xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    
    j = j + 2;
    outputdir_psd = [save_path filesep 'PSD']; %'D:\Dropbox (UFL)\0.0 Writing_postdoc\MiM_YA_paper\Figure\PSD';
    if ~exist(outputdir_psd,'dir')
        mkdir(outputdir_psd)
    end
    if g == 1
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_',num2str(k),'.pdf']),'ContentType','vector')
    elseif g == 2
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_',num2str(k),'_speed.pdf']),'ContentType','vector')
    end
end

%% Paper Figure: Violin plot for the average theta alpha betas power
measure_name_plot = {'theta_avg_power','alpha_avg_power','beta_avg_power'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'\theta','\alpha','\beta'};
% figure('color','white','position',[200 200 800 400]);
T_plot = table;
g = 1;
switch g
    case 1
        color_dark = color.terrain(2:end,:);
        color_light = color.terrain_shade(2:end,:);
        xtick_label_g = {'flat','low','med','high'};
    case 2
        color_dark = color.speed;
        color_light = color.speed_shade;
        xtick_label_g = {'0.25','0.5','0.75','1.0'};
end

%----------------------- PLOT ## By subject plot 
% figure('color','white','position',[100 300 600 180]);
j = 1;
for k = 4%3:length(fooof_diff_store{g}) % by each cluster
    switch k
        case {3,12} % sensorimotor area
            ylim_value_theta = [-1 2.5];
            ylim_value_alpha = [-2 18];
            ylim_value_beta  = [-1 10];
        case {7,9} % posterior area
            ylim_value_theta = [-1 2.5];
            ylim_value_alpha = [-2 12];
            ylim_value_beta  = [-1 8];
        case 14 % cingulate
            ylim_value_theta = [-1 6];
            ylim_value_alpha = [-2 10];
            ylim_value_beta  = [-1 7];
        case {6,13} % supplementary motor
            ylim_value_theta = [-0.5 6];
            ylim_value_alpha = [-2 4];
            ylim_value_beta  = [-2 6];
        case {11,10} % occipital
            ylim_value_theta = [-1 4];
            ylim_value_alpha = [-2 18];
            ylim_value_beta  = [-1 6.5];
        case 4 % caudate
            ylim_value_theta = [-1 4];
            ylim_value_alpha = [-2 7];
            ylim_value_beta  = [-0.1 4];
    end
    for i = 1:length(measure_name_plot)
        if mod(j,3) == 1 
            figure('color','white','position',[100 300 600 200]);
            j = 1;
        end
        measure_name = measure_name_plot{i};

        T_plot = C1_table(C1_table.study == num2str(g) & C1_table.cluster == num2str(k),:); % categorical variable can use ==

        %----------------------- 
        subplot(1,3,j)
        hold on;
        violinplot(T_plot.(measure_name),T_plot.cond,...
            'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
            'MarkerSize',15,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin','full','width',0.3,'QuartileStyle','shadow')
        if j == 1;ylabel('10*log(Flattened PSD)');else;ylabel('');end
        if g == 2;xlabel('m/s');end
        
        % axis padded
        fig_i = get(groot,'CurrentFigure');
        % fig_i.Position = [200,200,1820,920];
        box off
        title(title_plot{i});
        set(gca,'xticklabel', xtick_label_g,'fontsize',10);
        
        j = j + 1;
        T_stats_plot = T_stats(T_stats.study == num2str(g) & T_stats.cluster == num2str(k),:);
        switch i 
            case 1
                anova_stats = T_stats_plot.theta_anova;
                cond2_stats = T_stats_plot.theta_cond2;
                cond3_stats = T_stats_plot.theta_cond3;
                cond4_stats = T_stats_plot.theta_cond4;
                regress_sig = T_stats_plot.Th_num;
                regressline_stats = [T_stats_plot.Th_intercept T_stats_plot.Th_slope]; 
                R2 = T_stats_plot.Th_num_R2;
            case 2
                anova_stats = T_stats_plot.alpha_anova;
                cond2_stats = T_stats_plot.alpha_cond2;
                cond3_stats = T_stats_plot.alpha_cond3;
                cond4_stats = T_stats_plot.alpha_cond4;
                regress_sig = T_stats_plot.A_num;
                regressline_stats = [T_stats_plot.A_intercept T_stats_plot.A_slope];
                R2 = T_stats_plot.A_num_R2;
            case 3
                anova_stats = T_stats_plot.beta_anova;
                cond2_stats = T_stats_plot.beta_cond2;
                cond3_stats = T_stats_plot.beta_cond3;
                cond4_stats = T_stats_plot.beta_cond4;
                regress_sig = T_stats_plot.B_num;
                regressline_stats = [T_stats_plot.B_intercept T_stats_plot.B_slope];
                R2 = T_stats_plot.B_num_R2;
        end
        if g == 1 % terrain condition, categorical values
             if anova_stats < 0.05
                if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
                if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
                if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
             end
        elseif g == 2 % g == 2 % speed conditions, numeric values, regression line
            if anova_stats < 0.05
                % plot line
                x = 0:5;
                y = x*regressline_stats(2)*0.25 + regressline_stats(1);
                plot(x,y,'-','color','k','linewidth',1);
                xlim([0 5]);
                if regress_sig > 0.01 & regress_sig < 0.05
                    text(1,gety(gca)*1.2,{['* ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]});
                elseif regress_sig <= 0.01 & regress_sig > 0.001 
                    text(1,gety(gca)*1.2,{['** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]});
                else
                    text(1,gety(gca)*1.2,{['*** ','slope = ',num2str(round(regressline_stats(2),2))],['R^2 = ',num2str(round(R2,2))]});
                end
            end
        end        
        
         clear ydata
         ydata = [T_plot.alpha_avg_power;T_plot.beta_avg_power];
         max(ydata)
         switch i 
             case 1
                 ylim(ylim_value_theta);
             case 2
                 ylim(ylim_value_alpha);
             case 3
                 ylim(ylim_value_beta);
         end
%          ylim([-1 max(ydata)*1.5]);
    
%         ylim([0 Inf]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',measure_name)]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',measure_name)]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.jpg',measure_name)]);
    end
    if g == 1
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_Avg_Power',num2str(k),'.pdf']),'ContentType','vector')
    elseif g == 2
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_Avg_Power',num2str(k),'_speed_reg.pdf']),'ContentType','vector')
    end
end

%% Paper Figure: Violin plot for the Aperiodic fit slope and offset
measure_name_plot = {'aperiodic_exp','aperiodic_offset'}; % walking speed, stride duration, step variability, sacrum excursion variability for ML and AP
title_plot = {'aperiodic exponent','apriodic offset'};
% figure('color','white','position',[200 200 800 400]);
T_plot = table;
g = 2;
switch g
    case 1
        color_dark = color.terrain(2:end,:);
        color_light = color.terrain_shade(2:end,:);
        xtick_label_g = {'flat','low','med','high'};
    case 2
        color_dark = color.speed;
        color_light = color.speed_shade;
        xtick_label_g = {'0.25','0.5','0.75','1.0'};
end

%----------------------- PLOT ## By subject plot 
% figure('color','white','position',[100 300 600 180]);
j = 1;
for k = 6%[3:length(fooof_diff_store{g})]%3:length(fooof_diff_store{g}) % by each cluster
    switch k
        case {3,12} % sensorimotor area
            ylim_value_ap_exp = [0 2];
            ylim_value_ap_offset = [-3 1]*10;
        case {7,9} % posterior area
            ylim_value_ap_exp = [0 2];
            ylim_value_ap_offset = [-3 1]*10;
        case 14 % cingulate
            ylim_value_ap_exp = [0 2];
            ylim_value_ap_offset = [-3 1]*10;
        case {6,13} % supplementary motor
            ylim_value_ap_exp = [0 2];
            ylim_value_ap_offset = [-3 1]*10;
        case {11,10} % occipital
            ylim_value_ap_exp = [0 2];
            ylim_value_ap_offset = [-3 1]*10;
    end
    for i = 1:length(measure_name_plot)
        if mod(j,2) == 1 
            figure('color','white','position',[100 300 400 200]);
            j = 1;
        end
        measure_name = measure_name_plot{i};

        T_plot = C1_table(C1_table.study == num2str(g) & C1_table.cluster == num2str(k),:); % categorical variable can use ==

        %----------------------- 
        subplot(1,2,j)
        hold on;
        
        violinplot(T_plot.(measure_name)*10^(i-1),T_plot.cond,...
            'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
            'MarkerSize',15,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin','full','width',0.3,'QuartileStyle','shadow')
%         if j == 1;ylabel('10*log(Flattened PSD)');else;ylabel('');end
        if g == 2;xlabel('m/s');end
        hold off;
        % axis padded
        fig_i = get(groot,'CurrentFigure');
        % fig_i.Position = [200,200,1820,920];
        box off
        title(title_plot{i});
        set(gca,'xticklabel', xtick_label_g,'fontsize',10);
        
        switch i 
             case 1
                 ylim(ylim_value_ap_exp);
             case 2
                 ylim(ylim_value_ap_offset);
        end
        j = j + 1;
        T_stats_plot = T_stats(T_stats.study == num2str(g) & T_stats.cluster == num2str(k),:);
        switch i 
            case 1
                anova_stats = T_stats_plot.ap_exp_anova;
                cond2_stats = T_stats_plot.ap_exp2;
                cond3_stats = T_stats_plot.ap_exp3;
                cond4_stats = T_stats_plot.ap_exp4;
            case 2
                anova_stats = T_stats_plot.ap_offset_anova;
                cond2_stats = T_stats_plot.ap_offset2;
                cond3_stats = T_stats_plot.ap_offset3;
                cond4_stats = T_stats_plot.ap_offset4;           
        end
         if anova_stats < 0.05
            if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
            if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
            if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
         end
         clear ydata
         
         
%          ylim([-1 max(ydata)*1.5]);
    
%         ylim([0 Inf]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.fig',measure_name)]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Subjects_Fig_%s.jpg',measure_name)]);
    % saveas(fig_i,[save_dir filesep sprintf('Across_Trials_Fig_%s.jpg',measure_name)]);
    end
    if g == 1
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_AP_fit',num2str(k),'.pdf']),'ContentType','vector')
    elseif g == 2
        exportgraphics(gcf,fullfile(outputdir_psd ,['Cluster_PSD_AP_ft',num2str(k),'_speed.pdf']),'ContentType','vector')
    end
end

%% Make giant plot, with topography
g = 2;
switch g
    case 1
        color_dark = color.terrain(2:end,:);
        color_light = color.terrain_shade(2:end,:);
        xtick_label_g = {'flat','low','med','high'};
    case 2
        color_dark = color.speed;
        color_light = color.speed_shade;
        xtick_label_g = {'0p25','0p5','0p75','1p0'};
end
figure('color','white','position',[100 300 1200 1000]);
j = 1;
for k = [3 12 7 9 6 13 11 10 14] %3:length(fooof_group_results_org{g})
    if j == 31 
        figure('color','white','position',[100 300 1200 1000]);
        j = 1;
    end

    subplot(5,6,j)
    std_topoplot_CL(STUDY,k,'together');
    colormap(colormap_ersp)
    set(gcf,'color','w')

    subplot(5,6,j+1)
    for i = 1:length(spec_data_original{g}{k})
        data = spec_data_original{g}{k}{i}';
        JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
            color_dark(i,:),color_light(i,:));% input need to be a row vector
    end
    for i = 1:length(spec_data_original{g}{k})
        data = spec_data_original{g}{k}{i}';
        plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2);
    end  
   
    ax = gca;       
    axsignif = highlight_CL(ax, fooof_freq, pcond_org{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
    xlim([4 40]);ylim([-30 -5]);
    plot([0 40],[0 0],'--','color','black');
    xlabel('Frequency(Hz)');ylabel('10*log10(Power)');
    set(gca,'fontsize',10);
%     title(['Cluster ',num2str(k)]);
    xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    
    % ----------------------
    subplot(5,6,j+2)
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
            color_dark(i,:),color_light(i,:));% input need to be a row vector
    end
    for i = 1:length(fooof_diff_store{g}{k})
        data = fooof_diff_store{g}{k}{i}';
        plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',2);
    end  
   
    ax = gca;       
    axsignif = highlight_CL(ax, fooof_freq, pcond{g}{k}{1}(:,2), 'background', 'Frequency(Hz)');
    xlim([4 40]);
    plot([0 40],[0 0],'--','color','black');
    xlabel('Frequency(Hz)');ylabel('10*log10(Power)');
    set(gca,'fontsize',10);
%     title(['Cluster ',num2str(k)]);
    xline([3],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
    
    % -------------------
    for i = 1:length(measure_name_plot)
        measure_name = measure_name_plot{i};
        T_plot = C1_table(C1_table.study == num2str(g) & C1_table.cluster == num2str(k),:); % categorical variable can use ==

        subplot(5,6,j+i+2)
        hold on;
        violinplot(T_plot.(measure_name),T_plot.cond,...
            'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
            'MarkerSize',15,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter','HalfViolin',...
            'full','width',0.3,'QuartileStyle','shadow');
        ylabel('');
        xlabel('');
        hold off;
        % axis padded
        fig_i = get(groot,'CurrentFigure');
        % fig_i.Position = [200,200,1820,920];
        box off
        title(title_plot{i});
        set(gca,'xticklabel', xtick_label_g,'fontsize',10);
        T_stats_plot = T_stats(T_stats.study == num2str(g) & T_stats.cluster == num2str(k),:);
        switch i 
            case 1
                anova_stats = T_stats_plot.theta_anova;
                cond2_stats = T_stats_plot.theta_cond2;
                cond3_stats = T_stats_plot.theta_cond3;
                cond4_stats = T_stats_plot.theta_cond4;
            case 2
                anova_stats = T_stats_plot.alpha_anova;
                cond2_stats = T_stats_plot.alpha_cond2;
                cond3_stats = T_stats_plot.alpha_cond3;
                cond4_stats = T_stats_plot.alpha_cond4;
            case 3
                anova_stats = T_stats_plot.beta_anova;
                cond2_stats = T_stats_plot.beta_cond2;
                cond3_stats = T_stats_plot.beta_cond3;
                cond4_stats = T_stats_plot.beta_cond4;
        end
         if anova_stats < 0.05
            if cond2_stats < 0.05;sigline([1 2],'*',[],[],cond2_stats);end
            if cond3_stats < 0.05;sigline([1 3],'*',[],[],cond3_stats);end
            if cond4_stats < 0.05;sigline([1 4],'*',[],[],cond4_stats);end
         end
         clear ydata
         ydata = [T_plot.alpha_avg_power;T_plot.beta_avg_power];
         max(ydata)
%          ylim([-1 max(ydata)*1.5]);
    end
    
    j = j + 6;
end
%% 

study = 1;
for g = 2
    for k = 3:length(fooof_group_results_org{g})
        temp_table = fooof_group_results_org{g}{k}([fooof_group_results_org{g}{k}(:).study] == g);

        temp_table_C1 = temp_table([temp_table(:).cond] == 1);
        temp_table_C2 = temp_table([temp_table(:).cond] == 2);
        temp_table_C3 = temp_table([temp_table(:).cond] == 3);
        temp_table_C4 = temp_table([temp_table(:).cond] == 4);

        figure();set(gcf,'color','white','position',[200 200 1000 500]);
        subplot(2,3,1);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).theta_avg_power));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.Green);
        a2 = 10*(vertcat(temp_table_C2(:).theta_avg_power));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Yellow);
        a3 = 10*(vertcat(temp_table_C3(:).theta_avg_power));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Orange);
        a4 = 10*(vertcat(temp_table_C4(:).theta_avg_power));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Red);
        boxplot([a1 a2 a3 a4],'Colors','k');
        plot([1:4],[median(a1) median(a2) median(a3) median(a4)],'-','color','black','linewidth',1.5);
        ylabel('Theta Power');
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'flat','low','med','high'},'fontsize',10);

        subplot(2,3,2);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).alpha_avg_power));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.Green);
        a2 = 10*(vertcat(temp_table_C2(:).alpha_avg_power));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Yellow);
        a3 = 10*(vertcat(temp_table_C3(:).alpha_avg_power));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Orange);
        a4 = 10*(vertcat(temp_table_C4(:).alpha_avg_power));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Red);
        boxplot([a1 a2 a3 a4],'Colors','k');
        plot([1:4],[median(a1) median(a2) median(a3) median(a4) ],'-','color','black','linewidth',1.5);
        ylabel('Alpha Power');
        title(['Cluster ', num2str(k) , ' ',study]);
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'flat','low','med','high'},'fontsize',10);

        subplot(2,3,3);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).beta_avg_power));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.Green);
        a2 = 10*(vertcat(temp_table_C2(:).beta_avg_power));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Yellow);
        a3 = 10*(vertcat(temp_table_C3(:).beta_avg_power));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Orange);
        a4 = 10*(vertcat(temp_table_C4(:).beta_avg_power));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Red);
        boxplot([a1 a2 a3 a4 a5],'Colors','k');
        plot([1:4],[median(a1) median(a2) median(a3) median(a4)],'-','color','black','linewidth',1.5);
        ylabel('Beta Power');
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'flat','low','med','high'},'fontsize',10);

    end
end

%% aperiodic_exp, central_freq,r_squared across groups
% Looks like the peaks occur during wide range... not sure if the spectras
% are okay
for k = 4
    temp_table = fooof_group_results_org{k}([fooof_group_results_org{k}(:).cond] == 1);
    temp_table_H1 = temp_table([temp_table(:).group] == 1);
    temp_table_H2 = temp_table([temp_table(:).group] == 2);
    
    figure();set(gcf,'color','white');
    subplot(2,2,1);hold on;
    plot(ones(size(temp_table_H1,2),1),[temp_table_H1(:).aperiodic_exp],'bo');
    plot(ones(size(temp_table_H2,2),1),[temp_table_H2(:).aperiodic_exp],'ro');
    ylabel('Aperodic exponent');

    subplot(2,2,2);hold on;
    histogram(vertcat(temp_table_H1(:).central_freq));
    histogram(vertcat(temp_table_H2(:).central_freq));
    xlabel('Central frequency');ylabel('# Peaks')

    subplot(2,2,3);hold on;
    plot(ones(size(temp_table_H1,2),1),[temp_table_H1(:).r_squared],'bo');
    plot(ones(size(temp_table_H2,2),1),[temp_table_H2(:).r_squared],'ro');
    ylabel('R squared');

   
end

%%  Plot peaks occur within a specific band comparing young and old
figure();set(gcf,'color','white','position',[100 200 1000 1000]);
for k = 3:length(fooof_group_results_org)
    temp_table = fooof_group_results_org{k}([fooof_group_results_org{k}(:).cond] == 1 & strcmp({fooof_group_results_org{k}(:).study},'Imagined'));
    temp_table_H1 = temp_table([temp_table(:).group] == 1);
    temp_table_H2 = temp_table([temp_table(:).group] == 2);
    
    subplot(5,3,1+(k-3)*3);hold on;
    a1 = vertcat(temp_table_H1(:).alpha);
    if ~isempty(a1);plot(ones(length(a1),1),a1(:,2),'bo');end
    a2 = vertcat(temp_table_H2(:).alpha);
    if ~isempty(a2);plot(ones(length(a2),1)*2,a2(:,2),'ro');end
    ylabel('Alpha peak power');
    xlim([0 3]);
    set(gca,'xtick',linspace(1,2,2),'xticklabel',{'H1','H2'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    subplot(5,3,2+(k-3)*3);hold on;
    b1 = vertcat(temp_table_H1(:).beta);
    if ~isempty(b1);plot(ones(length(b1),1),b1(:,2),'bo');end
    b2 = vertcat(temp_table_H2(:).beta);
    if ~isempty(b2);plot(ones(length(b2),1)*2,b2(:,2),'ro');end
    ylabel('Beta peak power');
    xlim([0 3]);
    set(gca,'xtick',linspace(1,2,2),'xticklabel',{'H1','H2'},'fontsize',12);
    title(['Cluster ', num2str(k)]);
    
    subplot(5,3,3+(k-3)*3);hold on;
    t1 = vertcat(temp_table_H1(:).theta);
    if ~isempty(t1)
        plot(ones(length(t1),1),t1(:,2),'bo');
    end
    t2 = vertcat(temp_table_H2(:).theta);
    if ~isempty(t2)
        plot(ones(length(t2),1)*2,t2(:,2),'ro');
    end
    ylabel('Theta peak power');
    xlim([0 3]);
    set(gca,'xtick',linspace(1,2,2),'xticklabel',{'H1','H2'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

end

%%  Plot all components peaks occur within a specific band comparing young adults across conditions
% 
figure();set(gcf,'color','white','position',[100 300 1000 1000]);
for k = 3:7%length(fooof_group_results_org)
    temp_table = fooof_group_results_org{k}( [fooof_group_results_org{k}(:).group] == 1 & strcmp({fooof_group_results_org{k}(:).study},'Imagined'));
    temp_table_C1 = temp_table([temp_table(:).cond] == 1);
    temp_table_C2 = temp_table([temp_table(:).cond] == 2);
    temp_table_C3 = temp_table([temp_table(:).cond] == 3);
    temp_table_C4 = temp_table([temp_table(:).cond] == 4);
    temp_table_C5 = temp_table([temp_table(:).cond] == 5);
    
    if k <= 7
        
        subplot(5,3,1+(k-3)*3);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).theta));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
        a2 = 10*(vertcat(temp_table_C2(:).theta));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
        a3 = 10*(vertcat(temp_table_C3(:).theta));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
        a4 = 10*(vertcat(temp_table_C4(:).theta));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
        a5 = 10*(vertcat(temp_table_C5(:).theta));
        fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
        ylabel('Theta peak power');
        xlim([0 6]);
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
        title(['Cluster ', num2str(k)]);

        subplot(5,3,2+(k-3)*3);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).alpha));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
        a2 = 10*(vertcat(temp_table_C2(:).alpha));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
        a3 = 10*(vertcat(temp_table_C3(:).alpha));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
        a4 = 10*(vertcat(temp_table_C4(:).alpha));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
        a5 = 10*(vertcat(temp_table_C5(:).alpha));
        fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);

        ylabel('Alpha peak power');
        xlim([0 6]);
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
        title(['Cluster ', num2str(k)]);

        subplot(5,3,3+(k-3)*3);hold on;
        a1 = 10*(vertcat(temp_table_C1(:).beta));
        fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
        a2 = 10*(vertcat(temp_table_C2(:).beta));
        fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
        a3 = 10*(vertcat(temp_table_C3(:).beta));
        fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
        a4 = 10*(vertcat(temp_table_C4(:).beta));
        fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
        a5 = 10*(vertcat(temp_table_C5(:).beta));
        fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
        ylabel('Beta peak power');
        xlim([0 6]);
        set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
        title(['Cluster ', num2str(k)]);

    end
end

figure();set(gcf,'color','white','position',[100 300 1000 1000]);
i = 1;
for k = 8:length(fooof_group_results_org)
    temp_table = fooof_group_results_org{k}( [fooof_group_results_org{k}(:).group] == 1 & strcmp({fooof_group_results_org{k}(:).study},'Imagined'));
    temp_table_C1 = temp_table([temp_table(:).cond] == 1);
    temp_table_C2 = temp_table([temp_table(:).cond] == 2);
    temp_table_C3 = temp_table([temp_table(:).cond] == 3);
    temp_table_C4 = temp_table([temp_table(:).cond] == 4);
    temp_table_C5 = temp_table([temp_table(:).cond] == 5);
    
    subplot(5,3,1+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).theta));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).theta));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).theta));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).theta));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).theta));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
    ylabel('Theta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    
    subplot(5,3,2+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).alpha));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).alpha));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).alpha));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).alpha));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).alpha));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);

    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    subplot(5,3,3+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).beta));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).beta));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).beta));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).beta));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).beta));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    
    i = i + 1;
    
end

figure();set(gcf,'color','white','position',[100 300 1000 1000]);
i = 1;
for k = 13:length(fooof_group_results_org)
    temp_table = fooof_group_results_org{k}( [fooof_group_results_org{k}(:).group] == 1 & strcmp({fooof_group_results_org{k}(:).study},'Imagined'));
    temp_table_C1 = temp_table([temp_table(:).cond] == 1);
    temp_table_C2 = temp_table([temp_table(:).cond] == 2);
    temp_table_C3 = temp_table([temp_table(:).cond] == 3);
    temp_table_C4 = temp_table([temp_table(:).cond] == 4);
    temp_table_C5 = temp_table([temp_table(:).cond] == 5);
    
    subplot(5,3,1+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).theta));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).theta));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).theta));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).theta));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).theta));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
    ylabel('Theta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    
    subplot(5,3,2+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).alpha));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).alpha));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).alpha));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).alpha));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).alpha));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);

    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    subplot(5,3,3+(i-1)*3);hold on;
    a1 = 10*(vertcat(temp_table_C1(:).beta));
    fooof_scatter_plot(ones(size(a1,1),1),a1,color.darkGray);
    a2 = 10*(vertcat(temp_table_C2(:).beta));
    fooof_scatter_plot(ones(size(a2,1),1)*2,a2,color.Green);
    a3 = 10*(vertcat(temp_table_C3(:).beta));
    fooof_scatter_plot(ones(size(a3,1),1)*3,a3,color.Yellow);
    a4 = 10*(vertcat(temp_table_C4(:).beta));
    fooof_scatter_plot(ones(size(a4,1),1)*4,a4,color.Orange);
    a5 = 10*(vertcat(temp_table_C5(:).beta));
    fooof_scatter_plot(ones(size(a5,1),1)*5,a5,color.Red);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',12);
    title(['Cluster ', num2str(k)]);

    
    i = i + 1;
    
end
%% Take average of the identified peaks within participant
% Compare fooofed peaks within comps and participants
% only looking at alpha and beta band
color_jet = jet(15);
k = 12;
temp_table = fooof_group_results_org{k}( [fooof_group_results_org{k}(:).group] == 1 & strcmp({fooof_group_results_org{k}(:).study},'Walking'));
figure();set(gcf,'color','white','position',[200 300 500 800]);
cond_label  = {'rest','flat','low','med','high'};

for cond = 1:5
    temp_table_C1 = temp_table([temp_table(:).cond] == cond);

    subplot(5,2,cond*2-1); hold on;
    rand_vec = [-0.5:0.05:0.5];
    p = 1;        
    
    for i = 1:length(temp_table_C1)
        x1 = 10*vertcat(temp_table_C1(i).alpha);  
        fooof_scatter_plot(ones(size(x1,1),1)*rand_vec(p),x1,color_jet(temp_table_C1(i).subID,:),'connect',1);          
        if isempty(x1)
           scatter(rand_vec(p),0,[],'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        if (i<length(temp_table_C1)) & (temp_table_C1(i).subID ~= temp_table_C1(i+1).subID)
            p = p + 1;
        end
    end
    ylabel('Alpha peak power');
    xlim([-0.6 rand_vec(p)]);
    set(gca,'xtick',linspace(-0.6 ,(rand_vec(p)-0.6)/2,1),'xticklabel',cond_label(cond),'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 20])

    subplot(5,2,cond*2); hold on;
    rand_vec = [-0.5:0.05:0.5];
    p = 1;
    for i = 1:length(temp_table_C1)
        x1 = 10*vertcat(temp_table_C1(i).beta);            
        fooof_scatter_plot(ones(size(x1,1),1)*rand_vec(p),x1,color_jet(temp_table_C1(i).subID,:),'connect',1);         
        if isempty(x1)
           scatter(rand_vec(p),0,[],'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        if (i<length(temp_table_C1)) & (temp_table_C1(i).subID ~= temp_table_C1(i+1).subID)
            p = p + 1;
        end
    end
    ylabel('Beta peak power');
    xlim([-0.6 rand_vec(p)]);
    set(gca,'xtick',linspace(-0.6 ,(rand_vec(p)-0.6)/2,1),'xticklabel',cond_label(cond),'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 10])
end

%% Plot max of the peak within same comp and median within participants    
C1_table = table;
for k = 3:length(fooof_group_results_org)
    temp_table_C1 = fooof_group_results_org{k};
    for i = 1:length(temp_table_C1)
        if isempty(temp_table_C1(i).alpha)
            temp_table_C1(i).alpha = [NaN, NaN];
        end
        if isempty(temp_table_C1(i).beta)
            temp_table_C1(i).beta = [NaN, NaN];
        end
        if isempty(temp_table_C1(i).alpha_center)
            temp_table_C1(i).alpha_center = [NaN, NaN];
        end
        
        if isempty(temp_table_C1(i).beta_center)
            temp_table_C1(i).beta_center = [NaN, NaN];
        end

        [~,idx_a] = max(temp_table_C1(i).alpha(:,2));
        [~,idx_b] = max(temp_table_C1(i).beta(:,2));

        C1_table = vertcat(C1_table,table(temp_table_C1(i).subID,temp_table_C1(i).compID,{temp_table_C1(i).study},...
            temp_table_C1(i).cond,temp_table_C1(i).group,temp_table_C1(i).cluster,temp_table_C1(i).aperiodic_exp,...
            temp_table_C1(i).r_squared,temp_table_C1(i).alpha(idx_a,1),temp_table_C1(i).alpha(idx_a,2),temp_table_C1(i).beta(idx_b,1),temp_table_C1(i).beta(idx_b,2),...
            temp_table_C1(i).alpha_center(1,1),temp_table_C1(i).alpha_center(1,2),temp_table_C1(i).beta_center(1,1),temp_table_C1(i).beta_center(1,2),...
            'VariableNames',...
            {'subID','compID','study','cond','group','cluster','aperiodic_exp','r_squared','alpha_cf','alpha_p',...
            'beta_cf','beta_p','alpha_center','alpha_centerP','beta_center','beta_centerP'}));
    end
end

grp_C1_table = grpstats(C1_table,["subID","study","cond","group","cluster"],["nanmedian"],"DataVars",["r_squared","alpha_cf","alpha_p",...
            "beta_cf","beta_p","alpha_center","alpha_centerP","beta_center","beta_centerP"]);
grp_C1_table.Properties.VariableNames = {'subID','study','cond','group','cluster','GroupCount','med_r_squared','med_alpha_cf','med_alpha_p',...
            'med_beta_cf','med_beta_p','med_alpha_center','med_alpha_centerP','med_beta_center','med_beta_centerP'};
%% Plot med alpha p
figure();set(gcf,'color','white','position',[100 300 1000 1000]);
p = 1;study = 'Imagined';
for k = 3:8%length(fooof_group_results_org)
    subplot(6,2,p);hold on;
    clear med_alpha
    t1 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,study),:);
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond),10*t1.med_alpha_p(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_alpha_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_alpha(cond) = nanmedian(10*t1.med_alpha_p(t1.cond == cond ));
    end
    plot(1:5,med_alpha,'-','color','black','linewidth',2);
    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 20]);

    clear med_beta
    subplot(6,2,p+1);hold on;
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_beta_p(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_beta_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_beta_p(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_beta(cond) = nanmedian(10*t1.med_beta_p(t1.cond == cond ));
    end
    plot(1:5,med_beta,'-','color','black','linewidth',2);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 10]);
    p = p + 2;
end
figure();set(gcf,'color','white','position',[100 300 1000 1000]);
p = 1;
for k = 9:length(fooof_group_results_org)
    subplot(6,2,p);hold on;
    clear med_alpha
    t1 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,study),:);
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_alpha_p(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_alpha_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_alpha(cond) = nanmedian(10*t1.med_alpha_p(t1.cond == cond ));
    end
    plot(1:5,med_alpha,'-','color','black','linewidth',2);
    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 20]);

    clear med_beta
    subplot(6,2,p+1);hold on;
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_beta_p(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_beta_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.3,cond+0.3,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_beta(cond) = nanmedian(10*t1.med_beta_p(t1.cond == cond ));
    end
    plot(1:5,med_beta,'-','color','black','linewidth',2);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k)]);
    ylim([0 10]);
    p = p + 2;
end

%% Plot center alpha p and beta p
figure();set(gcf,'color','white','position',[100 300 1000 1000]);
p = 1;study = 'Walking';
for k = 3:8%length(fooof_group_results_org)
    subplot(6,2,p);hold on;
    clear med_alpha
    t1 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,study),:);
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond),10*t1.med_alpha_centerP(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_alpha_centerP(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_alpha(cond) = nanmedian(10*t1.med_alpha_centerP(t1.cond == cond ));
    end
    plot(1:5,med_alpha,'-','color','black','linewidth',2);
    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k),' Center Power']);
    ylim([0 20]);

    clear med_beta
    subplot(6,2,p+1);hold on;
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_beta_p(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_beta_centerP(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_beta_centerP(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_beta(cond) = nanmedian(10*t1.med_beta_centerP(t1.cond == cond ));
    end
    plot(1:5,med_beta,'-','color','black','linewidth',2);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k),' Center Power']);
    ylim([0 10]);
    p = p + 2;
end
figure();set(gcf,'color','white','position',[100 300 1000 1000]);
p = 1;
for k = 9:length(fooof_group_results_org)
    subplot(6,2,p);hold on;
    clear med_alpha
    t1 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,study),:);
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_alpha_centerP(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_alpha_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.2,cond+0.2,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_alpha(cond) = nanmedian(10*t1.med_alpha_centerP(t1.cond == cond ));
    end
    plot(1:5,med_alpha,'-','color','black','linewidth',2);
    ylabel('Alpha peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k),' Center Power']);
    ylim([0 20]);

    clear med_beta
    subplot(6,2,p+1);hold on;
    for cond = 1:5
        scatter(t1.cond(t1.cond == cond ),10*t1.med_beta_centerP(t1.cond == cond ),'MarkerFaceColor',color.terrain(cond ,:),'MarkerEdgeColor','black');
        if sum(isnan(t1.med_beta_p(t1.cond == cond )))
            numNaN = sum(isnan(t1.med_alpha_p(t1.cond == cond )));
            scatter(linspace(cond-0.3,cond+0.3,numNaN),zeros(numNaN,1),'MarkerFaceColor','white','MarkerEdgeColor','black');
        end
        med_beta(cond) = nanmedian(10*t1.med_beta_centerP(t1.cond == cond ));
    end
    plot(1:5,med_beta,'-','color','black','linewidth',2);
    ylabel('Beta peak power');
    xlim([0 6]);
    set(gca,'xtick',linspace(1,5,5),'xticklabel',{'rest','flat','low','med','high'},'fontsize',10);
    title(['Cluster ', num2str(k), ' Center Power']);
    ylim([0 10]);
    p = p + 2;
end
%% Run Prelim stats
k = 13
t1 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,'Imagined'),:);
t1.cond = categorical(t1.cond);
lme_alpha = fitlme(t1,'med_alpha_p ~ cond + (1|subID)');
anova(lme_alpha)
lme_beta = fitlme(t1,'med_beta_p ~ cond + (1|subID)');
anova(lme_beta)
lme_alpha_center = fitlme(t1,'med_alpha_centerP ~ cond + (1|subID)');
anova(lme_alpha_center )
lme_beta_center = fitlme(t1,'med_beta_centerP ~ cond + (1|subID)');
anova(lme_beta_center )

k = 13
t2 = grp_C1_table(grp_C1_table.cluster == k & grp_C1_table.group == 1 & strcmp(grp_C1_table.study,'Walking')& grp_C1_table.cond ~= 1,:);
t2.cond = categorical(t2.cond);
lme_alpha = fitlme(t2,'med_alpha_p ~ cond + (1|subID)');
anova(lme_alpha)
lme_beta = fitlme(t2,'med_beta_p ~ cond + (1|subID)');
anova(lme_beta)
lme_alpha_center = fitlme(t2,'med_alpha_centerP ~ cond + (1|subID)');
anova(lme_alpha_center )
lme_beta_center = fitlme(t2,'med_beta_centerP ~ cond + (1|subID)');
anova(lme_beta_center )
