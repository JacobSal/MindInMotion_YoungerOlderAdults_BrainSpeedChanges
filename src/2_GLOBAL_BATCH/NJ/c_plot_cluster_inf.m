%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/NJ/run_c_plot_cluster_inf.sh

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
%% DEFINE SOURCE DIRECTORY & CD ======================================== %%
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'NJ'];
%- cd to source directory
cd(source_dir)
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
    SLURM_POOL_SIZE = 1;
end
%% (PARAMETERS) ======================================================== %%
%## Hard Defines
%- spec generation parameters
RECOMPUTE_SPEC = true;
FREQ_FAC = 4;
FREQ_LIMITS = [1,70];
SPEC_MODE = 'psd'; %'fft'; %'psd'; 
%- spec plotting parameters
SPEC_FREQLIMITS = [3,50];
SPEC_NACCU = 2000;
% SPEC_XLIM = [1,100];
SPEC_YLIM = [-30,-8];
SPEC_CONDSTATS = 'on';        % ['on'|'off]
SPEC_STAT_METHOD = 'perm';    % ['param'|'perm'|'bootstrap']
SPEC_ALPHA = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
SPEC_MCORRECT = 'cluster'; %fdr
SPEC_GROUPSTATS = 'off';
SPEC_STAT_MODE = 'fieldtrip'; 
SPEC_STAT_FTMETHOD = 'montecarlo';
SPEC_SINGLETRIALS = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
%- 
TRIAL_TYPES = {'pre','post'};
COND_DESIGNS = {{'pre','post'}};
%- hardcode data_dir
DATA_SET = 'jacobsenN_dataset';
%- datetime override
dt = '06292023_NJ_Standing';
%## Soft Defines
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(TRIAL_TYPES)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_1 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
%-
%{
cluster_load_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
pick_cluster = 12;
clustering_weights.dipoles = 5;
clustering_weights.scalp = 5;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
cluster_alg = 'kmeans';
do_multivariate_data = 1;
evaluate_method = 'min_rv';
clustering_method = ['dipole_',num2str(clustering_weights.dipoles),...
    '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
    '_spec_',num2str(clustering_weights.spec)];
outputdir = [cluster_load_dir filesep 'clustering_solutions' filesep clustering_method filesep num2str(pick_cluster) filesep evaluate_method];
tmp = load(fullfile(outputdir,['cluster_update',num2str(pick_cluster),'.mat']));
STUDY.cluster = tmp.cluster_update;
%}
%-
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
    fprintf('Running Spectrum Calculation...');
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',FREQ_LIMITS,'logtrials','on'});
    end
end
%% PLOT DIPOLES, PSD'S, && TOPOPLOTS
% fprintf('==== Making Dipole Plots ====\n');
%     [~] = std_dipplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster),...
%             'mode','multicolor','figure','on');
%- main plot
%{
%## custom
inds = [2,3,4,7,8,9];

[~] = std_dipplot(STUDY,ALLEEG,'clusters',inds,...
            'mode','multicolor','figure','on');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_coronal.jpg')]);
%##
[~] = std_dipplot(STUDY,ALLEEG,'clusters',main_cl_inds(2:end),...
            'mode','multicolor','projlines','off');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_coronal.jpg')]);
%## outlier plot
[~] = std_dipplot(STUDY,ALLEEG,'clusters',outlier_cl_inds,...
            'mode','multicolor','figure','on');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_coronal.jpg')]);
%##
for cluster_i = inds
    std_dipplot(STUDY,ALLEEG,'clusters',cluster_i);
    fig_i = get(groot,'CurrentFigure');
%     saveas(fig_i,[save_dir filesep sprintf('DipPlot_0.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('cl%i_DipPlot_top.jpg',cluster_i)]);
    view([45,0,0])
%     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_1.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('cl%i_allDipPlot_sagittal.jpg',cluster_i)]);
    view([0,-45,0])
%     saveas(fig_i,[save_dir filesep sprintf('allDipPlot_2.fig')]);
    saveas(fig_i,[save_dir filesep sprintf('cl%i_allDipPlot_coronal.jpg',cluster_i)]);
end
%}
%- Spec plot
fprintf('==== Making Spectogram Plots ====\n');
std_specplot(STUDY,ALLEEG,'clusters',main_cl_inds,...
    'freqrange',[1,100]);
fig_i = get(groot,'CurrentFigure');
fig_i.Position = [500 300 1080 720];
saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allSpecPlot.jpg')]);
%- Topo plot
fprintf('==== Making Topograph Plots ====\n');
std_topoplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster));
fig_i = get(groot,'CurrentFigure');
fig_i.Position = [500 300 1080 720];
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.jpg')]);
close all
%{
%## POP VIEW PROPS
if ~exist([save_dir filesep 'component_props'],'dir') 
    mkdir([save_dir filesep 'component_props']);
    %## RECALCULATE ICAACT MATRICES
    ALLEEG = eeg_checkset(ALLEEG,'loaddata');
    for subj_i = 1:length(ALLEEG)
        if isempty(ALLEEG(subj_i).icaact)
            fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
            ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
            ALLEEG(subj_i).icaact    = reshape( ALLEEG(subj_i).icaact, size(ALLEEG(subj_i).icaact,1), ALLEEG(subj_i).pnts, ALLEEG(subj_i).trials);
        end
    end
    for cluster_i = 2:length(STUDY.cluster)
        sets_clust = STUDY.cluster(cluster_i).sets;
        for i = 1:length(sets_clust)
            subj_i = sets_clust(i);
            comps_clust = STUDY.cluster(cluster_i).comps(i);
            hold on;
            pop_prop_extended(ALLEEG(subj_i),0,comps_clust,NaN,...
            {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
            fig = get(groot,'CurrentFigure');
            fig.Name = sprintf('%s_cluster%i_ic%i',...
                                cluster_i,ALLEEG(subj_i).subject,comps_clust);
            fig.Position = [500 300 1280 720]; 
            hold off;

%                 saveas(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.fig',cluster_i,ALLEEG(subj_i).subject,comps_clust)]);
            saveas(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.png',cluster_i,ALLEEG(subj_i).subject,comps_clust)]);
%                 savefig(fig,[save_dir filesep 'component_props' filesep sprintf('cluster%i_%s_viewprops_ic%i.fig',cluster_i,ALLEEG(subj_i).subject,comps_clust)]);
            close all
        end
    end
end
close all
%}
%% ASSIGN GROUP NAME & STUDY STATS
STUDY = pop_statparams(STUDY, 'condstats', SPEC_CONDSTATS,...
                    'method',SPEC_STAT_METHOD,...
                    'singletrials',SPEC_SINGLETRIALS,'mode',SPEC_STAT_MODE,'fieldtripalpha',SPEC_ALPHA,...
                    'fieldtripmethod',SPEC_STAT_FTMETHOD,'fieldtripmcorrect',SPEC_MCORRECT,'fieldtripnaccu',SPEC_NACCU);
%% MAKEDESIGN & PLOT
for des_i = 1:length(COND_DESIGNS)
    %-
    [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
            'subjselect', {ALLEEG.subject},...
            'variable1','type',...
            'values1', COND_DESIGNS{des_i});
    for cluster_i = main_cl_inds
        %-
    %     [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
    %             'clusters',cluster_i,'comps','all','subject','','freqrange', FREQ_LIMITS,'subtractsubjectmean','on');
        %-
        [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
                'clusters',cluster_i,'freqrange',SPEC_FREQLIMITS,'plotmode','condensed',...
                'plotconditions','together','ylim',SPEC_YLIM); %'plotgroups','together'
        fig_i = get(groot,'CurrentFigure');
        %- set figure line colors
        cc = linspecer(length(COND_DESIGNS{des_i}));
        iter = 1;
        for i = 1:length(fig_i.Children(2).Children)
            set(fig_i.Children(2).Children(i),'LineWidth',1.5);
            set(fig_i.Children(2).Children(i),'Color',horzcat(cc(iter,:),0.6));
            if iter == size(cc,1)
                iter = 1;
            else
                iter = iter + 1;
            end                
        end
        %- tight layout
        % (06/22/2023) JS, seems like changing the Position property also
        % erases the statistics bar. Not sure how to fix?
%         set(fig_i.Children(2),'Position',[0.26,0.26,0.54,0.51]); %DEFAULT
%         set(fig_i.Children(3),'Position',[0.26,0.2345,0.54,0.0255]); %DEFAULT
        set(fig_i.Children(2),'Position',[0.15,0.15,0.8,0.8]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
        set(fig_i.Children(3),'Position',[0.15,0.15-0.0255,0.8,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
        set(fig_i.Children(1),'Location','southeast') %reset Legend
        drawnow;
        %- save
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        close(fig_i)
    end
end
%% Version History
%{
v1.0; (11/11/2022), JS: really need to consider updating bootstrap
    algorithm with parallel computing. Taking ~ 1 day per
    condition for all subjects and the bottle neck is entirely the
    bootstrap.

    Note: validateattributes and assert functions may be helpful
    in more clearly defining function inputs.
        e.g.  DO_PHASE_RND = true;
          errorMsg = 'Value must be (true/false). Determines whether a phase randomized distribution will be created.'; 
          validationFcn = @(x) assert(islogical(x),errorMsg);
v1.0; (12/5/2022) Need to adapt this to include all conditions
    within each SUBJ structure so connectivity can be calculated
    for the ALLEEG structure rather than the EEG structure.
    *** maybe try to ditch the SUBJ strucutre entirely for this
    round?
v1.0.01132023.0 : Initializing versioning for future iterations.
%}

