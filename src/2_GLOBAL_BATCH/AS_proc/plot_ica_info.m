%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/AS_proc/run_plot_ica_info.sh

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
%% REQUIRED SETUP 4 ALL SCRIPTS
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
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '2_GLOBAL_BATCH'  filesep 'AS_proc'];
%- cd to source directory
cd(source_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP
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
%% ===================================================================== %%
%## Data PATHS
fprintf('\nData PATHS\n');
%## DATASET SPECIFIC
SUBJ_PICS = {{'02','03','04','05','09','11','15','16','18','19','21','22',...
            '23','24','25','27','28','29','30','31','32','33','35','36','38'}};
SUBJ_ITERS = {(1:length(SUBJ_PICS{1}))};
%% (PARAMETERS) ======================================================== %%
fprintf('\nData Processing Parameters\n');
%## Hard Defines
%- dataset specific
DATA_SET = 'AS_dataset';
%- psd params
SPEC_YLIM = [-30,-5];
RECOMPUTE_SPEC = true;
FREQ_LIMITS = [3,50];
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
%- psd stats
condstats = 'on';        % ['on'|'off]
statsMethod = 'perm';    % ['param'|'perm'|'bootstrap']
Alpha = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
mcorrect = 'cluster'; %fdr
groupstats = 'off';
mode = 'fieldtrip';
singletrials = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
%- 
COND_CHARS = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'}; %'1Bounce_BM'
EVENT_CHARS = {'Subject_hit'}; 
test_1 = {'1Bounce_Human','2Bounce_Human'};
test_2 = {'2Bounce_Human','2Bounce_BM'};
test_3 = {'1Bounce_Human','2Bounce_Human','2Bounce_BM'};
COND_DESIGNS = {test_1,test_2,test_3};
%- datetime override
% dt = '05252023_bounces_1h2h2bm_JS';
% dt = '06122023_bounces_1h2h2bm_JS';
dt = '06152023_bounces_1h2h2bm_JS';
%## Soft Define
%- combinations of events and conditions
EVENT_COND_COMBOS = cell(length(COND_CHARS)*length(EVENT_CHARS),1);
cnt = 1;
for cond_i = 1:length(COND_CHARS)
    for event_i = 1:length(EVENT_CHARS)
        EVENT_COND_COMBOS{cnt} = sprintf('%s_%s',COND_CHARS{cond_i},EVENT_CHARS{event_i});
        cnt = cnt + 1;
    end
end
%- path for local data
DATA_DIR = [source_dir filesep '_data'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_1 = sprintf('%s_EPOCH_study',[EVENT_COND_COMBOS{:}]);
cluster_info_fpath = [STUDIES_DIR filesep 'as_cluster_info' filesep 'cluster_info.mat'];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(EVENT_COND_COMBOS)*4]);
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
% STUDY.cluster = STUDY.etc.cluster_save;
% disp(STUDY.cluster);
% disp(STUDY.urcluster);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
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
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',FREQ_LIMITS,'logtrials','on'});
    end
end
%% PLOT DIPOLES, PSD'S, && TOPOPLOTS
fprintf('==== Making Dipole Plots ====\n');
%     [~] = std_dipplot(STUDY,ALLEEG,'clusteras',2:length(STUDY.cluster),...
%             'mode','multicolor','figure','on');
%- main plot
%{
[~] = std_dipplot(STUDY,ALLEEG,'clusters',main_cl_inds(2:end),...
            'mode','multicolor','figure','on');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('allDipPlot_coronal.jpg')]);
%- outlier plot
[~] = std_dipplot(STUDY,ALLEEG,'clusters',outlier_cl_inds,...
            'mode','multicolor','figure','on');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('outlierDipPlot_coronal.jpg')]);
%}
%{
for cluster_i = 2:length(STUDY.cluster)
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
std_topoplot(STUDY,ALLEEG,'clusters',main_cl_inds);
fig_i = get(groot,'CurrentFigure');
fig_i.Position = [500 300 1080 720];
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.fig')]);
saveas(fig_i,[save_dir filesep sprintf('allTopoPlot.jpg')]);
close all
%## POP VIEW PROPS
%{
if ~exist([save_dir filesep 'component_props'],'dir') 
    mkdir([save_dir filesep 'component_props']);
end
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
close all
%}
%%
STUDY = pop_statparams(STUDY, 'condstats', condstats,...
                    'method',statsMethod,...
                    'singletrials',singletrials,'mode',mode,'fieldtripalpha',Alpha,...
                    'fieldtripmethod','montecarlo','fieldtripmcorrect',mcorrect,'fieldtripnaccu',10000);
%% MAKEDESIGN & PLOT
for des_i = 1:length(COND_DESIGNS)
    for cluster_i = main_cl_inds %inds
        %-
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
                'subjselect', {ALLEEG.subject},...
                'variable1','bounces',...
                'values1', COND_DESIGNS{des_i});
        %-
    %     [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
    %             'clusters',cluster_i,'comps','all','subject','','freqrange', FREQ_LIMITS,'subtractsubjectmean','on');
        %-
        [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG,...
                'clusters',cluster_i,'freqrange', FREQ_LIMITS,'plotmode','condensed',...
                'plotconditions','together','ylim',SPEC_YLIM); %'plotgroups','together'
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.fig',cluster_i,[COND_DESIGNS{des_i}{:}])]);
        saveas(fig_i,[save_dir filesep sprintf('%i_%s_specplot.jpg',cluster_i,[COND_DESIGNS{des_i}{:}])]);
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

