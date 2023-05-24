%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/gamma/MIM_OA_proc/run_plot_ica_info.sh

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
run_dir = [source_dir filesep '2_GLOBAL_BATCH' filesep 'MIM_OA_proc'];
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
%% (DATASET INFORMATION) =============================================== %%
%## (MIND IN MOTION) DATASET SPECIFIC (05/24/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H1008','H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};
SUBJ_NO_MRI = {'H2010', 'H2036', 'H2041', 'H2072', 'H3018','H3120'};
SUBJ_1YA = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
    'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1032','H1033','H1034','H1035',...
    'H1036','H1037','H1038','H1039','H1041','H1042','H1044','H1045','H1047','H1047'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2017', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2037', 'H2038',...
    'H2039', 'H2042', 'H2052', 'H2059', 'H2062', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- (OY) Subject Picks 
% SUBJ_PICS = {SUBJ_1YA}; 
% GROUP_NAMES = {'H1000''s'}; 
% SUBJ_ITERS = {1:length(SUBJ_1YA)}; 
%- (OA) Subject Picks 
SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
%% (PARAMETERS) ======================================================== %%
%## Hard Defines
%- datset name
DATA_SET = 'MIM_dataset';
%- datetime override
% dt = '04172023_MIM_OA_subset_N85_speed_terrain_merge';
% dt = '05012023_MIM_OA_subset_N85_speed_terrain_merge';
dt = '04172023_MIM_OA_subset_N85_speed_terrain_merge';
%- epoching params
% TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%## Soft Defines
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
study_fName_3 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(load_trials)*4]);
else
    POOL_SIZE = 1;
end
%- Create STUDY & ALLEEG structs
if ~exist([load_dir filesep study_fName_3 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_3 '_UNIX.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_3 '.study'],'filepath',load_dir);
    end
    [comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
end
% STUDY.cluster = STUDY.etc.cluster_save;
% disp(STUDY.cluster);
% disp(STUDY.urcluster);
%% PLOT DIPOLES, PSD'S, && TOPOPLOTS
% fprintf('==== Making Dipole Plots ====\n');
%     [~] = std_dipplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster),...
%             'mode','multicolor','figure','on');
%- main plot
%{
%## custom
inds = [2,3,4,7,8,9];
[~] = std_dipplot(STUDY,ALLEEG,'clusters',inds,...
            'mode','multicolor','figure','on',');
fig_i = get(groot,'CurrentFigure');
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_top.jpg')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_sagittal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('partDipPlot_coronal.jpg')]);]
%##
[~] = std_dipplot(STUDY,ALLEEG,'clusters',main_cl_inds(2:end),...
            'mode','multicolor','figure','on');
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
%%
%{
parfor subj_i = 1:length(ALLEEG)
    reject_struct = mim_reject_ics(ALLEEG(subj_i),ALLEEG(subj_i).filepath);
    tmp_bad = setdiff((1:size(ALLEEG(subj_i).icaweights,1)),find((reject_struct.IC_all_brain >= THRESH_BRAIN_SCORE & reject_struct.IC_all_brain ~= 9)));
    tmp_good = find(reject_struct.IC_all_brain >= THRESH_BRAIN_SCORE & reject_struct.IC_all_brain ~= 9);
    ALLEEG(subj_i).etc.urreject = [];
    ALLEEG(subj_i).etc.urreject.crit = [];
    ALLEEG(subj_i).etc.urreject.ic_keep = [];
    ALLEEG(subj_i).etc.urreject.ic_rej = [];
    ALLEEG(subj_i).etc.urreject.dipfit = [];
    if isempty(tmp_good)
        continue;
    end
    ALLEEG(subj_i).etc.urreject.crit = reject_struct;
    ALLEEG(subj_i).etc.urreject.ic_keep = tmp_good;
    ALLEEG(subj_i).etc.urreject.ic_rej = tmp_bad;
    ALLEEG(subj_i).etc.urreject.dipfit = ALLEEG(subj_i).dipfit;
    fprintf('** Subject %s has %i brain components\n',ALLEEG(subj_i).subject, length(tmp_good));
end
%}
%%
RECOMPUTE_SPEC = false;
%-
condstats = 'on';        % ['on'|'off]
statsMethod = 'perm';    % ['param'|'perm'|'bootstrap']
Alpha = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
mcorrect = 'cluster'; %fdr
groupstats = 'off';
mode = 'fieldtrip';
singletrials = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
%- 
speed_trials = {'0p25','0p5','0p75','1p0'};
terrain_trials = {'flat','low','med','high'};
COND_DESIGNS = {speed_trials,terrain_trials};
%-
inds = [2,3,4,7,8,9];
%## std_precomp.m params
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
FREQ_LIMITS = [3,40];
SPEC_YLIM = [-30,-5];
%% ASSIGN GROUP NAME & STUDY STATS
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = 'Older Adults';
    STUDY.datasetinfo(subj_i).group = 'Older Adults';
end
STUDY = pop_statparams(STUDY, 'condstats', condstats,...
                    'method',statsMethod,...
                    'singletrials',singletrials,'mode',mode,'fieldtripalpha',Alpha,...
                    'fieldtripmethod','montecarlo','fieldtripmcorrect',mcorrect,'fieldtripnaccu',10000);
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
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',FREQ_LIMITS,'logtrials','on'});
    end
end
%% PRECOMPUTE
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
%                                     'components',... 
%                                     'allcomps','on',...
%                                     'recompute','on',...
%                                     'spec','on',...
%                                     'specparams',...
%                                     {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
%                                     'freqrange',FREQ_LIMITS,'logtrials','on'});
%% MAKEDESIGN & PLOT
for des_i = 1:length(COND_DESIGNS)
    for cluster_i = inds % main_cl_inds
        %-
        [STUDY] = std_makedesign(STUDY, ALLEEG, des_i,...
                'subjselect', {ALLEEG.subject},...
                'variable1','cond',...
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

