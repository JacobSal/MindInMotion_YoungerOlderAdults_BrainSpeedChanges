f%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen, Chang Liu, Ryan Downey
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/_test/3_paper_MIM_HOA/run_alleeg_spectral_run.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
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
    DO_UNIX = false;
    PATH_EXT = 'M';
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
end
%## DEBUG: PATHROOT OVERRIDE
if DO_UNIX
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '_test' filesep '3_paper_MIM_HOA'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace
%% PARPOOL SETUP
if ~ispc
%     eeg_options;
    % see. eeg_optionsbackup.m for all eeglab options.
%     pop_editoptions('option_parallel',0,'option_storedisk',1,...
%         'option_saveversion6',0,'option_cachesize',1000,'option_savetwofiles',1);
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
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
%     pop_editoptions('option_parallel',0,'option_storedisk',1,...
%         'option_saveversion6',0,'option_cachesize',1000,'option_savetwofiles',1);
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% ================================================================= %%
%## PATHS
%- hardcode data_dir
DATA_SET = 'MIM_dataset';
DATA_DIR = [source_dir filesep '_data'];
%- path for local data
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
%## DATASET SPECIFIC
% SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
%             'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
%             'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
% SUBJ_MISSING_TRIAL_DATA = {'H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
%     'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};

SUBJ_2HMA = {'H2017', 'H2010', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2036', 'H2037', 'H2038',...
    'H2039', 'H2041', 'H2042', 'H2052', 'H2059', 'H2062', 'H2072', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
% SUBJ_3HMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
%     'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120'}; % JACOB,SAL(02/23/2023)
SUBJ_3NHMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120',...
    'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- Subject Picks
% SUBJ_PICS = {SUBJ_2HMA,SUBJ_3HMA,SUBJ_3NHMA};
SUBJ_PICS = {SUBJ_2HMA,SUBJ_3NHMA};
GROUP_NAMES = {'H2000''s','H3000''s'};
% SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3HMA),1:length(SUBJ_3NHMA)}; % JACOB,SAL(02/23/2023)
SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3NHMA)};
%% ===================================================================== %%
%## PROCESSING PARAMS
%- statistics
STAT_ALPHA = 0.05;
%- Spectral params
DO_TIMEWARP = true;
DO_BASELINE_CORRECTION = false;
RECOMPUTE_SPEC = false;
RECOMPUTE_ERSP = false;
FREQ_LIMITS = (1:100);
CYCLE_LIMITS = [3,0.8];
SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
FREQ_FAC = 4;
% PAD_RATIO = 2;
% SPEC_EPOCH_LIMS = [-1,3.4]; %[-2,2];
%- datetime override
% dt = '26122022';
% dt = '16032023_OA_subset';
% dt = '23032023_OA_subset';
% dt = '31032023_OA_subset_randi';
% dt = '04052023_OA_subset';
% dt = '04092023_MIM_OA_subset_N85_speed_terrain';
dt = '04172023_MIM_OA_subset_N85_speed_terrain_merge';
%- hard define
% load_trials = {'0p25'};
load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% load_trials = {'0p25','0p5','0p75','1p0'};
% load_trials = {'flat','low','med','high'}; 
study_fName_1 = sprintf('%s_all_comps_study',[load_trials{:}]);
study_fName_2 = sprintf('%s_reduced_comps_study',[load_trials{:}]);
study_fName_3 = sprintf('%s_MIM_study',[load_trials{:}]);
%- soft define
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
end
%%
if ~exist([load_dir filesep study_fName_2 '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [tmp,~] = pop_loadstudy('filename',[study_fName_2 '_UNIX.study'],'filepath',load_dir);
    else
        [tmp,~] = pop_loadstudy('filename',[study_fName_2 '.study'],'filepath',load_dir);
    end
end
% clust_1 = STUDY.cluster;
STUDY.cluster = tmp.cluster;
disp(STUDY.cluster);
disp(STUDY.urcluster);
% clear tmp
%% RECOMPUTE CLUSTERS?
%{
[tmpS]         = pop_clust(tmpS,tmpA,...
                'algorithm','kmeans',...
                'clus_num',CLUSTER_NUM);
%}
%% REORGANIZE DIPFITS
tmp_subjs = zeros(length(unique(STUDY.cluster(1).sets)),2);
tmp_subjs(:,1) = unique(STUDY.cluster(1).sets);
tmp_rmv_subjs = STUDY.etc.rmvd_subj.inds;
% set_inds = unique(tmpS.cluster(1).sets);
fprintf('==== removing subjects from cluster indices ====\n');
%- create an unscrambling array for removing subjects
iter = 1;
for subj_i = 1:length(tmp_subjs)
    if any(subj_i == tmp_rmv_subjs)
        continue;
    else
        tmp_subjs(subj_i,2) = iter;
        iter = iter + 1;
    end
end
for subj_i = 1:length(tmp_rmv_subjs)
    for cluster_i = 2:length(STUDY.cluster)
        inds = (STUDY.cluster(cluster_i).sets == tmp_rmv_subjs(subj_i));
        if any(inds)
            STUDY.cluster(cluster_i).sets(inds) = [];
            STUDY.cluster(cluster_i).comps(inds) = [];
        else
            continue;
        end
    end
end
for cluster_i = 2:length(STUDY.cluster)
    for comp_i = 1:length(STUDY.cluster(cluster_i).sets)
        inds = (tmp_subjs(:,1) == STUDY.cluster(cluster_i).sets(comp_i));
        if any(inds)
            STUDY.cluster(cluster_i).sets(comp_i) = tmp_subjs(inds,2);
        else
            continue;
        end
    end
end
%- parentcluster alterations
all_sets = [];
all_comps = [];
for clust_i = 2:length(STUDY.cluster)
    all_sets = [all_sets, STUDY.cluster(clust_i).sets];
    all_comps = [all_comps, STUDY.cluster(clust_i).comps];
end
STUDY.cluster(1).comps = all_comps;
STUDY.cluster(1).sets = all_sets;
%% PRECOMPUTE MEASURES
%## COMPUTE SPECTRUMS
if RECOMPUTE_SPEC
    [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG,...
                        'components',...
                        'recompute','on',...
                        'spec','on',...
                        'scalp','on',...
                        'specparams',...
                        {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                        'freqrange',FREQ_LIMITS});
end
%## COMPUTE ERSPs
if RECOMPUTE_ERSP
    allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    for subj_i = 1:length(ALLEEG)
        allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
    end
    % grandAvgWarpTo = median(allWarpTo);
    grandAvgWarpTo = mean(allWarpTo);
    disp(['Grand average (across all subj) warp to: ',num2str(grandAvgWarpTo)]);
    tmp = STUDY;
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i) = eeg_checkset(ALLEEG(subj_i),'loaddata');
        if isempty(ALLEEG(subj_i).icaact)
            fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
            ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        %- determine timewarping parameters
         if DO_TIMEWARP
            timewarp_param = ALLEEG(subj_i).timewarp.latencies;
            timewarpms_param = grandAvgWarpTo;
        else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,ALLEEG(subj),'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',80,'ntimesout',200,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,ALLEEG(subj_i),'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',80,'ntimesout',200,...
                    'baseline',nan,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    end
end
%%  Load cluster
for pick_cluster = [12]
    load(fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(pick_cluster),evaluate_method,['cluster_update',num2str(pick_cluster),'.mat']))
    STUDY.cluster = [];
    STUDY.cluster = cluster_update;
    % SUPER important, need to clean the cache; otherwise the study
    % will just load from previous cache and the results are wrong
    STUDY.cache = [];
    [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'subjselect', all_subjStr,'variable1','cond','values1', {'flat','low','med','high'});
    [STUDY] = std_makedesign(STUDY, ALLEEG, 2, 'subjselect', all_subjStr,'variable1','cond','values1', {'0p25','0p5','0p75','1p0'});

    %% Identify number clusters
    disp(['There are a total of ',num2str(length(STUDY.cluster)-2),' clusters']);
    for currentdesign = 1:2
        %- setup keywords for 
        if currentdesign == 1
            file_keyword = '_Terrain';
            study_walking = 2;
        elseif currentdesign == 2
            file_keyword = '_Speed';
            study_walking = 2;
        end
        %% SAVE all spec data
        if SAVE_SPEC
            freqrange = [1 100];
            STUDY.currentdesign = currentdesign;
            for i = 2:length(STUDY.cluster)
                clear specdata specfreqs specdata2 specfreqs2
                outputdir = fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(pick_cluster),evaluate_method,'ERSP_Plots',['Cluster_' num2str(i)]);
                if ~exist(outputdir);mkdir(outputdir);end

                disp('Computing specdata...')
                tic
                [STUDY, specdata, specfreqs, pgroup, pcond, pinter] = std_specplot(STUDY, ALLEEG, 'clusters',i,'comps','all','subject','','freqrange', freqrange);
                save(fullfile(outputdir,['readSPEC_', num2str(STUDY.currentdesign),file_keyword,'.mat']), 'specdata', 'specfreqs', 'pgroup', 'pcond', 'pinter');
                saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword,'_All_Comps_SPEC.fig']))
                toc 

                disp('Computing specdata with substract subject mean...')
                tic
                [STUDY2, specdata2, specfreqs2, pgroup2, pcond2, pinter2] = std_specplot(STUDY, ALLEEG, 'clusters',i,'comps','all','subject','','freqrange', freqrange,'subtractsubjectmean','on');
                save(fullfile(outputdir,['readSPEC_', num2str(STUDY.currentdesign),file_keyword,'_subSubjectMean.mat']), 'specdata2', 'specfreqs2', 'pgroup2', 'pcond2', 'pinter2');
                saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword,'_All_Comps_SPEC_subSubjectMean.fig']))
                toc 

                fprintf('All saved for %s \n', num2str(i))
                close

            end
        end
        %% SAVE all ERSP data for future plot
        if SAVE_ERSP
        % std_readesrp and std_readspec are obsolete
        % - Setup parameters
            condstats ='on';        % ['on'|'off]
            statsMethod ='perm';    % ['param'|'perm'|'bootstrap']
            Alpha = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
            mcorrect = 'cluster'; %fdr
            groupstats = 'off';
            mode = 'fieldtrip';
            singletrials = 'off' ;  %['on'|'off'] load single trials spectral data (if available). Default is 'off'.
            tmpSTUDY = pop_statparams(STUDY, 'condstats', condstats,...
                    'method',statsMethod,...
                    'singletrials',singletrials,'mode',mode,'fieldtripalpha',Alpha,...
                    'fieldtripmethod','montecarlo','fieldtripmcorrect',mcorrect,'fieldtripnaccu',10000);
            tmpSTUDY_commonbase = pop_erspparams(tmpSTUDY, 'subbaseline','on','timerange',[warpingvalues(1) warpingvalues(5)], 'ersplim',[-1.5 1.5]);  % 'subbaseline' - ['on'|'off'] subtract the same baseline across conditions for ERSP     
            STUDY.currentdesign = currentdesign;
            tmpSTUDY.currentdesign = currentdesign;
            tmpSTUDY_commonbase.currentdesign = currentdesign;
            
            for i = 3:length(STUDY.cluster)
                outputdir = fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(pick_cluster),evaluate_method,'ERSP_Plots',['Cluster_' num2str(i)]);
                if ~exist(outputdir);mkdir(outputdir);end
            end
            parpool(5)
            parfor i = 3:length(STUDY.cluster)
%                 clear allerspdata alltimes allfreqs
                outputdir = fullfile(save_study_folder,'clustering_solutions',clustering_method,num2str(pick_cluster),evaluate_method,'ERSP_Plots',['Cluster_' num2str(i)]);
%                 if ~exist(outputdir);mkdir(outputdir);end
                
                plotERSP_parfor(STUDY,tmpSTUDY,tmpSTUDY_commonbase, ALLEEG,myErspParams,myErspParams_trialbasefull,i,outputdir,file_keyword)
                
                disp(['SAVE ALL ERSP DATA successfully for cluster ', num2str(i)])
            end
        end
    end
end
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}