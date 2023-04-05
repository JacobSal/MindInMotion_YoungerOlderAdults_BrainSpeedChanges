%%% =================================================================== %%%
% IDENTIFY BRAIN-LIKE COMPONENTS
%%% =================================================================== %%%
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   For headmodel paper only

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
%- TIC
tic
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
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
run_dir = [source_dir filesep '3_ANALYZE' filesep '1_ic_validation'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = true;
setWorkspace
%% PARPOOL SETUP
if ~ispc
%     eeg_options;
    pop_editoptions('option_parallel',0,'option_storedisk',1);
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
    %- make meta data directory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions('option_parallel',0,'option_storedisk',1);
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
%- MIND IN MOTION (SUBSET (03/10/2023)
SUBJ_NORUN = {'H2012_FU', 'H2013_FU', 'H2018_FU', 'H2020_FU', 'H2021_FU',...
            'H3024_Case','H3029_FU','H3039_FU','H3063_FU','NH3021_Case', 'NH3023_Case','NH3025_Case', 'NH3030_FU',...
            'NH3068_FU', 'NH3036_FU', 'NH3058_FU'};
SUBJ_MISSING_TRIAL_DATA = {'H2012','H2018','H3024','NH3002', 'NH3004','NH3009',...
    'NH3023','NH3027', 'NH3028', 'NH3129', 'NH3040'};

SUBJ_2HMA = {'H2017', 'H2010', 'H2002', 'H2007', 'H2008', 'H2013', 'H2015',...
    'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2027', 'H2033', 'H2034', 'H2036', 'H2037', 'H2038',...
    'H2039', 'H2041', 'H2042', 'H2052', 'H2059', 'H2062', 'H2072', 'H2082',...
    'H2090', 'H2095', 'H2111', 'H2117'};
SUBJ_3HMA = {'H3018','H3029','H3034','H3039','H3042','H3046',...
    'H3047','H3053','H3063','H3072','H3073','H3077','H3092','H3103','H3107','H3120'}; % JACOB,SAL(02/23/2023)
SUBJ_3NHMA = {'NH3006', 'NH3007', 'NH3008', 'NH3010',...
    'NH3021', 'NH3025', 'NH3026',...
    'NH3030', 'NH3036',...
    'NH3041', 'NH3043', 'NH3051', 'NH3054', 'NH3055', 'NH3056', 'NH3058',...
    'NH3059', 'NH3066', 'NH3068', 'NH3069', 'NH3070', 'NH3071', 'NH3074',...
    'NH3076', 'NH3082', 'NH3086', 'NH3090', 'NH3102', 'NH3104', 'NH3105', 'NH3106',...
    'NH3108', 'NH3110', 'NH3112', 'NH3113', 'NH3114', 'NH3123', 'NH3128'}; % JACOB,SAL(02/23/2023)
%- trial types
% TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
%- Subject Picks
SUBJ_PICS = {SUBJ_2HMA,SUBJ_3HMA,SUBJ_3NHMA};
GROUP_NAMES = {'H2000''s','H3000''s','NH3000''s'};
SUBJ_ITERS = {1:length(SUBJ_2HMA),1:length(SUBJ_3HMA),1:length(SUBJ_3NHMA)}; % JACOB,SAL(02/23/2023)
%- Subject Picks
%- Subject Directory Information
OA_PREP_FPATH = '24022023_OA_prep'; % JACOB,SAL(02/23/2023)
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
%## CONVERT PATHS
if DO_UNIX
    OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR);
else
    OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR);
end
%% ===================================================================== %%
%## PROCESSING PARAMS
%- statistics
STAT_ALPHA = 0.05;
%- datetime override
dt = '16032023_OA_subset';
% dt = '23032023_OA_subset';
%- hard define
% load_trials = {'0p25'};
% load_trials = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
load_trials = {'0p25','0p5','0p75','1p0'};
% load_trials = {'flat','low','med','high'}; 
%- soft define
subjinfDir = [SUBJINF_DIR filesep sprintf('%s',dt)];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
if exist('SLURM_POOL_SIZE','var')
    POOL_SIZE = min([SLURM_POOL_SIZE,length(load_trials)*4]);
else
    POOL_SIZE = 1;
end
STUDIES = cell(1,length(load_trials));
ALLEEGS = cell(1,length(load_trials));
%- Create STUDY & ALLEEG structs
parfor (cond_i = 1:length(load_trials),POOL_SIZE)
    study_fName = sprintf('%s_MIM_study',load_trials{cond_i});
    if ~exist([load_dir filesep study_fName '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if ~ispc
            [tmpS,tmpA] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',load_dir);
        else
            [tmpS,tmpA] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);
        end
    end
end
%%
% clear all;close all;
% eeglab
% addpath 'M:\liu.chang1\scripts\MiM_HY'
% MiM_HY_config_params;
folder_name = 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10';
process_output_folder_Dipfit = ['M:\liu.chang1\STUDY-preprocess-HY_202212\',folder_name];

output_folder = [save_dir filesep 'ic_rejection_crit']; 
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end
%not include data with no custom electrode location
%%
for n = [1:length(all_subjStr)]
    runAllCriteria = 1;
    subjStr = all_subjStr{n}
    try
        filename_set = [subjStr,'_cleanEEG_',folder_name,'_ICA_ICLabel_dipfit_fem.set'];
        EEG = pop_loadset('filename',filename_set,'filepath',fullfile(process_output_folder_Dipfit,subjStr));
    catch 
        filename_set = [subjStr,'_cleanEEG_',folder_name,'ICA_ICLabel_dipfit_fem.set'];
        EEG = pop_loadset('filename',filename_set,'filepath',fullfile(process_output_folder_Dipfit,subjStr));
    end
    mkdir(fullfile(output_folder,subjStr,'Figures'))
    
    if runAllCriteria
        %% Crteria 1: Count brain components based on ICLabel
        % RUN IClabel Lite    
        iclabel_classifier = 'lite';
        EEG = iclabel(EEG,iclabel_classifier);
        EEG.etc.ic_classification.ICLabel.classifications_default = EEG.etc.ic_classification.ICLabel.classifications;
        %'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
        iclabel_classes = 1; 
        iclabel_threshold = 0.5;

        % RUN IClabel Lite
        classifications = EEG.etc.ic_classification.ICLabel.classifications;
        [classes ~] = find(transpose(EEG.etc.ic_classification.ICLabel.classifications == max(EEG.etc.ic_classification.ICLabel.classifications,[],2)));
        classes_to_keep = 1; % Brain components
        threshold_to_keep = 0.5;
        summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
        ICs_keep_brain = find(summed_scores_to_keep>threshold_to_keep);
        ICs_throw = find(summed_scores_to_keep<threshold_to_keep);
        ICs_keepsummed_scores_to_keep(ICs_throw) = 0;

        %-- 50% ICLabel
        ICs_keep_brain50 = find(summed_scores_to_keep>threshold_to_keep);
        threshold_to_keep = 0.5;
        summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
        ICs_keep_brain50 = find(summed_scores_to_keep>threshold_to_keep);

        %-- 75% ICLabel
        classes_to_keep = 1; % Brain components
        threshold_to_keep = 0.75;
        summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
        ICs_keep_brain75 = find(summed_scores_to_keep>threshold_to_keep);

        %-- 90% ICLabel
        classes_to_keep = 1; % Brain components
        threshold_to_keep = 0.9;
        summed_scores_to_keep = sum(classifications(:,classes_to_keep),2);
        ICs_keep_brain90 = find(summed_scores_to_keep>threshold_to_keep);

        % Store other classifications
        threshold_to_remove = 0.5;
        summed_scores_to_keep = sum(classifications(:,2),2);%Muscle
        ICs_muscle= find(summed_scores_to_keep>threshold_to_remove);

        summed_scores_to_keep = sum(classifications(:,3),2);%eye
        ICs_eye= find(summed_scores_to_keep>threshold_to_remove);

        %% Criteria 2 Spectrum graph: Plot PSD together for selected channels
        IC_potential = 1:size(EEG.icawinv,2);
        icaacttmp = EEG.icaact(IC_potential, :, :);
        spec_opt = {'freqrange',[0 80]};

        % compute the spectrum plots: TO DO: need to run it faster
        [spectra_psd,FREQ ]= spectopo( icaacttmp, EEG.pnts, EEG.srate,'percent',80,...
               'freqrange',[2 100], 'plot','off');
        % output unit from spectopo is in db

        % add a linear fit: note that the frequency used here is log(frequency)
        % to help characterize the 1/f structure. y axis should be log(power)
        fit_range = 2:40; % Original setting 2:100. now changed to 2:40 (2022-1-2)
        for i = 1:length(IC_potential)
            lsfit(i,:) = polyfit(log10(FREQ(fit_range,1)), spectra_psd(i,fit_range)', 1);
            spectra_psd_fit(i,:) = lsfit(i,1).*log10(FREQ)+lsfit(i,2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOTS
        % fit linear from 0-150Hz
        % fig_PSD = figure('color','w');
        % for i = 1:10
        %     plot(FREQ,spectra_psd(i,:),'color',c(floor((length(c)-20)/length(1:10))*i,:),'linewidth',1.2);hold on;
        % end
        % title(['IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        % ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
        % xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        % xlim([.5 80]);
        % legend(cellstr(num2str(1:10)),'Location','eastoutside');
        % legend box off

        % IC rejection criteria: only keep those with log fit < 0

        slope_thres = -0.2;  % Added on 1-3-2022 to avoid very flat spectrum that is not brain component

        ICs_keep_brain_spec = IC_potential(lsfit(:,1)<slope_thres)';% actually not a lot get discarded
        ICs_spec_dump = IC_potential(lsfit(:,1)>=slope_thres)';

        %% Inspect 
        numICs = 1:size(EEG.icawinv,2);
        figure();
        stem(numICs,lsfit(:,1)); hold on;
        stem(numICs(lsfit(:,1)> slope_thres),lsfit(lsfit(:,1)> slope_thres,1),'r');
        ylabel('slope')
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['spectral stem.jpg']))

        c = (othercolor('RdBu4'));%my personal color scheme
        % Log-log plot - retained
        fig_PSD2 = figure('color','w','position',[100 200 1200 500]);
        subplot(1,2,1)
        for i0 = 1:length(ICs_keep_brain_spec)
            % plot the fitted line
            p = ICs_keep_brain_spec(i0);
            semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',c(floor((length(c)-20)/length(ICs_keep_brain_spec))*i0,:),'linewidth',1.2);hold on;
            semilogx(FREQ,spectra_psd(p,:),'color',c(floor((length(c)-20)/length(ICs_keep_brain_spec))*i0,:),'linewidth',1.2);hold on;
            grid on;
        end
        title(['KEEP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
        xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        xlim([1 70]);
        legend(cellstr(num2str(ICs_keep_brain_spec)),'Location','eastoutside');
        legend box off
        subplot(1,2,2)
        for i0 = 1:length(ICs_keep_brain_spec)
            % plot the fitted line
            p = ICs_keep_brain_spec(i0);
            plot(FREQ,spectra_psd(p,:),'color',c(floor((length(c)-20)/length(ICs_keep_brain_spec))*i0,:),'linewidth',1.2);hold on;
            grid on;
        end
        title(['KEEP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
        xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        xlim([1 70]);
        legend(cellstr(num2str(ICs_keep_brain_spec)),'Location','eastoutside');
        legend box off
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['Keep potential brain components_spectral.jpg']))

        fig_PSD3 = figure('color','w');
        for i0 = 1:length(ICs_spec_dump)
            % plot the fitted line
            p = ICs_spec_dump(i0);
            semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',c(floor((length(c)-20)/length(ICs_spec_dump))*i0,:),'linewidth',1.2);hold on;
            semilogx(FREQ,spectra_psd(p,:),'color',c(floor((length(c)-20)/length(ICs_spec_dump))*i0,:),'linewidth',1.2);hold on;
            grid on;
        end
        title(['DUMP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
        xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        xlim([1 70]);
        legend(cellstr(num2str(ICs_spec_dump)),'Location','eastoutside');
        legend box off
    %     saveas(gcf,fullfile(save_IC_Rejection_folder,subjStr,'Figures',['potential brain components_spectral_',subDirNum,'.fig']))
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['Dump potential brain components_spectral.jpg']))
        %}

        %% Criteria 4: Scalp topographs and dipole location (if outside the brain or not)
        % Residual variance < 0.15
        IC_RV = vertcat(EEG.dipfit.model.rv);
        RV_threshold_brain = 15;
        ICs_RVthreshold_keep = find(IC_RV <= RV_threshold_brain/100);

        % How to find if the dipole location is inside brain or not

        %% Gather all IC rejection criteria
        All_IC_criteria.IClabel.brain50 = ICs_keep_brain50;
        All_IC_criteria.IClabel.brain75 = ICs_keep_brain75;
        All_IC_criteria.IClabel.brain90 = ICs_keep_brain90;

        All_IC_criteria.IClabel.brain = ICs_keep_brain50;
        All_IC_criteria.IClabel.muscle = ICs_muscle;
        All_IC_criteria.IClabel.muscle_threshold = threshold_to_remove;
        All_IC_criteria.IClabel.eye = ICs_eye;
        All_IC_criteria.IClabel.eye_threshold = threshold_to_remove;
        All_IC_criteria.IClabel.classification = EEG.etc.ic_classification.ICLabel.classifications;

        All_IC_criteria.Spectra.keep = ICs_keep_brain_spec;
        All_IC_criteria.Spectra.dump = ICs_spec_dump;

        All_IC_criteria.RV.keep = ICs_RVthreshold_keep;
        All_IC_criteria.RV.IC_RV = IC_RV;
        EEG.etc.IC_rej = All_IC_criteria;

        % Assign weights to each criteria 
        % If muscle % higher score remove
        IC_all_muscle = zeros(size(EEG.icawinv,2),1);
        IC_all_muscle(All_IC_criteria.IClabel.muscle) = IC_all_muscle(All_IC_criteria.IClabel.muscle)+2;
        IC_all_muscle(All_IC_criteria.Spectra.dump) = IC_all_muscle(All_IC_criteria.Spectra.dump)+1;
    %     IC_all_muscle(All_IC_criteria.Projection.EMG) = IC_all_muscle(All_IC_criteria.Projection.EMG)+2;

        % If eye. % higher score remove
        IC_all_eye = zeros(size(EEG.icawinv,2),1);
        IC_all_eye(All_IC_criteria.IClabel.eye) = IC_all_eye(All_IC_criteria.IClabel.eye)+2;
        IC_all_eye(All_IC_criteria.Spectra.dump) = IC_all_eye(All_IC_criteria.Spectra.dump)+1;

        % If brain % higher score keep
        IC_all_brain = zeros(size(EEG.icawinv,2),1);
        IC_all_brain(All_IC_criteria.IClabel.brain50) = IC_all_brain(All_IC_criteria.IClabel.brain50)+2;
        IC_all_brain(All_IC_criteria.IClabel.brain75) = IC_all_brain(All_IC_criteria.IClabel.brain75)+2;    
        IC_all_brain(All_IC_criteria.Spectra.keep) = IC_all_brain(All_IC_criteria.Spectra.keep)+1;
    %     IC_all_brain(All_IC_criteria.Projection.EMG) = IC_all_brain(All_IC_criteria.Projection.EMG)-3;
        IC_all_brain(All_IC_criteria.RV.keep) = IC_all_brain(All_IC_criteria.RV.keep)+5;
        IC_all_brain(size(EEG.icawinv,2)-5:size(EEG.icawinv,2)) = IC_all_brain(size(EEG.icawinv,2)-5:size(EEG.icawinv,2))-1;

        Output_ICRejection.All_IC_criteria = All_IC_criteria;
        Output_ICRejection.IC_all_brain = IC_all_brain;
        Output_ICRejection.IC_all_muscle = IC_all_muscle;
        Output_ICRejection.IC_all_eye = IC_all_eye;
        Output_ICRejection.Cleaning_Params = EEG.etc.Params;

        fileName = fullfile(output_folder,subjStr,[subjStr,'_ICRej.mat']);
        save(fileName,'Output_ICRejection');

        % Figure
        figure('color','white');
        subplot(3,1,1);
        stem(numICs,IC_all_muscle);title(['muscle: ',num2str(sum(IC_all_muscle == 3))]);ylabel('Dump: score');
        subplot(3,1,2);
        stem(numICs,IC_all_eye);title(['eye: ',num2str(sum(IC_all_eye == 2))]);ylabel('Dump: score');
        subplot(3,1,3);hold on;
        stem(numICs,IC_all_brain);title(['brain: ',num2str(sum(IC_all_brain == 8))]);ylabel('Keep: score');
        stem(numICs(IC_all_muscle>=2),IC_all_brain(IC_all_muscle>=2),'r');
        xlabel('IC number');legend('','flag muscle');
        mkdir(fullfile(output_folder,subjStr,'Figures'));
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['IC_score.jpg']))

        c = (othercolor('RdBu4'));%my personal color scheme
        % Log-log plot - retained
        IC_powpow = find(IC_all_brain >= 8);
        fig = figure('color','w','position',[100 200 1200 500]);
        subplot(1,2,1)
        for i0 = 1:length(IC_powpow)
            % plot the fitted line
            p = IC_powpow(i0);
            semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',c(floor((length(c)-20)/length(IC_powpow))*i0,:),'linewidth',1.2);hold on;
            semilogx(FREQ,spectra_psd(p,:),'color',c(floor((length(c)-20)/length(IC_powpow))*i0,:),'linewidth',1.2);hold on;
            grid on;
        end
        title(['KEEP: IC Activity PSD'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        ylabel('PSD (dB)', 'fontsize', 14); 
        xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        xlim([1 70]);
        legend(reshape(repmat(cellstr(num2str(IC_powpow)),2)',[],1),'Location','eastoutside');
        legend box off

        subplot(1,2,2)
        for i0 = 1:length(IC_powpow)
            % plot the fitted line
            p = IC_powpow(i0);
            plot(FREQ,spectra_psd(p,:),'color',c(floor((length(c)-20)/length(IC_powpow))*i0,:),'linewidth',1.2);hold on;
            grid on;
        end
        title(['KEEP: IC Activity PSD'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
        ylabel('PSD (dB)', 'fontsize', 14); 
        xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
        xlim([1 70]);
        legend(cellstr(num2str(IC_powpow)),'Location','eastoutside');
        legend box off
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['KEEP IC PSD.jpg']))
        
        summed_scores_to_keep = sum(classifications(:,1),2);
        titles_FEM = {};
        for i_title = 1:size(EEG.icawinv,2)
            titles_FEM{i_title} = [num2str(i_title) ': ' num2str(round(summed_scores_to_keep(i_title),3)),...
                ];
        end
        bemobil_plot_patterns_CL(EEG.icawinv,EEG.chanlocs,'chan_to_plot',find(IC_all_brain >= 8)','weights',summed_scores_to_keep,...
            'fixscale',0,'minweight',0.75,'titles',titles_FEM);
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['potential brain components_allcritera_customElectrode','.fig']))
        saveas(gcf,fullfile(output_folder,subjStr,'Figures',['potential brain components_allcritera_customElectrode','.jpg']))

    else
        load(fullfile(output_folder,subjStr,[subjStr,'_ICRej.mat']));
        IC_all_brain = Output_ICRejection.IC_all_brain;
    end
%     keyboard
    %% Last check Criteria 5: PowPow Cat Cross-Frequency Power-Power Coupling Analysis: A Useful Cross-Frequency Measure to Classify ICA-Decomposed EEG
    % It takes a long time to run. should pick only the ones that are
    % classified as 'brain'
    % %powpowcat parameters
    runPowPowCat = 1;
    if runPowPowCat
        upperFreqLimit = 100; %Frequency in Hz
        inputDataType = 2; %1, electrode data; 2, ICA time series
        methodType = 2;%1, Pearson's correlation; 2, Speaman's correlation
        numIterations = [];
        IC_powpow = find(IC_all_brain >= 8);
        fprintf('PowPowCAT parameters:\n upperFreqLimit= %i Hz\n inputDataType = ICs\n methodType= Spearman''s correlation (non-parametric)\n numIterations = %i\n',upperFreqLimit,numIterations);
        %run PowPowCAT
        % make a copy of EEG for powpowcat processing only
        if ~isempty(IC_powpow)
            EEG_powpow = EEG;
            EEG_powpow.icaact = EEG_powpow.icaact(IC_powpow,:);
            EEG_powpow = calc_PowPowCAT_CL(EEG_powpow, upperFreqLimit, inputDataType, methodType, numIterations);%CL version does not run stats
            EEG_powpow.setname = 'powpowcat';
            [ALLEEG, EEG_powpow, CURRENTSET] = eeg_store(ALLEEG, EEG_powpow , 0);
            eeglab redraw
            Plot35IcPushbutton_powpow(EEG_powpow,length(IC_powpow),IC_powpow)
            saveas(gcf,fullfile(output_folder,subjStr,['powpowcat.fig']))
        end
    end
    % eeglab redraw

    % SAVE IC_rejection 
    % EEG = pop_saveset( EEG, 'filepath', process_output_folder_IC_Rej, 'filename', strcat(subjStr,'_IC_rej'));
    

end

