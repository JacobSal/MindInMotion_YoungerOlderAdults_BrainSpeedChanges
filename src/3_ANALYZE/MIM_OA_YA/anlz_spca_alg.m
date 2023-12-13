%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_anlz_emg_artifact.sh

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
%- define the directory to the src folderd
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
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
%## hard define
%- datset name
DATA_SET = 'MIM_dataset';
OA_PREP_FPATH = 'EMG_ANALYSIS';
dt = 'tmp_emg_analysis';
DATA_DIR = [source_dir filesep '_data'];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep '_studies' filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)
study_fName = 'epoch_study';
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
if ~ispc
    [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',save_dir);
else
    [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',save_dir);
end     

%%
subject_chars = {ALLEEG.subject};
FREQS = (4:100);
N_FREQS = [];
WAVELET_STRUCT = struct('t',[0,1/EEG.srate],'f',[],'fc',1,'FWHM_tc',3,'squared','n');
ANALYSIS_TYPE = 'channel';
EVENT_CHAR = 'RHS';
EPOCH_MIN_MAX = [3,4.25];
N_RESAMPLES = 100;
TIMEWARP_EVENTS = {'RHS','LHS','LTO','RTO'};
CONDITION_BASE = 'rest';
% CONDITIONS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% CONDITIONS = {'flat','low','med','high'};
CONDITIONS = {'high'};
%##
cond_i = 1;
%%
for subj_i = 1:length(ALLEEG)
    %-
    fpath = [OUTSIDE_DATA_DIR filesep subject_chars{subj_i} filesep 'clean'];
    fname = dir([fpath filesep '*.set']);
    EEG = pop_loadset('filepath',fpath,'filename',fname(1).name);
    [EEG_chans,EMG_chans,Noise_chans] = getChannelTypes_func(EEG);
    EEG = pop_select(EEG, 'channel', [EEG_chans, EMG_chans]);
    %-
    inds1 = logical(strcmp({EEG.event.cond}, CONDITION_BASE));
    inds2 = logical(strcmp({EEG.event.type}, 'boundary'));
    val_inds = find(inds1 & ~inds2);
    FROM = [EEG.event(val_inds(1)).latency];
    TO = [EEG.event(val_inds(end)).latency];
    EEG_REST = pop_select(EEG, 'point', [FROM; TO]');
    baseline_data = permute(EEG_REST.data, [2,1]); % pnts x chans
    %- 
    inds = ismember({EEG.event.cond}, CONDITIONS);
    inds = find(inds);
    FROM = [EEG.event(inds(1)).latency];
    TO = [EEG.event(inds(end)).latency];
    EEG_GAIT = pop_select(EEG, 'point', [FROM; TO]');
%     EEG_GAIT = threshContinous(EEG_GAIT, cfg.V_tsh*2); % should this be extracted here?????
    EEG_GAIT.etc.valid_eeg = ones(size(EEG_GAIT.data,2),1);
    %- Maybe perform it on all dataset then perform timefreq decomp, use
    %specPCAdenoising plut derived coefficients on each gait cycle in
    %frequency domeain and project back to time domain after?
    [gait_avg,ERSP,GPM,output_struct] = spca_time_freq_decomp(EEG_GAIT,baseline_data);
    [ERSP_corr, GPM_corr, PSC1, ~,~] = specPCAdenoising(ERSP);
    %## PLOT
    figure(); set(gcf, 'position', [0 0 600 500]);
    plot(FREQS, squeeze(output_struct.baseline_ersp)', 'k-');
    ylabel('Amplitude (\muV)'), ylabel('Frequency (Hz)');
    grid on; box off
    title('Baseline ERSP (rest)');
    %- save all info together
    gait_ersp_struct = [];
    gait_ersp_struct.ID         = EEG_GAIT.subject;
    gait_ersp_struct.Noise_cov  = output_struct.baseline_cov;% noise cov for kernel computation
    gait_ersp_struct.F_Rest     = output_struct.baseline_ersp;
    gait_ersp_struct.TF         = gait_avg;
    gait_ersp_struct.ERSP_uncor = ERSP;
    gait_ersp_struct.GPM_uncor  = GPM;
    gait_ersp_struct.ERSP       = ERSP_corr;
    gait_ersp_struct.GPM        = GPM_corr;
    gait_ersp_struct.PSC1       = PSC1;
    gait_ersp_struct.numStrides = output_struct.cycle_cnt;
    gait_ersp_struct.numValidStrides = output_struct.valid_cycle_cnt;
    gait_ersp_struct.chanlocs   = EEG_GAIT.chanlocs;
    par_save([fpath filesep 'gait_ersp_spca.mat'],'gait_ersp_struct');
    
    %##
    chan_i = 1;
    DATA_TO_PLOT = {ERSP,GPM,ERSP_corr,GPM_corr};
    DATA_CHARS = {'original ERSP','orignial GPM','corrected ERSP','corrected GPM'};
    FIGURE_POSITION = [100,100,350,350];
%     EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    XTICK_LABEL = 'Gait Events';
    YTICK_LABEL = 'Frequency (Hz)';
    clim_ersp = [-2,2];
    FONT_SIZE = 12;
    %##
    for i = 1:length(DATA_TO_PLOT)
        in_dat = DATA_TO_PLOT{i};
        allersp = squeeze(in_dat(:,chan_i,:));
        alltimes = (1:size(allersp,1));
        allfreqs = FREQS;
        SUB_FREQ_LIMS = [min(allfreqs), max(allfreqs)];
        fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
        set(fig,'Units','inches','Position',[3 3 5 5])
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        horiz_shift = 0;
        tftopo(allersp,alltimes,allfreqs,'limits',... 
            [0 100 nan nan clim_ersp],...
            'logfreq','native');
        ax = gca;
        hold on;
        colormap(linspecer);
        %- adjust subplot position and height
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
    %         disp(get(ax,'Position'));
        %- set ylims
        ylim(log(SUB_FREQ_LIMS))
        if SUB_FREQ_LIMS(2) <= 50
            set(ax,'YTick',log([4.01,8,13,30,50])); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end  
        %- set color lims
        set(ax,'clim',clim_ersp);
        %- set x-axis & y-axis labels
        ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
        xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        colorbar
        hold off;
        fig = get(groot,'CurrentFigure');
    end
end
