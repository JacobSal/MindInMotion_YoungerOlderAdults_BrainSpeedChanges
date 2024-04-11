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
%% ===================================================================== %%
% SUBJ_1YA = {'H1002','H1004','H1007','H1009',...
%     'H1010','H1011','H1012','H1013','H1017',...
%     'H1018','H1019','H1020','H1022','H1024',...
%     'H1025','H1026','H1027','H1029','H1030','H1031',...
%     'H1032','H1033','H1034','H1035','H1036',...
%     'H1037','H1038','H1039','H1041','H1042',...
%     'H1044','H1045','H1046','H1047','H1048'}; % JACOB,SAL (04/18/2023)
SUBJ_2MA = {'H2002','H2007','H2008',...
    'H2013','H2015','H2017','H2020','H2021',...
    'H2022','H2023','H2025','H2026','H2027',...
    'H2033','H2034','H2037','H2038','H2039',...
    'H2042','H2052','H2059','H2062','H2082',...
    'H2090','H2095','H2111','H2117'};
SUBJ_3MA = {'H3029','H3034','H3039','H3053',...
    'H3063','H3072','H3077','H3103',...
    'H3107',...
    'NH3006','NH3007','NH3008','NH3010','NH3021',...
    'NH3026','NH3030','NH3036','NH3040',...
    'NH3041','NH3043','NH3054',...
    'NH3055','NH3058','NH3059','NH3066',...
    'NH3068','NH3069','NH3070','NH3074',...
    'NH3076','NH3086','NH3090','H3092','NH3102',...
    'NH3104','NH3105','NH3106','NH3108','NH3110',...
    'NH3112','NH3113','NH3114','NH3123','NH3128',...
    };
SUBJ_SLOW_WALKERS = {'H3042','H3046','H3047','H3073',...
    'H3092','NH3025','NH3051','NH3056','NH3071','NH3082'};
SUBJ_NO_MRI = {'H2010','H2012','H2018','H2036','H2041',...
    'H2072','H3018','H3120','NH3002','NH3009','NH3027','NH3129'};
SUBJ_MISSING_COND = {'H3024','NH3028'};
SUBJ_DONT_INC = {'NH3004','NH3023'};
SUBJ_PICS = {SUBJ_2MA,SUBJ_3MA};
GROUP_NAMES = {'H2000''s','H3000''s'}; 
SUBJ_ITERS = {1:length(SUBJ_2MA),1:length(SUBJ_3MA)};
subject_chars = {};
for i = 1:length(SUBJ_PICS)
    subject_chars = [subject_chars, SUBJ_PICS{i}];
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
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',save_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',save_dir);
% end 

%%
% subject_chars = {ALLEEG.subject};
FREQS = (4:100);
N_FREQS = [];
WAVELET_STRUCT = struct('t',[0,1/EEG.srate],'f',FREQS,'fc',1,'FWHM_tc',3,'squared','n');
ANALYSIS_TYPE = 'channel';
EVENT_CHAR = 'RHS';
EPOCH_MIN_MAX = [3,4.25];
N_RESAMPLES = 100;
TIMEWARP_EVENTS = {'RHS','LHS','LTO','RTO'};
CONDITION_BASE = 'rest';
% CONDITIONS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
% CONDITIONS = {'flat','low','med','high'};
ALL_CONDITIONS = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
CONDITIONS = {'high'};
EPOCH_MIN_MAX = [1,4.25];
SUB_CHAN_SELECT = {'Cz';'LSSCM';'LISCM';'LSTrap';'LITrap';'RISCM';'RSSCM';'RITrap';'RSTrap'};
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
    if ~isempty(SUB_CHAN_SELECT)
        EEG_REST = pop_select(EEG_REST, 'channel', SUB_CHAN_SELECT);
    end
    baseline_data = permute(EEG_REST.data, [2,1]); % pnts x chans
    %- 
    inds = ismember({EEG.event.cond}, CONDITIONS);
    inds = find(inds);
    FROM = [EEG.event(inds(1)).latency];
    TO = [EEG.event(inds(end)).latency];
    EEG_GAIT = pop_select(EEG, 'point', [FROM; TO]');
    if ~isempty(SUB_CHAN_SELECT)
        EEG_GAIT = pop_select(EEG_GAIT, 'channel', SUB_CHAN_SELECT);
    end
%     EEG_GAIT = threshContinous(EEG_GAIT, cfg.V_tsh*2); % should this be extracted here?????
    EEG_GAIT.etc.valid_eeg = ones(size(EEG_GAIT.data,2),1);
    %- Maybe perform it on all dataset then perform timefreq decomp, use
    %specPCAdenoising plut derived coefficients on each gait cycle in
    %frequency domeain and project back to time domain after?
    [gait_avg,ERSP,GPM,TF_new,output_struct] = spca_time_freq_decomp(EEG_GAIT,baseline_data);
    [ERSP_corr, GPM_corr, PSC1, ~,COEFFs] = specPCAdenoising(ERSP);
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
    par_save(gait_ersp_struct,fpath,'gait_ersp_spca.mat');
    
    %##
    chan_i = 1;
%     DATA_TO_PLOT = {ERSP,GPM,ERSP_corr,GPM_corr};
%     DATA_TO_PLOT = {ERSP,GPM};
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
        saveas(fig,[fpath filesep sprintf('ERSP_step%i.jpg',i)])
    end
    %% DENOISE
    %- 
%     inds = ismember({EEG.event.cond}, ALL_CONDITIONS);
    inds = ismember({EEG.event.cond}, CONDITIONS);
    inds = find(inds);
    FROM = [EEG.event(inds(1)).latency];
    TO = [EEG.event(inds(end)).latency];
    EEG_GAIT = pop_select(EEG, 'point', [FROM; TO]');
    fprintf('Using channel data as default...\n');
    if ~isempty(SUB_CHAN_SELECT)
        EEG_GAIT = pop_select(EEG_GAIT, 'channel', SUB_CHAN_SELECT);
    end
    data = permute(EEG_GAIT.data, [2,1]); % pnts x chans
    n_comps = EEG_GAIT.nbchan;
    hs_min_max = EPOCH_MIN_MAX*EEG.srate;
    data = bsxfun(@minus, data, mean(data,2));
    %- time frequency transform
    [TF_new, morlet_params] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared);
    TF = TF_new;
    %-
    idx_hs = find(strcmp({EEG_GAIT.event.type}, EVENT_CHAR));
    gait_tf_out = zeros(size(data,1),n_comps,N_FREQS); %strides/trials x pnts x chans x freqs
%     gait_t = zeros(size(data));
    %- step counter, increased for each valid step
    iter = 1;
    cnt = 1;
    %- resample each stride to the same legth (100 pnts)
    for cycle_cnt = 1:length(idx_hs)-1
        %- find first and last sample of stride
        cycle_edge = round([EEG_GAIT.event(idx_hs(cycle_cnt)).latency,...
            EEG_GAIT.event(idx_hs(cycle_cnt+1)).latency-1]); % first and last frame of gait cycle
        %- labels of all events within this cycle
        cycle_event = {EEG_GAIT.event([idx_hs(cycle_cnt):idx_hs(cycle_cnt+1)]).type};
        %- only keep labels of gait events to check their order:
        cycle_gaitEvent = cycle_event(contains(cycle_event,TIMEWARP_EVENTS));
        %-
        if hs_min_max(1) <= cycle_edge(2)-cycle_edge(1) &&... % check time until next HS
                cycle_edge(2)-cycle_edge(1) <= hs_min_max(2) &&...
                all(ismember(TIMEWARP_EVENTS,cycle_gaitEvent)) % oder of gait events correct
            %- 
            tf_cycle = TF_new(cycle_edge(1):cycle_edge(2),:,:); % extract data
%             tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*N_FREQS); % reshape to be able to use the resample function, skip but resample over different dimension?
            gait_tf = reshape(tf_cycle,cycle_edge(2)-cycle_edge(1)+1,n_comps,N_FREQS);
            [ERSP_corr, GPM_corr, PSC1, ~,COEFFs] = specPCAdenoising(gait_tf,COEFFs);
            TF_new(cycle_edge(1):cycle_edge(2),:,:) = ERSP_corr;
            cnt = cnt+1;
%             iter = iter + N_RESAMPLES;
        end
    end
    
    disp([num2str(round(cnt/cycle_cnt*100)) '% of the gait cycles are valid'])
    %## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
%     GPM = bsxfun(@minus,ERDS,mean(ERDS));
    inv_sig = icwt(squeeze(TF_new(:,1,:))');
    inv_orig = icwt(squeeze(TF(:,1,:))');
    %-
    cls = linspecer(3);
    fig = figure();
    hold on; 
    plot(data(:,1),'DisplayName', 'orig', 'Color', [cls(1,:), 0.75]);  
    % plot(inv_sig(:,1),'DisplayName', 'new'); 
    plot(inv_orig(:),'DisplayName', 'inverse_orig', 'Color', [cls(2,:), 0.3]); 
    plot(inv_sig(:),'DisplayName', 'clean', 'Color', [cls(3,:), 0.3]);  
    hold off;
    legend;
    saveas(fig,[fpath filesep 'sig_fig.jpg'])
    %- visualize
   
end

%%
%-
% inv_sig = ifft(squeeze(TF(:,1,:)));
inv_sig = icwt(squeeze(TF(:,1,:))','amor',flip(WAVELET_STRUCT.f),[WAVELET_STRUCT.f(1),WAVELET_STRUCT.f(end)],'WaveletParameters',[1,3]);
% inv_sig = inv_morlet_transform(TF,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc);
% inv_sig = (1/size(TF,3))*sum(squeeze(real(TF(:,1,:))),2);
figure; hold on; 
plot(data(:,1),'DisplayName', 'old'); 
plot(inv_sig(:,1),'DisplayName', 'new'); 
% plot(inv_sig(:),'DisplayName', 'new'); 
hold off;
legend;
%%
%{
TF = zeros(size(data,1),size(data,2),WAVELET_STRUCT.f(end)-WAVELET_STRUCT.f(1)+1);
inv_TF = zeros(size(data,1),size(data,2));
for chan_i = 1:size(data,2)
    [tmp,period,scale,coi,dj,para,k,W] = contwt(data(:,chan_i),WAVELET_STRUCT.t(2),1,1,WAVELET_STRUCT.f(1),WAVELET_STRUCT.f(end)-WAVELET_STRUCT.f(1),'Morlet',WAVELET_STRUCT.FWHM_tc);
    TF(:,chan_i,:) = tmp';
    inv_TF(:,chan_i) = invcwt(tmp,'Morlet',scale,para,k);
end
%-
for s = 1:length(scale)
    TF(:,:,s) = TF(:,:,s) * 2/sum(abs(W{s}));
end
TF = abs(TF);
allersp = squeeze(TF(1:100,1,:));
% allersp = squeeze(mtf_TF(1:100,1,:));
% tmp = p_struc.P_ifft;
% allersp = squeeze(tmp(1:100,1,:));
% allfreqs = f;
allfreqs = WAVELET_STRUCT.f(1):WAVELET_STRUCT.f(end);
figure;
tftopo(allersp,1:size(allersp,1),allfreqs,...
            'logfreq','native');
%}
