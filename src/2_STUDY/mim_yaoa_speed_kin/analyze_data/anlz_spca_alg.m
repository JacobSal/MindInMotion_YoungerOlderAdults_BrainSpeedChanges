%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_yaoa_speed_kin/.sh


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
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        STUDY_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',e)
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
        STUDY_DIR = fileparts(SCRIPT_DIR);
    end
end
%## Add Study & Script Paths
addpath(STUDY_DIR);
cd(STUDY_DIR);
fprintf(1,'Current folder: %s\n',STUDY_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
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
