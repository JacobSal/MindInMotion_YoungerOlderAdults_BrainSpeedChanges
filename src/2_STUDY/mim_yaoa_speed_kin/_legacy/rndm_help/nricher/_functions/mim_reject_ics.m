function [Output_ICRejection] = mim_reject_ics(EEG,save_dir,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- ICLabel version 
ICLABEL_VERSION = 'lite';
%- 
DO_POWPOWCAT = true;
DO_POWPOW_PLOT = false;
upperFreqLimit = 100; %Frequency in Hz
inputDataType = 2; %1, electrode data; 2, ICA time series
methodType = 2;%1, Pearson's correlation; 2, Speaman's correlation
numIterations = 1000;
%-
fit_range = 2:40; % Original setting 2:100. now changed to 2:40 (2022-1-2)
classes_to_keep = 1; % Brain components
thresholds_keep = [0.5,0.75,0.9]; % percentages to threshold the ICLabel brain label by for keeping
thresholds_remove = [0.5]; % percentages to threshold the ICLabel muscle and eye label by for removing
slope_thres = -0.2;  % Added on 1-3-2022 to avoid very flat spectrum that is not brain component
%-
SPEC_PERC = 80;
SPEC_FREQRANGE = [2 100];
%-
DO_PLOTTING = true;
%-
THRESHOLD_RV_BRAIN = 0.15;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,EEG,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## RECALCULATE ICAACT MATRICES
EEG = eeg_checkset(EEG,'loaddata');
if isempty(EEG.icaact)
    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

%## (Crteria 1) Count brain components based on ICLabel
try 
    EEG.etc.ic_classification.ICLabel.classifications;
catch 
    EEG = iclabel(EEG,ICLABEL_VERSION); % NOTE: iclabel classes: 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
end
EEG.etc.ic_classification.ICLabel.classifications_default = EEG.etc.ic_classification.ICLabel.classifications;
classifications = EEG.etc.ic_classification.ICLabel.classifications;
%## REJECT BASED ON ICLABELS
comps_keep = cell(1,length(thresholds_keep));
comps_remove = cell(2,length(thresholds_remove));
%- loop
for t_i = 1:length(thresholds_keep)
    tmp = sum(classifications(:,classes_to_keep),2);
    comps_keep{t_i} = find(tmp>thresholds_keep(t_i));
end
%- loop
for t_i = 1:length(thresholds_remove)
    tmp = sum(classifications(:,2),2);%Muscle
    comps_remove{1,t_i} = find(tmp>thresholds_remove(t_i));
    tmp = sum(classifications(:,3),2);%eye
    comps_remove{2,t_i} = find(tmp>thresholds_remove(t_i));
end
%% ===================================================================== %%
%## (Criteria 2) Spectrum graph: Plot PSD together for selected channels
IC_potential = 1:size(EEG.icawinv,2);
icaacttmp = EEG.icaact(IC_potential, :, :);
%- compute the spectrum plots: TO DO: need to run it faster
%output unit from spectopo is in db
[spectra_psd,FREQ] = spectopo( icaacttmp, EEG.pnts, EEG.srate,'percent',SPEC_PERC,...
       'freqrange',SPEC_FREQRANGE,'plot','off');
%- add a linear fit: note that the frequency used here is log(frequency)
%to help characterize the 1/f structure. y axis should be log(power)
lsfit = [];
spectra_psd_fit = [];
for i = 1:length(IC_potential)
    lsfit(i,:) = polyfit(log10(FREQ(fit_range,1)), spectra_psd(i,fit_range)', 1);
    spectra_psd_fit(i,:) = lsfit(i,1).*log10(FREQ)+lsfit(i,2);
end
%- IC rejection criteria: only keep those with log fit < 0
ICs_keep_brain_spec = IC_potential(lsfit(:,1)<slope_thres)';% actually not a lot get discarded
ICs_spec_dump = IC_potential(lsfit(:,1)>=slope_thres)';
numICs = 1:size(EEG.icawinv,2);
%% ===================================================================== %%
%## (Criteria 4) Scalp topographs and dipole location (if outside the brain or not)
% Residual variance < 0.15
IC_RV = vertcat(EEG.dipfit.model.rv);
IC_POSXYZ = vertcat(EEG.dipfit.model.posxyz);
ICs_RVthreshold_keep = find(IC_RV <= THRESHOLD_RV_BRAIN & all(~isnan(IC_POSXYZ),2));

%% ===================================================================== %%
%## Gather all IC rejection criteria
%-
All_IC_criteria.IClabel.brain50 = comps_keep{1};
All_IC_criteria.IClabel.brain75 = comps_keep{2};
All_IC_criteria.IClabel.brain90 = comps_keep{3};
%-
All_IC_criteria.IClabel.brain = comps_keep{1};
All_IC_criteria.IClabel.muscle = comps_remove{1,1};
All_IC_criteria.IClabel.muscle_threshold = thresholds_remove(1);
All_IC_criteria.IClabel.eye = comps_remove{2,1};
All_IC_criteria.IClabel.eye_threshold = thresholds_remove(1);
All_IC_criteria.IClabel.classification = EEG.etc.ic_classification.ICLabel.classifications;
%-
All_IC_criteria.Spectra.keep = ICs_keep_brain_spec;
All_IC_criteria.Spectra.dump = ICs_spec_dump;
%-
All_IC_criteria.RV.keep = ICs_RVthreshold_keep;
All_IC_criteria.RV.IC_RV = IC_RV;
EEG.etc.IC_rej = All_IC_criteria;

%- Assign weights to each criteria 
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

%%
%## (Criteria 5)  PowPowCat
% PowPow Cat Cross-Frequency Power-Power Coupling Analysis: A Useful Cross-Frequency Measure to Classify ICA-Decomposed EEG
% It takes a long time to run. should pick only the ones that are
% classified as 'brain'

if DO_POWPOWCAT
%     IC_powpow = find(Output_ICRejection.IC_all_brain >= 8);
    fprintf('PowPowCAT parameters:\n upperFreqLimit= %i Hz\n inputDataType = ICs\n methodType= Spearman''s correlation (non-parametric)\n numIterations = %i\n',upperFreqLimit,numIterations);
    %run PowPowCAT
    IC_powpow = find(IC_all_brain >= 7);
%     IC_powpow = find(IC_all_brain >= 5);
    % make a copy of EEG for powpowcat processing only
    if ~isempty(IC_powpow)
        EEG_powpow = EEG;
        EEG_powpow.icaact = EEG.icaact(IC_powpow,:);
%         EEG_powpow = pop_subcomp(EEG,IC_powpow,0,1);
        EEG_powpow = calc_PowPowCAT(EEG_powpow,upperFreqLimit,inputDataType,methodType,numIterations);
        EEG_powpow.setname = 'powpowcat';
        Plot35IcPushbutton_powpow(EEG_powpow,length(IC_powpow),IC_powpow)
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_dir filesep 'powpowcat.fig']);
        saveas(fig_i,[save_dir filesep 'powpowcat.jpg']);
        %- save EEG powpowcat?
        inds = PowPowCat_ICrej(EEG_powpow,DO_POWPOW_PLOT,save_dir);
        badPPC_IC = IC_powpow(inds);
    else
        badPPC_IC = [];
    end
else
    badPPC_IC = [];
end

%% SAVE OUTPUT
Output_ICRejection.All_IC_criteria = All_IC_criteria;
Output_ICRejection.IC_all_brain = IC_all_brain;
Output_ICRejection.IC_all_muscle = IC_all_muscle;
Output_ICRejection.IC_all_eye = IC_all_eye;
Output_ICRejection.IC_powpow_rej = badPPC_IC;
Output_ICRejection.Cleaning_Params = EEG.etc.Params;
%- save Output_ICRejection
fileName = [save_dir filesep sprintf('%s_ICRej.mat',EEG.subject)];
save(fileName,'Output_ICRejection');
%% ===================================================================== %%
%## (PLOTS) INSPECTION
if DO_PLOTTING
    %## (PLOT 1) Individual Independent Components
    %{
    N = 10;
    CMAP = linspecer(N);%my personal color scheme
    figure('color','w');
    hold on;
    for i = 1:N
        plot(FREQ,spectra_psd(i,:),'color',CMAP(i,:),'linewidth',1.2);hold on;
    end
    hold off;
    title('IC Activity Power Spectrum','units','normalized','fontsize',14,'FontWeight','Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([.5 80]);
    legend(cellstr(num2str(1:10)),'Location','eastoutside');
    legend box off
    fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = [500 300 1080 720];
%     saveas(fig_i,[save_dir filesep sprintf('%s_IC_powerspecs.fig',EEG.subject)]);
    saveas(fig_i,[save_dir filesep sprintf('%s_IC_powerspecs.jpg',EEG.subject)]);
    %}
    %## (PLOT 2) stemplot
    figure();
    hold on;
    stem(numICs,lsfit(:,1));
    if any(lsfit(:,1)> slope_thres)
        stem(numICs(lsfit(:,1)> slope_thres),lsfit(lsfit(:,1)>slope_thres,1),'r');
    end
    hold off;
    ylabel('slope')
    fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep sprintf('spectral_stem.jpg')]);
    %## (PLOT 3) POTENTIAL KEEP COMPONENTS
    % Log-log plot - retained
    CMAP = linspecer(length(ICs_keep_brain_spec));%my personal color scheme 
    figure('color','w');
    subplot(1,2,1)
    hold on;
    for i0 = 1:length(ICs_keep_brain_spec)
        % plot the fitted line
        p = ICs_keep_brain_spec(i0);
        semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',CMAP(i0,:),'linewidth',1.2);hold on;
        semilogx(FREQ,spectra_psd(p,:),'color',CMAP(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title(['KEEP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(ICs_keep_brain_spec)),'Location','eastoutside');
    legend box off
    subplot(1,2,2)
    hold on;
    for i0 = 1:length(ICs_keep_brain_spec)
        % plot the fitted line
        p = ICs_keep_brain_spec(i0);
        plot(FREQ,spectra_psd(p,:),'color',CMAP(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title(['KEEP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(ICs_keep_brain_spec)),'Location','eastoutside');
    legend box off
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep 'keep_potential_brain_components_spectral.jpg']);

    %## (PLOT 4) POTENTIAL REJECT COMPONENTS
    CMAP = linspecer(length(ICs_spec_dump));%my personal color scheme 
    figure('color','w');
    hold on;
    for i0 = 1:length(ICs_spec_dump)
        % plot the fitted line
        p = ICs_spec_dump(i0);
        semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',CMAP(i0,:),'linewidth',1.2);hold on;
        semilogx(FREQ,spectra_psd(p,:),'color',CMAP(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    hold off;
    title(['DUMP: IC Activity Power Spectrum'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('Power 10*log_{10}(uV^2/Hz)', 'fontsize', 14); 
    xlabel('Log(Frequency) (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(ICs_spec_dump)),'Location','eastoutside');
    legend box off
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    %     saveas(gcf,fullfile(save_IC_Rejection_folder,subjStr,'Figures',['potential brain components_spectral_',subDirNum,'.fig']))
    saveas(fig_i,[save_dir filesep 'reject_potential_brain_components_spectral.jpg'])
    
    %## (PLOT 5) STEM
    figure('color','white');
    subplot(3,1,1);
    stem(numICs,IC_all_muscle);title(['muscle: ',num2str(sum(IC_all_muscle == 3))]);ylabel('Dump: score');
    subplot(3,1,2);
    stem(numICs,IC_all_eye);title(['eye: ',num2str(sum(IC_all_eye == 2))]);ylabel('Dump: score');
    subplot(3,1,3);hold on;
    stem(numICs,IC_all_brain);title(['brain: ',num2str(sum(IC_all_brain == 8))]);ylabel('Keep: score');
    stem(numICs(IC_all_muscle>=2),IC_all_brain(IC_all_muscle>=2),'r');
    xlabel('IC number');legend('','flag muscle');
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep 'IC_score.jpg'])

    %## (PLOT 6) LOG-LOG PLOT 
    % Log-log plot - retained
    IC_powpow = find(IC_all_brain >= 8);
    CMAP = linspecer(length(IC_powpow));
    figure('color','w');
    subplot(1,2,1)
    for i0 = 1:length(IC_powpow)
        % plot the fitted line
        p = IC_powpow(i0);
        semilogx(FREQ,lsfit(p,1).*log10(FREQ)+lsfit(p,2),'--','color',CMAP(i0,:),'linewidth',1.2);hold on;
        semilogx(FREQ,spectra_psd(p,:),'color',CMAP(i0,:),'linewidth',1.2);hold on;
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
        plot(FREQ,spectra_psd(p,:),'color',CMAP(i0,:),'linewidth',1.2);hold on;
        grid on;
    end
    title(['KEEP: IC Activity PSD'],'units','normalized', 'fontsize', 14, 'FontWeight', 'Normal');
    ylabel('PSD (dB)', 'fontsize', 14); 
    xlabel('Frequency (Hz)', 'fontsize', 14, 'fontweight', 'normal'); 
    xlim([1 70]);
    legend(cellstr(num2str(IC_powpow)),'Location','eastoutside');
    legend box off
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = [500 300 1080 720];
    saveas(fig_i,[save_dir filesep 'KEEP_IC_PSD.jpg']);

    %## (PLOT 7) BEMOBIL PIPELINE PLOT
    summed_scores_to_keep = sum(classifications(:,1),2);
    titles_FEM = {};
    for i_title = 1:size(EEG.icawinv,2)
        titles_FEM{i_title} = [num2str(i_title) ': ' num2str(round(summed_scores_to_keep(i_title),3))];
    end
    bemobil_plot_patterns_CL(EEG.icawinv,EEG.chanlocs,'chan_to_plot',find(IC_all_brain >= 8)','weights',summed_scores_to_keep,...
        'fixscale',0,'minweight',0.75,'titles',titles_FEM);
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep 'potential_brain_components_allcritera_customElectrode.fig']);
    saveas(fig_i,[save_dir filesep 'potential_brain_components_allcritera_customElectrode.jpg']);
    
    %## (PLOT 8) POWPOWCAT
    

end
end

