function bad_ics_out = PowPowCat_ICrej(EEG,do_plot,save_dir,varargin)
% ADAPTED BY CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% PowPowCat for IC rejection
% Reject ICs that have moderate cross-freq coupling in low (<8Hz) and high
% (>30Hz) frequency windows
% Takes the median correlation coefficient value in both high and low
% frequency windows, not including the identity coeffient (which always
% ==1)
%
% % Usage:
% >> [badPCC_IC] = PowPowCat_ICrej(EEG, varargin);
%
% Required inputs:
%   EEG            - EEG dataset structures with PowPowCat precomputed
%
% Optional inputs:
%   plotstuff      - plot extended component properties of "bad" ICs and
%                    spectral covariance (default = 0 (off))
%   outputFigureFolder - main folder path to store figures, subfolder will
%                       be created. (default = EEG.filepath)
% Output:
%   badPPC_IC      - Independent components identified as having moderate
%                    correlation in low and high freq. windows (e.g. eye
%                    blinks and muscle artifacts, respectively)
% Code Designers: Noelle Jacobsen, Jacob Salminen
% Code Date: 10/17/2021, MATLAB 2020b
% Authored (10/17/2021): Noelle Jacobsen, University of Florida, 
% edited (07/15/2023): Jacob Salminen, University of Florida; cleaned up 

%## TIME
tic
%## DEFINE DEFAULTS
%- Cutoff for Correlations ot determine bad components
% NOTE: correlation coefficient threshold, 0.5-0.7= moderate correlation
BAD_IC_CUTOFF = 0.3; 
% NOTE (10/25/2023), Originally at 0.3; setting to 0.25 to be more
% aggressive on cleaning.
% NOTE (10/26/2023), Originally at 0.25; setting to 0.20 to be more
% aggressive on cleaning.
% NOTE (11/2/2023), Turning back to 0.3 for testing on new ICC
% parameterization.
% (11/7/2023), Turning to 0.1 & only performing on high frequency.
%-
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'do_plot',@islogical);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,EEG,do_plot,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%- extract powpowcat values
cov_matrix = EEG.etc.PowPowCAT.covMatrix;
freqs =  EEG.etc.PowPowCAT.freqs;
%- set params
lowfreq_coupling=cell(size(cov_matrix,3),1);
highfreq_coupling=cell(size(cov_matrix,3),1);
bad_ics_out=zeros(size(cov_matrix,3),1);
fprintf('\tIC#\t\tLowFreq Cov.\t\tHighFreq Cov.\n');
fprintf('________________________________________________');
for IC = 1:size(cov_matrix,3)
    %- extract correlations above and below alpha/beta frequencies (<8Hz,>30Hz)
    lowfreq_coupling{IC}=cov_matrix(freqs<8,freqs<8,IC); %below alpha
    highfreq_coupling{IC}=cov_matrix(freqs>30,freqs>30,IC); %above beta
    %## Low Freq Coupling
    CC = lowfreq_coupling{IC};
    %- remove diagnoal identity (coeff always equals one with itself)
    CC(logical(eye(size(CC)))) =[]; 
    CC=reshape(CC,size(lowfreq_coupling{IC},1)-1,size(lowfreq_coupling{IC},1));
    %- average (median)
    lowfreq_coupling_avg = median(median(CC));
    %## High Freq Coupling
    CC = highfreq_coupling{IC};
    %- remove diagnoal identity (coeff always equals one with itself)
    CC(logical(eye(size(CC)))) =[]; 
    CC=reshape(CC,size(highfreq_coupling{IC},1)-1,size(highfreq_coupling{IC},1));
    %- average (median)
    highfreq_coupling_avg = median(median(CC));
    %## TABLE
    fprintf('\n\t%i\t\t%.2f\t\t\t\t%.2f',IC,lowfreq_coupling_avg,highfreq_coupling_avg)
    %identify ICs with corr. coeff. above thresh in low and high freq windows
%     if lowfreq_coupling_avg>BAD_IC_CUTOFF || highfreq_coupling_avg>BAD_IC_CUTOFF
    if highfreq_coupling_avg>BAD_IC_CUTOFF
        bad_ics_out(IC) = IC;
        fprintf('\t**BAD**');
    end
end
bad_ics_out=bad_ics_out(bad_ics_out~=0);

%% (PLOTTING) ========================================================== %%
if do_plot
    save_dir_1 = [save_dir filesep 'summary'];
     if ~exist([save_dir filesep 'summary'],'dir')
        mkdir([save_dir filesep 'summary'])
     end
    save_dir_2 = [save_dir filesep 'powpowcat'];
    if ~exist([save_dir filesep 'powpowcat'],'dir')
        mkdir([save_dir filesep 'powpowcat'])
    end
    %## (PLOT 1) EXTENDED COMPONENT PROPERTIES
    for badIC = 1:length(bad_ics_out)
        fprintf('\nPlotting extended component properties of bad ICs');
        % plot extended comp properties
        pop_prop_extended(EEG,0,bad_ics_out(badIC),NaN,{'freqrange', [2 80]},{},1,'')
        fig_i = get(groot,'CurrentFigure');
        saveas(fig_i,[save_dir_2 filesep sprintf('%s_properties_%i.jpg',EEG.subject,bad_ics_out(badIC))]);
        close(fig_i);
    end
    %## (PLOT 2)
    %{
    myfig=figure;
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
    set(gcf,'Color','w');
    mytitle= sgtitle([EEG.subject,' Spectral Covariance ']);
    for IC = 1:size(cov_matrix,3)
        if IC<=36
            subplot(6,6,IC)
        elseif (IC == 37 || IC == 73 || IC== 109)
            saveas(myfig,[save_dir filesep 'PowPowCAT' filesep sprintf('%s_SpectralCovariance_%i.jpg',EEG.subject,IC-1)]);
            myfig= figure;%create a new figure so you don't have all ICs on one plot
            set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
            set(gcf,'Color','w')
            mytitle= sgtitle({'Component projection to inferior head and neck electrodes',...
                extractBefore(EEG.filename,'.set')});
            mytitle.Interpreter = 'none';
            close(myfig);
        end

        imagesc(cov_matrix(:,:,IC), [-0.8 0.8]);
        customColorMap = colormap(jet);
        colormap(customColorMap)
        currentAxesPosition = get(gca, 'position');
        colorbarHandle = colorbar;
        set(get(colorbarHandle, 'title'), 'String', 'Corr. Coef','fontsize',8)
        set(colorbarHandle, 'fontsize',8);
        axis xy
        axis square
        tickLabels = round(freqs(10:10:length(freqs))*10)/10;
        tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
        set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', tickLabels,...
            'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
            'fontsize', 5)
        xlabel('Frequency (Hz)', 'fontsize', 8)
        ylabel('Frequency (Hz)', 'fontsize', 8)
        if ismember(IC,bad_ics_out)
            title(['IC' int2str(IC)],'Color','r','fontsize',10)
        else
            title(['IC' int2str(IC)],'fontsize',10)
        end
    end
    saveas(myfig,[save_dir_2 filesep sprintf('%s_SpectralCovariance_%i.jpg',EEG.subject,IC-1)]);
    fprintf('\nFigures stored here: %s\n',[save_dir filesep 'PowPowCAT']);
    close(myfig);
    %}
    %## (PLOT 3)
    %{
    %- plot bad ICs spectral covariance subplot dimensions
    if length(bad_ics_out)>5
        M = 2;
        N = 5;
    else
        M = 1;
        N = length(bad_ics_out);
    end
    myfig = figure;
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
    set(gcf,'Color','w');
    mytitle= sgtitle([EEG.subject,' Spectral Covariance -Bad ICs']);
    for IC = 1:length(bad_ics_out)
        if length(bad_ics_out)>1
            subplot(M,N,IC)
        end
        imagesc(cov_matrix(:,:,bad_ics_out(IC)), [-0.8 0.8]);
        customColorMap = colormap(jet);
        colormap(customColorMap)
        currentAxesPosition = get(gca, 'position');
        colorbarHandle = colorbar;
        set(get(colorbarHandle, 'title'), 'String', 'Corr. Coef','fontsize',8)
        set(colorbarHandle, 'fontsize',12);
        axis xy
        axis square
        tickLabels = round(freqs(10:10:length(freqs))*10)/10;
        tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
        set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', tickLabels,...
            'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
            'fontsize', 8)
        xlabel('Frequency (Hz)', 'fontsize',12)
        ylabel('Frequency (Hz)', 'fontsize', 12)
        title(['**IC' int2str(IC),'**'],'Color','r','fontsize',16)
    end
    saveas(myfig,[save_dir_1 filesep sprintf('%s_BadSpectralCovariance_%i.jpg',EEG.subject,IC-1)]);
    fprintf('Summary of ICs with bad spectral covariance stored here: %s\n',[save_dir,'\Summary']);
    close(myfig);
    %}
end
%## TIME
toc
end

