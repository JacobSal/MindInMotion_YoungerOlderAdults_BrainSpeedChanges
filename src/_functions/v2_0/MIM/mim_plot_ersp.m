function [] = mim_plot_ersp(STUDY,ALLEEG,warping_times,save_dir,varargin)
%MIM_PLOT_ERSP Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%-
XTICK_LABEL = 'Gait Events';
YTICK_LABEL = 'Frequency (Hz)';
FONT_SIZE = 12;
SUBPLOT_WIDTH = 0.15;
SUBPLOT_HEIGHT = 0.7;
SHIFT_AMNT = 0.175;
STATS_TITLE = 'CUSTOM STATS';
%-
ersp_plot_type = 'tftopo';
ersp_pcond = [];

%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ersp_raw',@isnumeric);

addRequired(p,'ersp_raw',@isnumeric);
addRequired(p,'ersp_raw',@isnumeric);
addRequired(p,'ersp_raw',@isnumeric);
addRequired(p,'ersp_raw',@isnumeric);
addRequired(p,'ersp_raw',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'ersp_plot_type',ersp_plot_type,@ischar);
addParameter(p,'ersp_pcond',ersp_pcond,@isnumeric);
parse(p,STUDY,ALLEEG,warping_times,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
ersp_plot_type = p.Results.ersp_plot_type;
ersp_pcond = p.Results.ersp_pcond;
%% ===================================================================== %%
switch ersp_plot_type
    case 'tftopo'
        plot_tftopo(
    case 'contourf'
        
end

end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [fig] = plot_tftopo(allersp,alltimes,allfreqs,alltitles,allpcond,warping_times,clim_max,colormap_ersp,...
    EVENT_CHARS,SUB_FREQ_LIMS,FIGURE_POSITION)
    if length(clim_max) == 2
        clim_ersp = clim_max;
    else
        clim_ersp = [-clim_max,clim_max];
    end
    %%
    fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    set(fig,'Units','inches','Position',[3 3 14 5])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    horiz_shift = 0;
    hold on;
    for j = 1:length(allersp)
        subplot(1,length(allersp)+1,j); %,'position',[0.01+horiz_shift,0.1,0.5,0.5])
        ax = gca;
        tftopo(allersp{j},alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan clim_ersp],...
            'logfreq','native');
        hold on;
        colormap(colormap_ersp);
        %- adjust subplot position and height
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
%         disp(get(ax,'Position'));
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
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
        if j == 1
            ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
            xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        else
            xlabel('','FontSize',FONT_SIZE);
            ylabel('','fontsize',FONT_SIZE,'fontweight','bold');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
%         ylim()
        %- title
        title(alltitles{j});  
        horiz_shift = horiz_shift + SHIFT_AMNT;
    end
%     hold off;
    %%
    %## Add Stats To Plot
    if ~isempty(allpcond)
        subplot(1,length(allersp)+1,length(allersp)+1) % add one subplot for stats
        tftopo(double(allpcond),alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan  clim_ersp],...
            'logfreq','native')
        colormap(colormap_ersp);
        ax = gca;
        %-
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
        disp(get(ax,'Position'));
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.04;
        c.Limits = clim_ersp;
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',FONT_SIZE);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        set(hL,'Rotation',0);
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
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
        %- set y-axis labels
        xlabel('','FontSize',FONT_SIZE);
        ylabel(sprintf(''),'fontsize',FONT_SIZE,'fontweight','bold');
        %- set x-axis labels
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- title
        title(STATS_TITLE)
    else
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.05;
        c.Limits = clim_ersp;
    end
    hold off;
    fig = get(groot,'CurrentFigure');
end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [fig] = plot_contourf(ersp_raw,ersp_pcond,ersp_masked,alltimes,allfreqs,...
        alltitles,warping_times,clim_max,colormap_ersp,PANEL_OFFSET,EVENT_CHARS,SUB_FREQ_LIMS,FIGURE_POSITION)
    XTICK_LABEL = 'Gait Events';
    YTICK_LABEL = 'Frequency (Hz)';
    FONT_SIZE = 12;
    SUBPLOT_WIDTH = 0.15;
    SUBPLOT_HEIGHT = 0.7;
    SHIFT_AMNT = 0.175;
%     STATS_TITLE = 'CUSTOM STATS';
    COLORBAR_SHIFT = 0.05;
    ALPHA_MULTIPLE = 0.75;
    %%
    fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    set(fig,'Units','inches','Position',[3 3 14 5])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    horiz_shift = 0;
    hold on;
    for j = 1:length(ersp_raw)
        subplot(1,length(ersp_raw),j)
        alpha_mask = ones(size(ersp_pcond{j}))*ALPHA_MULTIPLE; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        alpha_mask(ersp_pcond{j} == 1) = 0; %0 is significant? 1 is not?
        contourf(alltimes, allfreqs, ersp_raw{j},200,...
                   'linecolor','none');
        hold on;
        imagesc(alltimes,allfreqs,ersp_masked{j},'AlphaData',alpha_mask);
        colormap(colormap_ersp);
        ax = gca;
        %- set figure size and position
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
%         disp(get(ax,'Position'));
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- color lim
        set(gca,'CLim',[-clim_max,clim_max],...
            'xlim',[warping_times(1) warping_times(end)],...
            'ydir','norm',...
            'ylim',[allfreqs(1) SUB_FREQ_LIMS(2)],...
            'yscale','log')
        %- set y-axis label
        if SUB_FREQ_LIMS(2) <= 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end
        %- set x-axis & y-axis labels
        if j == 1
            ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
            xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        else
            xlabel('','FontSize',FONT_SIZE);
            ylabel('','fontsize',FONT_SIZE,'fontweight','bold');
        end
        %- set x-axis ticks
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca;
        ax.XAxis.FontSize = FONT_SIZE;
        title(alltitles{j});  
        horiz_shift = horiz_shift + SHIFT_AMNT;
    end
    hold off;
    %- set color bar
    c = colorbar('XTick', -clim_max:0.2:clim_max);
    c.Position(1) = c.Position(1)+COLORBAR_SHIFT;
%     c.Limits = [-clim_max,clim_max];
%     caxis([-clim_max,clim_max]);
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',FONT_SIZE);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
end

