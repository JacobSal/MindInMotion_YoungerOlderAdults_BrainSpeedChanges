function [fig] = plot_txf_conds_tftopo(allersp,alltimes,allfreqs,alltitles,...
    allpcond,clim_max,colormap_ersp,...
    SUB_FREQ_LIMS,FIGURE_POSITION)
    %##
    XTICK_LABEL = 'Gait Events';
    YTICK_LABEL = 'Frequency (Hz)';
    FONT_SIZE = 12;
    SUBPLOT_WIDTH = 0.15;
    SUBPLOT_HEIGHT = 0.7;
    SHIFT_AMNT = 0.175;
    STATS_TITLE = 'CUSTOM STATS';
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
            [alltimes(1) alltimes(end) nan nan clim_ersp],...
            'logfreq','native');
        hold on;
        colormap(colormap_ersp);
        %- adjust subplot position and height
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]

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
%         %- set x-axis ticks
%         xrng = get(ax,'XLim');
%         if warping_times(1) < xrng(1)
%             warping_times(1) = xrng(1);
%         end
%         if warping_times(end) > xrng(end)
%             warping_times(end) = xrng(end);
%         end
%         set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
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
            [alltimes(1) alltimes(end) nan nan  clim_ersp],...
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
%         %- add vertical line
%         for i = 1:length(warping_times)
%             xline(ax,warping_times(i),'k--');
%         end
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
%         xrng = get(ax,'XLim');
%         if warping_times(1) < xrng(1)
%             warping_times(1) = xrng(1);
%         end
%         if warping_times(end) > xrng(end)
%             warping_times(end) = xrng(end);
%         end
%         set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
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