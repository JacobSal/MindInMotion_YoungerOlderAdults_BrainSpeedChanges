function fooof_scatter_plot(x,y,color,mode,plot_nan)
    if exist('mode') == 2
        mode = '';
    end
    if exist('plot_nan') == 0
        plot_nan = 0;
    end
    if ~isempty(y)
        try
            scatter(x,y(:,2),[],'MarkerFaceColor',color,'MarkerEdgeColor','black');
        catch
            scatter(x,y(:,1),[],'MarkerFaceColor',color,'MarkerEdgeColor','black');
        end
    end
    if plot_nan & (isempty(y))
        scatter(x,zeros(length(x)),[],'MarkerFaceColor','white','MarkerEdgeColor','black');
    end
    switch mode
        case 'connect'
            if ~isempty(y)
                plot(x,y(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
            end
    end
end