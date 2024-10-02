function [fig] = plot_multipane_txf(allersp,varargin)
alltimes = [];
allfreqs = [];
PLOT_STRUCT = struct('figure_position',[100,100,350,350],...
    'xtick_label','Time (s)',...
    'ytick_label','Frequency (Hz)',...
    'clim',[-2,2],...
    'font_size',12,...
    'freq_lims',[4,60],...
    'time_lims',[0,2],...
    'subplot_width',0.15,...
    'subplot_height',0.7,...
    'subplot_shift',0.175,...
    'colorbar_shift',0.05,...
    'plot_inchs',[3,3,14,5],...
    'colormap',linspecer,...
    'event_times',[0],...
    'event_chars',{{'0'}},...
    'title',['ersp'],...
    'subplot_titles',{{'this','subplot'}},...
    'cbar_intv',0.5,...
    'cbar_label',[{'\Delta Power'};{'(dB)'}]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@iscell);
%## OPTIONAL
addOptional(p,'alltimes',alltimes,@isnumeric);
addOptional(p,'allfreqs',allfreqs,@isnumeric);
%## PARAMETER
addParameter(p,'PLOT_STRUCT',PLOT_STRUCT,@(x) validate_struct(x,PLOT_STRUCT));
parse(p,allersp,varargin{:});
%## SET DEFAULTS
alltimes = p.Results.alltimes;
allfreqs = p.Results.allfreqs;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
%- ASSIGNED VALUES
if isempty(alltimes)
    fprintf('Using default assignment for alltimes. This may not be accurate to the actual times displayed...\n');
    alltimes = 0:size(allersp,1)-1;
end
if isempty(allfreqs)
    fprintf('Using default assignment for allfreqs. This may not be accurate to the actual frequency bands displayed...\n');
    allfreqs = 1:size(allersp,2);
end
if isempty(PLOT_STRUCT.freq_lims)
    PLOT_STRUCT.freq_lims = [min(allfreqs), max(allfreqs)];
end
if isempty(PLOT_STRUCT.time_lims)
    PLOT_STRUCT.time_lims = [min(alltimes), max(alltimes)];
end
%% ===================================================================== %%
%{
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',PLOT_STRUCT.plot_inchs)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
cnt = 1;
for j = 1:size(allersp,1)
    horiz_shift = 0;
    for i = 1:size(allersp,2)
        subplot(size(allersp,1),size(allersp,2),cnt)
%         alpha_mask = ones(size(ersp_pcond{j}))*ALPHA_MULTIPLE; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
%         alpha_mask(ersp_pcond{j} == 1) = 0; %0 is significant? 1 is not?
        contourf(alltimes,allfreqs,allersp{j,i},200,...
                   'linecolor','none','Fill_I','on');
%         imagesc(alltimes,allfreqs,ersp_masked{j},'AlphaData',alpha_mask);
        colormap(PLOT_STRUCT.colormap);
        ax = gca;
        %- set figure size and position
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]
    %         disp(get(ax,'Position'));
        %- add vertical line
        for k = 1:length(PLOT_STRUCT.event_times)
            xline(ax,PLOT_STRUCT.event_times(k),'k--');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if PLOT_STRUCT.event_times(1) < xrng(1)
            PLOT_STRUCT.event_times(1) = xrng(1);
        end
        if PLOT_STRUCT.event_times(end) > xrng(end)
            PLOT_STRUCT.event_times(end) = xrng(end);
        end
        set(ax,'XTick',PLOT_STRUCT.event_times,'XTickLabel',PLOT_STRUCT.event_chars);
        xtickangle(45)
        ax.XAxis.FontSize = PLOT_STRUCT.font_size;
        
        %- set y-axis label
        if PLOT_STRUCT.freq_lims(2) <= 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
        end
        %- set x-axis & y-axis labels
        if i == 1
            ylabel(PLOT_STRUCT.ytick_label,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
            xlabel(PLOT_STRUCT.xtick_label,'FontSize',PLOT_STRUCT.font_size);
        else
            xlabel('','FontSize',PLOT_STRUCT.font_size);
            ylabel('','fontsize',PLOT_STRUCT.font_size,'fontweight','bold');
        end
        if j == 1
            title(PLOT_STRUCT.subplot_titles{i});
        end
        %- color lim
        set(ax,'CLim',PLOT_STRUCT.clim,...
            'xlim',[PLOT_STRUCT.time_lims(1) PLOT_STRUCT.time_lims(end)],...
            'ydir','norm',...
            'ylim',PLOT_STRUCT.freq_lims,...
            'yscale','log')
%         hold off;
        horiz_shift = horiz_shift + PLOT_STRUCT.subplot_shift;
        cnt = cnt + 1;
        hold off;
    end
    if j == 1
        %- set color bar
        c = colorbar('XTick', PLOT_STRUCT.clim(1):PLOT_STRUCT.cbar_intv:PLOT_STRUCT.clim(2));
        c.Position(1) = c.Position(1)+PLOT_STRUCT.colorbar_shift;
        %- color bar label
        hL = ylabel(c,PLOT_STRUCT.cbar_label,'fontweight',...
            'bold','FontName','Arial','FontSize',PLOT_STRUCT.font_size);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        sgtitle(fig,PLOT_STRUCT.title);
        colormap(PLOT_STRUCT.colormap);
    end
end
hold off;
%}
%% ===================================================================== %%
%## TFTOPO
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',PLOT_STRUCT.plot_inchs)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
hold on;
cnt = 1;
for j = 1:size(allersp,1)
    horiz_shift = 0;
    for i = 1:size(allersp,2)
        subplot(size(allersp,1),size(allersp,2),cnt)
%         alpha_mask = ones(size(ersp_pcond{j}))*ALPHA_MULTIPLE; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
%         alpha_mask(ersp_pcond{j} == 1) = 0; %0 is significant? 1 is not?
        tftopo(double(allersp{j,i}),alltimes,allfreqs,'limits',... 
            [alltimes(1) alltimes(end) nan nan PLOT_STRUCT.clim],...
            'logfreq','native');
        hold on;
%         imagesc(alltimes,allfreqs,ersp_masked{j},'AlphaData',alpha_mask);
        colormap(PLOT_STRUCT.colormap);
        ax = gca;
        %- set figure size and position
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]
    %         disp(get(ax,'Position'));
        %- add vertical line
        for k = 1:length(PLOT_STRUCT.event_times)
            xline(ax,PLOT_STRUCT.event_times(k),'k--');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if PLOT_STRUCT.event_times(1) < xrng(1)
            PLOT_STRUCT.event_times(1) = xrng(1);
        end
        if PLOT_STRUCT.event_times(end) > xrng(end)
            PLOT_STRUCT.event_times(end) = xrng(end);
        end
        set(ax,'XTick',PLOT_STRUCT.event_times,'XTickLabel',PLOT_STRUCT.event_chars);
        xtickangle(45)
        ax.XAxis.FontSize = PLOT_STRUCT.font_size;
        
        %- set y-axis label
        if PLOT_STRUCT.freq_lims(2) <= 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
        end
        %- set x-axis & y-axis labels
        if i == 1
            ylabel(PLOT_STRUCT.ytick_label,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
            xlabel(PLOT_STRUCT.xtick_label,'FontSize',PLOT_STRUCT.font_size);
        else
            xlabel('','FontSize',PLOT_STRUCT.font_size);
            ylabel('','fontsize',PLOT_STRUCT.font_size,'fontweight','bold');
        end
        if j == 1
            title(PLOT_STRUCT.subplot_titles{i});
        end
        %- color lim
        set(ax,'CLim',PLOT_STRUCT.clim,...
            'xlim',[PLOT_STRUCT.time_lims(1) PLOT_STRUCT.time_lims(end)],...
            'ydir','norm',...
            'ylim',PLOT_STRUCT.freq_lims,...
            'yscale','log')
%         hold off;
        horiz_shift = horiz_shift + PLOT_STRUCT.subplot_shift;
        cnt = cnt + 1;
        hold off;
    end
    if j == 1
        %- set color bar
        c = colorbar('XTick', PLOT_STRUCT.clim(1):PLOT_STRUCT.cbar_intv:PLOT_STRUCT.clim(2));
        c.Position(1) = c.Position(1)+PLOT_STRUCT.colorbar_shift;
        %- color bar label
        hL = ylabel(c,PLOT_STRUCT.cbar_label,'fontweight',...
            'bold','FontName','Arial','FontSize',PLOT_STRUCT.font_size);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        sgtitle(fig,PLOT_STRUCT.title);
        colormap(PLOT_STRUCT.colormap);
    end
end
hold off;
end
%% ===================================================================== %%
function [b] = validate_struct(x,DEFAULT_STRUCT)
    b = false;
    struct_name = inputname(2);
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field names
    chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
    if ~all(chk)
        fprintf(2,'\nFields for struct do not match for %s\n',struct_name);
        return
    end
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nValue must be type %s, but is type %s\n',class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end