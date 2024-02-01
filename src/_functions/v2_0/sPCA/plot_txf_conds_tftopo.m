function [fig] = plot_txf_conds_tftopo(allersp,alltimes,allfreqs,...
    varargin)
%## INPUTS
%-
SUBPLOT_BOTTOM = 0.2;
SUBPLOT_INIT_SHIFT = 0.05;
%-
allpcond = [];
%-
alltitles = cell(length(allersp),1);
for i = 1:length(allersp)
    alltitles{i} = sprintf('cond%i',i);
end
PLOT_STRUCT = struct('figure_position_inch',[3,3,14,6],...
    'alltitles',{alltitles},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'clim',[-2,2],...
    'font_size',12,...
    'freq_lims',[],...
    'time_lims',[],...
    'subplot_width',0.15,...
    'subplot_height',0.65,...
    'shift_amnt',0.175,...
    'stats_title','F Statistic Mask',...
    'figure_title','');
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@iscell);
addRequired(p,'alltimes',@isnumeric);
addRequired(p,'allfreqs',@isnumeric);
%## OPTIONAL
addOptional(p,'allpcond',allpcond,@(x) isnumeric(x) || isempty(x) || islogical(x));
%## PARAMETER
addParameter(p,'PLOT_STRUCT',PLOT_STRUCT,@(x) validate_struct(x,PLOT_STRUCT));
parse(p,allersp,alltimes,allfreqs,varargin{:});
%## SET DEFAULTS
allpcond = p.Results.allpcond;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
if isempty(PLOT_STRUCT.freq_lims)
    PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
end
if isempty(PLOT_STRUCT.time_lims)
    PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
end  
%%
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',PLOT_STRUCT.figure_position_inch)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
horiz_shift = 0;
hold on;
for j = 1:length(allersp)
    if ~isempty(allpcond)
        subplot(1,length(allersp)+1,j);
    else
        subplot(1,length(allersp),j);
    end
    tftopo(allersp{j},alltimes,allfreqs,'limits',... 
    [PLOT_STRUCT.time_lims PLOT_STRUCT.freq_lims PLOT_STRUCT.clim],...
    'logfreq','native');
    % contourf(allersp);
    ax = gca;
    hold on;
    colormap(linspecer);
    %-
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[SUBPLOT_INIT_SHIFT+horiz_shift,SUBPLOT_BOTTOM,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]

    %- adjust subplot position and height
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
    %- set ylims
    ylim(log(PLOT_STRUCT.freq_lims))
    if PLOT_STRUCT.freq_lims(2) <= 50
        set(ax,'YTick',log([4.01,8,13,30,50])); 
        set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
    elseif PLOT_STRUCT.freq_lims(2) <= 100
        set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
    end  
    %- set color lims
    set(ax,'clim',PLOT_STRUCT.clim);
    %- set x-axis & y-axis labels
    if j == 1
        ylabel(PLOT_STRUCT.ylabel,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
        xlabel(PLOT_STRUCT.xlabel,'FontSize',PLOT_STRUCT.font_size);
    else
        xlabel('','FontSize',PLOT_STRUCT.font_size);
        ylabel('','fontsize',PLOT_STRUCT.font_size,'fontweight','bold');
    end
    %- set x-axis ticks
    if ~isempty(PLOT_STRUCT.xticklabel_times)
        xrng = get(ax,'XLim');
        if PLOT_STRUCT.xticklabel_times(1) < xrng(1)
            PLOT_STRUCT.xticklabel_times(1) = xrng(1);
        end
        if PLOT_STRUCT.xticklabel_times(end) > xrng(end)
            PLOT_STRUCT.xticklabel_times(end) = xrng(end);
        end
        set(ax,'XTick',PLOT_STRUCT.xticklabel_times);
        if ~isempty(PLOT_STRUCT.xticklabel_chars)
            set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
        end
        for i = 1:length(PLOT_STRUCT.xticklabel_times)
            xline(ax,PLOT_STRUCT.xticklabel_times(i),'k:')
        end
        xtickangle(45);
    end
    title(PLOT_STRUCT.alltitles{j});
    horiz_shift = horiz_shift + PLOT_STRUCT.shift_amnt;
end

%%
%## Add Stats To Plot
if ~isempty(allpcond)
    subplot(1,length(allersp)+1,length(allersp)+1)
    tftopo(allpcond,alltimes,allfreqs,'limits',... 
    [PLOT_STRUCT.time_lims PLOT_STRUCT.freq_lims PLOT_STRUCT.clim],...
    'logfreq','native');
    % contourf(allersp);
    ax = gca;
    hold on;
    colormap(linspecer);
    %-
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
    set(ax,'OuterPosition',[0 0 1 1]);
    set(ax,'Position',[SUBPLOT_INIT_SHIFT+horiz_shift,SUBPLOT_BOTTOM,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]
    
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.04;
    c.Limits = PLOT_STRUCT.clim;
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',PLOT_STRUCT.font_size);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust subplot position and height
    set(ax,'LineWidth',1)
    set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
    %- set ylims
    ylim(log(PLOT_STRUCT.freq_lims))
    if PLOT_STRUCT.freq_lims(2) <= 50
        set(ax,'YTick',log([4.01,8,13,30,50])); 
        set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
    elseif PLOT_STRUCT.freq_lims(2) <= 100
        set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
    end  
    %- set color lims
    set(ax,'clim',PLOT_STRUCT.clim);
    
    %- set x-axis & y-axis labels
    if j == 1
        ylabel(PLOT_STRUCT.ylabel,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
        xlabel(PLOT_STRUCT.xlabel,'FontSize',PLOT_STRUCT.font_size);
    else
        xlabel('','FontSize',PLOT_STRUCT.font_size);
        ylabel('','fontsize',PLOT_STRUCT.font_size,'fontweight','bold');
    end
    %- set x-axis ticks
    if ~isempty(PLOT_STRUCT.xticklabel_times)
        xrng = get(ax,'XLim');
        if PLOT_STRUCT.xticklabel_times(1) < xrng(1)
            PLOT_STRUCT.xticklabel_times(1) = xrng(1);
        end
        if PLOT_STRUCT.xticklabel_times(end) > xrng(end)
            PLOT_STRUCT.xticklabel_times(end) = xrng(end);
        end
        set(ax,'XTick',PLOT_STRUCT.xticklabel_times);
        for i = 1:length(PLOT_STRUCT.xticklabel_times)
            xline(ax,PLOT_STRUCT.xticklabel_times(i),'k:')
        end
        if ~isempty(PLOT_STRUCT.xticklabel_chars)
            set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
            xtickangle(45);
        end
    end
    title(PLOT_STRUCT.stats_title);
else
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = PLOT_STRUCT.clim;
end
if ~isempty(PLOT_STRUCT.figure_title)
    sgtitle(PLOT_STRUCT.figure_title);
end
hold off;
fig = get(groot,'CurrentFigure');
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
            fprintf(2,'\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end