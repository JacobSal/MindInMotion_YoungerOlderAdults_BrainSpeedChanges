function [fig] = plot_txf_conds_tftopo(allersp,alltimes,allfreqs,...
    varargin)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu
cat_logo();
%-
SUBPLOT_BOTTOM = 0.70;
SUBPLOT_INIT_SHIFT = 0.06;
% COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2];
% COLOR_LIM_ERR = 0.05;
COLORBAR_SHIFT = 0.07; %(02/17/2024) was 0.06
%-
allpcond = [];
allpgroup = [];
%-
alltitles = cell(length(allersp),1);
for i = 1:length(allersp)
    alltitles{i} = sprintf('cond%i',i);
end
DEFAULT_PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
    'alltitles',{alltitles},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[],...
    'time_lims',[],...
    'subplot_width',0.13,...
    'subplot_height',0.16,... %(02/17/2024) was 0.2
    'horiz_shift_amnt',0.17,...
    'vert_shift_amnt',0.22,...
    'group_titles',{{}},...
    'stats_title','F Statistic Mask',...
    'figure_title','');
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@iscell);
addRequired(p,'alltimes',@isnumeric);
addRequired(p,'allfreqs',@isnumeric);
%## OPTIONAL
addOptional(p,'allpcond',allpcond,@(x) iscell(x) || isempty(x) || islogical(x));
addOptional(p,'allpgroup',allpgroup,@(x) iscell(x) || isempty(x) || islogical(x));
%## PARAMETER
addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT));
parse(p,allersp,alltimes,allfreqs,varargin{:});
%## SET DEFAULTS
allpcond = p.Results.allpcond;
allpgroup = p.Results.allpgroup;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
if isempty(PLOT_STRUCT.freq_lims)
    PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
end
if isempty(PLOT_STRUCT.time_lims)
    PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
end
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT);
%% COLORLIMITS ALG
%## set color limits
if isempty(PLOT_STRUCT.clim)
    %##
    data = cat(3,allersp{:});
    bound = max([abs(prctile([data],5,'all')),abs(prctile([data],95,'all'))]);
    PLOT_STRUCT.clim = sort([-round(bound,1),round(bound,1)]);
end
%## GET INTERVALS
INTERVALS = round(linspace(PLOT_STRUCT.clim(1),PLOT_STRUCT.clim(2),7),2);
%%
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',PLOT_STRUCT.figure_position_inch)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
vert_shift = 0;
subp_dim1 = size(allersp,2)+double(~isempty(allpgroup));
supb_cnt = 1;
hold on;
for i = 1:size(allersp,2)
    horiz_shift = 0;
    for j = 1:size(allersp,1)
        if ~isempty(allpcond)
            subplot(subp_dim1,size(allersp,1)+double(~isempty(allpcond{1,i})),supb_cnt);
        else
            subplot(subp_dim1,size(allersp,1),supb_cnt);
        end
        tftopo(allersp{j,i},alltimes,allfreqs,'limits',... 
        [PLOT_STRUCT.time_lims PLOT_STRUCT.freq_lims PLOT_STRUCT.clim],...
        'logfreq','native');
        % contourf(allersp);
        ax = gca;
        hold on;
        colormap(linspecer);
        %-
        set(ax,'LineWidth',1)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[SUBPLOT_INIT_SHIFT+horiz_shift,SUBPLOT_BOTTOM+vert_shift,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]

        %- adjust subplot position and height
        set(ax,'LineWidth',1)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        %- set ylims
        ylim(log(PLOT_STRUCT.freq_lims))
        if PLOT_STRUCT.freq_lims(2) <= 50
            freqs = log([4.01,8,13,30,50]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            freqs = log([4.01,8,13,30,50,99.4843]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        end  
        %- set color lims
        set(ax,'clim',PLOT_STRUCT.clim);
        %- set x-axis & y-axis labels
        if j == 1 && i == subp_dim1
            xlabel(PLOT_STRUCT.xlabel,'FontSize',PLOT_STRUCT.font_size);
            ylabel(PLOT_STRUCT.ylabel,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
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
            for k = 1:length(PLOT_STRUCT.xticklabel_times)
                xline(ax,PLOT_STRUCT.xticklabel_times(k),'k--')
            end
        end
        if i == subp_dim1
            if ~isempty(PLOT_STRUCT.xticklabel_chars)
                set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
                xtickangle(45);
            end
        else
            set(ax,'XTickLabel',{});
        end
        title(PLOT_STRUCT.alltitles{j},'FontSize',8,'FontWeight','bold');
        horiz_shift = horiz_shift + PLOT_STRUCT.horiz_shift_amnt;
        supb_cnt = supb_cnt + 1;
    end
    %## Add Stats To Plot
    if ~isempty(allpcond)
        subplot(subp_dim1,length(allersp)+double(~isempty(allpcond{1,i})),supb_cnt)
        tftopo(allpcond{1,i}*1000,alltimes,allfreqs,'limits',... 
        [PLOT_STRUCT.time_lims PLOT_STRUCT.freq_lims PLOT_STRUCT.clim],...
        'logfreq','native');
        % (02/05/2024) adding *1000 to ensure high visibility of stats
        % contourf(allersp);
        ax = gca;
        hold on;
        colormap(linspecer);
        %-
        set(ax,'LineWidth',1)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[SUBPLOT_INIT_SHIFT+horiz_shift,SUBPLOT_BOTTOM+vert_shift,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]

        %- adjust subplot position and height
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        %- set ylims
        ylim(log(PLOT_STRUCT.freq_lims))
        if PLOT_STRUCT.freq_lims(2) <= 50
            freqs = log([4.01,8,13,30,50]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            freqs = log([4.01,8,13,30,50,99.4843]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        end  
        %- set color lims
        set(ax,'clim',PLOT_STRUCT.clim);
        %- set x-axis & y-axis labels
        if j == 1 && i == subp_dim1
            xlabel(PLOT_STRUCT.xlabel,'FontSize',PLOT_STRUCT.font_size);
            ylabel(PLOT_STRUCT.ylabel,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
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
            for k = 1:length(PLOT_STRUCT.xticklabel_times)
                xline(ax,PLOT_STRUCT.xticklabel_times(k),'k--')
            end
        end
        if i == subp_dim1
            if ~isempty(PLOT_STRUCT.xticklabel_chars)
                set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
                xtickangle(45);
            end
        else
            set(ax,'XTickLabel',{});
        end
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+COLORBAR_SHIFT;
        c.Limits = PLOT_STRUCT.clim;
        c.XTick = INTERVALS;
        %- color bar label
    %     hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
    %         'bold','FontName',PLOT_STRUCT.font_name,'FontSize',8);
        hL = ylabel(c,[{'\Delta Power (dB)'}],'fontweight',...
            'bold','FontName',PLOT_STRUCT.font_name,'FontSize',8);
    %     set(hL,'Rotation',0);
        set(hL,'Rotation',270);
        hL.Units = 'Inches';
        hL.Position = [SUBPLOT_INIT_SHIFT+horiz_shift-0.175,SUBPLOT_BOTTOM+0.0925,0];
        title(PLOT_STRUCT.stats_title,'FontSize',8,'FontWeight','bold');
        supb_cnt = supb_cnt+1;
    else
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+COLORBAR_SHIFT;
        c.Limits = PLOT_STRUCT.clim;
        c.XTick = INTERVALS;
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName',PLOT_STRUCT.font_name,'FontSize',8);
        set(hL,'Rotation',270);
        hL.Units = 'Inches';
        hL.Position = [SUBPLOT_INIT_SHIFT+horiz_shift-0.175,SUBPLOT_BOTTOM+0.0925,0];
    end
    if ~isempty(PLOT_STRUCT.group_titles)
        xx = PLOT_STRUCT.figure_position_inch(3)/2/PLOT_STRUCT.figure_position_inch(3);
        yy = PLOT_STRUCT.subplot_height+SUBPLOT_BOTTOM+vert_shift;
        a = annotation(fig,'textbox',[xx-0.05,yy-0.05,0.1,0.1],...
            'String',PLOT_STRUCT.group_titles{i},'LineStyle','none',...
            'FontName',PLOT_STRUCT.font_name,'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end
    vert_shift = vert_shift - PLOT_STRUCT.vert_shift_amnt;
end
if ~isempty(allpgroup)
    i = i + 1;
    horiz_shift = 0;
    for j = 1:size(allpgroup,2)
        if ~isempty(allpcond)
            subplot(subp_dim1,size(allersp,1)+double(~isempty(allpcond{1,1})),supb_cnt);
        else
            subplot(subp_dim1,size(allersp,1),supb_cnt);
        end
        tftopo(allpgroup{1,j},alltimes,allfreqs,'limits',... 
        [PLOT_STRUCT.time_lims PLOT_STRUCT.freq_lims PLOT_STRUCT.clim],...
        'logfreq','native');
        % contourf(allersp);
        ax = gca;
        hold on;
        colormap(linspecer);
        %-
        set(ax,'LineWidth',1)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[SUBPLOT_INIT_SHIFT+horiz_shift,SUBPLOT_BOTTOM+vert_shift,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]

        %- adjust subplot position and height
        set(ax,'LineWidth',1)
        set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size,'FontWeight','bold')
        %- set ylims
        ylim(log(PLOT_STRUCT.freq_lims))
        if PLOT_STRUCT.freq_lims(2) <= 50
            freqs = log([4.01,8,13,30,50]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            freqs = log([4.01,8,13,30,50,99.4843]);
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        end  
        %- set color lims
        set(ax,'clim',PLOT_STRUCT.clim);
        %- set x-axis & y-axis labels
        if j == 1 && i == subp_dim1
            xlabel(PLOT_STRUCT.xlabel,'FontSize',PLOT_STRUCT.font_size);
            ylabel(PLOT_STRUCT.ylabel,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
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
            for k = 1:length(PLOT_STRUCT.xticklabel_times)
                xline(ax,PLOT_STRUCT.xticklabel_times(k),'k--')
            end
        end
        if i == subp_dim1
            if ~isempty(PLOT_STRUCT.xticklabel_chars)
                set(ax,'XTickLabel',PLOT_STRUCT.xticklabel_chars);
                xtickangle(45);
            end
        else
            set(ax,'XTickLabel',{});
        end
        title(PLOT_STRUCT.alltitles{j},'FontSize',8,'FontWeight','bold');
        horiz_shift = horiz_shift + PLOT_STRUCT.horiz_shift_amnt;
        supb_cnt = supb_cnt + 1;
    end
    if ~isempty(PLOT_STRUCT.group_titles)
        xx = PLOT_STRUCT.figure_position_inch(3)/2/PLOT_STRUCT.figure_position_inch(3);
        yy = PLOT_STRUCT.subplot_height+SUBPLOT_BOTTOM+vert_shift;
        a = annotation(fig,'textbox',[xx-0.05,yy-0.05,0.1,0.1],...
            'String',PLOT_STRUCT.group_titles{i},'LineStyle','none',...
            'FontName',PLOT_STRUCT.font_name,'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end
end
if ~isempty(PLOT_STRUCT.figure_title)
    sgtitle(PLOT_STRUCT.figure_title);
end
hold off;
drawnow;
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
%##
function [struct_out] = set_defaults_struct(x,DEFAULT_STRUCT)
    struct_out = x;
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        if isempty(vals1{ind})
            struct_out.(fs1{ind}) = DEFAULT_STRUCT.(fs2{ind});
        end
    end
end