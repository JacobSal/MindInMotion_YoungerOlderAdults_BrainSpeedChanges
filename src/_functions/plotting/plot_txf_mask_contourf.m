function [fig] = plot_txf_mask_contourf(allersp,alltimes,allfreqs,...
    allersp_mask,allpcond,varargin)
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
GROUPTITLE_BOXSIZE = 0.1;
GROUPTITLE_FONTSIZE = 10;
GROUPTITLE_FONTWEIGHT = 'bold';
SUBPLOT_TITLE_FONTSIZE = 8;
SUBPLOT_TITLE_FONTWEIGHT = 'bold';
%-
pcond_chk = @(x) assert(iscell(x) && all((size(allersp)) == size(x)),...
    sprintf('allpcond is type %s and must be type cell, and its size is %s and must be %s',...
    class(x),mk_numarray_str(size(x)),mk_numarray_str(size(allersp)))); 
allersp_mask_chk = @(x) assert(iscell(x) && all((size(allersp)) == size(x)),...
    sprintf('allersp_mask is type %s and must be type cell, and its size is %s and must be %s',...
    class(x),mk_numarray_str(size(x)),mk_numarray_str(size(allersp)))); 

%-
alltitles = cell(length(allersp),1);
for i = 1:length(allersp)
    alltitles{i} = sprintf('cond%i',i);
end
%- 
DEFAULT_PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
    'alltitles',{alltitles},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',[],...
    'xticklabel_chars',{{}},...
    'xticklabel_angle',45,...
    'clim',[],...
    'font_size',8,...
    'font_name','Arial',...
    'freq_lims',[],...
    'time_lims',[],...
    'stats_title','F Stat (p<0.05)',...
    'figure_title','',...
    'contourf_grain',ceil((500/pi())),...
    'alpha_multiple',0.6,...
    'group_titles',{{}},...
    'group_titles_shift_x',0.0,...
    'group_titles_shift_y',0.65,...
    'subplot_width',0.13,...
    'subplot_height',0.16,... %(02/17/2024) was 0.2
    'subplot_shift_x',0.035,...
    'subplot_shift_y',0.05,...
    'subplot_init_y',0.7,...
    'subplot_init_x',0.06,...
    'colorbar_shift_x',0.145,...
    'colorbar_shift_y',0,...
    'colorbar_label_shift_x',0,...
    'colorbar_label_shift_y',0.005,...
    'colorbar_label','\Delta Power (dB)',...
    'colorbar_fontsize',8,...
    'colorbar_fontweight','bold',...
    'display_bandmarks',true,...
    'bandmarks',{{'\theta','\alpha','\beta','\gamma'}},...
    'bandmarks_shift_y',[-0.42,-0.2,0.05,0.325],...
    'bandmarks_shift_x',-0.03,...
    'bandmarks_fontsize',8,...
    'bandmarks_fontweight','bold');
    % 'colorbar_xshift',0.085,... % bigger = shift to right page
    % 'colorbar_yshift',0,... % bigger = shift to top page
    % 'colorbar_label_xshift',0.2,... % bigger = shift to left page (relative to colorbar)
    % 'colorbar_label_yshift',0.025); % bigger = shift to top page (relative to colorbar)
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%% Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@iscell);
addRequired(p,'alltimes',@isnumeric);
addRequired(p,'allfreqs',@isnumeric);
addRequired(p,'allersp_mask',allersp_mask_chk);
addRequired(p,'allpcond',pcond_chk);
%## OPTIONAL
addOptional(p,'allpgroup',@iscell);
%## PARAMETER
addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT));
parse(p,allersp,alltimes,allfreqs,allersp_mask,allpcond,varargin{:});
%## SET DEFAULTS
allpgroup = p.Results.allpgroup;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
if isempty(PLOT_STRUCT.freq_lims)
    PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
end
if isempty(PLOT_STRUCT.time_lims)
    PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
end
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT);
%-

%% COLORLIMITS ALG
%## set color limits
if isempty(PLOT_STRUCT.clim)
    %##
    data = cat(3,allersp{:});
    bound = max([abs(prctile([data],5,'all')),abs(prctile([data],95,'all'))]);
    PLOT_STRUCT.clim = sort([-round(bound,1),round(bound,1)]);
end
%## ASSIGN INTERVALS
colorbar_tick_intervals = round(linspace(PLOT_STRUCT.clim(1),PLOT_STRUCT.clim(2),7),2);
%%
inds1 = allfreqs>=PLOT_STRUCT.freq_lims(1) & allfreqs<=PLOT_STRUCT.freq_lims(2);
inds2 = alltimes>=PLOT_STRUCT.time_lims(1) & alltimes<=PLOT_STRUCT.time_lims(2);
tmp_freqs = allfreqs(inds1);
tmp_times = alltimes(inds2);
if size(allersp{1},1) ~= length(tmp_freqs) || size(allersp{1},2) ~= length(tmp_times)
    allersp = cellfun(@(x) x(inds1,inds2,:),allersp,'uniformoutput',false);
    allersp_mask = cellfun(@(x) x(inds1,inds2,:),allersp_mask,'uniformoutput',false);
    allpcond = cellfun(@(x) x(inds1,inds2,:),allpcond,'uniformoutput',false);
    allfreqs = tmp_freqs;
    alltimes = tmp_times;
end
PLOT_STRUCT.time_lims = [alltimes(1),alltimes(end)];
PLOT_STRUCT.freq_lims = [allfreqs(1),allfreqs(end)];
%%
fig = figure('color','white','renderer','Painters');
set(fig,'Units','inches','Position',PLOT_STRUCT.figure_position_inch)
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
subp_dim1 = size(allersp,2)+double(~isempty(allpgroup));
% supb_cnt = 1;
hold on;
set(gca,AXES_DEFAULT_PROPS{:})
vert_shift = PLOT_STRUCT.subplot_init_y;
for i = 1:size(allersp,2) % group
    horiz_shift = PLOT_STRUCT.subplot_init_x;
    for j = 1:size(allersp,1) % condition
        % subplot(subp_dim1,size(allersp,1),supb_cnt);
        ax = axes();
        %-
        alpha_mask = ones(size(allpcond{j,i}))*PLOT_STRUCT.alpha_multiple; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        alpha_mask(allpcond{j,i} == 1) = 0; %0 is significant? 1 is not?
        contourf(alltimes, allfreqs, allersp{j,i},PLOT_STRUCT.contourf_grain,...
                   'linecolor','none');
        hold on;
        imagesc(alltimes,allfreqs,allersp_mask{j,i},'AlphaData',alpha_mask);
        colormap(linspecer);
        % ax = gca;
        %-
        set(ax,'LineWidth',1,'FontName',PLOT_STRUCT.font_name,...
            'FontSize',PLOT_STRUCT.font_size,...
            'FontWeight','bold',...
            'OuterPosition',[0 0 1 1],...
            'Position',[horiz_shift,vert_shift,PLOT_STRUCT.subplot_width,PLOT_STRUCT.subplot_height]);  %[left bottom width height]
        set(ax,'CLim',PLOT_STRUCT.clim,...
            'xlim',PLOT_STRUCT.time_lims,...
            'ydir','norm',...
            'ylim',PLOT_STRUCT.freq_lims,...
            'yscale','log')
        %- yticks
        if PLOT_STRUCT.freq_lims(2) <= 50
            % freqs = log([4.01,8,13,30,50]);
            freqs = [4.01,8,13,30,50];
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        elseif PLOT_STRUCT.freq_lims(2) <= 100
            % freqs = log([4.01,8,13,30,50,99.4843]);
            freqs = [4.01,8,13,30,50];
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        elseif PLOT_STRUCT.freq_lims(2) <= 200
            freqs = [4.01,8,13,30,50,99.4843,195.5];
            set(ax,'YTick',freqs); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100','200'},'Fontsize',PLOT_STRUCT.font_size);
            for k = 1:length(freqs)
                yline(ax,freqs(k),'k:');
            end
        end  
        %- set color lims
        % set(ax,'clim',PLOT_STRUCT.clim);
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
                xtickangle(PLOT_STRUCT.xticklabel_angle);
            end
        else
            set(ax,'XTickLabel',{});
        end
        title(PLOT_STRUCT.alltitles{j},'FontSize',SUBPLOT_TITLE_FONTSIZE,...
            'FontWeight',SUBPLOT_TITLE_FONTWEIGHT);
        if j == size(allersp,1)
            %- set color bar
            c = colorbar();
            c.Position(1) = horiz_shift+PLOT_STRUCT.colorbar_shift_x;
            c.Position(2) = vert_shift+PLOT_STRUCT.colorbar_shift_y;
            c.Limits = PLOT_STRUCT.clim;
            c.XTick = colorbar_tick_intervals;
            %- color bar label
            hL = annotation(fig,'textbox',...
                [c.Position(1)+PLOT_STRUCT.colorbar_label_shift_x+(0.1/2),...
                c.Position(2)+PLOT_STRUCT.colorbar_label_shift_y+(0.1/2),0.1,0.1],...
                'String',PLOT_STRUCT.colorbar_label,'LineStyle','none',...
                'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.colorbar_fontsize,...
                'FontWeight',PLOT_STRUCT.colorbar_fontweight,...
                'HorizontalAlignment','left','VerticalAlignment','top',...
                'Units','normalized',...
                'Rotation',270);
        end
        horiz_shift = horiz_shift + PLOT_STRUCT.subplot_width + PLOT_STRUCT.subplot_shift_x;
        % supb_cnt = supb_cnt + 1;
    end
    
    if ~isempty(PLOT_STRUCT.group_titles)
        % GROUPTITLE_BOXSIZE = 0.1;
        xx = 0.5+(-GROUPTITLE_BOXSIZE/2)+PLOT_STRUCT.group_titles_shift_x;
        yy = vert_shift+PLOT_STRUCT.subplot_height*PLOT_STRUCT.group_titles_shift_y;
        a = annotation(fig,'textbox',[xx,yy,GROUPTITLE_BOXSIZE,GROUPTITLE_BOXSIZE],...
            'String',PLOT_STRUCT.group_titles{i},'LineStyle','none',...
            'FontName',PLOT_STRUCT.font_name,'FontSize',GROUPTITLE_FONTSIZE,...
            'FontWeight',GROUPTITLE_FONTWEIGHT,...
            'HorizontalAlignment','center','VerticalAlignment','top',...
            'Units','normalized');
    end
    if PLOT_STRUCT.display_bandmarks
        for ii = 1:length(PLOT_STRUCT.bandmarks)
            xx = PLOT_STRUCT.subplot_init_x+PLOT_STRUCT.bandmarks_shift_x;
            yy = vert_shift+PLOT_STRUCT.subplot_height*PLOT_STRUCT.bandmarks_shift_y(ii);
            a = annotation(fig,'textbox',[xx,yy,0.1,0.1],...
                'String',PLOT_STRUCT.bandmarks{ii},'LineStyle','none',...
                'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.bandmarks_fontsize,...
                'FontWeight',PLOT_STRUCT.bandmarks_fontweight,...
                'HorizontalAlignment','left','VerticalAlignment','top','Units','normalized');
        end
    end
    % supb_cnt = supb_cnt + 1;
    vert_shift = vert_shift - (PLOT_STRUCT.subplot_height + PLOT_STRUCT.subplot_shift_y);
end

if ~isempty(PLOT_STRUCT.figure_title)
    sgtitle(PLOT_STRUCT.figure_title,'FontWeight','bold',...
        'FontName',PLOT_STRUCT.font_name,'FontSize',14);
end
hold off;
drawnow;
fig = get(groot,'CurrentFigure');
end