function [ax] = group_violin(table_in,measure_char,cond_char,group_char,varargin)
%GROUP_VIOLIN Summary of this function goes here
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
%- FONTS
ax = [];
% GROUP_LAB_FONTSIZE = 10;
% GROUP_LAB_FONTWEIGHT = 'normal';
% XLAB_FONTSIZE = 10;
% YLAB_FONTSIZE = 10;
% XTICK_FONTSIZE = 10;
% XLAB_FONTWEIGHT = 'bold';
% YLAB_FONTWEIGHT = 'bold';
% TITLE_FONTSIZE = 10;
% TITLE_FONTWEIGHT = 'bold';
%-
% AXES_POS = [0.15 0.20 0.8 0.7];
% XLABEL_OFFSET = -0.225;
% GROUPLAB_POS = [0.225,0.775];
% GROUPLAB_POS = [0.125,0.475,0.812];

% GROUPLAB_YOFFSET = -0.285;
% REGRESS_TXT_XMULTI = 1; %-0.8;
% REGRESS_TXT_YMULTI = 1;
% REGRESS_TXT_SIZE = 6;
% REGRESS_TXT_SIZE = 10;
SCATTER_MARK_SIZE = 5;
STD_CUTOFF = 3;
DEFAULT_VIOLIN_PARAMS = {'Width',0.05,...
            'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
            'ShowMedian',true,'Bandwidth',0.1,'QuartileStyle','shadow',...
            'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
            'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
prc_ylim = [floor(prctile(table_in.(measure_char),1))-floor(std(table_in.(measure_char)))*1.5,...
            ceil(prctile(table_in.(measure_char),99))+ceil(std(table_in.(measure_char)))*1.5];
DEFAULT_PLOT_STRUCT = struct('color_map',linspecer(length(unique(table_in.(cond_char)))),...
    'cond_labels',{{string(unique(table_in.(cond_char)))}},...
    'cond_offsets',linspace(-0.3,0.3,length(unique(table_in.(cond_char)))),...
    'group_labels',{{string(unique(table_in.(group_char)))}},...
    'group_offsets',[0.125,0.475,0.812],...
    'group_lab_yoffset',-0.285,...
    'group_lab_fontweight','normal',...
    'group_lab_fontsize',12,...
    'y_label','',...
    'y_label_fontsize',12,...
    'y_label_fontweight','bold',...
    'ylim',prc_ylim,...
    'x_label','',...
    'x_label_fontsize',12,...
    'x_label_fontweight','bold',...
    'x_label_yoffset',-0.1,...
    'xlim',[],...
    'title',{{''}},...
    'title_fontsize',12,...
    'title_fontweight','normal',...
    'font_size',12,...
    'font_name','Arial',...
    'do_combine_groups',false,...
    'regresslab_txt_size',5,...
    'ax_position',[0,0,1,1],...
    'ax_line_width',1);

% STATS_STRUCT = struct('anova',{{}},...
%                       'pvals',{{}},...
%                       'pvals_pairs',{{}},...
%                       'regress_pval',{{}},...
%                       'regress_line',{{}},...
%                       'r2_coeff',0);
DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'line_type',{'best_fit'},... % ('best_fit' | 'means')
                          'regress_xvals',0,...
                          'subject_char',[],... % this option when filled prints removal of nan() info
                          'group_order',categorical({''}),...
                          'display_stats_char',false,... 
                          'stats_char',{{}},...
                          'bracket_conn_yshift',repmat(1.1,[length(unique(table_in.(group_char))),1]),...
                          'bracket_rawshifty_upper',0,...
                          'bracket_rawshifty_lower',0,...
                          'grp_sig_offset_x',zeros(length(unique(table_in.(group_char))),1),...
                          'grp_sig_offset_y',zeros(length(unique(table_in.(group_char))),1)); 

%- STATS_STRUCT EXAMPLE:
% STATS_STRUCT = struct('anova',{0.02,0.1},...
%                       'pvals',{[0.001,0.002,0.01],[0.2,0.3,0.05,0.1]},...
%                       'pvals_pairs,{[1,2],[1,3],[1,4]},...
%                       'regress_pval',{{}},...
%                       'regress_line',{{}},...
%                       'r2_ceoff',0);
%## CHECK FUNCTIONS
isnumscalar = @(x) (isnumeric(x)&isscalar(x));
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'table_in',@istable);
addRequired(p,'measure_char',@ischar);
addRequired(p,'cond_char',@ischar);
addRequired(p,'group_char',@ischar);
%## OPTIONAL
addOptional(p,'parent_axes',[],@(x) isobject(x) || isempty(x))
%## PARAMETER
addParameter(p,'VIOLIN_PARAMS',DEFAULT_VIOLIN_PARAMS,@iscell);
addParameter(p,'PLOT_STRUCT',DEFAULT_PLOT_STRUCT,@(x) validate_struct(x,DEFAULT_PLOT_STRUCT));
addParameter(p,'STATS_STRUCT',DEFAULT_STATS_STRUCT,@(x) validate_struct(x,DEFAULT_STATS_STRUCT));
%##
parse(p,table_in,measure_char,cond_char,group_char,varargin{:});
%## SET DEFAULTS
ax = p.Results.parent_axes;
VIOLIN_PARAMS = p.Results.VIOLIN_PARAMS;
PLOT_STRUCT = p.Results.PLOT_STRUCT;
STATS_STRUCT = p.Results.STATS_STRUCT;
% VIOLIN_PARAMS = set_defaults_struct(VIOLIN_PARAMS,DEFAULT_VIOLIN_PARAMS);
PLOT_STRUCT = set_defaults_struct(PLOT_STRUCT,DEFAULT_PLOT_STRUCT);
STATS_STRUCT = set_defaults_struct(STATS_STRUCT,DEFAULT_STATS_STRUCT);
%% ===================================================================== %%
%##
table_tmp = table_in;
%##
%- set groups
if any(isempty(STATS_STRUCT.group_order)) || any(isundefined(STATS_STRUCT.group_order))
    groups = unique(table_tmp.(group_char));
else
    warning('STATS_STRUCT.group_order is being used, make sure to change PLOT_STRUCT.group_labels to match...')
    groups = STATS_STRUCT.group_order;
end
if PLOT_STRUCT.do_combine_groups
    table_tmp.(group_char) = categorical(ones(height(table_tmp),1));
    PLOT_STRUCT.group_offsets = 0.5;
    groups = unique(table_tmp.(group_char));
    PLOT_STRUCT.group_labels = categorical({[sprintf('%s',PLOT_STRUCT.group_labels(1)), sprintf(' & %s',PLOT_STRUCT.group_labels(2:end))]});
else
    fprintf('Plotting groups seperately...\n')
end
%- set conditions
conds = unique(table_tmp.(cond_char));
%%
% if isempty(parent_axes)
%     ax = axes(figure);
% else
%     ax = axes(parent_axes);
% end
hold on;
cnt = 1;
g_o = 1; 
g_cnt = 1;
xticks = [];
xtick_labs = {};
gticks = [];
%- group primary, cond secondary
ymaxs = zeros(length(groups),length(conds));
sigline_ymax = zeros(length(groups),length(conds));
violins = cell(length(groups)*length(conds),1);
if ~isempty(PLOT_STRUCT.group_offsets)
    g_offset = PLOT_STRUCT.group_offsets(g_o);
    g_o = g_o + 1;
else
    g_offset = 0;
end
for i=1:length(groups)
    for j=1:length(conds)
        g_i = groups(i);
        c_i = conds(j);
        inds = table_tmp.(cond_char)==c_i & table_tmp.(group_char)==g_i;
        %- outliers
        tt = table_tmp(inds,:);
        vals = tt.(measure_char);
        sigline_ymax(i,j) = max(tt.(measure_char));
        chk = isnan(vals);
        %- remove NaNs from vals
        if any(chk)
            if ~isempty(STATS_STRUCT.subject_char)
                fprintf('\nRemoving NaN''s from data for subjects: \n');
                fprintf('%s, ',tt.(STATS_STRUCT.subject_char)(logical(chk)))
            end
            vals = vals(~chk);
        end
        %- calculate mean and std and crop out data beyond limits
        data_mean = mean(vals);
        data_std = std(vals);
        ind_crop = tt.(measure_char)>data_mean-STD_CUTOFF*data_std & tt.(measure_char)<data_mean+STD_CUTOFF*data_std;
        outliers = tt(~ind_crop,:);
        tt = tt(ind_crop,:);
        %- final table in
        offset = PLOT_STRUCT.cond_offsets(j)+g_offset;
        %- check input types
        if ~isnumscalar(g_i)
            tmp_gi = i; %double(g_i); %double(string(g_i));
        else
            tmp_gi = g_i;
        end
        %## MAIN MEASURES
        if ~isempty(tt)
            ymaxs(i,j) = max(tt.(measure_char));
            %- plot main violin for each condition & group
            bandwidth = range(vals)*0.1;
            ind = find(strcmp(VIOLIN_PARAMS,'Bandwidth'));
            if isempty(ind)
                VIOLIN_PARAMS = [VIOLIN_PARAMS, 'Bandwidth', bandwidth];
            else
                VIOLIN_PARAMS{ind+1} = bandwidth;
            end
            % violins(cnt) = Violin({table_in.(measure_char)},tmp_gi,VIOLIN_PARAMS{:});
            violins{cnt} = Violin({tt.(measure_char)},tmp_gi,VIOLIN_PARAMS{:});
            violins{cnt}.ScatterPlot.XData      = violins{cnt}.ScatterPlot.XData+offset;
            violins{cnt}.ViolinPlot.XData       = violins{cnt}.ViolinPlot.XData+offset;
            violins{cnt}.WhiskerPlot.XData      = violins{cnt}.WhiskerPlot.XData+offset;
            violins{cnt}.MedianPlot.XData       = violins{cnt}.MedianPlot.XData+offset;
            violins{cnt}.NotchPlots(1).XData    = violins{cnt}.NotchPlots(1).XData+offset;
            violins{cnt}.NotchPlots(2).XData    = violins{cnt}.NotchPlots(2).XData+offset;
            violins{cnt}.MeanPlot.XData         = violins{cnt}.MeanPlot.XData+offset;
            violins{cnt}.ViolinPlotQ.XData      = violins{cnt}.ViolinPlotQ.XData+offset;
            violins{cnt}.ViolinColor            = {PLOT_STRUCT.color_map(j,:)};
            xticks = [xticks,tmp_gi+offset];
            if iscell(PLOT_STRUCT.cond_labels)
                xtick_labs = [xtick_labs,{sprintf('%s',PLOT_STRUCT.cond_labels{j})}];
            else
                xtick_labs = [xtick_labs,{string(PLOT_STRUCT.cond_labels(j))}];
            end
        else
            fprintf('Condition %s & Group %s does not have entries\n',string(c_i),string(g_i))
        end
        %## OUTLIERS
        if ~isempty(outliers)
            g_ind = (outliers.(cond_char)==c_i) & (outliers.(group_char)==g_i);
            vals_y = outliers.(measure_char)(g_ind);
            if ~isempty(vals_y)
                in_x = repmat(xticks(cnt),size(vals_y));
                scatter(in_x,vals_y,SCATTER_MARK_SIZE,'k*','jitter','on', 'jitterAmount', 0.05)
            end
        else
            fprintf('Condition %s & Group %s does not have outliers\n',string(c_i),string(g_i))
        end
        cnt = cnt+1;
    end
    gticks(i) = mean(xticks((cnt-length(conds)):(cnt-1)));
    if ~isempty(PLOT_STRUCT.group_offsets) & i < length(groups)
        g_offset = PLOT_STRUCT.group_offsets(g_o);
        g_o = g_o + 1;
    else
        g_offset = g_offset + 0.2;
    end
end
hold on;
%##
% SIG_LINE_INCR = 0.1;
% ymax = max(table_tmp.(measure_char));
% yiter = get(gca,'ytick');
% yiter = yiter(2)-yiter(1);
cnt_g = 1;
hold_xlim = get(gca,'xlim');
annotes = [];
set_y = true;
% bracket_conn_y_init = STATS_STRUCT.bracket_y_perc;
bracket_rawshifty_upper = STATS_STRUCT.bracket_rawshifty_upper;
bracket_rawshifty_lower = STATS_STRUCT.bracket_rawshifty_lower;
% bracket_conn_yshift_perc = STATS_STRUCT.bracket_yshift_perc;
% bracket_conn_y_perc = STATS_STRUCT.bracket_y_perc;
for i=1:length(groups)
    %- GROUP STATISTICS
    if ~isempty(STATS_STRUCT.anova_grp)
        if STATS_STRUCT.anova_grp{i} < 0.05
            if STATS_STRUCT.pvals_grp{i} < 0.05
                x = zeros(1,2);
                bx1 = zeros(1,2);
                bx2 = zeros(1,2);
                by1 = zeros(1,2);
                by2 = zeros(1,2);
                g1 = STATS_STRUCT.pvals_grp_pairs{i}(1,1);
                g2 = STATS_STRUCT.pvals_grp_pairs{i}(1,2);
                % if i ~= 1
                    if g1 ~= 1
                        bx1(1) = (g1-1)*length(conds)+1;
                        bx1(2) = (g1-1)*length(conds)+length(conds);
                    else
                        bx1(1) = g1;
                        bx1(2) = g1+length(conds)-1;
                    end
                    if g2 ~= 1
                        bx2(1) = (g2-1)*length(conds)+1;
                        bx2(2) = (g2-1)*length(conds)+length(conds);
                    else
                        bx2(1) = g2;
                        bx2(2) = g2+length(conds)-1;
                    end
                    tmpy = [];
                    for tt = linspace(bx1(1),bx1(2),bx1(2)-bx1(1)+1)
                        tmpy = [tmpy, violins{tt}.ScatterPlot.YData];
                    end
                    by1(1) = max(tmpy)*1.025; 
                    by1(2) = max(tmpy)*1.05;
                    tmpy = [];
                    for tt = linspace(bx2(1),bx2(2),bx2(2)-bx2(1)+1)
                        tmpy = [tmpy, violins{tt}.ScatterPlot.YData];
                    end
                    by2(1) = max(tmpy)*1.025; 
                    by2(2) = max(tmpy)*1.05;
                    %-
                    for tt = 1:length(bx1)
                        bx1(tt) = violins{bx1(tt)}.MedianPlot.XData;
                        bx2(tt) = violins{bx2(tt)}.MedianPlot.XData;
                    end
                    %## USER ADJUSTMENTS
                    by1(1,2) = by1(1,2) + bracket_rawshifty_upper;
                    by2(1,2) = by2(1,2) + bracket_rawshifty_upper;
                    by1(1,1) = by1(1,1) + bracket_rawshifty_lower;
                    by2(1,1) = by2(1,1) + bracket_rawshifty_lower;
                    if length(STATS_STRUCT.grp_sig_offset_x) > 1 
                        sig_offset_x = STATS_STRUCT.grp_sig_offset_x(i);
                    else
                        sig_offset_x = STATS_STRUCT.grp_sig_offset_x;
                    end
                    if length(STATS_STRUCT.grp_sig_offset_y) > 1
                        sig_offset_y = STATS_STRUCT.grp_sig_offset_y(i);
                    else
                        sig_offset_y = STATS_STRUCT.grp_sig_offset_y;
                    end
                    if length(STATS_STRUCT.bracket_conn_yshift) > 1
                        bracket_conn_yshift = STATS_STRUCT.bracket_conn_yshift(i);
                    else
                        bracket_conn_yshift = STATS_STRUCT.bracket_conn_yshift;
                    end
                    %## PLOT BRACKET
                    pp = cus_sigbracket('+',STATS_STRUCT.pvals_grp{i},bx1,bx2,by1,by2,...
                        bracket_conn_yshift,sig_offset_x,sig_offset_y);
                    y = gety(pp);
                    set_y = false;
                    % bracket_conn_y_init = bracket_conn_y_init*(1+bracket_conn_yshift_perc)^2;
                % end
            end
        end
    end
    if set_y
        y = max(sigline_ymax(i,:));
    end
    %- CONDITION STATISTICS
    for j = 1:length(conds)
        if ~isempty(STATS_STRUCT.anova) && ~isempty(STATS_STRUCT.pvals)
            if STATS_STRUCT.anova{i} < 0.05
                if STATS_STRUCT.pvals{i}(j) < 0.05
                    % x = get(ax,'XTick');
                    x = zeros(1,2);
                    if i == 1
                        pair = STATS_STRUCT.pvals_pairs{i}{j};
                        % x = x(pair+1)-(x(2)-x(1))/2;
                        for tt = 1:length(pair)
                            x(tt) = violins{pair(tt)}.MedianPlot.XData;
                        end
                    else
                        pair = STATS_STRUCT.pvals_pairs{i}{j}+(i-1)*length(conds);
                        % x = x(pair+1)+(x(2)-x(1))/2;
                        for tt = 1:length(pair)
                            x(tt) = violins{pair(tt)}.MedianPlot.XData;
                        end
                    end
                    pp = cus_sigline(x,'*',[],y,STATS_STRUCT.pvals{i}(j));
                    % gety(pp)
                    y = gety(pp);
                    % y = y + yiter*0.4; %gety(gca); %SIG_LINE_INCR;
                end
            end
            
        end
    end
    %- REGRESSION-CONTINUOUS
    if ~all(isempty(STATS_STRUCT.regress_pval)) && ~all(isempty(STATS_STRUCT.regress_line)) %&& ~all(STATS_STRUCT.r2_coeff==0)
        if STATS_STRUCT.anova{i} < 0.05
            %- plot line of best fit
            if strcmp(STATS_STRUCT.line_type,'best_fit')
                x = xticks(cnt_g:cnt_g+length(conds)-1);
                incr = x(2)-x(1);
                x = [x(1)-incr, x];
                x = [x, x(end)+incr];
                y = STATS_STRUCT.regress_xvals*STATS_STRUCT.regress_line{i}(2) + STATS_STRUCT.regress_line{i}(1);
                plot(x,y,'-','color','k','linewidth',1);
            %-
            elseif strcmp(STATS_STRUCT.line_type,'means')
                x = xticks(cnt_g:cnt_g+length(conds)-1);
                y = STATS_STRUCT.regress_line{i}(:);
                plot(x,y,'-','color','k','linewidth',3);
            end
            x_txt = PLOT_STRUCT.group_offsets(i)-0.1;
            y_txt = 0.9;
            if STATS_STRUCT.display_stats_char
                if length(STATS_STRUCT.stats_char) ~= length(groups)
                    error('ERROR. STATS_STRUCT.stats_char must be length(unique(table.group_char)).')
                end
                str = STATS_STRUCT.stats_char{i};
                text(x_txt,y_txt,str,...
                        'FontSize',PLOT_STRUCT.regresslab_txt_size,...
                        'FontName',PLOT_STRUCT.font_name,...
                        'FontWeight','bold','Units','normalized');
            end
            % if STATS_STRUCT.display_one_regress && i == 1
            %     if STATS_STRUCT.do_include_intercept
            %         str = sprintf('b=%0.2g\nm=%0.2g\nR^2=%0.2g',STATS_STRUCT.regress_line{i}(1),STATS_STRUCT.regress_line{i}(2),STATS_STRUCT.r2_coeff(i));
            %     else
            %         str = sprintf('m=%0.2g\nR^2=%0.2g',STATS_STRUCT.regress_line{i}(2),STATS_STRUCT.r2_coeff(i));
            %     end
            %     if STATS_STRUCT.regress_pval{i} > 0.01 && STATS_STRUCT.regress_pval{i} < 0.05
            %         text(x_txt,y_txt,['* ',str],...
            %             'FontSize',PLOT_STRUCT.regresslab_txt_size,...
            %             'FontName',PLOT_STRUCT.font_name,...
            %             'FontWeight','bold','Units','normalized');
            %     elseif STATS_STRUCT.regress_pval{i} <= 0.01 && STATS_STRUCT.regress_pval{i} > 0.001
            %         text(x_txt,y_txt,['** ',str],...
            %             'FontSize',PLOT_STRUCT.regresslab_txt_size,...
            %             'FontName',PLOT_STRUCT.font_name,...
            %             'FontWeight','bold','Units','normalized');
            %     else
            %         text(x_txt,y_txt,['*** ',str],...
            %             'FontSize',PLOT_STRUCT.regresslab_txt_size,...
            %             'FontName',PLOT_STRUCT.font_name,...
            %             'FontWeight','bold','Units','normalized');
            %     end
            % end
        end
    end
    cnt_g = cnt_g + length(conds);
end

%## set figure color, units, and size
%- set axes units, and size
set(ax,'box','off')
set(ax,'LineWidth',PLOT_STRUCT.ax_line_width)
set(ax,'FontName',PLOT_STRUCT.font_name,'FontSize',PLOT_STRUCT.font_size)
set(ax,'OuterPosition',[0 0 1 1]);
set(gca,'Position',PLOT_STRUCT.ax_position);  %[left bottom width height] This I what I added, You need to play with this
%- xlabel
xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized',...
    'FontSize',PLOT_STRUCT.x_label_fontsize,'FontWeight',PLOT_STRUCT.x_label_fontweight);
pos1=get(xlh,'Position');
pos1(1,2)=pos1(1,2)+PLOT_STRUCT.x_label_yoffset;
set(xlh,'Position',pos1);
%- xticks
if ~isempty(PLOT_STRUCT.xlim)
    set(ax,'XLim',PLOT_STRUCT.xlim)
else
    set(ax,'XLim',hold_xlim)
end
set(ax,'XTick',sort(xticks));
set(ax,'XTickLabel',xtick_labs);
% if ischar(PLOT_STRUCT.xtick_labs{1})
%     if isempty(PLOT_STRUCT.xtick_labs{1})
%         set(ax,'XTickLabel',PLOT_STRUCT.xtick_labs{1});
%     else
%         set(ax,'XTickLabel',PLOT_STRUCT.xtick_labs);
%     end
% elseif PLOT_STRUCT.xtick_labs{1}
%     set(ax,'XTickLabel',xtick_labs);
% else
%     set(ax,'XTickLabel',xtick_labs);
% end
xtickangle(75);
%- set ylabel
ylabel(ax,PLOT_STRUCT.y_label,'FontSize',PLOT_STRUCT.y_label_fontsize,...
    'FontWeight',PLOT_STRUCT.y_label_fontweight);
ylim(ax,PLOT_STRUCT.ylim);
%- title
title(ax,PLOT_STRUCT.title,'FontSize',PLOT_STRUCT.title_fontsize,...
    'FontWeight',PLOT_STRUCT.title_fontweight);       
%- set group labels
% axpos = get(ax,'Position');
% axxlim = get(ax,'XLim');
axylim = get(ax,'YLim');
try
    % if ~isempty(PLOT_STRUCT.group_labels)
    % 3 = 4.5991
    % 2 = 5
    % 1 = 7.0769
    if ~isempty(PLOT_STRUCT.group_labels) %~isundefined(PLOT_STRUCT.group_labels)
        cnt_g = 1;
        for i = 1:length(groups)
            % x = gticks(i)/(axlim(1)+axlim(2));
            % x = gticks(i)/(axlim(2)+axpos(1)*axpos(3));
            % x = (gticks(i)/(axlim(2)+axlim(1)))*axpos(3)+axpos(1);%PLOT_STRUCT.group_offsets(i);
            x = gticks(i);
            % y = 0.3; %GROUPLAB_YOFFSET
            y = PLOT_STRUCT.group_lab_yoffset*abs(axylim(2)-axylim(1));
            % text(ax,x,y,0,char(PLOT_STRUCT.group_labels(i)),...
            %     'FontSize',GROUP_LAB_FONTSIZE,'FontWeight',GROUP_LAB_FONTWEIGHT,'HorizontalAlignment','center',...
            %     'Units','normalized');
            tt = text(ax,x,y,0,char(PLOT_STRUCT.group_labels(i)),...
                'FontSize',PLOT_STRUCT.group_lab_fontsize,'FontWeight',PLOT_STRUCT.group_lab_fontweight,...
                'FontName',PLOT_STRUCT.font_name,'HorizontalAlignment','center',...
                'Units','data');
            tt.Units = 'normalized';
            tt.Position(2) = PLOT_STRUCT.group_lab_yoffset;
            cnt_g = cnt_g + length(conds);
        end
    end
catch e
    error('Error. Plotting group labels failed...\n\n%s',getReport(e));
end
%%

% hold on;
% hold off;
end
%% ================================================================== %%
function pp = cus_sigbracket(lbl,pVal,bx1,bx2,by1,by2,connector_offset,sig_offset_x,sig_offset_y)
    %- bx1(1) are x points for group 1's first vertical line
    %- bx1(2) are x points for gorup 1's second vertical line (spanning)
    %- bx2(1) are x points for group 2's first vertical line
    %- bx2(2) are x points for group 2's second vertical line (spanning)
    %- by1(1) are the lower y points for the beginning of the horizontal"box", for group 1
    %- by1(2) are the uppper y points for the end of the horizontal "box", for group 1
    %- by2(1) & by2(2) are for group 2.
    %- bracket 1
    bxx1 = [bx1(1),bx1(1),bx1(2),bx1(2)]; %[bx1 is x 1st vert line, x for 2nd vert line]
    byy1 = [by1(1),by1(2),by1(2),by1(1)]; 
    %- bracket 2
    bxx2 = [bx2(1),bx2(1),bx2(2),bx2(2)]; 
    byy2 = [by2(1),by2(2),by2(2),by2(1)];
    %- connector 
    bm1 = mean(bx1);
    bm2 = mean(bx2);
    [byymax,~] = max([by1,by2]);
    boffset = abs(by2(2)-by2(1))*connector_offset;
    byconn = [max(by1),byymax+boffset,byymax+boffset,max(by2)];
    bxconn = [bm1,bm1,bm2,bm2];
    % Now plot the sig line on the current axis
    hold on;
    plot(bxx1,byy1,'-k','LineWidth',0.5);
    plot(bxx2,byy2,'-k','LineWidth',0.5);
    pp = plot(bxconn,byconn,'-k','LineWidth',0.5);
    %- add label
    if lbl
        if pVal < 0.001
            lbl = '+++';
        elseif pVal > 0.001 && pVal < 0.01
            lbl = '++';
        elseif pVal > 0.01 && pVal < 0.05
            lbl = '+';
        end
        text(mean(bxconn)+sig_offset_x, max(byconn)*1.05+sig_offset_y, lbl) % the sig star sign
    end
    hold on;
end
%% ================================================================== %%
function pp = cus_sigline(xs,lbl,h,yv,pVal)
    
    if nargin==1
        y=gety(gca);
        lbl=[];
    elseif nargin==2
        y=gety(gca);
    elseif nargin==3
        y=gety(h);
    elseif nargin>=4
        y=yv;
    end
    
    
    % Now plot the sig line on the current axis
    hold on
    xs=[xs(1);xs(2)];
    if y > 0
        pp=plot(xs,[1;1]*y*1.1,'-k', 'LineWidth',0.5);%line
    else
        pp=plot(xs,[1;1]*y + 2,'-k', 'LineWidth',0.5);%line
    end
    % plot(mean(xs), y*1.15, '*k','MarkerSize',5)% the sig star sign
    if lbl
        if pVal < 0.001
            lbl = '***';
        elseif pVal > 0.001 & pVal < 0.01
            lbl = '**';
        elseif pVal > 0.01 & pVal < 0.05
            lbl = '*'
        end
        if y > 0
            text(mean(xs)*1, y*1.12, lbl)% the sig star sign
        else
            text(mean(xs)*1, y + 2, lbl)% the sig star sign
        end
    end
    % plot([1;1]*xs(1),[y*1.05,y*1.1],'-k', 'LineWidth',1);%left edge drop
    % plot([1;1]*xs(2),[y*1.05,y*1.1],'-k', 'LineWidth',1);%right edge drop
    hold on;
end
%--------------------------------------------------------------------------
% Helper function that Returns the largest single value of ydata in a given
% graphic handle. It returns the given value if it is not a graphics
% handle. 
function y=gety(h)
    %Returns the largest single value of ydata in a given graphic handle
    %h= figure,axes,line. Note that y=h if h is not a graphics
    if isgraphics(h) 
        switch(get(h,'type'))
            case {'line','hggroup','patch'} %{'line','hggroup','patch','scatter'},
                y=max(get(h,'ydata'));
                return;
            otherwise
                ys=[];
                hs=get(h,'children');
                for n=1:length(hs)
                    ys=[ys,gety(hs(n))];
                end
                y=max(ys(:));
        end
    else
        y=h;
    end
end