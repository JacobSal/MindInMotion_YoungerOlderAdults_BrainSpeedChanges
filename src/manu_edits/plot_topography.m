function [ax1] = plot_topography(fig,STUDY,cl_n, ...
    groups,g_chars,g_chars_topo, ...
    ax_position,im_resize)
%PLOT_TOPOGRAPHY Summary of this function goes here
%   Detailed explanation goes here

%##
hold on;
%## PARAMS
AXES_DEFAULT_PROPS = {'box','off', ...
    'xtick',[], ...
    'ytick',[], ...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1]};
%--
PLOT_STRUCT = struct('font_name','Arial', ...
    'font_size',8, ...
    'font_weight','bold', ...
    'plot_title','A)');

%## TOPO PLOT
[~,h] = std_topoplot_CL(STUDY,cl_n,'together');
set(h,'color','w')
set(h,'Units','normalized');    
ax = get(h,'CurrentAxes');
colormap(h,linspecer);   
% pos = get(ax,'Position');
% opos = get(ax,'OuterPosition');
ac = get(ax,'Children');
ax1 = axes('Parent',fig,...
    'DataAspectRatio',get(ax,'DataAspectRatio'));
copyobj(ac,ax1);
set(ax1,AXES_DEFAULT_PROPS{:});
colormap(ax1,linspecer)
%-- group title counts
g_counts = cell(length(groups),1);
for g_i = 1:length(groups)
    g_inds = cellfun(@(x) strcmp(x,g_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cl_n).sets).group});
    if length(g_chars_topo{g_i}) == 1 || ischar(g_chars_topo{g_i})
        g_counts{g_i} =sprintf('%s N=%i',g_chars_topo{g_i},sum(g_inds));
    else
        g_counts{g_i} =sprintf('%s\n%s N=%i',g_chars_topo{g_i}{1},g_chars_topo{g_i}{2},sum(g_inds));
    end
end
ax1.Title.String = g_counts; %sprintf('N=%i',length(STUDY.cluster(cl_i).sets));
ax1.Title.Interpreter = 'none';
%-- axes edits
ax1.FontSize = PLOT_STRUCT.font_size; %PLOT_STRUCT.font_size;
ax1.FontName = PLOT_STRUCT.font_name;
ax1.FontWeight = PLOT_STRUCT.font_weight;
ax1.OuterPosition = [0,0,1,1];
ax1.Units = 'Normalized';
pp = get(ax1,'Position');    
ax1.Position = [ax_position(1),ax_position(2), ...
    ax_position(3)+pp(3)*im_resize,ax_position(4)+pp(4)*im_resize];  %[left bottom width height]
close(h)
end

