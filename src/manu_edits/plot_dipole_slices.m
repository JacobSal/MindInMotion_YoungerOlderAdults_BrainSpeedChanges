function [ax1,ax2,ax3] = plot_dipole_slices(fig,dip_fig_path, ...
    im_resize,paper_size,ax_position,label_position)
%PLOT_DIPOLE_SLICES Summary of this function goes here
%   Detailed explanation goes here
%## DIPOLE PLOTS =================================================== %%
hold on;
%--
PLOT_STRUCT = struct('font_name','Arial', ...
    'plot_title','A)');
%--
% AX_INIT_VERT_DIP = 0.83; %0.79; %0.8 --- %7.2; % inches for some reason, maybe a bug with normal units
% AX_INIT_HORIZ_DIP = 0.3846; %2.5; % inches for some reason, maybe a bug with normal units
% AX_INIT_HORIZ = 0.09;
% LAB_A_YOFFSET = -0.16;
% LAB_A_XOFFSET = -0.125;
% ax_position = [AX_INIT_HORIZ_DIP,AX_INIT_VERT_DIP,0,0];
% label_position = [AX_INIT_HORIZ+LAB_A_XOFFSET+(0.1/2),1+LAB_A_YOFFSET+(0.1/2),.1,.1];
%--
ax1_adjusts = [1.02,-0.04,-0.145];
ax2_adjusts = [1.55,0.15,-0.175];
ax3_adjusts = [1.2725,0.22,-0.11];
%--
DIP_SIZE = 7;
DIP_LINEWIDTH = 0.1;
DIP_BORDER_CDATA = [0,0,0];
DIP_MARKERFACEALPHA = 0.75;
%--
%## OPEN DIPOLE FIGURE & EXTRACT AXES
tmp = openfig(dip_fig_path);
cnt = length(tmp.Children(end).Children);
try
    while cnt > 3
        tmp.Children(end).Children(cnt).SizeData = DIP_SIZE;
        tmp.Children(end).Children(cnt).LineWidth = DIP_LINEWIDTH;
        tmp.Children(end).Children(cnt).CData = DIP_BORDER_CDATA; 
        tmp.Children(end).Children(cnt).MarkerFaceAlpha = DIP_MARKERFACEALPHA;
        cnt = cnt - 1;
    end
catch
    fprintf('not scatter 3d');
end
ax = get(tmp,'CurrentAxes');
view(ax,[90,0]);
daspect = get(ax,'DataAspectRatio');
ac = get(ax,'Children');
XLIM = [ac(1).XData(1,1),ac(1).XData(1,2)];
YLIM = [ac(1).YData(1,1),ac(1).YData(2,1)];
ZLIM = [ac(2).ZData(1,1),ac(2).ZData(2,1)];  

%## AXES 1
hold on;
ax1 = axes('Parent',fig, ...
    'DataAspectRatio',daspect, ...
    'Units','normalized');
copyobj(ac,ax1);
set(ax1,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1]);
ax1.XRuler.FirstCrossoverValue = 0;
ax1.YRuler.FirstCrossoverValue = 0;
ax1.ZRuler.FirstCrossoverValue = 0;
view([0,90])
camzoom(ax1,1.1^2)
set(ax1,'YLim',YLIM)
set(ax1,'XLim',XLIM)
drawnow;
pp = get(ax1,'Position');
%-
AX1_ADJ = ax1_adjusts(1); %1.02;
AX1_Y_O = ax1_adjusts(3); %-0.145; %-0.17;
AX1_X_O = ax1_adjusts(2); %-0.04; %-0.17;
ax1_x = pp(3)*im_resize*AX1_ADJ/paper_size(3);
ax1_y = pp(4)*im_resize*AX1_ADJ/paper_size(3);
set(ax1,'OuterPosition',[0,0,1,1], ...
    'Position',[ax_position(1)+ax1_x*AX1_X_O,ax_position(2)+ax1_y*AX1_Y_O,ax1_x,ax1_y]);

%## AXES 2
hold on;
ax2 = axes('Parent',fig, ...
    'DataAspectRatio',daspect, ...
    'Units','normalized');
copyobj(ac,ax2);
%-
set(ax2,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1]);
ax2.XRuler.FirstCrossoverValue = 0;
ax2.YRuler.FirstCrossoverValue = 0;
ax2.ZRuler.FirstCrossoverValue = 0;
set(ax2,'View',[90,0]);
camzoom(ax2,1.1^2)
set(ax2,'YLim',YLIM)
set(ax2,'ZLim',ZLIM)
drawnow;
pp2 = get(ax2,'Position');
AX2_ADJ = ax2_adjusts(1); %1.55;
AX2_Y_O = ax2_adjusts(3); %-0.175;
AX2_X_O = ax2_adjusts(2); %;0.15;
ax2_x = pp2(3)*im_resize*AX2_ADJ/paper_size(3);
ax2_y = pp2(4)*im_resize*AX2_ADJ/paper_size(4);    
set(ax2,'OuterPosition',[0,0,1,1], ...
    'Position',[ax_position(1)+ax1_x+(ax2_x*AX2_X_O),ax_position(2)+(ax2_y*AX2_Y_O),ax2_x,ax2_y]);

%## AXES 3
hold on;
ax3 = axes('Parent',fig, ...
    'DataAspectRatio',daspect, ...
    'Units','normalized');
copyobj(ac,ax3);
%-
set(ax3,'box','off', ...
    'xtick',[], ...
    'ytick',[],...
    'ztick',[], ...
    'xcolor',[1,1,1], ...
    'ycolor',[1,1,1], ...
    'zcolor',[1,1,1]);
ax3.XRuler.FirstCrossoverValue = 0;
ax3.YRuler.FirstCrossoverValue = 0;
ax3.ZRuler.FirstCrossoverValue = 0;
set(ax3,'view',[0,0])
set(ax3,'ZLim',ZLIM)
set(ax3,'XLim',XLIM)
camzoom(ax3,1.1^2)
pp3 = get(ax3,'Position');
AX3_ADJ = ax3_adjusts(1); %1.2725;
AX3_Y_O = ax3_adjusts(3); %-0.11;
AX3_X_O = ax3_adjusts(2); %0.22;
ax3_x = pp3(3)*im_resize*AX3_ADJ/paper_size(3);
ax3_y = pp3(4)*im_resize*AX3_ADJ/paper_size(4);    
set(ax3,'OuterPosition',[0,0,1,1], ...
    'Position',[ax_position(1)+ax1_x+ax2_x+(ax2_x*AX2_X_O)+(ax3_x*AX3_X_O),ax_position(2)+(ax3_y*AX3_Y_O),ax3_x,ax3_y]);
close(tmp)  

%## TITLE
annotation(fig,'textbox',label_position,...
    'String',PLOT_STRUCT.plot_title, ...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top', ...
    'LineStyle','none', ...
    'FontName',PLOT_STRUCT.font_name,...
    'FontSize',14, ...
    'FontWeight','Bold', ...
    'Units','normalized');

end

