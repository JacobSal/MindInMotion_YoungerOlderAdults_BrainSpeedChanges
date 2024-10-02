function [h_out] = cus_sigbracket(ax,p_value,bracket_1,bracket_2,varargin)
%CUS_SIGBRACKET Summary of this function goes here
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
axylim = get(ax,'YLim');
%-
DEF_LINE_STRUCT = struct('sig_sign','+',...
    'line_specs',{{'LineStyle','-','LineWidth',1,'Color','k'}},...
    'text_specs',{{'FontSize',10,'FontName','Arial','FontWeight','bold'}},...
    'bracket_conn',[],...
    'conn_offset_y_upper',[],...
    'bracket_offset_y_upper',0,...
    'bracket_offset_y_lower',0,...
    'sig_levels',[0.05,0.01,0.001],...
    'sig_offset_x',0,...
    'sig_offset_y',[]);
%## PARSER
p = inputParser;
%- REQUIRED
addRequired(p,'ax',@isgraphics)
addRequired(p,'p_value',@isnumeric);
addRequired(p,'bracket_1',@isnumeric); % e.g., [[1,1,2,2];[1,2,2,1]] (1,:) x-coords; (2,:) y-coords
addRequired(p,'bracket_2',@isnumeric); % e.g., [[3,3,4,4];[1,2,2,1]] (1,:) x-coords; (2,:) y-coords
%- PARAMETER
addParameter(p,'LINE_STRUCT',DEF_LINE_STRUCT,@(x) validate_struct(x,DEF_LINE_STRUCT));
parse(p,ax,p_value,bracket_1,bracket_2,varargin{:});
%- SET DEFAULTS
LINE_STRUCT = p.Results.LINE_STRUCT;
LINE_STRUCT = set_defaults_struct(LINE_STRUCT,DEF_LINE_STRUCT);
%% ===================================================================== %%
%- REDUNDANCY
LINE_STRUCT.sig_levels = sort(LINE_STRUCT.sig_levels,'descend');
%- BRACKET OFFSETS
% NOTE: bracket_offset_y_upper should be negative if bracket in neg.
% if (any(bracket_1(2,[2,3]) < 0) || any(bracket_2(2,[2,3])) && LINE_STRUCT.bracket_offset_y_upper > 0
%     LINE_STRUCT.bracket_offset_y_upper = -LINE_STRUCT.bracket_offset_y_upper;
% end
% if (any(bracket_1(2,[1,4]) < 0) || any(bracket_2(2,[1,4])) && LINE_STRUCT.bracket_offset_y_lower > 0
%     LINE_STRUCT.bracket_offset_y_lower = -LINE_STRUCT.bracket_offset_y_lower;
% end
bracket_1(2,[2,3]) = bracket_1(2,[2,3])+LINE_STRUCT.bracket_offset_y_upper; 
bracket_1(2,[1,4]) = bracket_1(2,[1,4])+LINE_STRUCT.bracket_offset_y_lower;
bracket_2(2,[2,3]) = bracket_2(2,[2,3])+LINE_STRUCT.bracket_offset_y_upper;
bracket_2(2,[1,4]) = bracket_2(2,[1,4])+LINE_STRUCT.bracket_offset_y_lower;
%- 
if isempty(LINE_STRUCT.conn_offset_y_upper)
    conn_offset_y_upper = 0.05*(abs(axylim(2)-axylim(1)));
else
    conn_offset_y_upper = LINE_STRUCT.conn_offset_y_upper;
end
%- BRACKET CONNECTOR
if isempty(LINE_STRUCT.bracket_conn)
    bm1 = mean(bracket_1(1,:));
    bm2 = mean(bracket_2(1,:));
    [byymax,~] = max([bracket_1(2,:),bracket_2(2,:)]);
    bracket_conn = [[bm1,bm1,bm2,bm2];...
        [max(bracket_1(2,:)),byymax+conn_offset_y_upper,byymax+conn_offset_y_upper,max(bracket_2(2,:))]];
else
    bracket_conn = LINE_STRUCT.bracket_conn;
end
%- SIGNIFICANCE SIGN OFFSET
if isempty(LINE_STRUCT.sig_offset_y)
    sig_offset_y = 0.05*(abs(axylim(2)-axylim(1)));
else
    sig_offset_y = LINE_STRUCT.sig_offset_y;
end
%## PLOT
p1 = plot(ax,bracket_1(1,:),bracket_1(2,:),LINE_STRUCT.line_specs{:});
p2 = plot(ax,bracket_2(1,:),bracket_2(2,:),LINE_STRUCT.line_specs{:});
p3 = plot(ax,bracket_conn(1,:),bracket_conn(2,:),LINE_STRUCT.line_specs{:});
h_out = [p1,p2,p3];
%- add label
chkl = find((p_value <= LINE_STRUCT.sig_levels),1,'last');
if ~isempty(chkl) && ~strcmp(LINE_STRUCT.sig_sign,'')
    sig_sign = repmat(LINE_STRUCT.sig_sign,[1,chkl]);
    tx = text(ax,mean(bracket_conn(1,:))+LINE_STRUCT.sig_offset_x, max(bracket_conn(2,:))+sig_offset_y,sig_sign,...
        LINE_STRUCT.text_specs{:}); % the sig star sign
    h_out = [h_out,tx];
end
hold on;
end
%##
%- bx1(1) are x points for group 1's first vertical line
%- bx1(2) are x points for gorup 1's second vertical line (spanning)
%- bx2(1) are x points for group 2's first vertical line
%- bx2(2) are x points for group 2's second vertical line (spanning)
%- by1(1) are the lower y points for the beginning of the horizontal"box", for group 1
%- by1(2) are the uppper y points for the end of the horizontal "box", for group 1
%- by2(1) & by2(2) are for group 2.
%- bracket 1
% bxx1 = [bx1(1),bx1(1),bx1(2),bx1(2)]; %[bx1 is x 1st vert line, x for 2nd vert line]
% byy1 = [by1(1),by1(2),by1(2),by1(1)]; 
%- bracket 2
% bxx2 = [bx2(1),bx2(1),bx2(2),bx2(2)]; 
% byy2 = [by2(1),by2(2),by2(2),by2(1)];
% Now plot the sig line on the current axis
