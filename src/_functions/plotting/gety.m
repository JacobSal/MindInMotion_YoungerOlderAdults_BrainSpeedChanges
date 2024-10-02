function y=gety(ax,varargin)
% Helper function that Returns the largest single value of ydata in a given
% graphic handle. It returns the given value if it is not a graphics
% handle. 
% Returns the largest single value of ydata in a given graphic handle
% h= figure,axes,line. Note that y=h if h is not a graphics
%## TIME
tic
%## DEFINE DEFAULTS
graphic_types = {'line','hggroup','patch'}; %,'scatter'
%## PARSER
p = inputParser;
%- REQUIRED
addRequired(p,'ax',@isgraphics)
addOptional(p,'graphic_types',graphic_types,@iscell)
parse(p,ax,varargin{:});
%% ===================================================================== %%
if isgraphics(ax) 
    switch get(ax,'type')
        case graphic_types %{'line','hggroup','patch','scatter'},
            y = max(get(ax,'ydata'));
            return;
        otherwise
            ys = [];
            hs = get(ax,'children');
            for n = 1:length(hs)
                ys = [ys,gety(hs(n))];
            end
            y = max(ys(:));
    end
else
    y = ax;
end
toc
end

