function [fig] = plot_txf(allersp,varargin)
%PLOT_TXF Summary of this function goes here
%   Detailed explanation goes here
alltimes = [];
allfreqs = [];
PLOT_STRUCT = struct('figure_position',[100,100,350,350],...
    'xtick_label','Percent Gait Cycle (%)',...
    'ytick_label','Frequency (Hz)',...
    'clim',[-2,2],...
    'font_size',12,...
    'freq_lims',[],...
    'time_lims',[],...
    'title',[]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@isnumeric);
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
%##
fig = figure('color','white','position',PLOT_STRUCT.figure_position,'renderer','Painters');
set(fig,'Units','inches','Position',[3 3 5 5])
set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
tftopo(allersp,alltimes,allfreqs,'limits',... 
    [PLOT_STRUCT.time_lims nan nan PLOT_STRUCT.clim],...
    'logfreq','native');
% contourf(allersp);
ax = gca;
hold on;
colormap(linspecer);
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
ylabel(PLOT_STRUCT.ytick_label,'FontSize',PLOT_STRUCT.font_size,'fontweight','bold');
xlabel(PLOT_STRUCT.xtick_label,'FontSize',PLOT_STRUCT.font_size);
% set(ax,'title',PLOT_STRUCT.title);
title(PLOT_STRUCT.title);
colorbar;
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
            fprintf(2,'\nValue must be type %s, but is type %s\n',class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end
