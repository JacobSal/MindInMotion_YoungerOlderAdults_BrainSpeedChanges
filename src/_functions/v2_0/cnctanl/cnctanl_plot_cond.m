function [handles] = cnctanl_plot_cond(to_from_conn,conditions,two_clusts,subj_clusts_comps,save_dir,varargin)
%UNTITLED Summary of this function goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%% Define Defaults
DO_LEGEND = true;
SAVE_FIG = true;
LEGEND_CHARS = {};
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'to_from_conn',@isnumeric);
addRequired(p,'conditions',@iscell);
addRequired(p,'two_clusts',@iscell);
addRequired(p,'subj_clusts_comps',@isnumeric);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
addOptional(p,'LEGEND_CHARS',LEGEND_CHARS,@iscell);

%## PARAMETER
parse(p,to_from_conn, conditions, two_clusts, subj_clusts_comps, save_dir, varargin{:});
%## SET DEFAULTS
LEGEND_CHARS = p.Results.LEGEND_CHARS;
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## COLORS
% cMap = jet;
% cMap = parula; %parula(size(subj_clusts_comps,1));
cMap = linspecer(size(to_from_conn,2));
color.lightGreen  = [186,228,179]/255;
color.Green = [44 162 95]/255;
color.lightRed = [252 174 145]/255;
color.Red = [222,45,38]/255;
%% Subject Specific traces
MARK_CHARS = {':o','-o','--o','-x','--x',':x','-^','--^',':^'};
%- parmas
cond_picks = 1:length(conditions);
X_LIM = [0.5,length(cond_picks)+0.5];
% Y_LIM = [-0.01,0.15];
Y_LIM = [0,0.15];
% Y_LIM = [-0.1,0.1];
% Y_LIM = [-0.4,0.4];

x = 1:length(cond_picks);
%- figure
fig = figure;
hold on;
h = [];
mark_iter = 1;
mark_char = 'o-';
for subj_i = 1:size(to_from_conn,2)
    if isempty(LEGEND_CHARS)
        dispname_h = sprintf('%i<->%i) subject %i',subj_clusts_comps(subj_i,1),subj_clusts_comps(subj_i,2),subj_i);
%         dispname_hh = sprintf('%i<->%i) subject %i',subj_clusts_comps(subj_i,2),subj_clusts_comps(subj_i,1),subj_i);
    else
        dispname_h = LEGEND_CHARS{subj_i};
%         dispname_hh = LEGEND_CHARS{subj_i};
    end
    h = [h, plot(x,squeeze(to_from_conn(1,subj_i,:))',mark_char,'DisplayName',dispname_h)];
    h(subj_i).Color = cMap(subj_i,:); %cMap(cc,:);
%     hh = plot(x,squeeze(to_from_conn(2,subj_i,:))','x-','DisplayName',dispname_hh);
%     hh.Color = cMap(cc,:);
    if mod(subj_i,2)
        mark_char = MARK_CHARS{mark_iter};
    else
        mark_char = MARK_CHARS{mark_iter};
    end
    if mark_iter < length(MARK_CHARS)
        mark_iter = mark_iter + 1;
    else
        mark_iter = 1;
    end
end
hold off;
xlim(X_LIM);
ylim(Y_LIM);
ylabel('Connectivity');
xlabel('Condition');
fig.CurrentAxes.XTick = 0:length(conditions(cond_picks));
fig.CurrentAxes.XTickLabel = {'',conditions{cond_picks},''};
title(sprintf('Connectivity of %s<->%s',two_clusts{1},two_clusts{2}))
fig.Position = [200,400,920,720];
fig_i = get(groot,'CurrentFigure');
if DO_LEGEND
    legend(h,'Location','northeast')
end
saveas(fig_i,[save_dir filesep sprintf('allsubj_%s-%s_%s.jpg',two_clusts{1},two_clusts{2},[conditions{:}])]);
if SAVE_FIG
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep sprintf('allsubj_%s-%s_%s.fig',two_clusts{1},two_clusts{2},[conditions{:}])]);
end
%% AVERAGES ACROSS SUBJECTS AND CONDITION
%- params
% avg_conn = squeeze(nanmean(to_from_conn,2));
avg_conn = squeeze(nanmedian(to_from_conn,2));
std_conn = squeeze(nanstd(to_from_conn,[],2));
X_LIM = [0.5,length(cond_picks)+0.5];
%- plot
x = 1:length(cond_picks);
h = [];
fig = figure; 
hold on;
h(1) = scatter(x,avg_conn(1,:),20,color.Green,'filled','DisplayName',sprintf('From %s',two_clusts{1})); 
h(2) = scatter(x,avg_conn(2,:),20,color.Red,'filled','DisplayName',sprintf('From %s',two_clusts{2})); 
h(3) = plot(avg_conn(1,:),'DisplayName','','Color',color.Green); 
h(4) = plot(avg_conn(2,:),'DisplayName','','Color',color.Red); 
h(5) = errorbar(avg_conn(1,:),std_conn(1,:),'DisplayName','','Color',color.Green);
h(6) = errorbar(avg_conn(2,:),std_conn(2,:),'DisplayName','','Color',color.Red);
xlim(X_LIM);
ylim(Y_LIM);
ylabel('Connectivity');
xlabel('Condition');
fig.CurrentAxes.XTick = 0:length(conditions);
fig.CurrentAxes.XTickLabel = {'',conditions{cond_picks},''};
title(sprintf('Connectivity of %s<->%s',two_clusts{1},two_clusts{2}))
legend(h(1:2),'Location','southeast')
hold off;
if SAVE_FIG
    fig_i = get(groot,'CurrentFigure');
    saveas(fig_i,[save_dir filesep sprintf('avgbtwn_%s-%s_%s.fig',two_clusts{1},two_clusts{2},[conditions{:}])]);
    saveas(fig_i,[save_dir filesep sprintf('avgbtwn_%s-%s_%s.jpg',two_clusts{1},two_clusts{2},[conditions{:}])]);
end
handles = h;
end

