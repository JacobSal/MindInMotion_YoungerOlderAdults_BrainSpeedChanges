function [] = mim_custom_spec_plots(STUDY,ALLEEG,save_dir,varargin)
%MIM_GEN_CLUSTER_FIGS Summary of this function goes here
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

%## TIME
tic
%## DEFINE DEFAULTS
%- ADMIN
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
SPEC_ALL_POS=[500 300 1080 920];
SPEC_SING_POS=[16 582 420 360];
% DO_SINGLE_CLUSTER_PLOTS = true;
%-
CLUSTERS_TO_PLOT = 1:length(STUDY.cluster);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'CLUSTERS_TO_PLOT',CLUSTERS_TO_PLOT,@isnumeric);
addParameter(p,'SPEC_PARAMS',SPEC_PARAMS,@isnumeric);

parse(p,STUDY,ALLEEG,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
CLUSTERS_TO_PLOT = p.Results.CLUSTERS_TO_PLOT;
SPEC_PARAMS = p.Results.SPEC_PARAMS;
%% ===================================================================== %%
c_names = {STUDY.cluster(CLUSTERS_TO_PLOT).name};
fprintf('Plotting cluster numbers:'); fprintf('%i,',CLUSTERS_TO_PLOT(1:end-1)); fprintf('%i',CLUSTERS_TO_PLOT(end)); fprintf('\n');
fprintf('Associated cluster names:'); cellfun(@(x) fprintf('%s,',x),c_names(1:end-1)); fprintf('%s',c_names{end}); fprintf('\n');
%% (SPEC) Spec plot conds for des_i and all groups
%{
fprintf('Plotting Spectograms for Conditions...\n');
for des_i = 1:length(STUDY.design)
    std_specplot(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,...
        'freqrange',SPEC_PARAMS.plot_freqrange,'design',des_i);
    fig_i = get(groot,'CurrentFigure');
    fig_i.Position = SPEC_ALL_POS;
    %- set figure line colors
    cc = linspecer(length(STUDY.design(des_i).variable.value));
    iter = 1;
    for ax_i = 2:length(fig_i.Children)
        for d = 1:length(fig_i.Children(ax_i).Children)
            %- pane 1
            set(fig_i.Children(ax_i).Children(d),'LineWidth',1.5);
            set(fig_i.Children(ax_i).Children(d),'Color',horzcat(cc(iter,:),0.6));
            if iter == size(cc,1)
                iter = 1;
            else
                iter = iter + 1;
            end                
        end
        %- zoom may not be need, but sort of cleans things up
%         camzoom(fig_i.Children(ax_i),1.1);
    end
    saveas(fig_i,[save_dir filesep sprintf('allSpecPlot_des%i.pdf',des_i)]);
end
%}
%% (SINGLE CLUSTSER PLOTS)
for i = 1:length(CLUSTERS_TO_PLOT)
    clust_i = CLUSTERS_TO_PLOT(i);
    if ~isempty(STUDY.cluster(clust_i).sets)
        %- (SPEC) Spec plot conds for des_i and all groups
        fprintf('Plotting Spectograms for Conditions...\n');
        for des_i = 1:length(STUDY.design)
            std_specplot(STUDY,ALLEEG,'clusters',clust_i,...
                'freqrange',SPEC_PARAMS.plot_freqrange,...
                'design',des_i);
            fig_i = get(groot,'CurrentFigure');
            hold on;
            set(fig_i,'Color','w')
            set(fig_i,'Units','inches','Position',[3 3 6 5])
            set(fig_i,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
%             fig_i.Position = %SPEC_SING_POS;
%             disp(get(fig_i,'Position'));
            %- set figure line colors
            cc = linspecer(length(STUDY.design(des_i).variable.value));
            iter = 1;
            for d = 1:length(fig_i.Children(2).Children)
                %- pane 1
                set(fig_i.Children(2).Children(d),'LineWidth',1.5);
                set(fig_i.Children(2).Children(d),'Color',horzcat(cc(iter,:),0.6));
                if iter == size(cc,1)
                    iter = 1;
                else
                    iter = iter + 1;
                end                
            end
            %- zoom may not be need, but sort of cleans things up
            set(fig_i.Children(2),'FontSize',13);
            set(fig_i.Children(3),'FontSize',13);
            set(fig_i.Children(2),'Position',[0.20,0.20,0.7,0.7]) %Default:[0.26,0.26,0.54,0.51]; Position::[left margin, lower margin, right margin, upper margin]
            set(fig_i.Children(3),'Position',[0.20,0.20-0.0255,0.7,0.0255]) %Default:[0.26,0.2345,0.54,0.0255]
            set(fig_i.Children(1),'Location','northeast') %reset Legend
            for c_i = 2:length(fig_i.Children)
                set(fig_i.Children(c_i),'FontName','Arial','FontSize',12,'FontWeight','bold')
%                 set(fig_i.Children(i),'OuterPosition',[0 0 1 1]);
            end
            hold off;
            drawnow;
            exportgraphics(fig_i,[save_dir filesep sprintf('allSpecPlot_des%i_cl%i.jpg',des_i,clust_i)],'Resolution',300);
        end
    else
        fprintf('\nNo Subjects in Cluster %s, number: %i\n',c_names{cllust_i},cllust_i);
    end
end

end

