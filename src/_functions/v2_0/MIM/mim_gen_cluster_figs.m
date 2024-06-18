function [] = mim_gen_cluster_figs(STUDY,ALLEEG,save_dir,varargin)
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
cat_logo();
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
% SPEC_ALL_POS=[500 300 1080 920];
% SPEC_SING_POS=[16 582 420 360];
DIP_ALL_POS=[16 100 1080 920];
DIP_SING_POS=[16 582 420 360];
TOPO_ALL_POS=[16 100 1240 920];
TOPO_SING_POS=[16 100 300 350];
DO_SINGLE_CLUSTER_PLOTS = true;
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
colors{1}  = [1 1 1];            % White
colors{2}  = [1 1 0];            % Yellow
colors{3}  = [221,52,151]/255;            % Pink
colors{4}  = [1 0 0];            % Red
colors{5}  = [250 140 0]/255;    % oragne
colors{6}  = [210 173 255]/255;   % purple
colors{7}  = [0.5 0.5 0];        % Olive
colors{8}  = [0.5 0 0.5];        % Purple
colors{9}  = [0.5 0 0];          % Maroon
colors{10} = [0 1 1];            % Aqua
colors{11} = [0 1 0];            % Lime
colors{12} = [0 0.5 0.5];        % Teal
colors{13} = [0 0.5 0];          % Green
colors{14} = [0 0 1];            % Blue
colors{15} = [0 0 0.5];          % Navy
colors{16} = [0.8 0.8 0.8];            % Gray
% Choosing and sorting 13 colors for clusters: Red, Green, Blue,
% Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, Gray, and White
colors = colors([4 11 14 2 13 10 5 6 15 16 1 7 9 3]);
%% ===================================================================== %%
c_names = {STUDY.cluster(CLUSTERS_TO_PLOT).name};
fprintf('Plotting cluster numbers:'); fprintf('%i,',CLUSTERS_TO_PLOT(1:end-1)); fprintf('%i',CLUSTERS_TO_PLOT(end)); fprintf('\n');
fprintf('Associated cluster names:'); cellfun(@(x) fprintf('%s,',x),c_names(1:end-1)); fprintf('%s',c_names{end}); fprintf('\n');
%## Update 7/20/2022 - Plot clustering results for AHA proposal
% Plot scalp topographs which also need to be averaged? 
if ~isfield(STUDY.cluster,'topo') 
    STUDY.cluster(1).topo = [];
end
for i = 1:length(CLUSTERS_TO_PLOT) % For each cluster requested
    clust_i = CLUSTERS_TO_PLOT(i);
    disp(clust_i)
    if isempty(STUDY.cluster(clust_i).topo)
        % Using this custom modified code to allow taking average within participant for each cluster
        STUDY = std_readtopoclust_CL(STUDY,ALLEEG,clust_i);
    end
end
%% (TOPO)
figure;
std_topoplot_CL(STUDY,CLUSTERS_TO_PLOT,'together');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',TOPO_ALL_POS,'color','w')
for c = 2:length(fig_i.Children)
    fig_i.Children(c).Title.Interpreter = 'none';
end
saveas(fig_i,fullfile(save_dir,'Cluster_topo_avg.jpg'));
saveas(fig_i,fullfile(save_dir,'Cluster_topo_avg.fig'));
%% (DIPOLE) Plot all dipole fit locations
% std_dipplot(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,'figure','off');
% fig_i = get(groot,'CurrentFigure');
% set(fig_i,'position',DIP_ALL_POS,'color','w')
% % saveas(fig_i,[save_dir filesep 'dipplot_seperatepanes.pdf']);
% saveas(fig_i,[save_dir filesep 'dipplot_seperatepanes.jpg']);
%% (DIPOLE) Plot dipole fit locations after averaged within participants
% std_dipplot_CL(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,'figure','off','mode','together_averaged');
% fig_i = get(groot,'CurrentFigure');
% set(fig_i,'position',DIP_SING_POS,'color','w')
% saveas(fig_i,[save_dir filesep 'Cluster_dipole_seperatepanes_avg.jpg']);
%% (DIPOLE) Plot dipole clusters 
STUDY.etc.dipparams.centrline = 'off';
std_dipplot_CL(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,'figure','off',...
    'mode','together_averaged_only','spheres','off','projlines','off',...
    'color_override',colors);
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',DIP_SING_POS,'color','w')
camzoom(1.1);
% camzoom(1.2^2);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.svg')],'ContentType','vector','Resolution',300);
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.fig')]);
view([45,0,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.svg')],'ContentType','vector','Resolution',300);
view([0,-45,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.svg')],'ContentType','vector','Resolution',300);
%{
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.jpg')]);
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')]);
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.fig')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal.jpg')]);
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal.jpg')]);
saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top.pdf')]);
%}
%% (DIPOLE) Plot dipole fit locations after averaged within participants and
% different clusters have different colors
STUDY.etc.dipparams.centrline = 'off';
std_dipplot_CL(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
fig_i = get(groot,'CurrentFigure');
set(fig_i,'position',DIP_SING_POS,'color','w')
% camzoom(1);
camzoom(1.2^2);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.svg')],'ContentType','vector','Resolution',300);
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
view([45,0,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.svg')],'ContentType','vector','Resolution',300);
view([0,-45,0])
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.jpg')],'Resolution',300);
exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.tiff')],'Resolution',900);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.pdf')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.eps')],'ContentType','vector','Resolution',300);
% exportgraphics(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.svg')],'ContentType','vector','Resolution',300);
%{
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.jpg')]);
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top.fig')]);
view([45,0,0])
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal.jpg')]);
view([0,-45,0])
saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal.jpg')]);
%}
%% (SPEC) Spec plot conds for des_i and all groups
% fprintf('Plotting Spectograms for Conditions...\n');
% for des_i = 1:length(STUDY.design)
%     std_specplot(STUDY,ALLEEG,'clusters',CLUSTERS_TO_PLOT,...
%         'freqrange',SPEC_PARAMS.plot_freqrange,'design',des_i);
%     fig_i = get(groot,'CurrentFigure');
%     fig_i.Position = SPEC_ALL_POS;
%     %- set figure line colors
%     cc = linspecer(length(STUDY.design(des_i).variable.value));
%     iter = 1;
%     for ax_i = 2:length(fig_i.Children)
%         for d = 1:length(fig_i.Children(ax_i).Children)
%             %- pane 1
%             set(fig_i.Children(ax_i).Children(d),'LineWidth',1.5);
%             set(fig_i.Children(ax_i).Children(d),'Color',horzcat(cc(iter,:),0.6));
%             if iter == size(cc,1)
%                 iter = 1;
%             else
%                 iter = iter + 1;
%             end                
%         end
%         %- zoom may not be need, but sort of cleans things up
%         camzoom(fig_i.Children(ax_i),1.1);
%     end
%     saveas(fig_i,[save_dir filesep sprintf('allSpecPlot_des%i.jpg',des_i)]);
% end
%% (SINGLE CLUSTSER PLOTS)
if DO_SINGLE_CLUSTER_PLOTS
    
%     COLOR_OPTS = colors; %linspecer(length(CLUSTERS_TO_PLOT)+1);
%     colors = cell(1,size(COLOR_OPTS,1));
%     for i = 1:size(COLOR_OPTS,1)
%         colors{i} = COLOR_OPTS(i,:);
%     end
    for i = 1:length(CLUSTERS_TO_PLOT)
        cluster_i = CLUSTERS_TO_PLOT(i);
        STUDY.etc.dipparams.centrline = 'off';
        std_dipplot_CL(STUDY,ALLEEG,'clusters',cluster_i,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'position',DIP_SING_POS,'color','w')
        for j = 1:length(fig_i.Children(2).Children)
            try
                fig_i.Children(2).Children(j).Color = colors{i};
            catch
                fprintf('Can''t change the color of this child\n')
            end
        end
        camzoom(1.2^2);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.jpg',cluster_i)],'Resolution',300);
        exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.tiff',cluster_i)],'Resolution',900);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.pdf',cluster_i)],'ContentType','vector','Resolution',300);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_top.svg',cluster_i)],'ContentType','vector','Resolution',300);
    %     saveas(fig_i,[save_dir filesep sprintf('new_dipplot_alldipspc_top.fig')]);
        view([45,0,0])
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.jpg',cluster_i)],'Resolution',300);
        exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.tiff',cluster_i)],'Resolution',900);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.pdf',cluster_i)],'ContentType','vector','Resolution',300);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_coronal.svg',cluster_i)],'ContentType','vector','Resolution',300);
        view([0,-45,0])
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.jpg',cluster_i)],'Resolution',300);
        exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.tiff',cluster_i)],'Resolution',900);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.pdf',cluster_i)],'ContentType','vector','Resolution',300);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_dipplot_alldipspc_sagittal.svg',cluster_i)],'ContentType','vector','Resolution',300);
        %- 
        %## TOPO PLOTS
        if ~isfield(STUDY.cluster,'topo') 
            STUDY.cluster(1).topo = [];
        end
        for j = 1:length(CLUSTERS_TO_PLOT) % For each cluster requested
            clust_i = CLUSTERS_TO_PLOT(j);
            disp(clust_i)
            if isempty(STUDY.cluster(clust_i).topo)
                % Using this custom modified code to allow taking average within participant for each cluster
                STUDY = std_readtopoclust_CL(STUDY,ALLEEG,clust_i);
            end
        end
        figure;
        std_topoplot_CL(STUDY,cluster_i,'together');
        fig_i = get(groot,'CurrentFigure');
        set(fig_i,'position',TOPO_SING_POS,'color','w')
        for c = 2:length(fig_i.Children)
            fig_i.Children(c).Title.Interpreter = 'none';
            fig_i.Children(c).FontSize = 13;
            fig_i.Children(c).FontName = 'Arial';
        end
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_cluster_topo_avg.jpg',cluster_i)],'Resolution',300);
        exportgraphics(fig_i,[save_dir filesep sprintf('%i_cluster_topo_avg.tiff',cluster_i)],'Resolution',900);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_cluster_topo_avg.pdf',cluster_i)],'ContentType','vector','Resolution',300);
%         exportgraphics(fig_i,[save_dir filesep sprintf('%i_cluster_topo_avg.svg',cluster_i)],'ContentType','vector','Resolution',300);
    end
    %     saveas(fig_i,fullfile(save_dir,sprintf('%i_cluster_topo_avg.fig',cluster_i)));
    %{
    for i = 1:length(CLUSTERS_TO_PLOT)
        clust_i = CLUSTERS_TO_PLOT(i);
        if ~isempty(STUDY.cluster(clust_i).sets)
            fprintf('\nPlotting Cluster %s, number: %i\n',c_names{i},clust_i);
            %## Load Scalp Topographies & Average across subjects 
            if ~isfield(STUDY.cluster,'topo') 
                STUDY.cluster(1).topo = [];
            end
            for c = 1:length(clust_i) % For each cluster requested
                cl_i = clust_i(c);
                disp(cl_i)
                if isempty(STUDY.cluster(cl_i).topo)
                    % Using this custom modified code to allow taking average within participant for each cluster
                    STUDY = std_readtopoclust_CL(STUDY,ALLEEG,cl_i);
                end
            end
            %## (TOPOPLOT) 
            figure;
            std_topoplot_CL(STUDY,clust_i,'together');
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'position',TOPO_SING_POS,'color','w');
            drawnow;
            for c = 2:length(fig_i.Children)
                fig_i.Children(c).Title.Interpreter = 'none';
                fig_i.Children(c).TitleFontSizeMultiplier = 1.4;
            end
            saveas(fig_i,fullfile(save_dir,sprintf('Cluster_topo_avg%i.jpg',clust_i)));
            %## (DIPOLE) Plot dipole clusters 
            STUDY.etc.dipparams.centrline = 'off';
            std_dipplot_CL(STUDY,ALLEEG,'clusters',clust_i,'figure','off','mode','together_averaged_only','spheres','off','projlines','off');
            fig_i = get(groot,'CurrentFigure');
            set(fig_i,'position',DIP_SING_POS,'color','w')
            set(fig_i, 'DefaultAxesTickLabelInterpreter', 'none')
            camzoom(1.2^2);
            saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top%i.jpg',clust_i)]);
            saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_top%i.fig',clust_i)]);
            view([45,0,0])
            saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_coronal%i.jpg',clust_i)]);
            view([0,-45,0])
            saveas(fig_i,[save_dir filesep sprintf('dipplot_avgdipspc_sagittal%i.jpg',clust_i)]);
            %## (DIPOLE) Plot dipole fit locations after averaged within participants and
            % different clusters have different colors
            STUDY.etc.dipparams.centrline = 'off';
            std_dipplot_CL(STUDY,ALLEEG,'clusters',clust_i,'figure','off','mode','together_averaged_multicolor','spheres','off','projlines','off');
            fig_i = get(groot,'CurrentFigure');
            set(fig_i, 'DefaultAxesTickLabelInterpreter', 'none')
            set(fig_i,'position',DIP_SING_POS,'color','w')
            camzoom(1.2^2);
            saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top%i.jpg',clust_i)]);
            saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_top%i.fig',clust_i)]);
            view([45,0,0])
            saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_coronal%i.jpg',clust_i)]);
            view([0,-45,0])
            saveas(fig_i,[save_dir filesep sprintf('dipplot_alldipspc_sagittal%i.jpg',clust_i)]);
        else
            fprintf('\nNo Subjects in Cluster %s, number: %i\n',c_names{cllust_i},cllust_i);
        end
    end
    %}
end
end

