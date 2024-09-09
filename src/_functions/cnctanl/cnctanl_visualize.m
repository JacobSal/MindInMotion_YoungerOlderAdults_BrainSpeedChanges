function [cfg,handles] = cnctanl_visualize(EEG,components,clusters,visual,varargin)
%VISUALCNCTANL Summary of this function goes here
%   Detailed explanation goes here

%## TIME
tic
%## DEFINE DEFAULTS
CLIM = [];
Defaults = {'nogui',...
    [],...
    {},...
    {}};
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'components',@isnumeric)
addRequired(p,'clusters',@iscell)
addRequired(p,'visual',@ischar)

%## OPTIONAL

%## PARAMETER
verify_path = (@(x) isempty(x) || isnumeric(x));
addParameter(p,'GUI_MODE',Defaults{1}, @ischar);
addParameter(p,'savePath',Defaults{2}, @ischar);
addParameter(p,'colorLimit',Defaults{2},verify_path);
addParameter(p,'connMeasures',Defaults{4},@iscell);
addParameter(p,'componentNames',Defaults{2});
addParameter(p,'CLIM',CLIM);
parse(p, EEG, components, clusters, visual, varargin{:});
%## SET DEFAULTS
GUI_MODE = p.Results.GUI_MODE;
savePath = p.Results.savePath;
colorLimit = p.Results.colorLimit;
connMeasures = p.Results.connMeasures;
componentNames = p.Results.componentNames;
CLIM = p.Results.CLIM;
%% ===================================================================== %%
if isempty(clusters)
    clusters = cell(1,length(components));
end

if isempty(componentNames)
    ComponentNames = cell(1,length(components));
    for j = 1:length(components)
        ComponentNames{j} =  sprintf('cl%s_co%i',clusters{j},components(j));
    end
else
    ComponentNames = componentNames;
end

if ~exist(savePath,'dir')
    mkdir(savePath)
end
if isempty(colorLimit)
    colorLimit = 95;
end
switch visual
    case 'TimeFreqChart'
        %% STEP 9: Visualize the Connectivity estimates in a Time-Frequency Grid
        
        %% ==== PARAMS ==== %%
        % THRESHOLDING
        % BASELINE
        %% ==== END: PARAMS ==== %%
        for cnd = 1:length(connMeasures)
            cfg = [];
            cfg.plotCondDiff = false;
            cfg.vismode = 'TimeXFrequency';
            cfg.msubset = 'all';
            cfg.MatrixLayout = [];
                cfg.MatrixLayout.arg_direct = 0;
                cfg.MatrixLayout.estimator = connMeasures{cnd};
                cfg.MatrixLayout.clim = colorLimit;
                cfg.MatrixLayout.arg_selection  = 'Full';
            cfg.clim = colorLimit; % percentage of data you want to capture in plots % unsure what this do. val = max(max(squeeze(prctile(abs(EEG.CAT.Conn.mCoh),100))))
            cfg.timeRange = [];
            cfg.freqValues = [1:60];
            cfg.windows = [];
            cfg.pcontour = [];
                cfg.pcontour.arg_direct = 0;
                cfg.pcontour.arg_selection = 0;
            cfg.baseline = [];
            cfg.fighandles  = [];
            cfg.smooth = 0;
            cfg.xord = [];
            cfg.yord = [];
            cfg.xloc = [];
            cfg.yloc = [];
            cfg.plotorder = [];
            cfg.topoplot = 'topoplot';
            cfg.nodelabels = ComponentNames;
            cfg.foilines = [7,15,30];
            cfg.foilinecolor = [0.7,0.7,0.7];
            cfg.events = {{0 'r' ':' 2}}; % (CELLs of CELLs), format: {position, color, line-type, line-thickness}
            cfg.freqscale = 'linear';
            cfg.transform = 'linear';
            cfg.yTickLoc = 'right';
            cfg.titleString = '';
            cfg.titleFontSize = 12;
            cfg.axesFontSize = 9;
            cfg.textColor = [1,1,1];
            cfg.linecolor = [1,1,1];
            cfg.patchcolor = [1,1,1];
            cfg.colormap = 'jet(300)';
            cfg.backgroundColor = [0,0,0];
            cfg.colorscheme = 'white';
            hold on;
            handles = feval(@vis_TimeFreqGrid,'ALLEEG',EEG,'Conn',EEG.CAT.Conn,cfg);
            fig = gcf;
            fig.Position = [0.05 0.1 0.9 0.8];       
            hold off;
            saveas(gcf,[savePath filesep sprintf('TimeFreqChart_%s.fig',connMeasures{cnd})]);
            saveas(gcf,[savePath filesep sprintf('TimeFreqChart_%s.png',connMeasures{cnd})]);
        end
        % This example plots a Time-Frequency Grid using "simple" percentile
        % statistics (this doesn't use the rigorous stats returned by
        % stat_surrogateStats).

        % For a single condition, we call pop_vis_TimeFreqGrid(EEG(cnd),...)
        % If we want to compare two conditions we can either use the dataset
        % returned by pop_stat_surrogateStat() with the 'Hab' statistics mode 
        % OR we can compare set1-set2 by calling 
        % pop_vis_TimeFreqGrid(EEG([set1 set2]), ... ) where set1,set2 are the 
        % indices of the datasets we want to compare.j
        
        %SEE. vis_TimeFreqGrid.m
%         [EEG,cfg,handles] = pop_vis_TimeFreqGrid(EEG,typeproc,cfg);
        
    case 'exTimeFreqChart'
        % MRI
%         path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
        mniMRI = fullfile(path2BEM, 'standard_mri.mat');
        % THRESHOLDING
        THRESH_CHAR = 'Statistics';
        % BASELINE
        %% ==== END: PARAMS ==== %%
        for cnd = 1:length(connMeasures)
            %% LOOP
            cfg = [];
            cfg.plotCondDiff = false;
            cfg.vismode = 'TimeXFrequency';
            cfg.msubset = 'all';
            cfg.thresholding = [];
            cfg.Stats = EEG.CAT.Stats;
            if ~strcmp(connMeasures{cnd},'mCoh')
                cfg.thresholding.arg_direct = 0;
                cfg.thresholding.arg_selection = 'Statistics';
                cfg.thresholding.plotci = false;
                cfg.thresholding.sigthreshmethod = 'pval';
                cfg.thresholding.alpha = 0.01;
                cfg.MatrixLayout = [];
                    cfg.MatrixLayout.arg_direct = 0;
                    cfg.MatrixLayout.estimator = connMeasures{cnd};
                    cfg.MatrixLayout.clim = colorLimit;
                    cfg.MatrixLayout.arg_selection  = 'Full';
            else
                continue;
            end
            cfg.clim = colorLimit; % percentage of data you want to capture in plots % unsure what this do. val = max(max(squeeze(prctile(abs(EEG.CAT.Conn.mCoh),100))))
            cfg.timeRange = [];
            cfg.freqValues = [1:100];
            cfg.windows = [];           
            cfg.baseline = [];
            cfg.topoplot = 'topoplot';
            cfg.fighandles  = [];
            cfg.smooth = 0;
            cfg.xord = [];
            cfg.yord = [];
            cfg.xloc = [];
            cfg.yloc = [];
            cfg.plotorder = [];
            cfg.nodelabels = ComponentNames;
            cfg.foilines = [7,15,30];
            cfg.foilinecolor = [0.7,0.7,0.7];
            cfg.events = {{0 'r' ':' 2}}; % (CELLs of CELLs), format: {position, color, line-type, line-thickness}
            cfg.freqscale = 'log';
            cfg.transform = 'linear';
            cfg.yTickLoc = 'right';
            cfg.titleString = sprintf('Subject %s',EEG.subject);
            cfg.titleFontSize = 10;
            cfg.axesFontSize = 9;
            cfg.textColor = [1,1,1];
            cfg.linecolor = [1,1,1];
            cfg.patchcolor = [1,1,1];
            cfg.colormap = 'jet(300)';
            cfg.backgroundColor = [0,0,0];
            cfg.colorscheme = 'white';
            hold on;
            handles = feval(@vis_TimeFreqGrid,'ALLEEG',EEG,'Conn',EEG.CAT.Conn,cfg);
%             handles = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,cfg);
%             fig = gcf;
%             fig.Position = [0.025 0.05 0.95 0.85];       
            hold off;
            saveas(gcf,[savePath filesep sprintf('exTimeFreqChart_%s.fig',connMeasures{cnd})]);
            saveas(gcf,[savePath filesep sprintf('exTimeFreqChart_%s.png',connMeasures{cnd})])
        end
    case 'BrainMovie'
        %% STEP 10: Visualize the Connectivity estimates in a 3D Brain-Movie
        pop_vis_causalBrainMovie3D(EEG,GUI_MODE,'stats',[],'connmethod','dDTF08','timeRange',[] ,'freqsToCollapse',[1:40] ,'collapsefun','max','resample',0,'subtractconds',false,'showNodeLabels',{'nodelabels' ComponentNames'},'nodesToExclude',{},'edgeColorMapping','PeakFreq','edgeSizeMapping','Connectivity','nodeColorMapping','Outflow','nodeSizeMapping','Power','baseline',[],'normalize',true,'useStats',[],'prcthresh',0.05,'absthresh',[],'footerPanelSpec',{'ICA_ERPenvelope' 'plottingmode' {'all' 'envelope'} 'envColor' [1 0 0] },'BMopts',{'size' [800 800]  'visible' 'on' 'latency' [] 'frames' [] 'figurehandle' [] 'cameraMenu' false 'rotationpath3d' {'AngleFactor' 1 'PhaseFactor' 0.75 'FramesPerCycle' []} 'view' [122 36]  'makeCompass' true 'project3d' 'off' 'theme' {'theme' 'classic'} 'Layers' {'scalp' {'scalpres' 'high' 'volumefile' [] 'scalptrans' 0.9 'scalpcolor' [1 0.75 0.65] } 'skull' [] 'csf' [] 'cortex' {'cortexres' 'mid' 'volumefile' [] 'cortextrans' 0.9 'cortexcolor' {'LONI_Atlas' 'colormapping' {'jet'}}} 'custom' []} 'facelighting' 'phong' 'opengl' 'on' 'events' {{0 'r' ':' 2}} 'flashEvents' false 'square' 'on' 'caption' true 'displayLegendLimitText' true 'showLatency' true 'dispRT' false 'backcolor' [0 0 0]  'graphColorAndScaling' {'nodeSizeLimits' [0.1 1]  'nodeColorLimits' [0 1]  'edgeSizeLimits' [0.1 0.8]  'edgeColorLimits' [0 1]  'nodeSizeDataRange' [] 'nodeColorDataRange' [] 'edgeSizeDataRange' [] 'edgeColorDataRange' [] 'centerDataRange' false 'edgeColormap' jet(64) 'diskscale' 0.2 'magnify' 1} 'outputFormat' {'framefolder' '' 'framesout' 'jpg' 'moviename' '' 'movieopts' {'videoname' ''} 'size' []} 'mri' 'standard_BEM_mri.mat' 'plotimgs' false 'coordformat' 'spherical' 'dipplotopt' {} 'bmopts_suppl' {} 'renderBrainMovie' true 'speedy' false 'mode' 'init_and_render' 'vars' []});
    case 'IndvCell_TxF'
        %## TIMEFREQ AVG PLOT
        for cnd = 1:length(connMeasures)
            if strcmp(connMeasures{cnd},'mCoh')
                continue;
            end
            clim = CLIM;
            done = [];
            for j = 1:length(ComponentNames)
                for i = 1:length(ComponentNames)
                    if j == i || any((i == done))
                        continue;
                    end
                    fprintf('%s) Plotting Time Frequency Plot...\n',EEG.subject)
                    %## Prepare the arguments for vis_TimeFreqCell()
                    %- 
                    subargs= [];
                    g.freqValues = EEG.CAT.Conn.freqs;
                    %- 
                    CEstimator = connMeasures{cnd};
                    OrigConnMatrix = EEG.CAT.Conn.(connMeasures{cnd});
                    ConnMatrix = EEG.CAT.Conn.(connMeasures{cnd});
                    erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes;
                    origFreqValues = EEG.CAT.Conn.freqs;
                    %- extract the stats matrix for this pair 
                    if isfield(EEG.CAT, 'Stats')
                        if ~isempty(EEG.CAT.Stats.(connMeasures{cnd}))
                            subargs.thresholding = [];
                            subargs.thresholding.arg_selection = 'Statistics';
                            subargs.thresholding.sigthreshmethod = 'pval';
                            subargs.thresholding.alpha.alpha = 0.05;
    %                         subargs.thresholding.alpha = 0.01;
                            subargs.thresholding.plotci = [];
                        else
                            subargs.thresholding    = 'None';
                        end             
                        StatsMatrix = EEG.CAT.Stats.(connMeasures{cnd}).(subargs.thresholding.sigthreshmethod);
                        %- Stats extract
                        if isequal(size(StatsMatrix),size(ConnMatrix))
                            Sji = squeeze(StatsMatrix(i,j,:,:));
                            Sij = squeeze(StatsMatrix(j,i,:,:));
                        elseif size(StatsMatrix,1)==2 && ndims(StatsMatrix)==5
                            Sji = permute(squeeze(squeeze(StatsMatrix(:,i,j,:,:))),[2 3 1]);
                            Sij = permute(squeeze(squeeze(StatsMatrix(:,j,i,:,:))),[2 3 1]);
                        else
                            Sji = StatsMatrix;
                            Sij = StatsMatrix;
                        end
                        %- statistics definition if available
                        if ~isempty(Sji)
                            subargs.StatsMatrix(1,:,:,:,:) = Sji;
                            subargs.StatsMatrix(2,:,:,:,:) = Sij;
                        end
                    else
                        subargs.StatsMatrix = [];
                        subargs.thresholding = 'None';
                    end
                    %- subargs definintion
                    subargs.elocs       = EEG.chanlocs;
                    subargs.chaninfo    = EEG.chaninfo;
                    subargs.alltimes    = erWinCenterTimes;
                    subargs.allfreqs    = origFreqValues;
                    %- connectivity extraction
                    subargs.ConnMatrix(1,:,:)  = squeeze(OrigConnMatrix(i,j,:,:));
                    subargs.ConnMatrix(2,:,:)  = squeeze(OrigConnMatrix(j,i,:,:));
                    %- color limits handle
                    tmp = mean(squeeze(subargs.ConnMatrix(1,:,:)),1);
                    tmp = subargs.ConnMatrix(:);
                    rnge = range(tmp);
                    mn = mean(tmp,1);
                    stdv = std(tmp);
                    med = median(tmp);
                    prc = prctile(tmp,[1,25,75,99]); 
                    if isempty(clim)
                        clim = [mn-(3*stdv),mn+(3*stdv)];
                    end
                    fprintf('Range: %0.2g\n',rnge);
                    fprintf('Mean: %0.2g\n',mn);
                    fprintf('Standard Deviation: %0.2g\n',stdv);
                    fprintf('Median: %0.2g\n',med);
                    fprintf('1st, 25th, 75th, and 99th Percentile: %0.2g, %0.2g, %0.2g, %0.3f\n',prc(1),prc(2),prc(3),prc(4));
                    fprintf('Min & Max: %0.2g & %0.2g\n',min(tmp),max(tmp));
                    %- topoplot plotting
                    subargs.topoplot    = 'topoplot';                    
                    if strcmpi(subargs.topoplot,'topoplot')
                        subargs.topovec     = squeeze(EEG.icawinv(:,EEG.CAT.curComps([j i])))';
                    elseif strcmpi(subargs.topoplot,'customtopo')
                        subargs.customTopoMatrix = g.customTopoMatrix([j i]);
                    else
                        subargs.topovec = [];
                        subargs.customTopoMatrix = {};
                    end
                    %- topoplot options
                    subargs.topoplot_opts = {};
                    %- plot formatting
                    subargs.baseline    = [];
                    subargs.freqscale   = 'log';
                    subargs.events      = {{0 'r' ':' 2}}; 
                    subargs.linewidth = 3;
                    subargs.titleString = sprintf('%s) IC %s->%s',EEG.subject,ComponentNames{i},ComponentNames{j});
                    subargs.titleFontSize   = 12;
                    subargs.axesFontSize    = 10;
                    subargs.textColor       = [1 1 1];
                    subargs.backgroundColor = [0 0 0];
                    subargs.clim            = clim;                    
                    subargs.bidir           = fastif(i==j,false,true);
                    subargs.connmethod      = CEstimator;
                    subargs.nodelabels      = ComponentNames([j i]);
                    subargs.dipplot         = [];
                    subargs.foilines        = [7 15 30];
                    subargs.foilinecolor    = [0.7 0.7 0.7];
                    subargs.smooth          = [];
                    subargs.colorscheme     = 'white';
                    subargs.colormap        = 'jet(300)';
                    out = feval(@vis_TimeFreqCell, subargs.ConnMatrix,...
                                subargs.alltimes ,subargs.allfreqs,...
                                subargs.StatsMatrix, [], subargs.elocs,...
                                subargs.chaninfo,subargs.connmethod,...
                                subargs.nodelabels,[],[],...
                                subargs.baseline,subargs.smooth,subargs.freqscale,subargs.linewidth,...
                                subargs.events,subargs.colormap,subargs.bidir,subargs.topoplot,...
                                subargs.topoplot_opts,{},subargs.topovec,...
                                subargs.dipplot,subargs.titleString,...
                                subargs.titleFontSize,subargs.axesFontSize,...
                                subargs.textColor,subargs.backgroundColor,...
                                subargs.colorscheme,subargs.foilines,...
                                subargs.foilinecolor,subargs.clim,...
                                subargs.thresholding);
                    fprintf('%s) Done.\n',EEG.subject)
                    fig = get(groot,'Currentfigure');
%                     fig.Position = [0.025 0.05 0.95 0.85];
                    if ~exist(savePath,'dir')
                        mkdir(savePath)
                    end
                    saveas(fig,[savePath filesep sprintf('cell_%s-%s_%s.fig',ComponentNames{i},ComponentNames{j},connMeasures{cnd})]);
                    saveas(fig,[savePath filesep sprintf('cell_%s-%s_%s.png',ComponentNames{i},ComponentNames{j},connMeasures{cnd})]);
                end
                done = [done j];
            end
        end
end

end

