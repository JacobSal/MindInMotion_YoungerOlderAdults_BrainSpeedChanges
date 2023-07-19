function [h] = cnctanl_vis_timefreq(EEG,comp_indxs,varargin)
%CNCTANL_VIS_TIMEFREQ Summary of this function goes here
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

%## TIME
tic
%## DEFINE DEFAULTS
%*
DISPLAY_NAMES = EEG.CAT.curComponentNames;
errorMsg = 'Value must be a CELL of CHARS (e.g., {''I'',''love'',''EEG''}). This will provide a name for the i''th component'; 
dn_validFcn = @(x) assert(iscell(x) && (length(x) == length(EEG.CAT.curComponentNames)),errorMsg);
%*
SAVE_DIR = [];
errorMsg = 'Value must be a CHAR. Path where figures will be saved.'; 
sd_validFcn = @(x) assert(ischar(x),errorMsg);
%*
CONN_MEASURES = EEG.CAT.configs.est_mvarConnectivity.connmethods;
errorMsg = 'Value must be a CELL of CHARS (e.g., {''I'',''love'',''EEG''}). Path where figures will be saved.'; 
cn_validFcn = @(x) assert(iscell(x),errorMsg);
%*
COLOR_LIM = [0,0.01];
errorMsg = 'Value must be formatted as [INT_min,INT_max]. Limits the the color scaling for plotting to.'; 
cl_validFcn = @(x) assert(isnumeric(x) && (length(x) == 2),errorMsg);
%## PARSE 
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'comp_indxs',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'DISPLAY_NAMES',DISPLAY_NAMES,dn_validFcn)
addParameter(p,'SAVE_DIR',SAVE_DIR,sd_validFcn);
addParameter(p,'CONN_MEASURES',CONN_MEASURES,cn_validFcn)
addParameter(p,'COLOR_LIM',COLOR_LIM,cl_validFcn)
parse(p, EEG, comp_indxs, varargin{:});
%## SET DEFAULTS
display_names = p.Results.DISPLAY_NAMES;
save_dir = p.Results.SAVE_DIR;
conn_measures = p.Results.CONN_MEASURES;
color_lim = p.Results.COLOR_LIM;
%% ===================================================================== %%
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
%## TIMEFREQ AVG PLOT
for conn_i = 1:length(conn_measures)
    fprintf('==== %s Plotting Time-Frequency Connectivity of %s ====\n',EEG.subject,conn_measures{conn_i})
    %- ignore measures that only perform within component coherencey
    if strcmp(conn_measures{conn_i},'mCoh')
        continue;
    end
    %- loop through unique connections
    done = [];
    for i = comp_indxs
        for j = comp_indxs
            if i == j || any((j == done))
                continue;
            end
            conn_meas = conn_measures{conn_i};
            h = mim_vis_TimeFreqCell(EEG,conn_meas,display_names,color_lim,j,i);
            fprintf('%s) Done.\n',EEG.subject)
            if ~isempty(save_dir)
                if ~exist(save_dir,'dir')
                    mkdir(save_dir)
                end
                fig = get(groot,'Currentfigure');
                saveas(fig,[save_dir filesep sprintf('cell_%s-%s_%s.fig',display_names{j},display_names{i},conn_measures{conn_i})]);
                saveas(fig,[save_dir filesep sprintf('cell_%s-%s_%s.png',display_names{j},display_names{i},conn_measures{conn_i})]);
            end
        end
        done = [done i];
    end
end
end
%% SUBFUNCTIONS
function [fig] = mim_vis_TimeFreqCell(EEG,conn_meas,comp_names,color_lim,i,j)
    %## Prepare the arguments for vis_TimeFreqCell()
    subargs= [];
    g.freqValues = EEG.CAT.Conn.freqs;
    %- 
    OrigConnMatrix = EEG.CAT.Conn.(conn_meas);
    ConnMatrix = EEG.CAT.Conn.(conn_meas);
    erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes;
    origFreqValues = EEG.CAT.Conn.freqs;
    %-
    if isfield(EEG.CAT,'masked_conn')
        connmethods = EEG.CAT.configs.est_mvarConnectivity.connmethods;
        for conn_i = 1:length(connmethods)
            EEG.CAT.Conn.(connmethods{conn_i})=EEG.CAT.masked_conn.(connmethods{conn_i});
        end
    end
    %- extract the stats matrix for this pair 
    if isfield(EEG.CAT, 'Stats')
        if ~isempty(EEG.CAT.Stats.(conn_meas))
            subargs.thresholding = [];
            subargs.thresholding.arg_selection = 'Statistics';
            subargs.thresholding.sigthreshmethod = 'pval'; %'ci';
%             subargs.thresholding.sigthreshmethod = 'ci';
            subargs.thresholding.alpha.alpha = 0.01;
%                         subargs.thresholding.alpha = 0.01;
            subargs.thresholding.plotci = [];
        else
            subargs.thresholding    = 'None';
        end             
        StatsMatrix = EEG.CAT.Stats.(conn_meas).(subargs.thresholding.sigthreshmethod);
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
    subargs.baseline        = [];
    subargs.freqscale       = 'log';
    subargs.events          = {{0 'r' ':' 2}}; 
    subargs.linewidth       = 3;
    subargs.titleString     = sprintf('%s) IC %s->%s',EEG.subject,comp_names{i},comp_names{j});
    subargs.titleFontSize   = 12;
    subargs.axesFontSize    = 10;
    subargs.textColor       = [1 1 1];
    subargs.backgroundColor = [0 0 0];
    subargs.clim            = color_lim;                    
    subargs.bidir           = fastif(i==j,false,true);
    subargs.connmethod      = conn_meas;
    subargs.nodelabels      = comp_names([j i]);
    subargs.dipplot         = [];
    subargs.foilines        = [7 15 30];
    subargs.foilinecolor    = [0.7 0.7 0.7];
    subargs.smooth          = [];
    subargs.colorscheme     = 'white';
    subargs.colormap        = 'jet(300)';
    fig = feval(@vis_TimeFreqCell, subargs.ConnMatrix,...
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
end

