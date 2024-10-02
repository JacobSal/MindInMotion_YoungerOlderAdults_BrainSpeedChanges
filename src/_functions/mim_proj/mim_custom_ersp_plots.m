function [] = mim_custom_ersp_plots(STUDY,cond_test,warping_times,cluster_ind,...
    cluster_load_ind,des_i,save_dir,varargin)
%MIM_GEN_ERSP_STATS Summary of this function goes here
% Description
% This MATLAB function generates and saves various Event-Related Spectral Perturbation (ERSP) plots for a given study using parallel processing. It supports multiple types of baseline correction and saves the output as MATLAB mat files.
% 
% Usage
% matlab
% Copy code
% mim_custom_ersp_plots(STUDY, cond_test, warping_times, cluster_ind, cluster_load_ind, des_i, save_dir, varargin)
% Inputs
% STUDY : Struct - EEG study structure containing information about the dataset and experiment.
% cond_test : Cell array - Cell array of conditions or tests to be analyzed.
% warping_times : Numeric array - Times used for warping ERSP data.
% cluster_ind : Numeric - Index representing the cluster to be processed.
% cluster_load_ind : Numeric - Index for loading the cluster-specific ERSP data.
% des_i : Numeric - Design index for the study.
% save_dir : String - Directory path where the generated plots will be saved.
% Optional Inputs
% 'ALLERSP' : Cell array (default: empty) - Preloaded ERSP data for all subjects.
% 'ALLTIMES' : Numeric array (default: empty) - Preloaded time values for ERSP data.
% 'ALLFREQS' : Numeric array (default: empty) - Preloaded frequency values for ERSP data.
% 'PCOND' : Cell array (default: empty) - Precomputed conditions for ERSP data.
% 'PGROUP' : Cell array (default: empty) - Precomputed groups for ERSP data.
% 'PINTER' : Cell array (default: empty) - Precomputed interactions for ERSP data.
% 'DO_SUBJ_PLOTS' : Logical (default: false) - Flag to enable subject-specific plots.
% 'CLUSTER_CLIM_MATCH' : Numeric array (default: empty) - Cluster color limits for matching.
% Outputs
% None. Generates and saves ERSP plots based on the specified parameters.
% Examples
% matlab
% Copy code
% % Example 1: Generate ERSP plots without subject-specific plots
% mim_custom_ersp_plots(STUDY, {'Condition1', 'Condition2'}, warping_times, 1, 2, 3, '/path/to/save');
% 
% % Example 2: Generate ERSP plots with subject-specific plots
% mim_custom_ersp_plots(STUDY, {'Condition1', 'Condition2'}, warping_times, 1, 2, 3, '/path/to/save', 'DO_SUBJ_PLOTS', true);

% Date
% Code Date: 07/17/2023
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
ALLERSP = [];
ALLTIMES = []; 
ALLFREQS = [];
PCOND = [];
PGROUP = [];
PINTER = [];
DO_SUBJ_PLOTS = false;
DO_BASELINE_CORRECT_1 = true;
DO_BASELINE_CORRECT_2 = true;
DO_BASELINE_CORRECT_3 = true;
% DO_BASELINE_CORRECT_4 = false;
ERSP_ALPHA = 0.05;
CLUSTER_CLIM_MATCH = [];
%- custom params
SUB_FREQ_LIMS = [4,60];
colormap_ersp = linspecer; %othercolor('RdYlBu11');
% colormap_ersp = colormap_ersp(end:-1:1,:);
colormap(colormap_ersp);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'cond_test',@iscell);
addRequired(p,'warping_times',@isnumeric);
addRequired(p,'cluster_ind',@isnumeric);
addRequired(p,'cluster_load_ind',@isnumeric);
addRequired(p,'des_i',@isnumeric);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'ALLERSP',ALLERSP,@iscell);
addParameter(p,'ALLTIMES',ALLTIMES,@isnumeric);
addParameter(p,'ALLFREQS',ALLFREQS,@isnumeric);
addParameter(p,'PCOND',PCOND,@iscell);
addParameter(p,'PGROUP',PGROUP,@iscell);
addParameter(p,'PINTER',PINTER,@iscell);
addParameter(p,'DO_SUBJ_PLOTS',DO_SUBJ_PLOTS,@islogical);
addParameter(p,'CLUSTER_CLIM_MATCH',CLUSTER_CLIM_MATCH,@(x) isempty(x)||isnumeric(x))
parse(p,STUDY,cond_test,...
    warping_times,cluster_ind,cluster_load_ind,des_i,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
DO_SUBJ_PLOTS = p.Results.DO_SUBJ_PLOTS;
CLUSTER_CLIM_MATCH = p.Results.CLUSTER_CLIM_MATCH;
allersp = p.Results.ALLERSP;
alltimes = p.Results.ALLTIMES;
allfreqs = p.Results.ALLFREQS;
pcond = p.Results.PCOND;
pgroup = p.Results.PGROUP;
pinter = p.Results.PINTER;
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
chk_in_1 = isempty(allersp) && isempty(alltimes) && isempty(allfreqs) && isempty(pcond) && isempty(pgroup) && isempty(pinter);
chk_in_2 = isempty(allersp) && isempty(alltimes) && isempty(allfreqs) && (isempty(pcond) || isempty(pgroup) || isempty(pinter));
if chk_in_1
    warning('Using default loadings for allersp, alltimes, allfreqs, and stats for each plot type')
    inside_load_flag = true;
elseif chk_in_2
    msg = [sprintf('allersp:\n'),evalc('disp(allersp)'),sprintf('alltimes:\n'),evalc('disp(alltimes)'),...
        sprintf('allfreqs:\n'),evalc('disp(allfreqs)'),sprintf('pcond:\n'),evalc('disp(pcond)'),...
        sprintf('pgroup:\n'),evalc('disp(pgroup)'),sprintf('pinter:\n'),evalc('disp(pinter)')];
    errID = 'mim_custom_ersp_plots:ImproperInputs';
    baseException = MException(errID,msg);
    throw(baseException);
else
    warning('Using inputs of allersp, alltimes, allfreqs, and stats for each plot type')
    inside_load_flag = false;
end
close all;
%## subject plots
if DO_SUBJ_PLOTS
    % plot each component with baseline to be the average of entire epoch
    %- raw ersp data no baselining
    if inside_load_flag
%         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_fpaths,'/');
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).ersp_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data = par_load(fpath,fname);
        allersp = ersp_data.allerspdata;
        alltimes = ersp_data.alltimes;
        allfreqs = ersp_data.allfreqs;
    end
    ss_save_dir = [save_dir filesep sprintf('%i',cluster_ind) filesep sprintf('%i',des_i)];
    if ~exist(ss_save_dir,'dir')
        mkdir(ss_save_dir)
    end
    fprintf('Using ERSP data:\n');
    disp(allersp)
%     cl_in = cluster_i + 1; % offset for parent cluster
    ersp_single_subj_plot(STUDY,allersp,alltimes,allfreqs,...
        warping_times,SUB_FREQ_LIMS,cluster_ind,colormap_ersp,ss_save_dir);
    
end
%% 1. Baseline correction = average of epoch within condition
if DO_BASELINE_CORRECT_1
    %- raw ersp data no baselining
    if inside_load_flag
%         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_fpaths,'/');
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).ersp_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data = par_load(fpath,fname);
        allersp = ersp_data.allerspdata;
        alltimes = ersp_data.alltimes;
        allfreqs = ersp_data.allfreqs;
    end
    ss_save_dir = [save_dir filesep sprintf('%i',cluster_ind) filesep sprintf('%i',des_i)];
    if ~exist(ss_save_dir,'dir')
        mkdir(ss_save_dir)
    end
    fprintf('Using ERSP data:\n');
    disp(allersp)
    ersp_baseline_plot_1(STUDY,allersp,alltimes,allfreqs,...
        warping_times,colormap_ersp,SUB_FREQ_LIMS,cluster_ind,ERSP_ALPHA,ss_save_dir);
end
%% 2. Common Baseline correction = average of epoch across condition (commonbaseline)
% Read allerspdata3, allerspdata3 should already be baselined
% allerspdata,  
if DO_BASELINE_CORRECT_2
    %- normalized ersp data
    if inside_load_flag
%         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_norm_fpaths,'/');
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).ersp_norm_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data_norm = par_load(fpath,fname);
        allersp = ersp_data_norm.allerspdata;
        alltimes = ersp_data_norm.alltimes;
        allfreqs = ersp_data_norm.allfreqs;
        pcond = ersp_data_norm.pcond;
    end
    ss_save_dir = [save_dir filesep sprintf('%i',cluster_ind) filesep sprintf('%i',des_i)];
    if ~exist(ss_save_dir,'dir')
        mkdir(ss_save_dir)
    end
%     [allersp,~,~] = mim_prune_ersp_trials(STUDY,allersp,cluster_i,des_i);
    
    fprintf('Using ERSP data:\n');
    disp(allersp)
%     fprintf('Using ERSP stats:\n');
%     disp(pcond)
    ersp_baseline_plots(STUDY,2,allersp,allfreqs,alltimes,pcond,...
        warping_times,colormap_ersp,SUB_FREQ_LIMS,cluster_ind,ss_save_dir);
%     ersp_baseline_plot_2(STUDY,allersp,allfreqs,alltimes,pcond,...
%         cond_test,warping_times,colormap_ersp,SUB_FREQ_LIMS,des_i,cluster_ind,ERSP_ALPHA,CLUSTER_CLIM_MATCH,ss_save_dir);
end
%% 3. Single trial full epoch correction + Common Baseline correction = average of epoch across condition (commonbaseline)
% Read allerspdata4, allerspdata4 should already be baselined
if DO_BASELINE_CORRECT_3
    if inside_load_flag
%         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_normcb_fpaths,'/');
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).ersp_normcb_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data_normcb = par_load(fpath,fname);
        allersp = ersp_data_normcb.allerspdata;
        alltimes = ersp_data_normcb.alltimes;
        allfreqs = ersp_data_normcb.allfreqs;
        pcond = ersp_data_normcb.pcond;
    end
    ss_save_dir = [save_dir filesep sprintf('%i',cluster_ind) filesep sprintf('%i',des_i)];
    if ~exist(ss_save_dir,'dir')
        mkdir(ss_save_dir)
    end
%     [allersp,~,~] = mim_prune_ersp_trials(STUDY,allersp,cluster_i,des_i)
    fprintf('Using ERSP data:\n');
    disp(allersp)
%     if inside_load_flag
%         fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_norm_fpaths,'/');
%         fpath = strjoin(fname(1:end-1),'/');
%         fname = fname{end};
%         ersp_data_norm = par_load(fpath,fname);
%         allersp = ersp_data_norm.allerspdata;
%         alltimes = ersp_data_norm.alltimes;
%         allfreqs = ersp_data_norm.allfreqs;
%         pcond = ersp_data_norm.pcond;
%     end
    ersp_baseline_plots(STUDY,3,allersp,allfreqs,alltimes,pcond,...
        warping_times,colormap_ersp,SUB_FREQ_LIMS,cluster_ind,ss_save_dir);
%     ersp_baseline_plot_3(STUDY,allersp,allfreqs,alltimes,...
%         pcond,warping_times,colormap_ersp,SUB_FREQ_LIMS,des_i,cluster_ind,ERSP_ALPHA,CLUSTER_CLIM_MATCH,ss_save_dir)
end
%% 4. 
% if DO_BASELINE_CORRECT_4
%     ersp_baseline_plot_4()
% end

end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [] = ersp_single_subj_plot(STUDY,allersp,alltimes,allfreqs,...
        warping_times,sub_freq_lims,cluster_i,colormap_ersp,save_dir)
    FREQ_BANDS = {4:8,9:11,8:13,13:30,30:60};
    FREQ_BANDS_CHARS = {'theta','mu','alpha','beta','gamma'};
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    COLOR_LIM_INTERVALS = [0.6,1.2,1.5,2,2.5];
    COLOR_LIM_ERR = 0.05;
    FIGURE_POSITION = [100,100,1480,350];
    PANEL_OFFSET = [-0.04,-0.03,0.02,-0.075]; % [LEFT,BOTTOM,WIDTH,HEIGHT];
    %- set titles
    if any(strcmp(STUDY.design(STUDY.currentdesign).variable(1).value,SPEED_REF_CHAR))
        condnames = SPEED_OVERRIDE_CHARS;
    else
        condnames = STUDY.design(STUDY.currentdesign).variable(1).value;
    end
    %-
    baseidx = (alltimes>=warping_times(1) & alltimes<=warping_times(end)); 
    freqidx = (allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
    allersp_cond_mean = cell(size(allersp));
    allersp_cond_mean_crop = cell(size(allersp));
    %## Calculate baselines
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            allersp_cond_mean{cond_i,group_i} = mean(allersp{cond_i,group_i},3);
            allersp_cond_mean_crop{cond_i,group_i} = allersp_cond_mean{cond_i,group_i}(freqidx,baseidx);
        end
    end
    %## PLOT
%     alltimes = alltimes(baseidx);
%     allfreqs = allfreqs(freqidx);
    subj_chars = cell(size(allersp{1},3),1);
    ic_val = zeros(size(allersp{1},3),1);
    freq_avgs = zeros(size(allersp{1},3),length(FREQ_BANDS)*size(allersp,1)*size(allersp,2));
    head_char = cell(length(FREQ_BANDS)*size(allersp,1)*size(allersp,2),1);
    for subj_i = 1:size(allersp{1},3)
        %- baseline subject to each condition's mean
        erspdata = cell(size(allersp));
        for group_i = 1:size(allersp,2)
            for cond_i = 1:size(allersp,1)
                %-
%                 erspdata{cond_i,group_i} = allersp{cond_i,group_i}(:,:,subj_i) - allersp_cond_mean{cond_i,group_i};
                %-
%                 erspdata{cond_i,group_i} = erspdata{cond_i,group_i} - allersp_cond_mean{cond_i,group_i};
%                 erspdata{cond_i,group_i} = allersp{cond_i,group_i}(freqidx,baseidx,subj_i);
%                 erspdata = erspdata - allersp_cond_mean_crop{cond_i,group_i};
                %-
                baseline = mean(allersp{cond_i,group_i}(:,baseidx,subj_i),2);
                erspdata{cond_i,group_i} = allersp{cond_i,group_i}(:,:,subj_i)-repmat(baseline,1,length(alltimes));
            end
        end
        %## set color limits
        climMat = [min(erspdata{4}(1:30,:),[],'all') max(erspdata{4}(1:30,:),[],'all')];
        clim_max = [];
        for i = 1:length(COLOR_LIM_INTERVALS)
            chk = any(climMat < (COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR)) && any(climMat > -(COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR));
            if chk
                clim_max = COLOR_LIM_INTERVALS(i);
                break
            end
        end
        if isempty(clim_max)
            clim_max = 1.5;
        end
        ic = STUDY.cluster(cluster_i).comps(subj_i);
        sub = STUDY.datasetinfo(STUDY.cluster(cluster_i).sets(subj_i)).subject;
        cl_ctr = STUDY.cluster(cluster_i).dipole.posxyz;
        all_dips = STUDY.cluster(cluster_i).all_diplocs;
        dists_to_ctr = sqrt(sum((all_dips-cl_ctr).^2,2));
        alltitles = std_figtitle('condnames',condnames, ...
         'clustname', sprintf('%s ic%i cl%i dist%0.1f',sub,ic,cluster_i,dists_to_ctr(subj_i)));
        %- plot time freq
        [fig] = plot_tftopo(erspdata,alltimes,allfreqs,alltitles,[],warping_times,clim_max,colormap_ersp,...
                        PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
        exportgraphics(fig,[save_dir filesep sprintf('%s_within_des%i_cl%i_stats.jpg',sub,STUDY.currentdesign,cluster_i)]);
%         exportgraphics(fig,[save_dir filesep sprintf('%s_within_des%i_cl%i_stats.pdf',sub,STUDY.currentdesign,cluster_i)],'ContentType','vector');
        close(fig)
        %- print averages
%         subj_chars{subj_i} = sub;
%         ic_val(subj_i) = ic;
%         cnt = 1;
%         for group_i = 1:size(allersp,2)
%             for cond_i = 1:size(allersp,1)
%                 for freq_i = 1:length(FREQ_BANDS)
%                     tmpidx = (allfreqs>=min(FREQ_BANDS{freq_i}) & allfreqs<=max(FREQ_BANDS{freq_i}));
%                     freq_avgs(subj_i,cnt) = mean(erspdata{cond_i,group_i}(tmpidx,baseidx),'all');
%                     head_char{cnt} = sprintf('group%i_cond%i_%s',group_i,cond_i,FREQ_BANDS_CHARS{freq_i});
%                     cnt = cnt+1;
%                 end
%             end
%         end
    end
%     T1 = table(subj_chars,'VariableNames',{'subject_char'});
%     T2 = table(ic_val,'VariableNames',{'independent_component'});
%     T3 = array2table(freq_avgs,'VariableNames',head_char);
%     T_out = [T1,T2,T3];
%     writetable(T_out,[save_dir filesep sprintf('freq_averages_des%i_cluster%i.csv',STUDY.currentdesign,cluster_i)])
end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [fig] = plot_tftopo(allersp,alltimes,allfreqs,alltitles,allpcond,warping_times,clim_max,colormap_ersp,...
    PANEL_OFFSET,EVENT_CHARS,SUB_FREQ_LIMS,FIGURE_POSITION)
    XTICK_LABEL = 'Gait Events';
    YTICK_LABEL = 'Frequency (Hz)';
    FONT_SIZE = 12;
    SUBPLOT_WIDTH = 0.15;
    SUBPLOT_HEIGHT = 0.7;
    SHIFT_AMNT = 0.175;
    STATS_TITLE = 'CUSTOM STATS';
    if length(clim_max) == 2
        clim_ersp = clim_max;
    else
        clim_ersp = [-clim_max,clim_max];
    end
    %%
    fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    set(fig,'Units','inches','Position',[3 3 14 5])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    horiz_shift = 0;
    hold on;
    for j = 1:length(allersp)
        subplot(1,length(allersp)+1,j); %,'position',[0.01+horiz_shift,0.1,0.5,0.5])
        ax = gca;
        tftopo(allersp{j},alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan clim_ersp],...
            'logfreq','native');
        hold on;
        colormap(colormap_ersp);
        %- adjust subplot position and height
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
%         disp(get(ax,'Position'));
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
        %- set ylims
        ylim(log(SUB_FREQ_LIMS))
        if SUB_FREQ_LIMS(2) <= 50
            set(ax,'YTick',log([4.01,8,13,30,50])); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end  
        %- set color lims
        set(ax,'clim',clim_ersp);
        %- set x-axis & y-axis labels
        if j == 1
            ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
            xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        else
            xlabel('','FontSize',FONT_SIZE);
            ylabel('','fontsize',FONT_SIZE,'fontweight','bold');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
%         ylim()
        %- title
        title(alltitles{j});  
        horiz_shift = horiz_shift + SHIFT_AMNT;
    end
%     hold off;
    %%
    %## Add Stats To Plot
    if ~isempty(allpcond)
        subplot(1,length(allersp)+1,length(allersp)+1) % add one subplot for stats
        tftopo(double(allpcond),alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan  clim_ersp],...
            'logfreq','native')
        colormap(colormap_ersp);
        ax = gca;
        %-
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
        disp(get(ax,'Position'));
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.04;
        c.Limits = clim_ersp;
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',FONT_SIZE);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        set(hL,'Rotation',0);
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
        %- set ylims
        ylim(log(SUB_FREQ_LIMS))
        if SUB_FREQ_LIMS(2) <= 50
            set(ax,'YTick',log([4.01,8,13,30,50])); 
            set(ax,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(ax,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(ax,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end  
        %- set color lims
        set(ax,'clim',clim_ersp);
        %- set y-axis labels
        xlabel('','FontSize',FONT_SIZE);
        ylabel(sprintf(''),'fontsize',FONT_SIZE,'fontweight','bold');
        %- set x-axis labels
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- title
        title(STATS_TITLE)
    else
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.05;
        c.Limits = clim_ersp;
    end
    hold off;
    fig = get(groot,'CurrentFigure');
end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [fig] = plot_contourf(ersp_raw,ersp_pcond,ersp_masked,alltimes,allfreqs,...
        alltitles,warping_times,clim_max,colormap_ersp,PANEL_OFFSET,EVENT_CHARS,SUB_FREQ_LIMS,FIGURE_POSITION)
    XTICK_LABEL = 'Gait Events';
    YTICK_LABEL = 'Frequency (Hz)';
    FONT_SIZE = 12;
    SUBPLOT_WIDTH = 0.15;
    SUBPLOT_HEIGHT = 0.7;
    SHIFT_AMNT = 0.175;
%     STATS_TITLE = 'CUSTOM STATS';
    COLORBAR_SHIFT = 0.05;
    ALPHA_MULTIPLE = 0.75;
    %%
    fig = figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    set(fig,'Units','inches','Position',[3 3 14 5])
    set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
    horiz_shift = 0;
    hold on;
    for j = 1:length(ersp_raw)
        subplot(1,length(ersp_raw),j)
        alpha_mask = ones(size(ersp_pcond{j}))*ALPHA_MULTIPLE; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        alpha_mask(ersp_pcond{j} == 1) = 0; %0 is significant? 1 is not?
        contourf(alltimes, allfreqs, ersp_raw{j},200,...
                   'linecolor','none');
        hold on;
        imagesc(alltimes,allfreqs,ersp_masked{j},'AlphaData',alpha_mask);
        colormap(colormap_ersp);
        ax = gca;
        %- set figure size and position
        %fig_i.CurrentAxes;
        set(ax,'LineWidth',1)
        set(ax,'FontName','Arial','FontSize',FONT_SIZE,'FontWeight','bold')
        set(ax,'OuterPosition',[0 0 1 1]);
        set(ax,'Position',[0.05+horiz_shift,0.2,SUBPLOT_WIDTH,SUBPLOT_HEIGHT]);  %[left bottom width height]
%         disp(get(ax,'Position'));
        %- add vertical line
        for i = 1:length(warping_times)
            xline(ax,warping_times(i),'k--');
        end
        %- set x-axis ticks
        xrng = get(ax,'XLim');
        if warping_times(1) < xrng(1)
            warping_times(1) = xrng(1);
        end
        if warping_times(end) > xrng(end)
            warping_times(end) = xrng(end);
        end
        set(ax,'XTick',warping_times,'XTickLabel',EVENT_CHARS);
        xtickangle(45)
        ax.XAxis.FontSize = FONT_SIZE;
        %- color lim
        set(gca,'CLim',[-clim_max,clim_max],...
            'xlim',[warping_times(1) warping_times(end)],...
            'ydir','norm',...
            'ylim',[allfreqs(1) SUB_FREQ_LIMS(2)],...
            'yscale','log')
        %- set y-axis label
        if SUB_FREQ_LIMS(2) <= 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',FONT_SIZE);
        elseif SUB_FREQ_LIMS(2) <= 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',FONT_SIZE);
        end
        %- set x-axis & y-axis labels
        if j == 1
            ylabel(YTICK_LABEL,'FontSize',FONT_SIZE,'fontweight','bold');
            xlabel(XTICK_LABEL,'FontSize',FONT_SIZE);
        else
            xlabel('','FontSize',FONT_SIZE);
            ylabel('','fontsize',FONT_SIZE,'fontweight','bold');
        end
        %- set x-axis ticks
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca;
        ax.XAxis.FontSize = FONT_SIZE;
        title(alltitles{j});
        horiz_shift = horiz_shift + SHIFT_AMNT;
    end
    hold off;
    %- set color bar
    c = colorbar('XTick', -clim_max:0.2:clim_max);
    c.Position(1) = c.Position(1)+COLORBAR_SHIFT;
%     c.Limits = [-clim_max,clim_max];
%     caxis([-clim_max,clim_max]);
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',FONT_SIZE);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
end
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plots(STUDY,plottype_i,allersp,allfreqs,alltimes,pcond_ersp,...
    warping_times,colormap_ersp,sub_freq_lims,cluster_i,save_dir)
%     SAVE_STATS = false;
%     fcn = @erspStats;
    DO_PAIRWISE_COMP = true;
%     TERRAIN_DES_INT = 1;
%     SPEED_DES_INT = 2;
    TERRAIN_REF_CHAR = 'flat';
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    FIGURE_POSITION = [100,100,1480,350];
    PANEL_OFFSET = [-0.04,-0.03,0.02,-0.075]; % [LEFT,BOTTOM,WIDTH,HEIGHT];
%     clim_max = 0.5; %[-2,2];
    EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    COLOR_LIM_INTERVALS = [0.6,1.2,1.5];
    COLOR_LIM_ERR = 0.05;
    %% ================================================================= %%
    baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(end)); 
    freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
    allersp_cond_mean = cell(size(allersp));
    allersp_cond_std = cell(size(allersp));
    allersp_std_crop = cell(size(allersp));
    allersp_cond_mean_crop = cell(size(allersp));
    allerspdata_crop = cell(size(allersp));
    allersp_subjmean = cell(size(allersp));
    %- this step may help when singletrials are loaded, but when not, it is
    %useless
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            allersp_subjmean{cond_i,group_i} = zeros(size(allersp{cond_i,group_i}));
            for subj_i = 1:size(allersp{cond_i,group_i},3)
                allersp_subjmean{cond_i,group_i}(:,:,subj_i) = nanmean(allersp{cond_i,group_i}(:,:,subj_i),3);
%                 allersp_subjmean{cond_i,group_i}(:,:,subj_i) = mean(allersp{cond_i,group_i}(:,:,subj_i),3);
            end
        end
    end
    %## Calculate Stats
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            allersp_cond_mean{cond_i,group_i} = mean(allersp_subjmean{cond_i,group_i},3);
            allersp_cond_std{cond_i,group_i} = std(allersp_subjmean{cond_i,group_i},[],3);
            allersp_std_crop{cond_i,group_i} = allersp_cond_std{cond_i,group_i}(freqidx,baseidx,:);
        
            allersp_cond_mean_crop{cond_i,group_i} = allersp_cond_mean{cond_i,group_i}(freqidx,baseidx);
            allerspdata_crop{cond_i,group_i} = allersp_subjmean{cond_i,group_i}(freqidx,baseidx,:);
        end
    end
%     [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp_subjmean,allfreqs,alltimes);
%     [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp,allfreqs,alltimes);
    [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp,allfreqs,alltimes);
    [pcond_ersp_crop, ~, ~] = erspStats(STUDY,allerspdata_crop,allfreqs(freqidx),alltimes(baseidx));
%     [pcond_ersp_crop, ~, ~] = erspStats(STUDY,allersp_cond_mean_crop,allfreqs,alltimes);
    %## set titles
    if any(strcmp(STUDY.design(STUDY.currentdesign).variable(1).value,SPEED_REF_CHAR))
        condnames = SPEED_OVERRIDE_CHARS;
    else
        condnames = STUDY.design(STUDY.currentdesign).variable(1).value;
    end
    alltitles = std_figtitle('condnames',condnames, ...
        'clustname', sprintf('CL%i',cluster_i));
    %## set color limits
    climMat = [min(allersp_cond_mean_crop{4}(1:30,:),[],'all') max(allersp_cond_mean_crop{4}(1:30,:),[],'all')];
    clim_max = [];
    for i = 1:length(COLOR_LIM_INTERVALS)
        chk = any(climMat < (COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR)) && any(climMat > -(COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR));
        if chk
            clim_max = COLOR_LIM_INTERVALS(i);
            break
        end
    end
    if isempty(clim_max)
        clim_max = 1.5;
    end
    %% ================================================================= %%
%    std_plottf(alltimes(baseidx),allfreqs(freqidx),allersp_std_crop, 'datatype','ersp', 'plotmode','normal',...
%         'titles',alltitles,'caxis',[0,4]);
    [fig] = plot_tftopo(allersp_std_crop,alltimes(baseidx),allfreqs(freqidx),alltitles,pcond_ersp_crop{1},... %pcond_ersp{1},...
        warping_times,[0,4],colormap_ersp,...
        PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
    exportgraphics(gcf,[save_dir filesep sprintf('erspplottype%i_stdplot_des%i_cl%i.pdf',plottype_i,STUDY.currentdesign,cluster_i)],'Resolution',300);
    exportgraphics(gcf,[save_dir filesep sprintf('erspplottype%i_stdplot_des%i_cl%i.jpg',plottype_i,STUDY.currentdesign,cluster_i)]);
    %% ================================================================= %%
%     %- plot
%     std_plottf(alltimes(baseidx),allfreqs(freqidx),allersp_cond_mean_crop, 'datatype','ersp', 'plotmode','normal',...
%         'titles',alltitles,'caxis',[-clim_max,clim_max])
%     %- save
%     saveas(gcf,[save_dir filesep sprintf('erspplottype2_orig_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
%     saveas(gcf,[save_dir filesep sprintf('erspplottype2_orig_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    %% ================================================================= %%
    %- plot time freq
    plot_allersp = cell(size(allersp_cond_mean));
    for i = 1:size(allersp_cond_mean,1)
        for j = 1:size(allersp_cond_mean,2)
            plot_allersp{i,j} = mean(allersp_cond_mean{i,j},3);
        end
    end
    [fig] = plot_tftopo(plot_allersp,alltimes,allfreqs,alltitles,pcond_ersp_nocrop{1},... %pcond_ersp{1},...
        warping_times,clim_max,colormap_ersp,...
        PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
    exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_des%i_cl%i.jpg',plottype_i,STUDY.currentdesign,cluster_i)],'Resolution',300);
    exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_des%i_cl%i.pdf',plottype_i,STUDY.currentdesign,cluster_i)],...
        'Resolution',300,'ContentType','vector');
    
    close(fig)
    %% ================================================================= %%
    %- plot time freq
    plot_allersp = cell(size(allersp_cond_mean_crop));
    for i = 1:size(allersp_cond_mean_crop,1)
        for j = 1:size(allersp_cond_mean_crop,2)
            plot_allersp{i,j} = mean(allersp_cond_mean_crop{i,j},3);
        end
    end
    plot_freqs = allfreqs(freqidx);
    plot_times = alltimes(baseidx);
    plot_cond = pcond_ersp_crop{1};
    [fig] = plot_tftopo(plot_allersp,plot_times,plot_freqs,alltitles,plot_cond,...
        warping_times,clim_max,colormap_ersp,...
        PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
    exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_notfull_des%i_cl%i.jpg',plottype_i,STUDY.currentdesign,cluster_i)],'Resolution',300);
    exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_notfull_des%i_cl%i.pdf',plottype_i,STUDY.currentdesign,cluster_i)],...
        'Resolution',300,'ContentType','vector');
    
    close(fig)
    %% ================================================================= %%
    %## Pairwise comparison
    % compute difference ersps
    if DO_PAIRWISE_COMP
        chk_1 = strcmp(TERRAIN_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        chk_2 = strcmp(SPEED_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        cond_chars = [STUDY.design(STUDY.currentdesign).variable(1).value];
        if any(chk_1)
            refErspCond = TERRAIN_REF_CHAR;
            refErspCond_fext = TERRAIN_REF_CHAR;
            refErspCond_ind = find(chk_1);
        elseif any(chk_2)
             %SPEED_REF_CHAR; %'1p0';
            refErspCond_ind = find(chk_2);
            refErspCond = SPEED_OVERRIDE_CHARS{refErspCond_ind};
            refErspCond_fext = SPEED_REF_CHAR;
        else
            error('Condition for reference ersp not found in STUDY design: %s',[STUDY.design(STUDY.currentdesign).variable(1).value])
        end
        inds_to_comp = setdiff(1:length(alltitles),refErspCond_ind);
        alltitles_pw = alltitles; %(inds_to_comp);
        if ~isempty(refErspCond)
            % mask differenec ersps- check that it's sig. different from zer
            erspDiff = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
            erspDiff_wind = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
            %- calculate pairwise statistics between conditions of interest
            for c = inds_to_comp
                %-
                fprintf('Computing Pair Stat for %s - %s...\n',refErspCond,cond_chars{c})
                curr_ersp = allersp{c,1};
                ref_ersp = allersp{refErspCond_ind,1};
                [tmp, ~, ~] = erspStats(STUDY,{curr_ersp;ref_ersp},allfreqs,alltimes);
%                 [pcond_ersp, ~, ~] = feval(fcn,STUDY,{curr_ersp;ref_ersp},allfreqs,alltimes);
                erspDiff(c).raw = mean(curr_ersp-ref_ersp,3);
                erspDiff(c).masked = erspDiff(c).raw.*tmp{1,1};
                erspDiff(c).pcond = tmp{1,1};
                %-
                curr_ersp_wind = allersp{c,1}(freqidx,baseidx,:);
                ref_ersp_wind = allersp{refErspCond_ind,1}(freqidx,baseidx,:);
                [tmp, ~, ~] = erspStats(STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs(freqidx),alltimes(baseidx));
%                 [pcond_ersp, ~, ~] = feval(fcn,STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
                erspDiff_wind(c).raw = mean(curr_ersp_wind-ref_ersp_wind,3);
                erspDiff_wind(c).masked = erspDiff_wind(c).raw.*tmp{1,1};
                erspDiff_wind(c).pcond = tmp{1,1};
            end
        end
%         if SAVE_STATS
%             mkdir(fullfile(save_dir,['cluster_',num2str(cluster_i)]));
%             save(fullfile(save_dir,['cluster_',num2str(cluster_i)],['Common_baseline_2',file_keyword,'.mat']),...
%                 'pcond_ersp_crop','pgroup_ersp_crop','pinter_ersp_crop','erspDiff','erspDiff_wind');
%         end
        %% ============================================================= %%
        ersp_raw = {erspDiff.raw};
        ersp_pcond = {erspDiff.pcond};
        ersp_masked = {erspDiff.masked};
        ersp_raw = ersp_raw(~cellfun(@isempty,ersp_raw));
        ersp_pcond = ersp_pcond(~cellfun(@isempty,ersp_pcond));
        ersp_masked = ersp_masked(~cellfun(@isempty,ersp_masked));
        plot_alltitles = cell(size(inds_to_comp));
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            plot_alltitles{j} = sprintf('%s - %s',alltitles_pw{cond_i},refErspCond);
        end
%         plot_freqs = allfreqs(freqidx);
%         plot_times = alltimes(baseidx);
        [fig] = plot_contourf(ersp_raw,ersp_pcond,ersp_masked,alltimes,allfreqs,plot_alltitles,...
        warping_times,clim_max,colormap_ersp,PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
        exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_compare-%s_des%i_cl%i.jpg',plottype_i,refErspCond_fext,STUDY.currentdesign,cluster_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_stats_compare-%s_des%i_cl%i.pdf',plottype_i,refErspCond_fext,STUDY.currentdesign,cluster_i)],...
            'Resolution',300,'ContentType','vector');
        
        close(fig)
        %% ============================================================= %%
        ersp_raw = {erspDiff_wind.raw};
        ersp_pcond = {erspDiff_wind.pcond};
        ersp_masked = {erspDiff_wind.masked};
        ersp_raw = ersp_raw(~cellfun(@isempty,ersp_raw));
        ersp_pcond = ersp_pcond(~cellfun(@isempty,ersp_pcond));
        ersp_masked = ersp_masked(~cellfun(@isempty,ersp_masked));
        plot_alltitles = cell(size(inds_to_comp));
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            plot_alltitles{j} = sprintf('%s - %s',alltitles_pw{cond_i},refErspCond);
        end
        plot_freqs = allfreqs(freqidx);
        plot_times = alltimes(baseidx);
        [fig] = plot_contourf(ersp_raw,ersp_pcond,ersp_masked,plot_times,plot_freqs,plot_alltitles,...
            warping_times,clim_max,colormap_ersp,PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
        exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_notfull_stats_compare-%s_des%i_cl%i.jpg',plottype_i,refErspCond_fext,STUDY.currentdesign,cluster_i)],'Resolution',300);
        exportgraphics(fig,[save_dir filesep sprintf('erspplottype%i_notfull_stats_compare-%s_des%i_cl%i.pdf',plottype_i,refErspCond_fext,STUDY.currentdesign,cluster_i)],...
            'Resolution',300,'ContentType','vector');
        
        close(fig)
    end
end
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plot_1(STUDY,allersp,alltimes,allfreqs,...
    warping_times,colormap_ersp,sub_freq_lims,cluster_i,alpha,save_dir)
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    FIGURE_POSITION = [100,100,1480,300];
    BOOT_NITERS = 2000;
    COLOR_LIM_INTERVALS = [0.6,1.2,1.5];
    COLOR_LIM_ERR = 0.05;
    CLUSTER_THRESHOLD = 1000;
    MASK_ALPHA = 0.9;
    PANEL_OFFSET = [];
    EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    %##
    baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(end)); 
    freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
%     allersp_cond_mean = cell(size(allersp));
%     allersp_cond_mean_crop = cell(size(allersp));
%     allerspdata_crop = cell(size(allersp));
    allersp_subjmean = cell(size(allersp));
    ersp_subj_mean_base_crop = cell(size(allersp));
    ersp_sub_base = cell(size(allersp));
    ersp_subj_mean_base = cell(size(allersp));
    %- this step may help when singletrials are loaded, but when not, it is
    %useless
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            allersp_subjmean{cond_i,group_i} = zeros(size(allersp{cond_i,group_i}));
            for subj_i = 1:size(allersp{cond_i,group_i},3)
                allersp_subjmean{cond_i,group_i}(:,:,subj_i) = nanmean(allersp{cond_i,group_i}(:,:,subj_i),3);
            end
        end
    end
    %## Calculate Stats
%     for group_i = 1:size(allersp,2)
%         for cond_i = 1:size(allersp,1)
%             allersp_cond_mean{cond_i,group_i} = mean(allersp_subjmean{cond_i,group_i},3);
%             allersp_cond_mean_crop{cond_i,group_i} = allersp_cond_mean{cond_i,group_i}(freqidx,baseidx);
%             allerspdata_crop{cond_i,group_i} = allersp_subjmean{cond_i,group_i}(freqidx,baseidx,:);
%         end
%     end
    %-
%     freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
%     subj_in_cluster = unique(STUDY.cluster(cluster_i).sets); %Subjects in this cluster
%     allersp_subjmean = cell(length(allersp),1);
%     for cond_i = 1:size(allersp,1)
%         allersp_subjmean{cond_i,1}(:,:,:) = zeros(size(allersp{cond_i},1),size(allersp{cond_i},2),length(subj_in_cluster));
%     end
%     %-
%     allerspdata_remove = allersp;        
%     subj_i = 1;
%     sub = unique(STUDY.cluster(cluster_i).sets);
%     for n = 1:length(unique(STUDY.cluster(cluster_i).sets)) %1:size(allerspdata_meanSubj{1},3) % 1:length(unique(STUDY.cluster(cluster_i).sets))    
%         comp_ind = STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
%         %- comp not using
%         if cluster_i == 1 %bad comp
%             for cond_i = 1:size(allersp,1)
%                 allerspdata_remove{cond_i}(:,:,subj_i) = nan(size(allersp{cond_i},1),size(allersp{cond_i},2),1);
%                 allersp_subjmean{cond_i}(:,:,n) = nanmean( allerspdata_remove{cond_i}(:,:,subj_i:subj_i + length(comp_ind)-1),3);     
%             end
%         else
%             for cond_i = 1:length(cond_test)
%                 allersp_subjmean{cond_i}(:,:,n) = nanmean( allerspdata_remove{cond_i}(:,:,subj_i:subj_i + length(comp_ind)-1),3);
%             end
%         end
%         subj_i = subj_i+length(comp_ind);
%     end
    
    %- reorganize allerspdata
    
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            tmp = allersp_subjmean{cond_i,group_i}; 
            %- calc mean power for each person
            base_time = mean(tmp(:,baseidx,:),2); 
            %- calc mean power across participant
%             base_time_subj = mean(baseline_allcomp,3);
            %- subtract baseline (mean across time) for each person
            ersp_sub_base{cond_i,group_i} = tmp-repmat(base_time,1,length(alltimes));
            %- subtract baseline (mean across time & subjects) for each condition
%             ersp_subj_mean_base{cond_i,group_i} = mean(tmp-repmat(base_time_subj,1,length(alltimes)),3);
            ersp_subj_mean_base_crop{cond_i,group_i} = ersp_sub_base{cond_i,group_i}(:,baseidx);
        end
    end
    %- set titles
    if any(strcmp(STUDY.design(STUDY.currentdesign).variable(1).value,SPEED_REF_CHAR))
        condnames = SPEED_OVERRIDE_CHARS;
    else
        condnames = STUDY.design(STUDY.currentdesign).variable(1).value;
    end
    alltitles = std_figtitle('condnames',condnames, ...
        'clustname', sprintf('CL%i',cluster_i));
    %## set color limits
    climMat = [min(ersp_subj_mean_base_crop{end}(1:30,:),[],'all') max(ersp_subj_mean_base_crop{end}(1:30,:),[],'all')];
    clim_max = [];
    for i = 1:length(COLOR_LIM_INTERVALS)
        chk = any(climMat < (COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR)) && any(climMat > -(COLOR_LIM_INTERVALS(i)+COLOR_LIM_ERR));
        if chk
            clim_max = COLOR_LIM_INTERVALS(i);
            break
        end
    end
    if isempty(clim_max)
        clim_max = 1.5;
    end 
    %% ================================================================= %%
    %## (PLOT) Paper Figure for YA paper and IEEE NER - significance masked ERSP for high terrain
    clust_ersp = cell(size(allersp,1),size(allersp,2));
    clust_maskedersp = cell(size(allersp,1),size(allersp,2));
    for cond_i = 1:size(allersp,1)
        fprintf('Performing Stats for Condition %i & Cluster %i\n',cond_i,cluster_i);
        %- Within condition signficance? (i.e., is a particular subject
        %very different from the average...?)
        if ~isnan(alpha)
            tmp = ersp_sub_base{cond_i,1}(freqidx,baseidx,:);% this is already sub baseline
            tmp_mean = mean(tmp,3);
            boot_freq = 1:size(tmp,1);
            boot_subj = 1:size(tmp,3);
            boot_surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
            surro = zeros(size(tmp,1),size(tmp,2),BOOT_NITERS);
            %- scramble time samples and calculate the average across
            %all times and all frequencies and store that value.
            for n = 1:BOOT_NITERS
                boot_time = randi(size(tmp,2),[size(tmp,2),1]); % random time samples
                tmpSurro = mean(tmp(boot_freq,boot_time,boot_subj),3);
                surro(:,:,n) = tmpSurro; % save 2000 iterations of surrogates 
            end
            %- Pull length(subject) surrogate averages from distribution then calc mean across
            %surrogates 
            for n = 1:BOOT_NITERS
                bootIdx  = randi(BOOT_NITERS,[size(tmp,3),1]);
                tmpSurro = mean(surro(:,:,bootIdx),3);
                boot_surro(:,:,n) = tmpSurro;
            end
            pvalMap = stat_surrogate_pvals(boot_surro,tmp_mean,'both');
            pvalMap(pvalMap>1)=1; 
            [p_masked, ~, ~, ~] = fdr_bh(pvalMap,alpha,'pdep',1);
            % debri removal
            [labelMap,~] = bwlabeln(p_masked);
            tmpDisp = sort(labelMap(:),'descend');
%             [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
            [occurrence,idx,~] = histcounts(tmpDisp,unique(tmpDisp));
            kMask = ismember(labelMap,idx((occurrence<CLUSTER_THRESHOLD)));
            finalMask = p_masked-kMask;
            clust_ersp{cond_i} = tmp_mean; 
            tmp = clust_ersp{cond_i}; 
            tmp(~finalMask) = 0;
            clust_maskedersp{cond_i} = tmp;
        else
            clust_ersp{cond_i} = mean(ersp_subj_mean_base{cond_i},3);
            clust_maskedersp{cond_i} = clust_ersp;
        end
    end
    %% ============================================================= %%
    ersp_raw = clust_ersp;
    ersp_pcond = clust_maskedersp;
    ersp_masked = clust_maskedersp;
    plot_freqs = allfreqs(freqidx);
    plot_times = alltimes(baseidx);
    [fig] = plot_contourf(ersp_raw,ersp_pcond,ersp_masked,plot_times,plot_freqs,alltitles,...
        warping_times,clim_max,colormap_ersp,PANEL_OFFSET,EVENT_CHARS,sub_freq_lims,FIGURE_POSITION);
    exportgraphics(fig,[save_dir filesep sprintf('erspplot_des%i_cl%i_bootstats.jpg',STUDY.currentdesign,cluster_i)],'Resolution',300);
    exportgraphics(fig,[save_dir filesep sprintf('erspplot_des%i_cl%i_bootstats.pdf',STUDY.currentdesign,cluster_i)],'Resolution',300,'ContentType','vector');
    
    close(fig)
end
%% (SUBFUNCTION) ======================================================= %%
%##
% exerpt from std_erspplot 
function [pcond, pgroup, pinter] = erspStats(STUDY,allersp,allfreqs,alltimes)
    %get stats parameters
    stats = STUDY.etc.statistics;
    stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
    if isempty(STUDY.design(STUDY.currentdesign).variable)
        stats.paired = { };
    else
        stats.paired = { STUDY.design(STUDY.currentdesign).variable(:).pairing };
    end

    %- get ersp params
    params = STUDY.etc.erspparams;
    params.plottf = [];
    % select specific time and freq
    % -----------------------------
    if ~isempty(params.plottf)
        if length(params.plottf) < 3
            params.plottf(3:4) = params.plottf(2);
            params.plottf(2)   = params.plottf(1);
        end
        [~, fi1] = min(abs(allfreqs-params.plottf(1)));
        [~, fi2] = min(abs(allfreqs-params.plottf(2)));
        [~, ti1] = min(abs(alltimes-params.plottf(3)));
        [~, ti2] = min(abs(alltimes-params.plottf(4)));
        for index = 1:length(allersp(:))
            allersp{index} = mean(mean(allersp{index}(fi1:fi2,ti1:ti2,:,:),1),2);
            allersp{index} = reshape(allersp{index}, [1 size(allersp{index},3) size(allersp{index},4) ]);
        end

        % prepare channel neighbor matrix for Fieldtrip
        statstruct = std_prepare_neighbors(STUDY, ALLEEG);
        stats.fieldtrip.channelneighbor = statstruct.etc.statistics.fieldtrip.channelneighbor;

        params.plottf = { params.plottf(1:2) params.plottf(3:4) };
        [pcond, pgroup, pinter] = std_stat(allersp, stats);
        if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end % single subject STUDY
    else
        [pcond, pgroup, pinter] = std_stat(allersp, stats);
        if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
                (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1))
            pcond = {}; pgroup = {}; pinter = {};
            disp('No statistics possible for single subject STUDY');
        end % single subject STUDY
    end
end