function [] = mim_custom_ersp_plots(STUDY,cond_test,warping_times,cluster_ind,...
    cluster_load_ind,des_i,save_dir,varargin)
%MIM_GEN_ERSP_STATS Summary of this function goes here
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
% alpha = 0.05; 
% plotComp = 0;
% YlimMax = 50;
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
colormap_ersp = colormap_ersp(end:-1:1,:);
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
%## subject plots
if DO_SUBJ_PLOTS
    % plot each component with baseline to be the average of entire epoch
    %- raw ersp data no baselining
    if inside_load_flag
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data = par_load(fpath,fname);
        allersp = ersp_data.allerspdata;
        alltimes = ersp_data.alltimes;
        allfreqs = ersp_data.allfreqs;
    end
    fprintf('Using ERSP data:\n');
    disp(allersp)
%     cl_in = cluster_i + 1; % offset for parent cluster
    ersp_single_subj_plot(STUDY,allersp,alltimes,allfreqs,...
        warping_times,SUB_FREQ_LIMS,cluster_ind,save_dir);
    
end
%% 1. Baseline correction = average of epoch within condition
if DO_BASELINE_CORRECT_1
    %- raw ersp data no baselining
    if inside_load_flag
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data = par_load(fpath,fname);
        allersp = ersp_data.allerspdata;
        alltimes = ersp_data.alltimes;
        allfreqs = ersp_data.allfreqs;
    end
    fprintf('Using ERSP data:\n');
    disp(allersp)
    ersp_baseline_plot_1(STUDY,allersp,alltimes,allfreqs,cond_test,...
        warping_times,colormap_ersp,SUB_FREQ_LIMS,cluster_ind,ERSP_ALPHA,CLUSTER_CLIM_MATCH,save_dir);
end
%% 2. Common Baseline correction = average of epoch across condition (commonbaseline)
% Read allerspdata3, allerspdata3 should already be baselined
% allerspdata,  
if DO_BASELINE_CORRECT_2
    %- normalized ersp data
    if inside_load_flag
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_norm_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data_norm = par_load(fpath,fname);
        allersp = ersp_data_norm.allerspdata;
        alltimes = ersp_data_norm.alltimes;
        allfreqs = ersp_data_norm.allfreqs;
        pcond = ersp_data_norm.pcond;
    end
    fprintf('Using ERSP data:\n');
    disp(allersp)
%     fprintf('Using ERSP stats:\n');
%     disp(pcond)
    ersp_baseline_plot_2(STUDY,allersp,allfreqs,alltimes,pcond,...
        cond_test,warping_times,colormap_ersp,SUB_FREQ_LIMS,des_i,cluster_ind,ERSP_ALPHA,CLUSTER_CLIM_MATCH,save_dir);
end
%% 3. Single trial full epoch correction + Common Baseline correction = average of epoch across condition (commonbaseline)
% Read allerspdata4, allerspdata4 should already be baselined
if DO_BASELINE_CORRECT_3
    if inside_load_flag
        fname = strsplit(STUDY.etc.mim_gen_ersp_data(des_i,cluster_load_ind).ersp_normcb_fpaths,'/');
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        ersp_data_normcb = par_load(fpath,fname);
        allersp = ersp_data_normcb.allerspdata;
        alltimes = ersp_data_normcb.alltimes;
        allfreqs = ersp_data_normcb.allfreqs;
        pcond = ersp_data_normcb.pcond;
    end
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
    ersp_baseline_plot_3(STUDY,allersp,allfreqs,alltimes,...
        pcond,warping_times,colormap_ersp,SUB_FREQ_LIMS,des_i,cluster_ind,ERSP_ALPHA,CLUSTER_CLIM_MATCH,save_dir)
end
%% 4. 
% if DO_BASELINE_CORRECT_4
%     ersp_baseline_plot_4()
% end

end
%% (SUBFUNCTIONS) ====================================================== %%
%## 
function [] = ersp_single_subj_plot(STUDY,allersp,alltimes,allfreqs,...
    warping_times,sub_freq_lims,cluster_i,save_dir)
%     FIGURE_POSITION = [100,100,1620,200];
    %## 
%     freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
    for i = 1:size(allersp{1},3) %1:length(STUDY.cluster(cluster_i).comps)
        ic = STUDY.cluster(cluster_i).comps(i);
        sub = STUDY.datasetinfo(STUDY.cluster(cluster_i).sets(i)).subject; 
        baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(5));
        erspdata = allersp{1}(:,:,i);
        baseline = mean(erspdata(:,baseidx,:),2);
        curr_ersp = erspdata(:,:,:)-repmat(baseline,1,length(alltimes));
        curr_ersp = mean(curr_ersp,3);
        curr_maskedersp = curr_ersp;
        %- plot
        figure('renderer','Painters');
        tftopo(curr_maskedersp,alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan nan nan],...
            'vert',warping_times(1:5),'logfreq','native');
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
        xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
        title(strcat({'Cluster '},num2str(cluster_i)));
        cbar('vert');
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, 'clustname', STUDY.cluster(cluster_i).name,...
            'subject', sub, 'compnames', num2str(ic));
        %- reorganize allerspdata
        single_comp_ersp = cell(length(allersp),1);
        single_comp_ersp_crop = cell(length(allersp),1);
        for cond_i = 1: length(allersp)
            erspdata = allersp{cond_i}(:,:,i);
            baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(5));                
            baseline = mean(erspdata(:,baseidx,:),2);                
            single_comp_ersp{cond_i,1} = mean(allersp{cond_i}(:,:,i)-repmat(baseline,1,length(alltimes)),3);
            single_comp_ersp_crop{cond_i,1} = single_comp_ersp{cond_i,1}(:,baseidx);
        end
        std_plottf(alltimes(baseidx),allfreqs, single_comp_ersp_crop, 'datatype','ersp', 'plotmode','normal','titles',alltitles)
        %- save
        saveas(gcf,[save_dir filesep sprintf('%s_within_%i_%i_stats_lim%i.jpg',sub,STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    end
    close all
end
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plot_1(STUDY,allersp,alltimes,allfreqs,cond_test,...
    warping_times,colormap_ersp,sub_freq_lims,cluster_i,alpha,cluster_match_chars,save_dir)
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    FIGURE_POSITION = [100,100,1480,300];
    BOOT_NITERS = 2000;
    COLOR_LIM_INTERVALS = [0.6,1.2,1.5];
    COLOR_LIM_ERR = 0.05;
    CLUSTER_THRESHOLD = 1000;
    MASK_ALPHA = 0.9;
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
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
%     for group_i = 1:size(allersp,2)
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
            clust_ersp = tmp_mean; 
            clust_maskedersp = clust_ersp; 
            clust_maskedersp(~finalMask) = 0;
        else
            clust_ersp = mean(ersp_subj_mean_base{cond_i},3);
            clust_maskedersp = clust_ersp;
        end
        %## (PLOT)
        subplot(1,length(allersp),cond_i)
        colormap(colormap_ersp);
        faceAlpha_mask = ones(size(clust_maskedersp))*MASK_ALPHA; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        faceAlpha_mask(clust_maskedersp ~=0 ) = 0; %0 is significant? 1 is not? 
        contourf(alltimes(baseidx), allfreqs(freqidx), clust_ersp,200,...
                   'linecolor','none');hold on;
        imagesc(alltimes(baseidx),allfreqs(freqidx),clust_maskedersp,'AlphaData',faceAlpha_mask);
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        set(gca,'clim',[-clim_max, clim_max],'xlim',[warping_times(1) warping_times(end)],...
            'ydir','norm','ylim',[allfreqs(1) sub_freq_lims(2)],'yscale','log')
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
        end
        if cond_i == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
        end
        xlabel('','Fontsize',12);
        title(alltitles{cond_i});
%         set(gca,'TitleFontSizeMultiplier',0.5);
        set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
    end
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-clim_max, clim_max];
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',90);
    hL.Position(1) = hL.Position(1)+1.2;
    hL.Position(2) = .13;
    saveas(gcf,[save_dir filesep sprintf('erspplot1_des%i_cl%i_stats_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    saveas(gcf,[save_dir filesep sprintf('erspplot1_des%i_cl%i_stats_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    close all
end
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plot_2(STUDY,allersp,allfreqs,alltimes,pcond_ersp,cond_test,...
    warping_times,colormap_ersp,sub_freq_lims,des_i,cluster_i,alpha,cluster_match_chars,save_dir)
%     SAVE_STATS = false;
    fcn = @erspStats;
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
    %## Subtracct out mean per condition
%     subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
%     allersp_subjmean = cell(length(cond_test),1);
%     for j = 1:size(allersp,1)
%         allersp_subjmean{j}(:,:,:) = zeros(size(allersp{j},1),size(allersp{j},2),length(subj_in_cluster));
%     end
%     allerspdata_remove = allersp;
%     subj_i = 1;
%     sub = unique(STUDY.cluster(cluster_i).sets);
%     for n = 1:length(unique(STUDY.cluster(cluster_i).sets))    
%         comp_ind =  STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
%         for j = 1:length(allersp_subjmean)
%             allersp_subjmean{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,subj_i:subj_i + length(comp_ind)-1),3);
%             fprintf('iter...\n');
%             disp(subj_i:subj_i + length(comp_ind)-1)
%             disp(comp_ind)
%             disp(n)
%         end     
%         subj_i = subj_i+length(comp_ind);
%     end
%     subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);
    baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(end)); 
    freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
    allersp_cond_mean = cell(size(allersp));
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
            allersp_cond_mean_crop{cond_i,group_i} = allersp_cond_mean{cond_i,group_i}(freqidx,baseidx);
            allerspdata_crop{cond_i,group_i} = allersp_subjmean{cond_i,group_i}(freqidx,baseidx,:);
        end
    end
%     [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp_subjmean,allfreqs,alltimes);
%     [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp,allfreqs,alltimes);
    
    [pcond_ersp_crop, ~, ~] = erspStats(STUDY,allerspdata_crop,allfreqs,alltimes);
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
     %- plot
    std_plottf(alltimes(baseidx),allfreqs(freqidx),allersp_cond_mean_crop, 'datatype','ersp', 'plotmode','normal',...
        'titles',alltitles,'caxis',[-clim_max,clim_max])
    %- save
    saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    
    %% ================================================================= %%
    %##  Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    for j = 1:length(allersp)
        clust_ersp = mean(allersp_cond_mean{j},3);
        subplot(1,length(allersp)+1,j)
        tftopo(clust_ersp,alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan -clim_max clim_max],...
            'logfreq','native');
        %- adjust panel height
        cur_pos = get(gca,'Position');
        cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
        cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
        cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1);
        cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
        set(gca,'Position',cur_pos)
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        %- set ylims
        ylim(log(sub_freq_lims))
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end  
        %- set color lims
        set(gca,'clim',[-clim_max, clim_max]);
        %- set y-axis labels
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        %- set x-axis labels
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
%         ylim()
        %- title
        title(alltitles{j});   
    end
    %## Add Stats To Plot
    subplot(1,length(allersp)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp{1}),alltimes,allfreqs,'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-clim_max, clim_max]],...
        'logfreq','native')
    %- title
    title('(eeglab) bootstrap stats');
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-clim_max, clim_max];
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust panel height
    cur_pos = get(gca,'Position');
    cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
    cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
    cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1)+0.025;
    cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
    set(gca,'Position',cur_pos);
    %- add vertical line
    xline(gca,warping_times(1),'k--');
    xline(gca,warping_times(2),'k--');
    xline(gca,warping_times(3),'k--');
    xline(gca,warping_times(4),'k--');
    xline(gca,warping_times(5),'k--');
    %- set ylims
    ylim(log(sub_freq_lims))
    if sub_freq_lims(2) <= 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif sub_freq_lims(2) <= 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    %- set y-axis labels
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    %- set x-axis labels
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
    xtickangle(45)
    ax = gca; 
    ax.XAxis.FontSize = 8;
    %- save
    if ~isnan(alpha)
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_des%i_cl%i.fig',STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_des%i_cl%i.jpg',STUDY.currentdesign,cluster_i)])
    else
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_des%i_cl%i.fig',STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_des%i_cl%i.jpg',STUDY.currentdesign,cluster_i)])
    end
    %% ================================================================= %%
    
    %##  Figure: tfplot with preset range
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    for j = 1:length(allersp)
        clust_ersp = mean(allersp_cond_mean_crop{j},3);
        subplot(1,length(allersp)+1,j)
        tftopo(clust_ersp,alltimes(baseidx),allfreqs(freqidx),'limits',... 
            [warping_times(1) warping_times(end) nan nan -clim_max clim_max],...
            'logfreq','native');
        %- adjust panel heights
        cur_pos = get(gca,'Position');
        cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
        cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
        cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1);
        cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
        set(gca,'Position',cur_pos)
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        %- set ylims
        ylim(log(sub_freq_lims))
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',log([4.01,8,12,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end  
        %- set color lims
        set(gca,'clim',[-clim_max, clim_max]);
        %- set y-axis labels
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        %- set x-axis labels
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        %- title
        title(alltitles{j});  
    end
    %## Add Stats To Plot
    subplot(1,length(allersp)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp_crop{1}),alltimes(baseidx),allfreqs(freqidx),'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-clim_max, clim_max]],...
        'logfreq','native');
    %- title
    title('(custom) bootstrap stats');
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-clim_max, clim_max];
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust panel height
    cur_pos = get(gca,'Position');
    cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
    cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
    cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1)+0.025;
    cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
    set(gca,'Position',cur_pos);
    %- add vertical line
    xline(gca,warping_times(1),'k--');
    xline(gca,warping_times(2),'k--');
    xline(gca,warping_times(3),'k--');
    xline(gca,warping_times(4),'k--');
    xline(gca,warping_times(5),'k--');
    %- set ylims
    ylim(log(sub_freq_lims))
    if sub_freq_lims(2) <= 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif sub_freq_lims(2) <= 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    %- set y-axis labels
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    %- set x-axis labels
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
    xtickangle(45)
    ax = gca; 
    ax.XAxis.FontSize = 8;
    if ~isnan(alpha)
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.fig',STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_notfull_des%i_cl%i.jpg',STUDY.currentdesign,cluster_i)])
    else
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_notfull_des%i_cl%i.fig',STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_notfull_des%i_cl%i.jpg',STUDY.currentdesign,cluster_i)])
    end
    %% ================================================================= %%
    %## Pairwise comparison
    % compute difference ersps
    if DO_PAIRWISE_COMP
        chk_1 = strcmp(TERRAIN_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        chk_2 = strcmp(SPEED_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        cond_chars = [STUDY.design(STUDY.currentdesign).variable(1).value];
        if any(chk_1)
            refErspCond = TERRAIN_REF_CHAR;
            refErspCond_ind = find(chk_1);
        elseif any(chk_2)
            refErspCond = SPEED_REF_CHAR; %'1p0';
            refErspCond_ind = find(chk_2);
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
                [tmp, ~, ~] = erspStats(STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
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
        figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            subplot(1,length(allersp),cond_i)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff(cond_i).pcond))*0.9; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff(cond_i).pcond == 1) = 0; %0 is significant? 1 is not?
            contourf(alltimes, allfreqs, erspDiff(cond_i).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes,allfreqs,erspDiff(cond_i).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            xline(gca,warping_times(1),'k--');
            xline(gca,warping_times(2),'k--');
            xline(gca,warping_times(3),'k--');
            xline(gca,warping_times(4),'k--');
            xline(gca,warping_times(5),'k--');
            %- color lim
            set(gca,'clim',[-clim_max,clim_max],'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) sub_freq_lims(2)],'yscale','log')
            %- set y-axis label
            if sub_freq_lims(2) <= 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif sub_freq_lims(2) <= 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            end
            %- set x-axis label
            xlabel('','Fontsize',10);
            set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(sprintf('%s - %s',alltitles_pw{cond_i},refErspCond));  
        end
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.05;
        c.Limits = [-clim_max,clim_max];
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        %- save
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_compare-%s_des%i_cl%i.fig',refErspCond,STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_stats_compare-%s_des%i_cl%i.jpg',refErspCond,STUDY.currentdesign,cluster_i)])
    %% ================================================================= %%
        %## Figure not full window
%         if ~isempty(refErspCond)
%             data = [];
%             for j = 1:3 %1:length(erspDiff_wind)
%                 data = [data, reshape(mean(erspDiff_wind(inds_to_comp(j)).raw,3).',1,[])];
%             end
%             IQR = iqr(data); %interquartile range
%             Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
%             myMin = -round(mean(data)+1.5*IQR,1);
%             erspDiff_clim = [myMin myMin*(-1)];
%         end
%         clim = erspDiff_clim;

        figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            subplot(1,length(allersp),cond_i)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff_wind(cond_i).pcond))*0.9; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff_wind(cond_i).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes, allfreqs, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes(baseidx), allfreqs(freqidx), erspDiff_wind(cond_i).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes(baseidx),allfreqs(freqidx),erspDiff_wind(cond_i).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            xline(gca,warping_times(1),'k--');
            xline(gca,warping_times(2),'k--');
            xline(gca,warping_times(3),'k--');
            xline(gca,warping_times(4),'k--');
            xline(gca,warping_times(5),'k--');
            set(gca,'clim',[-clim_max,clim_max],'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) sub_freq_lims(2)],'yscale','log')
            if sub_freq_lims(2) <= 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif sub_freq_lims(2) <= 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            end
            xlabel('','Fontsize',10);
            set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(sprintf('%s - %s',alltitles_pw{cond_i},refErspCond));  
        end
        %- set color bar
        c = colorbar();
        c.Position(1) = c.Position(1)+0.05;
        c.Limits = [-clim_max,clim_max];
        %- color bar label
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = .13;
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_notfull_stats_compare-%s_des%i_cl%i.fig',refErspCond,STUDY.currentdesign,cluster_i)])
        saveas(gcf,[save_dir filesep sprintf('erspplottype2_notfull_stats_compare-%s_des%i_cl%i.jpg',refErspCond,STUDY.currentdesign,cluster_i)])
    end
    close all
end
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plot_3(STUDY,allersp,allfreqs,alltimes,pcond_ersp,...
    warping_times,colormap_ersp,sub_freq_lims,des_i,cluster_i,alpha,cluster_match_chars,save_dir)
%     SAVE_STATS = false;
%     fcn = @erspStats;
    DO_PAIRWISE_COMP = true;
%     TERRAIN_DES_INT = 1;
%     SPEED_DES_INT = 2;    
    TERRAIN_REF_CHAR = 'flat';
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    FIGURE_POSITION = [100,100,1480,350];
    EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    PANEL_OFFSET = [-0.04,-0.03,0.02,-0.075]; % [LEFT,BOTTOM,WIDTH,HEIGHT];
    COLOR_LIM_INTERVALS = [0.6,1.2,1.5];
    COLOR_LIM_ERR = 0.05;
    %##
%     allerspdata_to_use = allersp;
%     subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
%     allerspdata_meanSubj = cell(length(allersp),1);
%     for j = 1:length(allersp)
%         allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata_to_use{j},1),size(allerspdata_to_use{j},2),length(subj_in_cluster));
%     end
%     allerspdata_remove = allerspdata_to_use;        
%     subj_i = 1;
%     sub = unique(STUDY.cluster(cluster_i).sets);
%     for n = 1:length(unique(STUDY.cluster(cluster_i).sets)) %1:size(allerspdata_meanSubj{1},3) %
%         comp_ind =  STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
%         for j = 1:size(allersp)
%             allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,subj_i:subj_i + length(comp_ind)-1),3);
%         end     
%         subj_i = subj_i+length(comp_ind);
%     end
    baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(end)); 
    freqidx = find(allfreqs>=sub_freq_lims(1) & allfreqs<=sub_freq_lims(2));
    allersp_cond_mean = cell(size(allersp));
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
            end
        end
    end
    %## Calculate Stats
    for group_i = 1:size(allersp,2)
        for cond_i = 1:size(allersp,1)
            allersp_cond_mean{cond_i,group_i} = mean(allersp_subjmean{cond_i,group_i},3);
            allersp_cond_mean_crop{cond_i,group_i} = allersp_cond_mean{cond_i,group_i}(freqidx,baseidx);
            allerspdata_crop{cond_i,group_i} = allersp_subjmean{cond_i,group_i}(freqidx,baseidx,:);
        end
    end
%     [pcond_ersp_nocrop, ~, ~] = erspStats(STUDY,allersp_subjmean,allfreqs,alltimes);
    [pcond_ersp_crop, ~, ~] = erspStats(STUDY,allerspdata_crop,allfreqs,alltimes);
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
    %- plot
    std_plottf(alltimes(baseidx),allfreqs(freqidx),allersp_cond_mean_crop, 'datatype','ersp', 'plotmode','normal',...
        'titles',alltitles,'caxis',[-clim_max,clim_max])
    %- save
    saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    %% (RECOMPUTE STATS) =============================================== %%
    %## (PLOT) Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    for j = 1:length(allersp)
        clust_ersp = mean(allersp_cond_mean{j},3);
        clust_maskedersp = clust_ersp;
        subplot(1,size(allersp,1)+1,j)
        tftopo(clust_maskedersp,alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan -clim_max clim_max],...
            'logfreq','native');
        %- adjust panel height
        cur_pos = get(gca,'Position');
        cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
        cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
        cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1);
        cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
        set(gca,'Position',cur_pos)
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        %- set ylims
        ylim(log(sub_freq_lims))
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',log([4.01,8,12,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end  
        %- set color lims
        set(gca,'clim',[-clim_max, clim_max]);
        %- set y-axis labels
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        %- set x-axis labels
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        %- title
        title(alltitles{j});
    end
    % --- add stats plot ------------
    subplot(1,size(allersp,1)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp{1}),alltimes,allfreqs,'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-clim_max, clim_max]],...
        'logfreq','native');
    %- title
    title('(eeglab) F stats');
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-clim_max, clim_max];
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust panel height
    cur_pos = get(gca,'Position');
    cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
    cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
    cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1)+0.025;
    cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
    set(gca,'Position',cur_pos);
    %- add vertical line
    xline(gca,warping_times(1),'k--');
    xline(gca,warping_times(2),'k--');
    xline(gca,warping_times(3),'k--');
    xline(gca,warping_times(4),'k--');
    xline(gca,warping_times(5),'k--');
    %- set ylims
    ylim(log(sub_freq_lims))
    if sub_freq_lims(2) <= 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif sub_freq_lims(2) <= 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    %- set y-axis labels
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    %- set x-axis labels
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
    xtickangle(45)
    ax = gca; 
    ax.XAxis.FontSize = 8;
    if ~isnan(alpha)
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    else
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    end
    %% ================================================================= %%
    %## (PLOT) Figure: tfplot with preset range
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    for j = 1:length(allersp)
        clust_ersp = mean(allersp_cond_mean_crop{j},3);
        clust_maskedersp = clust_ersp;
        subplot(1,length(allersp)+1,j)
        tftopo(clust_maskedersp,alltimes(baseidx),allfreqs(freqidx),'limits',... 
            [warping_times(1) warping_times(end) nan nan -clim_max clim_max],...
            'logfreq','native'); 
         %- adjust panel height
        cur_pos = get(gca,'Position');
        cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
        cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
        cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1);
        cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
        set(gca,'Position',cur_pos)
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        %- set ylims
        ylim(log(sub_freq_lims))
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',log([4.01,8,12,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end  
        %- set color lims
        set(gca,'clim',[-clim_max, clim_max]);
        %- set y-axis labels
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        %- set x-axis labels
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        %- title
        title(alltitles{j});
    end
    % --- add stats plot ------------
    subplot(1,length(allersp)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp_crop{1}),alltimes(baseidx),allfreqs(freqidx),'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-clim_max, clim_max]],...
        'logfreq','native');
    %- title
    title('F stats');
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-clim_max, clim_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust panel height
    cur_pos = get(gca,'Position');
    cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
    cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
    cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1)+0.025;
    cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
    set(gca,'Position',cur_pos);
    %- add vertical line
    xline(gca,warping_times(1),'k--');
    xline(gca,warping_times(2),'k--');
    xline(gca,warping_times(3),'k--');
    xline(gca,warping_times(4),'k--');
    xline(gca,warping_times(5),'k--');
    %- set ylims
    ylim(log(sub_freq_lims))
    if sub_freq_lims(2) <= 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif sub_freq_lims(2) <= 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    %- set y-axis labels
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    %- set x-axis labels
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
    xtickangle(45)
    ax = gca; 
    ax.XAxis.FontSize = 8;
    if ~isnan(alpha)
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_notfull_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_notfull_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    else
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_notfull_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_notfull_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    end
    %% ================================================================= %%
    %## (PLOT) Pairwise comparison
    % compute difference ersps
    if DO_PAIRWISE_COMP
        chk_1 = strcmp(TERRAIN_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        chk_2 = strcmp(SPEED_REF_CHAR,[STUDY.design(STUDY.currentdesign).variable(1).value]);
        if any(chk_1)
            refErspCond = TERRAIN_REF_CHAR;
            refErspCond_ind = find(chk_1);
        elseif any(chk_2)
            refErspCond = SPEED_REF_CHAR; %'1p0';
            refErspCond_ind = find(chk_2);
        else
            error('Condition for reference ersp not found in STUDY design: %s',[STUDY.design(STUDY.currentdesign).variable(1).value])
        end
        inds_to_comp = setdiff(1:length(alltitles),refErspCond_ind);
        alltitles_pw = alltitles; %(inds_to_comp);
        if ~isempty(refErspCond)
            % mask differenec ersps- check that it's sig. different from zero
            erspDiff = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
            erspDiff_wind = struct('raw',cell(length(inds_to_comp),1),'masked',cell(length(inds_to_comp),1),'pcond',cell(length(inds_to_comp),1));
            for c = inds_to_comp
                curr_ersp = allersp{c,1};
                ref_ersp = allersp{refErspCond_ind,1};
                [pcond_ersp, ~, ~] = erspStats(STUDY,{curr_ersp;ref_ersp},allfreqs,alltimes);
%                 [pcond, ~, ~] = feval(fcn, STUDY,{curr_ersp;ref_ersp},allfreqs,alltimes);
                erspDiff(c).raw = mean(curr_ersp-ref_ersp,3);
                erspDiff(c).masked = erspDiff(c).raw.*pcond_ersp{1,1};
                erspDiff(c).pcond = pcond_ersp{1,1};

                curr_ersp_wind = allersp{c,1}(freqidx,baseidx,:);
                ref_ersp_wind = allersp{refErspCond_ind,1}(freqidx,baseidx,:);
                [pcond_ersp, ~, ~] = erspStats(STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
%                 [pcond, ~, ~] = feval(fcn, STUDY,{curr_ersp_wind;ref_ersp_wind},allfreqs,alltimes);
                erspDiff_wind(c).raw = mean(curr_ersp_wind-ref_ersp_wind,3);
                erspDiff_wind(c).masked = erspDiff_wind(c).raw.*pcond_ersp{1,1};
                erspDiff_wind(c).pcond = pcond_ersp{1,1};
            end
        end
%         if SAVE_STATS
%             mkdir(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)]));
%             save(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)],['Common_baseline_3',file_keyword,'.mat']),...
%                 'pcond_ersp_crop','pgroup_ersp_crop','pinter_ersp_crop','erspDiff','erspDiff_wind');
%         end
        %% ================================================================= %%
        %## (PLOT) Figure Full comparing high terrain with flat terrain
        % - set clim for erspDiff
%         if ~isempty(refErspCond)
%             data = [];
%             for j = 1:length(erspDiff)
%                 data = [data, reshape(mean(erspDiff(j).raw,3).',1,[])];
%             end
%             IQR = iqr(data); %interquartile range
% %             Q1 = quantile(data,0.25);
%     %         myMin = round(Q1-1.5*IQR,1);
%             myMin = -round(mean(data)+1.5*IQR,1);
%             erspDiff_clim = [myMin myMin*(-1)];
%         end
%         clim = erspDiff_clim;
        %## figure
        figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            subplot(1,length(allersp),cond_i)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff(cond_i).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff(cond_i).pcond == 1) = 0; %0 is significant? 1 is not? 
            contourf(alltimes, allfreqs, erspDiff(cond_i).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes,allfreqs,erspDiff(cond_i).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            xline(gca,warping_times(1),'k--');
            xline(gca,warping_times(2),'k--');
            xline(gca,warping_times(3),'k--');
            xline(gca,warping_times(4),'k--');
            xline(gca,warping_times(5),'k--');
            set(gca,'clim',[-clim_max,clim_max],'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) sub_freq_lims(2)],'yscale','log')
            if sub_freq_lims(2) <= 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif sub_freq_lims(2) <= 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',10);
            title(sprintf('%s - %s',alltitles_pw{cond_i},refErspCond));  
        end
        %- set color bar
        c = colorbar();
        c.Limits = [-clim_max,clim_max];
        c.Position(1) = c.Position(1)+0.05;
        %- color bar title
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        set(hL,'Rotation',0);
    %     colorbar
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_compare-%s_des%i_cl%i_lim%i.fig',refErspCond,STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_compare-%s_des%i_cl%i_lim%i.jpg',refErspCond,STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        %% ================================================================= %%
        %## (PLOT)
        if ~isempty(refErspCond)
            data = [];
            for j = 1:length(erspDiff_wind)
                data = [data, reshape(mean(erspDiff_wind(j).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
%             Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+1.5*IQR,1);
            erspDiff_clim = [myMin myMin*(-1)];
        end
        clim = erspDiff_clim;
        figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
        for j = 1:length(inds_to_comp)
            cond_i = inds_to_comp(j);
            subplot(1,length(allersp),cond_i)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff_wind(cond_i).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff_wind(cond_i).pcond == cond_i) = 0; %0 is significant? 1 is not?
            contourf(alltimes(baseidx), allfreqs(freqidx), erspDiff_wind(cond_i).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes(baseidx),allfreqs(freqidx),erspDiff_wind(cond_i).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            xline(gca,warping_times(1),'k--');
            xline(gca,warping_times(2),'k--');
            xline(gca,warping_times(3),'k--');
            xline(gca,warping_times(4),'k--');
            xline(gca,warping_times(5),'k--');
            %- color lim
            set(gca,'clim',[-clim_max,clim_max],'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) sub_freq_lims(2)],'yscale','log')
            %- set y-axis label
            if sub_freq_lims(2) <= 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif sub_freq_lims(2) <= 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            end
            %- set x-axis label
            xlabel('','Fontsize',10);
            set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(sprintf('%s - %s',alltitles_pw{cond_i},refErspCond));
        end
        %- set color bar
        c = colorbar();
        c.Limits = [-clim_max,clim_max];
        c.Position(1) = c.Position(1)+0.05;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',0);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        set(hL,'Rotation',0);    
    	%- save
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_notfull_stats_compare-%s_des%i_cl%i_lim%i.fig',refErspCond,STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_notfull_stats_compare-%s_des%i_cl%i_lim%i.jpg',refErspCond,STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    end
    close all
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
%% (SUBFUNCTION) ======================================================= %%
%##
function [] = ersp_baseline_plot_4(STUDY,allersp,allfreqs,alltimes,pcond_ersp,...
    warping_times,colormap_ersp,sub_freq_lims,cluster_i,alpha,cluster_match_chars,save_dir)
%     SAVE_STATS = false;
%     fcn = @erspStats;
%     DO_PAIRWISE_COMP = true;
%     TERRAIN_DES_INT = 1;
%     SPEED_DES_INT = 2;    
%     TERRAIN_REF_CHAR = 'flat';
    SPEED_REF_CHAR = '1p0';
    SPEED_OVERRIDE_CHARS = {'0.25 m/s','0.5 m/s','0.75 m/s','1.0 m/s'};
    FIGURE_POSITION = [100,100,1480,350];
    EVENT_CHARS = {'RHS','LTO','LHS','RTO','RHS'};
    PANEL_OFFSET = [-0.04,-0.03,0.02,-0.075]; % [LEFT,BOTTOM,WIDTH,HEIGHT];
    CLIM_MAX = 2;
    %##
    %- hard code the color limit
    if ~isempty(cluster_match_chars)
        cluster_match_chars(cluster_i)
        switch cluster_i
            case {3,12} % sensorimotor area           
                clim_max = 1.2;
            case {7,9} % posterior area
                clim_max = 0.6;
            case 14 % cingulate
                clim_max = 0.6;
            case {6,13} % supplementary motor
                clim_max = 0.6;
            case {11,10} % occipital
                clim_max = 0.6;
            case 4 % caudate
                clim_max = 0.6;
        end
    else
        clim_max = 1.5;
    end
    %##
    allerspdata_to_use = allersp;
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
    allerspdata_meanSubj = cell(length(allersp),1);
    for j = 1:length(allersp)
        allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata_to_use{j},1),size(allerspdata_to_use{j},2),length(subj_in_cluster));
    end
    allerspdata_remove = allerspdata_to_use;        
    subj_i = 1;
    sub = unique(STUDY.cluster(cluster_i).sets);
    for n = 1:length(unique(STUDY.cluster(cluster_i).sets)) %1:size(allerspdata_meanSubj{1},3) %
        comp_ind =  STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
        for j = 1:size(allersp)
            allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,subj_i:subj_i + length(comp_ind)-1),3);
        end     
        subj_i = subj_i+length(comp_ind);
    end
    %- set titles
    if any(strcmp(STUDY.design(STUDY.currentdesign).variable(1).value,SPEED_REF_CHAR))
        condnames = SPEED_OVERRIDE_CHARS;
    else
        condnames = STUDY.design(STUDY.currentdesign).variable(1).value;
    end
    alltitles = std_figtitle('condnames',condnames, ...
        'clustname', sprintf('CL_%i',cluster_i));
    freqvalues = [4 100];
    baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(5)); 
    freqidx = find(allfreqs>=freqvalues(1) & allfreqs<=freqvalues(2));
    cluster_allcomp_ersp_mean = cell(length(allerspdata_meanSubj),1);
    cluster_allcomp_ersp_crop = cell(length(allerspdata_meanSubj),1);
    for j = 1:length(allerspdata_meanSubj)
        cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j},3);
        cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp_mean{j,1}(freqidx,baseidx);
    end
%     climMat = [min(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all') max(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all')];
%     clim_max = max(abs(climMat))+0.2;
%     std_plottf(alltimes(baseidx),allfreqs(freqidx), cluster_allcomp_ersp_crop, 'datatype','ersp', 'plotmode','normal',...
%         'titles',alltitles,'caxis',[-CLIM_MAX,CLIM_MAX])
%     %- save
%     saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
%     saveas(gcf,[save_dir filesep sprintf('erspplottype3_orig_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    %% (RECOMPUTE STATS) =============================================== %%
    % make allersp to be allersp_crop (contain full gait cycle but not full
    % epoch)
    allerspdata_crop = cell(length(allerspdata_meanSubj),1);
    for j = 1:length(allerspdata_meanSubj)
        allerspdata_crop{j,1} = allerspdata_meanSubj{j,1}(freqidx,baseidx,:);
    end
    [pcond_ersp_crop, ~, ~] = erspStats(STUDY,allerspdata_crop,allfreqs,alltimes);
    %## (PLOT) Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',FIGURE_POSITION,'renderer','Painters');
    for j = 1:length(allersp)
        clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
        clust_maskedersp = clust_ersp;
        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan -clim_max clim_max],...
            'logfreq','native');
        %- adjust panel height
        cur_pos = get(gca,'Position');
        cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
        cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
        cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1);
        cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
        set(gca,'Position',cur_pos)
        %- add vertical line
        xline(gca,warping_times(1),'k--');
        xline(gca,warping_times(2),'k--');
        xline(gca,warping_times(3),'k--');
        xline(gca,warping_times(4),'k--');
        xline(gca,warping_times(5),'k--');
        %- set ylims
        ylim(log(sub_freq_lims))
        if sub_freq_lims(2) <= 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif sub_freq_lims(2) <= 100
            set(gca,'YTick',log([4.01,8,12,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end  
        %- set color lims
        set(gca,'clim',[-CLIM_MAX, CLIM_MAX]);
        %- set y-axis labels
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        %- set x-axis labels
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        %- title
        title(alltitles{j});
    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp{1}),alltimes,allfreqs,'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-clim_max, clim_max]],...
        'logfreq','native');
    %- title
    title('(eeglab) F stats');
    %- set color bar
    c = colorbar();
    c.Position(1) = c.Position(1)+0.05;
    c.Limits = [-CLIM_MAX, CLIM_MAX];
    %- color bar label
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
    hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    set(hL,'Rotation',0);
    %- adjust panel height
    cur_pos = get(gca,'Position');
    cur_pos(4) = cur_pos(4)+PANEL_OFFSET(4);
    cur_pos(3) = cur_pos(3)+PANEL_OFFSET(3);
    cur_pos(1) = cur_pos(1)+PANEL_OFFSET(1)+0.025;
    cur_pos(2) = cur_pos(2)+PANEL_OFFSET(2);
    set(gca,'Position',cur_pos);
    %- add vertical line
    xline(gca,warping_times(1),'k--');
    xline(gca,warping_times(2),'k--');
    xline(gca,warping_times(3),'k--');
    xline(gca,warping_times(4),'k--');
    xline(gca,warping_times(5),'k--');
    %- set ylims
    ylim(log(sub_freq_lims))
    if sub_freq_lims(2) <= 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif sub_freq_lims(2) <= 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    %- set y-axis labels
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    %- set x-axis labels
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',EVENT_CHARS);
    xtickangle(45)
    ax = gca; 
    ax.XAxis.FontSize = 8;
    if ~isnan(alpha)
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_stats_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
    else
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_des%i_cl%i_lim%i.fig',STUDY.currentdesign,cluster_i,sub_freq_lims(2))]);
        saveas(gcf,[save_dir filesep sprintf('erspplottype3_des%i_cl%i_lim%i.jpg',STUDY.currentdesign,cluster_i,sub_freq_lims(2))])
    end
end
