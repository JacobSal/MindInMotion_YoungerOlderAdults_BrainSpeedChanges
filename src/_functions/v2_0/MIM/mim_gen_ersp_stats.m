function [] = mim_gen_ersp_stats(STUDY,ALLEEG,warping_times,save_dir,varargin)
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
alpha = 0.05; 
plotComp = 0;
YlimMax = 50;
performBaselineCorrect1 = 0;
performBaselineCorrect2 = 1;
performBaselineCorrect3 = 0;
% saveFig = 1;
saveStats = 1;
runPairwise = 1;
freqvalues = [4 YlimMax];
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'warping_times',@isnumeric);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,ALLEEG,warping_times,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
[~,main_cl_inds,~,valid_cluster] = eeglab_get_cluster_comps(STUDY);
fcn = @erspStats;
%% loop through clusters to plot
% for des_i = 1:length(STUDY.design)
%     for cluster_i = main_cl_inds
des_i = 1;
cluster_i = 3;
%- raw ersp data no baselining
ersp_data = par_load(STUDY.etc.mim_gen_ersp_data.ersp_fpaths{cluster_i,des_i});
allerspdata2 = ersp_data.allerspdata;
alltimes2 = ersp_data.alltimes;
allfreqs2 = ersp_data.allfreqs;
%- normalized ersp data
ersp_data_norm = par_load(STUDY.etc.mim_gen_ersp_data.ersp_norm_fpaths{cluster_i,des_i});
allerspdata3 = ersp_data_norm.allerspdata;
alltimes3 = ersp_data_norm.alltimes;
allfreqs3 = ersp_data_norm.allfreqs;
%- normalized ersp data wtih common baseline
%{
tmp = load(STUDY.etc.mim_gen_ersp_data.ersp_normcb_fpaths{cluster_i,des_i});
ersp_data_normcb = tmp.SAVE_VAR;
allerspdata3 = ersp_data_normcb.allerspdata;
alltimes3 = ersp_data_normcb.alltimes;
allfreqs3 = ersp_data_normcb.allfreqs;
%}

%%
if plotComp % plot each component with baseline to be the average of entire epoch
    alltimes = alltimes2;
    allfreqs = allfreqs2;
    freqidx = find(allfreqs>=freqvalues(1) & allfreqs<=freqvalues(2));
    for i = 1:length(STUDY.cluster(cluster_i).comps)
        ic = STUDY.cluster(cluster_i).comps(i);
        sub = STUDY.datasetinfo(STUDY.cluster(cluster_i).sets(i)).subject; 

        baseidx = find(alltimes2>=warping_times(1) & alltimes2<=warping_times(5));
        erspdata = allerspdata2{1}(:,:,i);
        baseline = mean(erspdata(:,baseidx,:),2);
        curr_ersp = erspdata(:,:,:)-repmat(baseline,1,length(alltimes2));
        curr_ersp = mean(curr_ersp,3);
        curr_maskedersp = curr_ersp;

        figure('renderer','Painters','print','-bestfit');
        tftopo(curr_maskedersp,alltimes,allfreqs,'limits',... 
            [warping_times(1) warping_times(end) nan nan nan nan],...
            'vert',warping_times(1:5),'logfreq','native');
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
        xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
        title(strcat({'Cluster '},num2str(cluster_i)));
        cbar('vert');
        vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
        alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, 'clustname', STUDY.cluster(cluster_i).name,...
            'subject', sub, 'compnames', num2str(ic));
        % reorganize allerspdata
        for j = 1: length(allerspdata2)
            erspdata = allerspdata2{j}(:,:,i);
            baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(5));                
            baseline = mean(erspdata(:,baseidx,:),2);                
            single_comp_ersp{j,1} = mean(allerspdata2{j}(:,:,i)-repmat(baseline,1,length(alltimes)),3);
            single_comp_ersp_crop{j,1} = single_comp_ersp{j,1}(:,baseidx);
        end
        std_plottf(alltimes(baseidx),allfreqs, single_comp_ersp_crop, 'datatype','ersp', 'plotmode','normal','titles',alltitles)
        if saveFig
            if ~exist(fullfile(save_dir,'comps'),'dir')
                mkdir(fullfile(save_dir,'comps'),'dir');
            end
            saveas(gcf,fullfile(save_dir,'comps',['Component_'  '_ERSP_' sub '_IC' num2str(ic) , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
        end
        close all
    end
end
%% 1. Baseline correction = average of epoch within condition
if performBaselineCorrect1
    alltimes = alltimes2;
    allfreqs = allfreqs2;
    freqidx = find(allfreqs>=freqvalues(1) & allfreqs<=freqvalues(2));
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
    for j = 1:4
        allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata2{j},1),size(allerspdata2{j},2),length(subj_in_cluster));
    end
    allerspdata_remove = allerspdata2;        
    p = 1;
    sub = unique(STUDY.cluster(cluster_i).sets);
    for n = 1:length(unique(STUDY.cluster(cluster_i).sets))    
        comp_ind = STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
        % comp not using
        if cluster_i == 1 %bad comp
            for j = 1:4
%                     keyboard
                allerspdata_remove{j}(:,:,p) = nan(size(allerspdata2{j},1),size(allerspdata2{j},2),1);
                allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,p:p + length(comp_ind)-1),3);     
            end
        else
            for j = 1:4
                allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,p:p + length(comp_ind)-1),3);
            end
        end
        p = p+length(comp_ind);
    end
    alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, ...
        'clustname', STUDY.cluster(cluster_i).name);
    % reorganize allerspdata
    clear erspdata baseline pboot
    for j = 1: length(allerspdata_meanSubj)
        erspdata = allerspdata_meanSubj{j}(:,:,:);
        baseidx = find(alltimes>=warping_times(1) & alltimes<=warping_times(5));                
        baseline_allcomp = mean(erspdata(:,baseidx,:),2); % mean power for each person
        baseline = mean(baseline_allcomp,3);%mean power across participant
        cluster_allcomp_ersp{j,1} = allerspdata_meanSubj{j}(:,:,:)-repmat(baseline_allcomp,1,length(alltimes));% subtract baseline for each person
%             cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline,1,length(alltimes)),3);
        cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp{j,1}(:,baseidx);
    end
    climMat = [min(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all') max(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all')];
    climMat_max = max(abs(climMat));

    %## (PLOT) Paper Figure for YA paper and IEEE NER - significance masked ERSP for high terrain
    freqvalues = [4 YlimMax];
    freqidx = find(allfreqs3>=freqvalues(1) & allfreqs3<=freqvalues(2));
    figure('color','white','position',[200 200 700 150],'renderer','Painters');
    for j = 1:length(allerspdata2)%1:length(allerspdata2)
        if ~isnan(alpha)
            curr_ersp_temp = cluster_allcomp_ersp{j,1}(freqidx,baseidx,:);% this is already sub baseline
            curr_ersp_temp_mean = mean(curr_ersp_temp,3);
            surro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
            for n = 1:2000
                bootLatency = randi(size(curr_ersp_temp,2),[size(curr_ersp_temp,2),1]); %random time sample
                bootFreq = 1:size(curr_ersp_temp,1);
                bootIc = 1:size(curr_ersp_temp,3); 
                tmpSurro = mean(curr_ersp_temp(bootFreq,bootLatency,bootIc),3);
                surro(:,:,n) = tmpSurro; %save 2000 iterations of surrogates 
            end
            bootSurro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
            for n = 1:2000
                bootIdx  = randi(2000,[size(curr_ersp_temp,3),1]);
                tmpSurro = mean(surro(:,:,bootIdx),3);
                bootSurro(:,:,n) = tmpSurro;
            end

            pvalMap = stat_surrogate_pvals(bootSurro,curr_ersp_temp_mean,'both');
            pvalMap(pvalMap>1)=1; 
            [p_masked, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalMap,0.05,'pdep',1);

            % debri removal
            [labelMap,uniqueLabelNum] = bwlabeln(p_masked);
            tmpDisp = sort(labelMap(:),'descend');
            [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
            sortOccurrence = sort(occurrence,'descend');
%             disp(num2str(sortOccurrence(2:10)));
            threshold = 1000;
            threshOccurrence = occurrence;
            threshIdx = find(threshOccurrence<threshold);
            kMask = ismember(labelMap,idx(threshIdx));
            finalMask = p_masked-kMask;

            clust_ersp = curr_ersp_temp_mean; 
            clust_maskedersp = clust_ersp; 
            clust_maskedersp(~finalMask) = 0;
%             curr_maskedersp(~p_masked) = 0; 

        else
            clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
            clust_maskedersp = clust_ersp;
        end   

        subplot(1,length(allerspdata2),j)
        colormap(colormap_ersp);
        faceAlpha_mask = ones(size(clust_maskedersp))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        faceAlpha_mask(clust_maskedersp ~=0 ) = 0; %0 is significant? 1 is not? 
    %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
    %             'linecolor','none');hold on;
        contourf(alltimes(baseidx), allfreqs(freqidx), clust_ersp,200,...
                   'linecolor','none');hold on;
        imagesc(alltimes(baseidx),allfreqs(freqidx),clust_maskedersp,'AlphaData',faceAlpha_mask);
        %- add vertical line
        vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
        set(gca,'clim',[-climMat_max, climMat_max],'xlim',[warping_times(1) warping_times(end)],...
            'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
        if YlimMax == 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
        elseif YlimMax == 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
        end
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
        end
        xlabel('','Fontsize',12);
%             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
        title(title_keyword{j});
        set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
    end
    hp4 = get(subplot(1,length(allerspdata2),4),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.008  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',90);
    hL.Position(1) = hL.Position(1)+1.2;
%             hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    c.Position(2) = .10;
    c.Position(4) = .5;    
%         colorbar

    saveas(gcf,fullfile(save_dir,'allerspdata',['allerspdata_within_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
    saveas(gcf,fullfile(save_dir,'allerspdata',['allerspdata_within_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
end
%% 2. Common Baseline correction = average of epoch across condition (commonbaseline)
% Read allerspdata3, allerspdata3 should already be baselined
% allerspdata,  
if performBaselineCorrect2

    clear baseline_allcomp baseline_allcond_mean baseline_allcond_median baseline_allcond_time baseline_allcond
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster

    allerspdata_to_use = allerspdata3;
    pcon_to_use = pcon3;

    clear allerspdata_meanSubj
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
    for j = 1:4
        allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata_to_use{j},1),size(allerspdata_to_use{j},2),length(subj_in_cluster));
    end
    allerspdata_remove = allerspdata_to_use;        
    p = 1;
    sub = unique(STUDY.cluster(cluster_i).sets);
    for n = 1:length(unique(STUDY.cluster(cluster_i).sets))    
        comp_ind =  STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
        for j = 1:4
            allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,p:p + length(comp_ind)-1),3);
        end     
        p = p+length(comp_ind);
    end
    alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, ...
        'clustname', STUDY.cluster(cluster_i).name);    

    freqvalues = [4 YlimMax];
    baseidx = find(alltimes3>=warping_times(1) & alltimes3<=warping_times(5)); 
    freqidx = find(allfreqs3>=freqvalues(1) & allfreqs3<=freqvalues(2));
    for j = 1:length(allerspdata_meanSubj)
        cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j},3);
        cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp_mean{j,1}(freqidx,baseidx);
    end
    climMat = [min(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all') max(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all')];
    climMat_max = max(abs(climMat))+0.2;
    for j = 1:length(allerspdata_meanSubj)
        allerspdata_crop{j,1} = allerspdata_meanSubj{j,1}(freqidx,baseidx,:);
    end
    [pcond_ersp_crop, pgroup_ersp_crop, pinter_ersp_crop] = erspStats(STUDY,allerspdata_crop);

    %##  Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',[200 200 700 150],'renderer','Painters');
    for j = 1:length(allerspdata3)
        clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes3,allfreqs3,'limits',... 
            [warping_times(1) warping_times(end) nan nan -climMat_max climMat_max],...
            'vert',warping_times(1:5),'logfreq','native');
    %             ylim(log([3 50]));    
        ylim(log([4 YlimMax]))
        if YlimMax == 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif YlimMax == 100
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end    
        set(gca,'clim',[-climMat_max, climMat_max]);
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        title(title_keyword{j});   
    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcon_to_use{1}),alltimes3,allfreqs3,'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warping_times(1:5),'logfreq','native');
    %             ylim(log([3 50]));
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    xlabel('','Fontsize',10);
    set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 8;

    %             cbar('vert',1:64,[-climMat_max, climMat_max])
    hp4 = get(subplot(1,length(allerspdata3)+1,5),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
            hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    c.Position(2) = .13;
    c.Position(4) = .73;

    if ~isnan(alpha)
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
    else
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),num2str(YlimMax),file_keyword,'.jpg']));
    end

    %##  Figure: tfplot with preset range
    figure('color','white','position',[200 200 850 150],'renderer','Painters');
    for j = 1:length(allerspdata_to_use)
        clust_ersp = mean(cluster_allcomp_ersp_crop{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes3(baseidx),allfreqs3(freqidx),'limits',... 
            [warping_times(1) warping_times(end) nan nan -climMat_max climMat_max],...
            'vert',warping_times(1:4),'logfreq','native');
    %             ylim(log([3 50]));    
        ylim(log([4 YlimMax]))
        if YlimMax == 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif YlimMax == 100
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end    
        set(gca,'clim',[-climMat_max, climMat_max]);
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        xlabel('','Fontsize',10);
        set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        title(title_keyword{j});   
    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp_crop{1}),alltimes3(baseidx),(freqidx),'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warping_times(1:4),'logfreq','native');
    %             ylim(log([3 50]));
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    

    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    xlabel('Time (ms)','Fontsize',10);
    title(strcat({'F stats Cluster '},num2str(cluster_i)));
    colormap(colormap_ersp);

    %             cbar('vert',1:64,[-climMat_max, climMat_max])
    hp4 = get(subplot(1,length(allerspdata_to_use)+1,5),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
            hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    c.Position(2) = .13;
    c.Position(4) = .73;

    if ~isnan(alpha)
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_stats_notfull_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_stats_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    else
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_notfull_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_eeglab_2_',num2str(cluster_i),'_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    end
    %## Pairwise comparison
    % compute difference ersps
    if runPairwise
        if STUDY.currentdesign == 1
            refErspCond = 'flat';
            refErspCond_ind = strmatch(refErspCond,[STUDY.design(1).variable(1).value]);
        else
            refErspCond = '1p0';
            refErspCond_ind = strmatch(refErspCond,[STUDY.design(2).variable(1).value]);
        end
        if isempty(refErspCond_ind)
            error('Condition for reference ersp not found in STUDY design: %s',refErspCond)
        end
        if ~isempty(refErspCond)
            r = refErspCond_ind; %reference ersp index
            ind = 1:length(allerspdata3);
            ind= setdiff(ind,r);
            % mask differenec ersps- check that it's sig. different from zero
            clear erspDiff
            for c = ind
                curr_ersp = allerspdata_to_use{c,1};
                ref_ersp = allerspdata_to_use{r,1};
                %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
                [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp;ref_ersp});
                [erspDiff(c).raw] = [mean(curr_ersp-ref_ersp,3)];
                [erspDiff(c).masked] = [erspDiff(c).raw.*pcond{1,1}];
                [erspDiff(c).pcond] = pcond{1,1};

                curr_ersp_wind = allerspdata_to_use{c,1}(freqidx,baseidx,:);
                ref_ersp_wind = allerspdata_to_use{r,1}(freqidx,baseidx,:);
                %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
                [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp_wind;ref_ersp_wind});
                [erspDiff_wind(c).raw] = [mean(curr_ersp_wind-ref_ersp_wind,3)];
                [erspDiff_wind(c).masked] = [erspDiff_wind(c).raw.*pcond{1,1}];
                [erspDiff_wind(c).pcond] = pcond{1,1};
            end
        end
        if saveStats
            mkdir(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)]));
            save(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)],['Common_baseline_2',file_keyword,'.mat']),...
                'pcond_ersp_crop','pgroup_ersp_crop','pinter_ersp_crop','erspDiff','erspDiff_wind');
        end
        %## Figure Full comparing high terrain with flat terrain
        % - set clim for erspDiff
        if ~isempty(refErspCond)
            data = [];
            for j = 1:3%1:length(~isempty(erspDiff(:).raw))
                data = [data, reshape(mean(erspDiff(ind(j)).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
            Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+1.5*IQR,1);
            erspDiff_clim = [myMin myMin*(-1)];
        end
        clim = erspDiff_clim;

        f1 = figure('color','white','position',[200 200 700 150],'renderer','Painters');
        for j = 1:length(allerspdata_to_use)-1
            subplot(1,length(allerspdata_to_use),j)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff(ind(j)).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff(ind(j)).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes3, allfreqs3, erspDiff(ind(j)).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes3,allfreqs3,erspDiff(ind(j)).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
            if YlimMax == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif YlimMax == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            end
            xlabel('','Fontsize',10);
            set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(title_keyword{j});   
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = clim;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.2;
    %             hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;    
    %     colorbar

        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_terrain_v_flat_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_terrain_v_flat_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
    %         saveas(gcf,fullfile(output_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_stats.jpg']));
        %## Figure not full window
        if ~isempty(refErspCond)
            data = [];
            for j = 1:3%1:length(erspDiff_wind)
                data = [data, reshape(mean(erspDiff_wind(ind(j)).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
            Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+1.5*IQR,1);
            erspDiff_clim = [myMin myMin*(-1)];
        end
        clim = erspDiff_clim;

        f1 = figure('color','white','position',[200 200 700 150],'renderer','Painters');
        for j = 1:length(allerspdata_to_use)-1
            subplot(1,length(allerspdata_to_use),j)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff_wind(ind(j)).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff_wind(ind(j)).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes3(baseidx), allfreqs3(freqidx), erspDiff_wind(ind(j)).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes3(baseidx),allfreqs3(freqidx),erspDiff_wind(ind(j)).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
            if YlimMax == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif YlimMax == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            end
            xlabel('','Fontsize',10);
            set(gca,'xtick',warping_times,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(title_keyword{j});   
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = clim;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.7;
        hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;

        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata3',['allerspdata3_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.jpg']));

    end
end
%% 3. Single trial full epoch correction + Common Baseline correction = average of epoch across condition (commonbaseline)
    % Read allerspdata4, allerspdata4 should already be baselined
if performBaselineCorrect3
    clear baseline_allcomp baseline_allcond_mean baseline_allcond_median baseline_allcond_time baseline_allcond
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster

    allerspdata_to_use = allerspdata4;
    pcon_to_use = pcon4;

    clear allerspdata_meanSubj
    subj_in_cluster = unique(STUDY.cluster(cluster_i).sets);%Subjects in this cluster
    for j = 1:4
        allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata_to_use{j},1),size(allerspdata_to_use{j},2),length(subj_in_cluster));
    end
    allerspdata_remove = allerspdata_to_use;        
    p = 1;
    sub = unique(STUDY.cluster(cluster_i).sets);
    for n = 1:length(unique(STUDY.cluster(cluster_i).sets))    
        comp_ind =  STUDY.cluster(cluster_i).comps(STUDY.cluster(cluster_i).sets == sub(n));
        for j = 1:4
            allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,p:p + length(comp_ind)-1),3);
        end     
        p = p+length(comp_ind);
    end
    alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, ...
        'clustname', STUDY.cluster(cluster_i).name);    
    freqvalues = [4 100];
    baseidx = find(alltimes3>=warping_times(1) & alltimes3<=warping_times(5)); 
    freqidx = find(allfreqs3>=freqvalues(1) & allfreqs3<=freqvalues(2));
    for j = 1:length(allerspdata_meanSubj)
        cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j},3);
        cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp_mean{j,1}(freqidx,baseidx);
    end
    climMat = [min(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all') max(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all')];
    climMat_max = max(abs(climMat))+0.2;
    std_plottf(alltimes3(baseidx),allfreqs3(freqidx), cluster_allcomp_ersp_crop, 'datatype','ersp', 'plotmode','normal',...
        'titles',alltitles,'caxis',[-1.5 1.5])
    %         saveas(gcf,fullfile(outputdir,['Component_'  '_ERSP_CLUSTER' , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
    % ----------------------------------------
    % Recompute stats
    % -----------------------------------------
    % make allersp to be allersp_crop (contain full gait cycle but not full
    % epoch)
    for j = 1:length(allerspdata_meanSubj)
        allerspdata_crop{j,1} = allerspdata_meanSubj{j,1}(freqidx,baseidx,:);
    end
    [pcond_ersp_crop, pgroup_ersp_crop, pinter_ersp_crop] = erspStats(STUDY,allerspdata_crop);
    %## (PLOT) Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',[200 200 800 250],'renderer','Painters');
    for j = 1:length(allerspdata3)
        clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes3,allfreqs3,'limits',... 
            [warping_times(1) warping_times(end) nan nan -climMat_max climMat_max],...
            'vert',warping_times(1:5),'logfreq','native');
    %             ylim(log([3 50]));    
        ylim(log([4 YlimMax]))
        if YlimMax == 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif YlimMax == 100
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end    
        set(gca,'clim',[-climMat_max, climMat_max]);
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        end
        xlabel('Time (ms)','Fontsize',10);
        title(strcat({'Cluster '},num2str(cluster_i),'-',title_keyword{j}));   

    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcon_to_use{1}),alltimes3,allfreqs3,'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warping_times(1:5),'logfreq','native');
    %             ylim(log([3 50]));
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    xlabel('Time (ms)','Fontsize',10);
    title(strcat({'F stats Cluster '},num2str(cluster_i)));
    colormap(colormap_ersp);

    %             cbar('vert',1:64,[-climMat_max, climMat_max])
    hp4 = get(subplot(1,length(allerspdata3)+1,5),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL.Position(2) = .13;
    c.Position(2) = .13;
    c.Position(4) = .73;

    if ~isnan(alpha)
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.jpg']));
    else
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_',num2str(YlimMax),file_keyword,'.jpg']));
    end

    %## (PLOT) Figure: tfplot with preset range
    figure('color','white','position',[200 200 800 250],'renderer','Painters');
    for j = 1:length(allerspdata_to_use)
        clust_ersp = mean(cluster_allcomp_ersp_crop{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes3(baseidx),allfreqs3(freqidx),'limits',... 
            [warping_times(1) warping_times(end) nan nan -climMat_max climMat_max],...
            'vert',warping_times(1:4),'logfreq','native');
    %             ylim(log([3 50]));    
        ylim(log([4 YlimMax]))
        if YlimMax == 50
            set(gca,'YTick',log([4.01,8,13,30,50])); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif YlimMax == 100
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end    
        set(gca,'clim',[-climMat_max, climMat_max]);
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        end
        xlabel('Time (ms)','Fontsize',10);
        title(strcat({'Cluster '},num2str(cluster_i),'-',title_keyword{j}));   

    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcond_ersp_crop{1}),alltimes3(baseidx),(freqidx),'limits',... 
        [warping_times(1) warping_times(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warping_times(1:4),'logfreq','native');
    %             ylim(log([3 50]));
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    

    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    xlabel('Time (ms)','Fontsize',10);
    title(strcat({'F stats Cluster '},num2str(cluster_i)));
    colormap(colormap_ersp);

    hp4 = get(subplot(1,length(allerspdata_to_use)+1,5),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.7;
            hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    c.Position(2) = .13;
    c.Position(4) = .73;

    if ~isnan(alpha)
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_stats_notfull_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_stats_notfull_',num2str(YlimMax),file_keyword,'.jpg']));
    else
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_notfull_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_eeglab_2_',num2str(cluster_i),'_notfull_',num2str(YlimMax),file_keyword,'.jpg']));
    end
    %## (PLOT) Pairwise comparison
    % compute difference ersps
    if runPairwise
        refErspCond = 'flat';
        refErspCond_ind = find(strmatch(refErspCond,[STUDY.design(1).variable(1).value]));
        if isempty(refErspCond_ind)
            error('Condition for reference ersp not found in STUDY design: %s',refErspCond)
        end
        if ~isempty(refErspCond)
            r = refErspCond_ind; %reference ersp index
            ind = 1:length(allerspdata3);
            ind = setdiff(ind,r);
            % mask differenec ersps- check that it's sig. different from zero
            clear erspDiff
            for c = ind
                curr_ersp = allerspdata_to_use{c,1};
                ref_ersp = allerspdata_to_use{r,1};
                %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
                [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp;ref_ersp});
                [erspDiff(c).raw] = [mean(curr_ersp-ref_ersp,3)];
                [erspDiff(c).masked] = [erspDiff(c).raw.*pcond{1,1}];
                [erspDiff(c).pcond] = pcond{1,1};

                curr_ersp_wind = allerspdata_to_use{c,1}(freqidx,baseidx,:);
                ref_ersp_wind = allerspdata_to_use{r,1}(freqidx,baseidx,:);
                %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
                [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp_wind;ref_ersp_wind});
                [erspDiff_wind(c).raw] = [mean(curr_ersp_wind-ref_ersp_wind,3)];
                [erspDiff_wind(c).masked] = [erspDiff_wind(c).raw.*pcond{1,1}];
                [erspDiff_wind(c).pcond] = pcond{1,1};
            end
        end
        if saveStats
            mkdir(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)]));
            save(fullfile(stats_output_folder,['cluster_',num2str(cluster_i)],['Common_baseline_3',file_keyword,'.mat']),...
                'pcond_ersp_crop','pgroup_ersp_crop','pinter_ersp_crop','erspDiff','erspDiff_wind');
        end
        %## (PLOT) Figure Full comparing high terrain with flat terrain
        % - set clim for erspDiff
        if ~isempty(refErspCond)
            data = [];
            for j = 1:length(erspDiff)
                data = [data, reshape(mean(erspDiff(j).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
            Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+1.5*IQR,1);
            erspDiff_clim = [myMin myMin*(-1)];
        end
        clim = erspDiff_clim;

        f1 = figure('color','white','position',[200 200 800 250],'renderer','Painters');
        for j = 2:length(allerspdata_to_use)
            subplot(1,length(allerspdata_to_use),j)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff(j).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff(j).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes3, allfreqs3, erspDiff(j).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes3,allfreqs3,erspDiff(j).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
            if YlimMax == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif YlimMax == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',10);
            title(strcat({'Cluster '},num2str(cluster_i),' ',title_keyword{j},'-','flat'));  
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = clim;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.7;
    %             hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;    
    %     colorbar

        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_stats_',num2str(YlimMax),file_keyword,'.jpg']));
        %## (PLOT)
        if ~isempty(refErspCond)
            data = [];
            for j = 1:length(erspDiff_wind)
                data = [data, reshape(mean(erspDiff_wind(j).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
            Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+1.5*IQR,1);
            erspDiff_clim = [myMin myMin*(-1)];
        end
        clim = erspDiff_clim;

        f1 = figure('color','white','position',[200 200 800 250],'renderer','Painters');
        for j = 2:length(allerspdata_to_use)
            subplot(1,length(allerspdata_to_use),j)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(erspDiff_wind(j).pcond))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(erspDiff_wind(j).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes3(baseidx), allfreqs3(freqidx), erspDiff_wind(j).raw,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes3(baseidx),allfreqs3(freqidx),erspDiff_wind(j).masked,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warping_times(2) warping_times(3) warping_times(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warping_times(1) warping_times(end)],...
                'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
            if YlimMax == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            elseif YlimMax == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            end    
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',10);
            title(strcat({'Cluster '},num2str(cluster_i),' ',title_keyword{j},'-','flat'));  
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = clim;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+1.7;
    %             hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;    
    %     colorbar

        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.jpg']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile(save_dir,'allerspdata4',['allerspdata4_terrain_v_flat_eeglab_',num2str(cluster_i),'_notfull_stats_',num2str(YlimMax),file_keyword,'.fig']));
    end
end
end
%% ===================================================================== %%
%## SUBFUNCTIONS
% exerpt from std_erspplot 
function [pcond, pgroup, pinter] = erspStats(STUDY,allersp)
    %get stats parameters
    stats = STUDY.etc.statistics;
    stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
    if isempty(STUDY.design(STUDY.currentdesign).variable)
        stats.paired = { };
    else
        stats.paired = { STUDY.design(STUDY.currentdesign).variable(:).pairing };
    end

    %get ersp params
    params = STUDY.etc.erspparams;
    params.plottf =[];
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
