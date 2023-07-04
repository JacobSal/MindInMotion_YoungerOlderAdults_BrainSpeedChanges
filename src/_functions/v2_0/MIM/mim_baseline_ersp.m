function [fig_handle] = mim_plot_ersp(ersp_data,ers_times,ersp_freqs,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Jacob Salminen
% Code Date: 06/09/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- timewarping_values
timewarping_values = [];
%-
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ersp_data',@isstruct);
addRequired(p,'ersp_times',@isstruct);
addRequired(p,'ersp_freqs',@isnumeric);
%## OPTIONAL
addOptional(p,'timewarping_values',timewarping_values,@isnumeric);
addOptional(p,'designnumber',designnumber,@isnumeric);
%## PARAMETER
parse(p,ersp_data,ersp_times,ersp_freqs,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
%## ASSIGN PARAMETERS FOR READING


            ic = STUDY.cluster(k).comps(i);
            sub = STUDY.datasetinfo(STUDY.cluster(k).sets(i)).subject; 

            baseidx = find(alltimes2>=warpingvalues(1) & alltimes2<=warpingvalues(5));
            erspdata = allerspdata2{1}(:,:,i);
            baseline = mean(erspdata(:,baseidx,:),2);
            curr_ersp = erspdata(:,:,:)-repmat(baseline,1,length(alltimes2));

            %Bootstrap and significance mask
%             if ~isnan(alpha)
%                 pboot = bootstat(curr_ersp,'mean(arg1,3);','boottype','shuffle',...
%                 'label','ERSP','bootside','both','naccu',200,...
%                 'basevect',[1:length(baseidx)],'alpha',0.05,'dimaccu',2);         
%                 curr_ersp = mean(curr_ersp,3);
%                 curr_maskedersp = curr_ersp;
%                 curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
%             else
                curr_ersp = mean(curr_ersp,3);
                curr_maskedersp = curr_ersp;
%             end   

            figure('renderer','Painters','print','-bestfit');
            tftopo(curr_maskedersp,alltimes,allfreqs,'limits',... 
                [warpingvalues(1) warpingvalues(end) nan nan nan nan],...
                'vert',warpingvalues(1:5),'logfreq','native');
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
            title(strcat({'Cluster '},num2str(k)));
    %         ylim(log([4 50]))
            cbar('vert');
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
            alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, 'clustname', STUDY.cluster(k).name,...
                'subject', sub, 'compnames', num2str(ic));
            % reorganize allerspdata
            for j = 1: length(allerspdata2)
                erspdata = allerspdata2{j}(:,:,i);
                baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));                
                baseline = mean(erspdata(:,baseidx,:),2);                
                single_comp_ersp{j,1} = mean(allerspdata2{j}(:,:,i)-repmat(baseline,1,length(alltimes)),3);
                single_comp_ersp_crop{j,1} = single_comp_ersp{j,1}(:,baseidx);
            end
            std_plottf(alltimes(baseidx),allfreqs, single_comp_ersp_crop, 'datatype','ersp', 'plotmode','normal','titles',alltitles)
            if saveFig
                if ~exist(fullfile(outputdir,'comps'))
                    mkdir(fullfile(outputdir,'comps'));
                end
                saveas(gcf,fullfile(outputdir,'comps',['Component_'  '_ERSP_' sub '_IC' num2str(ic) , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
            end
            close all
        end
    end
    %keyboard
    %% Summary results
    %% 1. Baseline correction = average of epoch within condition
    if performBaselineCorrect1
        clear allerspdata_meanSubj
        subj_in_cluster = unique(STUDY.cluster(k).sets);%Subjects in this cluster
        for j = 1:4
            allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata2{j},1),size(allerspdata2{j},2),length(subj_in_cluster));
        end
        allerspdata_remove = allerspdata2;        
        p = 1;
        sub = unique(STUDY.cluster(k).sets);
        for n = 1:length(unique(STUDY.cluster(k).sets))    
            comp_ind = STUDY.cluster(k).comps(STUDY.cluster(k).sets == sub(n));
            % comp not using
            if k == 1 %bad comp
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
            'clustname', STUDY.cluster(k).name);
        % reorganize allerspdata
        clear erspdata baseline pboot
        for j = 1: length(allerspdata_meanSubj)
            erspdata = allerspdata_meanSubj{j}(:,:,:);
            baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));                
            baseline_allcomp = mean(erspdata(:,baseidx,:),2); % mean power for each person
            baseline = mean(baseline_allcomp,3);%mean power across participant
            cluster_allcomp_ersp{j,1} = allerspdata_meanSubj{j}(:,:,:)-repmat(baseline_allcomp,1,length(alltimes));% subtract baseline for each person
%             cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline,1,length(alltimes)),3);
            cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp{j,1}(:,baseidx);
        end
        climMat = [min(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all') max(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all')];
        climMat_max = max(abs(climMat));
%         std_plottf(alltimes2(baseidx),allfreqs2, cluster_allcomp_ersp_crop, 'datatype','ersp', 'plotmode','normal',...
%             'titles',alltitles,'caxis',[-climMat_max climMat_max])
%         saveas(gcf,fullfile(outputdir,['Component_'  '_ERSP_CLUSTER' , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
%%      
%         for n = 1:size(cluster_allcomp_ersp{1},3)
%             figure();
%             tftopo(mean(cluster_allcomp_ersp{2}(:,:,n),3),alltimes,allfreqs,'limits',... 
%                 [warpingvalues(1) warpingvalues(end) nan nan nan nan],...
%                 'vert',warpingvalues(1:5),'logfreq','native');
%             set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
%             set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
%             ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
%             xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
%             title(strcat({'Cluster '},num2str(k)));
%     %         ylim(log([4 50]))
%             cbar('vert');
%         end
%%      
        %{
        figure('color','white','position',[200 200 800 250],'renderer','Painters');
        for j = 1:length(allerspdata2)
            if ~isnan(alpha)
    %                 temp = squeeze(mean(cluster_comp_ersp{j,1} ,3));
    %                 temp = squeeze(mean(cluster_allcomp_ersp{j,1}(:,baseidx,:),3));
                curr_ersp_temp = cluster_allcomp_ersp{j,1}(freqidx,baseidx,:);% this is already sub baseline
                curr_ersp_temp_mean = mean(curr_ersp_temp,3);
                % - Old way of performing bootstat 
                %{ 
                temp = cluster_allcomp_ersp{j,1}(:,baseidx,:);
                [pboot Pboottrialstmp, Pboottrials] = bootstat(temp,'mean(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',2000,...
                    'basevect',[1:length(baseidx)],'alpha',0.01,'dimaccu',2);         
                clust_ersp = mean(cluster_allcomp_ersp{j},3);
                clust_maskedersp = clust_ersp;
                clust_maskedersp(clust_ersp > repmat(pboot(:,1),[1 size(clust_ersp,2)]) & clust_ersp < repmat(pboot(:,2),[1 size(clust_ersp,2)])) = 0;
                %}
                % -----------------------------------
                % multiple comparison correction % problem with using fdr,
                % inconsistent threshold use across conditions, therefore,
                % it is probably better to set the alpha to be very
                % conservative - 20230416 note
%                 maskersp = [];
%                 maskitc  = [];                
% %                 if isempty(find(~isnan(pboot))) % if ERSP lims not provided
%                     if ndims(Pboottrials) < 3, Pboottrials = Pboottrials'; end
%                     exactp_ersp = compute_pvals(mean(temp,3), Pboottrials);
%                     alphafdr = fdr(exactp_ersp, 0.05);
%                     if alphafdr ~= 0
%                         fprintf('\n ERSP correction for multiple comparisons using FDR, alpha_fdr = %3.6f\n', alphafdr);
%                     else fprintf('\n ERSP correction for multiple comparisons using FDR, nothing significant\n', alphafdr);
%                     end
%                     maskersp = exactp_ersp >= alphafdr;  
%                     clust_maskedersp = clust_ersp;
%                     clust_maskedersp(maskersp) = 0;
% %                 end   
            
                % - New way of performing bootstat based on Amanda and Makoto's
                % code - 20230508 from https://sccn.ucsd.edu/pipermail/eeglablist/2014/008354.html
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
%                 p_masked = pvalMap ;p_masked(pvalMap > 0.05) = 0;p_masked(pvalMap < 0.05) = 1;
                [p_masked, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvalMap,0.05,'pdep',1);

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

                subplot(1,length(allerspdata2)+1,j) % add one subplot for stats
                tftopo(clust_maskedersp,alltimes(baseidx),allfreqs(freqidx),'limits',... 
                    [warpingvalues(1) warpingvalues(end) nan nan -climMat_max climMat_max],...
                    'vert',warpingvalues(1:4),'logfreq','native');
        %             ylim(log([3 50]));
                ylim(log([4 YlimMax]))
                if YlimMax == 50
                    set(gca,'YTick',log([4.01,8,13,30,50])); 
                    set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
                elseif YlimMax == 100
                    set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
                    set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
                end   
                if j == 1
                    ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
                else
                    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
                end
                xlabel('Time (ms)','Fontsize',10);
                title(strcat({'Cluster '},num2str(k),'-',terrain_keyword{j}));
    %             colormap(othercolor('RdBu11'));
        end
%         cbar('vert',1:64,[-climMat_max, climMat_max])
    %         saveas(gcf,fullfile(outputdir,['STATS_'  '_ERSP_CLUSTER' , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
    
        % --- add stats plot ------------
        subplot(1,length(allerspdata2)+1,5) % add one subplot for stats
        tftopo(double(pcon2{1}),alltimes2,allfreqs2,'limits',... 
            [warpingvalues(1) warpingvalues(end) nan nan  [-climMat_max, climMat_max]],...
            'vert',warpingvalues(1:5),'logfreq','native');
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
        title(strcat({'F stats Cluster '},num2str(k)));
        colormap(colormap_ersp)
    %             cbar('vert',1:64,[-climMat_max, climMat_max])
        hp4 = get(subplot(1,length(allerspdata2)+1,5),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = [-climMat_max, climMat_max];
%         hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
%             'bold','FontName','Arial','FontSize',9);
%         set(hL,'Rotation',0);
%         hL.Position(1) = hL.Position(1)+1.7;
    %             hL.Position(2) = hL.Position(2)+0.025;
%         hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;
        
        if ~isnan(alpha)
            saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_eeglab_1_',num2str(k),'_stats_',num2str(YlimMax),'.fig']));
            saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_eeglab_1_',num2str(k),'_stats_',num2str(YlimMax),'.jpg']));
        else
            saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_eeglab_1_',num2str(k),num2str(YlimMax),'.fig']));
            saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_eeglab_1_',num2str(k),num2str(YlimMax),'.jpg']));
        end
        %}
        %% Paper Figure for YA paper and IEEE NER - significance masked ERSP for high terrain
        freqvalues = [4 YlimMax];
        freqidx = find(allfreqs3>=freqvalues(1) & allfreqs3<=freqvalues(2));
        figure('color','white','position',[200 200 700 150],'renderer','Painters');
        for j = 1:length(allerspdata2)%1:length(allerspdata2)
            if ~isnan(alpha)
    %                 temp = squeeze(mean(cluster_comp_ersp{j,1} ,3));
    %                 temp = squeeze(mean(cluster_allcomp_ersp{j,1}(:,baseidx,:),3));
                %{
                temp = cluster_allcomp_ersp{j,1}(:,baseidx,:);
                [pboot Pboottrialstmp, Pboottrials] = bootstat(temp,'mean(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',2000,...
                    'basevect',[1:length(baseidx)],'alpha',0.01,'dimaccu',2);         
                clust_ersp = mean(cluster_allcomp_ersp{j}(:,:,:),3);
                clust_maskedersp = clust_ersp;
                clust_maskedersp(clust_ersp > repmat(pboot(:,1),[1 size(clust_ersp,2)]) & clust_ersp < repmat(pboot(:,2),[1 size(clust_ersp,2)])) = 0;
                %}
                % -----------------------------------
                % multiple comparison correction 
%                 maskersp = [];
%                 maskitc  = [];                
% %                 if isempty(find(~isnan(pboot))) % if ERSP lims not provided
%                     if ndims(Pboottrials) < 3, Pboottrials = Pboottrials'; end
%                     exactp_ersp = compute_pvals(mean(temp,3), Pboottrials);
%                     alphafdr = fdr(exactp_ersp, 0.05);
%                     if alphafdr ~= 0
%                         fprintf('\n ERSP correction for multiple comparisons using FDR, alpha_fdr = %3.6f\n', alphafdr);
%                     else fprintf('\n ERSP correction for multiple comparisons using FDR, nothing significant\n', alphafdr);
%                     end
%                     maskersp = exactp_ersp >= alphafdr;  
%                     clust_maskedersp = clust_ersp;
%                     clust_maskedersp(maskersp) = 0;
%                 end   
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
            faceAlpha_mask = ones(size(clust_maskedersp))*0.8; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(clust_maskedersp ~=0 ) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes(baseidx), allfreqs(freqidx), clust_ersp,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes(baseidx),allfreqs(freqidx),clust_maskedersp,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',[-climMat_max, climMat_max],'xlim',[warpingvalues(1) warpingvalues(end)],...
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
            set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
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
            
%         saveas(gcf,fullfile( summary_figure_folder,['allerspdata_within_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
%         saveas(gcf,fullfile( summary_figure_folder,['allerspdata_within_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.jpg']));
%         saveas(gcf,fullfile(output_folder,['allerspdata_within_high_eeglab_',num2str(k),'_stats.jpg']));
         exportgraphics(gcf,fullfile(
    end
    %% Summary results - 
  if performBaselineCorrect2
    %% 2. Common Baseline correction = average of epoch across condition (commonbaseline)
    % Read allerspdata3, allerspdata3 should already be baselined
    % allerspdata, 
    clear baseline_allcomp baseline_allcond_mean baseline_allcond_median baseline_allcond_time baseline_allcond
    subj_in_cluster = unique(STUDY.cluster(k).sets);%Subjects in this cluster
    
    allerspdata_to_use = allerspdata3;
    pcon_to_use = pcon3;
    
    clear allerspdata_meanSubj
    subj_in_cluster = unique(STUDY.cluster(k).sets);%Subjects in this cluster
    for j = 1:4
        allerspdata_meanSubj{j,1}(:,:,:) = zeros(size(allerspdata_to_use{j},1),size(allerspdata_to_use{j},2),length(subj_in_cluster));
    end
    allerspdata_remove = allerspdata_to_use;        
    p = 1;
    sub = unique(STUDY.cluster(k).sets);
    for n = 1:length(unique(STUDY.cluster(k).sets))    
        comp_ind =  STUDY.cluster(k).comps(STUDY.cluster(k).sets == sub(n));
        for j = 1:4
            allerspdata_meanSubj{j}(:,:,n) = nanmean( allerspdata_remove{j}(:,:,p:p + length(comp_ind)-1),3);
        end     
        p = p+length(comp_ind);
    end
    alltitles = std_figtitle('condnames',STUDY.design(STUDY.currentdesign). variable(1). value, ...
        'clustname', STUDY.cluster(k).name);    
    
    % reorganize allerspdata
%     clear erspdata baseline pboot
%     for j = 1: length(allerspdata_meanSubj)
%         erspdata = allerspdata_meanSubj{j}(:,:,:);
%         baseidx = find(alltimes3>=warpingvalues(1) & alltimes3<=warpingvalues(5));                
%         baseline_allcomp = erspdata(:,baseidx,:); 
%         for sub = 1:size(baseline_allcomp,3) % each participant
%             baseline_allcond{sub}(:,:,j) = baseline_allcomp(:,:,sub);
%         end
%     end
%     for sub = 1:size(baseline_allcomp,3)
%         baseline_allcond_mean{sub} = mean(mean(baseline_allcond{sub},2),3);
%         baseline_allcond_median{sub} = median(baseline_allcond{sub},3);
%         baseline_allcond_time(:,:,sub) = mean(mean(baseline_allcond{sub},2),3);
%     end
% %     baseline = mean(baseline_allcomp,3);
%     for j = 1:length(allerspdata_meanSubj)
%         cluster_allcomp_ersp_crop{j,1} = allerspdata_meanSubj{j}(:,:,:)-squeeze(repmat(baseline_allcond_time,1,length(alltimes3)));
%         cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline_allcond_time,1,length(alltimes3)),3);
%         cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp_mean{j,1}(:,baseidx);
%     end
% %     climMat = [min(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all') max(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all')];
% %     climMat_max = max(abs(climMat))+0.2;
    freqvalues = [4 YlimMax];
    baseidx = find(alltimes3>=warpingvalues(1) & alltimes3<=warpingvalues(5)); 
    freqidx = find(allfreqs3>=freqvalues(1) & allfreqs3<=freqvalues(2));
    for j = 1:length(allerspdata_meanSubj)
        cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j},3);
        cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp_mean{j,1}(freqidx,baseidx);
    end
    climMat = [min(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all') max(cluster_allcomp_ersp_crop{4}(1:30,:),[],'all')];
    climMat_max = max(abs(climMat))+0.2;
%     std_plottf(alltimes3(baseidx),allfreqs3(freqidx), cluster_allcomp_ersp_crop, 'datatype','ersp', 'plotmode','normal',...
%         'titles',alltitles,'caxis',[-1.5 1.5])
%         saveas(gcf,fullfile(outputdir,['Component_'  '_ERSP_CLUSTER' , ' Design', file_keyword,num2str(STUDY.currentdesign)]));
    % ----------------------------------------
    % Recompute stats
    % -----------------------------------------
    % make allersp to be allersp_crop (contain full gait cycle but not full
    % epoch)
    for j = 1:length(allerspdata_meanSubj)
        allerspdata_crop{j,1} = allerspdata_meanSubj{j,1}(freqidx,baseidx,:);
    end
%     [pcond_ersp_crop, pgroup_ersp_crop, pinter_ersp_crop] = erspStats(STUDY,allerspdata_crop);
    
%%  Figure: full tfplot and set limits, stats results use the precomputed results  
    figure('color','white','position',[200 200 700 150],'renderer','Painters');
    for j = 1:length(allerspdata3)
        clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j)
        tftopo(clust_maskedersp,alltimes3,allfreqs3,'limits',... 
            [warpingvalues(1) warpingvalues(end) nan nan -climMat_max climMat_max],...
            'vert',warpingvalues(1:5),'logfreq','native');
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
        set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        title(title_keyword{j});   
    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    tftopo(double(pcon_to_use{1}),alltimes3,allfreqs3,'limits',... 
        [warpingvalues(1) warpingvalues(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warpingvalues(1:5),'logfreq','native');
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
    set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 8;

%             cbar('vert',1:64,[-climMat_max, climMat_max])
    hp4 = get(subplot(1,length(allerspdata3)+1,5),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',0);
    hL.Position(1) = hL.Position(1)+1.2;
            hL.Position(2) = hL.Position(2)+0.025;
    hL.Position(2) = .13;
    c.Position(2) = .13;
    c.Position(4) = .73;
    colormap(colormap_ersp);
    
    if ~isnan(alpha)
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
    else
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),num2str(YlimMax),file_keyword,'.jpg']));
    end
    
    %%  Figure: tfplot with preset range
    figure('color','white','position',[200 200 850 150],'renderer','Painters');
    for j = 1:length(allerspdata_to_use)
        clust_ersp = mean(cluster_allcomp_ersp_crop{j},3);
        clust_maskedersp = clust_ersp;

        subplot(1,length(allerspdata_to_use)+1,j);hold on;
        contourf(alltimes3(baseidx),allfreqs3(freqidx),clust_maskedersp,200,...
         'linecolor','none');hold on;
            %- add vertical line
        vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
        vline(warpingvalues(1),'k-');
        hline(50,'k-');
        set(gca,'xlim',[warpingvalues(1) warpingvalues(end)],...
            'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
%             ylim(log([3 50]));    
        if YlimMax == 50
            set(gca,'YTick',[4,8,13,30,50]); 
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
        xlim([warpingvalues(1) warpingvalues(end)]);
        set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        title(title_keyword{j});   
        box on;
    end
    % --- add stats plot ------------
    subplot(1,length(allerspdata_to_use)+1,5) % add one subplot for stats
    contourf(alltimes3(baseidx),allfreqs3(freqidx),double(pcond_ersp_crop{1}),200,...
         'linecolor','none');     
%             ylim(log([3 50]));
    vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
    set(gca,'xlim',[warpingvalues(1) warpingvalues(end)],...
        'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
    set(gca,'clim',[-climMat_max, climMat_max]);
    if YlimMax == 50
        set(gca,'YTick',[4,8,13,30,50]); 
        set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        set(gca,'linewidth',1);
    elseif YlimMax == 100
        set(gca,'YTick',[4,8,13,30,50,99.4843]); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
    end    
    xlabel('','Fontsize',10);
    xlim([warpingvalues(1) warpingvalues(end)]);
    set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 8;
    colormap(colormap_ersp);
    ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
    title(strcat({'F stats '}));
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
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_stats_notfull_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_stats_notfull_',num2str(YlimMax),file_keyword,'.jpg']));
    else
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_notfull_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_eeglab_2_',num2str(k),'_notfull_',num2str(YlimMax),file_keyword,'.jpg']));
    end
    %% Pairwise comparison
    % compute difference ersps
    if runPairwise
        clear erspDiff_wind erspDiff
        if STUDY.currentdesign == 1
            refErspCond = 'flat';
            refErspCond_ind = strmatch(refErspCond,[STUDY.design(1).variable(1).value]);
        else
            refErspCond = '1p0';
            refErspCond_ind = strmatch(refErspCond,[STUDY.design(2).variable(1).value]);
        end
        if ~isempty(refErspCond)
            r = refErspCond_ind; %reference ersp index
            ind = 1:length(allerspdata3);
            ind = setdiff(ind,r);
        end
        
        if isempty(refErspCond_ind)
            error('Condition for reference ersp not found in STUDY design: %s',refErspCond)
        end
        %{
        if ~isempty(refErspCond)
            r = refErspCond_ind; %reference ersp index
            ind = 1:length(allerspdata3);
            ind= setdiff(ind,r);
            % mask differenec ersps- check that it's sig. different from zero
            clear erspDiff
            for c = ind
                clear pcond
                curr_ersp = allerspdata_to_use{c,1};
                ref_ersp = allerspdata_to_use{r,1};
                %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
%                 [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp;ref_ersp});
%                 [erspDiff(c).raw] = [mean(curr_ersp-ref_ersp,3)];
%                 [erspDiff(c).masked] = [erspDiff(c).raw.*pcond{1,1}];
%                 [erspDiff(c).pcond] = pcond{1,1};
                
                clear pcond
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
            mkdir(fullfile(stats_output_folder,['cluster_',num2str(k)]));
            save(fullfile(stats_output_folder,['cluster_',num2str(k)],['Common_baseline_2',file_keyword,'.mat']),...
                'pcond_ersp_crop','pgroup_ersp_crop','pinter_ersp_crop','erspDiff','erspDiff_wind');
        end
        %}
        %% Figure Full comparing high terrain with flat terrain
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
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warpingvalues(1) warpingvalues(end)],...
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
            set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(title_keyword{ind(j)});   
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
        c.Limits = clim;
        hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
            'bold','FontName','Arial','FontSize',9);
        set(hL,'Rotation',90);
        hL.Position(1) = hL.Position(1)+0.5;
    %             hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = .13;
        c.Position(2) = .13;
        c.Position(4) = .73;    
    %     colorbar

%         saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
%         saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
%         saveas(gcf,fullfile(output_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_stats.jpg']));
        %% Figure not full window
        if ~isempty(refErspCond)
            data = [];
            for j = 1:length(ind)
                data = [data, reshape(mean(erspDiff_wind(ind(j)).raw,3).',1,[])];
            end
            IQR = iqr(data); %interquartile range
            Q1 = quantile(data,0.25);
    %         myMin = round(Q1-1.5*IQR,1);
            myMin = -round(mean(data)+2*IQR,1);
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
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',clim,'xlim',[warpingvalues(1) warpingvalues(end)],...
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
            set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            title(title_keyword{ind(j)});   
        end
        hp4 = get(subplot(1,length(allerspdata_to_use),3),'Position');
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

        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_notfull_stats_',num2str(YlimMax),file_keyword,'.fig']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_notfull_stats_',num2str(YlimMax),file_keyword,'.pdf']));
        saveas(gcf,fullfile( summary_figure_folder,['allerspdata3_terrain_v_flat_eeglab_',num2str(k),'_notfull_stats_',num2str(YlimMax),file_keyword,'.jpg']));

    
  end
  end
end
%## TIME
toc
end

