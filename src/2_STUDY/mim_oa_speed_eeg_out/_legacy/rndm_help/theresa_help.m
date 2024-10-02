

TIME_BOUND = [-500,500];
TIME_BOUND_BASE = [-500,0];
FREQ_BOUND = [3,60];
time_inds = alltimes >= TIME_BOUND(1) & alltimes <= TIME_BOUND(2);
freq_inds = allfreqs >= FREQ_BOUND(1) & allfreqs <= FREQ_BOUND(2);
% sub_times = alltimes(time_inds);
%##
tmp_dat = load('C:\Users\jsalminen\Downloads\readESRP_1_Single_subBase.mat');
allersp = tmp_dat.allerspdata2;
allfreqs = tmp_dat.allfreqs2;
alltimes = tmp_dat.alltimes2;
% tmp_dat = load('C:\Users\jsalminen\Downloads\readESRP_1_Single.mat');
% allersp = tmp_dat.allerspdata;
% allfreqs = tmp_dat.allfreqs;
% alltimes = tmp_dat.alltimes;
%-
tmp = load('-mat','C:\Users\jsalminen\Downloads\STUDY-preprocess-P3_20240125_S_SubEpoch_numTrials.study');
STUDY = tmp.STUDY;
STUDY = std_makedesign(STUDY,[],1,'subjselect',{},...
    'variable1','cond','values1', {'Pong','TableTennis'}); %updating STUDY to make sure it is clear that there are condition values
STUDY = pop_statparams(STUDY,'condstats','on',...
    'method','perm',...
    'singletrials','off',...
    'mode','fieldtrip',...
    'fieldtripalpha',0.05,...
    'fieldtripmethod','montecarlo',...
    'fieldtripmcorrect','cluster',...
    'fieldtripnaccu',2000);
STUDY = pop_erspparams(STUDY, 'subbaseline','off',...
    'timerange',[TIME_BOUND(1),TIME_BOUND(2)],...
    'ersplim', [-1,1]);
%%
% [allersp_b1_f,allersp_b1,~,~] = eeglab_baseln(allersp,alltimes,allfreqs,...
%             TIME_BOUND_1,FREQ_BOUND,...
%             'DO_COMMON_BASE',false,...
%             'DO_SUBJ_BASE',true);
[allersp_b1_f,allersp_b1,~,~] = eeglab_baseln(allersp,alltimes,allfreqs,...
            TIME_BOUND_BASE,FREQ_BOUND,...
            'DO_COMMON_BASE',false,...
            'DO_SUBJ_BASE',true);
for cond_i = 1:length(allersp)
    allersp_b1_f{cond_i} = allersp_b1_f{cond_i}(freq_inds,time_inds,:);
end
[pcond_ersp_crop, ~, ~] = ersp_stats_conds(STUDY,allersp,allfreqs,alltimes); 
[pcond_erspbase_crop, ~, ~] = ersp_stats_conds(STUDY,allersp_b1_f,allfreqs(freq_inds),alltimes(time_inds)); 

%%
% p1 = cell(length(allersp),1);
p2 = cell(length(allersp),1);
p3 = cell(length(allersp),1);
for cond_i = 1:length(allersp)
    % p1{cond_i} = squeeze(mean(allersp_b1{cond_i},3));
    p2{cond_i} = squeeze(mean(allersp{cond_i},3));
    p3{cond_i} = squeeze(mean(allersp_b1_f{cond_i},3));
end
%##
PLOT_STRUCT = struct('figure_position_inch',[0.5,0.5,6.5,9],...
    'alltitles',{{}},...
    'xlabel','Gait Cycle Time (ms)',...
    'ylabel','Frequency (Hz)',...
    'xticklabel_times',{[]},...
    'xticklabel_chars',{{}},...
    'clim',[],...
    'font_size',[],...
    'font_name','Arial',...
    'freq_lims',FREQ_BOUND,...
    'time_lims',TIME_BOUND,...
    'subplot_width',[],...
    'subplot_height',[],... %(02/17/2024) was 0.2
    'horiz_shift_amnt',[],...
    'vert_shift_amnt',[],...
    'group_titles',{{}},...
    'stats_title','F Statistic Mask',...
    'figure_title','');
PLOT_STRUCT.alltitles = {'1','2','none'};
PLOT_STRUCT.clim = [-0.25,0.25];
% fig = plot_txf_conds_tftopo(p1,sub_times,allfreqs,{},{},...
%     'PLOT_STRUCT',PLOT_STRUCT);
% drawnow;
PLOT_STRUCT.figure_title = 'original';
fig = plot_txf_conds_tftopo(p2,alltimes,allfreqs,pcond_ersp_crop,{},...
    'PLOT_STRUCT',PLOT_STRUCT);
drawnow;
PLOT_STRUCT.figure_title = 'baselined';
fig = plot_txf_conds_tftopo(p3,alltimes(time_inds),allfreqs(freq_inds),pcond_erspbase_crop,{},...
    'PLOT_STRUCT',PLOT_STRUCT);
drawnow;
% fig = plot_txf_conds_tftopo(p3,alltimes,allfreqs,pcond_ersp_crop,{},...
%     'PLOT_STRUCT',PLOT_STRUCT);
%%  ===================================================================  %%
TIME_BOUND_BASE = [-1000,0];
BOOT_NITERS = 2000;
BOOT_ALPHA = 0.05;
BOOT_CLUST_THRESH = 100;
[allersp_b1_f,allersp_b1,~,~] = eeglab_baseln(allersp,alltimes,allfreqs,...
            TIME_BOUND_BASE,FREQ_BOUND,...
            'DO_COMMON_BASE',true,...
            'DO_SUBJ_BASE',true);
clust_ersp = cell(size(allersp_b1_f,1),size(allersp_b1_f,2));
clust_maskedersp = cell(size(allersp_b1_f,1),size(allersp_b1_f,2));
group_i = 1;
for cond_i = 1:size(allersp_b1_f,1)
    fprintf('Performing Stats for Condition %i\n',cond_i);
    tmp = allersp_b1_f{cond_i,group_i};
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
    [p_masked, ~, ~, ~] = fdr_bh(pvalMap,BOOT_ALPHA,'pdep',1);
    % debri removal
    [labelMap,~] = bwlabeln(p_masked);
    tmpDisp = sort(labelMap(:),'descend');
%             [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
    [occurrence,idx,~] = histcounts(tmpDisp,unique(tmpDisp));
    kMask = ismember(labelMap,idx((occurrence<BOOT_CLUST_THRESH)));
    finalMask = p_masked-kMask;
    clust_ersp{cond_i} = tmp_mean; 
    tmp = clust_ersp{cond_i}; 
    tmp(~finalMask) = 0;
    clust_maskedersp{cond_i} = tmp;
end
PLOT_STRUCT.subplot_width = 0.13;
PLOT_STRUCT.subplot_height = 0.16;
PLOT_STRUCT.horiz_shift_amnt = 0.17;
PLOT_STRUCT.vert_shift_amnt = 0.22;
PLOT_STRUCT.alltitles = {'Pong','TableTennis'};
PLOT_STRUCT.clim = [-0.25,0.25];
[fig] = plot_txf_mask_contourf(clust_ersp,alltimes,allfreqs,clust_maskedersp,clust_maskedersp,{},...
    'PLOT_STRUCT',PLOT_STRUCT);
drawnow;
% exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_%s_bootstraps_ersp_sb.tiff',cl_i,des_i,groups{groups_ind(group_i)})],'Resolution',1000);
% exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_bootstraps_ersp_sb.jpg',cl_i,des_i)],'Resolution',300);
% close(fig);