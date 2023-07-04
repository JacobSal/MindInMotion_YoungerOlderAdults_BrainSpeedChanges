%% TESTING ICA vs. PCA CLUSTERING
[tmp_study,tmp_alleeg,comps_out] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
[tmp_study2] = cluster_ica_reduce(MAIN_STUDY);
%##COMPARISON PLOTS
std_dipplot(tmp_study2,MAIN_ALLEEG,'clusters',2:length(tmp_study2.cluster),'mode','multicolor');
std_dipplot(tmp_study,tmp_alleeg,'clusters',2:length(tmp_study.cluster),'mode','multicolor');
% view([45,-45,0])
% view([0,-45,0])
% view([0,0,45])
%- Spec plot ICA
specMin = 15;
specMax = 40;
std_specplot(tmp_study2,MAIN_ALLEEG,'clusters',1:length(tmp_study2.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Spec plot PCA
std_specplot(tmp_study,tmp_alleeg,'clusters',1:length(tmp_study.cluster),'ylim',[specMin,specMax],'freqrange',[1,65]);
%- Topo plot ICA
std_topoplot(tmp_study2,MAIN_ALLEEG,'clusters', 'all');
%- Topo plot PCA
std_topoplot(tmp_study,tmp_alleeg,'clusters', 'all');
%## POP VIEW PROPS (ICA)
if ~exist([save_dir filesep 'component_props'],'dir')
    mkdir([save_dir filesep 'component_props']);
end
for cluster_i = 2:length(tmp_study2.cluster)
    sets_clust = tmp_study2.cluster(cluster_i).sets;
    for i = 1:length(sets_clust)
        subj_i = sets_clust(i);
        comps_clust = tmp_study2.cluster(cluster_i).comps(i);
        hold on;
        pop_prop_extended(MAIN_ALLEEG(subj_i),0,comps_clust,NaN,...
        {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
        fig = gcf;
        fig.Position = [500 300 920 480]; 
        hold off;
%         savefig(fig,[save_dir filesep 'component_props' filesep sprintf('%i_%s_viewprops_co%i_cl%i.fig',cluster_i,MAIN_ALLEEG(subj_i).subject,comps_clust,cluster_i)]);
    end
end
%## POP VIEW PROPS (PCA)
if ~exist([save_dir filesep 'component_props'],'dir')
    mkdir([save_dir filesep 'component_props']);
end
for cluster_i = 2:length(tmp_study.cluster)
    sets_clust = tmp_study.cluster(cluster_i).sets;
    for i = 1:length(sets_clust)
        subj_i = sets_clust(i);
        comps_clust = tmp_study.cluster(cluster_i).comps(i);
        hold on;
        pop_prop_extended(tmp_alleeg(subj_i),0,comps_clust,NaN,...
        {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]});
        fig = gcf;
        fig.Position = [500 300 920 480]; 
        hold off;
%         savefig(fig,[save_dir filesep 'component_props' filesep sprintf('%i_%s_viewprops_co%i_cl%i.fig',cluster_i,tmp_alleeg(subj_i).subject,comps_clust,cluster_i)]);
    end
end