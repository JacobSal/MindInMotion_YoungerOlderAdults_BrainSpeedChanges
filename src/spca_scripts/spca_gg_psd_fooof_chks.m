%   Project Title: MIM YOUNGER AND OLDER ADULTS KINEMATICS-EEG ANALYSIS
%
%   Code Designer: Jacob salminen
%## SBATCH (SLURM KICKOFF SCRIPT)
% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_STUDY/mim_oa_speed_eeg_out/.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
clear java;
%}
%% SET WORKSPACE ======================================================= %%
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic
ADD_CLEANING_SUBMODS = false;
%## Determine Working Directories
if ~ispc
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
        SRC_DIR = fileparts(SRC_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = getenv('SCRIPT_DIR');
        SRC_DIR = getenv('SRC_DIR');
    end
else
    try
        SCRIPT_DIR = matlab.desktop.editor.getActiveFilename;
        SCRIPT_DIR = fileparts(SCRIPT_DIR);
    catch e
        fprintf('ERROR. PWD_DIR couldn''t be set...\n%s',getReport(e))
        SCRIPT_DIR = dir(['.' filesep]);
        SCRIPT_DIR = SCRIPT_DIR(1).folder;
    end
    SRC_DIR = fileparts(SCRIPT_DIR);    
end
%## Add Study, Src, & Script Paths
addpath(SRC_DIR);
cd(SRC_DIR);
fprintf(1,'Current folder: %s\n',SRC_DIR);
%## Set PWD_DIR, EEGLAB path, _functions path, and others...
set_workspace
%% (DATASET INFORMATION) =============================================== %%
[SUBJ_PICS,GROUP_NAMES,SUBJ_ITERS,~,~,~,~] = mim_dataset_information('yaoa_spca');
%% (PARAMETERS) ======================================================== %%
fprintf('Assigning Params\n');
%## Hard Define
%- statistics & conditions
SPEED_CUTOFF = 0.1;
SPEED_VALS = {'0.25','0.5','0.75','1.0';
              '0p25','0p5','0p75','1p0'};
TERRAIN_VALS = {'flat','low','med','high'};
COLOR_MAPS_TERRAIN = linspecer(4);
custom_yellow = [254,223,0]/255;
COLOR_MAPS_TERRAIN = [COLOR_MAPS_TERRAIN(3,:);custom_yellow;COLOR_MAPS_TERRAIN(4,:);COLOR_MAPS_TERRAIN(2,:)];
COLOR_MAPS_SPEED = linspecer(4*3);
COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all');
%## hard define
%- FOOOF
% settings = struct('peak_width_limits',[1,8],...
%     'min_peak_height',0.05,...
%     'max_n_peaks',3);
settings = struct('peak_width_limits',[1,8],...
    'min_peak_height',0.05,...
    'max_n_peaks',5);
f_range = [3, 40];
theta_band_lims = [4, 8];
alpha_band_lims = [8 12];
beta_band_lims  = [12 30];
alpha1_band_lims = [8,10.5];
alpha2_band_lims = [10.5,13];
beta1_band_lims = [13,20];
beta2_band_lims = [20,30];
%% (PATHS) ============================================================= %%
%- datset name
DATA_SET = 'MIM_dataset';
%- study name
STUDY_DNAME = '10172024_MIM_YAOAN89_antsnorm_dipfix_iccREMG0p4_powpow0p3_skull0p01_15mmrej_speed';
studies_fpath = [PATHS.src_dir filesep '_data' filesep DATA_SET filesep '_studies'];
%- load study file
study_fpath = [studies_fpath filesep sprintf('%s',STUDY_DNAME)];
%- load cluster
CLUSTER_K = 11;
CLUSTER_STUDY_NAME = 'temp_study_rejics5';
cluster_fpath = [study_fpath filesep '__iclabel_cluster_kmeansalt_rb5'];
cluster_study_fpath = [cluster_fpath filesep 'icrej_5'];
cluster_k_dir = [cluster_study_fpath filesep sprintf('%i',CLUSTER_K)];
%- save directory 
ANALYSIS_DNAME = 'spca_fooof_psd_anl';
save_dir = [cluster_k_dir filesep ANALYSIS_DNAME];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% LOAD STUDY
if ~ispc
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s_UNIX.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
else
    tmp = load('-mat',[cluster_study_fpath filesep sprintf('%s.study',CLUSTER_STUDY_NAME)]);
    STUDY = tmp.STUDY;
end
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '_UNIX.study'],'filepath',cluster_study_fpath);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[CLUSTER_STUDY_NAME '.study'],'filepath',cluster_study_fpath);
% end
%-
cl_struct = par_load(cluster_k_dir,sprintf('cl_inf_%i.mat',CLUSTER_K));
STUDY.cluster = cl_struct;
[comps_out,main_cl_inds,outlier_cl_inds] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICS = main_cl_inds;
%%
%## RE-POP PARAMS
STUDY_DESI_PARAMS = {{'subjselect',{},...
            'variable2','cond','values2',{'flat','low','med','high'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}},...
            {'subjselect',{},...
            'variable2','cond','values2',{'0p25','0p5','0p75','1p0'},...
            'variable1','group','values1',{'H1000''s','H2000''s','H3000''s'}}};
STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY.cache = [];
for des_i = 1:length(STUDY_DESI_PARAMS)
    [STUDY] = std_makedesign(STUDY,ALLEEG,des_i,STUDY_DESI_PARAMS{des_i}{:});
end
%% LOAD FOOOF DATA
tmp = load([save_dir filesep 'psd_feature_table.mat']);
FOOOF_TABLE = tmp.FOOOF_TABLE;
%% MULTI-CLUSTER PLOT OF ALL SUBJECTS ================================== %%
%{
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
conditions = unique(FOOOF_TABLE.cond_char);
groups = unique(FOOOF_TABLE.group_char);
%-
% des_i = 1;
% cond_i = 1;
% group_i = 1;
IM_RESIZE = 0.5;
HZ_DIM = 4;
VERTICAL_SHIFT =  0.2;
HORIZONTAL_SHIFT = 0.2;
HORIZ_START = 0.08;
VERTICAL_START = 0.75;
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
%##
for j = 1:length(designs)
    des_i = double(string(designs(j)));
    %-
    switch des_i
        case 1
            color_dark = COLOR_MAPS_TERRAIN;
            color_light = COLOR_MAPS_TERRAIN;
            GROUP_CMAP_OFFSET = [0,0.1,0.1];
            xtick_label_g = {'flat','low','med','high'};
        case 2
            color_dark = COLOR_MAPS_SPEED;
            color_light = COLOR_MAPS_SPEED+0.15;
            GROUP_CMAP_OFFSET = [0.15,0,0];
            xtick_label_g = {'0.25','0.50','0.75','1.0'};
    end
    for cond_i = 1:length(xtick_label_g)
        for group_i = 1:length(groups)
            fig = figure('color','white','renderer','Painters');
            sgtitle(sprintf('Design %i, Condition %s, Group %s',des_i,xtick_label_g{cond_i},groups(group_i)),'FontName','Arial','FontSize',14,'FontWeight','bold','Interpreter','none');
            set(fig,'Units','inches','Position',[0.5,0.5,6,9.5])
            set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
            hold on;
            set(gca,AXES_DEFAULT_PROPS{:})
            
            %-
            horiz_shift = HORIZ_START;
            hz = 0;
            vert_shift = 0;
            for i = 1:length(clusters)
                hz = hz + 1;
                cl_i = double(string(clusters(i)));
                axes();
                hold on;
                %-
                fooof_psd = fooof_diff_store{des_i}{cl_i}{cond_i,group_i}';
                fooof_psd_mean = mean(fooof_psd);
                subjs = plot(fooof_freq,fooof_psd,'color',[0,0,0,0.15],'linestyle','-','linewidth',2,'displayname','orig. subj psd');
                mean_plot = plot(fooof_freq,fooof_psd_mean,'color',color_dark(cond_i,:),'linestyle','-','linewidth',4,'displayname','orig. subj psd');
                
                %-
                ax = gca;
                plot([0 40],[0 0],'--','color','black');
                xlim([4 40]);
                ylim([-2 10]);
                xlabel('Frequency(Hz)');
                if hz ~= 1
                    ylabel('');
                else
                    ylabel('10*log_{10}(Power)');
                end
                set(ax,'FontName','Arial',...
                    'FontSize',12,...
                    'FontWeight','bold')
                xline(3,'--'); xline(8,'--'); xline(13,'--'); xline(30,'--');
                set(ax,'FontName','Arial','FontSize',10,...
                    'FontWeight','bold')
                title(sprintf('CL%i',cl_i))
                set(ax,'OuterPosition',[0,0,1,1]);
                set(ax,'Position',[horiz_shift,VERTICAL_START-vert_shift,0.3*IM_RESIZE,0.25*IM_RESIZE]);  %[left bottom width height]
                if hz < HZ_DIM
                    horiz_shift = horiz_shift + HORIZONTAL_SHIFT;
                    
                else
                    vert_shift = vert_shift + VERTICAL_SHIFT;
                    horiz_shift = HORIZ_START;
                    hz = 0;
                end
            end
            hold off;
            exportgraphics(fig,[save_dir filesep sprintf('cl%i_des%i_cond%i_group%i_cluster_fooofs_subjs_means.tiff',cl_i,des_i,cond_i,group_i)],'Resolution',600)
            close(fig);
        end
    end
end
%}
%% SANITY CHECK: APERIODIC EXP., TERRAIN WALKING SPEED, APERIODIC OFFSET
%{
%## (STATS STRUCT) ====================================================== %%
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'measure_tag',categorical({''}),...
    'design_tag',categorical({''}),...
    'mod_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{''},...
    'anova_preds_p',{[]},...
    'anova_preds_stat',{[]},...
    'anova_preds_df',{[]},...
    'mod_preds_terms',{{''}},...
    'mod_preds_p',[],...
    'mod_preds_stat',[],...
    'mod_preds_coeff',[],...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',[],...
    'multi_comp_p',[],...
    'multi_comp_coeff',[],...
    'multi_comp_lci_uci',[],...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_struct = DEF_STATS_TRACK_STRUCT;
cnts = 1;
%##
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
g_chars_subp = {'YA','OHMA','OLMA'};
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
FIGURE_POSITION =[1,0,6.5,9];
IM_RESIZE = 0.85;
AX_H  = 0.2;
AX_W = 0.275;
AX_HORIZ_SHIFT = 0.06;
AX_VERT_SHIFT = 0.11;
AX_INIT_HORIZ = 0.06;
AX_INIT_VERT_VIO = 0.8;
HZ_DIM = 3;
MEASURES_PLOT = {'aperiodic_exp','speed_ms','aperiodic_offset'};
% measure_plot = 'aperiodic_exp';
%-
for k = 1:length(MEASURES_PLOT)
    measure_plot = MEASURES_PLOT{k};
    for i = 1:length(designs)
        des_i = double(string(designs(i)));
        switch des_i
            case 1
                color_dark = COLOR_MAPS_TERRAIN;
                color_light = COLOR_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
        end
        %## AXES LIMITS
        fig = figure('color','white','renderer','Painters');
        TITLE_XSHIFT = 0.4;
        TITLE_YSHIFT = 0.975;
        TITLE_BOX_SZ = [0.4,0.4];
        set(fig,'Units','inches','Position',FIGURE_POSITION)
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        set(gca,AXES_DEFAULT_PROPS{:});
        %##
        horiz_shift = 0;
        vert_shift = 0;
        hz = 1;
        %- ylim autoset
        temp_table = FOOOF_TABLE(FOOOF_TABLE.design_id == num2str(des_i),:);
        YLIMS = [prctile(temp_table.(measure_plot),1)-std(temp_table.(measure_plot))*1.5,...
            prctile(temp_table.(measure_plot),99)+std(temp_table.(measure_plot))*3];
        for j = 1:length(clusters)
            %## STATS\
            cl_i = double(string(clusters(j)));
            temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{[]},...
                          'regress_xvals',0,...
                          'subject_char',[],... % this option when filled prints removal of nan() info
                          'group_order',categorical({''}),...
                          'do_include_intercept',false);
            STATS_STRUCT = DEFAULT_STATS_STRUCT;
            cnt = 1;
            switch des_i
                case 1
                    tmptmp_table = table(categorical(string(temp_table.subj_char)),double(temp_table.(measure_plot)),...
                        categorical(string(temp_table.cond_id)),categorical(string(temp_table.cond_char)),...
                        categorical(string(temp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_id','cond_char','group_id'});
                case 2
                    tmptmp_table = table(categorical(string(temp_table.subj_char)),double(temp_table.(measure_plot)),...
                        double(string(temp_table.cond_char)),double(string(temp_table.cond_char)),...
                        categorical(string(temp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_id','cond_char','group_id'});
            end
            % for g_i = 1:length(unique(temp_table.group_id))
            %     switch des_i
            %         case 1
            %             %## LINEAR MODEL
            %             % t1.log_avg_power= log(t1.(MEASURE_NAMES{meas_i})+5);
            %             mod_lme = sprintf('%s ~ 1 + cond_id + (1|subj_char)',measure_plot);
            %             % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
            %             stats_out = fitlme(tmptmp_table,mod_lme);
            %             anova_out = anova(stats_out);
            %             %## GATHER STATS
            %             %- test normality
            %             [norm_h,norm_p] = lillietest(stats_out.residuals);
            %             %- get effects
            %             [~,bnames,~] = stats_out.fixedEffects();
            %             [~,brnames,bretable] = stats_out.randomEffects();
            %             %- intercept only model
            %             % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            %             % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            %             % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            %             R2 = stats_out.Rsquared.Adjusted;
            %             %- intercept only model
            %             altmod_out = sprintf('%s ~ 1 + (1|subj_char)',measure_plot);
            %             altstats_out = fitlme(tmptmp_table,altmod_out);
            %             %- alternative f2?
            %             R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
            %             R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
            %             alt_f2 = (R22-R21)/(1-R22);
            %             %- populate struct
            %             stats_struct(cnts).stat_test_mod = mod_lme;
            %             stats_struct(cnts).measure_tag = categorical({measure_plot});
            %             stats_struct(cnts).design_tag = categorical(des_i);
            %             stats_struct(cnts).mod_tag = categorical(2);
            %             % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
            %             stats_struct(cnts).mod_resp_terms = measure_plot;
            %             stats_struct(cnts).anova_preds_terms = t(:,1)';
            %             tmp = t(:,7)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_p = tmp;
            %             tmp = t(:,6)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_stat = tmp;
            %             tmp = t(:,3)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_df =tmp;
            %             stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
            %             stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
            %             stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
            %             stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
            %             stats_struct(cnts).multi_comp_terms = gnames';
            %             stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
            %             [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
            %             stats_struct(cnts).multi_comp_p = adj_p;
            %             stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
            %             stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
            %             stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
            %             stats_struct(cnts).norm_test_p = norm_p;
            %             stats_struct(cnts).norm_test_h = norm_h;
            %             stats_struct(cnts).effect_size = alt_f2;
            %             stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
            %             cnts = cnts + 1;
            %             stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
            %             %## PLOT =============================================================== %%
            %             %##
            %             aa = anova_out.pValue(strcmp(anova_out.Term,'cond_id'));
            %             c2s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_2'));
            %             c3s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_3'));
            %             c4s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_4'));
            %             STATS_STRUCT(cnt).anova{group_i}=aa;
            %             STATS_STRUCT(cnt).pvals{group_i}=[1,c2s,c3s,c4s];
            %             STATS_STRUCT(cnt).pvals_pairs{group_i}={[1,1],[1,2],[1,3],[1,4]};
            %         case 2
            %             %## LINEAR MODEL
            %             % t1.log_avg_power= log(t1.(measure_plot)+5);
            %             % mod_lme = sprintf('%s ~ 1 + cond_id + (1|subj_char)',measure_plot);
            %             mod_lme = sprintf('%s ~ 1 + cond_id + group_id + (1|subj_char)',measure_plot);
            %             % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
            %             stats_out = fitlme(tmptmp_table,mod_lme);
            %             anova_out = anova(stats_out);
            %             %## GATHER STATS
            %             %- test normality
            %             [norm_h,norm_p] = lillietest(stats_out.residuals);
            %             %- get effects
            %             [~,bnames,~] = stats_out.fixedEffects();
            %             [~,brnames,bretable] = stats_out.randomEffects();
            %             %- intercept only model
            %             % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
            %             % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
            %             % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
            %             R2 = stats_out.Rsquared.Adjusted;
            %             %- intercept only model
            %             altmod_out = sprintf('%s ~ 1 + (1|subj_char)',measure_plot);
            %             altstats_out = fitlme(tmptmp_table,altmod_out);
            %             %- alternative f2?
            %             R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
            %             R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
            %             alt_f2 = (R22-R21)/(1-R22);
            %             %- populate struct
            %             stats_struct(cnts).stat_test_mod = mod_lme;
            %             stats_struct(cnts).measure_tag = categorical({measure_plot});
            %             stats_struct(cnts).design_tag = categorical(des_i);
            %             stats_struct(cnts).mod_tag = categorical(2);
            %             % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
            %             stats_struct(cnts).mod_resp_terms = measure_plot;
            %             stats_struct(cnts).anova_preds_terms = t(:,1)';
            %             tmp = t(:,7)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_p = tmp;
            %             tmp = t(:,6)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_stat = tmp;
            %             tmp = t(:,3)';
            %             tmp = tmp(~cellfun(@isempty,tmp));
            %             stats_struct(cnts).anova_preds_df =tmp;
            %             stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
            %             stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
            %             stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
            %             stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
            %             stats_struct(cnts).multi_comp_terms = gnames';
            %             stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
            %             [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
            %             stats_struct(cnts).multi_comp_p = adj_p;
            %             stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
            %             stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
            %             stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
            %             stats_struct(cnts).norm_test_p = norm_p;
            %             stats_struct(cnts).norm_test_h = norm_h;
            %             stats_struct(cnts).effect_size = alt_f2;
            %             stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
            %             cnts = cnts + 1;
            %             stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
            %             %## PLOT =============================================================== %%
            %             %##
            %             aa =  anova_out.pValue(strcmp(anova_out.Term,'cond_id'));
            %             c2s = [];
            %             c3s = [];
            %             c4s = [];
            %             rs = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id')));
            %             rls = [double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'))),...
            %                 double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_id')))]; 
            %             r2 = R2;
            %             STATS_STRUCT(cnt).anova{group_i}=anova_out.pValue(strcmp(anova_out.Term,'cond_id'));
            %             STATS_STRUCT(cnt).regress_pval{group_i}=rs;
            %             STATS_STRUCT(cnt).regress_line{group_i}=rls;
            %             STATS_STRUCT(cnt).r2_coeff(group_i)=r2;
            %             STATS_STRUCT(cnt).regress_xvals=(0:5)*0.25;
            %     end
            %     cnt = cnt + 1;
            % end
            %## GROUPWISE
            switch des_i
                case 1
                    % mod_lme = sprintf('%s~1+group_id+cond_id',measure_plot);
                    mod_out = sprintf('%s ~ 1 + cond_id + group_id + (1|subj_char)',measure_plot);
                    stats_out = fitlme(tmptmp_table,mod_out);
                    pred_terms = stats_out.CoefficientNames;
                    % anova_out = anova(stats_out);
                    [p,t,anova_out,terms] = anovan(tmptmp_table.(measure_plot),{tmptmp_table.cond_id, tmptmp_table.group_id},...
                        'sstype',3,'varnames',{'trial_char','group_id'},'model','linear','Display','off');
                    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[2],...
                        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
                    disp(stats_out)
                    %- test normality
                    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
                    %- intercept only model
                    altmod_out = sprintf('%s ~ 1',measure_plot);
                    altstats_out = fitlm(tmptmp_table,altmod_out);
                    %- alternative f2?
                    R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                    R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                    alt_f2 = (R22-R21)/(1-R22);
                    %##
                    %- populate struct
                    stats_struct(cnts).stat_test_mod = mod_out;
                    stats_struct(cnts).measure_tag = categorical({measure_plot});
                    stats_struct(cnts).design_tag = categorical(des_i);
                    stats_struct(cnts).mod_tag = categorical(2);
                    % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
                    stats_struct(cnts).mod_resp_terms = measure_plot;
                    stats_struct(cnts).anova_preds_terms = t(:,1)';
                    tmp = t(:,7)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_p = tmp;
                    tmp = t(:,6)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_stat = tmp;
                    tmp = t(:,3)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_df =tmp;
                    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
                    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
                    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
                    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
                    stats_struct(cnts).multi_comp_terms = gnames';
                    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
                    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                    stats_struct(cnts).multi_comp_p = adj_p;
                    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
                    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
                    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
                    stats_struct(cnts).norm_test_p = norm_p;
                    stats_struct(cnts).norm_test_h = norm_h;
                    stats_struct(cnts).effect_size = alt_f2;
                    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                    cnts = cnts + 1;
                    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                    %## PLOT =============================================================== %%
                    multi_p = adj_p;
                    aa_p = t(:,7)';
                    mod_coeff = stats_out.Coefficients.Estimate;
                    anova_p = aa_p{2};
                    anova_grp_p = aa_p{3};
                    terr_p = [1,stats_out.Coefficients.pValue(2:4)']';
                    grp_p = [0,stats_out.Coefficients.pValue(3:4)']';
                    speed_r2 = stats_out.Rsquared.Adjusted;
                    speed_xvals = (0:5)*0.25;
                    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
                        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
                        'pvals',{{(terr_p),(terr_p),(terr_p)}},...
                        'pvals_pairs',{{{[1,1],[1,2],[1,3],[1,4]},...
                            {[1,1],[1,2],[1,3],[1,4]},...
                            {[1,1],[1,2],[1,3],[1,4]}}},...
                        'pvals_grp',{num2cell(adj_p)},...
                        'pvals_grp_pairs',{num2cell(comparisons(:,1:2),2)},...
                        'regress_pval',{{}},...
                        'regress_line',{{}},...
                        'line_type',{'best_fit'},... % ('best_fit' | 'means')
                        'regress_xvals',speed_xvals,...
                        'subject_char',[],... % this option when filled prints removal of nan() info
                        'group_order',categorical({''}),...
                        'display_stats_char',true,...
                        'stats_char',{{}},...
                        'bracket_conn_yshift',[1,1,1],...
                        'bracket_rawshifty_upper',0.2,...
                        'bracket_rawshifty_lower',0,...
                        'grp_sig_offset_x',[0,0,0],... %zeros(length(unique(tmp_table.(GROUP_TABLE_VAR)))),...
                        'grp_sig_offset_y',[0,0,0]);
                case 2
                    % mod_lme = sprintf('%s~1+group_id+cond_id',measure_plot);
                    mod_out = sprintf('%s ~ 1 + cond_id + group_id + (1|subj_char)',measure_plot);
                    stats_out = fitlme(tmptmp_table,mod_out);
                    pred_terms = stats_out.CoefficientNames;
                    % anova_out = anova(stats_out);
                    [p,t,anova_out,terms] = anovan(tmptmp_table.(measure_plot),{tmptmp_table.cond_id, tmptmp_table.group_id},...
                        'sstype',3,'varnames',{'trial_char','group_id'},'model','linear','Display','off');
                    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[2],...
                        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
                    disp(stats_out)
                    %- test normality
                    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
                    %- intercept only model
                    altmod_out = sprintf('%s ~ 1',measure_plot);
                    altstats_out = fitlm(tmptmp_table,altmod_out);
                    %- alternative f2?
                    R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                    R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                    alt_f2 = (R22-R21)/(1-R22);
                    %##
                    %- populate struct
                    stats_struct(cnts).stat_test_mod = mod_out;
                    stats_struct(cnts).measure_tag = categorical({measure_plot});
                    stats_struct(cnts).design_tag = categorical(des_i);
                    stats_struct(cnts).mod_tag = categorical(2);
                    % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
                    stats_struct(cnts).mod_resp_terms = measure_plot;
                    stats_struct(cnts).anova_preds_terms = t(:,1)';
                    tmp = t(:,7)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_p = tmp;
                    tmp = t(:,6)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_stat = tmp;
                    tmp = t(:,3)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_df =tmp;
                    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
                    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
                    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
                    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
                    stats_struct(cnts).multi_comp_terms = gnames';
                    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
                    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                    stats_struct(cnts).multi_comp_p = adj_p;
                    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
                    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
                    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
                    stats_struct(cnts).norm_test_p = norm_p;
                    stats_struct(cnts).norm_test_h = norm_h;
                    stats_struct(cnts).effect_size = alt_f2;
                    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                    cnts = cnts + 1;
                    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                    %## PLOT =============================================================== %%
                    multi_p = adj_p;
                    aa_p = t(:,7)';
                    mod_coeff = stats_out.Coefficients.Estimate;
                    anova_p = aa_p{2};
                    anova_grp_p = aa_p{3};
                    speed_p = stats_out.Coefficients.pValue(2);
                    grp_p = [0,stats_out.Coefficients.pValue(3:4)']';
                    speed_r2 = stats_out.Rsquared.Adjusted;
                    speed_xvals = (0:5)*0.25;
                    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
                        'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
                        'pvals',{{}},...
                        'pvals_pairs',{{}},...
                        'pvals_grp',{num2cell(adj_p)},...
                        'pvals_grp_pairs',{num2cell(comparisons(:,1:2),2)},...
                        'regress_pval',{{speed_p,speed_p,speed_p}},...
                        'regress_line',{{[mod_coeff(1), mod_coeff(2)],...
                            [mod_coeff(1)+mod_coeff(3), mod_coeff(2)],...
                            [mod_coeff(1)+mod_coeff(4), mod_coeff(2)]}},...
                        'line_type',{'best_fit'},... % ('best_fit' | 'means')
                        'regress_xvals',speed_xvals,...
                        'subject_char',[],... % this option when filled prints removal of nan() info
                        'group_order',categorical({''}),...
                        'display_stats_char',true,...
                        'stats_char',{{}},...
                        'bracket_conn_yshift',[1,1,1],...
                        'bracket_rawshifty_upper',0.2,...
                        'bracket_rawshifty_lower',0,...
                        'grp_sig_offset_x',[0,0,0],... %zeros(length(unique(tmp_table.(GROUP_TABLE_VAR)))),...
                        'grp_sig_offset_y',[0,0,0]);
            end
            %##
            % cl_i = double(string(clusters(j)));
            % temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            %- ylim autoset
            % % YLIMS = [min(temp_table.(measure_plot))-0.5,max(temp_table.(measure_plot))+1];
            %## PLOT
            ax = axes();
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_STRUCT = struct('color_map',color_dark,...
                'cond_labels',{cellstr(string(unique(temp_table.cond_char)))},...
                'cond_offsets',[-0.35,-0.1,0.15,0.40],...
                'group_labels',{g_chars_subp},...
                'group_offsets',[0.125,0.475,0.812],...
                'group_lab_yoffset',-0.23,...
                'group_lab_fontweight','normal',...
                'group_lab_fontsize',10,...
                'y_label',measure_plot,...
                'y_label_fontsize',10,...
                'y_label_fontweight','bold',...
                'ylim',YLIMS,...
                'x_label',x_label,...
                'x_label_fontsize',10,...
                'x_label_fontweight','bold',...
                'x_label_yoffset',-.14,...
                'xlim',[],...
                'title',{{sprintf('Cluster %i',cl_i)}},...
                'title_fontsize',10,...
                'title_fontweight','bold',...
                'font_size',9,...
                'font_name','Arial',...
                'do_combine_groups',false,...
                'regresslab_txt_size',6,...
                'ax_position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE],...
                'ax_line_width',1);
            ax = group_violin(tmptmp_table,measure_plot,'cond_id','group_id',...
                ax,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_STRUCT,...
                'STATS_STRUCT',STATS_STRUCT);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            
            %##
            if hz < HZ_DIM
                horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            else
                vert_shift = vert_shift - (AX_H*IM_RESIZE + AX_VERT_SHIFT*IM_RESIZE);
                horiz_shift = 0;
                hz = 0;
            end
            hz = hz + 1;
        end
        drawnow;
        %## SAVE
        hold off;
        savefig(fig,[save_dir filesep sprintf('%s_violin_group_plot_des%i.fig',measure_plot,des_i)])
        exportgraphics(fig,[save_dir filesep sprintf('%s_violin_group_plot_des%i.png',measure_plot,des_i)],'Resolution',300)
        close(fig);
    end
end
%}
%% (STATISTICS CALCULATIONS) =========================================== %%
%## STRUCT IMPLEMENTATION
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{''},...
    'measure_tag',categorical({''}),...
    'design_tag',categorical({''}),...
    'group_tag',categorical({''}),...
    'cluster_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{{''}},...
    'anova_preds_p',{[]},...
    'anova_preds_stat',{[]},...
    'anova_preds_df1',{[]},...
    'anova_preds_df2',{[]},...
    'mod_preds_terms',{{''}},...
    'mod_preds_p',[],...
    'mod_preds_stat',[],...
    'mod_preds_coeff',[],...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',[],...
    'multi_comp_p',[],...
    'multi_comp_coeff',[],...
    'multi_comp_lci_uci',[],...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_struct = DEF_STATS_TRACK_STRUCT;
DEF_STATS_TRACK_XLSX = struct('stat_test_mod',{''},...
    'measure_tag',categorical({''}),...
    'design_tag',categorical({''}),...
    'group_tag',categorical({''}),...
    'cluster_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{{''}},...
    'anova_preds_p',{''},...
    'anova_preds_stat',{''},...
    'anova_preds_df1',{''},...
    'anova_preds_df2',{''},...
    'mod_preds_terms',{''},...
    'mod_preds_p',{''},...
    'mod_preds_stat',{''},...
    'mod_preds_coeff',{''},...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',{''},...
    'multi_comp_p',{''},...
    'multi_comp_coeff',{''},...
    'multi_comp_lci_uci',{''},...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_xlsx = DEF_STATS_TRACK_XLSX;
cnts = 1;
%-
TERRAIN_DES_ID = 1;
SPEED_DES_ID = 2;
meas_names = {'theta_avg_power','alpha_avg_power','beta_avg_power'};
% meas_names = {'theta_avg_power','alpha_avg_power','beta_avg_power','alpha1_avg_power','alpha2_avg_power','beta1_avg_power','beta2_avg_power'};
% LOG_MEASURE_NAMES = {'log_theta_avg_power','log_alpha_avg_power','log_beta_avg_power'};
%% ===================================================================== %%
SPEED_MEAS_ANL = 'cond_char';
% SPEED_MEAS_ANL = 'speed_div_stat';
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
cntsts = 1;
%-
for des_i = 1:length(designs)
    for cl_i = 1:length(clusters)
        for g_i = 1:length(groups)
            %##
            %- extract table for cluster and design
            inds = FOOOF_TABLE.cluster_id == clusters(cl_i) &...
                FOOOF_TABLE.design_id == designs(des_i) & FOOOF_TABLE.group_id == groups(g_i);
            tmptmp = FOOOF_TABLE(inds,:);
            switch double(string(des_i))
                case TERRAIN_DES_ID
                    for meas_i = 1:length(meas_names)
                        %- NOTE: need to reassign to new table because of how
                        %categorical variables will hold onto removed
                        %entries causing rank defiecencies.
                        tmp = table(categorical(string(tmptmp.subj_char)),double(tmptmp.(meas_names{meas_i})),...
                            categorical(string(tmptmp.cond_id)),'VariableNames',{'subj_char',meas_names{meas_i},'cond_id'});
                        %## LINEAR MODEL
                        % t1.log_avg_power= log(t1.(MEASURE_NAMES{meas_i})+5);
                        mod_out = sprintf('%s ~ 1 + cond_id + (1|subj_char)',meas_names{meas_i});
                        % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
                        % stats_out = fitlme(tmp,mod_out,'FitMethod','ML','DummyVarCoding','effects');
                        stats_out = fitlme(tmp,mod_out,'FitMethod','REML','DummyVarCoding','effects');
                        anova_out = anova(stats_out);
                        %## GATHER STATS
                        %- test normality
                        [norm_h,norm_p] = lillietest(stats_out.residuals);
                        %- get effects
                        [~,bnames,~] = stats_out.fixedEffects();
                        [~,brnames,bretable] = stats_out.randomEffects();
                        %- intercept only model
                        % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
                        % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
                        % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
                        R2 = stats_out.Rsquared.Adjusted;
                        %- intercept only model
                        altmod_out = sprintf('%s ~ 1 + (1|subj_char)',meas_names{meas_i});
                        % altstats_out = fitlme(tmp,altmod_out);
                        altstats_out = fitlme(tmp,altmod_out,'FitMethod','REML','DummyVarCoding','effects');
                        %- alternative f2?
                        R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                        R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                        alt_f2 = (R22-R21)/(1-R22);
                        % disp(stats_out)
                        %- populate xlsx
                        stats_xlsx(cnts).stat_test_mod = mod_out;
                        stats_xlsx(cnts).measure_tag = categorical(meas_names(meas_i));
                        stats_xlsx(cnts).design_tag = designs(des_i);
                        stats_xlsx(cnts).group_tag = groups(g_i);
                        stats_xlsx(cnts).cluster_tag = clusters(cl_i);
                        % stats_xlsx(cnts).resp_terms = meas_names(meas_i);
                        stats_xlsx(cnts).mod_resp_terms = meas_names{meas_i};
                        %-
                        stats_xlsx(cnts).anova_preds_terms = cell2csv_util(anova_out.Term); %anova_out.Term';
                        stats_xlsx(cnts).anova_preds_p = cell2csv_util(anova_out.pValue); %anova_out.pValue';
                        stats_xlsx(cnts).anova_preds_stat = cell2csv_util(anova_out.FStat); %anova_out.FStat';
                        stats_xlsx(cnts).anova_preds_df1 = cell2csv_util(anova_out.DF1); %anova_out.DF1';
                        stats_xlsx(cnts).anova_preds_df2 = cell2csv_util(anova_out.DF2); %anova_out.DF2';
                        stats_xlsx(cnts).mod_preds_p = cell2csv_util(stats_out.Coefficients.pValue);
                        stats_xlsx(cnts).mod_preds_terms = cell2csv_util(stats_out.Coefficients.Name);
                        stats_xlsx(cnts).mod_preds_stat = cell2csv_util(stats_out.Coefficients.tStat);
                        stats_xlsx(cnts).mod_preds_coeff = cell2csv_util(stats_out.Coefficients.Estimate);
                        % stats_xlsx(cnts).multi_comp_terms = cell2csv_util(gnames);
                        % stats_xlsx(cnts).multi_comp_t1_t2 = cell2csv_util(comparisons(:,1:2));
                        % [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                        % stats_xlsx(cnts).multi_comp_p = cell2csv_util(adj_p);
                        % stats_xlsx(cnts).multi_comp_coeff = cell2csv_util(comparisons(:,4));
                        % stats_xlsx(cnts).multi_comp_lci_uci = cell2csv_util([comparisons(:,3),comparisons(:,5)]);
                        stats_xlsx(cnts).mod_r2 = cell2csv_util(stats_out.Rsquared.Adjusted);
                        stats_xlsx(cnts).norm_test_p = norm_p;
                        stats_xlsx(cnts).norm_test_h = norm_h;
                        stats_xlsx(cnts).effect_size = alt_f2;
                        stats_xlsx(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        %- populate struct (2)
                        stats_struct(cnts).stat_test_mod = mod_out;
                        stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
                        stats_struct(cnts).design_tag = designs(des_i);
                        stats_struct(cnts).group_tag = groups(g_i);
                        stats_struct(cnts).cluster_tag = clusters(cl_i);
                        % stats_struct(cnts).resp_terms = meas_names(meas_i);
                        stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
                        %-
                        stats_struct(cnts).anova_preds_terms = (anova_out.Term'); %anova_out.Term';
                        stats_struct(cnts).anova_preds_p = (anova_out.pValue'); %anova_out.pValue';
                        stats_struct(cnts).anova_preds_stat = (anova_out.FStat'); %anova_out.FStat';
                        stats_struct(cnts).anova_preds_df1 = (anova_out.DF1'); %anova_out.DF1';
                        stats_struct(cnts).anova_preds_df2 = (anova_out.DF2'); %anova_out.DF2';
                        stats_struct(cnts).mod_preds_p = (stats_out.Coefficients.pValue');
                        stats_struct(cnts).mod_preds_terms = (stats_out.Coefficients.Name');
                        stats_struct(cnts).mod_preds_stat = (stats_out.Coefficients.tStat');
                        stats_struct(cnts).mod_preds_coeff = (stats_out.Coefficients.Estimate');
                        % stats_struct(cnts).multi_comp_terms = (gnames);
                        % stats_struct(cnts).multi_comp_t1_t2 = (comparisons(:,1:2));
                        % [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                        % stats_struct(cnts).multi_comp_p = (adj_p);
                        % stats_struct(cnts).multi_comp_coeff = (comparisons(:,4));
                        % stats_struct(cnts).multi_comp_lci_uci = ([comparisons(:,3),comparisons(:,5)]);
                        stats_struct(cnts).mod_r2 = (stats_out.Rsquared.Adjusted);
                        stats_struct(cnts).norm_test_p = norm_p;
                        stats_struct(cnts).norm_test_h = norm_h;
                        stats_struct(cnts).effect_size = alt_f2;
                        stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        cnts = cnts + 1;
                        stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                        stats_xlsx(cnts) = DEF_STATS_TRACK_XLSX;
                    end
                case SPEED_DES_ID
                    for meas_i = 1:length(meas_names)
                        %- NOTE: need to reassign to new table because of how
                        %categorical variables will hold onto removed
                        %entries causing rank defiecencies.
                        %- speed factors
                        % tmp = table(categorical(string(tmptmp.subj_char)),double(tmptmp.(MEASURE_NAMES{meas_i})),...
                        %     categorical(string(tmptmp.cond_id)),'VariableNames',{'subj_char',MEASURE_NAMES{meas_i},'cond_id'});
                        %- speed cont.
                        % tmp = table(categorical(string(tmptmp.subj_char)),double(tmptmp.(MEASURE_NAMES{meas_i})),...
                        %     double(string(tmptmp.cond_id)),'VariableNames',{'subj_char',MEASURE_NAMES{meas_i},'cond_id'});
                        %- speed_diff_stat.
                        tmp = table(categorical(string(tmptmp.subj_char)),double(tmptmp.(meas_names{meas_i})),...
                            double(string(tmptmp.(SPEED_MEAS_ANL))),'VariableNames',{'subj_char',meas_names{meas_i},'cond_id'});
                        %## LINEAR MODEL
                        % t1.log_avg_power= log(t1.(MEASURE_NAMES{meas_i})+5);
                        mod_out = sprintf('%s ~ 1 + cond_id + (1|subj_char)',meas_names{meas_i});
                        % mod_lme = sprintf('%s ~ 1 + cond_id^2 + cond_id + (1|subj_char)',MEASURE_NAMES{meas_i});
                        % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
                        % stats_out = fitlme(tmp,mod_out,'FitMethod','ML','DummyVarCoding','effects');
                        stats_out = fitlme(tmp,mod_out,'FitMethod','REML','DummyVarCoding','effects');
                        anova_out = anova(stats_out);
                        % [p,t,anova_out,terms] = anovan(tmp.(meas_names{meas_i}),{tmp.cond_id},...
                        %     'sstype',3,'varnames',{'cond_id'},'model','linear','Display','off');
                        % [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[1,2],...
                        %     'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
                        % disp(stats_out)
                        %## GATHER STATS
                        %- test normality
                        [norm_h,norm_p] = lillietest(stats_out.residuals);
                        %- get effects
                        [~,bnames,~] = stats_out.fixedEffects();
                        [~,brnames,bretable] = stats_out.randomEffects();
                        %- intercept only model
                        % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
                        % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
                        % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
                        R2 = stats_out.Rsquared.Adjusted;
                        %- intercept only model
                        altmod_out = sprintf('%s ~ 1 + (1|subj_char)',meas_names{meas_i});
                        % altstats_out = fitlme(tmp,altmod_out);
                        altstats_out = fitlme(tmp,altmod_out,'FitMethod','REML','DummyVarCoding','effects');
                        %- alternative f2?
                        R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                        R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                        alt_f2 = (R22-R21)/(1-R22);
                        %## LINEAR
                        %- populate xlsx
                        stats_xlsx(cnts).stat_test_mod = mod_out;
                        stats_xlsx(cnts).measure_tag = categorical(meas_names(meas_i));
                        % stats_xlsx(cnts).design_tag = categorical(des_i);
                        stats_xlsx(cnts).design_tag = designs(des_i);
                        stats_xlsx(cnts).group_tag = groups(g_i);
                        stats_xlsx(cnts).cluster_tag = clusters(cl_i);
                        stats_xlsx(cnts).mod_resp_terms = meas_names{meas_i};
                        stats_xlsx(cnts).anova_preds_terms = cell2csv_util(anova_out.Term); %anova_out.Term';
                        stats_xlsx(cnts).anova_preds_p = cell2csv_util(anova_out.pValue); %anova_out.pValue';
                        stats_xlsx(cnts).anova_preds_stat = cell2csv_util(anova_out.FStat); %anova_out.FStat';
                        stats_xlsx(cnts).anova_preds_df1 = cell2csv_util(anova_out.DF1); %anova_out.DF1';
                        stats_xlsx(cnts).anova_preds_df2 = cell2csv_util(anova_out.DF2); %anova_out.DF2';
                        stats_xlsx(cnts).mod_preds_p = cell2csv_util(stats_out.Coefficients.pValue);
                        stats_xlsx(cnts).mod_preds_terms = cell2csv_util(stats_out.Coefficients.Name);
                        stats_xlsx(cnts).mod_preds_stat = cell2csv_util(stats_out.Coefficients.tStat);
                        stats_xlsx(cnts).mod_preds_coeff = cell2csv_util(stats_out.Coefficients.Estimate);
                        % stats_xlsx(cnts).multi_comp_terms = cell2csv_util(gnames);
                        % stats_xlsx(cnts).multi_comp_t1_t2 = cell2csv_util(comparisons(:,1:2));
                        % [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                        % stats_xlsx(cnts).multi_comp_p = cell2csv_util(adj_p);
                        % stats_xlsx(cnts).multi_comp_coeff = cell2csv_util(comparisons(:,4));
                        % stats_xlsx(cnts).multi_comp_lci_uci = cell2csv_util([comparisons(:,3),comparisons(:,5)]);
                        stats_xlsx(cnts).mod_r2 = cell2csv_util(stats_out.Rsquared.Adjusted);
                        stats_xlsx(cnts).norm_test_p = norm_p;
                        stats_xlsx(cnts).norm_test_h = norm_h;
                        stats_xlsx(cnts).effect_size = alt_f2;
                        stats_xlsx(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        %- populate struct (2)
                        stats_struct(cnts).stat_test_mod = mod_out;
                        stats_struct(cnts).measure_tag = categorical(meas_names(meas_i));
                        stats_struct(cnts).design_tag = designs(des_i);
                        stats_struct(cnts).group_tag = groups(g_i);
                        stats_struct(cnts).cluster_tag = clusters(cl_i);
                        stats_struct(cnts).mod_resp_terms = meas_names{meas_i};
                        stats_struct(cnts).anova_preds_terms = (anova_out.Term'); %anova_out.Term';
                        stats_struct(cnts).anova_preds_p = (anova_out.pValue'); %anova_out.pValue';
                        stats_struct(cnts).anova_preds_stat = (anova_out.FStat'); %anova_out.FStat';
                        stats_struct(cnts).anova_preds_df1 = (anova_out.DF1'); %anova_out.DF1';
                        stats_struct(cnts).anova_preds_df2 = (anova_out.DF2'); %anova_out.DF2';
                        stats_struct(cnts).mod_preds_p = (stats_out.Coefficients.pValue');
                        stats_struct(cnts).mod_preds_terms = (stats_out.Coefficients.Name');
                        stats_struct(cnts).mod_preds_stat = (stats_out.Coefficients.tStat');
                        stats_struct(cnts).mod_preds_coeff = (stats_out.Coefficients.Estimate');
                        % stats_struct(cnts).multi_comp_terms = (gnames);
                        % stats_struct(cnts).multi_comp_t1_t2 = (comparisons(:,1:2));
                        % [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                        % stats_struct(cnts).multi_comp_p = (adj_p);
                        % stats_struct(cnts).multi_comp_coeff = (comparisons(:,4));
                        % stats_struct(cnts).multi_comp_lci_uci = ([comparisons(:,3),comparisons(:,5)]);
                        stats_struct(cnts).mod_r2 = (stats_out.Rsquared.Adjusted);
                        stats_struct(cnts).norm_test_p = norm_p;
                        stats_struct(cnts).norm_test_h = norm_h;
                        stats_struct(cnts).effect_size = alt_f2;
                        stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        cnts = cnts + 1;
                        stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                        stats_xlsx(cnts) = DEF_STATS_TRACK_XLSX;
                    end
                otherwise
                    error('Design %i not defined...\n',des_i);
            end
            
        end
    end
end
stats_struct = struct2table(stats_struct);
writetable(stats_struct,[save_dir filesep 'STATS_TRACK_STRUCT_TABLE_speedlin_alt.xlsx']);
% writetable(stats_struct,['./' filesep 'STATS_TRACK_STRUCT_TABLE_speedlin_alt.xlsx']);
save([save_dir filesep 'STATS_TRACK_STRUCT_speedlin_alt.mat'],'stats_struct');
%-
stats_xlsx = struct2table(stats_xlsx);
writetable(stats_xlsx,[save_dir filesep 'STATS_TRACK_STRUCT_TABLE_speedlin_csv_alt.xlsx']);
% writetable(stats_xlsx,['./' filesep 'STATS_TRACK_STRUCT_TABLE_speedlin_csv_alt.xlsx']);
save([save_dir filesep 'STATS_TRACK_STRUCT_speedlin_csv_alt.mat'],'stats_xlsx');

%% CONDENSED ANOVA PLOTS FOR EEG POWER
%{
%## (STATS STRUCT) ====================================================== %%
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'measure_tag',categorical({''}),...
    'measure_char',{''},...
    'design_tag',categorical({''}),...
    'design_num',[],...
    'cluster_num',[],...
    'mod_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{''},...
    'anova_preds_p',{[]},...
    'anova_preds_stat',{[]},...
    'anova_preds_df',{[]},...
    'mod_preds_terms',{{''}},...
    'mod_preds_p',[],...
    'mod_preds_stat',[],...
    'mod_preds_coeff',[],...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',[],...
    'multi_comp_p',[],...
    'multi_comp_coeff',[],...
    'multi_comp_lci_uci',[],...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_struct = DEF_STATS_TRACK_STRUCT;
cnts = 1;
%##
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
g_chars_subp = {'YA','OHMA','OLMA'};
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
FIGURE_POSITION =[1,1,6.5,9];
%-
FONT_SIZE_VIO_REG = 8;
FONT_SIZE_VIO = 6;
IM_RESIZE = 0.85;
AX_H  = 0.18;
AX_W = 0.27;
AX_HORIZ_SHIFT = 0.06;
AX_VERT_SHIFT = 0.11;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT_VIO = 0.8;
HZ_DIM = 3;
YLIMS = [-0.5,3];
AXES_FONT_SIZE_VIO = 10;
GROUP_LAB_FONTSIZE = 10;
GROUP_LAB_FONTWEIGHT = 'bold ';
XLAB_FONTSIZE = 10;
YLAB_FONTSIZE = 10;
XTICK_FONTSIZE = 10;
XLAB_FONTWEIGHT = 'bold';
YLAB_FONTWEIGHT = 'bold';
TITLE_FONTSIZE = 10;
TITLE_FONTWEIGHT = 'bold';
XLABEL_OFFSET = -.05;
GROUP_LAB_YOFFSET = -0.275;
MEASURES_PLOT = {'theta_avg_power','alpha_avg_power','beta_avg_power','alpha1_avg_power','alpha2_avg_power','beta1_avg_power','beta2_avg_power'};
% measure_plot = 'aperiodic_exp';
for k = 1:length(MEASURES_PLOT)
    measure_plot = MEASURES_PLOT{k};
    for i = 2 %1:length(designs)
        des_i = double(string(designs(i)));
        switch des_i
            case 1
                color_dark = COLOR_MAPS_TERRAIN;
                color_light = COLOR_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
        end
        %## AXES LIMITS
        fig = figure('color','white','renderer','Painters');
        TITLE_XSHIFT = 0.4;
        TITLE_YSHIFT = 0.975;
        TITLE_BOX_SZ = [0.4,0.4];
        set(fig,'Units','inches','Position',FIGURE_POSITION)
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        set(gca,AXES_DEFAULT_PROPS{:});
        %##
        horiz_shift = 0;
        vert_shift = 0;
        hz = 1;
        %- ylim autoset
        temp_table = FOOOF_TABLE(FOOOF_TABLE.design_id == num2str(des_i),:);
        YLIMS = [prctile(temp_table.(measure_plot),1)-std(temp_table.(measure_plot))*1.5,prctile(temp_table.(measure_plot),99)+std(temp_table.(measure_plot))*3];
        for j = 1:length(clusters)
            %## STATS
            cl_i = double(string(clusters(j)));
            temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                'anova_grp',{{}},...
                'pvals',{{}},...
                'pvals_pairs',{{}},...
                'pvals_grp',{{}},...
                'pvals_grp_pairs',{{}},...
                'regress_pval',{{}},...
                'regress_line',{{}},...
                'line_type',{'best_fit'},... % ('best_fit' | 'means')
                'regress_xvals',speed_xvals,...
                'subject_char',[],... % this option when filled prints removal of nan() info
                'group_order',categorical({''}),...
                'display_stats_char',true,...
                'stats_char',{{}},...
                'bracket_conn_yshift',[1,1,1],...
                'bracket_rawshifty_upper',0.2,...
                'bracket_rawshifty_lower',0,...
                'grp_sig_offset_x',[0,0,0],... %zeros(length(unique(tmp_table.(GROUP_TABLE_VAR)))),...
                'grp_sig_offset_y',[0,0,0]);
            STATS_STRUCT = DEFAULT_STATS_STRUCT;
            % cnt = 1;
            for g_i = 1:length(unique(temp_table.group_id))
                tmptmp_table = FOOOF_TABLE(FOOOF_TABLE.group_id == string(g_i) & FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
                
                switch des_i
                    case 1
                        tmptmp_table = table(categorical(string(tmptmp_table.subj_char)),double(tmptmp_table.(measure_plot)),...
                            categorical(string(tmptmp_table.cond_id)),categorical(string(tmptmp_table.cond_char)),...
                            categorical(string(tmptmp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_id','cond_char','group_id'});
                    case 2
                        tmptmp_table = table(categorical(string(tmptmp_table.subj_char)),double(tmptmp_table.(measure_plot)),...
                            double(string(tmptmp_table.cond_char)),...
                            categorical(string(tmptmp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_char','group_id'});
                end
                switch des_i
                    case 1
                        %## LINEAR MODEL
                        % t1.log_avg_power= log(t1.(MEASURE_NAMES{meas_i})+5);
                        mod_out = sprintf('%s ~ 1 + cond_id + (1|subj_char)',measure_plot);
                        % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
                        stats_out = fitlme(tmptmp_table,mod_out,'FitMethod','REML');
                        anova_out = anova(stats_out);
                        %## GATHER STATS
                        %- test normality
                        [norm_h,norm_p] = lillietest(stats_out.residuals);
                        %- get effects
                        [~,bnames,~] = stats_out.fixedEffects();
                        [~,brnames,bretable] = stats_out.randomEffects();
                        %- intercept only model
                        % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
                        % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
                        % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
                        R2 = stats_out.Rsquared.Adjusted;
                        %- intercept only model
                        altmod_out = sprintf('%s ~ 1 + (1|subj_char)',measure_plot);
                        altstats_out = fitlme(tmptmp_table,altmod_out);
                        %- alternative f2?
                        R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                        R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                        alt_f2 = (R22-R21)/(1-R22);
                        %- populate struct
                        stats_struct(cnts).stat_test_mod = mod_out;
                        stats_struct(cnts).measure_tag = categorical({measure_plot});
                        stats_struct(cnts).design_tag = categorical(des_i);
                        stats_struct(cnts).mod_tag = categorical(2);
                        % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
                        stats_struct(cnts).mod_resp_terms = measure_plot;
                        stats_struct(cnts).anova_preds_terms = t(:,1)';
                        tmp = t(:,7)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_p = tmp;
                        tmp = t(:,6)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_stat = tmp;
                        tmp = t(:,3)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_df =tmp;
                        stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
                        stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
                        stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
                        stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
                        stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
                        stats_struct(cnts).norm_test_p = norm_p;
                        stats_struct(cnts).norm_test_h = norm_h;
                        stats_struct(cnts).effect_size = alt_f2;
                        stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        cnts = cnts + 1;
                        stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                        %## PLOT =============================================================== %%
                        %##
                        aa = anova_out.pValue(strcmp(anova_out.Term,'cond_id'));
                        c2s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_2'));
                        c3s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_3'));
                        c4s = stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_id_4'));
                        STATS_STRUCT.anova{g_i}=aa;
                        STATS_STRUCT.pvals{g_i}=[1,c2s,c3s,c4s];
                        STATS_STRUCT.pvals_pairs{g_i}={[1,1],[1,2],[1,3],[1,4]};
                        
                    case 2
                       
                        %## LINEAR MODEL
                        % t1.log_avg_power= log(t1.(measure_plot)+5);
                        % mod_lme = sprintf('%s ~ 1 + cond_id + (1|subj_char)',measure_plot);
                        mod_out = sprintf('%s ~ 1 + cond_char + group_id + (1|subj_char)',measure_plot);
                        % mod_lme = 'theta_avg_power ~ 1 + cond + (1|speed_ms)';
                        stats_out = fitlme(tmptmp_table,mod_out,'FitMethod','REML');
                        anova_out = anova(stats_out);
                        %## GATHER STATS
                        %- test normality
                        [norm_h,norm_p] = lillietest(stats_out.residuals);
                        %- get effects
                        [~,bnames,~] = stats_out.fixedEffects();
                        [~,brnames,bretable] = stats_out.randomEffects();
                        %- intercept only model
                        % altmod_lme = sprintf('%s ~ 1 + (1|subj_char)',measure_name);
                        % altstats_out = fitlme(T_vals_plot,altmod_lme,'DummyVarCoding','effects');
                        % R2 = 1-(altstats_out.LogLikelihood/stats_out.LogLikelihood)^(2/stats_out.NumVariables);
                        R2 = stats_out.Rsquared.Adjusted;
                        %- intercept only model
                        altmod_out = sprintf('%s ~ 1 + (1|subj_char)',measure_plot);
                        altstats_out = fitlme(tmptmp_table,altmod_out);
                        %- alternative f2?
                        R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                        R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                        alt_f2 = (R22-R21)/(1-R22);
                        %- populate struct
                        stats_struct(cnts).stat_test_mod = mod_out;
                        stats_struct(cnts).measure_tag = categorical({measure_plot});
                        stats_struct(cnts).design_tag = categorical(des_i);
                        stats_struct(cnts).mod_tag = categorical(2);
                        % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
                        stats_struct(cnts).mod_resp_terms = measure_plot;
                        stats_struct(cnts).anova_preds_terms = t(:,1)';
                        tmp = t(:,7)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_p = tmp;
                        tmp = t(:,6)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_stat = tmp;
                        tmp = t(:,3)';
                        tmp = tmp(~cellfun(@isempty,tmp));
                        stats_struct(cnts).anova_preds_df =tmp;
                        stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
                        stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
                        stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
                        stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
                        stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
                        stats_struct(cnts).norm_test_p = norm_p;
                        stats_struct(cnts).norm_test_h = norm_h;
                        stats_struct(cnts).effect_size = alt_f2;
                        stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                        cnts = cnts + 1;
                        stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                        %## PLOT =============================================================== %%
                        %##
                        aa =  anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
                        c2s = [];
                        c3s = [];
                        c4s = [];
                        rs = double(stats_out.Coefficients.pValue(strcmp(stats_out.Coefficients.Name,'cond_char')));
                        rls = [double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'(Intercept)'))),...
                            double(stats_out.Coefficients.Estimate(strcmp(stats_out.Coefficients.Name,'cond_char')))]; 
                        r2 = R2;
                        STATS_STRUCT.anova{g_i}=anova_out.pValue(strcmp(anova_out.Term,'cond_char'));
                        STATS_STRUCT.regress_pval{g_i}=rs;
                        STATS_STRUCT.regress_line{g_i}=rls;
                        STATS_STRUCT.regress_xvals=(0:5)*0.25;
                        STATS_STRUCT.display_stats_char = true;
                        if aa > 0.01 && aa < 0.05
                            str = sprintf('* m=%0.2g\nb=%0.2g\nR^2=%0.2g',rls(2),rls(1),r2);
                        elseif aa <= 0.01 && aa > 0.001
                            str = sprintf('** m=%0.2g\nb=%0.2g\nR^2=%0.2g',rls(2),rls(1),r2);
                        elseif aa <= 0.001
                            str = sprintf('** m=%0.2g\nb=%0.2g\nR^2=%0.2g',rls(2),rls(1),r2);
                        else
                            str = '';
                        end
                        STATS_STRUCT.stats_char{g_i}=str;
                end
            end
            %##
            % cl_i = double(string(clusters(j)));
            % temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            %- ylim autoset
            % % YLIMS = [min(temp_table.(measure_plot))-0.5,max(temp_table.(measure_plot))+1];
            %## PLOT
            ax = axes();
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_STRUCT = struct('color_map',color_dark,...
                'cond_labels',{cellstr(string(unique(temp_table.cond_char)))},...
                'cond_offsets',[-0.35,-0.1,0.15,0.40],...
                'group_labels',{g_chars_subp},...
                'group_offsets',[0.125,0.475,0.812],...
                'group_lab_yoffset',-0.23,...
                'group_lab_fontweight','normal',...
                'group_lab_fontsize',10,...
                'y_label',measure_plot,...
                'y_label_fontsize',10,...
                'y_label_fontweight','bold',...
                'ylim',YLIMS,...
                'x_label',x_label,...
                'x_label_fontsize',10,...
                'x_label_fontweight','bold',...
                'x_label_yoffset',-.14,...
                'xlim',[],...
                'title',{{sprintf('Cluster %i',cl_i)}},...
                'title_fontsize',10,...
                'title_fontweight','bold',...
                'font_size',9,...
                'font_name','Arial',...
                'do_combine_groups',false,...
                'regresslab_txt_size',6,...
                'ax_position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE],...
                'ax_line_width',1);
            ax = group_violin(temp_table,measure_plot,'cond_char','group_id',...
                ax,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_STRUCT,...
                'STATS_STRUCT',STATS_STRUCT);
            %##
            if hz < HZ_DIM
                horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            else
                vert_shift = vert_shift - (AX_H*IM_RESIZE + AX_VERT_SHIFT*IM_RESIZE);
                horiz_shift = 0;
                hz = 0;
            end
            hz = hz + 1;
        end
        drawnow;
        %## SAVE
        hold off;
        savefig(fig,[save_dir filesep sprintf('%s_violin_group_plot_nongroup_des%i.fig',measure_plot,des_i)])
        exportgraphics(fig,[save_dir filesep sprintf('%s_violin_group_plot_nongroup_des%i.png',measure_plot,des_i)],'Resolution',300)
        close(fig);
    end
end
stats_struct_table = struct2table(stats_struct);
%}
%% RANOVAS FOR EEG POWER
%{
%## (STATS STRUCT) ====================================================== %%
DEF_STATS_TRACK_STRUCT = struct('stat_test_mod',{{''}},...
    'measure_tag',categorical({''}),...
    'measure_char',{''},...
    'design_tag',categorical({''}),...
    'design_num',[],...
    'cluster_num',[],...
    'mod_tag',categorical({''}),...
    'mod_resp_terms',{''},...
    'rnd_terms',{''},...
    'anova_preds_terms',{''},...
    'anova_preds_p',{[]},...
    'anova_preds_stat',{[]},...
    'anova_preds_df',{[]},...
    'mod_preds_terms',{{''}},...
    'mod_preds_p',[],...
    'mod_preds_stat',[],...
    'mod_preds_coeff',[],...
    'mod_r2',[],...
    'multi_comp_terms',{''},...
    'multi_comp_t1_t2',[],...
    'multi_comp_p',[],...
    'multi_comp_coeff',[],...
    'multi_comp_lci_uci',[],...
    'norm_test_p',[],...
    'norm_test_h',[],...
    'effect_size',[],...
    'effect_size_calc',{''});
stats_struct = DEF_STATS_TRACK_STRUCT;
cnts = 1;
%##
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
g_chars_subp = {'YA','OHMA','OLMA'};
%-
AXES_DEFAULT_PROPS = {'box','off','xtick',[],'ytick',[],...
        'ztick',[],'xcolor',[1,1,1],'ycolor',[1,1,1]};
FIGURE_POSITION =[1,1,6.5,9];
%-
FONT_SIZE_VIO_REG = 10;
FONT_SIZE_VIO = 6;
IM_RESIZE = 0.85;
AX_H  = 0.18;
AX_W = 0.27;
AX_HORIZ_SHIFT = 0.06;
AX_VERT_SHIFT = 0.11;
AX_INIT_HORIZ = 0.07;
AX_INIT_VERT_VIO = 0.8;
HZ_DIM = 3;
YLIMS = [-0.5,3];
AXES_FONT_SIZE_VIO = 10;
GROUP_LAB_FONTSIZE = 10;
GROUP_LAB_FONTWEIGHT = 'bold ';
XLAB_FONTSIZE = 10;
YLAB_FONTSIZE = 10;
XTICK_FONTSIZE = 10;
XLAB_FONTWEIGHT = 'bold';
YLAB_FONTWEIGHT = 'bold';
TITLE_FONTSIZE = 10;
TITLE_FONTWEIGHT = 'bold';
XLABEL_OFFSET = -.05;
GROUP_LAB_YOFFSET = -0.275;
MEASURES_PLOT = {'theta_avg_power','alpha_avg_power','beta_avg_power','alpha1_avg_power','alpha2_avg_power','beta1_avg_power','beta2_avg_power'};
% measure_plot = 'aperiodic_exp';
%%
for k = 1:length(MEASURES_PLOT)
    measure_plot = MEASURES_PLOT{k};
    for i = 2 %1:length(designs)
        %%
        des_i = double(string(designs(i)));
        switch des_i
            case 1
                color_dark = COLORS_MAPS_TERRAIN;
                color_light = COLORS_MAPS_TERRAIN;
                xtick_label_g = {'flat','low','med','high'};
                x_label = 'terrain';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.50','0.75','1.0'};
                x_label = 'speed (m/s)';
                cond_offsets = [-0.35,-0.1,0.15,0.40];
        end
        %## AXES LIMITS
        fig = figure('color','white','renderer','Painters');
        TITLE_XSHIFT = 0.4;
        TITLE_YSHIFT = 0.975;
        TITLE_BOX_SZ = [0.4,0.4];
        set(fig,'Units','inches','Position',FIGURE_POSITION)
        set(fig,'PaperUnits','inches','PaperSize',[1 1],'PaperPosition',[0 0 1 1])
        set(gca,AXES_DEFAULT_PROPS{:});
        %##
        horiz_shift = 0;
        vert_shift = 0;
        hz = 1;
        %- ylim autoset
        temp_table = FOOOF_TABLE(FOOOF_TABLE.design_id == num2str(des_i),:);
        YLIMS = [prctile(temp_table.(measure_plot),1)-std(temp_table.(measure_plot))*1.5,prctile(temp_table.(measure_plot),99)+std(temp_table.(measure_plot))*3];
        for j = 1:length(clusters)
            %## STATS
            cl_i = double(string(clusters(j)));
            temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            DEFAULT_STATS_STRUCT = struct('anova',{{}},...
                          'anova_grp',{{}},...
                          'pvals',{{}},...
                          'pvals_pairs',{{}},...
                          'pvals_grp',{{}},...
                          'pvals_grp_pairs',{{}},...
                          'regress_pval',{{}},...
                          'regress_line',{{}},...
                          'r2_coeff',{[]},...
                          'regress_xvals',0,...
                          'subject_char',[],... % this option when filled prints removal of nan() info
                          'group_order',categorical({''}),...
                          'do_include_intercept',false);
            STATS_STRUCT = DEFAULT_STATS_STRUCT;
            cnt = 1;
            switch des_i
                case 1
                    tmptmp_table = table(categorical(string(temp_table.subj_char)),double(temp_table.(measure_plot)),...
                        categorical(string(temp_table.cond_id)),categorical(string(temp_table.cond_char)),...
                        categorical(string(temp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_id','cond_char','group_id'});
                case 2
                    tmptmp_table = table(categorical(string(temp_table.subj_char)),double(temp_table.(measure_plot)),...
                        double(string(temp_table.cond_char)),...
                        categorical(string(temp_table.group_id)),'VariableNames',{'subj_char',measure_plot,'cond_char','group_id'});
            end
            %## GROUPWISE
            switch des_i
                case 1
                    % mod_lme = sprintf('%s~1+group_id+cond_id',measure_plot);
                    mod_lme = sprintf('%s ~ 1 + cond_id + group_id + (1|subj_char)',measure_plot);
                    stats_out = fitlme(tmptmp_table,mod_lme);
                    pred_terms = stats_out.CoefficientNames;
                    % anova_out = anova(stats_out);
                    [p,t,anova_out,terms] = anovan(tmptmp_table.(measure_plot),{tmptmp_table.cond_id, tmptmp_table.group_id},...
                        'sstype',3,'varnames',{'trial_char','group_id'},'model','linear','Display','off');
                    [comparisons,means,~,gnames] = multcompare(anova_out,'Dimension',[2],...
                        'display','off','Alpha',0.05); % comparisons columns: [pred1,pred2,lowerCI,estimate,upperCI,p-val]
                    disp(stats_out)
                    % %## RANOVA
                    % tt_out = unstack(tmp_table,meas_names{meas_i},'trial_char');
                    % % mod_out = sprintf('x0_25,x0_5,x0_75,x1~group_id+1');
                    % mod_out = sprintf('x0_25,x0_5,x0_75,x1_0~group_id+1');
                    % tmptmp = table([1,2,3,4]','VariableNames',{'speed'});
                    % %-
                    % rm_out = fitrm(tt_out,mod_out,'WithinDesign',tmptmp);
                    % [ran_out,A,C,D] = ranova(rm_out);
                    % [tbl] = rm_out.multcompare('group_id');
                    % %##
                    %- test normality
                    [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
                    %- intercept only model
                    altmod_out = sprintf('%s ~ 1',measure_plot);
                    altstats_out = fitlm(tmptmp_table,altmod_out);
                    %- alternative f2?
                    R21 = altstats_out.SSR/altstats_out.SST; % coefficient of determination
                    R22 = stats_out.SSR/stats_out.SST; % coefficient of determination
                    alt_f2 = (R22-R21)/(1-R22);
                    %##
                    %- populate struct
                    stats_struct(cnts).stat_test_mod = mod_lme;
                    stats_struct(cnts).measure_tag = categorical({measure_plot});
                    stats_struct(cnts).design_tag = categorical(des_i);
                    stats_struct(cnts).mod_tag = categorical(2);
                    % stats_struct(cnts).resp_terms = MEASURES_PLOT(meas_i);
                    stats_struct(cnts).mod_resp_terms = measure_plot;
                    stats_struct(cnts).anova_preds_terms = t(:,1)';
                    tmp = t(:,7)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_p = tmp;
                    tmp = t(:,6)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_stat = tmp;
                    tmp = t(:,3)';
                    tmp = tmp(~cellfun(@isempty,tmp));
                    stats_struct(cnts).anova_preds_df =tmp;
                    stats_struct(cnts).mod_preds_p = stats_out.Coefficients.pValue;
                    stats_struct(cnts).mod_preds_terms = stats_out.Coefficients.Name';
                    stats_struct(cnts).mod_preds_stat = stats_out.Coefficients.tStat;
                    stats_struct(cnts).mod_preds_coeff = stats_out.Coefficients.Estimate;
                    stats_struct(cnts).multi_comp_terms = gnames';
                    stats_struct(cnts).multi_comp_t1_t2 = comparisons(:,1:2);
                    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(comparisons(:,6));
                    stats_struct(cnts).multi_comp_p = adj_p;
                    stats_struct(cnts).multi_comp_coeff = comparisons(:,4);
                    stats_struct(cnts).multi_comp_lci_uci = [comparisons(:,3),comparisons(:,5)];
                    stats_struct(cnts).mod_r2 = stats_out.Rsquared.Adjusted;
                    stats_struct(cnts).norm_test_p = norm_p;
                    stats_struct(cnts).norm_test_h = norm_h;
                    stats_struct(cnts).effect_size = alt_f2;
                    stats_struct(cnts).effect_size_calc = '(R2_full-R2_null)/(1-R2_full)';
                    cnts = cnts + 1;
                    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                    %## PLOT =============================================================== %%
                    multi_p = adj_p;
                    aa_p = t(:,7)';
                    mod_coeff = stats_out.Coefficients.Estimate;
                    anova_p = aa_p{2};
                    anova_grp_p = aa_p{3};
                    terr_p = [1,stats_out.Coefficients.pValue(2:4)']';
                    grp_p = [0,stats_out.Coefficients.pValue(3:4)']';
                    speed_r2 = stats_out.Rsquared.Adjusted;
                    speed_xvals = (0:5)*0.25;
                    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
                            'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
                            'pvals',{{(terr_p),(terr_p),(terr_p)}},...
                            'pvals_pairs',{{{[1,1],[1,2],[1,3],[1,4]},...
                                {[1,1],[1,2],[1,3],[1,4]},...
                                {[1,1],[1,2],[1,3],[1,4]}}},...
                            'pvals_grp',{num2cell(adj_p)},...
                            'pvals_grp_pairs',{num2cell(comparisons(:,1:2),2)},...
                            'regress_pval',{{}},...
                            'regress_line',{{}},...
                            'r2_coeff',{[]},...
                            'regress_xvals',speed_xvals,...
                            'subject_char',[],... % this option when filled prints removal of nan() info
                            'group_order',categorical({''}),...
                            'do_include_intercept',false,...
                            'bracket_yshift_perc',1,...
                            'bracket_y_perc',1,...
                            'bracket_rawshifty_upper',0.2,...
                            'bracket_rawshifty_lower',0,...
                            'grp_sig_offset_x',0,...
                            'grp_sig_offset_y',0.1); 
                case 2
                    %## RANOVA
                    tt_out = unstack(tmptmp_table,measure_plot,'cond_char');
                    mod_out = sprintf('x0_25,x0_5,x0_75,x1~group_id+1');
                    tmptmp = table([1,2,3,4]','VariableNames',{'speed'});
                    %- run test
                    rm_out = fitrm(tt_out,mod_out,'WithinDesign',tmptmp);
                    [ran_out,A,C,D] = ranova(rm_out);
                    %- test coeffs
                    disp(cat(1,A{:})*table2array(rm_out.Coefficients)*C)
                    [ran_mult] = rm_out.multcompare('group_id');
                    %- test normality
                    % [norm_h,norm_p] = lillietest(stats_out.Residuals.Raw);
                    sphr = rm_out.mauchly(C);
                    sphr_corr = rm_out.epsilon(C);
                    %## RANOVA
                    %- populate struct
                    stats_struct(cnts).stat_test_mod = rm_out.BetweenModel;
                    stats_struct(cnts).measure_tag = categorical({measure_plot});
                    stats_struct(cnts).measure_char = measure_plot;
                    stats_struct(cnts).design_tag = categorical(des_i);
                    stats_struct(cnts).design_num = des_i;
                    stats_struct(cnts).cluster_num = cl_i;
                    stats_struct(cnts).mod_tag = categorical({'ranova_2'});
                    stats_struct(cnts).mod_resp_terms = rm_out.ResponseNames;
                    %- anova
                    stats_struct(cnts).anova_preds_terms = ran_out.Properties.RowNames';
                    stats_struct(cnts).anova_preds_df = ran_out.DF';
                    %- liberal
                    % tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
                    % stats_struct(cnts).anova_preds_p = tmp(1:2)';
                    %- conservative
                    tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
                    stats_struct(cnts).anova_preds_p = tmp(1:2)';
                    %-
                    stats_struct(cnts).anova_preds_stat = ran_out.F;
                    %- model
                    stats_struct(cnts).mod_preds_coeff = rm_out.Coefficients; % remove from table before xlsx
                    stats_struct(cnts).mod_preds_terms = rm_out.BetweenFactorNames';
                    stats_struct(cnts).multi_comp_terms = ran_mult.Properties.VariableNames(1:2);
                    stats_struct(cnts).multi_comp_p = ran_mult.pValue;
                    stats_struct(cnts).multi_comp_t1_t2 = num2cell(ran_mult{:,1:2},2)';
                    stats_struct(cnts).multi_comp_lci_uci = [ran_mult.Lower,ran_mult.Upper];
                    stats_struct(cnts).norm_test_p = norm_p;
                    stats_struct(cnts).norm_test_h = norm_h;
                    cnts = cnts + 1;
                    stats_struct(cnts) = DEF_STATS_TRACK_STRUCT;
                    %## PLOT =============================================================== %%
                    %## RANOVA VISUALIZE
                    % tmp = max([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
                    tmp = min([ran_out.pValue,ran_out.pValueGG,ran_out.pValueHF,ran_out.pValueLB],[],2);
                    comparisons = [1,2;1,3;2,3];
                    [h,crit_p,adj_ci_cvrg,adj_p] = fdr_bh(ran_mult.pValue);
                    adj_p = [adj_p(1),adj_p(2),adj_p(4)];
                    anova_p = 0;
                    % anova_grp_p = ran_out.pValue(2);
                    anova_grp_p = tmp(2);
                    % ran_means_1 = [rm_out.Coefficients{1,:}].*[0.25,0.5,0.75,1];
                    % ran_means_2 = [rm_out.Coefficients{2,:}].*[0.25,0.5,0.75,1]+ran_means_1;
                    % ran_means_3 = [rm_out.Coefficients{3,:}].*[0.25,0.5,0.75,1]+ran_means_1;
                    ran_means_1 = [rm_out.Coefficients{1,:}];
                    ran_means_2 = [rm_out.Coefficients{2,:}]+ran_means_1;
                    ran_means_3 = [rm_out.Coefficients{3,:}]+ran_means_1;
                    if anova_grp_p > 0.01 && anova_grp_p < 0.05
                        str = {sprintf('* pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
                    elseif anova_grp_p <= 0.01 && anova_grp_p > 0.001
                        str = {sprintf('** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
                    elseif anova_grp_p <= 0.001
                        str = {sprintf('*** pVal=%0.3g\nMS_{error}=%0.2g',tmp(2),ran_out.MeanSq(3)),'',''};
                    else
                        str = {'','',''};
                    end
                    STATS_STRUCT = struct('anova',{{anova_p,anova_p,anova_p}},...
                            'anova_grp',{{anova_grp_p,anova_grp_p,anova_grp_p}},...
                            'pvals',{{}},...
                            'pvals_pairs',{{}},...
                            'pvals_grp',{num2cell(adj_p)},...
                            'pvals_grp_pairs',{num2cell(comparisons,2)},...
                            'regress_pval',{{0,0,0}},...
                            'regress_line',{{ran_means_2,...
                                ran_means_3,...
                                ran_means_1}},...
                            'line_type',{'means'},...
                            'stats_char',{str},...
                            'display_stats_char',true,...
                            'bracket_yshift_perc',0.5,...
                            'bracket_y_perc',0.5,...
                            'bracket_rawshifty_upper',0,...
                            'bracket_rawshifty_lower',0,...
                            'grp_sig_offset_x',0,...
                            'grp_sig_offset_y',0);
            end
            %##
            % cl_i = double(string(clusters(j)));
            % temp_table = FOOOF_TABLE(FOOOF_TABLE.cluster_id == num2str(cl_i) & FOOOF_TABLE.design_id == num2str(des_i),:);
            %- ylim autoset
            % % YLIMS = [min(temp_table.(measure_plot))-0.5,max(temp_table.(measure_plot))+1];
            %## PLOT
            ax = axes();
            VIOLIN_PARAMS = {'width',0.1,...
                'ShowWhiskers',false,'ShowNotches',false,'ShowBox',true,...
                'ShowMedian',true,'Bandwidth',0.15,'QuartileStyle','shadow',...
                'HalfViolin','full','DataStyle','scatter','MarkerSize',8,...
                'EdgeColor',[0.5,0.5,0.5],'ViolinAlpha',{0.2 0.3}};
            PLOT_STRUCT = struct('color_map',color_dark,...
                'cond_labels',unique(tmptmp_table.cond_char),'group_labels',categorical(g_chars_subp),...
                'cond_offsets',cond_offsets,...
                'group_offsets',[0.125,0.475,0.812],...
                'group_lab_yoffset',-0.285,...
                'y_label',measure_plot,...
                'title',sprintf('Cluster %i',cl_i),'font_size',FONT_SIZE_VIO,'ylim',YLIMS,...
                'font_name','Arial','x_label',x_label,'do_combine_groups',false,...
                'regresslab_txt_size',FONT_SIZE_VIO_REG);
            ax = group_violin(tmptmp_table,measure_plot,'cond_char','group_id',...
                ax,...
                'VIOLIN_PARAMS',VIOLIN_PARAMS,...
                'PLOT_STRUCT',PLOT_STRUCT,...
                'STATS_STRUCT',STATS_STRUCT);
            set(ax,'OuterPosition',[0,0,1,1]);
            set(ax,'Position',[AX_INIT_HORIZ+horiz_shift,AX_INIT_VERT_VIO+vert_shift,AX_W*IM_RESIZE,AX_H*IM_RESIZE]);  %[left bottom width height]
            ax.Children(1).FontSize = GROUP_LAB_FONTSIZE;
            ax.Children(2).FontSize = GROUP_LAB_FONTSIZE;
            ax.Children(3).FontSize = GROUP_LAB_FONTSIZE;
            % ax.Children(1).Position(2) = GROUP_LAB_YOFFSET;
            % ax.Children(2).Position(2) = GROUP_LAB_YOFFSET;
            % ax.Children(3).Position(2) = GROUP_LAB_YOFFSET;
            yt = yticks(ax);
            xlh = xlabel(ax,PLOT_STRUCT.x_label,'Units','normalized','FontSize',XLAB_FONTSIZE,'FontWeight',XLAB_FONTWEIGHT);
            pos1=get(xlh,'Position');
            pos1(1,2)=pos1(1,2)+XLABEL_OFFSET;
            set(xlh,'Position',pos1);
            if hz ~= 1
                ylabel('');
            end
            %##
            if hz < HZ_DIM
                horiz_shift = horiz_shift + AX_W*IM_RESIZE + AX_HORIZ_SHIFT*IM_RESIZE;
            else
                vert_shift = vert_shift - (AX_H*IM_RESIZE + AX_VERT_SHIFT*IM_RESIZE);
                horiz_shift = 0;
                hz = 0;
            end
            hz = hz + 1;
        end
        drawnow;
        %## SAVE
        hold off;
        savefig(fig,[save_dir filesep sprintf('%s_violin_group_plot_ranova_des%i.fig',measure_plot,des_i)])
        exportgraphics(fig,[save_dir filesep sprintf('%s_violin_group_plot_ranova_des%i.png',measure_plot,des_i)],'Resolution',300)
        close(fig);
    end
end
stats_struct_table = struct2table(stats_struct);
%}
%% Perform time series stats on the flattened curve
%## STATS
try
    STUDY.etc = rmfield(STUDY.etc,'statistics');
end
STUDY = pop_statparams(STUDY,...
    'groupstats','off',...
    'condstats','on',...
    'method','perm',...
    'singletrials','off',...
    'mode','fieldtrip',...
    'fieldtripalpha',NaN,...
    'fieldtripmethod','montecarlo',...
    'fieldtripmcorrect','fdr',...
    'fieldtripnaccu',4000);
stats = STUDY.etc.statistics;
stats.paired{1} = 'on'; % Condition stats
stats.paired{2} = 'off'; % Group stats
%% ===================================================================== %%
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%- fooof_diff_store needs to be freq x subject, and condition by row
design_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
%-
pcond = cell(length(designs),1);
pgroup = cell(length(designs),1);
pinter = cell(length(designs),1);
statcond = cell(length(designs),1);
statgroup = cell(length(designs),1);
statinter = cell(length(designs),1);
cnt = 1;

for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        %## Ho : all samples come from the same distribution
        %## Ha : all samples come from different distributions
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(fooof_diff_store{des_i}{cl_i}, stats);
        %- 
        clust_t{cnt} = cl_i;
        design_t{cnt} = des_i;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter;
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond{des_i}{cl_i} = temp_pcond;
        pgroup{des_i}{cl_i} = temp_pcond;
        pinter{des_i}{cl_i} = temp_pinter;
        statcond{des_i}{cl_i} = temp_statcond;
        statgroup{des_i}{cl_i} = temp_statgroup;
        statinter{des_i}{cl_i} = temp_statinter;
        %%%%%
        for k0 = 1:length(pcond{des_i}{cl_i})
            pcond{des_i}{cl_i}{k0}(:,2) = pcond{des_i}{cl_i}{k0}(:,1)<0.05;    
            pcond_test_t{cnt} = pcond{des_i}{cl_i}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup{des_i}{cl_i})
            if ~isempty(pgroup{des_i}{cl_i}{k0})
                pgroup{des_i}{cl_i}{k0}(:,2) = pgroup{des_i}{cl_i}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup{des_i}{cl_i}{k0}(:,1)<0.05;    
            end
        end
        for k0 = 1:length(pinter{des_i}{cl_i})
            if ~isempty(pinter{des_i}{cl_i}{k0})
                pinter{des_i}{cl_i}{k0}(:,2) = pinter{des_i}{cl_i}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter{des_i}{cl_i}{k0}(:,1)<0.05;    
            end
        end
        cnt = cnt + 1;
    end
end
% table_out = table(design_t,clust_t,pcond_t,pcond_test_t,pgroup_t,pgroup_test_t,pinter_t,pinter_test_t,statcond_t,statgroup_t,statinter_t);
% save([save_dir filesep 'fooof_psd_stats.mat'],'table_out');
save([save_dir filesep 'fooof_pcond.mat'],'pcond');
%% ===================================================================== %%
% fooof_diff_store needs to be freq x subject, and condition by row
freq_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
% subj_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
% cond_t = cell(2*(length(fooof_diff_store{g})-3+1),1);
design_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
clust_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pcond_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pgroup_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
pinter_test_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statcond_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statgroup_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
statinter_t = cell(2*(length(fooof_diff_store{1})-3+1),1);
%-
designs = unique(FOOOF_TABLE.design_id);
clusters = unique(FOOOF_TABLE.cluster_id);
groups = unique(FOOOF_TABLE.group_id);
%-
pcond_org = cell(length(designs),1);
pgroup_org = cell(length(designs),1);
pinter_org = cell(length(designs),1);
statcond_org = cell(length(designs),1);
statgroup_org = cell(length(designs),1);
statinter_org = cell(length(designs),1);
cnt = 1;
for i = 1:length(designs)
    des_i = double(string(designs(i)));
    for j = 1:length(clusters)
        cl_i = double(string(clusters(j)));
        [temp_pcond, temp_pgroup, temp_pinter, temp_statcond, temp_statgroup, temp_statinter] = std_stat(spec_data_original{des_i}{cl_i}, stats);
        %- 
        clust_t{cnt} = cl_i;
        design_t{cnt} = des_i;
        pcond_t{cnt} = temp_pcond;
        pgroup_t{cnt} = temp_pgroup;
        pinter_t{cnt} = temp_pinter; %{temp_pinter};
        statcond_t{cnt} = temp_statcond;
        statgroup_t{cnt} = temp_statgroup;
        statinter_t{cnt} = temp_statinter;
        %-
        pcond_org{des_i}{cl_i} = temp_pcond;
        pgroup_org{des_i}{cl_i} = temp_pcond;
        pinter_org{des_i}{cl_i} = temp_pinter;
        statcond_org{des_i}{cl_i} = temp_statcond;
        statgroup_org{des_i}{cl_i} = temp_statgroup;
        statinter_org{des_i}{cl_i} = temp_statinter;
        for k0 = 1:length(pcond_org{des_i}{cl_i})
            pcond_org{des_i}{cl_i}{k0}(:,2) = pcond_org{des_i}{cl_i}{k0}(:,1)<0.05; 
            pcond_test_t{cnt} = pcond_org{des_i}{cl_i}{k0}(:,1)<0.05;    
        end
        for k0 = 1:length(pgroup_org{des_i}{cl_i})
            if ~isempty(pgroup_org{des_i}{cl_i}{k0})
                pgroup_org{des_i}{cl_i}{k0}(:,2) = pgroup_org{des_i}{cl_i}{k0}(:,1)<0.05;  
                pgroup_test_t{cnt} = pgroup_org{des_i}{cl_i}{k0}(:,1)<0.05;   
            end
        end
        for k0 = 1:length(pinter_org{des_i}{cl_i})
            if ~isempty(pinter_org{des_i}{cl_i}{k0})
                pinter_org{des_i}{cl_i}{k0}(:,2) = pinter_org{des_i}{cl_i}{k0}(:,1)<0.05;
                pinter_test_t{cnt} = pinter_org{des_i}{cl_i}{k0}(:,1)<0.05; %{pinter_org{g}{k}{k0}(:,1)<0.05};  
            end
        end            
        cnt=cnt+1;
    end
end
save([save_dir filesep 'org_pcond.mat'],'pcond_org');
%% Sanity check - time series plots from aperiodic subtraction
% TERRAIN_DES_ID = 1;
% SPEED_VALS = {'0.25','0.5','0.75','1.0';
%               '0p25','0p5','0p75','1p0'};
% TERRAIN_VALS = {'flat','low','med','high'};
% COLORS_MAPS_TERRAIN = linspecer(4);
% custom_yellow = [254,223,0]/255;
% COLORS_MAPS_TERRAIN = [COLORS_MAPS_TERRAIN(3,:);custom_yellow;COLORS_MAPS_TERRAIN(4,:);COLORS_MAPS_TERRAIN(2,:)];
% COLOR_MAPS_SPEED = linspecer(4*3);
% COLOR_MAPS_SPEED = [COLOR_MAPS_SPEED(1,:);COLOR_MAPS_SPEED(2,:);COLOR_MAPS_SPEED(3,:);COLOR_MAPS_SPEED(4,:)];
% %##
% %-
% designs = unique(FOOOF_TABLE.design_id);
% clusters = unique(FOOOF_TABLE.cluster_id);
% groups = unique(FOOOF_TABLE.group_id);
% %-
% for i = 1:length(designs)
%     des_i = double(string(designs(i)));
%     for j = 1:length(clusters)
%         cl_i = double(string(clusters(j)));
%         data_min = min([fooof_diff_store{des_i}{cl_i}{:}],[],'all');
%         data_max = max([fooof_diff_store{des_i}{cl_i}{:}],[],'all');
%         cnt = 1;
%         figure('color','white');
%         for k = 1:size(fooof_diff_store{des_i}{cl_i},1)
%             for l = 1:size(fooof_diff_store{des_i}{cl_i},2)
%                 hold on;
%                 subplot(size(fooof_diff_store{des_i}{cl_i},1),size(fooof_diff_store{des_i}{cl_i},2),cnt)
%                 data = fooof_diff_store{des_i}{cl_i}{k,l};
%                 switch des_i
%                     case 1
%                         color_dark = COLORS_MAPS_TERRAIN;
%                         color_light = COLORS_MAPS_TERRAIN;
%                         GROUP_CMAP_OFFSET = [0,0.1,0.1];
%                         xtick_label_g = {'flat','low','med','high'};
%                     case 2
%                         color_dark = COLOR_MAPS_SPEED;
%                         color_light = COLOR_MAPS_SPEED+0.15;
%                         GROUP_CMAP_OFFSET = [0.15,0,0];
%                         xtick_label_g = c_chars; %{'0.25','0.50','0.75','1.0'};
%                 end
%                 plot(fooof_freq,data,'color',color_dark(k,:));
%                 ylabel('log10(Power)')
%                 ylim([data_min data_max]);
%                 xlabel('Frequency(Hz)');
%                 % title(['Cluster ',num2str(cl_i)]);
%                 cnt = cnt + 1;
%             end
%         end
%         hold off;
%         drawnow;
%     end
% end
