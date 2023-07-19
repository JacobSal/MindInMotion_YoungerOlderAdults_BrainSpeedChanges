function [spec_fpaths,spec_ss_fpaths,ersp_fpaths,ersp_norm_fpaths,ersp_normcb_fpaths] = mim_gen_ersp_data(STUDY,ALLEEG,warping_times,...
    ersp_load_params,save_dir,varargin)
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
DO_SPEC_CALC = true;
DO_ERSP_CALC = true;
STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
SPEC_PARAMS = struct('freqrange',[1,100],...
    'subject','',...
    'comps','all');
ERSP_PARAMS = struct('subbaseline','on',...
    'timerange',[warping_times(1) warping_times(5)],...
    'ersplim',[-2,2],...
    'freqrange',[1,200]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'warping_times',@isnumeric);
addRequired(p,'ersp_load_params',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'DO_SPEC_CALC',DO_SPEC_CALC,@islogical);
addParameter(p,'DO_ERSP_CALC',DO_ERSP_CALC,@islogical);
addParameter(p,'ERSP_PARAMS',ERSP_PARAMS,@isstruct);
addParameter(p,'SPEC_PARAMS',SPEC_PARAMS,@isstruct);
addParameter(p,'STAT_PARAMS',STAT_PARAMS,@isstruct);
parse(p,STUDY,ALLEEG,warping_times,ersp_load_params,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
ERSP_PARAMS = p.Results.ERSP_PARAMS;
SPEC_PARAMS = p.Results.SPEC_PARAMS;
STAT_PARAMS = p.Results.STAT_PARAMS;
DO_SPEC_CALC = p.Results.DO_SPEC_CALC;
DO_ERSP_CALC = p.Results.DO_ERSP_CALC;
%% ===================================================================== %%
%## PARAMS SETUP
tmpSTUDY = pop_statparams(STUDY, 'condstats', STAT_PARAMS.condstats,...
        'method',STAT_PARAMS.method,...
        'singletrials',STAT_PARAMS.singletrials,'mode',STAT_PARAMS.mode,...
        'fieldtripalpha',STAT_PARAMS.fieldtripalpha,'fieldtripmethod',STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',STAT_PARAMS.fieldtripnaccu);
tmpSTUDY_commonbase = pop_erspparams(tmpSTUDY, 'subbaseline',ERSP_PARAMS.subbaseline,...
        'timerange',ERSP_PARAMS.timerange, 'ersplim',ERSP_PARAMS.ersplim);  % 'subbaseline' - ['on'|'off'] subtract the same baseline across conditions for ERSP     
[~,main_cl_inds,~,valid_clusters] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICKS = main_cl_inds;
%-
spec_fpaths = cell(length(STUDY.design),length(CLUSTER_PICKS));
spec_ss_fpaths = cell(length(STUDY.design),length(CLUSTER_PICKS));
ersp_fpaths = cell(length(STUDY.design),length(CLUSTER_PICKS));
ersp_norm_fpaths = cell(length(STUDY.design),length(CLUSTER_PICKS));
ersp_normcb_fpaths = cell(length(STUDY.design),length(CLUSTER_PICKS));
%% SPECTRUM CALCULATIONS
if DO_SPEC_CALC
    for des_i = 1:length(STUDY.design)
        store_1 = cell(1,length(CLUSTER_PICKS));
        store_2 = cell(1,length(CLUSTER_PICKS));
        STUDY.currentdesign = des_i;
        design_char = [];
        for i = 1:length(STUDY.design(des_i).variable)
            if i==1
                design_char = [STUDY.design(des_i).variable(1).value{:}];
            else
                design_char = [design_char '_' STUDY.design(des_i).variable(i).value{:}];
            end
        end
        %## Compute Specs
        parfor i = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(i);
            %##
            fprintf('Computing specdata for cluster %i...\n',cluster_i)
            tic
            try
                [~,specdata,specfreqs,pgroup,pcond,pinter] = std_specplot(STUDY, ALLEEG, 'clusters',cluster_i,...
                    'comps',SPEC_PARAMS.comps,'subject',SPEC_PARAMS.subject,'freqrange',SPEC_PARAMS.freqrange);
                %- save dat
                spec_data = struct('specdata',specdata,'specfreqs',specfreqs,...
                    'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
                par_save(spec_data,save_dir,sprintf('spec_data_cl%i_%s.mat',cluster_i,design_char));
                store_1{i} = [save_dir filesep sprintf('spec_data_cl%i_%s.mat',cluster_i,design_char)];
                %- save fig
                fig_i = get(groot,'CurrentFigure');
    %             fig_i.Position = [500 300 1480 920];
                saveas(fig_i,[save_dir filesep sprintf('spec_plot_cl%i_%s_allcomps',cluster_i,design_char)])
            catch e
                fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'stack. %s\n'],e.identifier,e.message,getReport(e));
            end
            toc
            %##
            fprintf('Computing specdata with substracted subject mean for cluster %i...\n',cluster_i)
            tic
            try
                [~,specdata,specfreqs,pgroup,pcond,pinter] = std_specplot(STUDY, ALLEEG, 'clusters',cluster_i,...
                    'comps',SPEC_PARAMS.comps,'subject',SPEC_PARAMS.subject,'freqrange',SPEC_PARAMS.freqrange,'subtractsubjectmean','on');
                %- save dat
                spec_data = struct('specdata',specdata,'specfreqs',specfreqs,...
                    'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
                par_save(spec_data,save_dir,sprintf('spec_data_cl%i_%s_subtractmean.mat',cluster_i,design_char));
                store_2{i} = [save_dir filesep sprintf('spec_data_cl%i_%s_subtractmean.mat',cluster_i,design_char)];
                %- save fig
                fig_i = get(groot,'CurrentFigure');
    %             fig_i.Position = [500 300 1480 920];
                saveas(fig_i,[save_dir filesep sprintf('spec_plot_cl%i_%s_allcomps_subtractsubjectmean',cluster_i,design_char)])
            catch e
                fprintf(['error. identifier: %s\n',...
                     'error. %s\n',...
                     'stack. %s\n'],e.identifier,e.message,getReport(e));
            end
            toc 
            fprintf('done.\n')
            close all
        end
        spec_fpaths(des_i,:) = store_1;
        spec_ss_fpaths(des_i,:) = store_2;
    end
end
%% ERSP CALCULATIONS
% NOTE (07/18/2023) std_ersplot_customParams.m has edits (find:'Chang') as
% well as a call to another edited function std_readata_customParams.m
if DO_ERSP_CALC
    for des_i = 1:length(STUDY.design)
        store_1 = cell(1,length(CLUSTER_PICKS));
        store_2 = cell(1,length(CLUSTER_PICKS));
        store_3 = cell(1,length(CLUSTER_PICKS));
        STUDY.currentdesign = des_i;
        design_char = [];
        for i = 1:length(STUDY.design(des_i).variable)
            if i==1
                design_char = [STUDY.design(des_i).variable(1).value{:}];
            else
                design_char = [design_char '_' STUDY.design(des_i).variable(i).value{:}];
            end
        end
        parfor i = 1:length(CLUSTER_PICKS)
            cluster_i = CLUSTER_PICKS(i);
            % ----------------------------------------------------------- %
            %## ERSP calculation for no normalization 
            fprintf('Gathering ERSP without any normalization for cluster %i\n',cluster_i)
            tic
            [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot(STUDY,ALLEEG,'clusters',cluster_i,...
                'freqrange',ERSP_PARAMS.freqrange); %PLOTS THE ENTIRE cluster_i'th CLUSTER ,% did I apply subtractsubjectmean???, no correction at all with the current setting
            %- save dat
            ersp_data = struct('allerspdata',allerspdata,'alltimes',alltimes,'allfreqs',allfreqs,...
                'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
            par_save(ersp_data,save_dir,sprintf('ersp_data_cl%i_%s.mat',cluster_i,STUDY.design));
            store_1{i} = [save_dir filesep sprintf('ersp_data_cl%i_%s.mat',cluster_i,STUDY.design)];
            %- save fig
            fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
            saveas(fig_i,[save_dir filesep sprintf('allcomps_ersp_plot_cl%i_%s',cluster_i,design_char)])
            close(fig_i)
            toc
            % ----------------------------------------------------------- %
            %## ERSP calculation for normalization 
            fprintf('Gathering ERSP after baseline correction using times {%0.2f,%0.2f] for cluster %i\n',warping_times(1),warping_times(5),cluster_i)
            tic
            disp(tmpSTUDY.etc.statistics)
            disp(ersp_load_params.common_base);
            [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot_customParams(tmpSTUDY,ALLEEG,ersp_load_params.common_base,'clusters',cluster_i,...
                'freqrange',ERSP_PARAMS.freqrange);
            %- save dat
            ersp_data = struct('allerspdata',allerspdata,'alltimes',alltimes,'allfreqs',allfreqs,...
                'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
            par_save(ersp_data,save_dir,sprintf('ersp_data_cl%i_%s_subbaselined.mat',cluster_i,design_char));
            store_2{i} = [save_dir filesep sprintf('ersp_data_cl%i_%s_subbaselined.mat',cluster_i,design_char)];
            %- save fig
            fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
            saveas(fig_i,[save_dir filesep sprintf('ersp_plot_cl%i_%s_allcomps_subbaselined',cluster_i,design_char)])
            close(fig_i)
            toc
            % ----------------------------------------------------------- %
            disp('Gathering esrp after baseline correction using warp(1) to warp(5) and apply commonbase')
            tic
            disp(tmpSTUDY_commonbase.etc.statistics)
            disp(tmpSTUDY_commonbase.etc.erspparams)
            [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot_customParams(tmpSTUDY_commonbase,ALLEEG,ersp_load_params.common_base,'clusters',cluster_i,...
                'freqrange',ERSP_PARAMS.freqrange); 
            %- save dat
            ersp_data = struct('allerspdata',allerspdata,'alltimes',alltimes,'allfreqs',allfreqs,...
                'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
            par_save(ersp_data,save_dir,sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',cluster_i,design_char));
            store_3{i} = [save_dir filesep sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',cluster_i,design_char)];
            ersp
            %- save fig
            fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
            saveas(fig_i,[save_dir filesep sprintf('ersp_plot_cl%i_%s_allcomps_subbaselined_commonbase',cluster_i,design_char)])
            close(fig_i)
        end
        ersp_fpaths(des_i,:) = store_1;
        ersp_norm_fpaths(des_i,:) = store_2;
        ersp_normcb_fpaths(des_i,:) = store_3;
    end
end
STUDY.etc.mim_gen_ersp_data = struct('ersp_fpaths',ersp_fpaths,...
    'spec_fpaths',spec_fpaths,...
    'spec_ss_fpaths',spec_ss_fpaths,...
    'ersp_norm_fpaths',ersp_norm_fpaths,...
    'ersp_normcb_fpaths',ersp_normcb_fpaths);
end

