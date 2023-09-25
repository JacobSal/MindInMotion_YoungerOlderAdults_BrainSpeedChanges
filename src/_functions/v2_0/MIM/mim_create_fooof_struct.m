function [fooof_group_results_org,fooof_diff_store] = mim_create_fooof_struct(STUDY,save_dir,varargin)
%MIM_CREATE_FOOOF_STRUCT Summary of this function goes here
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
%-
SPEC_FPATHS = cell(length(STUDY.design),1);
for i = 1:length(STUDY.design)
    SPEC_FPATHS{i,:} = {STUDY.etc.mim_gen_ersp_data([STUDY.etc.mim_gen_ersp_data.des_ind] ==i).spec_fpaths};
end
SPEC_FPATHS = cat(1,SPEC_FPATHS{:});
%-
CLUSTER_PICKS = unique([STUDY.etc.mim_gen_ersp_data.clust_ind_cl]);
%-
SETTINGS = struct();  % Use defaults
SETTINGS.peak_width_limits = [1 8];%default [1 6] % Amanda used [1 8] in her paper
SETTINGS.min_peak_height = 0.05;
% settings.peak_threshold = 2;
SETTINGS.max_n_peaks = 3; % originally set to be 2 - 2023-06-07
% the settings are consitent with fooof on github
f_range = [3, 40];
theta_band = [4, 8];
alpha_band = [8 12];
beta_band  = [12 30];
%- custom params
colormap_ersp = linspecer; %othercolor('RdYlBu11');
% colormap_ersp = colormap_ersp(end:-1:1,:);
colormap(colormap_ersp);
close all;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'SETTINGS',SETTINGS,@isstruct);
addParameter(p,'CLUSTER_PICKS',CLUSTER_PICKS,@isnumeric);
addParameter(p,'SPEC_FPATHS',SPEC_FPATHS,@iscell);
parse(p,STUDY,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%## MAKE DIRS
if ~exist(save_dir,'dir')
    warning('''save_dir'' is not a directory... fooof data will not be save');
    save_dir = [];
end
SETTINGS = p.Results.SETTINGS;
CLUSTER_PICKS = p.Results.CLUSTER_PICKS;
SPEC_FPATHS = p.Results.SPEC_FPATHS;
%% ===================================================================== %%
%## Calculate scalp maps 
% if ~isfield(STUDY.cluster,'topo'), STUDY.cluster(1).topo = []; end
% for clus = CLUSTER_PICKS % For each cluster requested
%     if isempty(STUDY.cluster(clus).topo)
%         STUDY = std_readtopoclust_CL(STUDY,ALLEEG, clus);% Using this custom modified code to allow taking average within participant for each cluster
%     end
% end
%## VALIDATION PRINTOUTS
cl_inds = [STUDY.etc.mim_gen_ersp_data.clust_ind_cl];
des_inds = [STUDY.etc.mim_gen_ersp_data.des_ind];
fprintf('Clusters to plot: '); fprintf('%i,', CLUSTER_PICKS); fprintf('\n');
fprintf('Designs to plot: '); fprintf('%i,', unique(des_inds)); fprintf('\n');
fprintf('Statistics Parameters:\n');
disp(STUDY.etc.statistics)
fprintf('Statistics Fieldtrip Parameters:\n');
disp(STUDY.etc.statistics.fieldtrip)
fprintf('Statistics EEGLAB Parameters:\n');
disp(STUDY.etc.statistics.eeglab)
fprintf('ERSP Parameters:\n');
disp(STUDY.etc.erspparams)
%## FOOOF SETUP & LOOP
fooof_results = cell(length(STUDY.design),1);
fooof_group_results_org = cell(length(STUDY.design),1);
fooof_apfit_store = cell(length(STUDY.design),1);
fooof_diff_store = cell(length(STUDY.design),1);
spec_data_original = cell(length(STUDY.design),1);
for des_i = 1:length(STUDY.design)
    cond_test = STUDY.design(des_i).variable(1).value;
    fprintf('Current design: %i\n',STUDY.currentdesign);
    fprintf('Running Design: '); fprintf('%s,',cond_test{1:end-1}); fprintf('%s',cond_test{end}); fprintf('\n');
%     STUDY.currentdesign = des_i;
    for j = 1:length(CLUSTER_PICKS)
        cluster_i = CLUSTER_PICKS(j);
        fname = SPEC_FPATHS{des_i,j};
        fpath = strjoin(fname(1:end-1),'/');
        fname = fname{end};
        tmp = par_load(fpath,fname);
        specdata = {tmp.specdata}';
        specfreqs = {tmp.specfreqs}';
        specfreqs = specfreqs{1};
        %## (DEBUG) ALTERNATIVE
%         if inside_load_flag
%             cluster_load_ind = find(logical(cl_inds == cluster_i) & logical(des_inds == des_i));
%             fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).spec_fpaths,'/');
% %             fname = strsplit(STUDY.etc.mim_gen_ersp_data(cluster_load_ind).spec_ss_fpaths,'/');
%             fpath = strjoin(fname(1:end-1),'/');
%             fname = fname{end};
%             tmp = par_load(fpath,fname);
%             specdata = {tmp.specdata}';
%             specfreqs = {tmp.specfreqs}';
%             specfreqs = specfreqs{1};
%         end
        %% Run fooof
        %- Input should be in linear spacing   
        i_ind = 0;
        for group = 1:size(specdata,2) % in case there is young and old adult group
            for cond = 1:size(specdata,1) % different level of terrains
                specdata_nolog = 10.^(specdata{cond,group}/10);
                % Run FOOOF
                return_model = true;
                for subj_i = 1:size(specdata{cond,group},2)
                    fooof_results{des_i}{cond,group}{subj_i} = fooof(specfreqs, specdata_nolog(:,subj_i), f_range, SETTINGS, return_model);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).subID = STUDY.cluster(cluster_i).sets(subj_i);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).compID = STUDY.cluster(cluster_i).comps(subj_i);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).study = des_i;%1 = terrain, 2 = speed
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).cond = cond;
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).group = group;
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).cluster = cluster_i;
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).aperiodic_exp = fooof_results{des_i}{cond,group}{subj_i}.aperiodic_params(2);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).aperiodic_offset = fooof_results{des_i}{cond,group}{subj_i}.aperiodic_params(1);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).central_freq = fooof_results{des_i}{cond,group}{subj_i}.peak_params(:,1);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).power = fooof_results{des_i}{cond,group}{subj_i}.peak_params(:,2);
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).r_squared = fooof_results{des_i}{cond,group}{subj_i}.r_squared;

                    %## Compute Average Power After Flatten Curve
%                     fooof_diff = fooof_results{g}{cond,group}{i}.power_spectrum - fooof_results{g}{cond,group}{i}.ap_fit;
                    % Super important, the output is already logged, the
                    % only difference is the magnitude by 10
                    fooof_diff = 10*(fooof_results{des_i}{cond,group}{subj_i}.power_spectrum) - 10*(fooof_results{des_i}{cond,group}{subj_i}.ap_fit);
                    fooof_freq = fooof_results{des_i}{cond,group}{subj_i}.freqs;
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).theta_avg_power = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).alpha_avg_power = mean(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{des_i}{cluster_i}(i_ind + subj_i).beta_avg_power = mean(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));

                    %- data structure needs to be freq x subject
                    fooof_diff_store{des_i}{cluster_i}{cond}(:,subj_i) = fooof_diff';
                    fooof_apfit_store{des_i}{cluster_i}{cond}(:,subj_i) = 10*(fooof_results{des_i}{cond,group}{subj_i}.ap_fit);
                    
                    %- store original spec data
                    spec_data_original{des_i}{cluster_i}{cond} = specdata{cond,group}(specfreqs >= f_range(1) & specfreqs <= f_range(2),:);
                end
                i_ind = i_ind + size(specdata{cond,group},2);
            end
        end
    end
end
for g = 1:2
    for k = 3:length(fooof_group_results_org{g})
        for i = 1:length(fooof_group_results_org{g}{k})
            if ~isempty(fooof_group_results_org{g}{k}(i).central_freq)
                fooof_group_results_org{g}{k}(i).theta = [];
                fooof_group_results_org{g}{k}(i).alpha = [];
                fooof_group_results_org{g}{k}(i).beta = [];
                for j = 1:length(fooof_group_results_org{g}{k}(i).central_freq)
                    cf = fooof_group_results_org{g}{k}(i).central_freq(j);
                    if cf > theta_band(1) & cf <= theta_band(2)
                        fooof_group_results_org{g}{k}(i).theta = [fooof_group_results_org{g}{k}(i).theta; cf fooof_group_results_org{g}{k}(i).power(j)];
                    elseif cf > alpha_band(1) & cf <= alpha_band(2)
                        fooof_group_results_org{g}{k}(i).alpha = [fooof_group_results_org{g}{k}(i).alpha; cf fooof_group_results_org{g}{k}(i).power(j)];
                    elseif cf > beta_band(1) & cf <= beta_band(2)
                        fooof_group_results_org{g}{k}(i).beta = [fooof_group_results_org{g}{k}(i).beta; cf fooof_group_results_org{g}{k}(i).power(j)];
                    end
                end
                if length(fooof_group_results_org{g}{k}(i).theta) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).theta(:,1)-6));
                    temp_power = fooof_group_results_org{g}{k}(i).theta(:,2);
                    fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).theta(indx,1) temp_power(indx)];
                end
                if length(fooof_group_results_org{g}{k}(i).alpha) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).alpha(:,1)-10));
                    temp_power = fooof_group_results_org{g}{k}(i).alpha(:,2);
                    fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).alpha(indx,1) temp_power(indx)];
                end
                if length(fooof_group_results_org{g}{k}(i).beta) > 1
                    [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).beta(:,1)-20));
                    temp_power = fooof_group_results_org{g}{k}(i).beta(:,2);
                    fooof_group_results_org{g}{k}(i).beta_center = [fooof_group_results_org{g}{k}(i).beta(indx,1) temp_power(indx)];
                end
            end
        end
    end
end
if ~isempty(save_dir)
    par_save(fooof_group_results_org,[save_dir filesep 'fooof_results_summary.mat']);
    par_save(fooof_diff_store,[save_dir filesep 'fooof_diff_store.mat']);
end
% keyboard
end