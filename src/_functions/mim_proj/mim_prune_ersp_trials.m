function [allersp_pruned_avg,allersp_pruned_all,trial_nums] = mim_prune_ersp_trials(STUDY,allerspdata,cluster_i,des_i,varargin)
%MIM_PRUNE_ERSP_TRIALS Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic;
%## DEFINE DEFAULTS
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'allerspdata',@iscell);
addRequired(p,'des_i',@isnumeric);
addRequired(p,'cluster_i',@isnumeric);
%## OPTIONAL
%## PARAMETER
parse(p,STUDY,allerspdata,cluster_i,des_i,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
%% ===================================================================== %%
PERCENTILE_VALS = [10,25,50,75,90];
PERCENTILE_ITER = 5;
COND_CHARS = STUDY.design(des_i).variable.value;
BAD_TRIAL_PERC_CUTOFF = 0.25;
BAD_SUBJ_PERC_CUTOFF = 0.30;
avg_fcn = @mean;
%-
allersp_pruned_all = cell(size(allerspdata));
allersp_pruned_avg = cell(size(allerspdata));
for group_i = 1:size(allerspdata,2)
    for cond_i = 1:size(allerspdata,1)
        trial_nums = zeros(length(unique(STUDY.cluster(cluster_i).sets)),1);
        subjs = unique(STUDY.cluster(cluster_i).sets);
        for i = 1:length(subjs)
            subj_i = subjs(i);
            subj_epoch_n = strcmp({STUDY.datasetinfo(subj_i).trialinfo.cond},COND_CHARS{cond_i});
            trial_nums(i) = sum(subj_epoch_n); %find(subj_epoch_n,1,'last') - find(subj_epoch_n,1,'last');
        end
        cnt1 = 1;
        cnt2 = trial_nums(1);
        %- percentile across all subjects in condition.
        prct = prctile(allerspdata{cond_i,group_i},PERCENTILE_VALS,[3,1,2]);
        time_freq_area = size(allerspdata{cond_i,group_i},1)*size(allerspdata{cond_i,group_i},2);
        prunes = cell(length(trial_nums),1);
        prunes_avg = cell(length(trial_nums),1);
        bad_subjs = zeros(length(trial_nums),1);
        for i = 2:length(trial_nums)+1
            tmp = allerspdata{cond_i,group_i}(:,:,cnt1:cnt2);
            %- percentile across all trials in a subject.
%             prct = prctile(tmp,PERCENTILE_VALS,[3,1,2]);
            logi = zeros(size(tmp,3),2);
            bad_area_p = zeros(size(tmp,3),2);
%             bad_trials = zeros(size(tmp,3),1);
            for trial_i = 1:size(tmp,3)
                logi(trial_i,1) = sum(tmp(:,:,trial_i) < -prct(PERCENTILE_ITER),'all');
                logi(trial_i,2) = sum(tmp(:,:,trial_i) > prct(PERCENTILE_ITER),'all');
                bad_area_p(trial_i,1) = logi(trial_i,1)/time_freq_area;
                bad_area_p(trial_i,2) = logi(trial_i,2)/time_freq_area;
            end
            bad_trials = any(bad_area_p > BAD_TRIAL_PERC_CUTOFF,2);
            bad_subjs(i-1) = sum(bad_trials)/length(bad_trials) > BAD_SUBJ_PERC_CUTOFF;
            prunes{i-1} = tmp(:,:,~bad_trials);
            prunes_avg{i-1} = feval(avg_fcn,prunes{i-1},3);
            if i <= length(trial_nums)
                cnt1 = cnt2+1;
                cnt2 = cnt2+trial_nums(i);
                disp(cnt1);
                disp(cnt2);
            end
        end
        allersp_pruned_all{cond_i,group_i} = cat(3,prunes{:});
        allersp_pruned_avg{cond_i,group_i} = cat(3,prunes_avg{:});
    end
end
%## TIME
toc;
end

