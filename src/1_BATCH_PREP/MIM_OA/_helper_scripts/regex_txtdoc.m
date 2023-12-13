
%%
fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\1_BATCH_PREP\MIM_OA\_hpg_logs';
fname = '13711340_MIM_prep.log';
fid = fileread([fpath filesep fname]);
reg_out = regexp(fid,'[+-]?\d*\.?\d* channels were rejected','match');
nums_out = cellfun(@(x) str2double(regexp(x,'[+-]?\d*\.?\d*','match')),reg_out);
nums_out = abs(nums_out);
%##
% x_in = 0:max(nums_out);
% dist_out = fitdist(nums_out','poisson');
% y_out = poisson(x_in,dist_out.lambda);
%##
fig = figure;
hold on;
hh = histogram(nums_out);
% plot(x_in,y_out*max(hh.Values)*3);
title(sprintf('ICanClean Rejected Components N=%i (EEG R^2=0.65, EMG R^2=0.3)',length(nums_out)));
xlabel('Number of Components Rejected');
ylabel('Counts');
ylim([0,20])
vline(mean(nums_out),'-',sprintf('Mu=%0.1f',mean(nums_out)),[-0.01,0.7]);
hold off;
%%
% fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\1_BATCH_PREP\MIM_OA\_hpg_logs';
% fname = '11904335_MIM_prep.log';
% fid = fileread([fpath filesep fname]);
% reg_out = regexp(fid,'[+-]?\d*\.?\d* channels were rejected','match');
% nums_out = cellfun(@(x) str2double(regexp(x,'[+-]?\d*\.?\d*','match')),reg_out);
% nums_out = abs(nums_out);
% %##
% out = poisson(2.842,nums_out);
% figure;
% hh = histogram(nums_out);
% title(sprintf('ICanClean Rejected Components N=%i (EEG R^2=0.65, EMG R^2=0.4)',length(nums_out)));
% xlabel('Number of Components Rejected');
% ylabel('Counts');
% ylim([0,25])
% vline(mean(nums_out),'-',sprintf('Mu=%0.1f',mean(nums_out)),[-0.01,0.7]);
%%
% fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\1_BATCH_PREP\MIM_OA\_hpg_logs';
% fname = '7397262_MIM_prep.log';
% fid = fileread([fpath filesep fname]);
% reg_out = regexp(fid,'[+-]?\d*\.?\d* channels were rejected','match');
% nums_out = cellfun(@(x) str2double(regexp(x,'[+-]?\d*\.?\d*','match')),reg_out);
% nums_out = abs(nums_out);
% %##
% figure;
% hh = histogram(nums_out);
% title(sprintf('ICanClean Rejected Components N=%i (EEG R^2=0.65, EMG R^2=0.4)',length(nums_out)));
% xlabel('Number of Components Rejected');
% ylabel('Counts');
% ylim([0,20])
% vline(mean(nums_out),'-',sprintf('Mu=%0.1f',mean(nums_out)),[-0.01,0.7]);
%%
fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\1_BATCH_PREP\MIM_OA\_hpg_logs';
fname = '13909746_MIM_prep.log';
fid = fileread([fpath filesep fname]);
reg_out = regexp(fid,'[+-]?\d*\.?\d* channels were rejected','match');
nums_out = cellfun(@(x) str2double(regexp(x,'[+-]?\d*\.?\d*','match')),reg_out);
nums_out = abs(nums_out);
%##
figure;
hh = histogram(nums_out);
title(sprintf('ICanClean Rejected Components N=%i (EEG R^2=0.65, EMG R^2=0.3)',length(nums_out)));
xlabel('Number of Components Rejected');
ylabel('Counts');
ylim([0,20])
vline(mean(nums_out),'-',sprintf('Mu=%0.1f',mean(nums_out)),[-0.01,0.7]);
%%
fpath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\1_BATCH_PREP\MIM_OA\_hpg_logs';
fname = '13909746_MIM_prep.log';
fid = fileread([fpath filesep fname]);
reg_out = regexp(fid,'[+-]?\d*\.?\d* channels were rejected','match');
nums_out = cellfun(@(x) str2double(regexp(x,'[+-]?\d*\.?\d*','match')),reg_out);
nums_out = abs(nums_out);
%##
figure;
hh = histogram(nums_out);
title(sprintf('ICanClean Rejected Components N=%i (EEG R^2=0.60, EMG R^2=0.4)',length(nums_out)));
xlabel('Number of Components Rejected');
ylabel('Counts');
ylim([0,20])
vline(mean(nums_out),'-',sprintf('Mu=%0.1f',mean(nums_out)),[-0.01,0.7]);