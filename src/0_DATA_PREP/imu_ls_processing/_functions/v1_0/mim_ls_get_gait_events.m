function [EEG_LS] = mim_ls_get_gait_events(EEG_LS,save_dir,varargin)
%MIM_LS_GET_GAIT_EVENTS Summary of this function goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/13/2023, MATLAB 2019a
% Copyright (C) Ryan Downey
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%% Define Defaults
GRF_THRESHOLD = 20; %units = N 
BAD_STEPS_THRESHOLD = 850;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG_LS');
addRequired(p,'save_dir');
%## OPTIONAL
%## PARAMETER
parse(p, EEG_LS, save_dir, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS
%% ===================================================================== %%
%## finding gait events
fprintf('==== Finding Gait Events from Loadsol ====\n');
leftInd = strcmpi({EEG_LS.chanlocs.labels},'GRF_L'); %strcmpi for case insensitive
rightInd = strcmpi({EEG_LS.chanlocs.labels},'GRF_R');

LeftFP_Fz = EEG_LS.data(leftInd,:);
RightFP_Fz = EEG_LS.data(rightInd,:);

left_logical=LeftFP_Fz > GRF_THRESHOLD; %makes logical of contact (1 or 0) based on exceeding threshold
right_logical=RightFP_Fz > GRF_THRESHOLD;
left_heel_strike=find(diff(left_logical)>0); %finds the indices when GRF exceeds threshold
right_heel_strike=find(diff(right_logical)>0);
left_toe_off=find(diff(left_logical)<0); %finds the indices when GRF exceeds threshold
right_toe_off=find(diff(right_logical)<0);

% myHSfig = figure; plot(abs(RightFP_Fz)); hold on; plot(right_heel_strike,abs(RightFP_Fz(right_heel_strike)),'r^');
% plot(abs(LeftFP_Fz),'k'); plot(left_heel_strike,abs(LeftFP_Fz(left_heel_strike)),'ro');
% plot(right_toe_off,abs(RightFP_Fz(right_toe_off)),'g^');
% plot(left_toe_off,abs(LeftFP_Fz(left_toe_off)),'go');
% 
% drawnow;
% %saveas(gcf,['plots/' csv_filename(firstletter:lastletter) '_FZ_FP_fixed'])
%%
diff_rightHS=diff(right_heel_strike);
diff_leftHS=diff(left_heel_strike);
std_heelstrike=std([diff_rightHS'; diff_leftHS']);
mean_heelstrike=trimmean([diff_rightHS'; diff_leftHS'],10);

diff_rightTO=diff(right_toe_off);
diff_leftTO=diff(left_toe_off);
std_toeoff=std([diff_rightTO'; diff_leftTO']);
mean_toeoff=trimmean([diff_rightTO'; diff_leftTO'],10);

%ryan note: thresholds below are in terms of samples of time. May want to
%adjust scale or change to percentage based. Originally everything below
%was > 150
badsteps_rightHS=find(abs(diff_rightHS-mean_heelstrike) > BAD_STEPS_THRESHOLD);  %if time point is 100 above mean, mark as bad
badsteps_leftHS=find(abs(diff_leftHS-mean_heelstrike) > BAD_STEPS_THRESHOLD);  %if time point is 100 above mean, mark as bad
badsteps_rightTO=find(abs(diff_rightTO-mean_toeoff) > BAD_STEPS_THRESHOLD);  %if time point is 100 above mean, mark as bad
badsteps_leftTO=find(abs(diff_leftTO-mean_toeoff) > BAD_STEPS_THRESHOLD);  %if time point is 100 above mean, mark as bad

disp(['bad right heel strikes: ' num2str((length(badsteps_rightHS)/length(right_heel_strike))*100) '% ('  num2str(length(badsteps_rightHS)) '/'  num2str(length(right_heel_strike)) ')' ]);
disp(['bad left heel strikes: '  num2str((length(badsteps_leftHS)/length(left_heel_strike))*100) '% ('  num2str(length(badsteps_leftHS)) '/'  num2str(length(left_heel_strike)) ')' ]);
disp(['bad right toe offs: ' num2str((length(badsteps_rightTO)/length(right_toe_off))*100) '% ('  num2str(length(badsteps_rightTO)) '/'  num2str(length(right_toe_off)) ')' ]);
disp(['bad left toe offs: '  num2str((length(badsteps_leftTO)/length(left_toe_off))*100) '% ('  num2str(length(badsteps_leftTO)) '/'  num2str(length(left_toe_off)) ')' ]);

% figure(5); title('difference in heel strike time'); hold on; plot(diff_rightHS); plot(diff_leftHS,'r');
% plot(1:length(diff_rightHS),mean_heelstrike+std_heelstrike*2);
% plot(1:length(diff_rightHS),mean_heelstrike-std_heelstrike*2);
% plot(badsteps_rightHS,diff_rightHS(badsteps_rightHS),'k^');
% plot(badsteps_leftHS,diff_leftHS(badsteps_leftHS),'k^');
% set(gcf,'Position',[100 525 560 420])
% %saveas(gcf,['plots/' csv_filename(firstletter:lastletter) '_heelstrike_fixed'])
%
% figure(6); title('difference in toe off time'); hold on; plot(diff_rightTO); plot(diff_leftTO,'r');
% plot(1:length(diff_rightTO),mean_toeoff+std_toeoff*2);
% plot(1:length(diff_rightTO),mean_toeoff-std_toeoff*2);
% plot(badsteps_rightTO,diff_rightTO(badsteps_rightTO),'k^');
% plot(badsteps_leftTO,diff_leftTO(badsteps_leftTO),'k^');
% set(gcf,'Position',[900 525 560 420])
% % %saveas(gcf,['plots/' csv_filename(firstletter:lastletter) '_toeoff_fixed'])
%%
deletesteps_rightHS=find(diff_rightHS(badsteps_rightHS)-200<0);
deletesteps_leftHS=find(diff_leftHS(badsteps_leftHS)-200<0);
deletesteps_rightTO=find(diff_rightTO(badsteps_rightTO)-200<0);
deletesteps_leftTO=find(diff_leftTO(badsteps_leftTO)-200<0);

right_heel_strike([1 badsteps_rightHS((deletesteps_rightHS))+1])=NaN;
left_heel_strike([1 badsteps_leftHS((deletesteps_leftHS))+1])=NaN;
right_toe_off([1 badsteps_rightTO((deletesteps_rightTO))+1])=NaN;
left_toe_off([1 badsteps_leftTO((deletesteps_leftTO))+1])=NaN;

LHS=left_heel_strike(~isnan(left_heel_strike));
RHS=right_heel_strike(~isnan(right_heel_strike));
LTO=left_toe_off(~isnan(left_toe_off));
RTO=right_toe_off(~isnan(right_toe_off));
badsteps=[length(badsteps_rightHS) length(badsteps_leftHS) length(badsteps_rightTO) length(badsteps_leftTO)];
badstepsname=['badsteps_rightHS, ', 'badsteps_leftHS, ', 'badsteps_rightTO, ', 'badsteps_leftTO '];

%% update EEG structure with gait events
eventStrings = {'LHS','LTO','RHS','RTO'};
eventSamples = {LHS, LTO, RHS, RTO};

for eventType_i = 1:length(eventStrings)
    eventStr = eventStrings{eventType_i};
    tempEvents_sample = eventSamples{eventType_i};
    for event_j = 1:length(tempEvents_sample)
        EEG_LS.event(end+1).latency = tempEvents_sample(event_j);
        EEG_LS.event(end).duration = 1;
        EEG_LS.event(end).duration = 1; %note that we use 'end' here instead of end+1 since the length of EEG.event has been updated
        EEG_LS.event(end).channel = 0;
        EEG_LS.event(end).type = eventStr;
        EEG_LS.event(end).code = eventStr;
        EEG_LS.event(end).datetime = EEG_LS.data(4,:);
    end
end

%let EEG rearrange the event order
EEG_LS = eeg_checkset(EEG_LS,'eventconsistency');

%make urevent structure
EEG_LS.urevent = EEG_LS.event;

%go back to event structure and add in corresponding urevent numbers
for event_i = 1:length(EEG_LS.event)
    EEG_LS.event(event_i).urevent = event_i;
end

end

