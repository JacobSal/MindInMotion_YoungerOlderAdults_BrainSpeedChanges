function [EEG_LS] = create_loadsol_struct(loadsol_table,rec_start_struct,save_dir,varargin)
%CREATE_LOADSOL_EVENTS Summary of this function goes here
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
% Code Date: 01/10/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu

%% Define Defaults
GRF_THRESHOLD = 20; %units = N 
BAD_STEPS_THRESHOLD = 850;
BAD_STEPS_SAMPS_DELETE = 200;
SYNC_PULSE_WIDTH = 0.05; 
SYNC_THRESHOLD = 6000; %default 6000 N/s. Note @200Hz a thres of 6000 N/s is equivalent to a 30 N jump from one sample point to the next %2021-04-29 RJD temp changed to 2500 for subject NH3002
                   %used 2000N for H1046
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'loadsol_table',@istable);
addRequired(p,'rec_start_struct',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p, loadsol_table, rec_start_struct, save_dir, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS
%% ===================================================================== %%
%## Extract data from loadsol table
datetime_data = loadsol_table.datetime;
left_foot_data = loadsol_table.left_foot_N;
right_foot_data = loadsol_table.right_foot_N;
sync_foot_data = loadsol_table.sync_foot_N;
sample_rate = 1/seconds(loadsol_table.time_seconds(2)-loadsol_table.time_seconds(1));
%## Find sync events
%- edit this function if you want to temporarily change the threshold for 
%  finding sync events
[sync_events_t, sync_events_samp] = findLoadSolSyncEvents(sync_foot_data,sample_rate,save_dir,SYNC_THRESHOLD,SYNC_PULSE_WIDTH); 
sync_datetimes = loadsol_table.datetime(sync_events_samp);
medianDiffSync_t = median(diff(sync_events_t));
disp(['Sync events are ',num2str(medianDiffSync_t),' seconds apart on average.']);
%- concatenate loadsol data
concatData =  [left_foot_data; right_foot_data; sync_foot_data];
%- create event structure for loadsol
EEG_LS = pop_importdata('dataformat','array','nbchan',size(concatData,1),'data',concatData,'setname','LoadSol','srate',sample_rate,'pnts',size(concatData,2),'xmin',0);
%- define labels/types for channels
EEG_LS.chanlocs(1).labels = 'GRF_L';
EEG_LS.chanlocs(2).labels = 'GRF_R';
EEG_LS.chanlocs(3).labels = 'LoadSolSync';
EEG_LS.chanlocs(4).labels = 'datetime';

EEG_LS.chanlocs(1).type = 'BioM';
EEG_LS.chanlocs(2).type = 'BioM';
EEG_LS.chanlocs(3).type = 'Sync';
EEG_LS.chanlocs(4).type = 'datetime';

%- Add recording time structure to EEG structure
EEG_LS.recordingStart = loadsol_table.datetime;

for event_i = 1:length(sync_events_samp)
    EEG_LS.event(event_i).type = 'Sync';
    EEG_LS.event(event_i).latency = sync_events_samp(event_i);
    EEG_LS.event(event_i).datetime = sync_datetimes(event_i);
end
%% GAIT CYCLE EVENT FINDER
fprintf('==== Finding Gait Events from Loadsol ====\n');
leftInd = strcmpi({EEG_LS.chanlocs.labels},'GRF_L'); %strcmpi for case insensitive
rightInd = strcmpi({EEG_LS.chanlocs.labels},'GRF_R');

LeftFP_Fz = EEG_LS.data(leftInd,:);
RightFP_Fz = EEG_LS.data(rightInd,:);

%- find points where ground reaction force indicates foot contact
left_logical=LeftFP_Fz > GRF_THRESHOLD; %makes logical of contact (1 or 0) based on exceeding threshold
right_logical=RightFP_Fz > GRF_THRESHOLD;
%- extract events using logicals
left_heel_strike=find(diff(left_logical)>0); %finds the indices when GRF exceeds threshold
right_heel_strike=find(diff(right_logical)>0);
left_toe_off=find(diff(left_logical)<0); %finds the indices when GRF exceeds threshold
right_toe_off=find(diff(right_logical)<0);
%- recording grf at event
left_grf_hs=LeftFP_Fz((diff(left_logical)>0));
right_grf_hs=RightFP_Fz(diff(right_logical)>0);
left_grf_to=LeftFP_Fz(diff(left_logical)<0);
right_grf_to=RightFP_Fz(diff(right_logical)<0);

%## PLOT
cMaps = linspecer(4);
figure('Name',sprintf('Subject-%s_Loadsol',EEG_LS.subject));
hold on;
%- raw signals
plot(EEG_LS.times/1000,LeftFP_Fz,'-','Color',[252 174 145]/255,'DisplayName','Left Foot')
plot(EEG_LS.times/1000,RightFP_Fz,'-','Color',[189 215 231]/255,'DisplayName','Right Foot')
%- extracted events
stem(EEG_LS.times(left_heel_strike)/1000,20*ones(1,length(left_heel_strike)),'Color',cMaps(1,:),'DisplayName','LHS')
stem(EEG_LS.times(right_heel_strike)/1000,20*ones(1,length(right_heel_strike)),'Color',cMaps(2,:),'DisplayName','RHS')
stem(EEG_LS.times(left_toe_off)/1000,20*ones(1,length(left_toe_off)),'Color',cMaps(3,:),'DisplayName','LTO')
stem(EEG_LS.times(right_toe_off)/1000,20*ones(1,length(right_toe_off)),'Color',cMaps(4,:),'DisplayName','RTO')
xlabel('Time (s)')
ylabel('Ground Reaction Force (N)')
hold off;
legend();
drawnow;
saveas(gcf,fullfile(save_dir, 'validation_figure_ls.jpg'))
%% Reject Bad Foot Strikes
diff_rightHS=diff(right_heel_strike);
diff_leftHS=diff(left_heel_strike);
% std_heelstrike=std([diff_rightHS'; diff_leftHS']);
mean_heelstrike=trimmean([diff_rightHS'; diff_leftHS'],10);

diff_rightTO=diff(right_toe_off);
diff_leftTO=diff(left_toe_off);
% std_toeoff=std([diff_rightTO'; diff_leftTO']);
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

%%
deletesteps_rightHS=find(diff_rightHS(badsteps_rightHS)-BAD_STEPS_SAMPS_DELETE<0);
deletesteps_leftHS=find(diff_leftHS(badsteps_leftHS)-BAD_STEPS_SAMPS_DELETE<0);
deletesteps_rightTO=find(diff_rightTO(badsteps_rightTO)-BAD_STEPS_SAMPS_DELETE<0);
deletesteps_leftTO=find(diff_leftTO(badsteps_leftTO)-BAD_STEPS_SAMPS_DELETE<0);

right_heel_strike([1 badsteps_rightHS((deletesteps_rightHS))+1])=NaN;
left_heel_strike([1 badsteps_leftHS((deletesteps_leftHS))+1])=NaN;
right_toe_off([1 badsteps_rightTO((deletesteps_rightTO))+1])=NaN;
left_toe_off([1 badsteps_leftTO((deletesteps_leftTO))+1])=NaN;

LHS=left_heel_strike(~isnan(left_heel_strike));
RHS=right_heel_strike(~isnan(right_heel_strike));
LTO=left_toe_off(~isnan(left_toe_off));
RTO=right_toe_off(~isnan(right_toe_off));
GRF_LHS = left_grf_hs(~isnan(left_heel_strike));
GRF_RHS = right_grf_hs(~isnan(right_heel_strike));
GRF_LTO = left_grf_to(~isnan(left_toe_off));
GRF_RTO = right_grf_to(~isnan(right_toe_off));
badsteps=[length(badsteps_rightHS) length(badsteps_leftHS) length(badsteps_rightTO) length(badsteps_leftTO)];
badstepsname=['badsteps_rightHS, ', 'badsteps_leftHS, ', 'badsteps_rightTO, ', 'badsteps_leftTO '];

%% update EEG structure with gait events
eventStrings = {'LHS','LTO','RHS','RTO'};
eventSamples = {LHS, LTO, RHS, RTO};
eventGRF = {GRF_LHS, GRF_LTO, GRF_RHS, GRF_RTO};
for eventType_i = 1:length(eventStrings)
    eventStr = eventStrings{eventType_i};
    tempEvents_sample = eventSamples{eventType_i};
    tempGrf = eventGRF{eventType_i};
    for event_j = 1:length(tempEvents_sample)
        EEG_LS.event(end+1).latency = tempEvents_sample(event_j);
        EEG_LS.event(end).duration = 1;
        EEG_LS.event(end).duration = 1; %note that we use 'end' here instead of end+1 since the length of EEG.event has been updated
        EEG_LS.event(end).channel = 0;
        EEG_LS.event(end).type = eventStr;
        EEG_LS.event(end).code = eventStr;
        EEG_LS.event(end).datetime = loadsol_table.datetime(tempEvents_sample(event_j));
        EEG_LS.event(end).grf = tempGrf(event_j);
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

%% SUBFUNCTIONS
function [ syncEvents_t, syncEvents_samp ] = findLoadSolSyncEvents( sync_foot_data, sample_rate, save_dir, SYNC_THRESHOLD, SYNC_PULSE_WIDTH )
%FINDLOADSOLSYNCEVENTS Summary of this function goes here
%   Detailed explanation goes here
%## NOTES:
%       2021-04-29 RJD fixed error in equation for thresholding sync signal
%       (multiply vs divide). This didn't have any negative consequences since it
%       would only kick in once we started analyzing data that was recorded at a
%       sample frequency other than 200Hz which the old threshold was originally
%       calculated for. 

deltaT = 1/sample_rate; %sample period in seconds
                   
disp(['Using threshold of ',num2str(SYNC_THRESHOLD),' Newtons per second to find rising edges of sync signal'])

event_ind = find((diff(sync_foot_data)/deltaT)>=SYNC_THRESHOLD ); %2021-04-29 fixed typo in equation (see other note above for more info)
% event_ind = find(diff(analogSyncSig/Fs)>thresh); 
event_ind = event_ind + 1; %so rising edge happens at timestamp where signal actually starts to go up rather than right before it (to match how normal hardware would mark an event)

sampleWindow = round(SYNC_PULSE_WIDTH*sample_rate); %Default was originally round(0.025 seconds * Fs) which is 5 samples at 200 Hz %2021-04-29 RJD changed to 0.05*Fs (roughly 10 samples) primarily for subj NH3002 but it could be good for others too
newEvent_ind = event_ind(1); %grab first
%remove extra
for event_i = 2:length(event_ind)
    if event_ind(event_i)-newEvent_ind(end) > sampleWindow
        newEvent_ind(end+1) = event_ind(event_i);
    end
end
numRemoved = length(event_ind)-length(newEvent_ind);
disp(['Removed ',num2str(numRemoved),' extra events after simple thresholding']);
event_ind = newEvent_ind;
disp(['Found ',num2str(length(event_ind)),' sync events on loadsol data']);

%## Create figure to visualize loadsol events
myfig = figure; hold on;
stem(diff(event_ind)/sample_rate,'Color',[252 174 145]/255);
plot([1,length(event_ind)-1], median(diff(event_ind/sample_rate))*ones(1,2));
title('Relative time between events');
xlabel('Event num');
ylabel('Duration (s)');
hold off;
drawnow;
saveas(myfig,fullfile(save_dir, 'timeBetweenSyncEvents.fig'));
disp(['Median diff between events = ',num2str(median(diff(event_ind))/sample_rate)]);

%## Create figure to visualize loadsol events
time = 0:1/sample_rate:length(sync_foot_data)/sample_rate-1/sample_rate;
myfig = figure; hold on;
plot(time, sync_foot_data);
stem(time(event_ind),max(sync_foot_data)*ones(size(event_ind)),'Color',[252 174 145 100]/255);
title('Events plotted on top of analog sync signal');
xlabel('Local time (s)');
ylabel('Fake foot force (N)');
hold off;
drawnow;
saveas(myfig,fullfile(save_dir, 'syncEventsMarked.fig'));

syncEvents_samp = event_ind;
syncEvents_t = time(syncEvents_samp);
end

