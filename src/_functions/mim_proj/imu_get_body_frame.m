function [NEW_EEG_IMU] = imu_get_body_frame(EEG_IMU, save_dir, varargin)
%EXTRACT_GAIT_EVENTS Summary of this function goes here
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
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%##
tic
%% DEFINE DEFAULTS
if ~which('quaternRotate')
    error('You need to add the Gait-Tracking-With-x-IMU-master free library to get quaternion math');
end
HPfilter_passBand   = 0.2; %Hi Colton. Ryan here. I decided to with 0.2 Hz instead of 0.25. Let me know how it looks
LPfilter_passBand   = 20; %Hz
% MIM_R_FOLDER        = ['R:' filesep 'Ferris-Lab' filesep 'share' filesep 'MindInMotion'];
% MIM_R_FOLDER = [filesep 'blue' filesep 'jsalminen' filesep 'GitHub' filesep 'par_EEGProcessing',...
%     filesep 'src' filesep '_data' filesep 'MIM_dataset' filesep '_studies' filesep 'subject_mgmt'];
% TRIAL_CROPPING_FNAME = 'Trial_Cropping_V2_test.xlsx';
TRIAL_CROPPING_XLSX = ['R:' filesep 'Ferris-Lab' filesep 'share' filesep 'MindInMotion' filesep 'Trial_Cropping_V2_test.xlsx'];
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'EEG_IMU',@isstruct);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
addParameter(p,'TRIAL_CROPPING_XLSX',TRIAL_CROPPING_XLSX,@ischar);
%## PARAMETER
parse(p, EEG_IMU, save_dir, varargin{:});
%## SET DEFAULTS
%- OPTIONALS

%- PARAMETER
TRIAL_CROPPING_XLSX = p.Results.TRIAL_CROPPING_XLSX;
%- PERMS
%% ===================================================================== %%
%## Cropping parameters from problematic IMU excel
fprintf('==== Cropping MIM Trials For Bad Timeframes ====\n');
tic
do_crop = cell(1,length(unique({EEG_IMU.event.trialName})));
crop_idxs = cell(1,length(unique({EEG_IMU.event.trialName})));
for trial_i = 1:length(unique({EEG_IMU.event.trialName}))
    [~,idx_end] = regexp(EEG_IMU.event(trial_i).trialName,sprintf('%s_',EEG_IMU.subject));
    trial_name = EEG_IMU.event(trial_i).trialName; trial_name = trial_name(idx_end+1:end);
    num_check = regexp(trial_name,'\d');
    if ~isempty(num_check)
        trial_num = str2double(trial_name(num_check(end)));
    else
        trial_num = 1;
    end
    if trial_num <= 2
        [ tmp_dc, tmp_ci, ~, ~ ] = mim_check_trials(EEG_IMU.subject,trial_name,TRIAL_CROPPING_XLSX); %grabbing cropping times if needed
        % do_crop{trial_i} = false;
        % crop_idxs{trial_i} = [];
        if isnan(tmp_ci)
            fprintf('%s) trial %s not found. Not including...\n',EEG_IMU.subject,trial_name);
        else
            do_crop{trial_i}=tmp_dc;
            crop_idxs{trial_i}=tmp_ci;
            if do_crop{trial_i}
                fprintf('Cropping data for subject %s and trial %s\n',EEG_IMU.subject,trial_name);
            end
        end
        
    else
        fprintf('Trial number is %i, not performing cropping for %s...\n',trial_num,trial_name);
        do_crop{trial_i} = false;
        crop_idxs{trial_i} = [];
    end
end
toc
%% Grab accelerometer channel data & APDM auto-generated quanternions
acc_chan = regexp({EEG_IMU.chanlocs.labels},'acc'); acc_chan = cellfun(@(x) ~isempty(x),acc_chan);
gyro_chan = regexp({EEG_IMU.chanlocs.labels},'gyro'); gyro_chan = cellfun(@(x) ~isempty(x),gyro_chan);
mag_chan = regexp({EEG_IMU.chanlocs.labels},'mag'); mag_chan = cellfun(@(x) ~isempty(x),mag_chan);
orient_chan = regexp({EEG_IMU.chanlocs.labels},'ori'); orient_chan = cellfun(@(x) ~isempty(x),orient_chan);

%- converting to column data FYI
sensor_frame_gyro = EEG_IMU.data(gyro_chan,:)'; % 3xtime_pnts 
sensor_frame_mag = EEG_IMU.data(mag_chan,:)'; % 3xtime_pnts
sensor_frame_acc = EEG_IMU.data(acc_chan,:)'; % 3xtime_pnts
%- Dataset containing the orientation quaternio
sensor_frame_ori = EEG_IMU.data(orient_chan,:)'; % 4xtime_pnts
%% convert acc data to world frame with quaternions
world_frame_acc = zeros(size(sensor_frame_acc));
%- go through each time point
for samp_i =1:size(sensor_frame_acc,1)
    %- grab the quaternion that tells us how our local IMU frame relates to the world frame (sensor fusion already done)
    qtemp = sensor_frame_ori(samp_i,:); 
    %- rotate the local accelerations to world frame
    world_frame_acc(samp_i,:)= quaternRotate(sensor_frame_acc(samp_i,:), qtemp); 
end
%- 
% world_frame_acc(isnan(world_frame_acc)) = deal(0);
%## remove the constant effect of gravity and/or constant accelerations in any other direction 
gConstant = mean(world_frame_acc(:,3)); %units m/s^2 %should be positive or negative????
fprintf('Gravity = %s m/s^2\n',num2str(gConstant));
%acc_worldFrame_noGravity = acc_worldFrame - [0 0 gConstant];
world_frame_accNoGrav = world_frame_acc - mean(world_frame_acc);

%Optionally filter the accelerometer %Note: it may cause some weird results to
%filter acc/vel/pos rather than just vel/pos
world_frame_accNoGravFilter = FilterIMUfunction(world_frame_accNoGrav,HPfilter_passBand,LPfilter_passBand, EEG_IMU.srate);

%## Integrate acceleration to obtain velocity
world_frame_vel = cumtrapz((world_frame_accNoGrav))*(1/EEG_IMU.srate); %m/s

% vel_est = detrend(vel_est); %2021-03-24 RD: note may not be needed to detrend given we are now BP filtering the velocity

%FILTERING velocity data colton 3/24/21
world_frame_vel = FilterIMUfunction(world_frame_vel,HPfilter_passBand,LPfilter_passBand, EEG_IMU.srate);

%Integrate velocity to obtain position
world_frame_pos = ( cumtrapz((world_frame_vel))*(1/EEG_IMU.srate) ); %m

%Filter position data
world_frame_pos = FilterIMUfunction(world_frame_pos,HPfilter_passBand,LPfilter_passBand, EEG_IMU.srate);

%% PLOTS
%-
% trial_name = EEG_IMU.event(TRIAL_IND).trialName;
TRIAL_CHAR_TO_TEST = 'low'; % flat, low, med, high, 0p25, 0p5, 0p75, 1p0
SAMPLES_TO_PLOT = 500;
%- prep
inds = cellfun(@(x) contains(x,TRIAL_CHAR_TO_TEST),{EEG_IMU.event.trialName});
[~,TRIAL_IND] = find(inds,1,'first');
trial_idxs = [EEG_IMU.event(TRIAL_IND).latency, EEG_IMU.event(TRIAL_IND+1).latency];
ind_1 = trial_idxs(1)+500;
ind_2 = trial_idxs(1)+SAMPLES_TO_PLOT+500;
time = EEG_IMU.times(ind_1:ind_2)/1000;
% idx_1 = find([EEG_IMU.event.latency] <= ind_1,1,'last');
% idx_2 = find([EEG_IMU.event.latency] >= ind_2,1,'first');
trial_name = EEG_IMU.event(TRIAL_IND).trialName;
% plot 1
gyrX = sensor_frame_gyro(ind_1:ind_2,1);
gyrY = sensor_frame_gyro(ind_1:ind_2,2);
gyrZ = sensor_frame_gyro(ind_1:ind_2,3);
s_accX = sensor_frame_acc(ind_1:ind_2,1);
s_accY = sensor_frame_acc(ind_1:ind_2,2);
s_accZ = sensor_frame_acc(ind_1:ind_2,3);
w_accX = world_frame_acc(ind_1:ind_2,1);
w_accY = world_frame_acc(ind_1:ind_2,2);
w_accZ = world_frame_acc(ind_1:ind_2,3);
acc_fX = world_frame_accNoGravFilter(ind_1:ind_2,1);
% plot 2
w_velX = world_frame_vel(ind_1:ind_2,1);
w_velY = world_frame_vel(ind_1:ind_2,2);
w_velZ = world_frame_vel(ind_1:ind_2,3);
w_posX = world_frame_pos(ind_1:ind_2,1);
w_posY = world_frame_pos(ind_1:ind_2,2);
w_posZ = world_frame_pos(ind_1:ind_2,3);
%## graph for gyroscope and accelerometer
%{
figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', sprintf('Sensor Data: %s',trial_name));
ax(1) = subplot(2,1,1);
    hold on;
    plot(time, gyrX, 'r');
    plot(time, gyrY, 'g');
    plot(time, gyrZ, 'b');
    title('Gyroscope');
    xlabel('Time (s)');
    ylabel('Angular velocity (^\circ/s)');
    legend('X', 'Y', 'Z');
    hold off;
ax(2) = subplot(2,1,2);
    hold on;
    plot(time, s_accX, 'r');
    plot(time, s_accY, 'g');
    plot(time, s_accZ, 'b');
    plot(time, w_accX, '--r');
    plot(time, w_accY, '--g');
    plot(time, w_accZ, '--b');
    plot(time, acc_fX, ':r', 'LineWidth', 2);
    title('Accelerometer');
    xlabel('Time (s)');
    ylabel('Acceleration (g)');
    legend('sensor_X', 'sensor_Y', 'sensor_Z', 'world_X', 'world_Y', 'world_Z', 'filtered world_X');
    hold off;
linkaxes(ax,'x');
fig_i = get(groot,'CurrentFigure');
% saveas(fig_i,[save_dir filesep sprintf('accValidation_%s.fig',trial_name)]);
saveas(fig_i,[save_dir filesep sprintf('accValidation_%s.jpg',trial_name)]);
%-
figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', sprintf('Sensor Estimates: %s',trial_name));
ax(1) = subplot(2,1,1);
    hold on;
    plot(time, w_velX, 'r');
    plot(time, w_velY, 'g');
    plot(time, w_velZ, 'b');
    title('Velocity Estimate');
    xlabel('Time (s)');
    ylabel('Velocity  (m/s)');
    legend('X', 'Y', 'Z');
    hold off;
ax(2) = subplot(2,1,2);
    hold on;
    plot(time, w_posX, 'r');
    plot(time, w_posY, 'g');
    plot(time, w_posZ, 'b');
    title('Position Estimate');
    xlabel('Time (s)');
    ylabel('Position (m)');
    legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
    hold off;
linkaxes(ax,'x');
fig_i = get(groot,'CurrentFigure');
% saveas(fig_i,[save_dir filesep sprintf('velposValidation_%s.fig',trial_name)]);
saveas(fig_i,[save_dir filesep sprintf('velposValidation_%s.jpg',trial_name)]);
%}
%% find avg  body frame directions (front/left/up) in world frame and then rotate all data to avg body frame of reference
%smartest dummy method
imu_to_world_trans = mean(sensor_frame_ori); %goes from IMU to world from
% qtemp_World2IMU = quaternConj(qtemp_IMU2World); % world frame to sensor frame (not sensor to world)

samplePoint_IMUforward_inNWU    = quaternRotate([0 0 -1],imu_to_world_trans);
samplePoint_IMUleft_inNWU       = quaternRotate([0 -1 0],imu_to_world_trans);
tempRotMat                      = [samplePoint_IMUforward_inNWU(1:2)/sqrt(sum(samplePoint_IMUforward_inNWU(1:2).^2)); samplePoint_IMUleft_inNWU(1:2)./sqrt(sum(samplePoint_IMUleft_inNWU(1:2).^2))]; %2x2 %pca rotation matrix. score = centeredData*coeff; centeredData = SCORE*COEFF'
tempRotMat(3,3)                 = 1; %make it a full 3-d roation matrix

%note: tempRotMat now takes us from bodys frame to world frame thus
%inverse should give us rotation necessary to go from world frame to
%body frame (keeping front and left directions fixed to the transverse
%plane)
tempRotMat = pinv(tempRotMat); %ryan update in future 2021-03-16

% zRotAngle = acos(tempRotMat(1,1)); %analytic solution
% tempQ = rotMat2quatern( euler2rotMat(0,0,-zRotAngle) ); %negative to flip reference frame?  %Needs Gait-Tracking-.... library and subfolders to be added to path

%rotate things from global frame to something more like A/P and mediolateral
pos_est_VAF = world_frame_pos*tempRotMat; %Front, Left, Up
vel_est_VAF = world_frame_vel*tempRotMat; %Front, Left, Up
acc_est_VAF = world_frame_accNoGrav*tempRotMat; %Front, Left, Up

%% make static image of entire timecourse of trial from 3 orthogonal angles
%{
cnt = 1;
for trial_i = 1:length(unique({EEG_IMU.event.trialName}))
    trial_idxs = [EEG_IMU.event(cnt).latency, EEG_IMU.event(cnt+1).latency];
    figure('Name',sprintf('Subject-%s_Trial-%s',EEG_IMU.subject,EEG_IMU.event(cnt).trialName));
    hold on;
    if do_crop{trial_i}
        ExactCropLatency = round(EEG_IMU.srate*crop_idxs{trial_i}+1);
        nPairsTemp = size(crop_idxs{trial_i});
        nPairsReal = nPairsTemp(1); %roundabout way but length does work for 2x2 matrix as well
        for pair_i = 1:1:nPairsReal
            UsefulCropLatencies = ExactCropLatency(pair_i,1):ExactCropLatency(pair_i,2);
            subplot(1,3,1); plot(-pos_est_VAF(UsefulCropLatencies,2),pos_est_VAF(UsefulCropLatencies,1),':'); if pair_i < nPairsReal; hold on; end; xlabel('Right'); ylabel('Foward'); title('Transverse'); %top down 
            subplot(1,3,2); plot(pos_est_VAF(UsefulCropLatencies,1),pos_est_VAF(UsefulCropLatencies,3),':'); if pair_i < nPairsReal; hold on; end; xlabel('Forward'); ylabel('Up'); title('Sagittal'); %side view
            subplot(1,3,3); plot(-pos_est_VAF(UsefulCropLatencies,2),pos_est_VAF(UsefulCropLatencies,3),':'); if pair_i < nPairsReal; hold on; end; xlabel('Right'); ylabel('Up'); title('Coronal'); %from behind or in front
        end
    else 
        subplot(1,3,1); plot(-pos_est_VAF(trial_idxs(1):trial_idxs(2),2),pos_est_VAF(trial_idxs(1):trial_idxs(2),1),':'); xlabel('Right'); ylabel('Foward'); title('Transverse'); %top down 
        subplot(1,3,2); plot(pos_est_VAF(trial_idxs(1):trial_idxs(2),1),pos_est_VAF(trial_idxs(1):trial_idxs(2),3),':'); xlabel('Forward'); ylabel('Up'); title('Sagittal'); %side view
        subplot(1,3,3); plot(-pos_est_VAF(trial_idxs(1):trial_idxs(2),2),pos_est_VAF(trial_idxs(1):trial_idxs(2),3),':'); xlabel('Right'); ylabel('Up'); title('Coronal'); %from behind or in front
    end
    hold off;
    fig_i = get(groot,'CurrentFigure');
    % saveas(fig_i,[save_dir filesep sprintf('stateSpaceValid_%s.fig',EEG_IMU.event(cnt).trialName)]);
    saveas(fig_i,[save_dir filesep sprintf('stateSpaceValid_%s.jpg',EEG_IMU.event(cnt).trialName)]);
    cnt = cnt + 2;
end
close all;
%}
%% IMU 6DOF MOVIE (takes a lot of time to render)
%{
TRIAL_CHAR_TO_TEST = 'low'; % flat, low, med, high, 0p25, 0p5, 0p75, 1p0
% SAMPLES_TO_PLOT = 500;
%- prep
inds = cellfun(@(x) contains(x,TRIAL_CHAR_TO_TEST),{EEG_IMU.event.trialName});
[~,TRIAL_IND] = find(inds,1,'first');
trial_idxs = [EEG_IMU.event(TRIAL_IND).latency, EEG_IMU.event(TRIAL_IND+1).latency];
trial_idxs = [trial_idxs(1)+500,trial_idxs(1)+SAMPLES_TO_PLOT+500];
trial_name = EEG_IMU.event(TRIAL_IND).trialName;
pos_in = pos_est_VAF(trial_idxs(1):trial_idxs(2),:);
rot_in = quatern2rotMat(quaternConj(sensor_frame_ori(trial_idxs(1):trial_idxs(2),:)));
disp_imu_movie(pos_in,rot_in,[save_dir filesep sprintf('6DOF_movie_%s',trial_name)]);
%}
%% make new EEG set with A/P M/L vert positions
% chans_keep = {'gyro','mag','barom','acc'};
NEW_EEG_IMU = EEG_IMU;
n_chans = size(EEG_IMU.data,1);
NEW_EEG_IMU.pnts = size(NEW_EEG_IMU.data,2);
NEW_EEG_IMU.etc.imu_to_world_trans = imu_to_world_trans;
NEW_EEG_IMU.etc.world_to_ap_trans = tempRotMat;
%no fancy pca. use estimate of average Front, Left, Up frame
data = cat(1,EEG_IMU.data,pos_est_VAF',vel_est_VAF',acc_est_VAF',world_frame_pos',(sensor_frame_ori)');
NEW_EEG_IMU.data = data;
NEW_EEG_IMU.chanlocs(n_chans+1).labels  = 'body_xPos';
NEW_EEG_IMU.chanlocs(n_chans+2).labels  = 'body_yPos';
NEW_EEG_IMU.chanlocs(n_chans+3).labels  = 'body_zPos';
NEW_EEG_IMU.chanlocs(n_chans+4).labels  = 'body_xVel';
NEW_EEG_IMU.chanlocs(n_chans+5).labels  = 'body_yVel';
NEW_EEG_IMU.chanlocs(n_chans+6).labels  = 'body_zVel';
NEW_EEG_IMU.chanlocs(n_chans+7).labels  = 'body_xAcc';
NEW_EEG_IMU.chanlocs(n_chans+8).labels  = 'body_yAcc';
NEW_EEG_IMU.chanlocs(n_chans+9).labels  = 'body_zAcc';
NEW_EEG_IMU.chanlocs(n_chans+10).labels  = 'world_xPos';
NEW_EEG_IMU.chanlocs(n_chans+11).labels  = 'world_yPos';
NEW_EEG_IMU.chanlocs(n_chans+12).labels  = 'world_zPos';
NEW_EEG_IMU.chanlocs(n_chans+13).labels = 'q1';
NEW_EEG_IMU.chanlocs(n_chans+14).labels = 'qi';
NEW_EEG_IMU.chanlocs(n_chans+15).labels = 'qj';
NEW_EEG_IMU.chanlocs(n_chans+16).labels = 'qk';
for i = 1:length(NEW_EEG_IMU.chanlocs)
    if contains(NEW_EEG_IMU.chanlocs(i).labels,'Back')
        out = strsplit(NEW_EEG_IMU.chanlocs(i).labels,'_');
        out = strjoin(['orig',out(2:end)],'_');
        NEW_EEG_IMU.chanlocs(i).labels = out;
    end
    NEW_EEG_IMU.chanlocs(i).type = 'IMU_Opal';
end
%##
toc
end

%% SUBFUNCTIONS
function [ xfilt ] = FilterIMUfunction(x,HPfilter_passBand,LPfilter_passBand, srate)
%This functino expects column vectors
%It applies 4th order butterworth filtering (forward and reverse with filtfilt)
%You can do BP filtering with this function. 
%You can also do HP only if you set the LP input to inf.
%This function does some special version of padding prior to filtering to prevent edge effects 

%## filter position data
%note: we are doing fancy copying of fake data to the left and right of
%original dataset so we can filter the data without edge effects. Our
%original IMU data is sharply cut into 180-sec long trials. We want to use
%all of that data, but we also need to filter. Since filters use windows,
%we need extra padded data.

%temporarily add extra data to edges to help avoid edge effects
%borrowing code found online for minimizing edge effects
R=0.1; % 10% of signal
Nr=5000;
N=size(x,1);
NR=min(round(N*R),Nr); % At most xxx points
for i=1:size(x,2)
    x1(:,i)=2*x(1,i)-flipud(x(2:NR+1,i));  % maintain continuity in level and slope
    x2(:,i)=2*x(end,i)-flipud(x(end-NR:end-1,i));
end
x=[x1;x;x2]; %append extra data to each size of the original (x)

% Do filtering
disp('We are doing Butterworth bandpass filtering');
disp(['[',num2str([HPfilter_passBand, LPfilter_passBand]),'] Hz']);
if isinf(LPfilter_passBand) %highpass only
    dig_filt = designfilt('highpassiir', 'FilterOrder', 4,...
                            'HalfPowerFrequency', HPfilter_passBand, 'SampleRate',srate); %hp only
%     x = filtfilt(dig_filt.Coefficients(1,1:3),dig_filt.Coefficients(2,1:3),x);
    x = filter(dig_filt,x); % JS (02/09/2023), changed to using filter as recommended by the designfilt function
else %if HP and LP cutoffs provided (non-inf)
    dig_filt = designfilt('bandpassiir', 'FilterOrder', 4,...
                            'HalfPowerFrequency1', HPfilter_passBand, 'HalfPowerFrequency2', LPfilter_passBand, 'SampleRate', srate);%bandpass
%     x = filtfilt(dig_filt.Coefficients(1,1:3),dig_filt.Coefficients(2,1:3),x);
    x = filter(dig_filt,x); % JS (02/09/2023), changed to using filter as recommended by the designfilt function
end

%crop back to original data (temporarily added extra to the front and back)
x=x(NR+1:end-NR,:); %go back to original time window

%define output filtered data
xfilt = x; 
end
