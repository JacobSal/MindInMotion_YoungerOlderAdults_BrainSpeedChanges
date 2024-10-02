function [EEG_IMU] = create_imu_struct(imu_table,rec_start_struct,sensor_names,varargin)
%CREATE_IMU_STRUCT Summary of this function goes here
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

%% DEFINE DEFAULTS
%- Estimate Sample Rate
% convert from microseconds unix time to time since start in seconds
IMU_SAMPLE_RATE = round(1/median(diff((imu_table{1}.time-imu_table{1}.time(1))/(10^6))));
SIGS_PER_IMU = 14; % 1 baro, 3 acc (x/y/z), 3 mag (x/y/z), 3 gyro (x/y/z), 4 orientations (s/x/y/z), 
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'imu_table',@iscell);
addRequired(p,'rec_start_struct',@isstruct);
addRequired(p,'sensor_names',@iscell);
%## OPTIONAL
%## PARAMETER
parse(p, imu_table, rec_start_struct, sensor_names, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS

%% ===================================================================== %%
%## use all data except time and remove headers
chan_names = imu_table{1}.Properties.VariableNames;
%- extract trial indices and concatentate data
imu_idxs = zeros(2,length(imu_table));
cnt = 1;
imu_data = [];
imu_idxs(1,1) = 1;
imu_idxs(2,1) = size(imu_table{1},1);
sensor_start_chans = zeros(1,length(chan_names));
for trial_i = 1:length(imu_table)
    time_chan = regexp(chan_names,'time');
    time_chan = cellfun(@(x) ~isempty(x),time_chan);
    for s_i = 1:length(sensor_names)
        sensor_name = sensor_names{s_i};
        tmp = regexp(chan_names,sensor_name);
        tmp = cellfun(@(x) ~isempty(x),tmp);
        sensor_start_chans = sensor_start_chans + tmp;
    end
    if any(sensor_start_chans)
        sensor_start_chans = logical(sensor_start_chans + time_chan);
        tmp = imu_table{trial_i};
        if cnt == 1
            imu_data = tmp(:,sensor_start_chans);
            imu_idxs(1,trial_i) = 1;
            imu_idxs(2,trial_i) = size(imu_table{1},1);
        else
            imu_data = cat(1,imu_data,tmp(:,sensor_start_chans));
            imu_idxs(1,trial_i) = imu_idxs(2,trial_i-1)+1;
            imu_idxs(2,trial_i) = imu_idxs(2,trial_i-1)+size(tmp,1);
        end
        cnt = cnt + 1;
    else
        warning('Sensor position ''%s'' not found',sensor_name);
    end
end
%- extract data from table && transpose
imu_data = imu_data{1:end,:}';
%- number of channels 
n_sigs = size(imu_data,1);
%- total length of data
n_pnts = size(imu_data,2);
%- create IMU event struct
EEG_IMU = pop_importdata('dataformat','array','nbchan',n_sigs,'data',imu_data,'setname','IMU','srate',IMU_SAMPLE_RATE,'pnts',n_pnts,'xmin',0);

%- update chanlocs info which channel names and type
for ch_i = 1:n_sigs
    EEG_IMU.chanlocs(ch_i).labels = chan_names{ch_i};
    EEG_IMU.chanlocs(ch_i).type = 'BioM';
    EEG_IMU.urchanlocs(ch_i).labels = chan_names{ch_i};
    EEG_IMU.urchanlocs(ch_i).type = 'BioM';
end

cnt = 1;
for trial_i = 1:length(imu_table)
    %- add events at trial start and trial end (to help synchronize events
    %  with EEG stream later)
    EEG_IMU.event(cnt).type	= 'TrialStart';
    EEG_IMU.event(cnt).latency = imu_idxs(1,trial_i); 
    EEG_IMU.event(cnt).trialName = rec_start_struct(trial_i).trialName;
    %- add recording time structure to EEG.event structure
    EEG_IMU.event(cnt).datetime = rec_start_struct(trial_i).datetime;

    EEG_IMU.event(cnt+1).type = 'TrialEnd';
    EEG_IMU.event(cnt+1).latency = imu_idxs(2,trial_i);
    EEG_IMU.event(cnt+1).trialName = rec_start_struct(trial_i).trialName;
    %- add recording time structure to EEG.event structure
    EEG_IMU.event(cnt+1).datetime = rec_start_struct(trial_i).datetime+seconds((imu_idxs(2,trial_i)-imu_idxs(1,trial_i))/IMU_SAMPLE_RATE);
    cnt = cnt + 2;
end
EEG_IMU.subject = rec_start_struct(1).subjectName;
EEG_IMU.setname = sprintf('%s_allTrials_IMU',rec_start_struct(1).subjectName); %e.g. H1012_high_1
end

