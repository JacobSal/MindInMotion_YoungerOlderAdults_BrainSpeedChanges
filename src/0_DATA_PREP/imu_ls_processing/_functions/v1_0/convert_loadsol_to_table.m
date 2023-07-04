function [loadsol_table,rec_start_struct] = convert_loadsol_to_table(txt_file_path,varargin)
%CONVERT_LOADSOL_TO_MAT Summary of this function goes here
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
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'txt_file_path',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p, txt_file_path, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%% ===================================================================== %%
files_dir = dir([txt_file_path filesep '*.txt']);
files_name = cell(length(files_dir),1);
rec_start_struct = struct('yearStr',[],...
                        'monthStr',[],...
                        'dayStr',[],...
                        'hourStr',[],...
                        'minStr',[],...
                        'secStr',[],...
                        'fracSecStr',[],...
                        'datetime',[]);
loadsol_table = cell(1,length(files_dir));
for f_i = 1:length(files_dir)
    files_name{f_i} = files_dir(f_i).name;
    %## Check file naming format
    us_idx = strfind(files_name{f_i},'_');
    if isempty(us_idx)
        error('File format not as expected, check underscores')
    elseif length(us_idx) > 1
        disp('We found more than one underscore, we are assuming the last one is better');
        us_idx = us_idx(end);
    end
    %## Get time & data informatino from file name
    trial_name = files_name{f_i}(1:us_idx-1);
    date_str = files_name{f_i}(us_idx+1:us_idx+21); %example 19-09-16 03-01-00-269 for sept 16th 2019 at 3:01:00.269 PM
    rec_start_struct(f_i).yearStr = ['20',date_str(1:2)];
    rec_start_struct(f_i).monthStr = date_str(4:5);
    rec_start_struct(f_i).datStr = date_str(7:8);
    rec_start_struct(f_i).hourStr = date_str(10:11);
    tmp = str2num(rec_start_struct(f_i).hourStr);
    if tmp < 7 %file claims recording started before 7 am (no way man)
        tmp = tmp + 12; %adjust values to military time
        rec_start_struct(f_i).hourStr = num2str(tmp); %note: don't need to worry about ensuring the hour is 2 digits long to match original string format because tempHr will be between 12 and 19 after converting to military time
    end
    rec_start_struct(f_i).minStr = date_str(13:14);
    rec_start_struct(f_i).secStr = date_str(16:17);
    rec_start_struct(f_i).fracSecStr = date_str(19:21);
    rec_start_struct(f_i).datetime = datetime(date_str,'Format','yy-MM-dd HH-mm-ss-SSS');
    if hour(rec_start_struct(f_i).datetime) < 7 %file claims recording started before 7 am (no way man)
        rec_start_struct(f_i).datetime = rec_start_struct(f_i).datetime + hours(12); %adjust values to military time
    end
    rec_start_struct(f_i).datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    %## Load loadsol data
    fprintf('Loading %s from...\n%s\n',files_name{f_i},txt_file_path);
    tmp = importdata([txt_file_path filesep files_name{f_i}]);
    %## Get information from ASCII file
    %- common params
%     sample_rate = round(1/median(diff(tmp.data(:,1)')));
    world_time      = rec_start_struct(f_i).datetime + seconds(tmp.data(:,1)');
    device_time     = seconds(tmp.data(:,1)');
    tempLabels      = tmp.textdata{3,1};
    spaceInd        = strfind(tempLabels,' ');
    startNamesInd   = spaceInd((diff(spaceInd) ~= 1))+1;
    endNamesInd     = spaceInd(startNamesInd)-1;
    sensorList      = cell(1,length(startNamesInd));
    %- get sensors from txt file
    for name_i = 1:length(startNamesInd)
        sensorList{name_i} = tempLabels((startNamesInd(name_i):endNamesInd(name_i)));
    end
    %## Determine left, right, and sync foot
    leftLegs = find(contains(sensorList,'-L'));
    rightLegs = find(contains(sensorList,'-R'));
    
    try %try auto figuring out way
        if (length(leftLegs) == 2) && (length(rightLegs) == 1) %sync signal is a left leg
                rightFootNum = rightLegs;
                rightFootPartialName = sensorList{rightFootNum}(1:3); %e.g. P1U for a cetain size, P1V for another size
                leftFootNum = leftLegs( contains(sensorList(leftLegs),rightFootPartialName));
                syncFootNum = leftLegs( ~contains(sensorList(leftLegs),rightFootPartialName));
        elseif (length(leftLegs) == 1) && (length(rightLegs) == 2) %sync signal is a right leg
                leftFootNum = leftLegs;
                leftFootPartialName = sensorList{leftFootNum}(1:3); %e.g. P1U for a cetain size, P1V for another size
                rightFootNum = rightLegs( contains(sensorList(rightLegs),leftFootPartialName));
                syncFootNum = rightLegs( ~contains(sensorList(rightLegs),leftFootPartialName));
        else %Roehl
%                 error('I do not know what is going on');
                rightFootNum = rightLegs;
                leftFootNum = leftLegs;
                syncFootNum = setdiff(1:length(sensorList),[rightFootNum leftFootNum]);
        end
    catch 
        error('Check loadsol naming convention in .txt file')
    end
    
    %## Get data for each loadsol
    %- NOTE: Multiply by 2 because data set comes as the column vectors 
    %        [sensor1Time, sensor1Y, sensor2Time, sensor2Y, ...] where all 
    %        the sensor Times are the same

    leftFootData = tmp.data(:,leftFootNum*2)'; 
    rightFootData = tmp.data(:,rightFootNum*2)'; 
    syncFootData = tmp.data(:,syncFootNum*2)'; 
    
    %## Plot data
    %{
    fig_i = figure; 
    hold on;
    title('Raw Signals');
    worldTime = rec_start_struct.datetime + seconds(tmp.data(:,1)');
    plot(worldTime,leftFootData,'b');
    plot(worldTime,rightFootData,'r');
    plot(worldTime,syncFootData,'k');
    ylabel('Vertical force (N)')
    legend('Left','Right','Sync');
    hold off;
    saveas(fig_i,fullfile(save_dir, 'rawSignals.fig'));
    %movefile('rawSignals.fig', save_dir); %Moving to Reports
    %}
    loadsol_table{f_i} = table(device_time,world_time,leftFootData,rightFootData,syncFootData,'VariableNames',{'time_seconds','datetime','left_foot_N','right_foot_N','sync_foot_N'});
end


end


