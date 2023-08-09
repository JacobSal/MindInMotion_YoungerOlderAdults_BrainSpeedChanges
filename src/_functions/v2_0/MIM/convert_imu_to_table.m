function [imu_table,rec_start_struct] = convert_imu_to_table(imu_file_path,varargin)
%CONVERT_IMU_TO_STRUCT Summary of this function goes here
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
addRequired(p,'imu_file_path',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p, imu_file_path, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- PERMS
SIGS_PER_IMU = 14;
%% ===================================================================== %%
file_list = dir([imu_file_path filesep '*.csv']);
imu_table = cell(1,length(file_list));
% numFiles = length(tempFileList); 
% disp(['We found ',num2str(numFiles),' files.']);
rec_start_struct = struct('yearStr',[],...
                        'monthStr',[],...
                        'dayStr',[],...
                        'hourStr',[],...
                        'minStr',[],...
                        'secStr',[],...
                        'fracSecStr',[],...
                        'datetime',[],...
                        'trialName',[],...
                        'filename',[],...
                        'filepath',[]);
for f_i = 1:length(file_list)
    fName = file_list(f_i).name;
    fPath = imu_file_path;
    
    %## Get time & data informatino from file name
    dashInd = strfind(fName,'-');
    fullDateStr = fName(1:dashInd(2)-1); %example 19-09-16 03-01-00-269 for sept 16th 2019 at 3:01:00.269 PM
    trialName = fName(dashInd(2)+1:end-4);
    %- store information
    rec_start_struct(f_i).trialName = trialName;
    rec_start_struct(f_i).yearStr = fullDateStr(1:4);
    rec_start_struct(f_i).monthStr = fullDateStr(5:6);
    rec_start_struct(f_i).dayStr = fullDateStr(7:8);
    
    %- make sure to convert to military time (e.g. 1pm = 13). 
    %  Otherwise there is confusion between AM/PM for the LoadSol files
    rec_start_struct(f_i).hourStr = fullDateStr(10:11); 
    tempHr = str2num(rec_start_struct(f_i).hourStr);
    
    %- file claims recording started before 7 am (no way man)
    if tempHr < 7 
        tempHr = tempHr + 12; %adjust values to military time
        %- NOTE: don't need to worry about ensuring the hour is 2 digits 
        %        long to match original string format because tempHr will 
        %        be between 12 and 19 after converting to military time
        rec_start_struct(f_i).hourStr = num2str(tempHr); 
    end
    rec_start_struct(f_i).minStr = fullDateStr(12:13);
    rec_start_struct(f_i).secStr = fullDateStr(14:15);
    
    %- probably not needed. done to match liveamp
    %  ryan test new datetime code
    rec_start_struct(f_i).fracSecStr = '000';
    rec_start_struct(f_i).datetime = datetime(fullDateStr,'Format','yyyyMMdd-HHmmss');
    if hour(rec_start_struct(f_i).datetime) < 7 %file claims recording started before 7 am (no way man)
        rec_start_struct(f_i).datetime = rec_start_struct(f_i).datetime + hours(12); %adjust values to military time
    end
    rec_start_struct(f_i).datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
%     disp(recordingStart);

    %## Load file & determine a precise start time
    tmp = importdata([fPath filesep fName]);
    
    %- figure out start time to millisecond precision (i.e. using data 
    %  containedin the file itself rather than just the filename which 
    %  only goes to seconds precision)
    timeRaw = tmp.data(:,1)'; %microseconds since Jan 1 1970? 
    tStartGlobal = datetime(timeRaw(1)*1E-6, 'ConvertFrom', 'posixtime');
    tStartGlobal.Format =  'yyyy-MM-dd HH:mm:ss.SSS';
    
    %## Determine start time of IMU file
    tempDiff = tStartGlobal - rec_start_struct(f_i).datetime;
    
    %- mod 3600 accounts for some integer number of hours difference 
    secondsShift = mod(seconds(tempDiff),3600); 
    
    %- the two ways of calculating the start time should generally agree 
    %  within a second but sometimes things just arent perfect in life so 
    %  lets give him some leeway with 2 secondsNote though that the hour 
    %  gets shifted because of UTC time)
    
    %- make recordingStart more accurate by incorporating the posixtime 
    %  (UTC) data from first column of the excell sheet

    %- 2021-04-30 RJD increased secondsShift threshold from 2 to 5 seconds to
    %                 accomodate for H2025_med_1

    if secondsShift  <= 5 
        rec_start_struct(f_i).datetime = rec_start_struct(f_i).datetime + seconds(secondsShift); %shift by the tiny offset but ignore the shift of hours (because one was UTC and another was Eastern time)
        rec_start_struct(f_i).datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
        disp(['secondsShift = ', num2str(secondsShift)]);
    end
    
    %- 2021-07-22 added this logic b/c special situation with 'H3034-SP_p75_3'
    %  where the correct offset between the posixtime and the file name 
    %  timestamp was very tiny but negative so mod() resulted in 3599.9979 
    %  which is really -0.0021. Surprised we didn't encounter this before
    if secondsShift >= 3595 && secondsShift < 3600 
        disp('Applying special correction for small negative shift');
        disp(['secondsShift was ', num2str(secondsShift)]);
        secondsShift = secondsShift-3600;
        disp(['secondsShift is now ', num2str(secondsShift)]);
        rec_start_struct(f_i).datetime = rec_start_struct(f_i).datetime + seconds(secondsShift); 
        rec_start_struct(f_i).datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    end
    
    %- 2022-10-25 added this logic b/c special situation with 'H1037-SP_p25_1'
    %  where the correct offset between the posixtime and the file name 
    %  timestamp was very tiny but negative so mod() resulted in 3507 
    if secondsShift >= 3504 && secondsShift < 3507 
        disp('Applying special correction for small negative shift');
        disp(['secondsShift was ', num2str(secondsShift)]);
        secondsShift = secondsShift-3507;
        disp(['secondsShift is now ', num2str(secondsShift)]);
        rec_start_struct(f_i).datetime = rec_start_struct(f_i).datetime + seconds(secondsShift); 
        rec_start_struct(f_i).datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';    
    end
    
    if secondsShift > 3600
        warndlg('There is a disagreement between the start time as it appears in the first column of the raw CSV (posixtime) and what is listed in the filename (datestring)');
        warning('There is a disagreement between the start time as it appears in the first column of the raw CSV (posixtime) and what is listed in the filename (datestring)');
        warning(['secondsShift = ', num2str(secondsShift)]);
        warning('Skipping importing this file because we are scared about the disagreement being so big');
        return;
    end

    %## Fix text data
    %- newer versions of matlab can use split function directly 
    strA = split(tmp.textdata(1,1),',')'; 
    
    %- However, if your matlab version is stupid and the split function only 
    %  returns back string arrays not also cell arrays of chars then...
    if isstring(strA) 
        for str_i = 1:length(strA)
            tmp.textdata(1,str_i) = {char(strA(str_i))}; %turn each individual string element to a char cell
        end
    elseif iscell(strA)
        tmp.textdata(1,:) = split(tmp.textdata(1,1),',')'; %fix text issues
    else
        disp('I could not figure out how to parse the IMU text headers');
    end
    
    %- newer versions of matlab can use split function directly
    strB = split(tmp.textdata(2,1),',')';
    
    %- However, if your matlab version is stupid and the split function 
    %  only returns back string arrays not also cell arrays of chars then...
    if isstring(strB) 
        for str_i = 1:length(strB)
            tmp.textdata(2,str_i) = {char(strB(str_i))}; %turn each individual string element to a char cell
        end
    elseif iscell(strB)
        tmp.textdata(2,:) = split(tmp.textdata(2,1),',')'; %fix text issues
    else
        disp('I could not figure out how to parse the IMU text headers');
    end
    %Find out how many signals and IMUs we are dealing with
    numSigs = size(tmp.data,2); %first entry is time FYI
    
    %- each imu has acc_x/y/z gyro_x/y/z mag_x/y/z barom orient_x/y/z. 
    %  numSigs-1 to account for first entry being time
    numIMU = (numSigs-1)/SIGS_PER_IMU; 
%     numPts = size(tmp.data,1);

    %- Go thru each IMU (more for first participant, and only one for everyone
    %  else as of 4/13/2020
    tempChanNames = cell(1,numSigs);
    t_ind = 1;
    tempChanNames{t_ind} = 'time';

    acc_all = [];
    gyro_all = [];
    mag_all = [];
    barom_all = [];
    orient_all = [];

    for IMU_i = 1:numIMU
        %- account for time being first entry and each imu having 14 signals.
        indOffset = 1 + (IMU_i-1)*14;  
        
        %- generate table variable names
        imu_name =  tmp.textdata{1,4 + indOffset}; 
        imu_name = join(strsplit(imu_name,' '),'_'); imu_name = imu_name{1};
        
        %- assign indices for each IMU available
        acc_ind = (1:3) + indOffset;
        gyro_ind = (4:6) + indOffset;
        mag_ind = (7:9) + indOffset;
        barom_ind = (10) + indOffset; 
        orient_ind = (11:14) + indOffset; 

        acc_all = [acc_all, acc_ind];
        gyro_all = [gyro_all, gyro_ind];
        mag_all = [mag_all, mag_ind];
        barom_all = [barom_all, barom_ind];
        orient_all = [orient_all, orient_ind];

        tempChanNames(indOffset*ones(size(1:14)) + (1:14)) = {[imu_name,'_acc_x'], [imu_name,'_acc_y'], [imu_name,'_acc_z'], ...
                                                                [imu_name,'_gyro_x'], [imu_name,'_gyro_y'], [imu_name,'_gyro_z'], ...
                                                                [imu_name,'_mag_x'], [imu_name,'_mag_y'], [imu_name,'_mag_z'], ...
                                                                [imu_name,'_barom'],...
                                                                [imu_name,'_ori_s'], [imu_name,'_ori_x'], [imu_name,'_ori_y'], [imu_name,'_ori_z']};
    end
    
    imu_table{f_i} = array2table(tmp.data,'VariableNames',tempChanNames);
end

end

%% SUBFUNCTIONS
