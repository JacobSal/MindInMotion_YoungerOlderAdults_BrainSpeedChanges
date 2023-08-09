
function [ do_crop_imu, exact_crop_imu, do_crop_ls, exact_crop_ls ] = mim_check_trials(subject_name,trial_name,input_dir,varargin)
%MIM_CHECK_TRIALS - looks through a master Excel sheet to see if we
% need to ignore the current trial. Crops trials based on 'TrialCrop'
% sheet, and removes gait events based on 'Loadsol' sheet.
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
% Copyright (C) Ryan Downey, Chang Liu, Jacob Salminen, jsalminen@ufl.edu
% Created by Ryan Downey - 01/01/2021
% Modified by Chang Liu - 07/14/2022
% Modified by Jacob Salminen - 02/09/2022

%% DEFINE DEFAULTS
% MasterTable = [];
% LoadsolTable = [];
MIM_R_FOLDER = ['R:' filesep 'Ferris-Lab' filesep 'share' filesep 'MindInMotion'];
%## PARSE
p = inputParser;
%## REQUIRED
addRequired(p,'subject_name',@ischar);
addRequired(p,'trial_name',@ischar);
addRequired(p,'input_dir',@ischar);
%## OPTIONAL
addParameter(p,'MIM_R_FOLDER',MIM_R_FOLDER,@ischar);
%## PARAMETER
parse(p, subject_name, trial_name, input_dir, varargin{:});
%## SET DEFAULTS
%- OPTIONALS

%- PARAMETER

%- PERMS
%% ===================================================================== %%
persistent MasterTable
persistent LoadsolTable

if isempty(input_dir)
    disp('couldn''t find the excel doc so the function can''t work')
end
if isempty(MasterTable)
    MasterTable = readtable(input_dir,'Sheet','trialCrop','Range','B2:V127'); %loads it up
    MasterTable.Properties.VariableNames(1:21) = {'SortingNum','SubjCode','SP_0p5_1',...
        'SP_0p5_2',	'SP_0p25_1','SP_0p25_2','SP_0p75_1','SP_0p75_2','SP_1p0_1',	'SP_1p0_2',	...
        'TM_flat_1','TM_flat_2','TM_high_1','TM_high_2','TM_low_1',	'TM_low_2',	'TM_med_1',	'TM_med_2',	'Rest','MotorImagery_1','MotorImagery_2'};
end
SubjectMatch = contains(MasterTable{:,'SubjCode'},subject_name); %look thru the subj to find match
SubjectMatchInd = find(SubjectMatch); %grab corresponding index

if isempty(SubjectMatchInd)
    do_crop_imu = false; %haven't ran the trial so there is not even a match
    exact_crop_imu = []; %2021-12-07 RJD added to avoid error in analyzeLoadsolManySubj script
else
    Result = MasterTable{SubjectMatchInd,trial_name}; %result as a cell
    ResultString = Result{1}; %the result as a string
    exact_crop_imu = str2num(ResultString); %result as a vector. the exact seconds to crop it by
    %should get something in format of [boundary boundary]

    %2021-06-16 RJD replaced DoCrop  logic because it could not handle
    %excel having an empty cell (i.e. ResultString = '', ExactCrop = [])
    if isempty(exact_crop_imu)
        do_crop_imu = false;
    else
        do_crop_imu = true;
    end

%         if ResultString(1) == '[' && ResultString(length(ResultString)) == ']' %kinda bad logic but also kinda good 
%             DoCrop = true;
%         else
%             DoCrop = false;
%         end

end
%% 
if isempty(LoadsolTable)
%     LoadsolTable = readtable(input_dir,'Sheet','loadsolCrop'); %loads it up
    LoadsolTable = readtable(input_dir,'Sheet','trialCrop','Range','B2:V127'); %loads it up
    LoadsolTable.Properties.VariableNames(1:21) = {'SortingNum','SubjCode','SP_0p5_1',...
        'SP_0p5_2',	'SP_0p25_1','SP_0p25_2','SP_0p75_1','SP_0p75_2','SP_1p0_1',	'SP_1p0_2',	...
        'TM_flat_1','TM_flat_2','TM_high_1','TM_high_2','TM_low_1',	'TM_low_2',	'TM_med_1',	'TM_med_2',	'Rest','MotorImagery_1','MotorImagery_2'};
end
SubjectMatch = contains(LoadsolTable{:,'SubjCode'},subject_name); %look thru the subj to find match
SubjectMatchInd = find(SubjectMatch); %grab corresponding index

if isempty(SubjectMatchInd)
    do_crop_ls = false; %haven't ran the trial so there is not even a match
    exact_crop_ls = []; %2021-12-07 RJD added to avoid error in analyzeLoadsolManySubj script
else
    Result = LoadsolTable{SubjectMatchInd,trial_name}; %result as a cell
    ResultString = Result{1}; %the result as a string
    exact_crop_ls = str2num(ResultString); %result as a vector. the exact seconds to crop it by
    %should get something in format of [boundary boundary]

    %2021-06-16 RJD replaced DoCrop  logic because it could not handle
    %excel having an empty cell (i.e. ResultString = '', ExactCrop = [])
    if isempty(exact_crop_ls)
        do_crop_ls = false;
    else
        do_crop_ls = true;
    end

%         if ResultString(1) == '[' && ResultString(length(ResultString)) == ']' %kinda bad logic but also kinda good 
%             DoCrop = true;
%         else
%             DoCrop = false;
%         end

end
end

