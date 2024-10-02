function [MasterTable] = mim_read_master_sheet(varargin)
%MIM_READ_MASTER_SHEET - looks through a master Excel sheet to see if we
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
% Modified by Jacob Salminen - 08/22/2023

%% DEFINE DEFAULTS
input_fPath = 'M:\jsalminen\GitHub\par_EEGProcessing\src\_data\MIM_dataset\_studies\subject_mgmt\subject_logging.xlsx';
if ~ispc
    input_fPath = convertPath2UNIX(input_fPath);
else
    input_fPath = convertPath2Drive(input_fPath);
end
SHEET_NAME = 'table_reduced_inf';
%## PARSER OBJ
p = inputParser;
%## REQUIRED

%## OPTIONAL
addOptional(p,'input_fPath',input_fPath,@ischar);
%## PARAMETER
addParameter(p,'SHEET_NAME',SHEET_NAME,@ischar);
%## PARSE
parse(p, varargin{:});
%## SET DEFAULTS
%- OPTIONALS
input_fPath = p.Results.input_fPath;
%- PARAMETER
SHEET_NAME = p.Results.SHEET_NAME;
%- PERMS
%% ===================================================================== %%
if ~exist(input_fPath,'file')
    error('''input_dir'' does not exist: %s',input_fPath)
end
MasterTable = readtable(input_fPath,'Sheet',SHEET_NAME); %loads it up
end

