%GETDATA4SUBJN Loads the ICA & EEG struct for the subject designated by
%subjStr & subDirNum.
%   IN:
%       subjStr, CHAR
%           String of subject you are trying to load
%       subDirNum, NUMERIC
%           The sub-directory number for the AMICA/ICA analysis you ran and
%           would like to analyze
%       ALLEEG, STRUCT
%           Struct containing any previously loaded EEG data
%       CURRENTSET, NUMERIC
%           Number representing the current dataset you are working on in
%           ALLEEG
%       MIMDataFolder, PATH/CHAR
%           Path to the Mind-In-Motion data folder located either on your
%           computer (local temp folder) or on the R:/ drive.
%       overWriteCurrent, LOGICAL
%           logical that determines if the current set will overwrite the
%           previously loaded dataset. false, will not overwrite. true,
%           will overwrite.
%   OUT:
%       ALLEEG, STRUCT
%           Struct containing the newly loaded EEG dataset(s)
%       EEG, STRUCT
%           Struct containing the newly loaded EEG dataset
%       CURRENTSET, NUMERIC
%           Number that represents the dataset that is loaded.
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220112.0
%   Previous Version: n/a
%   Summary:  
%
function [ALLEEG,EEG,CURRENTSET] = subj_loadICA(SUBJ,processNum,varargin)
% SUBJ_LOADICA this loads
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {[]};
p = inputParser;
%## REQUIRED
addRequired(p,'SUBJ',@isstruct) ;
addRequired(p,'processNum',@isnumeric);
%## OPTIONAL
addOptional(p,'fName',Defaults{1},@ischar);
addOptional(p,'fPath',Defaults{1},@ischar);
%## PARAMETER
parse(p, SUBJ, processNum, varargin{:});
%## SET DEFAULTS
fName = p.Results.fName;
fPath = p.Results.fPath;
%-
OVER_WRITE_CURRENT_SET = 1;
ALLEEG = [];
CURRENTSET = 1;
if isempty(fName)
    fName = SUBJ.pathsEEG.ica(processNum).filename;
    warning('subj_loadICA.m) Using filename: %s',fName)
end
if isempty(fPath)
    fPath = SUBJ.pathsEEG.ica(processNum).filepath;
    warning('subj_loadICA.m) Using filepath: %s',fPath)
end
%% ===================================================================== %%
%## load EEG file contianing AMICA results
EEG = pop_loadset('filename',fName,'filepath',fPath);
if isempty(EEG.icaweights)
    %## CONCATENATE DATA | OVERWRITE DATA
    if OVER_WRITE_CURRENT_SET    
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
    else
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    end

    %## REMOVE ICA RESULTS FROM PREVIOUS ANALYSIS
    EEG = rmfield(EEG,'icaweights');
    EEG.icaweights = [];
    EEG = rmfield(EEG,'icawinv');
    EEG.icawinv = [];
    EEG = rmfield(EEG,'icaact');
    EEG.icaact = [];
    EEG = rmfield(EEG,'icasphere');
    EEG.icasphere = [];
    EEG = rmfield(EEG,'icachansind');
    EEG.icachansind = []; %also remove icachanind?

    %## LOAD CURRENT AMICA RESULTS FOR subjStr & subDirNum
    EEG = pop_loadmodout(EEG,fPath); 
    % make sure everything is ok and calculate the activation matrix
    pop_editoptions('option_computeica',1); %randomly needed for some people to run even if you manually selected this checkbox with the GUI
    EEG = eeg_checkset(EEG, 'ica'); %Note: I have no idea why it says "LLt not set" or what LLt is

    disp(['HEY YOU! Yeah you! You got ',num2str(size(EEG.etc.amica.W,3)),' AMICA model(s) to choose from!']);
else
    fprintf('ica already available.\n');
end
%% store updated EEG structure
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );

%## TIME
toc

end
%% Version History
%{
v1.0.20210112.0, JS: initial configuration. Script adapted from
addAMICAresults

%}