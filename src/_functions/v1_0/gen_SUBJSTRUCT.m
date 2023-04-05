function [SUBJSTRUCT] = getSubjsDirData(PATHS, subjNums, varargin)
%GETSUBJDIRDATA Summary of this function goes here
%   Detailed explanation goes here
%   
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {false, false,...
    [], ...
    'TmpStudy', ...
    [], ...
    'HModel_JS',...
    false};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct)
addRequired(p,'subjNums')
%## OPTIONAL
addOptional(p,'studyName',Defaults{4},@ischar)
%## PARAMETER
validPath = @(x) ischar(x) || isempty(x);
addParameter(p,'saveData',Defaults{2},@islogical) % logical determining whether data will be saved for certain steps
addParameter(p,'savePath',Defaults{3},validPath)
addParameter(p,'modelFile',Defaults{6})
addParameter(p,'loadData',Defaults{7})
parse(p, PATHS, subjNums, varargin{:});
%## SET DEFAULTS
saveData = p.Results.saveData;
savePath = p.Results.savePath;
studyName = p.Results.studyName; 
modelFile = p.Results.modelFile;
loadData = p.Results.loadData;

%## SAVE FOLDER HANDLER (HIGHEST ORDER FUNCTION ONLY)
if isempty(savePath)
    savePath = [PATHS.localStorage 'subjinf'];
else
    savePath = [savePath filesep 'subjinf'];
end
tmpSave = fullfile(savePath,sprintf('SUBJSTRUCT_%s.mat',studyName));
% ----------------------------------------------------------------------- %
%## Get all subject path data and organize into health and non-healthy
[allSubj,~,idxval] = getDirInfMIM(PATHS);
[Hidx, NHidx, idxA1, idxA2, idxA3] = idxval{:};

%## Extract IMU and LoadSol information for participants
disp('Assigning patient data to subject strucutre');
IMUstats = struct('cond', [], ...
                  'trialOM', [], ...
                  'trailSbS', [], ...
                  'condPath', '', ...
                  'trialOMPath', '', ...
                  'trialSbSPath', '');
LoadsolStats = struct('cond', [], ...
                      'trialOM', [], ...
                      'trailSbS', [], ...
                      'condPath', '', ...
                      'trialOMPath', '', ...
                      'trialSbSPath', '');
EEGstats = struct('icaDirNum',[],...
                  'icaPath','', ...
                  'icaSETname',[],...
                  'rawEEGpath','',...
                  'trialEEGpath','',...
                  'headModelPath','', ...
                  'elecPath',''); %...
%                   'hpgCleanPath',[]); % may implement at a later date (07/05/2022)
icaSETname = [];
icaDirNum = [];
SourceLocalize = struct('componentsAll',[],...
                        'componentsKeep',[],...
                        'mrifile',[],...
                        'hasAcpc',0,...
                        'hasFid',0,...
                        'acpcNames',{{'ac','pc','xzpoint','rpoint'}},...
                        'fidNames',{{'nasion','lhj','rhj','xzpoint'}},...
                        'acpcCoords',{{[0,0,0],[0,0,0],[0,0,0],[0,0,0]}},...
                        'fidCoords',{{[0,0,0],[0,0,0],[0,0,0],[0,0,0]}});
CnctAnl = struct('filepath',[],...
                 'filename',[],...
                 'cluster',[],...
                 'components',[]);
SUBJSTRUCT = struct('SubjectStr',"",...
                       'Nback',0, ...
                       'ImaginedWalking',0, ...
                       'IMUstats', IMUstats, ...
                       'LoadsolStats', LoadsolStats, ...
                       'Health', '', ...
                       'AgeInt', 0, ...
                       'EEGstats',EEGstats,...
                       'SourceLocalize',SourceLocalize,...
                       'CnctAnl',CnctAnl);
%## check if subject numbers are provided
if isempty(subjNums)
    subjNums = (1:length(allSubj));
end

for i = 1:length(subjNums)
    subjStr = allSubj(subjNums(i)).name;
    subjP = allSubj(subjNums(i)).folder;
    SUBJSTRUCT(i).SubjectStr = subjStr;
    fprintf('%i) Processing subject: %s\n',subjNums(i),subjStr)
    % Does the subject have N-back?
    tmp = dir(fullfile(subjP,subjStr,'EPrime','N-back'));
    SUBJSTRUCT(i).Nback = sum([tmp.bytes]) > 0;
    % Does the subject have Imagined?
    tmp = dir(fullfile(subjP,subjStr,'EPrime','Imagined'));
    SUBJSTRUCT(i).ImaginedWalking = sum([tmp.bytes]) > 0;
    % See what ICA's were performed
    SUBJSTRUCT(i).EEGstats = EEGstats;
    tmp = dir(fullfile(subjP,subjStr,'EEG','AMICA'));
    if ~isempty(tmp)
        icaNums = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
        SUBJSTRUCT(i).EEGstats.icaDirNum = icaNums;
        SUBJSTRUCT(i).EEGstats.icaPath = tmp.folder;
        icaSname = cell(length(icaNums),1);
        % Gather the names of all ICA decomps to find a common one for the STUDY
        % group
        for j = 1:length(icaNums)
            tmpi = dir([tmp(j).folder filesep icaNums{j}]);
            tmpo = cellfun(@(x) regexp(x,'.set'), {tmpi.name}, 'UniformOutput', false);
            tmpo = cellfun(@(x) x > 0, tmpo, 'UniformOutput', false);
            index = false(1, numel(tmpo));
            for k = 1:numel(tmpo)
              index(k) = ~isempty((tmpo{k}));
            end
            tmpo = strsplit(tmpi(index).name,'.'); tmpo = strsplit(tmpo{1},'_');
            if length(tmpo) > 2
                outi = strjoin(tmpo(2:end),'_');
            elseif length(tmpo) > 1 
                outi = tmpo{2};
            else
                outi = tmpo{:};
            end
            icaSname{j} = outi;        
        end
        SUBJSTRUCT(i).EEGstats.icaSETname = icaSname;
    elseif ~strcmp(studyName,'ALL')
        try
            mkdir(fullfile(subjP,subjStr,'EEG','AMICA'))
        catch
            error(['Subject ''%s'' does not have a viable ICA decomposition, please',...
                '\ngenerate ICA file then run this script again.'],subjStr)
        end
    end
    % add paths for RAW and TRIAL EEG data.
    tmp = dir(fullfile(subjP,subjStr,'EEG','Raw'));
    try
        SUBJSTRUCT(i).EEGstats.rawEEGpath = tmp(1).folder;
    catch
        warning('Raw data not collected for this subject');
    end
    tmp = dir(fullfile(subjP,subjStr,'EEG','Trials'));
    try
        SUBJSTRUCT(i).EEGstats.trialEEGpath = tmp(1).folder;
    catch
        warning(['Trial data not processed for this subject...\n',...
            'use 0_preProcessRaw/11_extractTrials before running ICA']);
    end
    
    % find if there is a headscan/headmodel available for subject
    tmp = dir(fullfile(subjP,subjStr,'HeadScan',modelFile));
    if sum([tmp.bytes]) > 0
        SUBJSTRUCT(i).EEGstats.headModelPath = fullfile(tmp(strcmp({tmp.name},'vol.mat')).folder,'vol.mat');
        SUBJSTRUCT(i).EEGstats.elecPath = fullfile(tmp(strcmp({tmp.name},'elec_aligned_init.mat')).folder,'elec_aligned_init.mat');
    else
        warning('A valid Headmodel and electrode .mat file is unavailable for participant ''%s''', subjStr)
    end
    % attach conditions/trial information  and path to conditions
    % IMU
    SUBJSTRUCT(i).IMUstats = IMUstats;
    tmp = dir(fullfile(subjP,subjStr,'IMU','Processed_Conditions'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).IMUstats.cond = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
        SUBJSTRUCT(i).IMUstats.condPath = tmp.folder;
    end
    tmp = dir(fullfile(subjP,subjStr,'IMU','Processed_Trials','Stride_by_Stride_Structs'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).IMUstats.trialSbSPath = tmp.folder;
        SUBJSTRUCT(i).IMUstats.trialSbS = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
    end
    tmp = dir(fullfile(subjP,subjStr,'IMU','Processed_Trials','Outcome_Measures'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).IMUstats.trialOM = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
        SUBJSTRUCT(i).IMUstats.trialOMPath = tmp.folder;
    end
    % LOADSOL
    SUBJSTRUCT(i).LoadsolStats = LoadsolStats;
    tmp = dir(fullfile(subjP,subjStr,'Loadsol','Processed_Conditions'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).LoadsolStats.cond = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
        SUBJSTRUCT(i).LoadsolStats.condPath = tmp.folder;
    end
    tmp = dir(fullfile(subjP,subjStr,'Loadsol','Processed_Trials','Stride_by_Stride_Gait_Structs'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).LoadsolStats.trialSbSPath = tmp.folder;
        SUBJSTRUCT(i).LoadsolStats.trialSbS = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
    end
    tmp = dir(fullfile(subjP,subjStr,'Loadsol','Processed_Trials','Outcome_Measures'));
    if ~isempty(tmp)
        SUBJSTRUCT(i).LoadsolStats.trialOM = cellfun(@(x) x,{tmp(3:end).name},'UniformOutput',false);
        SUBJSTRUCT(i).LoadsolStats.trialOMPath = tmp.folder;
    end 
    % HEALTH AND AGE DESIGNATION
    SUBJSTRUCT(i).Health = char(Hidx(i+3)*'He' + NHidx(i+3)*'NH');
    SUBJSTRUCT(i).AgeInt = idxA1(i+3)*1 + idxA2(i+3)*2 + idxA3(i+3)*3;
    
    % ASSIGN commonICAs empty struct for later use
    SUBJSTRUCT(i).CommonICAs = CommonICAs;
    
    % ASSIGN sourceLocalize empty struct for later use
    SUBJSTRUCT(i).SourceLocalize = SourceLocalize;
    tmp = dir(fullfile(subjP,subjStr,'MRI','Raw'));
    if ~isempty(tmp) && length(tmp) > 2  % added 'length(tmp)>2' (06/29/22) because of hidden folders
        SUBJSTRUCT(i).SourceLocalize.mrifile = fullfile(tmp(1).folder,'T1.nii');
    else
        SUBJSTRUCT(i).SourceLocalize.mrifile = [];
    end
    % ASSIGN path for subject struct     
    SUBJSTRUCT(i).PATHS.studyPath = tmpSave;
    SUBJSTRUCT(i).PATHS.PATHS = PATHS;
    SUBJSTRUCT(i).PATHS.tmpICA = [];
    SUBJSTRUCT(i).PATHS.tmpTrial = [];
    
    % ASSIGN cnctAnl for subject struct
    SUBJSTRUCT(i).CnctAnl = CnctAnl;
end
%## CHECK 'ALL' SUBJSTRUCT FOR ALREADY STORED ANALYSES
if ~strcmp(studyName,'ALL')    
    [SUBJSTRUCT] = checkSubjDirData(PATHS,SUBJSTRUCT,'dataDir',savePath,'saveData',false);
end
%## GET COMMON ICAs
if ~strcmp(studyName,'ALL')
    [SUBJSTRUCT] = findCommonICA(SUBJSTRUCT); 
end
%## CHECK FOR EXISTING SUBJSTRUCT
if ~isfolder(fileparts(tmpSave))
    mkdir(fileparts(tmpSave))
end
%## WRAP UP & SAVE
if saveData % && ~isfile(tmpSave) % Making funcitonallity more specific (07/01/2022), JS
    save(tmpSave,'SUBJSTRUCT');
elseif loadData && isfile(tmpSave)
    fprintf('Study ''%s'' already exists. Loading...\n',studyName);
    tmpL = load(tmpSave,'SUBJSTRUCT');
    SUBJSTRUCT = tmpL.SUBJSTRUCT;
else
    warning('Study ''%s'' does not exist, and saving was not specified.\nNo file saved.\n',studyName)
end

%## TIME
toc
end

