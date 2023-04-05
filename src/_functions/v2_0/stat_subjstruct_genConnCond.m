function [SUBJSTRUCT] = stat_subjstruct_genConnCond(PATHS,SUBJSTRUCT,processNum_1,processNum_2,studyName,studyDir,varargin)
%STAT_SUBJSTRUCT_GENCONNCOND Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {true};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum_1',@isnumeric);
addRequired(p,'processNum_2',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,PATHS,SUBJSTRUCT,processNum_1,processNum_2,studyName,studyDir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%## PARAMS
ASSIGN_BOOTSTRAP_MEAN = false;
DO_BOOTSTRAP = true;
DO_PHASE_RND = false;
DO_NONZERO_TEST = false;
%## MAKE DIRS
mkdir(studyDir);
%## DETECT OS
try
    if strncmp(computer,'PC',2)
        DO_UNIX = false;
    else
        DO_UNIX = true;
    end
catch
    error('OSError:unknownOS','ERROR. You are working in an unknown Operating System.');
end
if DO_UNIX
    studyDir = convertPath2UNIX(studyDir,'dferris');
else
    studyDir = convertPath2Drive(studyDir,'M');
end
%% ===================================================================== %%
%## DIARY
diaryPath = [PATHS.path4diaries filesep studyName filesep '4_diary_genConnStats.txt'];
if ~exist([PATHS.path4diaries filesep studyName],'dir')
    mkdir([PATHS.path4diaries filesep studyName]);
end
diary(diaryPath)
%## EEGLAB
eeglab;
%% ASSIGN LABELS FOR PLOTTING
%- condition/trial name
ANALYSIS_NAME = SUBJSTRUCT(1).META.epoched.trialType;
%- cells to store components and clusters for each subject
components = cell(length(SUBJSTRUCT),1);
fs = fields(SUBJSTRUCT(1).cnctAnl.(ANALYSIS_NAME));
clusts = cell(length(SUBJSTRUCT),1);
comps = cell(length(SUBJSTRUCT),1);
clusters = cell(length(SUBJSTRUCT),1);
%-
for cnd = 1:length(SUBJSTRUCT)
    components{cnd} = SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).components;
    fprintf('==== Creating Labels for Subject %s ====\n',SUBJSTRUCT(cnd).subjStr);
    for j = 1:(length(fs)-1)
        clusts{cnd} = [ clusts{cnd,1} SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ClusterNumber'];
        comps{cnd} = [ comps{cnd,1} SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ComponentNumber'];
    end
    val = components{cnd};
    comps{cnd} = unique(comps{cnd});
    for i = 1:length(val)
        tmp = clusts{cnd}(val(i) == comps{cnd});
        clusters{cnd} = [clusters{cnd}, {sprintf('%i',tmp)}];
    end
    %- PRINT OUTS
    fprintf('Plotting Components: '); fprintf('%i,',comps{cnd}'); fprintf('\n');
    fprintf('Associated Clusters: '); 
    for i = 1:length(clusters{cnd})
        if isempty(clusters{cnd}{i})
            clusters{cnd}{i} = 'NaN';
        end
        fprintf('%s,',clusters{cnd}{i}');
    end
    fprintf('\n');
end
%% CONNECTIVITY
parfor cnd = 1:length(SUBJSTRUCT)
    ALLEEG = [];
    %## PARAMS
    %## SETUP
    %## LOOP PRINTS
    fprintf('==== Processing %s ====', SUBJSTRUCT(cnd).subjStr);
    %- Load condition 1  
    fprintf('Loading %s',SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_1).label);
    fPath = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_1).filepath;
    fName = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_1).filename;
    if DO_UNIX
        fPath = convertPath2UNIX(fPath,'dferris');
    else
        fPath = convertPath2Drive(fPath,'M');
    end
    %## LOAD EEG
    EEG = par_load(fPath,fName,[]);
    EEG = cnctanl_loadCAT(EEG,SUBJSTRUCT,processNum_1,'Bootstrap')
    ALLEEG(1) = EEG;
    %- Load condition 2
    fprintf('Loading %s',SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_2).label);
    fPath = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_2).filepath;
    fName = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum_2).filename;
    if DO_UNIX
        fPath = convertPath2UNIX(fPath,'dferris');
    else
        fPath = convertPath2Drive(fPath,'M');
    end
    %## LOAD EEG
    EEG = par_load(fPath,fName,[]);
    EEG = cnctanl_loadCAT(EEG,SUBJSTRUCT,processNum_1,'Bootstrap')
    ALLEEG(2) = EEG;
    %% GENERATE BETWEEN CONDITION STATS
    ALLEEG = stat_cnctanl_connCond(PATHS,ALLEEG,studyDir,studyName,connMethods,components,clusters);
    %% ==== COMPLEX VISUALS ==== %%
    %- assign saveDir
    saveDir = [studyDir filesep EEG.subject];
    if ~exist(saveDir,'dir')
        fprintf(1,'Making Directory: %s',saveDir);
        mkdir(saveDir)
    end    
    %% EX-TIMEFREQCHART (STATS)   
    [~,~] = visualCnctAnl(ALLEEG,PATHS,components,clusters,...
        'exTimeFreqChart',...
        'connMeasures',connMethods,...
        'bounds',{},...
        'savePath',saveDir);
    fprintf('%s) Done.',SUBJSTRUCT(cnd).subjStr)
end
diary('off');
end

