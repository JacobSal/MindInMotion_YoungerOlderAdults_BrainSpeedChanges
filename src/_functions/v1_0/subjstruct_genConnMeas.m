function [SUBJSTRUCT] = subjstruct_genConnMeas(PATHS,SUBJSTRUCT,processNum,studyName,studyDir,connMethods,brainChars,brainCoords,varargin)
%GENCONNMEASSUBJS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {10};
p = inputParser; 
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
addRequired(p,'connMethods',@iscell);
addRequired(p,'brainChars',@iscell);
addRequired(p,'brainCoords',@iscell);
%## 
addOptional(p,'anatThresh',Defaults{1},@isnumeric);
%## PARAMETER

%## PARSE
parse(p,PATHS,SUBJSTRUCT,processNum,studyName,studyDir,connMethods,brainChars,brainCoords,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
ANAT_THRESH = p.Results.anatThresh; 
%- PARAMETER
%- Define Defaults
%## PARAMS 
COMP_THRESH = 1; % minimum 1, maximum inf 
REGEN_FILES = false;
DO_BOOTSTRAP = true;
WINDOW_LENGTH = 0.5;
WINDOW_STEP_SIZE = 0.025;
NEW_SAMPLE_RATE = [];
RV_THRESH = 0.15;
BRAIN_THRESH = 0.5;
CONSTRAIN_TO_ANATOMY = false;
REJECT_SUBJECTS = true;
ASSIGN_BOOTSTRAP_MEAN = false;
BRAIN_STRUCT = struct('AnatTable',[]);
components = cell(length(SUBJSTRUCT),1);
ANALYSIS_NAME = SUBJSTRUCT(1).META.epoched.trialType;
CONN_METHODS = connMethods;
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
%% ===================================================================== %%
%## DIARY
diaryPath = [PATHS.path4diaries filesep studyName filesep '3_diary_genConnMeas.txt'];
if ~exist([PATHS.path4diaries filesep studyName] ,'dir')
    mkdir([PATHS.path4diaries filesep studyName]);
end
diary(diaryPath)
%## EEGLAB
eeglab;
%% CONNECTIVITY
%Load STUDY (assuming it exists)
fprintf('\n==== LOADING STUDY ====\n')
if DO_UNIX
    studyDir = convertPath2UNIX(studyDir,'dferris');
    [STUDY, ALLEEG] = pop_loadstudy('filename',[studyName,'_UNIX','.study'],'filepath',studyDir);
else
    studyDir = convertPath2Drive(studyDir,'M');
    [STUDY, ALLEEG] = pop_loadstudy('filename',[studyName,'.study'],'filepath',studyDir);
end
%% 
fprintf('\n==== ASSIGNING ANATOMY ====\n')
% ==== DEPRECTAED ==== %
% distOut = zeros(length(STUDY.cluster)-1,length(BRAIN_COORDS));
% setNum = cell(length(STUDY.cluster)-1,1);
% for i = 1:length(BRAIN_COORDS)    
%     cnt = 1;
%     for j = 2:length(STUDY.cluster) % skip parentcluster
%         [~,cntr] = std_centroid(STUDY,ALLEEG,j,'dipole');
%         cntr = cntr{1};
%         cntrCoord = cntr.dipole.posxyz;
%         distOut(cnt,i) = sqrt(sum((cntrCoord-BRAIN_COORDS{i}).^2));
%         setNum{cnt} = unique(STUDY.cluster(j).sets);
%         cnt = cnt + 1;        
%     end
% end
% %## Subject Rejection
% [M,I] = min(distOut);
% I = I + 1;
% sOut = [STUDY.cluster(I).sets];
% cOut = [STUDY.cluster(I).comps];
% ==== DEPRECTAED ==== %
%% CONSTRAIN
if CONSTRAIN_TO_ANATOMY
    for cnd = 1:length(SUBJSTRUCT)        
        subjStr = SUBJSTRUCT(cnd).subjStr;
        fprintf('Processing %s...\n',subjStr);
        SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME) = [];
        compsOut = cell(length(brainChars),1);
        for j = 1:length(brainChars)
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}) = BRAIN_STRUCT;
            [compsOut{j},~,~,~,AnatTable] = get_Dist2Anatomy(STUDY,ALLEEG,cnd,brainCoords{j},brainChars{j},subjStr,'threshold',ANAT_THRESH);
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}).AnatTable = AnatTable;
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}).AnatLoc = brainCoords{j};
        end
        components{cnd} = unique(compsOut{cnd});
        SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).components = components{cnd};
    end
else
    for cnd = 1:length(SUBJSTRUCT)        
        subjStr = SUBJSTRUCT(cnd).subjStr;
        fprintf('Processing %s...\n',subjStr);
        SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME) = [];
        compsOut = cell(length(brainChars),1);
        for j = 1:length(brainChars)
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}) = BRAIN_STRUCT;
            [compsOut{j},~,~,~,AnatTable] = get_Dist2Anatomy(STUDY,ALLEEG,cnd,brainCoords{j},brainChars{j},subjStr,'threshold',ANAT_THRESH);
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}).AnatTable = AnatTable;
            SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(brainChars{j}).AnatLoc = brainCoords{j};
        end
    end
    for cnd = 1:length(SUBJSTRUCT) %#ok<UNRCH>
        comps = SUBJSTRUCT(cnd).sourceLocalize.dipoleTable.ComponentNumber;
        rv = SUBJSTRUCT(cnd).sourceLocalize.dipoleTable.ResidualVariance;
        pBrain = SUBJSTRUCT(cnd).sourceLocalize.dipoleTable.ICLabel_pBrain;
        idx_rv = logical([rv{:}] <= RV_THRESH)';
        idx_pBrain = logical(pBrain >= BRAIN_THRESH);
        idx = idx_rv & idx_pBrain;
        components{cnd} = comps(idx);
        SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).components = components{cnd};
    end
end
%% SUBJECT REJECTION
keep = zeros(length(SUBJSTRUCT),1);
%- reject if they have too few components
for cnd = 1:length(SUBJSTRUCT)        
    if length(components{cnd})>COMP_THRESH
        keep(cnd) = 1;
    end        
end
disp(keep)
keep = logical(keep);
%- Reject subjects
if REJECT_SUBJECTS
    ALLEEG = ALLEEG(keep);
    SUBJSTRUCT = SUBJSTRUCT(keep);
    components = components(keep);
end

if all(keep == 0)
    error('ALL SUBJECTS REJETED. STOPPING ANLAYSIS.');
end
%% MAIN CONNECTIVITY PROCESSING (see. indvCnctAnl.m)
% generate connectivity measures
fprintf('\n==== CALCULATING CONNECTIVITY MEASURES ====\n')
parfor cnd = 1:length(SUBJSTRUCT)
    %## LOOP PARAMS
    SUBJ = SUBJSTRUCT(cnd);
    fPath = SUBJ.pathsEEG.cnctAnl(processNum).filepath;
    fName = SUBJ.pathsEEG.cnctAnl(processNum).filename;
    if DO_UNIX
        fPath = convertPath2UNIX(fPath,'dferris');
    else
        fPath = convertPath2Drive(fPath,'M');
    end
    subjStr = SUBJSTRUCT(cnd).subjStr;
    %## LOAD EEG
    EEG = pop_loadset('filename',ALLEEG(cnd).filename,'filepath', ALLEEG(cnd).filepath);
    %- Prints
    fprintf('%s) Processing componets:\n',subjStr)
    fprintf('%i,',components{cnd}'); fprintf('\n');
    %## GENERATE CONNECTIVITY MEASURES
    if ~exist([fPath filesep fName],'file') || REGEN_FILES
        tmpEEG = cnctanl_connMeas(EEG,PATHS,components{cnd},CONN_METHODS,fPath,...
            'WindowLengthSec',WINDOW_LENGTH,'WindowStepSizeSec', WINDOW_STEP_SIZE,...
            'NewSamplingRate',NEW_SAMPLE_RATE);
    else
        tmpEEG = par_load(fPath,fName);
    end
    %% (BOOTSTRAPPING) GROUP STATISTICS 
    % (09/22/2022) JS, Might want to try and speed up bootstrap by
    % addapting stat_surrogateGen.m to use parfor for bootstrapping... If
    % possible? doesn't seem built well in the first place, so maybe?
    if DO_BOOTSTRAP
        fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')        
        chk = strsplit(fName,'.');
        if ~exist([fPath filesep [chk{1}, '_BootStrap.mat']],'file') || REGEN_FILES
            tmpEEG = cnctanl_groupStats(tmpEEG,PATHS,'BootStrap');
            %- save BootStrap distribution
            tmp = tmpEEG.CAT.PConn;
            par_save(tmp,fPath,fName,[],'_BootStrap');
        else
            EEG.CAT.PConn = par_load(fPath,[chk{1}, '_BootStrap.mat'],[]);
            tmpEEG = EEG;
        end
        %- assign mean of bootstrap as Conn value
        if ASSIGN_BOOTSTRAP_MEAN
            tmpEEG.CAT.Conn = stat_getDistribMean(tmpEEG.CAT.PConn); %#ok<UNRCH>
        end
        %- clear bootstrap calculation
        tmpEEG.CAT.PConn = []        
    end
     % save EEG 
    par_save(tmpEEG,fPath,fName,[]);   
end
diary('off');
end

