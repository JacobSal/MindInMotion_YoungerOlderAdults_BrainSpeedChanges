function [SUBJSTRUCT] = subjstruct_reject(PATHS,SUBJSTRUCT,processNum,studyName,studyDir,connMethods,brainChars,brainCoords,varargin)
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
RV_THRESH = 0.15;
BRAIN_THRESH = 0.5;
CONSTRAIN_TO_ANATOMY = false;
REJECT_SUBJECTS = true;
BRAIN_STRUCT = struct('AnatTable',[]);
components = cell(length(SUBJSTRUCT),1);
ANALYSIS_NAME = SUBJSTRUCT(1).META.epoched.trialType;
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
% diaryPath = [PATHS.path4diaries filesep studyName filesep '3_diary_genConnMeas.txt'];
% if ~exist([PATHS.path4diaries filesep studyName] ,'dir')
%     mkdir([PATHS.path4diaries filesep studyName]);
% end
% diary(diaryPath)
%## EEGLAB
% eeglab;
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
    res = cell(1,length(SUBJSTRUCT));
    parfor subj_i = 1:length(SUBJSTRUCT)
        SUBJ = SUBJSTRUCT(subj_i);
        res{subj_i} = reject_by_iclabel(SUBJ,brainChars,brainCoords,STUDY,ALLEEG,...
                                    subj_i,ANALYSIS_NAME,ANAT_THRESH,RV_THRESH,...
                                    BRAIN_THRESH,BRAIN_STRUCT);
    end
    for subj_i = 1:length(SUBJSTRUCT)
        SUBJSTRUCT(subj_i) = res{subj_i};
        components{subj_i} = SUBJSTRUCT(subj_i).cnctAnl.(ANALYSIS_NAME).components;
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
    SUBJSTRUCT = SUBJSTRUCT(keep);
end

if all(keep == 0)
    error('ALL SUBJECTS REJETED. STOPPING ANLAYSIS.');
end
% SUBJSTRUCT = subj_saveStruct(SUBJSTRUCT);
% STUDY = pop_savestudy(STUDY,ALLEEG);
end

function SUBJ = reject_by_iclabel(SUBJ,brainChars,brainCoords,STUDY,ALLEEG,...
                                    INDEX,ANALYSIS_NAME,ANAT_THRESH,RV_THRESH,...
                                    BRAIN_THRESH,BRAIN_STRUCT)
    subjStr = SUBJ.subjStr;
    fprintf('Processing %s...\n',subjStr);
    
    SUBJ.cnctAnl.(ANALYSIS_NAME) = [];
    compsOut = cell(1,length(brainChars));
    for j = 1:length(brainChars)
        BCoord = brainCoords{j};
        BChar = brainChars{j};
        [SUBJ,compsOut{j}] = anatomy_calc(SUBJ,BCoord,BChar,STUDY,ALLEEG,INDEX,ANALYSIS_NAME,ANAT_THRESH,BRAIN_STRUCT);
    end
    comps = SUBJ.sourceLocalize.dipoleTable.ComponentNumber;
    rv = SUBJ.sourceLocalize.dipoleTable.ResidualVariance;
    pBrain = SUBJ.sourceLocalize.dipoleTable.ICLabel_pBrain;
    idx_rv = logical([rv{:}] <= RV_THRESH)';
    idx_pBrain = logical(pBrain >= BRAIN_THRESH);
    idx = idx_rv & idx_pBrain;
    components = comps(idx);
    SUBJ.cnctAnl.(ANALYSIS_NAME).components = components;
end

function [SUBJ,compsOut] = anatomy_calc(SUBJ,BCoord,BChar,STUDY,ALLEEG,INDEX,ANALYSIS_NAME,ANAT_THRESH,BRAIN_STRUCT)
    SUBJ.cnctAnl.(ANALYSIS_NAME).(BChar) = BRAIN_STRUCT;
    [compsOut,~,~,~,AnatTable] = get_Dist2Anatomy(STUDY,ALLEEG,INDEX,BCoord,...
                                        BChar,sprintf('%i',INDEX),'threshold',ANAT_THRESH);
    SUBJ.cnctAnl.(ANALYSIS_NAME).(BChar).AnatTable = AnatTable;
    SUBJ.cnctAnl.(ANALYSIS_NAME).(BChar).AnatLoc = BCoord;
end