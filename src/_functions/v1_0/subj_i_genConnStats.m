function [SUBJ] = subj_i_genConnStats(PATHS,SUBJ,processNum,studyName,studyDir,varargin)
%GENCONNSTATSSUBJS Summary of this function goes here
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
addRequired(p,'SUBJ',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,PATHS,SUBJ,processNum,studyName,studyDir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Set Defaults
%## PARAMS
REGEN_FILES = true;
DO_PHASE_RND = true;
DO_NONZERO_TEST = true;
RESAMPLE_DATA = false;
FREQUENCY_BAND_TO_PLOT = [7:50];
ANALYSIS_NAME = SUBJ.META.epoched.trialType;
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
% (09/10/2022) was having error with ft_version (resulting from a call to
% ft_defaults) by hard coding the matlab search path the error seems to be
% fixed. No promises... (09/11/2022) Seems that I needed to also include
% some fixes to my batch pathing by setting as such:
% job = batch(pp,@subjstruct_genConnStats,1,{PATHS,SUBJSTRUCT,PROCESSNUM,studyName,studyDir},...
%                         'Pool',POOL_SIZE,...
%                         'AutoAttachFiles',false,...
%                         'AutoAddClientPath',false,...
%                         'CurrentFolder',PATHS.path4src,...
%                         'AttachedFiles',PATHS.path4functions,...
%                         'AdditionalPaths',struct2cell(PATHS));
%## DIARY
% diaryPath = [PATHS.path4diaries filesep studyName filesep '4_diary_genConnStats.txt'];
% if ~exist([PATHS.path4diaries filesep studyName],'dir')
%     mkdir([PATHS.path4diaries filesep studyName]);
% end
% diary(diaryPath)
%## EEGLAB
% eeglab;
%% ASSIGN LABELS FOR PLOTTING
%- condition/trial name

%- cells to store components and clusters for each subject
% components = cell(length(SUBJSTRUCT),1);
fs = fields(SUBJ.cnctAnl.(ANALYSIS_NAME));
% clusts = cell(length(SUBJSTRUCT),1);
% comps = cell(length(SUBJSTRUCT),1);
% clusters = cell(length(SUBJSTRUCT),1);
clusters = {};
%-
components = SUBJ.cnctAnl.(ANALYSIS_NAME).components;
fprintf('==== Creating Labels for Subject %s ====\n',SUBJ.subjStr);
for j = 1:(length(fs)-1)
    clusts = SUBJ.cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ClusterNumber';
    comps = SUBJ.cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ComponentNumber';
end
val = components;
comps = unique(comps);
for i = 1:length(val)
    tmp = clusts(val(i) == comps);
    clusters = [clusters, {sprintf('%i',tmp)}];
end
%- PRINT OUTS
fprintf('Plotting Components: '); fprintf('%i,',comps'); fprintf('\n');
fprintf('Associated Clusters: '); 
for i = 1:length(clusters)
    if isempty(clusters{i})
        clusters{i} = 'NaN';
    end
    fprintf('%s,',clusters{i}');
end
fprintf('\n');
%% CONNECTIVITY
TIME_LIMITS = SUBJ.META.epoched.epochTimeLimits;
%## PARAMS
timeLim = TIME_LIMITS;
%## SETUP
%## LOOP PRINTS
fprintf('==== Processing %s ====', SUBJ.subjStr);
fPath = SUBJ.pathsEEG.cnctAnl(processNum).filepath;
fName = SUBJ.pathsEEG.cnctAnl(processNum).filename;
if DO_UNIX
    fPath = convertPath2UNIX(fPath,'dferris');
else
    fPath = convertPath2Drive(fPath,'M');
end
%## LOAD RESULTS
EEG = par_load(fPath,fName,[]);
%- Resample if you like
if RESAMPLE_DATA
    NEW_HZ = 500;
    EEG = pop_resample(EEG, NEW_HZ, 0.8, 0.4); 
end
%% GENERATE PHASE RANDOMIZED DISTRIBUTION    
% see. stat_surrogateGen
% see. stat_surrogateStats
%- Generate Phase Randomized Distribution
if DO_PHASE_RND
    fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_PhaseRnd.mat']],'file') || REGEN_FILES
        [tmpEEG,~] = cnctanl_groupStats(EEG,PATHS,'PhaseRnd');
        %- Save Phase randomized distribution
        tmp = tmpEEG.CAT.PConn;
        par_save(tmp,fPath,fName,[],'_PhaseRnd');
    else
        EEG.CAT.PConn = par_load(fPath,[chk{1}, '_PhaseRnd.mat'],[]);
        tmpEEG = EEG;
    end
    fprintf('done.\n')
end
%- Conduct Nonzero Test (Need Phase Randomized Distribution)
if DO_NONZERO_TEST
    fprintf('\n==== NONZERO STATISTICS ====\n')
    chk = strsplit(fName,'.');
    if ~exist([fPath filesep [chk{1}, '_NonZero.mat']],'file') || REGEN_FILES
        [tmpEEG,~] = cnctanl_groupStats(tmpEEG,PATHS,'NonZero');
        %- Save Nonzero Stats
        tmp = tmpEEG.CAT.Stats;
        par_save(tmp,fPath,fName,[],'_NonZero');
    else
        tmpEEG.CAT.Stats = par_load(fPath,[chk{1}, '_NonZero.mat'],[]);
    end
    fprintf('done.\n')
end
%- assign saveDir
saveDir = [studyDir filesep EEG.subject];
if ~exist(saveDir,'dir')
    fprintf(1,'Making Directory: %s',saveDir);
    mkdir(saveDir)
end
%% ==== COMPLEX VISUALS ==== %%  
%## Generate Connectivity Matrix
CEstimator = 'dDTF08'; %SUBJ.cnctAnl.connMethods;
stat_alpha = 0.05;
catConn = tmpEEG.CAT.Conn.(CEstimator);
catConnMask = tmpEEG.CAT.Stats.(CEstimator).pval < stat_alpha;
%-
[SUBJ,meanMat,stdvMat] = gen_connMatrix(SUBJ,catConn,catConnMask,(7:40));
SUBJ.cnctAnl.connMatrix.comps = components;
SUBJ.cnctAnl.connMatrix.clust = clusts;
%-
I = eye(size(meanMat));
I = (I == 0);
meanMat = meanMat.*I;
mn = nanmean(meanMat(:));
stdv = nanmean(stdvMat(:));
clim = [mn-1*stdv,mn+1*stdv];
% TIMEFREQCHART (NO STATS)
%## VISUALIZE TIME FREQ GRID (components to component)
fprintf('%s) Plotting Time Frequency Plot...\n',SUBJ.subjStr)
[~,~] = cnctanl_visualize(tmpEEG,PATHS,components,clusters,...
    'TimeFreqChart',...
    'connMeasures',{CEstimator},...%SUBJ.cnctAnl.connMethods,...
    'colorLimit',clim,...
    'savePath',saveDir);
fprintf('%s) Done.',SUBJ.subjStr)
%% EX-TIMEFREQCHART (STATS)   
fprintf('%s) Plotting Time Frequency Plot Pval-Masked...\n',SUBJ.subjStr)
[~,~] = cnctanl_visualize(tmpEEG,PATHS,components,clusters,...
    'exTimeFreqChart',...
    'connMeasures',{CEstimator},...%SUBJ.cnctAnl.connMethods,...
    'colorLimit',clim,...
    'savePath',saveDir);
fprintf('%s) Done.',SUBJ.subjStr)
%% POP VIEW PROPS
fprintf('%s) Plotting Component Properties Plot...\n',SUBJ.subjStr)
for i = 1:size(tmpEEG.icaact,1)
    pop_prop_extended(tmpEEG,0,i,NaN,...
        {'freqrange',[1 65],'freqfac',4,'limits',[1,60,-30,0]},...
        {'erpstd','on','specaxis','log','limits',[timeLim(1)*1000,timeLim(2)*1000,-1,1]});
    fig = get(groot,'CurrentFigure');
    fig.Position = [15 15 1080 920];
    savefig(fig,[saveDir filesep sprintf('viewprops_%i.fig',i)]);
end
%- clear STATS & PCONN for space saving.
tmpEEG.CAT.Stats = [];
tmpEEG.CAT.PConn = [];
%- save EEG
par_save(tmpEEG,fPath,fName,[]);
end

