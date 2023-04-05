function [SUBJSTRUCT] = subjstruct_genConnStats(PATHS,SUBJSTRUCT,processNum,studyName,studyDir,varargin)
%GENCONNSTATSSUBJS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

diary 'subjstruct_genConnStats'
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {true};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,PATHS,SUBJSTRUCT,processNum,studyName,studyDir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Define Defaults
%## PARAMS
DO_PHASE_RND = true;
DO_NONZERO_TEST = true;
RESAMPLE_DATA = false;
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
diaryPath = [PATHS.path4diaries filesep '4_diary_genConnStats.txt'];
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
    for j = 1:(length(fs)-1)
        clusts{cnd} = [ clusts{cnd,1} SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ClusterNumber'];
        comps{cnd} = [ comps{cnd,1} SUBJSTRUCT(cnd).cnctAnl.(ANALYSIS_NAME).(fs{j}).AnatTable.ComponentNumber'];
    end
    val = unique(comps{cnd});
    for i = 1:length(val)
        components{cnd} = [components{cnd} val(i)];
        tmp = unique(clusts{cnd}(val(i) == comps{cnd}));
        clusters{cnd} = [clusters{cnd}, {sprintf('%i',tmp)}];
    end
end

% for cnd = 1:length(components)
%     for j = 1:length(components{cnd})
%         idx = (components{cnd}(j) == comps);
%         vals = unique(clusts(idx));            
%         if length(vals) > 1
%             str = [sprintf('%i-',vals(1:end-1)) sprintf('%i',vals(end))];
%         else
%             str = sprintf('%i',vals);
%         end
%         clusters{cnd} = [clusters{cnd}, {str}];
%     end
% end
%% CONNECTIVITY
%%     
parfor cnd = 1:length(SUBJSTRUCT)
    %## PARAMS
    %## SETUP
    %## LOOP PRINTS
    fprintf('==== Processing %s ====', SUBJSTRUCT(cnd).subjStr);
    fPath = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum).filepath;
    fName = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum).filename;
    %## LOAD RESULTS
    EEG = par_load(fPath,fName,[]);
    %% GENERATE PHASE RANDOMIZED DISTRIBUTION
    fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
    % see. stat_surrogateGen
    % see. stat_surrogateStats
    %- Generate Phase Randomized Distribution
    if DO_PHASE_RND
        fPath = SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum).filepath;
        fName = [SUBJSTRUCT(cnd).pathsEEG.cnctAnl(processNum).filename '_PhaseRnd'];
        EEG.CAT.PConn = par_load(fPath,fName,[]);
    end
    %- Conduct Nonzero Test (Need Phase Randomized Distribution)
    if DO_NONZERO_TEST
        [tmpEEG,~] = cnctanl_groupStats(tmpEEG,PATHS,'NonZero');
    end
    %- save EEG
    par_save(tmpEEG,fPath,fName,[]);
    %- assign saveDir
    saveDir = [studyDir filesep EEG.subject];
    if exist(saveDir,'dir')
        mkdir(saveDir)
    end
    %% ==== COMPLEX VISUALS ==== %%    
    % TIMEFREQCHART (NO STATS)
    %## VISUALIZE TIME FREQ GRID (components to component)
    [~,~] = cnctanl_visualize(tmpEEG,PATHS,components{cnd},clusters{cnd},...
        'TimeFreqChart',...
        'connMeasures',SUBJSTRUCT(cnd).cnctAnl.connMethods,...
        'bounds',{},...
        'savePath',saveDir);

    %% EX-TIMEFREQCHART (STATS)   
    [~,~] = cnctanl_visualize(tmpEEG,PATHS,components{cnd},clusters{cnd},...
        'exTimeFreqChart',...
        'connMeasures',SUBJSTRUCT(cnd).cnctAnl.connMethods,...
        'bounds',{},...
        'savePath',saveDir);

    %% POP VIEW PROPS
    for i = 1:size(tmpEEG.icaact,1)
        hold on;
        pop_prop_extended(tmpEEG,0,i,NaN,{'freqrange',[1 65],'freqfac',4})
        fig = gcf;
        fig.Position = [500 300 920 480]; 
        hold off;
        savefig(fig,[saveDir sprintf('viewprops_%i.fig',i)]);
    end
end
end

