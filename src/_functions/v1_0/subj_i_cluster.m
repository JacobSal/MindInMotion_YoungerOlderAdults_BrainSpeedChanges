function [SUBJSTRUCT] = subj_i_cluster(PATHS,SUBJSTRUCT,processNum,studyName,...
                                        studyDir,varargin)
%CLUSTERSUBJS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {true,...
            []};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct);
addRequired(p,'SUBJSTRUCT',@isstruct);
addRequired(p,'processNum',@isnumeric);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
%## OPTIONAL
%## PARAMETER
veriLims = (@(x) (isnumeric(x) && length(x) == 2) || isempty(x));
veriEpoch = (@(x) (isnumeric(x) && length(x) == 2) || isempty(x));
addParameter(p,'freqLim',Defaults{2},veriLims); 
addParameter(p,'cycleLim',Defaults{2},veriLims); 
addParameter(p,'epochTimeLimits',Defaults{2},veriEpoch);
parse(p,PATHS,SUBJSTRUCT,processNum,studyName,studyDir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
freqLim         = p.Results.freqLim;
cycleLim        = p.Results.cycleLim;
epochTimeLimits = p.Results.epochTimeLimits;
%- Define Defaults
if isempty(epochTimeLimits)
    epochTimeLimits = SUBJSTRUCT(1).META.epoched.epochTimeLimits;
end
if isempty(freqLim)
    freqLim = [1,100]; % Power Spectrum [frequency min, frequency max] in (Hz)
end
if isempty(cycleLim)
    cycleLim = [0.8,3]; % Cycle limits for power spectrum calculation
end
%%
%## PARAMS
% precluster components for study
DO_DIPOLE_REJ = true;
CREATE_STUDY = true;
%- specparams
SPEC_MODE = 'psd'; %options: 'psd','fft','pburg','pmtm'
LOG_TRIALS = 'off'; %options: 'on'
FREQ_FAC = 4;
SPEC_CONTINUOUS = 'on';
ANALYSIS_SESS = 'first';
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
%## DIARY
diaryPath = [PATHS.path4diaries filesep studyName '2_diary_cluster.txt'];
if ~exist([PATHS.path4diaries filesep studyName],'dir')
    mkdir([PATHS.path4diaries filesep studyName]);
end
diary(diaryPath)
%## EEGLAB
eeglab;
%## FIELDTRIP
ft_defaults;
%% CREATE STUDY
disp('Creating study.')
if CREATE_STUDY
    ANALYSIS_COND = string(SUBJSTRUCT(1).META.epoched.trialType);
    ANALYSIS_GROUP = string(SUBJSTRUCT(1).SAVE.label);
    %* empty ALLEEG structure for repopulating
    ALLEEG = cell(1,length(SUBJSTRUCT)); 
    %* empty STUDY structure for repopulating
    STUDY = []; 
    %* conditions (e.g., ["rest","rest","rest","rest","rest"])
    conds = repmat(string(ANALYSIS_COND),1,length(SUBJSTRUCT));
    %* sessions (e.g., [1,1,1,1,1])
    sessions = repmat(1,1,length(SUBJSTRUCT));
    %* group (e.g., ["1","1","1","1","1"])
    group = repmat(ANALYSIS_GROUP,1,length(SUBJSTRUCT));
    fsPrev = {};
    % populate ALLEEG struct
    for cnd=1:length(SUBJSTRUCT)
        fPath = SUBJSTRUCT(cnd).pathsEEG.epoched(processNum).filepath;
        fName = SUBJSTRUCT(cnd).pathsEEG.epoched(processNum).filename;
        subjStr = SUBJSTRUCT(cnd).subjStr;
        if DO_UNIX
            fPath = convertPath2UNIX(fPath,'dferris');
        else
            fPath = convertPath2Drive(fPath,'M');
        end
        fprintf(1,'Loading Subject %s\n',subjStr)
        EEG = pop_loadset('filepath',fPath,'filename',fName);
        EEG.filepath = fPath;
        EEG.filename = fName;
        EEG.subject = subjStr;
        EEG.group = char(group(cnd));
        EEG.condition = char(conds(cnd));
        EEG.session = sessions(cnd);
        fs = fields(EEG);
        % delete fields not present in other structs.
        out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); out = [out{:}];
        delFs = fs(~out);
        if ~isempty(fsPrev) && any(~out)
            for j = 1:length(delFs)
                EEG = rmfield(EEG,delFs{j});
                fprintf("%s) Removing fields %s",EEG.subject,delFs{j})
            end
        else
            fsPrev = fs;
        end 
        ALLEEG{cnd} = EEG;          
    end
    %- CONCATENATE ALLEEG
    val = [];
    ALLEEG = cellfun(@(x) [val; x], ALLEEG);
    %% CREATE STUDY
    % initiailize study
    fprintf('Initializing STUDY.\n')
    [STUDY, ALLEEG] = std_editset( STUDY, ALLEEG,...
                                    'name',studyName);
    % make sure all .mat files have a .fdt file associated with it.
    % (08/28/22) updatedat turnned off 
    % (08/28/22) savedat turned off
    [STUDY, ALLEEG] = std_editset( STUDY, ALLEEG,...
                                    'updatedat','off','savedat','off','resave','off',...
                                    'filename',studyName,...
                                    'filepath',studyDir);
%     [STUDY, ALLEEG] = std_editset( STUDY, ALLEEG,...
%                                     'commands',{'dipselect',0.2,'inbrain','on'});

    %% STUDY CONFIGURATION MODIFICATION
    %Make STUDY design
    % (04/13/2022): play around with these funcitons a bit, may be fun for
    % quick permutations of conditions and sessions
    if ~isempty(ANALYSIS_COND)
        fprintf('attaching conditions\n');
        [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'condition',ANALYSIS_COND);
    end
    if ~isempty(ANALYSIS_SESS)
        fprintf('attaching sessions\n');
        [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'session',ANALYSIS_SESS);
    end
    %% SAVE STUDY
    if DO_UNIX
        [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',[studyName '_UNIX'],'filepath',studyDir);
    else
        [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',studyName,'filepath',studyDir);
    end
else
    %Load STUDY (assuming it exists)
    if DO_UNIX
        [STUDY, ALLEEG] = pop_loadstudy('filename',[studyName,'_UNIX','.study'],'filepath',studyDir);        
    else
        [STUDY, ALLEEG] = pop_loadstudy('filename',[studyName,'.study'],'filepath',studyDir); 
    end
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = (1:length(ALLEEG));
end

%% DIPOLE REJECTION
if DO_DIPOLE_REJ
    compList = [];
    setList  = [];
    for subj_i = 1:length(SUBJSTRUCT)
        if ~isfield(ALLEEG(subj_i).etc,'ic_classification')
           error('please run ICLabel and write down the code for it');
        else
           goodICLabel = ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications(:,1)' >= 0.5;
        end
        goodRV = [ALLEEG(subj_i).dipfit.model.rv] <= 0.15;
        tempComps = find(goodRV & goodICLabel);
        if length(tempComps)<2
           [val, sortedComps] = sort([ALLEEG(subj_i).dipfit.model.rv],'ascend');
           tempComps = sortedComps(1:2); 
        end
        compList = [compList, tempComps];
        setList  = [setList, repmat(subj_i,1,length(tempComps))];
        STUDY.datasetinfo(subj_i).comps = tempComps;
        fprintf('%s) Removing components: ',SUBJSTRUCT(subj_i).subjStr); fprintf('%i, ', find(~(goodRV & goodICLabel))); fprintf('\n');
        % save results
        compNum = (1:length(ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications(:,1)))';
        ICLclass = ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications;
        rvOut = {ALLEEG(subj_i).dipfit.model.rv}';
        t = table(compNum,rvOut, ICLclass(:,1),ICLclass(:,2),ICLclass(:,3),ICLclass(:,4),...
            ICLclass(:,5),ICLclass(:,6),ICLclass(:,7),...
            'VariableNames',{'ComponentNumber','ResidualVariance','ICLabel_pBrain','ICLabel_pMuscle',...
            'ICLabel_pEye','ICLabel_pHeart','ICLabel_pLineNoise','ICLabel_pChannelNoise','ICLabel_pOther'});
        SUBJSTRUCT(subj_i).sourceLocalize.dipoleTable = t;
    end
    STUDY.preclust = [];
    STUDY.cluster.sets = setList;
    STUDY.cluster.comps = compList;
    uS = unique(setList);
    suum = 0;
    for cnd = 1:length(uS)    
        fprintf('** Subject %s has %i brain components\n',SUBJSTRUCT(cnd).subjStr, length(compList(setList == cnd)));
        suum = suum + length(compList(setList == cnd));
        SUBJSTRUCT(subj_i).sourceLocalize.study.components = compList(setList == cnd);
    end
    avgComps = ceil(suum/length(uS));
else
    avgComps = 12;
end
%% PRECLUSTERING PROCESS
if true
    % make sure ALLEEG & STUDY structs are consistent
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
    % remove dipoles from analysis.
    fprintf('\n');
    % [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG,...
    %                                 'components',...
    %                                 'rmicacomps','on');
    % this computes power spectral density for designated study design
    fprintf('\n');
    [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG,...
                                    'components',...                               
                                    'recompute','on',...
                                    'spec','on',...
                                    'erp','on',...
                                    'scalp','on',...
                                    'specparams',...
                                    {'continuous',SPEC_CONTINUOUS,'specmode',SPEC_MODE,'logtrials',LOG_TRIALS,'freqfac',FREQ_FAC,'freqrange',freqLim});
    cyc1 = 3;
    cyc2 = 0.8;
    [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG,...
                                'components',...
                                'ersp','on',...
                                'erspparams',...
                                {'cycles' [cyc1 cyc2],'freqs',freqLim,...
                                'padratio',2});


    %% PRECLUSTERING
    fprintf('\n');
    % CLUST METHOD 1
    % [STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
    %                        { 'spec', 'npca', 15, 'norm' 1 'weight' 1 } , ...
    %                        { 'dipoles', 'norm' 1, 'weight' 3 });
    % CLUST METHOD 2
    [STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
                            { 'dipoles', 'norm' 1, 'weight' 3 });

    % this uses the designated algorithm to cluster components
    % Try: sweeping from 5 to 15 cluster numbers and see how they compare
    % Try: using average number of components for all subjects.
    fprintf('\n');
    [STUDY]         = pop_clust(STUDY, ALLEEG,...
                        'algorithm','kmeans',...
                        'clus_num',  avgComps); % , 'outliers',  1 );
    % check the STUDY & ALLEEG for consistency, and save.
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); 
    if DO_UNIX
        [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',[studyName '_UNIX'],'filepath',studyDir);
    else
        [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',studyName,'filepath',studyDir);
    end
    % STUDY = pop_clustedit(STUDY,ALLEEG) %## DEBUG/GUI
    %% SUBOUT CUSTOM MRI FOR MNI
    % fprintf('Reassigning MRI for MNI plotting...\n');
    % % pop_clustedit(STUDY, ALLEEG, clusts); %## DEBUG
    % path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
    % mniMRI = fullfile(path2BEM, 'standard_mri.mat');
    % for i = 1:length(SUBJSTRUCT)
    %      ALLEEG(i).dipfit.mrifile = mniMRI;
    % end
    % mri = ft_read_mri(mniMRI);
    % disp('done')
    %% VISUALIZATION (part 1)
    if ~SUBJSTRUCT(1).META.batch.clusterVisualizationDone
        clusts = (1:length(STUDY.cluster))'; % remove the parent cluster?
        clustMax = length(STUDY.cluster);
        specMin = -50;
        specMax = 1;
        % look into "no component indices, making incremental ones"
        %     STUDY = pop_dipparams(STUDY, 'axistight', 'off',...
        %                             'projimg', 'off',...
        %                             'projlines', 'off',...
        %                             'density', 'off',...
        %                             'centrline', 'on');
        STUDY = pop_specparams(STUDY,'freqrange',[1,60]);
        fig = cell(length(clustMax),1);
        for clus = 1:clustMax
            fprintf('Plotting cluster %i\n',clus)
            %## DIPOLE PLOT
            fig{clus} = figure();
            std_dipplot(STUDY,ALLEEG,'clusters',clus,'figure','off');
            fig{clus}.Position = [100 100 920 720];
            view([45,-45,45])
            saveas(fig{clus},[studyDir filesep sprintf('DipPlot_%i.fig',clus)]);
            saveas(fig{clus},[studyDir filesep sprintf('DipPlot_%i.jpg',clus)]);
        end    
        for clus = 1:clustMax
            %## ERP PLOT
            std_erpplot(STUDY,ALLEEG,'clusters',clus,...
                                'ylim',[-0.5,0.5]);
            fig_i = get(groot,'CurrentFigure');
            fig_i.Position = [15 15 1080 750];        
            saveas(fig_i,[studyDir filesep sprintf('ErpPlot_%i.fig',clus)]);
            saveas(fig_i,[studyDir filesep sprintf('ErpPlot_%i.jpg',clus)]);
        end
        for clus = 1:clustMax
            %## SPEC PLOT
            STUDY = std_specplot(STUDY,ALLEEG,'clusters',clus, 'ylim',[specMin,specMax]);
            fig_i = get(groot,'CurrentFigure');
            fig_i.Position = [15 15 1080 750];
            saveas(fig_i,[studyDir filesep sprintf('SpecPlot_%i.fig',clus)]);
            saveas(fig_i,[studyDir filesep sprintf('SpecPlot_%i.jpg',clus)]);
        end
        for clus = 1:clustMax
            %## TOPO PLOT
            STUDY = std_topoplot(STUDY,ALLEEG,'clusters', clus);
            fig_i = get(groot,'CurrentFigure');
            fig_i.Position = [15 15 1080 750];
            saveas(fig_i,[studyDir filesep sprintf('TopoPlot_%i.fig',clus)]);
            saveas(fig_i,[studyDir filesep sprintf('TopoPlot_%i.jpg',clus)]);
        end
        %'dipcolor',{'red','green','blue','purple','orange'})
        %% VISUALIZATION (part 2)
        % see. pop_erspparams()
        % STUDY = pop_erspparams(STUDY,...
        %                         'freqrange',[1,40],'ersplim',[-2,2],...
        %                         'averagemode','rms','maskdata','off',...
        %                         'averagechan','off','subbaseline','off');

        % PLOTS
        %## DIPOLE PLOT
        [STUDY] = std_dipplot(STUDY,ALLEEG,'clusters',clusts,'figure','off');
        fig_j = get(groot,'CurrentFigure');
        fig_j.Position = [100 100 920 720];
        view([45,-45,45])
        saveas(fig_j,[studyDir filesep sprintf('allDipPlot.fig')]);
        saveas(fig_j,[studyDir filesep sprintf('allDipPlot.jpg')]);
        %## ERP PLOT
        STUDY = std_erpplot(STUDY,ALLEEG,'clusters',clusts,...
                            'mode', 'together',...
                            'plotmode', 'condensed',...
                            'ylim',[-0.5,0.5]);
        fig_i = get(groot,'CurrentFigure');
        fig_i.Position = [15 15 1080 750];    
        saveas(fig_i,[studyDir filesep sprintf('allErpPlot.fig')]);
        saveas(fig_i,[studyDir filesep sprintf('allErpPlot.jpg')]);
        %## SPEC PLOT
        STUDY = std_specplot(STUDY,ALLEEG,'clusters',clusts, 'ylim',[specMin,specMax]);
        fig_i = get(groot,'CurrentFigure');
        fig_i.Position = [15 15 1080 750];
        saveas(fig_i,[studyDir filesep sprintf('allSpecPlot.fig')]);
        saveas(fig_i,[studyDir filesep sprintf('allSpecPlot.jpg')]);
        %## TOPO PLOT
        STUDY = std_topoplot(STUDY,ALLEEG,'clusters', clusts);
        fig_i = get(groot,'CurrentFigure');
        fig_i.Position = [15 15 1080 750];
        saveas(fig_i,[studyDir filesep sprintf('allTopoPlot.fig')]);
        saveas(fig_i,[studyDir filesep sprintf('allTopoPlot.jpg')]);
        %## ERSP PLOT
        % STUDY = std_erspplot(STUDY,ALLEEG,'clusters',clusts);
        % saveas(gcf,[STUDY_DIR filesep sprintf('allErspPlot.fig')]);
        % saveas(gcf,[STUDY_DIR filesep sprintf('allErspPlot.jpg')]);
        SUBJSTRUCT(1).META.batch.clusterVisualizationDone = true;
    end
end
%% SAVE STUDY
if DO_UNIX
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',[studyName '_UNIX'],'filepath',studyDir);
else
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename',studyName,'filepath',studyDir);
end
%% CONVERT ALLEEG PATHS TO UNIX!
% STUDY.filepath = convertPath2UNIX(STUDY.filepath,'dferris');
tmpSTUDY = STUDY;
if DO_UNIX
    for cnd = 1:length({STUDY.datasetinfo.filepath})
        tmpSTUDY.datasetinfo(cnd).filepath = convertPath2Drive(tmpSTUDY.datasetinfo(cnd).filepath,'M');
    end
    [~]         = pop_savestudy( tmpSTUDY, ALLEEG,...
                        'filename',studyName,'filepath',studyDir);
else
    for cnd = 1:length({STUDY.datasetinfo.filepath})
        tmpSTUDY.datasetinfo(cnd).filepath = convertPath2UNIX(tmpSTUDY.datasetinfo(cnd).filepath,'dferris');
    end
    [~]         = pop_savestudy( tmpSTUDY, ALLEEG,...
                        'filename',[studyName '_UNIX'],'filepath',studyDir);
end
diary('off');
end
%% ===================================================================== %%
%## FUNCTIONS
