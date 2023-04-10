function [ALLEEG,STUDY] = eeglab_cluster(ALLEEG,studyName,studyDir,varargin)
%SUBJ_I_CLUSTER Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%## DEFINE DEFAULTS
DO_DIPOLE_REJECTION = true;
%- spec params
SPEC_MODE = 'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
LOG_TRIALS = 'off'; %options: 'on'
FREQ_FAC = 4;
SPEC_CONTINUOUS = 'off';
PAD_RATIO = 2;
SPEC_EPOCH_LIMS = [0,4];
%* 
FREQ_LIMITS = [1,100];
errorMsg = 'Value must be of format [INT1,INT2] where INT2 > INT1. Sets value limits for calculated spectral measure'; 
fl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%* 
CYCLE_LIMITS = [3,0.8];
errorMsg = 'Value must be of format [DOUBLE1,DOUBLE2]. Varies the length of the window used to determine the power of a certain frequency band';
cl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%*
DO_ERSP_CALC = false;
%- cluster params
SET_CLUSTER_N = 10; %8;
DO_VISUALIZATION = false;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'studyName',@ischar);
addRequired(p,'studyDir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'FREQ_LIMITS',FREQ_LIMITS,fl_validFcn); 
addParameter(p,'CYCLE_LIMITS',CYCLE_LIMITS,cl_validFcn); 
addParameter(p,'SPEC_MODE',SPEC_MODE,@ischar); 
addParameter(p,'LOG_TRIALS',LOG_TRIALS,@ischar); 
addParameter(p,'FREQ_FAC',FREQ_FAC,@isnumeric); 
addParameter(p,'SPEC_CONTINUOUS',SPEC_CONTINUOUS,@ischar); 
addParameter(p,'PAD_RATIO',PAD_RATIO,@isnumeric); 
addParameter(p,'SPEC_EPOCH_LIMS',SPEC_EPOCH_LIMS,@isnumeric); 
addParameter(p,'SET_CLUSTER_N',SET_CLUSTER_N,@isnumeric); 
addParameter(p,'DO_VISUALIZATION',DO_VISUALIZATION,@islogical); 
addParameter(p,'DO_ERSP_CALC',DO_ERSP_CALC,@islogical);
parse(p,ALLEEG,studyName,studyDir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
FREQ_LIMITS         = p.Results.FREQ_LIMITS;
CYCLE_LIMITS        = p.Results.CYCLE_LIMITS;
SPEC_MODE           = p.Results.SPEC_MODE;
LOG_TRIALS          = p.Results.LOG_TRIALS;
FREQ_FAC            = p.Results.FREQ_FAC;
SPEC_CONTINUOUS     = p.Results.SPEC_CONTINUOUS;
PAD_RATIO           = p.Results.PAD_RATIO;
SET_CLUSTER_N       = p.Results.SET_CLUSTER_N;
SPEC_EPOCH_LIMS     = p.Results.SPEC_EPOCH_LIMS;
%%
% DIP_NUM = 1;
% DIP_PLOT = 'off';
%## MAKE DIRS
if ~exist(studyDir,'dir')
    mkdir(studyDir);
end
%% ===================================================================== %%
%* empty STUDY structure for repopulating
STUDY = [];
%% CREATE STUDY
% initiailize study
fprintf('\n==== Making Study Modifications ====\n')
[STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,...
                                'name',studyName);
% make sure all .mat files have a .fdt file associated with it.
% (03/08/23) why were these turned off? <-- to save memory!! (03/16/2023),
% use eeglab_options to set memory options so it doesn't conflict.
% (08/28/22) updatedat turnned off 
% (08/28/22) savedat turned off
[STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,...
                                'updatedat','off',...
                                'savedat','off','resave','on',...
                                'commands', {'dipselect',0.5},...
                                'filename',studyName,...
                                'filepath',studyDir);
%     [STUDY, ALLEEG] = std_editset( STUDY, ALLEEG,...
%                                     'commands',{'dipselect',0.2,'inbrain','on'});

%% STUDY CONFIGURATION MODIFICATION
%Make STUDY design
% (04/13/2022): play around with these funcitons a bit, may be fun for
% quick permutations of conditions and sessions
% if ~isempty(ANALYSIS_COND)
%     fprintf('attaching conditions\n');
%     [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'condition',ANALYSIS_COND);
% end
% if ~isempty(ANALYSIS_SESS)
%     fprintf('attaching sessions\n');
%     [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'session',ANALYSIS_SESS);
% end
%% SAVE STUDY
for subj_i = 1:length(ALLEEG)
    %- check iclabel 
    if ~isfield(ALLEEG(subj_i).etc,'ic_classification')
       error('For some reason ICLABELS were not calculated for subject %s',ALLEEG(subj_i).subject)
    end
    %- check dipfit
    if ~isfield(ALLEEG(subj_i).dipfit,'model')
        disp(ALLEEG(subj_i).dipfit)
        error('For some reason DIPFITS were not calculated for subject %s',ALLEEG(subj_i).subject)
    end
end
if ~ispc
    [STUDY, ALLEEG] = pop_savestudy(STUDY,ALLEEG,'filename',[studyName '_UNIX'],'filepath',studyDir);
else
    [STUDY, ALLEEG] = pop_savestudy(STUDY,ALLEEG,'filename',studyName,'filepath',studyDir);
end

%% DIPOLE REJECTION'
ICLABEL_BRAIN_CRIT = 0.5;
RV_BRAIN_CRIT = 0.15;
if DO_DIPOLE_REJECTION
%     DipoleTable = cell(1,length(ALLEEG));
    comp_list_keep = [];
    set_list_keep  = [];
    comp_list_rej = [];
    set_list_rej  = [];
    for subj_i = 1:length(ALLEEG)
        %- check iclabel 
        goodICLabel = ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications(:,1)' >= ICLABEL_BRAIN_CRIT;
        %- check dipfit
        goodRV = [ALLEEG(subj_i).dipfit.model.rv] <= RV_BRAIN_CRIT;
        %- find independent components with a good iclabel and dipfit
        %residudal variance
        tmp_ic = find(goodRV & goodICLabel);
        tmp_rej = find(~(goodRV & goodICLabel));
        if length(tmp_ic)<2
           [val, sortedComps] = sort([ALLEEG(subj_i).dipfit.model.rv],'ascend');
           tmp_ic = sortedComps(1:2); 
        end
        comp_list_keep = [comp_list_keep, tmp_ic];
        set_list_keep  = [set_list_keep, repmat(subj_i,1,length(tmp_ic))];
        comp_list_rej = [comp_list_rej, tmp_rej];
        set_list_rej = [set_list_rej, repmat(subj_i,1,length(tmp_rej))];
%         STUDY.datasetinfod(subj_i).etc.ic_classification.ICLabel.classifications(:,1)))';
%         ICLclass = ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications;
%         rvOut = {ALLEEG(subj_i).dipfit.model.rv}';
%         t = table(compNum,rvOut, ICLclass(:,1),ICLclass(:,2),ICLclass(:,3),ICLclass(:,4),...
%             ICLclass(:,5),ICLclass(:,6),ICLclass(:,7),...
%             'VariableNames',{'ComponentNumber','ResidualVariance','ICLabel_pBrain','ICLabel_pMuscle',...
%             'ICLabel_pEye','ICLabel_pHeart','ICLabel_pLineNoise','ICLabel_pChannelNoise','ICLabel_pOther'});
         ALLEEG(subj_i).reject.gcompreject = ~(goodRV & goodICLabel);
         ALLEEG(subj_i).reject.rejmanual = 'yes';
         ALLEEG(subj_i) = pop_subcomp(ALLEEG(subj_i),tmp_rej,0,0);
    end
%     STUDY.preclust = [];
%     STUDY.cluster.sets = set_list_keep;
%     STUDY.cluster.comps = comp_list_keep;
    uS = unique(set_list_keep);
    suum = 0;
    for subj_i = 1:length(uS)    
        fprintf('** Subject %s has %i brain components\n',ALLEEG(subj_i).subject, length(comp_list_keep(set_list_keep == subj_i)));
        suum = suum + length(comp_list_keep(set_list_keep == subj_i));
    end
end
if isempty(SET_CLUSTER_N) && exist('uS','var')
    avgComps = ceil(suum/length(uS));
else
    avgComps = SET_CLUSTER_N;
end
%% PRECLUSTERING PROCESS
% make sure ALLEEG & STUDY structs are consistent
[STUDY, ALLEEG] = std_checkset(STUDY,ALLEEG);
% remove dipoles from analysis.
fprintf('\n');
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                'components',...
                                'rmicacomps','on');
% this computes power spectral density for designated study design
fprintf('==== Performing Precomputing ====\n');
if strcmp(SPEC_MODE,'psd')
%     SPEC_EPOCHLIM = [-1,1];
%     BOUNDARIES = (1:length(SPEC_EPOCHLIM):ALLEEG(1).xmax);
    BOUNDARIES = [];
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                    'components',...                               
                                    'recompute','on',...
                                    'spec','on',...
                                    'scalp','on',...
                                    'specparams',...
                                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                                    'freqrange',FREQ_LIMITS,'boundaries',BOUNDARIES});
else
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                    'components',...                               
                                    'recompute','on',...
                                    'spec','on',...
                                    'scalp','on',...
                                    'specparams',...
                                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                                    'freqrange',FREQ_LIMITS});
end
if DO_ERSP_CALC
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                'components',...
                                'ersp','on',...
                                'erspparams',...
                                {'cycles' CYCLE_LIMITS,'freqs',FREQ_LIMITS,...
                                'padratio',PAD_RATIO});
end
fprintf('done.\n');

%% PRECLUSTERING
fprintf('\n');
% CLUST METHOD 1
% [STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
%                        { 'spec', 'npca', 15, 'norm' 1 'weight' 1 } , ...
%                        { 'dipoles', 'norm' 1, 'weight' 3 });
% CLUST METHOD 2
[STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
                        { 'dipoles', 'norm' 1, 'weight' 5 },...
                        { 'scalp', 'norm' 1, 'weight' 5 });

% this uses the designated algorithm to cluster components
% Try: sweeping from 5 to 15 cluster numbers and see how they compare
% Try: using average number of components for all subjects.
fprintf('\n');
[STUDY]         = pop_clust(STUDY,ALLEEG,...
                    'algorithm','kmeans',...
                    'clus_num',avgComps); %,'outliers',1);

% STUDY = std_movecomp(STUDY,ALLEEG, from_cluster, to_cluster, comps);
%## KEEP TRACK OF REJECTED COMPS
STUDY.urcluster = STUDY.cluster;
STUDY.urcluster(end+1).name = sprintf('rvcrit%0.2f_iclabelcrit_Brain%i_rejected_comps',RV_BRAIN_CRIT,ICLABEL_BRAIN_CRIT);
STUDY.urcluster(end).parent = '';
STUDY.urcluster(end).sets = set_list_rej;
STUDY.urcluster(end).comps = comp_list_rej;
% check the STUDY & ALLEEG for consistency, and save.
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG); 
%% SAVE STUDY && CONVERT ALLEEG PATHS TO UNIX
%## ROBUST SAVE STUDY
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).etc.full_setfile.filename = ALLEEG(subj_i).filename;
    ALLEEG(subj_i).etc.full_setfile.filepath = ALLEEG(subj_i).filepath;
    ALLEEG(subj_i).filename = sprintf('%s_%s_ICA_TMPEEG',ALLEEG(subj_i).subject,'reducedcomps');
    ALLEEG(subj_i) = pop_saveset('filename',ALLEEG(subj_i).filename,'filepath',ALLEEG(subj_i).filepath);
end
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG); 
[STUDY,ALLEEG] = eeglab_save_study(STUDY,ALLEEG,...
                                    studyName,studyDir,...
                                    'RESAVE_DATASETS','on');
end
%% ===================================================================== %%
%## FUNCTIONS
