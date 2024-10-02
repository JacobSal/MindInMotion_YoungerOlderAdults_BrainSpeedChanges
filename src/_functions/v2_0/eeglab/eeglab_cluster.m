function [STUDY,ALLEEG] = eeglab_cluster(STUDY,ALLEEG,varargin)
%SUBJ_I_CLUSTER Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu,
%## TIME
tic
%## DEFINE DEFAULTS
%- switches
DO_ERSP_CALC = false;
DO_SPEC_CALC = true;
%- spec params
SPEC_MODE = 'psd'; %'psd'; %options: 'psd','fft','pburg','pmtm'
LOG_TRIALS = 'on'; %options: 'on'
FREQ_FAC = 4;
PAD_RATIO = 2;
%* 
FREQ_LIMITS = [1,100];
errorMsg = 'Value must be of format [INT1,INT2] where INT2 > INT1. Sets value limits for calculated spectral measure'; 
fl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%* 
CYCLE_LIMITS = [3,0.8];
errorMsg = 'Value must be of format [DOUBLE1,DOUBLE2]. Varies the length of the window used to determine the power of a certain frequency band';
cl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%- cluster params
% (06/14/2023) JS, Trying K=9 for Older Adults
KMEANS_CLUSTER_N = 9; %10; %8; 
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'FREQ_LIMITS',FREQ_LIMITS,fl_validFcn); 
addParameter(p,'CYCLE_LIMITS',CYCLE_LIMITS,cl_validFcn); 
addParameter(p,'SPEC_MODE',SPEC_MODE,@ischar); 
addParameter(p,'LOG_TRIALS',LOG_TRIALS,@ischar); 
addParameter(p,'FREQ_FAC',FREQ_FAC,@isnumeric); 
addParameter(p,'PAD_RATIO',PAD_RATIO,@isnumeric);
addParameter(p,'KMEANS_CLUSTER_N',KMEANS_CLUSTER_N,@isnumeric); 
addParameter(p,'DO_ERSP_CALC',DO_ERSP_CALC,@islogical);
addParameter(p,'DO_SPEC_CALC',DO_SPEC_CALC,@islogical); 
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
FREQ_LIMITS         = p.Results.FREQ_LIMITS;
CYCLE_LIMITS        = p.Results.CYCLE_LIMITS; 
SPEC_MODE           = p.Results.SPEC_MODE;
LOG_TRIALS          = p.Results.LOG_TRIALS; 
FREQ_FAC            = p.Results.FREQ_FAC;
PAD_RATIO           = p.Results.PAD_RATIO;
KMEANS_CLUSTER_N    = p.Results.KMEANS_CLUSTER_N;
DO_SPEC_CALC        = p.Results.DO_SPEC_CALC;
DO_ERSP_CALC        = p.Results.DO_ERSP_CALC;
%% ===================================================================== %%
%## SAVE STUDY
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
%%
% make sure ALLEEG & STUDY structs are consistent
[STUDY, ALLEEG] = std_checkset(STUDY,ALLEEG);
% remove dipoles from analysis.
fprintf('\n');
% this computes power spectral density for designated study design
if DO_SPEC_CALC
    fprintf('==== Performing PSD & Topo Precomputing ====\n');
    parfor subj_i = 1:length(ALLEEG)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        [~, ~] = std_precomp(tmp, EEG,...
                    'components',...
                    'recompute','on',...
                    'spec','on',...
                    'scalp','on',...
                    'savetrials','on',...
                    'specparams',...
                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                    'freqrange',FREQ_LIMITS,'logtrials',LOG_TRIALS});
    end
end
%%
if DO_ERSP_CALC
    fprintf('==== Performing ERSP Precomputing ====\n');
    for subj_i = 1:length(ALLEEG)
        ALLEEG(subj_i).group = '1';
        STUDY.datasetinfo(subj_i).group = '1';
    end
    for subj_i = 1:length(ALLEEG)
        %- extract per condition timewarp from event struct
    %     new_warp = struct('latencies',[],'epochs',[],'eventSequence',[],'warpto',[]);
    %     new_warp.latencies = cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.latencies);
    %     new_warp.epochs = 1:length(new_warp.latencies);
    %     new_warp.eventSequence = ALLEEG(subj_i).etc.timewarp_by_cond(1).eventSequence;
    %     new_warp.warpto = ALLEEG(subj_i).timewarp.warpto; %median(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto)); % could use mean here, using median for now
        all_warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto)); % could use mean here, using median for now
        %- reassign per condition timewarping
        ALLEEG(subj_i).timewarp.warpto = all_warpto;
    end
    % allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    all_warpto = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
    for subj_i = 1:length(ALLEEG)
        all_warpto(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
    end
    % grandAvgWarpTo = median(allWarpTo);
    grandAvgWarpTo = nanmedian(all_warpto);
    disp(['Grand average (across all subj) warp to: ',num2str(grandAvgWarpTo)]);
    parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        end
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        %- determine timewarping parameters
         if DO_TIMEWARP
            timewarp_param = EEG.timewarp.latencies;
            timewarpms_param = grandAvgWarpTo;
         else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length(FREQ_LIMITS),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',[grandAvgWarpTo(1) grandAvgWarpTo(end)],...
                    'commonbase','on','trialbase','off','basenorm','on','padratio',PAD_RATIO}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel','off','cycles',CYCLE_LIMITS,...
                    'nfreqs',length(FREQ_LIMITS),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'padratio',PAD_RATIO}); %ERSP
        end
    end
end
fprintf('done.\n');
%% PRECLUSTERING
fprintf('Clustering...\n');
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                'components');
% make sure ALLEEG & STUDY structs are consistent
[STUDY, ALLEEG] = std_checkset(STUDY,ALLEEG);
fprintf('\n');
% CLUST METHOD (MIM 04/23/2023)
[STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
                        { 'dipoles', 'norm' 1, 'weight' 5 },...
                        { 'scalp', 'norm' 1, 'weight' 5 });
% [STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
%                         { 'dipoles', 'norm' 1, 'weight' 5 });

% this uses the designated algorithm to cluster components
% Try: sweeping from 5 to 15 cluster numbers and see how they compare
% Try: using average number of components for all subjects.
fprintf('\n');
[STUDY]         = pop_clust(STUDY,ALLEEG,...
                    'algorithm','kmeans',...
                    'clus_num',KMEANS_CLUSTER_N); %,'outliers',1);

%% SAVE STUDY && CONVERT ALLEEG PATHS TO UNIX
STUDY.urcluster = STUDY.cluster;
% check the STUDY & ALLEEG for consistency, and save.
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
end

function custom_cluster(STUDY,ALLEEG,kmeans_cluster_n)

end
