function [] = mim_precomp_ersp(STUDY,ALLEEG,warping_times,varargin)
%MIM_GEN_ERSP Summary of this function goes here
%   This is a CUSTOM function
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%- hard defines (for now 07/16/2023)
%- ERSP PARAMS
DO_PARALLEL = 'off';
ERSP_SUBBASELINE = 'on';
ERSP_TIMERANGE = [warping_times(1) warping_times(5)];
ERSP_LIM = [-2,2];
ERSP_FREQRANGE = [1,200];
%- soft defines
%## TIME
tic
%## INITIATE PARSER
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@islogical);
addRequired(p,'warping_times',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'DO_BASELINE_CORRECTION',DO_BASELINE_CORRECTION,@islogical);
addParameter(p,'DO_PARALLEL',DO_PARALLEL,@ischar
addParameter(p,'ERSP_SUBBASELINE',ERSP_SUBBASELINE,@ischar);
addParameter(p,'ERSP_TIMERANGE',ERSP_TIMERANGE,@isnumeric);
addParameter(p,'ERSP_LIM',ERSP_LIM,@isnumeric);
addParameter(p,'ERSP_FREQRANGE',ERSP_FREQRANGE,@isnumeric);
parse(p,STUDY,ALLEEG,warping_times,varargin{:});
%## SET DEFAULTS
%- PARAMETER
DO_BASELINE_CORRECTION = p.Results.DO_BASELINE_CORRECTION;
DO_PARALLEL = p.Results.DO_PARALLEL;
ERSP_SUBBASELINE = p.Results.ERSP_SUBBASELINE;
ERSP_TIMERANGE = p.Results.ERSP_TIMERANGE;
ERSP_LIM = p.Results.ERSP_LIM;
ERSP_FREQRANGE = p.Results.ERSP_FREQRANGE;

%% ===================================================================== %%
STUDY = pop_statparams(STUDY, 'condstats', STAT_PARAMS.condstats,...
        'method',STAT_PARAMS.method,...
        'singletrials',STAT_PARAMS.singletrials,'mode',STAT_PARAMS.mode,...
        'fieldtripalpha',STAT_PARAMS.fieldtripalpha,'fieldtripmethod',STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',STAT_PARAMS.fieldtripnaccu);

STUDY = pop_erspparams(tmpSTUDY, 'subbaseline',ERSP_PARAMS.subbaseline,...
        'timerange',ERSP_PARAMS.timerange, 'ersplim',ERSP_PARAMS.ersplim);  % 'subbaseline' - ['on'|'off'] subtract the same baseline across conditions for ERSP     
[~,main_cl_inds,~,valid_clusters] = eeglab_get_cluster_comps(STUDY);
CLUSTER_PICKS = main_cl_inds;
%##
parfor subj_i=1:length(ALLEEG)
    exitcode = 0;
    try
        EEG = ALLEEG(subj_i);
        tmp = STUDY;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        %- overrride datasetinfo to trick std_precomp to run.
        tmp.datasetinfo = STUDY.datasetinfo(subj_i);
        tmp.datasetinfo(1).index = 1;
        %- determine timewarping parameters
         if DO_TIMEWARP
            timewarp_param = EEG.timewarp.latencies;
            timewarpms_param = warping_times;
         else
             timewarp_param = [];
             timewarpms_param = [];
        end
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel',DO_PARALLEL,'cycles',CYCLE_LIMITS,...
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param,'baseline',b_lims,...
                    'commonbase','on','trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
            [~, ~] = std_precomp(tmp,EEG,'components','savetrials','on',...
                    'recompute','on','ersp','on','itc','off',...
                    'erspparams',{'parallel',DO_PARALLEL,'cycles',CYCLE_LIMITS,...
                    'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,...
                    'baseline',nan(),'timewarp',timewarp_param,...
                    'timewarpms',timewarpms_param}); %ERSP
        end
    catch e
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
        quit(exitcode);
    end
exitcode=1;
quit(exitcode)
end
[~, ersp_data, ersp_times, ersp_freqs, ~, params_ersp] = std_readdata(STUDY,ALLEEG,...
                            'clusters',cluster_i,'singletrials',ERSP_SINGLETRIALS,... 
                            'datatype','ersp','freqrange',ERSP_FREQLIMITS,...
                            'design',des_i);
%% ALTERNATIVE (not subject specific for timewarping...
%{
if DO_BASELINE_CORRECTION
    % Baseline correction
    [~, ~] = std_precomp(STUDY,ALLEEG,'components','savetrials','on',...
            'recompute','on','ersp','on','itc','off',...
            'erspparams',{'parallel','on','cycles',CYCLE_LIMITS,...
            'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
            'timewarpms',timewarpms_param,'baseline',b_lims,...
            'commonbase','on','trialbase','off','basenorm','on'}); %ERSP
else
    % No baseline correction
    [~, ~] = std_precomp(STUDY,ALLEEG,'components','savetrials','on',...
            'recompute','on','ersp','on','itc','off',...
            'erspparams',{'parallel','on','cycles',CYCLE_LIMITS,...
            'nfreqs',length((ERSP_FREQLIMITS(1):ERSP_FREQLIMITS(2))),'ntimesout',TIMEWARP_NTIMES,...
            'baseline',nan(),'timewarp',timewarp_param,...
            'timewarpms',timewarpms_param}); %ERSP
end
%}
end

