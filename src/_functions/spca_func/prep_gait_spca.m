function [gait_avg,ERDS,GPM,TF,output_struct] = prep_gait_spca(EEG_gait,base_txf_mean,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here

WAVELET_STRUCT = struct('t',[0,1/EEG_gait.srate],...
    'f',(4:100),...
    'fc',1,...
    'FWHM_tc',3,...
    'squared','n',...
    'data_type','double');
SPCA_PARAMS = struct('analysis_type','component',...
    'event_char','RHS',...
    'epoch_min_max',[1,4.25],...
    'n_resamples',100,...
    'timewarp_events',{{'RHS','LHS','LTO','RTO'}},...
    'condition_base','rest',...
    'condition_gait',{{''}});
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG_gait',@isstruct);
addRequired(p,'base_txf_data',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'SPCA_PARAMS',SPCA_PARAMS,@(x) validate_struct(x,SPCA_PARAMS));
addParameter(p,'WAVELET_STRUCT',WAVELET_STRUCT,@(x) validate_struct(x,WAVELET_STRUCT));
parse(p,EEG_gait,base_txf_mean,varargin{:});
%## SET DEFAULTS
SPCA_PARAMS = p.Results.SPCA_PARAMS;
WAVELET_STRUCT = p.Results.WAVELET_STRUCT;

%% ===================================================================== %
[TF,noise_cov,params] = eeg_txf_decomp(EEG_gait,SPCA_PARAMS.analysis_type,'WAVELET_STRUCT',WAVELET_STRUCT);
%- take magnitude (not power), pnts x chans x freqs
TF = abs(TF);
%% MARTIN S RESAMPLING CODE (ACTS AS TIMEWARPING) ====================== %%
%- ASSIGNED VALUES
n_freqs = size(TF,3);
n_comps = size(TF,2);
hs_min_max = SPCA_PARAMS.epoch_min_max*EEG_gait.srate; % time of next RHS in s
idx_hs = find(strcmp({EEG_gait.event.type}, SPCA_PARAMS.event_char));
gait_tf = zeros(length(idx_hs)-1,SPCA_PARAMS.n_resamples,n_comps*n_freqs,WAVELET_STRUCT.data_type); %strides/trials x pnts x chans x freqs
%- step counter, increased for each valid step
cnt = 1; 
%- resample each stride to the same legth (100 pnts)
for cycle_cnt = 1:length(idx_hs)-1
    %- find first and last sample of stride
    cycle_edge = round([EEG_gait.event(idx_hs(cycle_cnt)).latency,...
        EEG_gait.event(idx_hs(cycle_cnt+1)).latency-1]); % first and last frame of gait cycle
    %- labels of all events within this cycle
    cycle_event = {EEG_gait.event([idx_hs(cycle_cnt):idx_hs(cycle_cnt+1)]).type};
    %- only keep labels of gait events to check their order:
    cycle_gaitEvent = cycle_event(contains(cycle_event,SPCA_PARAMS.timewarp_events));
    %-
    if hs_min_max(1) <= cycle_edge(2)-cycle_edge(1) &&... % check time until next HS
            cycle_edge(2)-cycle_edge(1) <= hs_min_max(2) &&...
            all(ismember(SPCA_PARAMS.timewarp_events,cycle_gaitEvent)) %&& ...% oder of gait events correct
%             all(EEG_gait.etc.valid_eeg(cycle_edge(1):cycle_edge(2))) % no high amplitude samples
        %- 
        tf_cycle = TF(cycle_edge(1):cycle_edge(2),:,:); % extract data
        tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*n_freqs); % reshape to be able to use the resample function, skip but resample over different dimension?
        gait_tf(cnt,:,:) = resample(tf_cycle,SPCA_PARAMS.n_resamples,cycle_edge(2)-cycle_edge(1)+1,0); % resample and store
        cnt = cnt+1;
    end
end
disp([num2str(round(cnt/cycle_cnt*100)) '% of the gait cycles are valid'])
%- Gait_TF now: strides/trials x pnts + chans x freqs
gait_tf = reshape(gait_tf,size(gait_tf,1),SPCA_PARAMS.n_resamples,n_comps,n_freqs); % reshape to trials x pnts x chans x freqs
gait_avg = squeeze(mean(gait_tf)); % average over trials
ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_txf_mean));
%## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
GPM = bsxfun(@minus,ERDS,mean(ERDS));
output_struct = struct('cycle_cnt',cycle_cnt,'valid_cycle_cnt',cnt,'baseline_cov',noise_cov,'morlet_params',params);
%% (ALTERNATIVE) USE TIMEWARPPING EEGLAB =================================== %%
% functions: newtimef.m, timefreq.m, timewarp.m
% (12/9/2023) JS, the timewarp function is the real meat that is needed to
% warp gait events to a grand average. Could use the morlet_transform_fast
% after that. 
end
%% ===================================================================== %%
function [b] = validate_struct(x,DEFAULT_STRUCT)
    b = false;
    struct_name = inputname(2);
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field names
    chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
    if ~all(chk)
        fprintf(2,'\nFields for struct do not match for %s\n',struct_name);
        return
    end
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nValue must be type %s, but is type %s\n',class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end

