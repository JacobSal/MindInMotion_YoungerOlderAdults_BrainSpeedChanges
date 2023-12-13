function [gait_avg,ERDS,GPM,TF,output_struct] = spca_time_freq_decomp(EEG_gait,EEG_baseline,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here

WAVELET_STRUCT = struct('t',[0,1/EEG_gait.srate],...
    'f',(4:100),...
    'fc',1,...
    'FWHM_tc',3,...
    'squared','n');
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
addRequired(p,'EEG_baseline',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'SPCA_PARAMS',SPCA_PARAMS,@(x) validate_struct(x,SPCA_PARAMS));
addParameter(p,'WAVELET_STRUCT',WAVELET_STRUCT,@(x) validate_struct(x,WAVELET_STRUCT));
parse(p,EEG_gait,EEG_baseline,varargin{:});
%## SET DEFAULTS
SPCA_PARAMS = p.Results.SPCA_PARAMS;
WAVELET_STRUCT = p.Results.WAVELET_STRUCT;
%- ASSIGNED VALUES
n_freqs = length(WAVELET_STRUCT.f);
%% ===================================================================== %%
% (12/9/2023) JS, probably a bug with eeg_checkset with current pipeline.
% It won't load the data and deletes the icaact.
switch SPCA_PARAMS.analysis_type
    case 'channel'
        if isempty(EEG_gait.data)
            EEG_gait = eeg_checkset(EEG_gait,'loaddata');
        end
        data = permute(EEG_gait.data, [2,1]); % pnts x chans
        n_comps = EEG_gait.nbchan;
    case 'component'
        if isempty(EEG_gait.icaact)
            EEG_gait = eeg_checkset(EEG_gait,'loaddata');
            fprintf('%s) Recalculating ICA activations\n',EEG_gait.subject);
            EEG_gait.icaact = (EEG_gait.icaweights*EEG_gait.icasphere)*EEG_gait.data(EEG_gait.icachansind,:);
            EEG_gait.icaact = reshape( EEG_gait.icaact, size(EEG_gait.icaact,1), EEG_gait.pnts, EEG_gait.trials);
        end
        data = permute(EEG_gait.icaact, [2,1]); % pnts x chans! --> BS way?
        n_comps = size(data, 2);
    otherwise
        fprintf('Using channel data as default...\n');
        data = permute(EEG_gait.data, [2,1]); % pnts x chans
        n_comps = EEG_gait.nbchan;
end
%-
hs_min_max = SPCA_PARAMS.epoch_min_max*EEG_gait.srate; % time of next RHS in s
%- CAR (common average refrence), make sure you have 'clean' data before
data = bsxfun(@minus, data, mean(data,2));
%- time frequency transform
fprintf('\nRunning Time-Frequency Decomposition Using Morlet Wavelets...\n');
tt = tic;
[TF,morlet_params] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared);
fprintf('Time: %0.2f\n',toc(tt));
%## (ALTERNATIVE) MORLET WAVELET PACKAGE
%{
TF = zeros(size(data,1),size(data,2),WAVELET_STRUCT.f(end)-WAVELET_STRUCT.f(1)+1);
for chan_i = 1:size(data,2)
    [tmp,period,scale,coi,dj,param,k] = contwt(data(:,chan_i),WAVELET_STRUCT.t(2),1,1,WAVELET_STRUCT.f(1),WAVELET_STRUCT.f(end)-WAVELET_STRUCT.f(1),'Morlet',WAVELET_STRUCT.FWHM_tc);
    TF(:,chan_i,:) = tmp';
end
contwt_struct = [];
contwt_struct.period = period;
contwt_struct.scale = scale;
contwt_struct.coi = coi;
contwt_struct.dj = dj;
contwt_struct.param = param;
contwt_struct.k = k;
%}
%- take magnitude (not power), pnts x chans x freqs
TF = abs(TF);
%% MARTIN S RESAMPLING CODE (ACTS AS TIMEWARPING) ====================== %%
idx_hs = find(strcmp({EEG_gait.event.type}, SPCA_PARAMS.event_char));
gait_tf = zeros(length(idx_hs)-1,SPCA_PARAMS.n_resamples,n_comps*n_freqs); %strides/trials x pnts x chans x freqs
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
%## baseline correct to dB power change to standing baseline (also called ERSP)
[tf_rest,noise_cov] = txf_baseline(EEG_baseline,SPCA_PARAMS.analysis_type,WAVELET_STRUCT);
%- average over time, keep magnitude (not power -> would amplify outliers)
base_mean(1,:,:) = squeeze(mean(abs(tf_rest(10+1:end-10,:,:)),1));
%-
ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_mean));
%## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
GPM = bsxfun(@minus,ERDS,mean(ERDS));
% output_struct = struct('cycle_cnt',cycle_cnt,'valid_cycle_cnt',cnt,'baseline_ersp',tf_rest,'baseline_cov',noise_cov,'contwt_struct',contwt_struct);
output_struct = struct('cycle_cnt',cycle_cnt,'valid_cycle_cnt',cnt,'baseline_ersp',base_mean,'baseline_cov',noise_cov,'morlet_params',morlet_params);
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

