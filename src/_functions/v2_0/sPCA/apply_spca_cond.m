function [ERSP_corr,GPM_corr,output_struct] = apply_spca_cond(EEG_cond,base_mean,COEFFS,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here

WAVELET_STRUCT = struct('t',[0,1/EEG_cond.srate],...
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
addRequired(p,'EEG_cond',@isstruct);
addRequired(p,'base_mean',@isnumeric);
addRequired(p,'COEFFS',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'SPCA_PARAMS',SPCA_PARAMS,@(x) validate_struct(x,SPCA_PARAMS));
addParameter(p,'WAVELET_STRUCT',WAVELET_STRUCT,@(x) validate_struct(x,WAVELET_STRUCT));
parse(p,EEG_cond,base_mean,COEFFS,varargin{:});
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
        if isempty(EEG_cond.data)
            EEG_cond = eeg_checkset(EEG_cond,'loaddata');
        end
        data = permute(EEG_cond.data, [2,1]); % pnts x chans
        n_comps = EEG_cond.nbchan;
    case 'component'
        if isempty(EEG_cond.icaact)
            EEG_cond = eeg_checkset(EEG_cond,'loaddata');
            fprintf('%s) Recalculating ICA activations\n',EEG_cond.subject);
            EEG_cond.icaact = (EEG_cond.icaweights*EEG_cond.icasphere)*EEG_cond.data(EEG_cond.icachansind,:);
            EEG_cond.icaact = reshape( EEG_cond.icaact, size(EEG_cond.icaact,1), EEG_cond.pnts, EEG_cond.trials);
        end
        data = permute(EEG_cond.icaact, [2,1]); % pnts x chans! --> BS way?
        n_comps = size(data, 2);
    otherwise
        fprintf('Using channel data as default...\n');
        data = permute(EEG_cond.data, [2,1]); % pnts x chans
        n_comps = EEG_cond.nbchan;
end
%% ===================================================================== %%
[TF,noise_cov,params] = eeg_txf_decomp(EEG_cond,SPCA_PARAMS.analysis_type,'WAVELET_STRUCT',WAVELET_STRUCT);
%- take magnitude (not power), pnts x chans x freqs
TF = abs(TF);
%% MARTIN S RESAMPLING CODE (ACTS AS TIMEWARPING) ====================== %%
hs_min_max = SPCA_PARAMS.epoch_min_max*EEG_cond.srate; % time of next RHS in s
idx_hs = find(strcmp({EEG_cond.event.type}, SPCA_PARAMS.event_char));
gait_tf = zeros(length(idx_hs)-1,SPCA_PARAMS.n_resamples,n_comps*n_freqs); %strides/trials x pnts x chans x freqs
%## METHOD 1
%{
gait_ersp_out = zeros(output_struct.cycle_cnt,SPCA_PARAMS.n_resamples,n_comps,n_freqs); %strides/trials x pnts x chans x freqs
gait_gpm_out = zeros(output_struct.cycle_cnt,SPCA_PARAMS.n_resamples,n_comps,n_freqs); %strides/trials x pnts x chans x freqs
%}
%## METHOD 2
%- step counter, increased for each valid step
% iter = 1;
cnt = 1;
%- resample each stride to the same legth (100 pnts)
for cycle_cnt = 1:length(idx_hs)-1
    %- find first and last sample of stride
    cycle_edge = round([EEG_cond.event(idx_hs(cycle_cnt)).latency,...
        EEG_cond.event(idx_hs(cycle_cnt+1)).latency-1]); % first and last frame of gait cycle
    %- labels of all events within this cycle
    cycle_event = {EEG_cond.event([idx_hs(cycle_cnt):idx_hs(cycle_cnt+1)]).type};
    %- only keep labels of gait events to check their order:
    cycle_gaitEvent = cycle_event(contains(cycle_event,SPCA_PARAMS.timewarp_events));
    %-
    if hs_min_max(1) <= cycle_edge(2)-cycle_edge(1) &&... % check time until next HS
            cycle_edge(2)-cycle_edge(1) <= hs_min_max(2) &&...
            all(ismember(SPCA_PARAMS.timewarp_events,cycle_gaitEvent)) % oder of gait events correct
        %## METHOD 1
        %{
        tf_cycle = TF_new(cycle_edge(1):cycle_edge(2),:,:); % extract data
        tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*N_FREQS); % reshape to be able to use the resample function, skip but resample over different dimension?
        gait_tf = resample(tf_cycle,SPCA_PARAMS.n_resamples,cycle_edge(2)-cycle_edge(1)+1,0);
        gait_tf = reshape(gait_tf,SPCA_PARAMS.n_resamples,n_comps,n_freqs);
        gait_tf = 20*bsxfun(@minus,log10(gait_tf), log10(base_mean));
        [ERSP_corr, GPM_corr, PSC1,~,~] = specPCAdenoising(gait_tf,COEFFs);
        gait_ersp_out(cnt,:,:,:) = ERSP_corr;
        gait_gpm_out(cnt,:,:,:) = GPM_corr;
        %##
        fig = plot_txf(squeeze(ERSP_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_trial%i_ersp_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i},cnt)]);
        %##
        fig = plot_txf(squeeze(GPM_corr(:,CHAN_INT,:)),[],WAVELET_STRUCT.f);
        exportgraphics(fig,[fpath filesep sprintf('chan%i_cond%s_trial%i_gpm_corr.jpg',CHAN_INT,SPCA_PARAMS.condition_gait{cond_i},cnt)]);
        %}
        %## METHOD 2
        tf_cycle = TF(cycle_edge(1):cycle_edge(2),:,:); % extract data
        tf_cycle = reshape(tf_cycle,size(tf_cycle,1),n_comps*n_freqs); % reshape to be able to use the resample function, skip but resample over different dimension?
        gait_tf(cnt,:,:) = resample(tf_cycle,SPCA_PARAMS.n_resamples,cycle_edge(2)-cycle_edge(1)+1,0); % resample and store
        cnt = cnt+1;
    end
end
%## METHOD 1
%{
gait_avg = squeeze(mean(gait_ersp_out)); 
gait_gpm_avg = squeeze(mean(gait_gpm_out));
par_save(gait_avg,fpath,sprintf('cond%s_spca_ersp_corr.mat',SPCA_PARAMS.condition_gait{cond_i}));
par_save(gait_gpm_avg,fpath,sprintf('cond%s_spca_gpm_corr.mat',SPCA_PARAMS.condition_gait{cond_i}));
par_save(gait_ersp_out,fpath,sprintf('cond%s_spca_ersp_alltrials.mat',SPCA_PARAMS.condition_gait{cond_i}));
par_save(gait_gpm_out,fpath,sprintf('cond%s_spca_gpm_alltrials.mat',SPCA_PARAMS.condition_gait{cond_i}));
%}
%## METHOD 2
gait_tf = reshape(gait_tf,size(gait_tf,1),SPCA_PARAMS.n_resamples,n_comps,n_freqs); % reshape to trials x pnts x chans x freqs
gait_avg = squeeze(mean(gait_tf)); % average over trials
%-
ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_mean));
%## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
GPM = bsxfun(@minus,ERDS,mean(ERDS));
%## SPCA denoising
[ERSP_corr, GPM_corr, PSC1,~,~] = specPCAdenoising(ERDS,COEFFS);
fprintf('\n%s) Plotting validations...\n',EEG_cond.subject);
disp([num2str(round(cnt/cycle_cnt*100)) '% of the gait cycles are valid'])

% output_struct = struct('cycle_cnt',cycle_cnt,'valid_cycle_cnt',cnt,'baseline_ersp',tf_rest,'baseline_cov',noise_cov,'contwt_struct',contwt_struct);
output_struct = struct('cycle_cnt',cycle_cnt,...
    'valid_cycle_cnt',cnt,...
    'baseline_ersp',base_mean,...
    'baseline_cov',noise_cov,...
    'morlet_params',params,...
    'psc1',PSC1,...
    'erds_orig',ERDS,...
    'gpm_orig',GPM);
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

