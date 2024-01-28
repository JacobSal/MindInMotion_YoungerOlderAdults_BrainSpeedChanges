function [ersp_out,gpm_out,gait_avg,output_struct] = apply_spca_cond_timewarp(txf_mat_struct,base_txf_mean,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## Define Parser
p = inputParser;
COEFFS = [];
COND_STR = [];
%## REQUIRED
addRequired(p,'txf_mat_struct',@isstruct);
addRequired(p,'base_txf_mean',@isnumeric);
%## OPTIONAL
addParameter(p,'COEFFS',COEFFS,@isnumeric);
addParameter(p,'COND_STR',COND_STR,@ischar);
%## PARAMETER
parse(p,txf_mat_struct,base_txf_mean,varargin{:});
%## SET DEFAULTS
COEFFS = p.Results.COEFFS;
COND_STR = p.Results.COND_STR;
%% (ALTERNATIVE) USE TIMEWARPPING EEGLAB =================================== %%
% functions: newtimef.m, timefreq.m, timewarp.m
% (12/9/2023) JS, the timewarp function is the real meat that is needed to
% warp gait events to a grand average. Could use the morlet_transform_fast
% after that.
%##
%- reshape data [trials x pnts x chans x freq]
fn = fieldnames(txf_mat_struct);
inds = find(contains(fieldnames(txf_mat_struct),'comp'));
test = txf_mat_struct.(fn{inds(1)});
%##
%{
% gait_txf = zeros(size(test,3),size(test,2),length(inds),size(test,1),'single');
gait_txf = zeros(size(test,3),size(test,2),length(inds),size(test,1),'double');

for i = 1:length(inds)
    gait_txf(:,:,i,:) = reshape(txf_mat_struct.(fn{inds(i)}),size(test,3),size(test,2),1,size(test,1));
end
%}
%##
% gait_txf = zeros(size(test,3),size(test,2),length(inds),size(test,1),'single');
gait_txf = zeros(size(test,1),size(test,2),size(test,3),length(inds),'double');

for i = 1:length(inds)
    gait_txf(:,:,:,i) = txf_mat_struct.(fn{inds(i)}); % freq x time x epoch x chan
end

if ~isempty(COND_STR)
    inds = strcmp({txf_mat_struct.trialinfo.cond},COND_STR);
    gait_txf = squeeze(gait_txf(:,:,inds,:));
end
%- average across trials
gait_avg = squeeze(mean(abs(gait_txf),3)); % average over trials and take abs() to ensure amplitude calc
gait_avg = permute(gait_avg,[2,3,1]);
%-
ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_txf_mean)); % decibel calculation
%## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
GPM = bsxfun(@minus,ERDS,mean(ERDS));
%## SPCA denoising
if ~isempty(COEFFS)
    [ersp_out, gpm_out, PSC1,~,~] = specPCAdenoising(ERDS,COEFFS);
    output_struct = struct('baseline_ersp',base_txf_mean,...
        'psc1',PSC1,...
        'erds_orig',ERDS,...
        'gpm_orig',GPM);
else
    ersp_out = ERDS;
    gpm_out = GPM;
    output_struct = struct('baseline_ersp',base_txf_mean,...
        'psc1',[],...
        'erds_orig',[],...
        'gpm_orig',[]);
end

end

