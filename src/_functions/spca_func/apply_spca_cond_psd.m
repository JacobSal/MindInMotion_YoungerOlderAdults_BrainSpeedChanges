function [psd_out,psd_avg,output_struct] = apply_spca_cond_psd(f_mat_struct,base_f_mean,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%## Define Parser
p = inputParser;
COEFFS = [];
COND_STR = [];
%## REQUIRED
addRequired(p,'f_mat_struct',@isstruct);
addRequired(p,'base_f_mean',@isnumeric);
%## OPTIONAL
addParameter(p,'COEFFS',COEFFS,@isnumeric);
addParameter(p,'COND_STR',COND_STR,@ischar);
%## PARAMETER
parse(p,f_mat_struct,base_f_mean,varargin{:});
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
fn = fieldnames(f_mat_struct);
inds = find(contains(fieldnames(f_mat_struct),'comp'));
test = f_mat_struct.(fn{inds(1)});
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
psd_f = zeros(size(test,1),size(test,2),length(inds),'double');
for i = 1:length(inds)
    psd_f(:,:,i) = f_mat_struct.(fn{inds(i)}); % freq x epoch x chan
end

if ~isempty(COND_STR)
    inds = strcmp({f_mat_struct.trialinfo.cond},COND_STR);
    psd_f = squeeze(psd_f(:,inds,:));
end
%- average across trials
psd_f = permute(psd_f,[2,1,3]);
psd_avg = mean(psd_f); % average over trials and take abs() to ensure amplitude calc
psd_avg = permute(psd_avg,[1,3,2]);
% gait_avg = permute(gait_avg,[]);
%-
psd_baselined = bsxfun(@minus,psd_avg,base_f_mean);
% ERDS = 20*bsxfun(@minus,log10(gait_avg), log10(base_f_mean)); % decibel calculation
%## further baseline correct to dB change to mean gait cycle baseline (aka gait power modulation)
% GPM = bsxfun(@minus,psd_baselined,mean(psd_baselined));
% GPM = [];
%## SPCA denoising
if ~isempty(COEFFS)
    [psd_out, ~, psc1,~,~] = specPCAdenoising(psd_baselined,COEFFS);
    output_struct = struct('baseline_psd',base_f_mean,...
        'psc1',psc1,...
        'psd_orig_avg',psd_avg,...
        'psd_orig_baselined',psd_baselined);
else
    psd_out = psd_baselined;
    output_struct = struct('baseline_psd',base_f_mean,...
        'psc1',[],...
        'psd_orig_avg',psd_avg,...
        'psd_orig_baselined',psd_baselined);
end

end

