function [ersps_baseline, noise_cov] = txf_baseline(EEG,analysis_type,WAVELET_STRUCT)
%## v1.2 Nadine Jacobsen, April 2022: added IC option
%% ===================================================================== %%
% (12/9/2023) JS, probably a bug with eeg_checkset with current pipeline.
% It won't load the data and deletes the icaact.
switch analysis_type
    case 'channel'
        if isempty(EEG.data)
            EEG = eeg_checkset(EEG,'loaddata');
        end
        data = permute(EEG.data, [2,1]); % pnts x chans
%         n_comps = EEG_gait.nbchan;
    case 'component'
        if isempty(EEG.icaact)
            EEG = eeg_checkset(EEG,'loaddata');
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        data = permute(EEG.icaact, [2,1]); % pnts x chans! --> BS way?
%         n_comps = size(EEG_gait.icasphere, 1);
    otherwise
        fprintf('Using channel data as default...\n');
        data = permute(EEG.data, [2,1]); % pnts x chans
%         n_comps = EEG_gait.nbchan;
end
%- CAR (common average refrence), make sure you have 'clean' data before
data = bsxfun(@minus, data, mean(data,2));
%- compute covariance matrix for inverse kernel computations (brainstorm)
noise_cov = cov(data);
%- Time-frequency analysis (function adapted by Seeber from brainstorm)
[ersps_baseline,~] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared);
end