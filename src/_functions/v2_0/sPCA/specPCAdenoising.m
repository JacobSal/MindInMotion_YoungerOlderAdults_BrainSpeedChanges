function [ERSP_corr, GPM_corr, ERSP_psc1, PSC, V] = specPCAdenoising(ERSP, varargin)
% spectral PCA denoising. PCA of gait ERSP and removal of first first principal 
% component as describend in:
% Seeber M, Scherer R, Wagner J, Solis-Escalante T, Müller-Putz GR. 
% High and low gamma EEG oscillations in central sensorimotor areas are 
% conversely modulated during the human gait cycle. 
% Neuroimage [Internet]. 2015 May 15 [cited 2018 Mar 8];112:318–26. 
% Available from: https://linkinghub.elsevier.com/retrieve/pii/S105381191500227X
%
% INPUT
% - ERSP        [matrix: times x chans x freqs] mean gait ERSP, power as dB change to standing BL
% optional:
% - V           [matrix:  freqs x freqs] unmixing matrix = sorted spectral PCA eigenvectors obtained
%               from previous call to this function. E.g. to correct condition-specific
%               data with  eigenvectors from the grand average
%
% OUTPUT
% - ERSP_corr   [matrix: times x chans x freqs] corrected (first principal 
%               component removed) mean gait ERSP, power as dB change to standing BL
% - GPM_corr    [matrix: times x chans x freqs] corrected (first principal 
%               component removed) mean gait phase related power modulation, 
%               power as change to mean gait cycle BL
% - ERSP_psc1   [matrix: times x chans x freqs] only first principal 
%               component of the mean gait ERSP, power as dB change to
%               standing BL. I.e. data that was removed
% - PSC         [cell: freqs; each cell --> matrix: times x chans x freqs] all kept
%               principal components of the mean gait ERSP, power as dB change to
%               standing BL. Could be used to identify whether further PC
%               should be removed.
% - V           [matrix:  freqs x freqs] unmixing matrix = sorted spectral PCA eigenvectors.
%               Can be used to call to this function, e.g. to correct 
%               condition-specific data with eigenvectors from grand average
%
% written by Martin Seeber,FBM Lab, UniGe, 2019
% adapted by Nadine Jacobsen, University of Oldenburg, May 2022
% v1.0 last changed May-30-2022

if isempty(varargin) % check whether eigenvalues are provided, if not calculate from data
    % PCA of ERSP --> dim = num freqs
    if size(ERSP,1) == 1
        ERSP_cc = squeeze(ERSP);
        % ERSP_cc = permute(ERSP_cc,[2,1]);
        CC = cov(squeeze(ERSP_cc)); %covariance matrix: freq x freq
    else
        CC = cov(squeeze(mean(ERSP))); %covariance matrix: freq x freq
    end
    [v, d] = eig(CC); % obtain eigenvalues
    
    % find component with greatest eigenvalue
    [D, sort_ix] = sort(diag(d),'descend');
    V = v(:,sort_ix);
else
    V = varargin{1};
end

Ai = V.'; % take inverse of unmixing matrix V to get mixing matrix Ai
% transpose is quicker than inv() and sufficient since V vectors are orthogonal

% set first component to 0
ERSP_psc1 = zeros(size(ERSP));
ERSP_corr = zeros(size(ERSP));
PSC = cell(1,length(V));

% loop through all samples and correct them
for t_cnt = 1:size(ERSP,1)
    F_psc = squeeze(ERSP(t_cnt,:,:)) * V;
    ERSP_psc1(t_cnt,:,:) = F_psc(:,1) * Ai(1,:); % times x chans x freqs
    ERSP_corr(t_cnt,:,:) = F_psc(:,2:end) * Ai(2:end,:);
    for p_cnt = 1:size(F_psc,2)
        if t_cnt==1
          PSC{p_cnt} = zeros(size(ERSP));
        end
        PSC{p_cnt}(t_cnt,:,:)  = F_psc(:,p_cnt) * Ai(p_cnt,:);
    end 
end

GPM_corr = bsxfun(@minus,ERSP_corr,mean(ERSP_corr));

end
