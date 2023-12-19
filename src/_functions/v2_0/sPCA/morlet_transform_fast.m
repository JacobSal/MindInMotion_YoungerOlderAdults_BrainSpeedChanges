function [P,param_struct] = morlet_transform_fast(x,t,f,fc,FWHM_tc,squared,data_type)
% MORLET_TRANSFORM: 
%     Applies complex Morlet wavelet transform to the timeseries stored in the 
%     matrix x with size (ntimes x ntimeseries). It returns a wavelet coefficient map 
%     (by default squared)
%
% INPUTS:
%    - x       : (ntimes x ntimeseries) a vector of the timeseries
%    - t       : (ntimes x 1) a vector of times (in secs)
%    - f       : (nfreqs x 1) a vector of frequencies in which to estimate the
%                wavelet transform (in Hz). Default is 1:60Hz
%    - fc      : (default is 1) central frequency of complex Morlet wavelet in Hz
%    - FWHM_tc : (default is 3) FWHM of complex Morlet wavelet in
%                time. Also see morlet_design.m. 
%    - squared : 'y' (default) or 'n'. Flag that decided whether the function returns the
%                squared coefficients (y) or not (n). Squaring represents neural power in
%                the corresponding frequency.
%
% OUTPUT:
%    - P: (x ntimes x ntimeseries x nfreqs), or
%    - P: (nTimes x nFreqs), if ntimeseries = 1
%         A matrix of wavelet coefficients (by default squared)
%
% EXAMPLE:
%    t = [0:.01:1].';
%    f = [5:20].';
%    P = morlet_transform(sin(2*pi*10*t),t,f);
%    Coefs = morlet_transform(sin(2*pi*10*t),t,f,[],[],'n');

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Dimitrios Pantazis, 2010
% adapted by Martin Seeber, 2015/17

FWHM_tc(FWHM_tc < 1) = 1;

%- signal parameters
Ts = t(2)-t(1); %sampling period of signal
Fs = 1/Ts; %sampling frequency of signal

%- complex morlet wavelet parameters
scales = f ./ fc; %scales for wavelet
nscales = length(scales);

sigma_tc(1:nscales) = FWHM_tc / sqrt(8*log(2));

%- compute wavelet kernels for each scale
W = cell(nscales,1);
precision = 3;

for s = 1:nscales
    time = scales(s)* [-precision*sigma_tc(s)/scales(s): 1/Fs : precision*sigma_tc(s)/scales(s)].';
    %-
    W{s} = morlet_wavelet(time,fc,sigma_tc(s));
    %## VALIDATION PLOT
%     W{s} = sqrt(scales(s)) * morlet_wavelet(time,fc,sigma_tc(s));
%     if ~mod(log2(s),1)
%         figure;plot(time,real(W{s}));
%     end
end
%compute wavelet coefficients
nx = size(x,2); %number of timeseries
ntimes = size(x,1); %number of timepoints

if nx > 1
    P = zeros(ntimes,nx,nscales,data_type);
    for s = 1:nscales
        P(:,:,s) = conv2(x,W{s},'same') * 2/(sum(abs(W{s})));
    end

elseif nx ==1
    P = zeros(ntimes,nscales,data_type);
    for s = 1:nscales
        P(:,s) = conv2(x,W{s},'same') * 2/(sum(abs(W{s}))); 
    end
end

%-if return squared coefficients
if strcmp(squared,'y')
    P = abs(P).^2; %return neural power
end
param_struct = struct('W',W,'scales',scales,'sigma_tc',sigma_tc,'precision',precision);
end

% % Convert dimensions in the format we want:
% % [nTimeseries x nFreqs x nTimes] => [nTimeseries x nTimes x nFreqs]
% P = permute(P, [1,3,2]);



    