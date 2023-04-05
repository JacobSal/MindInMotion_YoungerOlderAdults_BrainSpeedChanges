function [SUBJ,meanMat,stdvMat,medMat,rngeMat] = gen_connMatrix(SUBJ,catConn,catConnMask,varargin)
%GEN_CONNMATRIX Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {[]};
p = inputParser;
%## REQUIRED
addRequired(p,'SUBJ',@isstruct);
addRequired(p,'catConn',@isnumeric);
addRequired(p,'catConnMask',@islogical);
%## OPTIONAL
addOptional(p,'freqs',Defaults{1},@isnumeric);
%## PARAMETER
parse(p,SUBJ,catConn,catConnMask,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%- Set Defaults
freqs = p.Results.freqs;
if isempty(freqs)
    freqs = (1:size(catConn,3));
end
%% ===================================================================== %%
NUM_STATS = 4;
PAR_CPU = 20;
subjStr = SUBJ.subjStr;
connValMat = cell(size(catConn,1),size(catConn,2));
tic
ticBytes(gcp)
parfor (i = 1:size(catConn,1),PAR_CPU)
    extract_stat = zeros(NUM_STATS,size(catConn,2));
    for j = 1:size(catConn,2)
        %- connectivity extraction
        ConnMatrix  = squeeze(catConn(i,j,freqs,:)).*catConnMask(i,j,freqs,:);
        %- color limits handle
        tmp = ConnMatrix(:);
        tmp(tmp == 0) = nan();
        rnge = range(tmp);
        mn = nanmean(tmp,1);
        stdv = nanstd(tmp);
        med = nanmedian(tmp);        
        msg = [...
        sprintf('==== %s) Component %i to Component %i ====\n',subjStr,i,j),...
        sprintf('Range: %0.5f\n',rnge),...
        sprintf('Mean: %0.5f\n',mn),...
        sprintf('Standard Deviation: %0.5f\n',stdv),...
        sprintf('Median: %0.5f\n',med),...
        sprintf('Min & Max: %0.5f & %0.5f\n',min(tmp),max(tmp))]
        disp(msg);
        extract_stat(1,j) = mn;
        extract_stat(2,j) = stdv;
        extract_stat(3,j) = med;
        extract_stat(4,j) = rnge;
    end
    connValMat{i} = extract_stat;
end
toc
tocBytes(gcp)
%- extract stats
statMat = cat(3,connValMat{:});
%- seperate stats
meanMat = squeeze(statMat(1,:,:));
stdvMat = squeeze(statMat(2,:,:));
medMat = squeeze(statMat(3,:,:));
rngeMat = squeeze(statMat(4,:,:));
%- save to SUBJSTRUCT
SUBJ.cnctAnl.connMatrix.statMat = statMat;
SUBJ.cnctAnl.connMatrix.dimLabs = {{'mean','standard dev','median','range'},'componentTO','componentFROM'};
end

