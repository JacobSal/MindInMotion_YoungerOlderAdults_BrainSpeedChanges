function [statMat,extract_sig] = gen_connMatrix(EEG,conn_name,varargin)
%GEN_CONNMATRIX Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
freqs = EEG.CAT.Conn.freqs; %(1:size(EEG.CAT.Conn.(conn_name),3));
stat_alpha = 0.05;
%- 
NUM_STATS = 4;
DO_MASKING = false;
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'conn_name',@ischar);
%## OPTIONAL
addOptional(p,'freqs',freqs,@isnumeric);
addOptional(p,'stat_alpha',stat_alpha,@isnumeric);
%## PARAMETER
addParameter(p,'DO_MASKING',DO_MASKING,@islogical);
addParameter(p,'NUM_STATS',NUM_STATS,@islogical);
parse(p,EEG,conn_name,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
freqs = p.Results.freqs;
stat_alpha = p.Results.stat_alpha;
%- PARAMETER
DO_MASKING = p.Results.DO_MASKING;
NUM_STATS = p.Results.NUM_STATS;

%% ===================================================================== %%
subjStr = EEG.subject;
connValMat = cell(size(EEG.CAT.Conn.(conn_name),1),size(EEG.CAT.Conn.(conn_name),2));
extract_sig = zeros(size(EEG.CAT.Conn.(conn_name),1),size(EEG.CAT.Conn.(conn_name),2),size(EEG.CAT.Conn.(conn_name),4));
catConnMask = logical(ones(size(EEG.CAT.Conn.(conn_name))));
%- 
if isfield(EEG.CAT,'Stats') && DO_MASKING
    catConnMask = EEG.CAT.Stats.(conn_name).pval < stat_alpha;
end
%- 
tic
% parfor (i = 1:size(catConn,1),PAR_CPU)
for i = 1:size(EEG.CAT.Conn.(conn_name),1)
    extract_stat = zeros(NUM_STATS,size(EEG.CAT.Conn.(conn_name),2));
    for j = 1:size(EEG.CAT.Conn.(conn_name),2)
        %- connectivity extraction
        ConnMatrix  = squeeze(EEG.CAT.Conn.(conn_name)(i,j,freqs,:)).*squeeze(catConnMask(i,j,freqs,:));
        %- color limits handle
        tmp = sum(ConnMatrix,1);
        tmp(tmp == 0) = nan();
        
        rnge = range(tmp);
        mn = nanmean(tmp);
        stdv = nanstd(tmp);
        med = nanmedian(tmp);        
        msg = [...
        sprintf('==== %s) Component %i to Component %i ====\n',subjStr,i,j),...
        sprintf('Range: %0.5f\n',rnge),...
        sprintf('Mean: %0.5f\n',mn),...
        sprintf('Standard Deviation: %0.5f\n',stdv),...
        sprintf('Median: %0.5f\n',med),...
        sprintf('Min & Max: %0.5f & %0.5f\n',min(tmp),max(tmp))];
        disp(msg);
        extract_sig(i,j,:) = tmp;
        extract_stat(1,j) = mn;
        extract_stat(2,j) = stdv;
        extract_stat(3,j) = med;
        extract_stat(4,j) = rnge;
    end
    connValMat{i} = extract_stat;
end
%- extract stats
statMat = cat(3,connValMat{:});
%- seperate stats
% meanMat = squeeze(statMat(1,:,:));
% stdvMat = squeeze(statMat(2,:,:));
% medMat = squeeze(statMat(3,:,:));
% rngeMat = squeeze(statMat(4,:,:));
%## TIME
toc
end

