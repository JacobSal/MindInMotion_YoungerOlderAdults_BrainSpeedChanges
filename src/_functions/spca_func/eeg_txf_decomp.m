function [TF,noise_cov,params] = eeg_txf_decomp(EEG,analysis_type,varargin)
%SPCA_TIME_FREQ_DECOMP Summary of this function goes here
%   Detailed explanation goes here

WAVELET_STRUCT = struct('t',[0,1/EEG.srate],...
    'f',(4:100),...
    'fc',1,...
    'FWHM_tc',3,...
    'squared','n',...
    'data_type','double');
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct);
addRequired(p,'analysis_type',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'WAVELET_STRUCT',WAVELET_STRUCT,@(x) validate_struct(x,WAVELET_STRUCT));
parse(p,EEG,analysis_type,varargin{:});
%## SET DEFAULTS
WAVELET_STRUCT = p.Results.WAVELET_STRUCT;
%% ===================================================================== %%
% (12/9/2023) JS, probably a bug with eeg_checkset with current pipeline.
% It won't load the data and deletes the icaact.
switch analysis_type
    case 'channel'
        if isempty(EEG.data)
            EEG = eeg_checkset(EEG,'loaddata');
        end
        data = permute(EEG.data, [2,1]); % pnts x chans
    case 'component'
        if isempty(EEG.icaact)
            EEG = eeg_checkset(EEG,'loaddata');
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end
        data = permute(EEG.icaact, [2,1]); % pnts x chans! --> BS way?
    otherwise
        fprintf('Using channel data as default...\n');
        data = permute(EEG.data, [2,1]); % pnts x chans
end
%- CAR (common average refrence), make sure you have 'clean' data before
data = bsxfun(@minus, data, mean(data,2));
%- compute covariance matrix for inverse kernel computations (brainstorm)
noise_cov = cov(data);
%- time frequency transform
fprintf('\nRunning Time-Frequency Decomposition Using Morlet Wavelets...\n');
tt = tic;
try
    if ispc
        user = memory;
        fprintf('\nAvailable memory: %0.2g GB\n',user.MemAvailableAllArrays/(10e9));
    end
    [TF,params] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared,WAVELET_STRUCT.data_type);
catch MExc
    fprintf('\n%s\n',getReport(MExc));
    fprintf('Trying segmentating data...\n');
    %## VARIABLE PROCEDURE (lower precision)
    fprintf('Using a lower precision array...\n');
    WAVELET_STRUCT.data_type = 'single';
%     tf_out = zeros(dim1,dim2,length(WAVELET_STRUCT.f),'uint16');
%     tf_out = zeros(dim1,dim2,length(WAVELET_STRUCT.f),'single');
    [TF,params] = morlet_transform_fast(data,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared,WAVELET_STRUCT.data_type);
    TF = abs(TF);
    %## FILE WRITE PROCEDURE
    %{
    %-
    dim1 = size(data,1);
    dim2 = size(data,2);
    user = memory;
    fprintf('\nAvailable memory: %0.2g GB\n',user.MemAvailableAllArrays/(10e9));
    %- subtract off 2GB for system
    mem_alloc = user.MemAvailableAllArrays - 2*10e9;
    %- estimate the increments needed to perform calculation
    spacings = ceil((length(WAVELET_STRUCT.f)*dim2*dim1*8^3)/(mem_alloc));
    %-
    params = [];
    fcnt = 1;
    fname = sprintf('txf_data_%i.bin',fcnt);
    while exist(fname,'file')
        fcnt = fcnt + 1;
        fname = sprintf('txf_data_%i.bin',fcnt);
    end
    fileID = fopen(fname,'w');
    inds = (1:floor(dim1/spacings));
    padn = floor(EEG.srate/pi);
    for sp = 1:spacings
        dat_in = data(inds,:);
        dat_in = padarray(dat_in,[padn,0],'both');
        fprintf('Percent Remaining: %0.1f ...\n',(spacings-sp)/spacings*100);
        [tmp,params] = morlet_transform_fast(dat_in,WAVELET_STRUCT.t,WAVELET_STRUCT.f,WAVELET_STRUCT.fc,WAVELET_STRUCT.FWHM_tc,WAVELET_STRUCT.squared);
        tmp(1:padn,:,:) = [];
        tmp(end-padn+1:end,:,:) = [];
        tmp = abs(tmp);
        tmp = floor(tmp*1000);
        if sp == 1
            fwrite(fileID,tmp,'uint16');
            fclose(fileID);
            fileID = fopen(fname,'a');
        else
            fwrite(fileID,tmp,'uint16');
        end
        inds = inds + sp*floor(dim1/spacings);
        fprintf('done\n');
    end
    fclose(fileID);
    %}
end
fprintf('Time: %0.2f\n',toc(tt));
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

