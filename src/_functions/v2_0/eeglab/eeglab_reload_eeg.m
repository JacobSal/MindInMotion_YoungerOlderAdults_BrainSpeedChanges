function [EEG] =  eeglab_reload_eeg(fpath,fname)
%LOAD_SUBJECT Summary of this function goes here
%
%   IN:
%       fpath,
%       fname,
%   OUT:
%       struct_out, STRUCT ARRAY
%
%   Version History --> See details at the end of the script.
%   Previous Version: n/a
%   Summary:  
%
cat_logo();
%## TIME
t = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'fpath',@ischar);
addRequired(p,'fname',@ischar);
%## PARSE
parse(p, fpath, fname);
%## SET DEFAULTS
%% ===================================================================== %%
EEG = pop_loadset('filepath',fpath,...
    'filename',fname);
%-
EEG = eeg_checkset(EEG,'loaddata');
if isempty(EEG.icaact)
    fprintf('%s) Recalculating ICA activations\n',EEG.subject);
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end
fprintf('eeglab_load_subject completed in %0.2fs',toc(t));
end

