function [EEG] = eeglab_pop_subcomp(EEG,comps,rmv_bool)
%EEGLAB_POP_SUBCOMP Summary of this function goes here
%   eeglabs pop_subcomp, but it works.
if logical(rmv_bool)
    comps_remove = comps;
    comps_keep = setdiff_bc(1:size(EEG.icaweights,1), comps);
else
    comps_keep = comps;
    comps_remove = setdiff_bc(1:size(EEG.icaweights,1), comps);
end
fprintf('Computing projection and removing %d components ....\n', length(comps_remove));
compproj = EEG.icawinv(:, comps_keep)*eeg_getdatact(EEG, 'component', comps_keep, 'reshape', '2d');
compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
EEG.data(EEG.icachansind,:,:) = compproj;
EEG.setname = [ EEG.setname ' pruned with ICA'];
EEG.icaact  = [];
goodinds    = setdiff_bc(1:size(EEG.icaweights,1), comps_remove);
EEG.icawinv     = EEG.icawinv(:,goodinds);
EEG.icaweights  = EEG.icaweights(goodinds,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];
% EEG.icachansind = [];
%- iclabel mods
fprintf('making iclabel modifications\n');
if isfield(EEG.etc, 'ic_classification')
    if isfield(EEG.etc.ic_classification, 'ICLabel') 
        if isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            if ~isempty(EEG.etc.ic_classification.ICLabel.classifications)
                EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodinds,:);
            end
        end
    end
end
try
    EEG.dipfit.model = EEG.dipfit.model(goodinds);
catch e
    fprintf(['error. identifier: %s\n',...
         'error. %s\n',...
         'error. on subject %s\n',...
         'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
end
end

