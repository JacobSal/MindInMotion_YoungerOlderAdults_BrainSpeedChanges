function [EEG] = rmv_ics(EEG)
%RMV_ICS Summary of this function goes here
%   Detailed explanation goes here
%%
% (07/06/2023) JS, make sure to checkset before editing
% EEG.icachansind otherwise dialogue box will appear and make
% things not work on hpg.
%- (07/06/2023) JS, code taken from EEGLAB
components = EEG.etc.urreject.ic_rej;
fprintf('Computing projection and removing %d components ....\n', length(components));
component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
EEG.data(EEG.icachansind,:,:) = compproj;
EEG.setname = [ EEG.setname ' pruned with ICA'];
EEG.icaact  = [];
goodinds    = setdiff_bc(1:size(EEG.icaweights,1), components);
EEG.icawinv     = EEG.icawinv(:,goodinds);
EEG.icaweights  = EEG.icaweights(goodinds,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];
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


