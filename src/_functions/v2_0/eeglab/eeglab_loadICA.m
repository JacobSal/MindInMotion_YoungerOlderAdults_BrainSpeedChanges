function [ALLEEG,EEG,CURRENTSET] = eeglab_loadICA(fName,fPath,varargin)
%EEGLAB_LOADICA Loads the ICA & EEG struct for the subject designated by
%subjStr & subDirNum.

%   IN:
%   OUT:
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220112.0
%   Previous Version: n/a
%   Summary:  
%
%## TIME
tic
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'fName',@ischar);
addRequired(p,'fPath',@ischar);
%## PARAMETER
parse(p, fName, fPath, varargin{:});
%## SET DEFAULTS
%-
OVER_WRITE_CURRENT_SET = 1;
ALLEEG = [];
CURRENTSET = 1;
%% ===================================================================== %%
fprintf(1,'==== Loading ICA ====\n');
%## load EEG file contianing AMICA results
EEG = pop_loadset('filename',fName,'filepath',fPath);
if isempty(EEG.icaweights)
    %## CONCATENATE DATA | OVERWRITE DATA
    if OVER_WRITE_CURRENT_SET    
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
    else
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    end

    %## REMOVE ICA RESULTS FROM PREVIOUS ANALYSIS
    EEG = rmfield(EEG,'icaweights');
    EEG.icaweights = [];
    EEG = rmfield(EEG,'icawinv');
    EEG.icawinv = [];
    EEG = rmfield(EEG,'icaact');
    EEG.icaact = [];
    EEG = rmfield(EEG,'icasphere');
    EEG.icasphere = [];
    EEG = rmfield(EEG,'icachansind');
    EEG.icachansind = []; %also remove icachanind?

    %## LOAD CURRENT AMICA RESULTS FOR subjStr & subDirNum
    EEG = pop_loadmodout(EEG,fPath);
    %- sometimes need to reset variables when using pop_subcomp();
%     if size(EEG.icaact,1) ~= EEG.nbchan
%         EEG = eeg_checkset(EEG, 'ica'); %Note: I have no idea why it says "LLt not set" or what LLt is
%         EEG.icaweights = EEG.etc.amica.W(:,:,1);
%         EEG.icasphere  = EEG.etc.amica.S(1:size(EEG.etc.amica.W,1),:);
%         EEG.icawinv    = EEG.etc.amica.A(:,:,1);
%         EEG.icaact     = [];
%         EEG            = eeg_checkset(EEG);
%         EEG.icaact     = eeg_getica(EEG);
%         EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications_default;
%     end
    % make sure everything is ok and calculate the activation matrix
%     pop_editoptions('option_computeica',1); %randomly needed for some people to run even if you manually selected this checkbox with the GUI
else
    fprintf('ica already available.\n');
end
%% store updated EEG structure
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );

%## TIME
toc

end
%% Version History
%{
v1.0.20210112.0, JS: initial configuration. Script adapted from
addAMICAresults

%}