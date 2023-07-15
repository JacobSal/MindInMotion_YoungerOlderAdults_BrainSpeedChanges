function [STUDY,ALLEEG] = nj_create_study(ALLEEG,study_fName,study_fPath,varargin)
%MIM_CREATE_STUDY Summary of this function goes here
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
THRESH_BRAIN_SCORE = 8;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'study_fName',@ischar);
addRequired(p,'study_fPath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,ALLEEG,study_fName,study_fPath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
% DIP_NUM = 1;
% DIP_PLOT = 'off';
%## MAKE DIRS
if ~exist(study_fPath,'dir')
    mkdir(study_fPath);
end
%% ===================================================================== %%
%## DIPOLE REJECTION

tmp_rmv_subjs = zeros(1,length(ALLEEG));
%## DIPOLE REJECTION
parfor subj_i = 1:length(ALLEEG)
    %- use rejection criteria from NJ to determine bad ic's
    %- log good & bad components
    tmp_bad = find(ALLEEG(subj_i).reject.gcompreject);
    tmp_good = find(~ALLEEG(subj_i).reject.gcompreject);
    ALLEEG(subj_i).etc.urreject = [];
    ALLEEG(subj_i).etc.urreject.crit = [];
    ALLEEG(subj_i).etc.urreject.ic_keep = [];
    ALLEEG(subj_i).etc.urreject.ic_rej = [];
    ALLEEG(subj_i).etc.urreject.dipfit = [];
    if isempty(tmp_good)
        fprintf('** Subject %s has 0 brain components\n',ALLEEG(subj_i).subject);
    else
%         ALLEEG(subj_i).etc.urreject.crit = reject_struct;
        ALLEEG(subj_i).etc.urreject.ic_keep = tmp_good;
        ALLEEG(subj_i).etc.urreject.ic_rej = tmp_bad;
        ALLEEG(subj_i).etc.urreject.dipfit = ALLEEG(subj_i).dipfit;
        fprintf('** Subject %s has %i brain components\n',ALLEEG(subj_i).subject, length(tmp_good));
    end
    %- remove IC's from struct
    if length(ALLEEG(subj_i).etc.urreject.ic_keep) < 2 %|| isempty(ALLEEG(subj_i).etc.urreject)
        fprintf('** Subject %s rejected.\n',ALLEEG(subj_i).subject);
        tmp_rmv_subjs(subj_i) = 1;
    else
        % (07/06/2023) JS, make sure to checkset before editing
        % EEG.icachansind otherwise dialogue box will appear and make
        % things not work on hpg.
%         ALLEEG(subj_i) = eeg_checkset(ALLEEG(subj_i),'loaddata');
%         ALLEEG(subj_i).icachansind = ALLEEG(subj_i).etc.urreject.ic_keep;
%         if isempty(ALLEEG(subj_i).icaact)
%             fprintf('%s) Recalculating ICA activations\n',ALLEEG(subj_i).subject);
%             ALLEEG(subj_i).icaact = (ALLEEG(subj_i).icaweights*ALLEEG(subj_i).icasphere)*ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:);
% %             ALLEEG(subj_i).icaact = reshape(ALLEEG(subj_i).icaact,size(ALLEEG(subj_i).icaact,1),ALLEEG(subj_i).pnts,ALLEEG(subj_i).trials);
%         end
%         ica_weights = (ALLEEG(subj_i).icaact/ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:))/ALLEEG(subj_i).icasphere;
%         ica_winv = (ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:)/ALLEEG(subj_i).icaact);
%         tmp = sum(sqrt((ALLEEG(subj_i).icaweights-ica_weights).^2),[1,2]);
%         fprintf('%7ssum(sqrt((icaweights_new-icaweights_old).^2)) = %0.3f\n','',tmp);
%         %## Update ALLEEG & comps_out
%         ALLEEG(subj_i).icaweights = ica_weights;
%         ALLEEG(subj_i).icawinv = ica_winv;
        %- (07/06/2023) JS, code taken from EEGLAB
        components = ALLEEG(subj_i).etc.urreject.ic_rej;
        fprintf('Computing projection and removing %d components ....\n', length(components));
        component_keep = setdiff_bc(1:size(ALLEEG(subj_i).icaweights,1), components);
        compproj = ALLEEG(subj_i).icawinv(:, component_keep)*eeg_getdatact(ALLEEG(subj_i), 'component', component_keep, 'reshape', '2d');
        compproj = reshape(compproj, size(compproj,1), ALLEEG(subj_i).pnts, ALLEEG(subj_i).trials);
        ALLEEG(subj_i).data(ALLEEG(subj_i).icachansind,:,:) = compproj;
        ALLEEG(subj_i).setname = [ ALLEEG(subj_i).setname ' pruned with ICA'];
        ALLEEG(subj_i).icaact  = [];
        goodinds    = setdiff_bc(1:size(ALLEEG(subj_i).icaweights,1), components);
        ALLEEG(subj_i).icawinv     = ALLEEG(subj_i).icawinv(:,goodinds);
        ALLEEG(subj_i).icaweights  = ALLEEG(subj_i).icaweights(goodinds,:);
        ALLEEG(subj_i).specicaact  = [];
        ALLEEG(subj_i).specdata    = [];
        ALLEEG(subj_i).reject      = [];

        %- iclabel mdos
        if isfield(ALLEEG(subj_i).etc, 'ic_classification')
            if isfield(ALLEEG(subj_i).etc.ic_classification, 'ICLabel') 
                if isfield(ALLEEG(subj_i).etc.ic_classification.ICLabel, 'classifications')
                    if ~isempty(ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications)
                        ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications = ALLEEG(subj_i).etc.ic_classification.ICLabel.classifications(goodinds,:);
                    end
                end
            end
        end
        %- dipfit mods
        try
            ALLEEG(subj_i).dipfit.model = ALLEEG(subj_i).dipfit.model(goodinds);
        catch e
            fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,ALLEEG(subj_i).subject,getReport(e));
        end
    end
end
ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
%% CREATE STUDY
% initiailize study
fprintf('\n==== Making Study Modifications ====\n')
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName,...
                                'filename',study_fName,...
                                'filepath',study_fPath);
%## (NJ) DIPOLE STRUCT
% STUDY.urcluster = tmp_cluster;
% STUDY.cluster = tmp_cluster;
% STUDY.etc.rmvd_subj.inds = tmp_rmv_subjs;
%## SAVE
% parfor (subj_i = 1:length(ALLEEG),POOL_SIZE)
parfor subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).etc.full_setfile.filename = ALLEEG(subj_i).filename;
    ALLEEG(subj_i).etc.full_setfile.filepath = ALLEEG(subj_i).filepath;
    ALLEEG(subj_i).filename = sprintf('%s_%s_ICA_TMPEEG',ALLEEG(subj_i).subject,'reducedcomps');
    ALLEEG(subj_i) = pop_saveset(ALLEEG(subj_i),'filename',ALLEEG(subj_i).filename,'filepath',ALLEEG(subj_i).filepath);
end
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
end
%% SUBFUNCTIONS