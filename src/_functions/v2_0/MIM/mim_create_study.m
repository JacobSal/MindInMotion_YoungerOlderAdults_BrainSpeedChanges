function [STUDY,ALLEEG] = mim_create_study(ALLEEG,study_fName,study_fPath,varargin)
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
STUDY = [];
%% ===================================================================== %%
tmp_rmv_subjs = zeros(1,length(ALLEEG));
%## DIPOLE REJECTION
brain_ics_count = zeros(length(ALLEEG),1);
brain_ics_scores = cell(length(ALLEEG),1);
powpow_rej_count = zeros(length(ALLEEG),1);
all_rej_count = zeros(length(ALLEEG),1);
powpow_rej_nums = cell(length(ALLEEG),1);
brain_ic_nums = cell(length(ALLEEG),1);
all_rej_nums = cell(length(ALLEEG),1);
subj_str = cell(length(ALLEEG),1);
good_diplocs_orig = cell(length(ALLEEG),1);
bad_diplocs_orig = cell(length(ALLEEG),1);
bad_diplocs_mni = cell(length(ALLEEG),1);
good_diplocs_mni = cell(length(ALLEEG),1);
parfor subj_i = 1:length(ALLEEG)
% for subj_i = 1:length(ALLEEG)
    fprintf('Rejecting IC''s for subject %s...\n',ALLEEG(subj_i).subject);
    %- use rejection criteria to determine bad ic's
    reject_struct = mim_reject_ics(ALLEEG(subj_i),ALLEEG(subj_i).filepath);
    %- log good & bad components
    chk = (reject_struct.IC_all_brain >= THRESH_BRAIN_SCORE & reject_struct.IC_all_brain ~= 9);
    chk_w_powpow = unique(cat(1,find(chk),reject_struct.IC_powpow_rej));
    tmp_bad = setdiff(find((1:size(ALLEEG(subj_i).icaweights,1))),chk_w_powpow);
    tmp_good = chk_w_powpow';
    ALLEEG(subj_i).etc.urreject = [];
    ALLEEG(subj_i).etc.urreject.crit = [];
    ALLEEG(subj_i).etc.urreject.ic_keep = [];
    ALLEEG(subj_i).etc.urreject.ic_rej = [];
    ALLEEG(subj_i).etc.urreject.dipfit = [];
    if isempty(tmp_good)
        tmp_good = 0;
        fprintf('** Subject %s has 0 brain components\n',ALLEEG(subj_i).subject);
    else
        ALLEEG(subj_i).etc.urreject.crit = reject_struct;
        ALLEEG(subj_i).etc.urreject.ic_keep = tmp_good;
        ALLEEG(subj_i).etc.urreject.ic_rej = tmp_bad;
        ALLEEG(subj_i).etc.urreject.dipfit = ALLEEG(subj_i).dipfit;
        fprintf('** Subject %s has %i brain components\n',ALLEEG(subj_i).subject, length(tmp_good));
    end
    %- table printouts
    brain_ics_count(subj_i) = length(tmp_good);
    powpow_rej_count(subj_i) = length(reject_struct.IC_powpow_rej);
    all_rej_count(subj_i) = length(tmp_bad);
    tmp = reject_struct.IC_powpow_rej;
    if ~isempty(tmp)
        powpow_rej_nums{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %reject_struct.IC_powpow_rej;
    else
        powpow_rej_nums{subj_i} = '';
    end
    if ~isempty(tmp_good)
        brain_ic_nums{subj_i} = ['[' sprintf('%i,',tmp_good(1:end-1)) sprintf('%i',tmp_good(end)) ']']; %sprintf('%i,',tmp_good); %tmp_good;
    else
        brain_ic_nums{subj_i} = '';
    end
    if ~isempty(tmp_bad)
        all_rej_nums{subj_i} = ['[' sprintf('%i,',tmp_bad(1:end-1)) sprintf('%i',tmp_bad(end)) ']']; %sprintf('%i,',tmp_bad); %tmp_bad
    else
        all_rej_nums{subj_i} = '';
    end
    tmp = reject_struct.IC_all_brain;
    tmp = tmp(chk);
    if ~isempty(tmp)
        brain_ics_scores{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %sprintf('%i,',tmp(chk)); %tmp(chk);
    else
        brain_ics_scores{subj_i} = '';
    end
    subj_str{subj_i} = ALLEEG(subj_i).subject;
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
        %- iclabel mods
        fprintf('making iclabel modifications\n');
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
        fprintf('making dipfit modifications\n');
        tmp = {ALLEEG(subj_i).dipfit.model(tmp_good).pos_old};
%             tmp = tmp(chk);
        if ~isempty(tmp)
            print_tmp = [];
            for i = 1:length(tmp)
                print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
            end
            good_diplocs_orig{subj_i} = print_tmp; %sprintf('%i,',tmp(chk)); %tmp(chk);
        else
            good_diplocs_orig{subj_i} = '';
        end
        %-
        tmp = {ALLEEG(subj_i).dipfit.model(tmp_bad).pos_old};
%             tmp = tmp(chk);
        if ~isempty(tmp)
            print_tmp = [];
            for i = 1:length(tmp)
                print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
            end
            bad_diplocs_orig{subj_i} = print_tmp; %sprintf('%i,',tmp(chk)); %tmp(chk);
%                 bad_diplocs{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %sprintf('%i,',tmp(chk)); %tmp(chk);
        else
            bad_diplocs_orig{subj_i} = '';
        end
        %-
        tmp = {ALLEEG(subj_i).dipfit.model(tmp_bad).mnipos};
        if ~isempty(tmp)
            print_tmp = [];
            for i = 1:length(tmp)
                print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
            end
            bad_diplocs_mni{subj_i} = print_tmp; %sprintf('%i,',tmp(chk)); %tmp(chk);
%                 bad_diplocs{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %sprintf('%i,',tmp(chk)); %tmp(chk);
        else
            bad_diplocs_mni{subj_i} = '';
        end
        tmp = {ALLEEG(subj_i).dipfit.model(tmp_good).mnipos};
        if ~isempty(tmp)
            print_tmp = [];
            for i = 1:length(tmp)
                print_tmp = [print_tmp,sprintf('[%0.2f,%0.2f,%0.2f],',tmp{i})];
            end
            good_diplocs_mni{subj_i} = print_tmp; %sprintf('%i,',tmp(chk)); %tmp(chk);
%                 bad_diplocs{subj_i} = ['[' sprintf('%i,',tmp(1:end-1)) sprintf('%i',tmp(end)) ']']; %sprintf('%i,',tmp(chk)); %tmp(chk);
        else
            good_diplocs_mni{subj_i} = '';
        end
        try
            ALLEEG(subj_i).dipfit.model = ALLEEG(subj_i).dipfit.model(goodinds);
        catch e
            fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,ALLEEG(subj_i).subject,getReport(e));
        end
    end
    fprintf('DONE. rejecting IC''s for subject %s.\n',ALLEEG(subj_i).subject);
end
ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
%% WRITE TO XLSX
fprintf('Writing subject rejection criteria table...\n');
tmp_table = table(brain_ics_count,powpow_rej_count,all_rej_count,powpow_rej_nums,brain_ic_nums,...
    all_rej_nums,brain_ics_scores,good_diplocs_orig,bad_diplocs_orig,good_diplocs_mni,bad_diplocs_mni,'RowNames',subj_str);
writetable(tmp_table,[study_fPath filesep 'rejection_crit.xlsx'],'WriteRowNames',true,'WriteVariableNames',true) 
fprintf('DONE. Writing subject rejection criteria table.\n');
%% REMOVE COMPS (version 1)
% (06/17/2023) JS, changing line 70 from < 3 to <= 3 (losing subjects w/ 3
% brain comps)
% tmp_rmv_subjs = zeros(1,length(ALLEEG));
% for subj_i = 1:length(ALLEEG)
%     if length(ALLEEG(subj_i).etc.urreject.ic_keep) <= 3 || isempty(ALLEEG(subj_i).etc.urreject)
%         fprintf('** Subject %s rejected.\n',ALLEEG(subj_i).subject);
%         tmp_rmv_subjs(subj_i) = 1;
%         continue;
%     end
%     ALLEEG(subj_i) = pop_subcomp(ALLEEG(subj_i),ALLEEG(subj_i).etc.urreject.ic_rej,0,0);
% end
% ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
%% REMOVE COMPS (version 2)
% (06/17/2023) JS, changing line 70 from < 3 to <= 3 (losing subjects w/ 3
% brain comps)
% (06/20/2023) JS, removing pop_subcomp from pipeline as it removes ~4
% subjects who have 2 brain components. changing to < 2.
% tmp_rmv_subjs = zeros(1,length(ALLEEG));
% TMP_ALLEEG = cell(1,length(ALLEEG));
% parfor subj_i = 1:length(ALLEEG)
%     if length(ALLEEG(subj_i).etc.urreject.ic_keep) < 2 || isempty(ALLEEG(subj_i).etc.urreject)
%         fprintf('** Subject %s rejected.\n',ALLEEG(subj_i).subject);
% %         tmp_rmv_subjs(subj_i) = 1;
%     else
%         ALLEEG(subj_i).icachansind = ALLEEG(subj_i).etc.urreject.ic_keep;
%         ALLEEG(subj_i) = eeg_checkset(ALLEEG(subj_i),'loaddata');
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
%     end
% end
% ALLEEG = ALLEEG(~logical(tmp_rmv_subjs));
% ALLEEG = TMP_ALLEEG(~cellfun(@isempty,TMP_ALLEEG));
%% CREATE STUDY
% (11/22/2023) is this step really needed? takes up a lot of space...
% initiailize study
fprintf('\n==== Making Study Modifications ====\n');
[STUDY, ALLEEG] = std_editset([],ALLEEG,...
                                'updatedat','off',...
                                'savedat','off',...
                                'name',study_fName,...
                                'filename',study_fName,...
                                'filepath',study_fPath);
% make sure all .mat files have a .fdt file associated with it.
% (03/08/23) why were these turned off? <-- to save memory!! (03/16/2023),
% use eeglab_options to set memory options so it doesn't conflict.
% (08/28/22) updatedat turnned off 
% (08/28/22) savedat turned off
parfor subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).etc.full_setfile.filename = ALLEEG(subj_i).filename;
    ALLEEG(subj_i).etc.full_setfile.filepath = ALLEEG(subj_i).filepath;
    % ALLEEG(subj_i).filename = sprintf('%s_allcond_ICA_TMPEEG.set',ALLEEG(subj_i).subject); %sprintf('%s_%s_ICA_TMPEEG',ALLEEG(subj_i).subject,'reducedcomps');
    ALLEEG(subj_i) = pop_saveset(ALLEEG(subj_i),'filename',ALLEEG(subj_i).filename,'filepath',ALLEEG(subj_i).filepath);
end
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG); 

end

