%% RECOOPING STUDY FILES AFTER ERRROR

dt = '16022023';
study_fName = sprintf('copy_study');
save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
PARSE_TYPE = 'Constant';
SUFFIX_PATH_EPOCHED = 'Epoched';
path2BEM  = [PATHS.path4EEGlab filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
mniMRI = fullfile(path2BEM, 'standard_mri.mat');
mniVol = fullfile(path2BEM, 'standard_vol.mat');
mniChan1005 = fullfile(path2BEM,'elec','standard_1005.elc');
TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];
%%
fprintf(1,'\n==== LOADING CLUSTE0R STUDY DATA ====\n');
if ~ispc
    [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',save_dir);
else
    [MAIN_STUDY,MAIN_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',save_dir);
end
MAIN_ALLEEG = eeg_checkset(MAIN_ALLEEG,'loaddata');
for subj_i = 1:length(MAIN_ALLEEG)
    if isempty(MAIN_ALLEEG(subj_i).icaact)
        fprintf('%s) Recalculating ICA activations\n',MAIN_ALLEEG(subj_i).subject);
        MAIN_ALLEEG(subj_i).icaact = (MAIN_ALLEEG(subj_i).icaweights*MAIN_ALLEEG(subj_i).icasphere)*MAIN_ALLEEG(subj_i).data(MAIN_ALLEEG(subj_i).icachansind,:);
    end
end

%%
%- assiging important dipfit model information for later recall
fprintf('==== Reassigning MRI for MNI plotting ====\n');
for subj_i = 1:length(MAIN_ALLEEG)
    try
        MAIN_ALLEEG(subj_i).dipfit.coord_transform;
        MAIN_ALLEEG(subj_i).dipfit.mrifile;
        MAIN_ALLEEG(subj_i).dipfit.hdmfile;
        MAIN_ALLEEG(subj_i).dipfit.coordformat;
    catch
        fprintf('MNI pop_dipfit_settings...\n');
        MAIN_ALLEEG(subj_i) = pop_dipfit_settings(MAIN_ALLEEG(subj_i),'coordformat','MNI','coord_transform',COORD_TRANSFORM_MNI,...
                'hdmfile',MNI_VOL,'mrifile',MNI_MRI,'chanfile',MNI_CHAN_1005);
    end
end
%% PCA REDUCTION
[MAIN_STUDY,MAIN_ALLEEG,comps_out,outliers] = cluster_pca_reduce(MAIN_STUDY,MAIN_ALLEEG);
%% CREATE STUDY FILE FOR EACH CONDITION
fprintf('==== SAVING STUDIES ====\n');
for cond_i = 1:length(TRIAL_TYPES)
    tmp_alleeg = cell(1,length(tmp));
    study_fName = sprintf('%s_MIM_study',TRIAL_TYPES{cond_i});
    
    %- extract each condition for each subject
    for subj_i = 1:length(tmp)
        if length(tmp{subj_i}) == 1
            EEG = tmp{subj_i};
        else
            EEG = tmp{subj_i}(cond_i);
        end
        tmp_alleeg{subj_i} = EEG;
    end
    tmp_alleeg = cellfun(@(x) [[]; x], tmp_alleeg);
    %## CREATE NEW STUDY STRUCTURED
    [tmp_study, tmp_alleeg] = std_editset(MAIN_STUDY,tmp_alleeg,...
                                'rmclust','off',...
                                'addchannellabels','on',...
                                'name',study_fName,...
                                'commands',{'remove',find(rmv_subj)});
    %- study modifications
    tmp_study.urcluster = MAIN_STUDY.cluster;
    tmp_study.rmvd_subj_inds = find(rmv_subj);
    tmp_study.filename = [study_fName '.study'];
    tmp_study.name = study_fName;
    tmp_study.condition = TRIAL_TYPES{cond_i};
    %## ROBUST SAVE
    [tmp_study,tmp_alleeg] = eeglab_save_study(tmp_study,tmp_alleeg,...
                                        study_fName,save_dir,...
                                        'STUDY_COND',TRIAL_TYPES{cond_i});
end
