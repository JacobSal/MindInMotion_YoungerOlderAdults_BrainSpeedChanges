%{
%## v1.0.01132023.0 : Initializing versioning for future iterations. Previous
    Params:
        %## PATHS
        %- hardcode data_dir
        DATA_SET = 'jacobsenN_dataset';
        DATA_DIR = [source_dir filesep '_data'];
        %- path for local data
        STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
        SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
        %## DATASET SPECIFIC
        SUBJ_PICS = {'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34',...
                             'S35','S36','S37','S38','S39','S40','S41','S42','S43','S44',...
                             'S46','S47','S48'};
        DATASET_PROCESSNUM = 2;
        TRIAL_TYPES = {'pre','post'};
        %- hardcode data_dir
        DATA_SET = 'jacobsenN_dataset';
        DATA_DIR = [source_dir filesep '_data'];
        %- path for local data
        STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
        SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
        %## DATASET SPECIFIC
        SUBJ_PICS = {'S25','S26','S27','S28','S29','S30','S31','S32','S33','S34',...
                             'S35','S36','S37','S38','S39','S40','S41','S42','S43','S44',...
                             'S46','S47','S48'};
        DATASET_PROCESSNUM = 2;
        TRIAL_TYPES = {'pre','post'};
        %## PROCESSING PARAMS
        %- subjstruct and saving
        SAVE_EEG = false;
        REGEN_ALL_SUBJSTRUCT = false;
        %## ==== SCRIPT PARAMS ==== ##%
        %- component rejection
        DO_DIPOLE_REJ = true;
        THRESHOLD_DIPFIT_RV = 0.15;
        THRESHOLD_BRAIN_ICLABEL = 0.65;
        %- subj_i_epoch
        PARSE_TYPE = 'Constant'; %
        EPOCH_TIME_LIMITS = [-2,2]; 
        TRIAL_LENGTH = 3*60; % trial length in seconds
        PER_OVERLAP = 0.0; % percent overlap between epochs
        TRIAL_BEGIN_STR = '';
        TRIAL_END_STR = '';
        EVENT_TRIAL_PARSER = '';
        EVENT_COND_PARSER = '';
        %- connectivity process
        FREQS = (3:60);
        CONN_METHODS = {'dDTF','GGC','ffDTF'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
        % ANAT_THRESH = 30; % distance to anatomical location in (mm)
        % subj_i_genConnMeas
        CNCTANL_TOOLBOX = 'sift'; %'bsmart'
        DO_BOOTSTRAP = false;
        WINDOW_LENGTH = 0.5;
        WINDOW_STEP_SIZE = 0.025;
        NEW_SAMPLE_RATE = [];
        ASSIGN_BOOTSTRAP_MEAN = false;
        SAVE_CONN_BOOTSTRAP = false;
        % subj_i_genConnStats
        DO_PHASE_RND = false;
%## v1.1.
%}
