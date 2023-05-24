%{
%## v1.2.02172023.0 : YA WITH TIME WARPING (MORE SUBJECTS)
    %## PATHS
    %- hardcode data_dir
    DATA_SET = 'MIM_dataset';
    DATA_DIR = [source_dir filesep '_data'];
    %- path for local data
    STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
    SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
    %## DATASET SPECIFIC
    %- MIND IN MOTION
    % SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
    %     'H1022','H1024','H1026','H1027','H1033','H1034'}; % CHANG,LIU(12/26/2022)
    SUBJ_YNG = {'H1002','H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1018','H1019','H1020',...
        'H1022','H1024','H1026','H1027','H1029','H1030','H1031','H1033','H1034','H1035',...
        'H1036','H1037','H1038','H1039','H1041','H1042','H1045','H1047','H1048'}; % CHANG,LIU(02/15/2023)
    TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
    %- Subject Picks
    SUBJ_PICS = {SUBJ_YNG};
    SUBJ_ITERS = {1:length(SUBJ_YNG)}; % CHANG,LIU(12/26/2023)
    % SUBJ_ITERS = {[],1:length(SUBJ_HMA),1:length(SUBJ_NMA)}; % CHANG,LIU(02/15/2023)
    %- Subject Directory Information
    PREPROCESS_NAME = 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10'; % CHANG,LIU(02/15/2023)
    OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-HY_202212';  % CHANG,LIU(02/15/2023)
    % PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % CHANG,LIU(12/26/2022)
    % OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % CHANG,LIU(12/26/2022)
    if DO_UNIX
        OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR);
    else
        OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR);
    end
    % OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
    OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping
    %% ===================================================================== %%
    %## PROCESSING PARAMS
    % pop_editoptions('option_parallel',1);
    %- study group and saving
    SAVE_EEG = false; % saves EEG structures throughout processing
    %- component rejection crit
    THRESHOLD_DIPFIT_RV = 0.15;
    THRESHOLD_BRAIN_ICLABEL = 0.50;
    %- MIM specific epoching
    EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
    %- subj_i_epoch
    PARSE_TYPE = 'Constant'; %
    EPOCH_TIME_LIMITS = [-1,1];
    TRIAL_LENGTH = 3*60; % trial length in seconds
    PER_OVERLAP = 0.0; % percent overlap between epochs
    TRIAL_BEGIN_STR = 'TrialStart';
    TRIAL_END_STR = 'TrialEnd';
    EVENT_TRIAL_PARSER = 'type';
    EVENT_COND_PARSER = 'cond';
    %- connectivity process
    CONN_FREQS = (1:100);
    CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
    %- subj_i_genConnMeas
    CNCTANL_TOOLBOX = 'sift'; %'bsmart'
    WINDOW_LENGTH = 0.5;
    WINDOW_STEP_SIZE = 0.025;
    NEW_SAMPLE_RATE = [];
    DO_BOOTSTRAP = true;
    % ASSIGN_BOOTSTRAP_MEAN = false;
    % SAVE_CONN_BOOTSTRAP = false;
    %- subj_i_genConnStats
    DO_PHASE_RND = true;
    %- eeglab_cluster.m spectral params
    FREQ_LIMITS = [1,100];
    CYCLE_LIMITS = [3,0.8];
    SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
    FREQ_FAC = 4;
    PAD_RATIO = 2;
%## v1.1.02152023.0 : YA WTIH TIME WARPING
    %## PATHS
    %- hardcode data_dir
    DATA_SET = 'MIM_dataset';
    DATA_DIR = [source_dir filesep '_data'];
    %- path for local data
    STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
    SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
    %## DATASET SPECIFIC
    %- MIND IN MOTION
    SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
        'H1022','H1024','H1026','H1027','H1033','H1034'}; % CHANG,LIU(12/26/2022)
    TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
    %- Subject Picks
    SUBJ_ITERS = {1:length(SUBJ_YNG),[],[]}; % CHANG,LIU(12/26/2023)
    %- Subject Directory Information
    PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % CHANG,LIU(12/26/2022)
    OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % CHANG,LIU(12/26/2022)
    if DO_UNIX
        OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR);
    else
        OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR);
    end
    % OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
    OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping
    %% ===================================================================== %%
    %## PROCESSING PARAMS
    %- study group and saving
    SAVE_EEG = false; % saves EEG structures throughout processing
    %- component rejection crit
    THRESHOLD_DIPFIT_RV = 0.15;
    THRESHOLD_BRAIN_ICLABEL = 0.50;
    %- MIM specific epoching
    EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
    %- subj_i_epoch
    PARSE_TYPE = 'Constant'; %
    EPOCH_TIME_LIMITS = [-1,1];
    TRIAL_LENGTH = 3*60; % trial length in seconds
    PER_OVERLAP = 0.0; % percent overlap between epochs
    TRIAL_BEGIN_STR = 'TrialStart';
    TRIAL_END_STR = 'TrialEnd';
    EVENT_TRIAL_PARSER = 'type';
    EVENT_COND_PARSER = 'cond';
    %- connectivity process
    CONN_FREQS = (1:100);
    CONN_METHODS = {'dDTF','GGC','dDTF08'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
    %- subj_i_genConnMeas
    CNCTANL_TOOLBOX = 'sift'; %'bsmart'
    WINDOW_LENGTH = 0.5;
    WINDOW_STEP_SIZE = 0.025;
    NEW_SAMPLE_RATE = [];
    DO_BOOTSTRAP = false;
    % ASSIGN_BOOTSTRAP_MEAN = false;
    % SAVE_CONN_BOOTSTRAP = false;
    %- subj_i_genConnStats
    DO_PHASE_RND = true;
    %- eeglab_cluster.m spectral params
    FREQ_LIMITS = [1,100];
    CYCLE_LIMITS = [3,0.8];
    SPEC_MODE = 'fft'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
    FREQ_FAC = 4;
    PAD_RATIO = 2;
%## v1.0.01132023.0 : Initializing versioning for future iterations. Previous
    %## PATHS
    %- hardcode data_dir
    DATA_SET = 'MIM_dataset';
    DATA_DIR = [source_dir filesep '_data'];
    %- path for local data
    STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
    SUBJINF_DIR = [DATA_DIR filesep DATA_SET filesep '_subjinf'];
    %- MIND IN MOTION
    SUBJ_YNG = {'H1004','H1007','H1009','H1010','H1011','H1012','H1013','H1017','H1020',...
        'H1022','H1024','H1026','H1027','H1033','H1034'};
    SUBJ_HMA = {'H2002', 'H2010', 'H2015', 'H2017', 'H2020', 'H2021', 'H2022', 'H2023',...
        'H2025', 'H2026', 'H2034', 'H2059', 'H2062', 'H2082', 'H2095'};
    SUBJ_NMA = {'NH3008', 'NH3043', 'NH3055', 'NH3059', 'NH3069', ...
        'NH3070', 'NH3074', 'NH3086', 'NH3090', 'NH3104', 'NH3105', 'NH3106', 'NH3112', 'NH3114'};
    TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
    %- Subject Picks
    SUBJ_PICS = {SUBJ_YNG,SUBJ_HMA,SUBJ_NMA};
    SUBJ_ITERS = {1:length(SUBJ_YNG),[],[]};
    %- Subject Directory Information
    PREPROCESS_NAME = '8std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'; % value to help with file looping
    OUTSIDE_DATA_DIR = 'M:\liu.chang1\STUDY-preprocess-MiM_20220725'; % value to help with file looping
    if DO_UNIX
        OUTSIDE_DATA_DIR = convertPath2UNIX(OUTSIDE_DATA_DIR,'dferris');
    else
        OUTSIDE_DATA_DIR = convertPath2Drive(OUTSIDE_DATA_DIR,'M');
    end
    % OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_templateElec';
    OUTSIDE_DATA_SUFFIX = 'ICA_ICLabel_dipfit_fem'; % value to help with file looping
    %## Params
    dt = '26122022';
    %## PROCESSING PARAMS
    % pop_editoptions('option_parallel',1);
    %- subjstruct and saving
    SAVE_EEG = false; % saves EEG structures throughout processing
    %- component rejection crit
    THRESHOLD_DIPFIT_RV = 0.15;
    THRESHOLD_BRAIN_ICLABEL = 0.50;
    %- subj_i_epoch
    PARSE_TYPE = 'Constant'; %
    EPOCH_TIME_LIMITS = [-2,2];
    TRIAL_LENGTH = 3*60; % trial length in seconds
    PER_OVERLAP = 0.0; % percent overlap between epochs
    TRIAL_BEGIN_STR = 'TrialStart';
    TRIAL_END_STR = 'TrialEnd';
    EVENT_TRIAL_PARSER = 'type';
    EVENT_COND_PARSER = 'cond';
    %- connectivity process
    FREQS = (3:60);
    CONN_METHODS = {'dDTF','GGC','ffDTF'}; % Options: 'S', 'dDTF08', 'GGC', 'mCoh', 'iCoh'
    % ANAT_THRESH = 30; % distance to anatomical location in (mm)
    %- subj_i_genConnMeas
    CNCTANL_TOOLBOX = 'sift'; %'bsmart'
    WINDOW_LENGTH = 0.5;
    WINDOW_STEP_SIZE = 0.025;
    NEW_SAMPLE_RATE = [];
    DO_BOOTSTRAP = false;
    ASSIGN_BOOTSTRAP_MEAN = false;
    SAVE_CONN_BOOTSTRAP = false;
    %- subj_i_genConnStats
    DO_PHASE_RND = true;
    %- eeglab_cluster.m spectral params
    FREQ_LIMITS = [1,100];
    CYCLE_LIMITS = [3,0.8];
    SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
    FREQ_FAC = 4;
    PAD_RATIO = 2;
%}
