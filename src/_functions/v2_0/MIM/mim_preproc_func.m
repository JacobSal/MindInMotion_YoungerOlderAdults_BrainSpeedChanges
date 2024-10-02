function [EEG,amica_cmd,mim_config] = mim_preproc_func(subj_char,save_fPath,fOutput,fDataset,varargin)
%MAIN_FUNC Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%
% NOTES:
%       Code designed for Mind In Motion Study (NIHU01) ran from the
%       University of Florida.
%       prep step 1: High pass filter, cleanline, and merge
%       
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 02/23/2023, MATLAB 2019a
% Copyright (C) Chang Liu, Designer & Creator (20210820)
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu, Code Dweeb
% 
% Modified - 20230205 - add EMG high filter for ICC (Chang Liu)
% Modified - 20230223 - function wrapper (Jacob Salminen)

%## TIME
tic
%## DEFINE DEFAULTS
%- IC_Rejection excel
MiM_config = [];
STUDY_FOLDER = [];
if ~isempty(STUDY_FOLDER)
    disp('Loading Powpowcat...');
    % Specify sheet and range
    switch Study_folder
        case '5std_iCC0p65_ChanRej0p5_TimeRej0p4_winTol10'
            sheet_name = "Sheet1";
        case '5std_iCC0p9_ChanRej0p5_TimeRej0p4_winTol10'
            sheet_name = "Sheet2";
        case 'EMG_HP5std_iCC0p65_iCCEMG0p4_ChanRej0p7_TimeRej0p4_winTol10'
            sheet_name = "Sheet3";    
    end
    opts = spreadsheetImportOptions("NumVariables", 8);
    opts.Sheet = char(sheet_name);
    opts.DataRange = "A2:H36";
    % Specify column names and types
    opts.VariableNames = ["num", "SubjCode", "PowpowcatIC_powpowDefinitelyBad", "PowpowcatIC_powpowMaybeBad", "PowpowcatIC_numberDefinitelyBad", "PowpowcatIC_numberMaybeBad", "Status","Number_of_components"];
    opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string","string"];
    % Specify variable properties
    opts = setvaropts(opts, ["SubjCode", "PowpowcatIC_powpowDefinitelyBad", "PowpowcatIC_powpowMaybeBad", "PowpowcatIC_numberDefinitelyBad", "PowpowcatIC_numberMaybeBad", "Status","Number_of_components"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["SubjCode", "PowpowcatIC_powpowDefinitelyBad", "PowpowcatIC_powpowMaybeBad", "PowpowcatIC_numberDefinitelyBad", "PowpowcatIC_numberMaybeBad", "Status","Number_of_components"], "EmptyFieldRule", "auto");
    % Import the data
    PowpowCat_Rej = readtable([fDataset filesep 'subject_mgmt' filesep 'PowpowCat_Rejection_HY.xlsx'], opts, "UseExcel", false);
    PowpowCat_Rej.Properties.VariableNames = {'Num','SubjStr','No_use_bad','No_use_maybebad','IC_bad','IC_maybebad','Status','Number_components'};
    MiM_config.PowpowCat_Rej = PowpowCat_Rej;
end
%##
fullRankAvRefBool = false;
doCleanline = 1;
% overwrite_file = 0;                              
do_processEMG = 1;
%## folder settings
mim_config.shareFolder      = 'R:\Ferris-Lab\share\MindInMotion\Data\';
% mim_config.EEG_input        = save_fPath;
mim_config.subjMgmtFolder   = [fDataset filesep 'subject_mgmt' filesep 'Trial_Cropping_V2_test.xlsx'];
%## Pre-processing parameters
%- 
mim_config.fullRankAvRefBool = false; %default false for regular avg ref, true for full rank version
mim_config.avgRefPCAReduction = 1; %1 for just eeg, 3 for eeg+noise+emg % what is the difference?
mim_config.std_threshold = 3;
%- ICanClean parameters
mim_config.ICC.params=[];
mim_config.ICC.params.statsEEG_stande = 4;
mim_config.ICC.params.cleanWindow = 2;
mim_config.ICC.params.RTBool = false;
% mim_config.ICC.params.windowLength = 2; % 2023-1-10 Ryan said 2 is better % previously used 1; % sliding window that gets cleaned 
mim_config.ICC.params.rhoSqThres_source = 0.65; % 2023-1-10 Use 0.65 for 4sec window %Default 0.9, NH3066 try 0.8, H3063 = 0.95, NH3068, H3072 try 0.85
%* (07/11/2023) JS, trying rhoSqThres_source 0.65
%* (07/12/2023) JS, trying rhoSqThres_source 0.9 as a less aggressive cutoff 
%* (10/30/2023) JS, trying rhoSqThres_source 0.6 as a  more aggressive cutoff
mim_config.ICC.params.cleanXwith = 'X'; 
%-- old parameterization (05/05/2023)
% mim_config.ICC.params.extraTime_pre = 1; % 2023-1-10 Ryan said use 1 (2-MiM_config.ICC.params.windowLength)/2; %wider window for stats
% mim_config.ICC.params.extraTime_post = 1; %MiM_config.ICC.params.extraTime_pre; %wider window for stats
%- ICanClean Muscle parameters
mim_config.ICC_muscle.params = mim_config.ICC.params;
mim_config.ICC_muscle.params.rhoSqThres_source = 0.4; %with 4sec window, 0.9 is too high, 0 bad sources get removed...- -(EMG without highpass filter)
%* (07/11/2023) JS, trying rhoSqThres_source 0.4
%* (10/26/2023) JS, trying rhoSqThres_source 0.3  as a more aggressive cutoff 
%* (10/30/2023) JS, trying rhoSqThres_source 0.25  as a more aggressive cutoff 
%- ICanClean Canonical Correlation Analysis parameters
lagAmount_samples = 1;
CCA_Rsq_thres = 0.2; %default 0.2
mim_config.CCA.lagAmount_samples = lagAmount_samples;
mim_config.CCA.CCA_Rsq_thres = CCA_Rsq_thres;
%- Artifact Subspace Rejection parameters
mim_config.ASR_correct.burstCritASR = 10;
mim_config.ASR_correct.useExternalCalibASR = 1;
%- CleanLine Parameters
mim_config.chan_crit1 = 0.70; %used 0.5 for young participants, 0.7 for H1044 %default = 0.85; Some notes: used 0.7 for all MiM participant after iCC (!), except for H1026,H1030 used 0.5 %the default 0.85 seems to be too aggressive to my data. 
mim_config.wind_crit = 0.4; %default = 0.25
mim_config.chan_max_broken_time = 0.6;
mim_config.chan_detected_fraction_threshold = 0.6;
mim_config.flatline_crit = 'off';
mim_config.line_noise_crit = 'off';
mim_config.kurt_crit = 5;%default = 5
mim_config.chan_detect_num_iter = 5;
mim_config.window_crit_tolerances = [-Inf, 10]; %default = 7
mim_config.kurt_prob_std_rm = 0;
%- further filtering settings
mim_config.filter_lowCutoffFreqAMICA = 1.75;% 1.75 is 1.5Hz cutoff!
mim_config.filter_AMICA_highPassOrder = 1650;
mim_config.filter_highCutoffFreqAMICA = [];
mim_config.filter_AMICA_lowPassOrder = [];
%- AMICA Parameters
mim_config.emailStr = 'jsalminen@ufl.edu';
%## Script Specific Params
cushion_sec = 5; %sec
%## PREPROCESSING PARAMS
EMG_highpass = 1; %use .set with highpass filter EMG
do_iCC = 0;
do_iCC_and_ChanRej_TimeRej = 0; %This is the pipeline we decide to use - 2022-09-02
do_iCC_muscle = 0;
do_iCC_and_ChanRej_TimeRej_iCCMuscle = 1; %Testing this pipeline now - 2023-02-05
do_CCA = 0;
do_ASR = 0;
do_ChanRej = 0;
do_TimeRej = 0; %TimeRej include ChanRej. Pick ChanRej or TimeRej
do_HP = {}; %{'8std_ChanRej'};%{'8std_ChanRej_TimeRej'};
do_ChanRej_iCC = 0;
do_postASR_ChanRej = {}; %
do_postASR_ChanRej_TimeRej = {};
autoChRejMethod = [num2str(mim_config.std_threshold),'std'];
cleaningMethod = '';
finalChToKeepForICA = {'EEG'};
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'subj_char',@ischar);
addRequired(p,'fPath',@ischar);
addRequired(p,'fOutput',@ischar);
addRequired(p,'fDataset',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,subj_char,save_fPath,fOutput,fDataset,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETERS
%- MIMCONFIG
mim_config.CleanEEG_output     = [fOutput filesep 'clean'];
mim_config.amicaEEG_output     = [fOutput filesep 'clean'];
%% ===================================================================== %%
%## (SCRIPT PART 1)
EEG = [];
ALLEEG = [];
%## Find All Trials of Interest
fileList_TM = dir([save_fPath filesep 'TM*.set']);
fileList_SP = dir([save_fPath filesep 'SP*.set']);
fileList_Rest = dir([save_fPath filesep 'Rest.set']);
fileList_Imagined = dir([save_fPath filesep 'MotorImagery*.set']);
fileList = [fileList_Rest; fileList_TM; fileList_SP];% for this analysis, only include Rest, TM, and SP
if length(fileList) < 10
    error('Subject %s missing condition',subj_char)
end
nTrials = 1:size(fileList,1);
%% Trial Level Processing
% The processing includess a highpass filter cutoff at 1Hz, and then use
% cleanline at 60Hz for all channels
% eeglab; %initialize EEG and ALLEEG as empty
for trial_i = nTrials
    %Load trial
    EEG = pop_loadset('filename',fileList(trial_i).name,'filepath',save_fPath);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); %append to new EEG to ALLEEG

    %% Deal with trial cropping
    fprintf('Cropping trial %i for subject %s',trial_i,subj_char)
%     [DoCrop, ExactCrop, DoCrop_loadsol, ExactCrop_loadsol ] = CropTrialCheckFunc_checkLoadsol(subj_char,fileList(trial_i).name(1:end-4),...
%         mim_config.subjMgmtFolder);
    [ DoCrop, ExactCrop, DoCrop_loadsol, ExactCrop_loadsol ] = mim_check_trials(subj_char,fileList(trial_i).name(1:end-4),...
        mim_config.subjMgmtFolder);
    if DoCrop_loadsol % 
        % remove gait events not in the ExactCrop_loadsol range
        % first convert the time --> latency
%             temp = [EEG.event.latency];
        disp('You are removing bad gait events');
        cushion_latency = EEG.srate*cushion_sec;
        crop_latency_range = EEG.srate*ExactCrop_loadsol + cushion_latency + 1;
        total_gaitevent_idx = find(strcmp('LHS',{EEG.event.type}) | strcmp('LTO',{EEG.event.type})|...
            strcmp('RHS',{EEG.event.type}) | strcmp('RTO',{EEG.event.type}));
        gaitevent_idx_inrange = [];
        for p = 1:size(crop_latency_range,1)
            gaitevent_idx_inrange = horzcat(gaitevent_idx_inrange,find( [EEG.event.latency] > ...
                crop_latency_range(p,1) & [EEG.event.latency] < crop_latency_range(p,2)));
        end
        [gaitevent_idx_outrange] = setdiff(total_gaitevent_idx,gaitevent_idx_inrange);
        EEG.event(gaitevent_idx_outrange) = [];

%             figure();plot(temp(gaitevent_idx_outrange)./500,ones(length(gaitevent_idx_outrange)) ,'o');hold on;plot(temp./500,zeros(length(temp)),'o');
%             figure();plot(temp(gaitevent_idx_inrange)./500,ones(length(gaitevent_idx_inrange)) ,'o');
    end
    if DoCrop %only deal with shortened trials for now (cannot ged rid of other events that happen during trial)
        disp('You are cropping the trial');
        startInd = find(strcmpi('TrialStart',{EEG.event.type}));
        startTime = EEG.times( round(EEG.event(startInd).latency) )/1000;
        endInd = find(strcmpi('TrialEnd',{EEG.event.type}));
        ExactCropLatencies = EEG.srate*ExactCrop(end)+EEG.event(startInd).latency;        
        endTime = EEG.times( round(ExactCropLatencies))/1000;
        %- add trial end
        EEG.event(end+1) = EEG.event(endInd);
        EEG.event(end).latency = ExactCropLatencies;
        EEG.event(end).datetime = [];

        cushion_sec = 5;
        ExactCrop_update = ExactCrop+cushion_sec;
        ExactCrop_update(1) = 0;
        ExactCrop_update(end,end) = endTime+cushion_sec;
        EEG = pop_select( EEG, 'time',ExactCrop_update); % This is wrong in the original code from Roehl. That code doesn't take into account the cushion time added for each trial
    end

    %## HP filter EEG, EMG, and Noise
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'chantype',{'EEG' 'EMG' 'Noise'}); %1 Hz HP filter (results in -6bB cutoff at 0.5 Hz)
    if do_processEMG
        EEG = pop_eegfiltnew(EEG, 'locutoff',20,'chantype',{ 'EMG' }); %high pass filter 20 Hz
    end
    %- Keep only a subset of chans
    EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %define channels
    EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
    EEG = pop_select( EEG, 'channel', sort([EEG_chans, EMG_chans,Noise_chans])); %select subset of channels
    %     EEG = pop_select( EEG, 'channel', sort([EEG_chans])); %select just EEG channels

    %- update channel labels again just in case?
    EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %redefine channels (numbering changed since we took subset of channels)
    EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
    Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));

    %## Cleanline
    if doCleanline % removal of line noise (set to remove 60 Hz below)
        % Reference before cleanline filter
        EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
        EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[EEG_chans EMG_chans Noise_chans],'computepower',1,...
            'linefreqs',[60 120],'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',...
            0,'scanforlines',1,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
        [ALLEEG,~,~] = eeg_store(ALLEEG, EEG, CURRENTSET); %These params followed PREP pipeline
    else
        %Update ALLEEG with whatever changes you made to EEG variable before moving on to next trial
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET); %overwrites current set    
    end
end
%## Session level Processing
EEG = pop_mergeset(ALLEEG,1:length(ALLEEG), 0);
EEG.setname = sprintf('%s_merged_eeg',subj_char);
EEG.filename = sprintf('%s_merged_raw',subj_char);
[~, EEG, ~] = eeg_store( ALLEEG, EEG, 0 ); %append to new EEG to ALLEEG

% eeglab redraw
% % TO DO: SAVE the MERGED EEG so we don't need to start from beginning
% if doCleanline && do_processEMG
%     fName = strcat(subj_char,'_HP1hz_cleanline_merge_EMG');
%     [EEG] = pop_saveset(EEG,'filename',fName,...
%         'filepath',fOutput,...
%         'savemode','twofiles');
% elseif doCleanline && do_processEMG == 0
%     fName = strcat(subj_char,'_HP1hz_cleanline_merge');
%     % TO DO: SAVE the MERGED EEG so we don't need to start from beginning
%     [EEG] = pop_saveset(EEG,'filename',fName,...
%         'filepath',fOutput,...
%         'savemode','twofiles');
% else
%     fName = strcat(subj_char,'_HP1hz_merge');
%     [EEG] = pop_saveset(EEG,'filename',fName,...
%         'filepath',fOutput,...
%         'savemode','twofiles');        
% end
toc
%%

%% ===================================================================== %%
%## SCRIPT PART 2
tic
fprintf('==== COMPLETING %s PREPROCESSING ====\n\n',subj_char) 
if ischar(finalChToKeepForICA) %e.g. just 'EEG', not {'EEG'} or {'EEG','Noise'}
    finalChToKeepForICA = {finalChToKeepForICA}; %turn 'EEG' into {'EEG'} to make code below more robust
end
if fullRankAvRefBool
    avg_ref_pca_reduction = 0;
else
    avg_ref_pca_reduction = length(finalChToKeepForICA); %1 for just {'EEG'} and 3 for {'EEG','Noise','EMG'}
end

%% Load subject filtered merged file
if EMG_highpass
    fileName = [subj_char,'_HP1hz_cleanline_merge_EMG.set']; % edit here for file name
    cleaningMethod = horzcat(cleaningMethod,'EMG_HP');
else
    fileName = [subj_char,'_HP1hz_cleanline_merge.set']; % edit here for file name
end
% EEG = pop_loadset('filepath',mim_config.EEG_input,'filename',fileName);
EEG.urchanlocs = EEG.chanlocs; % keep one copy of old channel info; needed for chan_rej
%% Session Level Processing
EEG.setname = sprintf('%s_merged_eeg',subj_char);
EEG.filename = sprintf('%s_merged_raw',subj_char);
% eeglab redraw;
EEG.subject = subj_char;
%## preprocessing steps
%Re-ref EEG, EMG, and Noise to themselves
fprintf('==== %s 1.) REREFERENCING ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
%Reject bad channels
fprintf('==== %s 2.) AUTOREJECTION OF BAD CHANS ====\n',subj_char)
EEG = autoRejCh_func_CL(EEG,mim_config.std_threshold);
%Re-ref again (since remnants of rejected channels still exist from
%last average referencing
fprintf('==== %s 3.) REREFERENCING pt.2 ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
% reject and re-ref again
fprintf('==== %s 4.) AUTOREJECTION OF BAD CHANS pt.2 ====\n',subj_char)
EEG = autoRejCh_func_CL(EEG,mim_config.std_threshold);
% re-ref
fprintf('==== %s 5.) REREFERENCING pt.3 ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
cleaningMethod = horzcat(cleaningMethod,autoChRejMethod);
%update channel types b/c that is always a good idea
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type})); %redefine channels (numbering changed since we took subset of channels)
EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
%% CLEANING    
%## ICanClean Base
if do_iCC
    rho = num2str(mim_config.ICC.params.rhoSqThres_source);idx_rho = strfind(rho,'.');
    cleaningMethod = horzcat(cleaningMethod,'_iCC','0p',rho(idx_rho+1:end));
    fprintf('==== %s STARTING: ICANCLEAN ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,[EEG_chans], [Noise_chans], 0, mim_config.ICC.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end
%## ICanClean w/ noise chans && Chan Rej and Time Rej
if do_iCC_and_ChanRej_TimeRej % Roehl look here
    rho = num2str(mim_config.ICC.params.rhoSqThres_source);idx_rho = strfind(rho,'.');
    cleaningMethod = horzcat(cleaningMethod,'_iCC','0p',rho(idx_rho+1:end));
    fprintf('==== %s STARTING: ICANCLEAN ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,[EEG_chans], [Noise_chans], 0, mim_config.ICC.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
    if mim_config.kurt_prob_std_rm == 0
        kurt_prob_std_rm = '';
    else
        kurt_prob_std_rm = '_MoreChanRej';
    end
    winParam = num2str(mim_config.wind_crit);chanParam = num2str(mim_config.chan_crit1);winTol = num2str(mim_config.window_crit_tolerances(2));
    fprintf('==== %s STARTING: CHAN TIME REJ ====\n',subj_char);       
    [EEG_temp_clean,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,mim_config); 
    %- Time Rejection
    EEG = EEG_temp_clean_timerej;zthreshold = EEG.etc.clean_artifacts.window_zthreshold;
    cleaningMethod = horzcat(cleaningMethod,'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_TimeRej','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol,kurt_prob_std_rm);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
    preprocess_pipeline = [cleaningMethod];
    %- save output
    mkdir([mim_config.CleanEEG_output])
    fid = fopen([mim_config.CleanEEG_output filesep 'info.txt'],'w');
    fprintf(fid,'\n %.2f percent of frames were rejected\n', p_frames_rej);
    fprintf(fid,'\n %.2f channels were rejected\n', p_chan_rej);
    fclose(fid);
    save([mim_config.CleanEEG_output filesep 'rejection_info.mat'],'p_frames_rej','p_chan_rej','zthreshold');
end

%##
if do_iCC_muscle
    rho_iCC_muscle = num2str(mim_config.ICC_muscle.params.rhoSqThres_source);
    cleaningMethod = horzcat(cleaningMethod,'_iCCEMG','0p',rho_iCC_muscle(strfind(rho_iCC_muscle,'.')+1:end));
    fprintf('==== %s STARTING: ICANCLEAN MUSCLE ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,[EEG_chans], [EMG_chans], 0, mim_config.ICC_muscle.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end

%## ICanClean w/ noise Sensors and EMG sensors && Channel Rejection with Time Rejection
if do_iCC_and_ChanRej_TimeRej_iCCMuscle
    rho = num2str(mim_config.ICC.params.rhoSqThres_source);idx_rho = strfind(rho,'.');
    cleaningMethod = horzcat(cleaningMethod,'_iCC','0p',rho(idx_rho+1:end));
    fprintf('==== %s STARTING: ICANCLEAN EEG with Noise ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,[EEG_chans], [Noise_chans], 0, mim_config.ICC.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);

    rho_iCC_muscle = num2str(mim_config.ICC_muscle.params.rhoSqThres_source);
    cleaningMethod = horzcat(cleaningMethod,'_iCCEMG','0p',rho_iCC_muscle(strfind(rho_iCC_muscle,'.')+1:end));
    fprintf('==== %s STARTING: ICANCLEAN EEG with EMG ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG,[EEG_chans], [EMG_chans], 0, mim_config.ICC_muscle.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);

    if mim_config.kurt_prob_std_rm == 0 ; kurt_prob_std_rm = '';else;kurt_prob_std_rm = '_MoreChanRej';end;
    winParam = num2str(mim_config.wind_crit);chanParam = num2str(mim_config.chan_crit1);winTol = num2str(mim_config.window_crit_tolerances(2));
    fprintf('==== %s STARTING: CHAN TIME REJ ====\n',subj_char);       
    [EEG_temp_clean,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,mim_config); 
    %Time Rej
    EEG = EEG_temp_clean_timerej;zthreshold = EEG.etc.clean_artifacts.window_zthreshold;
    cleaningMethod = horzcat(cleaningMethod,'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_TimeRej','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol,kurt_prob_std_rm);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
    preprocess_pipeline = [cleaningMethod];
    %- save output
    mkdir([mim_config.CleanEEG_output])
    fid = fopen([mim_config.CleanEEG_output filesep 'info.txt'],'w');
    fprintf(fid,'\n %.2f percent of frames were rejected\n', p_frames_rej);
    fprintf(fid,'\n %.2f channels were rejected\n', p_chan_rej);
    fclose(fid);
    save([mim_config.CleanEEG_output filesep 'rejection_info.mat'],'p_frames_rej','p_chan_rej','zthreshold');
end

%## CCA
if do_CCA
    rho_CCA = num2str(mim_config.CCA.CCA_Rsq_thres);idx_rho_CCA = strfind(rho_CCA,'.');
    cleaningMethod = horzcat(cleaningMethod,'_CCA','0p',rho_CCA(idx_rho_CCA+1:end));
    fprintf('==== %s STARTING: CCA ====\n',subj_char);
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG)        
    EEG = autoLagCCA_wrap(EEG,mim_config.CCA);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end

%## ASR
if do_ASR
    BurstCriteria = num2str(mim_config.ASR_correct.burstCritASR);
    cleaningMethod = horzcat(cleaningMethod,'_ASRcorr',BurstCriteria);
    fprintf('==== %s STARTING: AUTOMATIC SUBSPACE RECONSTRUCTION ====\n',subj_char);
    [EEG] = ASR_correction(EEG,mim_config.ASR_correct);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end

%ChanRej
if do_ChanRej || do_TimeRej
    if mim_config.kurt_prob_std_rm == 0 ; kurt_prob_std_rm = '';else;kurt_prob_std_rm = '_MoreChanRej';end;
    winParam = num2str(mim_config.wind_crit);chanParam = num2str(mim_config.chan_crit1);winTol = num2str(mim_config.window_crit_tolerances(2));
    fprintf('==== %s STARTING: CHAN TIME REJ ====\n',subj_char);       
    [EEG_temp_clean,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,mim_config); 
    %Time Rej
    if do_ChanRej
        cleaningMethod = horzcat(cleaningMethod,'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol,kurt_prob_std_rm);
        EEG = EEG_temp_clean;
        EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
    end
    if do_TimeRej
        if p_frames_rej < 10 % if time window rejected under 10%, rejection
            Rej_time_window = 1;
        else
            Rej_time_window = 0;
        end
        if Rej_time_window 
            EEG = EEG_temp_clean_timerej;zthreshold = EEG.etc.clean_artifacts.window_zthreshold;
            cleaningMethod = horzcat(cleaningMethod,'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_TimeRej','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol,kurt_prob_std_rm);
            EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
        else
            EEG = EEG_temp_clean;
            cleaningMethod = horzcat(cleaningMethod,'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_TimeRej','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol,kurt_prob_std_rm);
            subj_char = horzcat(subj_char,'_noTimeRej');
            EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
        end
    end
    preprocess_pipeline = [cleaningMethod];
    %- save output
    mkdir([mim_config.CleanEEG_output])
    fid = fopen([mim_config.CleanEEG_output filesep 'info.txt'],'w');
    fprintf(fid,'\n %.2f percent of frames were rejected\n', p_frames_rej);
    fprintf(fid,'\n %.2f channels were rejected\n', p_chan_rej);
    fclose(fid);
    save([mim_config.CleanEEG_output filesep 'rejection_info.mat'],'p_frames_rej','p_chan_rej','zthreshold');
end

if ~isempty(do_HP)
    fprintf('==== %s STARTING: HIGH PASS FILTERING ====\n',subj_char);       
    cleaningMethod = horzcat(do_HP{1},'_HP');
    EEGpreprocess_fileName = [subj_char,'_cleanEEG_',do_HP{1},'.set'];
    EEGpreprocess_filePath = fullfile(mim_config.CleanEEG_output,do_HP{1},subj_char);
    EEG_preprocessed = pop_loadset('filepath',EEGpreprocess_filePath,'filename',EEGpreprocess_fileName);
    [EEG] = bemobil_filter_CL(EEG_preprocessed,...
        mim_config.filter_lowCutoffFreqAMICA, mim_config.filter_highCutoffFreqAMICA,...
        mim_config.filter_AMICA_highPassOrder, mim_config.filter_AMICA_lowPassOrder);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end

if do_ChanRej_iCC 
    cleaningMethod = horzcat(cleaningMethod,'_ChanRej_iCC');
    [EEG_temp_clean,EEG_temp_clean_timerej,p_frames_rej,p_chan_rej] = channelrejection_wrap(EEG,mim_config); 
    [EEG_chans, EMG_chans, Noise_chans] = getChannelTypes_func(EEG) ;
    [EEG] = iCanClean(EEG_temp_clean,[EEG_chans], [Noise_chans], 0, mim_config.ICC.params);
    EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
end

if ~isempty(do_postASR_ChanRej)
    winParam = num2str(mim_config.wind_crit);chanParam = num2str(mim_config.chan_crit1);winTol = num2str(mim_config.window_crit_tolerances(2));
    cleaningMethod = horzcat(do_postASR_ChanRej{1},'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol);
    [EEG] = do_postASR_ChanRej_wrap(EEG,subj_char,cleaningMethod, do_postASR_ChanRej,mim_config,fullRankAvRefBool);
end

if ~isempty(do_postASR_ChanRej_TimeRej)
    winParam = num2str(mim_config.wind_crit);chanParam = num2str(mim_config.chan_crit1);winTol = num2str(mim_config.window_crit_tolerances(2));
    cleaningMethod = horzcat(do_postASR_ChanRej_TimeRej{1},'_ChanRej','0p',chanParam(strfind(chanParam,'.')+1:end),'_','TimeRej','0p',winParam(strfind(winParam,'.')+1:end),'_winTol',winTol);
    [~,EEG] = do_postASR_ChanRej_wrap(EEG,subj_char,cleaningMethod,do_postASR_ChanRej_TimeRej,mim_config,fullRankAvRefBool);
end

%% Save EEG and log
preprocess_pipeline = [cleaningMethod];
fprintf('Saving %s EEG to %s\n',subj_char,fullfile(mim_config.CleanEEG_output,preprocess_pipeline));
EEG.etc.CleanType  = preprocess_pipeline;
EEG.etc.Params     = mim_config;
mkdir(mim_config.CleanEEG_output)

[EEG_chans,EMG_chans,Noise_chans] = getChannelTypes_func(EEG) ;
% re-ref
fprintf('==== %s FINAL REREFERENCING ====\n',subj_char)
EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankAvRefBool);
%- save output
EEG = pop_select(EEG,'channel',sort([EEG_chans]));%2022-5-13 not use EMG
EEG = pop_saveset(EEG,'filepath',mim_config.CleanEEG_output,...
                    'filename',sprintf('%s_cleanEEG_%s.set',subj_char,preprocess_pipeline),...
                    'savemode','twofiles');
%- saving EMG output
% EEG = pop_saveset(EEG,'filepath',mim_config.CleanEEG_output,...
%                     'filename',sprintf('%s_cleanEEG_%s_wEMG_wNoise',subj_char,preprocess_pipeline),...
%                     'savemode','twofiles');
%% HPG IMPLEMENTATION ONLY
disp([num2str(length(EEG_chans)),' remained']);
float_fPath = [EEG.filepath filesep sprintf('%s_cleanEEG_%s.fdt',subj_char,preprocess_pipeline)];
[EEG,amica_cmd] = mim_prep_hpg_amica(EEG,float_fPath,mim_config.amicaEEG_output,mim_config.emailStr,avg_ref_pca_reduction);
EEG = pop_saveset(EEG,'filepath',mim_config.CleanEEG_output,...
                    'filename',sprintf('%s_cleanEEG_%s.set',subj_char,preprocess_pipeline),...
                    'savemode','twofiles');
toc
end

