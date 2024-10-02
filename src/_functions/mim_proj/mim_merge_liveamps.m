function [EEG] = mim_merge_liveamps(liveamp_fpaths,varargin)
%MIM_MERGE_LIVEAMPS Summary of this function goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Jacob Salminen
% Code Date: 01/10/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%## LOCAL PERMS
LIVEAMP_SUBDIR_NAMES = {'Cort_GY','Cort_WR','Noise_GY','Noise_WR'}; %Define folders to look for
LIVEAMP_SHORT_NAMES = {'CGY','CWR','NGY','NWR'}; %shorthand notation to use when renaming sets and events later
LIVEAMP_CHAN_STRUCT = struct('chan_names',{{{'Fp1';'AFp1';'AFz';'AF3';'AF7';'AFF5h';'AFF1h';'F1';'F3';'F5'; ...
                    'F7';'F9';'FFT9h';'FFT7h';'FFC5h';'FFC3h';'FFC1h';'FCz';'FC1';'FC3';...
                    'FC5';'FT7';'FT9';'FTT9h';'FTT7h';'FCC5h';'FCC3h';'FCC1h';'C1';'C3';...
                    'C5';'T7';'TTP7h';'CCP5h';'CCP3h';'CCP1h';'CP1';'CP3';'CP5';'TP7';'LISCM';'TPP9h';...
                    'TPP7h';'CPP5h';'CPP3h';'CPP1h';'P1';'P3';'P5';'P7';'PPO5h';'PPO1h';...
                    'POz';'PO3';'PO7';'PPO9h';'LSSCM';'LSTrap';'POO9h';'O1';'POO1';'LITrap';...
                    'OI1h';'Iz';'CGY-x';'CGY-y';'CGY-z'},...
                    {'Cz';'C2';'C4';'C6';'T8';'FTT10h';'FTT8h';'FCC6h';'FCC4h';'FCC2h';...
                    'FC2';'FC4';'FC6';'FT8';'FT10';'FFT10h';'FFT8h';'FFC6h';'FFC4h';'FFC2h';...
                    'Fz';'F2';'F4';'F6';'F8';'F10';'AFF6h';'AFF2h';'AF4';'AF8';...
                    'Fp2';'AFp2';'Oz';'OI2h';'RITrap';'RISCM';'POO10h';'O2';'POO2';'PO4';'PO8';'PPO10h';...
                    'RSSCM';'PPO6h';'PPO2h';'Pz';'P2';'P4';'P6';'P8';'TPP10h';'TPP8h';...
                    'CPP6h';'CPP4h';'CPP2h';'CP2';'CP4';'CP6';'TP8';'RSTrap';'TTP8h';'CCP6h';...
                    'CCP4h';'CCP2h';'CWR-x';'CWR-y';'CWR-z'},...
                    {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-FTT10h';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';...
                    'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';...
                    'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';...
                    'N-Fp2';'N-AFp2';'N-Oz';'N-OI2h';'N-RITrap';'N-RISCM';'N-POO10h';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';...
                    'N-RSSCM';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-TPP10h';'N-TPP8h';...
                    'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-RSTrap';'N-TTP8h';'N-CCP6h';...
                    'N-CCP4h';'N-CCP2h';'NWR-x';'NWR-y';'NWR-z'},...
                    {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-FTT10h';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';...
                    'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';...
                    'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';...
                    'N-Fp2';'N-AFp2';'N-Oz';'N-OI2h';'N-RITrap';'N-RISCM';'N-POO10h';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';...
                    'N-RSSCM';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-TPP10h';'N-TPP8h';...
                    'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-RSTrap';'N-TTP8h';'N-CCP6h';...
                    'N-CCP4h';'N-CCP2h';'NWR-x';'NWR-y';'NWR-z'}}},...
    'liveamp_name',{{'cort_gy',...
                    'cort_wr',...
                    'noise_gy',...
                    'noise_wr'}},...
    'ref_chan',{{'CPz',...
                'CPz',...
                'N-CPz',...
                'N-CPz'}});
%##
cat_logo();
%## TIME
t = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'liveamp_fpaths',@ischar);
%## PARSE
parse(p, liveamp_fpaths, varargin{:});
%## SET DEFAULTS
%% ===================================================================== %%
%## DEFINE LOCAL PARAMS
EEG = [];
%## MAIN FUNCTION MEAT

%% Actually Gather and Merge Data
inputFolder = fullfile(MIMDataFolder,subjStr,'EEG','Raw');
%% IMPORT LIVEAMP FILES
%## State your business and clock in
disp('Gathering and importing all LiveAmpFiles');
tic
%## Select subject directory
chk = any(cellfun(@(x) ~exist(x,'dir'),liveamp_fpaths));
if chk
    fprintf('One of the LiveAmp folders doesn''t exist...\n')
    liveamp_fpaths = cell(length(LIVEAMP_SUBDIR_NAMES),1);
    tmpdir = uigetdir(startingDir, 'Pick folder which contains subfolders of all separate LiveAmp64 data (Cort_GY, Noise_WR, etc)');
    for i = 1:length(LIVEAMP_SUBDIR_NAMES)
        liveamp_fpaths{i} = [tmpdir filesep LIVEAMP_SUBDIR_NAMES{i}];
    end
end

%## process all subfolders
chk = any(cellfun(@(x) ~exist(x,'dir'),liveamp_fpaths));
% tmp_fname_noext = [];
for la_i = 1:length(liveamp_fpaths) %for each expected file
    subfolderName = liveamp_fpaths{la_i};
    
    if chk(la_i) %make sure subdirectory exists
        tmp_fdir = dir([liveamp_fpaths filehsep '*.vhdr']);
        switch length(tmp_fdir)
            case 0 
                error('No Brain Vision files (.vhdr) in %s...\n',liveamp_fpaths{la_i});
            case 1
                tmp_fname = tmp_fdir.name;
                [~,tmp_fname_noext,~] = fileparts(tmp_fname);
            otherwise
                error('Too many Brain Vision files (.vhdr) in %s...\n',liveamp_fpaths{la_i});
        end
        
        %## Load file
        liveamp_fpaths{la_i} = liveamp_fpaths{la_i}; %fullfile(subjDir,subfolderName);
        EEG = pop_loadbv(liveamp_fpaths{la_i}, tmp_fname); %load data using special plugin for BV files

        %## Get start time of recording
        vmrkFile = fullfile(liveamp_fpaths{la_i},[tmp_fname_noext,'.vmrk']);
        [EEG.recordingStart] = getBVRecordingStartTime(vmrkFile);
        disp(EEG.recordingStart)
        
        %## select correct channel info
        eegCapV2_criticalDateTime = datetime(2021,5,4,0,0,0); %made new EEG cap (V2) and used it first for MIM on May 4th (Note Amanda used it earlier for her own study, not MIM)
        if EEG.recordingStart.datetime >= eegCapV2_criticalDateTime
            %use v2 info *CHRISTINA*
            disp('I see you are using version 2 of our custom EEG cap (May 4th 2021 and beyond for MIM)');
            chSelection = lower(subfolderName);
            ind = strcmp(LIVEAMP_CHAN_STRUCT.liveamp_name,'cort_gy');
            chNames = LIVEAMP_CHAN_STRUCT(ind).chan_names;
            refCh = LIVEAMP_CHAN_STRUCT(ind).ref_chan;

            switch lower(chSelection)
                
                case 'cort_gy'
                    ind = strcmp(LIVEAMP_CHAN_STRUCT.liveamp_name,'cort_gy');
                    chNames = {'Fp1';'AFp1';'AFz';'AF3';'AF7';'AFF5h';'AFF1h';'F1';'F3';'F5'; ...
                        'F7';'F9';'FFT9h';'FFT7h';'FFC5h';'FFC3h';'FFC1h';'FCz';'FC1';'FC3';...
                        'FC5';'FT7';'FT9';'FTT9h';'FTT7h';'FCC5h';'FCC3h';'FCC1h';'C1';'C3';...
                        'C5';'T7';...
                        'TTP7h';'CCP5h';'CCP3h';'CCP1h';'CP1';'CP3';'CP5';'TP7';'LISCM';'TPP9h';...
                        'TPP7h';'CPP5h';'CPP3h';'CPP1h';'P1';'P3';'P5';'P7';'PPO5h';'PPO1h';...
                        'POz';'PO3';'PO7';'PPO9h';'LSSCM';'LSTrap';'POO9h';'O1';'POO1';'LITrap';...
                        'OI1h';'Iz';...
                        'CGY-x';'CGY-y';'CGY-z'};
                    refCh = 'CPz';
                case 'cort_wr'
                    chNames = {'Cz';'C2';'C4';'C6';'T8';'FTT10h';'FTT8h';'FCC6h';'FCC4h';'FCC2h';...
                        'FC2';'FC4';'FC6';'FT8';'FT10';'FFT10h';'FFT8h';'FFC6h';'FFC4h';'FFC2h';...
                        'Fz';'F2';'F4';'F6';'F8';'F10';'AFF6h';'AFF2h';'AF4';'AF8';...
                        'Fp2';'AFp2';...
                        'Oz';'OI2h';'RITrap';'RISCM';'POO10h';'O2';'POO2';'PO4';'PO8';'PPO10h';...
                        'RSSCM';'PPO6h';'PPO2h';'Pz';'P2';'P4';'P6';'P8';'TPP10h';'TPP8h';...
                        'CPP6h';'CPP4h';'CPP2h';'CP2';'CP4';'CP6';'TP8';'RSTrap';'TTP8h';'CCP6h';...
                        'CCP4h';'CCP2h';...
                        'CWR-x';'CWR-y';'CWR-z'};
                    refCh = 'CPz';
                case 'noise_gy'
                    chNames = {'N-Fp1';'N-AFp1';'N-AFz';'N-AF3';'N-AF7';'N-AFF5h';'N-AFF1h';'N-F1';'N-F3';'N-F5'; ...
                        'N-F7';'N-F9';'N-FFT9h';'N-FFT7h';'N-FFC5h';'N-FFC3h';'N-FFC1h';'N-FCz';'N-FC1';'N-FC3';...
                        'N-FC5';'N-FT7';'N-FT9';'N-FTT9h';'N-FTT7h';'N-FCC5h';'N-FCC3h';'N-FCC1h';'N-C1';'N-C3';...
                        'N-C5';'N-T7';...
                        'N-TTP7h';'N-CCP5h';'N-CCP3h';'N-CCP1h';'N-CP1';'N-CP3';'N-CP5';'N-TP7';'N-LISCM';'N-TPP9h';...
                        'N-TPP7h';'N-CPP5h';'N-CPP3h';'N-CPP1h';'N-P1';'N-P3';'N-P5';'N-P7';'N-PPO5h';'N-PPO1h';...
                        'N-POz';'N-PO3';'N-PO7';'N-PPO9h';'N-LSSCM';'N-LSTrap';'N-POO9h';'N-O1';'N-POO1';'N-LITrap';...
                        'N-OI1h';'N-Iz';...
                        'NGY-x';'NGY-y';'NGY-z'};
                    refCh = 'N-CPz';
                case 'noise_wr'
                    chNames = {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-FTT10h';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';...
                        'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';...
                        'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';...
                        'N-Fp2';'N-AFp2';...
                        'N-Oz';'N-OI2h';'N-RITrap';'N-RISCM';'N-POO10h';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';...
                        'N-RSSCM';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-TPP10h';'N-TPP8h';...
                        'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-RSTrap';'N-TTP8h';'N-CCP6h';...
                        'N-CCP4h';'N-CCP2h';...
                        'NWR-x';'NWR-y';'NWR-z'};
                    refCh = 'N-CPz';
            end
            
        else %use v1 cap info
            disp('It looks like you are using the original version (v1) of our custom EEG cap (any recording prior to May 4th 2021 for the MIM study)');
            chSelection = subfolderName;
            switch lower(chSelection)
                case 'cort_gy'
                    chNames = {'LSSCM';'AFp1';'AFz';'AF3';'AF7';'AFF5h';'AFF1h';'F1';'F3';'F5';'F7';'F9';'FFT9h';'FFT7h';'FFC5h';'FFC3h';'FFC1h';'FCz';'FC1';'FC3';'FC5';'FT7';'FT9';'LISCM';'FTT7h';'FCC5h';'FCC3h';'FCC1h';'C1';'C3';'C5';'T7';'TTP7h';'CCP5h';'CCP3h';'CCP1h';'CP1';'CP3';'CP5';'TP7';'TP9';'LSTrap';'TPP7h';'CPP5h';'CPP3h';'CPP1h';'P1';'P3';'P5';'P7';'PPO5h';'PPO1h';'POz';'PO3';'PO7';'PPO9h';'P9';'PO9';'LITrap';'O1';'POO1';'O9';'OI1h';'Iz';'CGY-x';'CGY-y';'CGY-z'};
                    refCh = 'CPz';
                case 'cort_wr'
                    chNames = {'Cz';'C2';'C4';'C6';'T8';'RISCM';'FTT8h';'FCC6h';'FCC4h';'FCC2h';'FC2';'FC4';'FC6';'FT8';'FT10';'FFT10h';'FFT8h';'FFC6h';'FFC4h';'FFC2h';'Fz';'F2';'F4';'F6';'F8';'F10';'AFF6h';'AFF2h';'AF4';'AF8';'RSSCM';'AFp2';'Oz';'OI2h';'O10';'PO10';'RITrap';'O2';'POO2';'PO4';'PO8';'PPO10h';'P10';'PPO6h';'PPO2h';'Pz';'P2';'P4';'P6';'P8';'RSTrap';'TPP8h';'CPP6h';'CPP4h';'CPP2h';'CP2';'CP4';'CP6';'TP8';'TP10';'TTP8h';'CCP6h';'CCP4h';'CCP2h';'CWR-x';'CWR-y';'CWR-z'};
                    refCh = 'CPz';
                case 'noise_gy'
                    chNames = {'N-LSEOG';'N-AFp1';'N-AFz';'N-AF3';'N-AF7';'N-AFF5h';'N-AFF1h';'N-F1';'N-F3';'N-F5';'N-F7';'N-F9';'N-FFT9h';'N-FFT7h';'N-FFC5h';'N-FFC3h';'N-FFC1h';'N-FCz';'N-FC1';'N-FC3';'N-FC5';'N-FT7';'N-FT9';'N-LIEOG';'N-FTT7h';'N-FCC5h';'N-FCC3h';'N-FCC1h';'N-C1';'N-C3';'N-C5';'N-T7';'N-TTP7h';'N-CCP5h';'N-CCP3h';'N-CCP1h';'N-CP1';'N-CP3';'N-CP5';'N-TP7';'N-TP9';'N-LSJaw';'N-TPP7h';'N-CPP5h';'N-CPP3h';'N-CPP1h';'N-P1';'N-P3';'N-P5';'N-P7';'N-PPO5h';'N-PPO1h';'N-POz';'N-PO3';'N-PO7';'N-PPO9h';'N-P9';'N-PO9';'N-LIJaw';'N-O1';'N-POO1';'N-O9';'N-OI1h';'N-Iz';'NGY-x';'NGY-y';'NGY-z'};
                    refCh = 'N-CPz';
                case 'noise_wr'
                    chNames = {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-RIEOG';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';'N-RSEOG';'N-AFp2';'N-Oz';'N-OI2h';'N-O10';'N-PO10';'N-RIJaw';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';'N-P10';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-RSJaw';'N-TPP8h';'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-TP10';'N-TTP8h';'N-CCP6h';'N-CCP4h';'N-CCP2h';'NWR-x';'NWR-y';'NWR-z'};
                    refCh = 'N-CPz';
            end
        end

        


        %## update set name
%         EEG.setname = subfolderName; %e.g. cort_GY
        EEG.setname = shortHandNames{la_i}; %e.g. CGY
        
        %## update channel info
        chType = {};
        for ch_i = 1:EEG.nbchan %Note: EEG.nbchan accomodates 64 and 67 channel case (whether accelerometers recorded or not) whereas had we used length(chNames) then it might throw an error when there are fewer channels than expected
            tempLabel = chNames{ch_i};
            EEG.chanlocs(ch_i).labels = chNames{ch_i};
            
            if contains(tempLabel,{'N-LSEOG', 'N-LIEOG', 'N-RSEOG', 'N-RIEOG'},'IgnoreCase',true) 
                chType{ch_i} = 'N-EOG'; %these are junk (idea to repurose noise sensors for extra EOG or EMG did not work out)
            elseif contains(tempLabel,{'N-LSJaw','N-LIJaw','N-RSJaw','N-RIJaw', ...
                                    'N-LSSCM', 'N-LISCM', 'N-LSTrap', 'N-LITrap', ...
                                    'N-RSSCM', 'N-RISCM', 'N-RSTrap', 'N-RITrap'},'IgnoreCase',true)
                chType{ch_i} = 'N-EMG';%these are junk. 
                %Note N-EOG and N-Jaw muscles we tried to record with V1 of
                %cap but it did not work out. Meanwhile, N-NeckMuscles
                %(SCM,Trap) are just electrodes connected to nothing at
                %all. That's just how we chose to name them for v2 version
                %of the cap to make them easy to delete later (N-LSSCM is
                %the noise electrode that would otherwise be paired to
                %whatever cortical electrode was repurposed to be a neck
                %EMG sensor but N-LSSCM is not recording anything of
                %value)
            elseif contains(tempLabel,'N-','IgnoreCase',true)
                chType{ch_i} = 'Noise';
            elseif contains(tempLabel,{'-x','-y','-z'},'IgnoreCase',true)
                chType{ch_i} = 'Acc';
            elseif contains(tempLabel,{'Trap','SCM'},'IgnoreCase',true)
                chType{ch_i} = 'EMG';
            else
                chType{ch_i} = 'EEG';
            end
        end
        
        EEG_chans = find(strcmpi(chType,'EEG')); %normal electrodes
        Noise_chans = find(strcmpi(chType,'Noise')); %noise electrodes
        NEOG_chans = find(strcmpi(chType,'N-EOG')); %EOG electrodes measured via noise system (sometimes works, sometimes doesn't)
        EMG_chans = find(strcmpi(chType,'EMG')); %neck EMG
        NEMG_chans = find(strcmpi(chType,'N-EMG')); %EMG electrodes measured via noise system (sometimes works, sometimes doesn't)
        Acc_chans = find(strcmpi(chType,'Acc')); %accelerometers
%         Bad_chans = find(strcmpi(chType,'Bad')); %channels we know to never both with in the first place
% %         EEG = pop_select(EEG,'nochannel',Bad_chans);
        
        EEG=pop_chanedit(EEG, 'settype',{num2str(EEG_chans), 'EEG'});
        EEG=pop_chanedit(EEG, 'settype',{num2str(Noise_chans), 'Noise'});
        EEG=pop_chanedit(EEG, 'settype',{num2str(NEOG_chans), 'N-EOG'});
        EEG=pop_chanedit(EEG, 'settype',{num2str(EMG_chans), 'EMG'});
        EEG=pop_chanedit(EEG, 'settype',{num2str(NEMG_chans), 'N-EMG'});
        EEG=pop_chanedit(EEG, 'settype',{num2str(Acc_chans), 'Acc'});
        
        
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, la_i);
        
        %add generic location info
        [EEGLABDirectory,~,~] = fileparts(which('eeglab'));
        %     fileLoc_BEM = fullfile(EEGLABDirectory,'plugins\dipfit2.3\standard_BEM\elec\standard_1005.elc');
%                 fileLoc_BEM = fullfile(EEGLABDirectory,'plugins\dipfit2.3\standard_BEM\elec\standard_1005_FixedO9andO10.elc'); %get modified file from Ryan
        
        %CHRISTINA BROOKS LOOK HERE 2021-10-13
        dipfitDir = dir(fullfile(EEGLABDirectory,'plugins','dipfit*')); %2021-10-13 not currently robust to there being .zip files alonside active dipfit directory
        dipfitDir = dipfitDir.name;
        
        try
                fileLoc_BEM = fullfile(EEGLABDirectory,'plugins',dipfitDir,'standard_BEM\elec\standard_1005_FixedO9andO10.elc'); %get modified file from Ryan
                EEG=pop_chanedit(EEG, 'nosedir','+Y', 'lookup', fileLoc_BEM);
                disp('Successfully loaded custom-made file for default locations that Ryan made: standard_1005_FixedO9andO10.elc');
                [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        catch
                disp('Could not find the custom .elc file Ryan made that contains default locations for channels O9 and O10. I am resorting to the premade version that comes with EEGLAB which lacks O9 and O10');
                fileLoc_BEM = fullfile(EEGLABDirectory,'plugins',dipfitDir,'standard_BEM\elec\standard_1005.elc'); %original (not as good) file. You will probably be missing location info for a couple of channels but at least you have something
                EEG=pop_chanedit(EEG, 'nosedir','+Y', 'lookup', fileLoc_BEM);
                [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        end

        %set reference
        EEG=pop_chanedit(EEG, 'setref',{num2str([EEG_chans EMG_chans]) refCh}); %only EEG and neck EMG chans have CPz reference
        EEG=pop_chanedit(EEG, 'setref',{num2str([Noise_chans NEOG_chans NEMG_chans]) refCh}); %only EEG and neck EMG chans have CPz reference
        
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
       
        
    else %tell user if subdirectory doesn't exist
        warning(['Cannot find data for ', subfolderName])
    end
    
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, la_i); 
end

%## Tell them you finished and that you hope they didn't have to wait too long
disp('Done importing all individual live amp files.');
%%
%more memory efficient
TMP = gatherAndImportLiveAmpFiles(liveamp_fpaths{la_i});
EEG = AlignLiveAmps_new(,DiaryOutputFolder);

%Add global times to events to make it easier to look at manually
EEG = addGlobalTimesToEvents(EEG);

%% Save data
EEG.setname = [subjStr,'_EEG']; 
end

%% ===================================================================== %%
%##
function [ALLEEG] = gatherAndImportLiveAmpFiles(inputFolder)
    %GATHERANDIMPORTLIVEAMPFILES Summary of this function goes here
    %   Detailed explanation goes here

    %## State your business and clock in
    disp('Gathering and importing all LiveAmpFiles');
    tic
    %## Select subject directory
    
    subjDir = inputFolder; 
    if isempty(subjDir)
        startingDir = cd;
        subjDir = uigetdir(startingDir, 'Pick folder which contains subfolders of all separate LiveAmp64 data (Cort_GY, Noise_WR, etc)');
    end
    
    %## process all subfolders
    individualDirs = {'Cort_GY', 'Cort_WR', 'Noise_GY', 'Noise_WR'}; %Define folders to look for
    LIVEAMP_SHORT_NAMES = {'CGY',    'CWR',      'NGY',      'NWR'}; %shorthand notation to use when renaming sets and events later
    fileNameNoExt = [];
    for LiveAmp_i = 1:length(individualDirs) %for each expected file
        subfolderName = individualDirs{LiveAmp_i};
        
        if isdir(fullfile(subjDir,subfolderName)) %make sure subdirectory exists
            disp(['Found ', subfolderName,' subfolder!']);
            tempFileSearch = dir(fullfile(subjDir,subfolderName,'*.vhdr'));
            switch length(tempFileSearch)
                case 0 
                    error('Could not find a Brain Vision (.vhdr) file in that subdirectory!!');
                    return;
                case 1
                    fileName = tempFileSearch.name;
                    [~,fileNameNoExt,EXT] = fileparts(fileName);
                otherwise
                    error('Too many Brain Vision files (.vhdr) found in that subdirectory!');
                    return;
            end
            
    
            %## Load file
            inputFolder = fullfile(subjDir,subfolderName);
    %         EEG = pop_loadbv_nodata(inputFolder, fileName); %load data using special plugin for BV files
            EEG = pop_loadbv(inputFolder, fileName); %load data using special plugin for BV files
            %## Get start time of recording
            vmrkFile = fullfile(inputFolder,[fileNameNoExt,'.vmrk']);
            [ EEG.recordingStart ] = getBVRecordingStartTime( vmrkFile );
            disp(EEG.recordingStart)
            
            %## select correct channel info
            eegCapV2_criticalDateTime = datetime(2021,5,4,0,0,0); %made new EEG cap (V2) and used it first for MIM on May 4th (Note Amanda used it earlier for her own study, not MIM)
            if EEG.recordingStart.datetime >= eegCapV2_criticalDateTime
                %use v2 info *CHRISTINA*
                disp('I see you are using version 2 of our custom EEG cap (May 4th 2021 and beyond for MIM)');
                chSelection = subfolderName;
                switch lower(chSelection)
                    case 'cort_gy'
                        chNames = {'Fp1';'AFp1';'AFz';'AF3';'AF7';'AFF5h';'AFF1h';'F1';'F3';'F5'; ...
                            'F7';'F9';'FFT9h';'FFT7h';'FFC5h';'FFC3h';'FFC1h';'FCz';'FC1';'FC3';...
                            'FC5';'FT7';'FT9';'FTT9h';'FTT7h';'FCC5h';'FCC3h';'FCC1h';'C1';'C3';...
                            'C5';'T7';...
                            'TTP7h';'CCP5h';'CCP3h';'CCP1h';'CP1';'CP3';'CP5';'TP7';'LISCM';'TPP9h';...
                            'TPP7h';'CPP5h';'CPP3h';'CPP1h';'P1';'P3';'P5';'P7';'PPO5h';'PPO1h';...
                            'POz';'PO3';'PO7';'PPO9h';'LSSCM';'LSTrap';'POO9h';'O1';'POO1';'LITrap';...
                            'OI1h';'Iz';...
                            'CGY-x';'CGY-y';'CGY-z'};
                        refCh = 'CPz';
                    case 'cort_wr'
                        chNames = {'Cz';'C2';'C4';'C6';'T8';'FTT10h';'FTT8h';'FCC6h';'FCC4h';'FCC2h';...
                            'FC2';'FC4';'FC6';'FT8';'FT10';'FFT10h';'FFT8h';'FFC6h';'FFC4h';'FFC2h';...
                            'Fz';'F2';'F4';'F6';'F8';'F10';'AFF6h';'AFF2h';'AF4';'AF8';...
                            'Fp2';'AFp2';...
                            'Oz';'OI2h';'RITrap';'RISCM';'POO10h';'O2';'POO2';'PO4';'PO8';'PPO10h';...
                            'RSSCM';'PPO6h';'PPO2h';'Pz';'P2';'P4';'P6';'P8';'TPP10h';'TPP8h';...
                            'CPP6h';'CPP4h';'CPP2h';'CP2';'CP4';'CP6';'TP8';'RSTrap';'TTP8h';'CCP6h';...
                            'CCP4h';'CCP2h';...
                            'CWR-x';'CWR-y';'CWR-z'};
                        refCh = 'CPz';
                    case 'noise_gy'
                        chNames = {'N-Fp1';'N-AFp1';'N-AFz';'N-AF3';'N-AF7';'N-AFF5h';'N-AFF1h';'N-F1';'N-F3';'N-F5'; ...
                            'N-F7';'N-F9';'N-FFT9h';'N-FFT7h';'N-FFC5h';'N-FFC3h';'N-FFC1h';'N-FCz';'N-FC1';'N-FC3';...
                            'N-FC5';'N-FT7';'N-FT9';'N-FTT9h';'N-FTT7h';'N-FCC5h';'N-FCC3h';'N-FCC1h';'N-C1';'N-C3';...
                            'N-C5';'N-T7';...
                            'N-TTP7h';'N-CCP5h';'N-CCP3h';'N-CCP1h';'N-CP1';'N-CP3';'N-CP5';'N-TP7';'N-LISCM';'N-TPP9h';...
                            'N-TPP7h';'N-CPP5h';'N-CPP3h';'N-CPP1h';'N-P1';'N-P3';'N-P5';'N-P7';'N-PPO5h';'N-PPO1h';...
                            'N-POz';'N-PO3';'N-PO7';'N-PPO9h';'N-LSSCM';'N-LSTrap';'N-POO9h';'N-O1';'N-POO1';'N-LITrap';...
                            'N-OI1h';'N-Iz';...
                            'NGY-x';'NGY-y';'NGY-z'};
                        refCh = 'N-CPz';
                    case 'noise_wr'
                        chNames = {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-FTT10h';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';...
                            'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';...
                            'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';...
                            'N-Fp2';'N-AFp2';...
                            'N-Oz';'N-OI2h';'N-RITrap';'N-RISCM';'N-POO10h';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';...
                            'N-RSSCM';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-TPP10h';'N-TPP8h';...
                            'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-RSTrap';'N-TTP8h';'N-CCP6h';...
                            'N-CCP4h';'N-CCP2h';...
                            'NWR-x';'NWR-y';'NWR-z'};
                        refCh = 'N-CPz';
                end
                
            else %use v1 cap info
                disp('It looks like you are using the original version (v1) of our custom EEG cap (any recording prior to May 4th 2021 for the MIM study)');
                chSelection = subfolderName;
                switch lower(chSelection)
                    case 'cort_gy'
                        chNames = {'LSSCM';'AFp1';'AFz';'AF3';'AF7';'AFF5h';'AFF1h';'F1';'F3';'F5';'F7';'F9';'FFT9h';'FFT7h';'FFC5h';'FFC3h';'FFC1h';'FCz';'FC1';'FC3';'FC5';'FT7';'FT9';'LISCM';'FTT7h';'FCC5h';'FCC3h';'FCC1h';'C1';'C3';'C5';'T7';'TTP7h';'CCP5h';'CCP3h';'CCP1h';'CP1';'CP3';'CP5';'TP7';'TP9';'LSTrap';'TPP7h';'CPP5h';'CPP3h';'CPP1h';'P1';'P3';'P5';'P7';'PPO5h';'PPO1h';'POz';'PO3';'PO7';'PPO9h';'P9';'PO9';'LITrap';'O1';'POO1';'O9';'OI1h';'Iz';'CGY-x';'CGY-y';'CGY-z'};
                        refCh = 'CPz';
                    case 'cort_wr'
                        chNames = {'Cz';'C2';'C4';'C6';'T8';'RISCM';'FTT8h';'FCC6h';'FCC4h';'FCC2h';'FC2';'FC4';'FC6';'FT8';'FT10';'FFT10h';'FFT8h';'FFC6h';'FFC4h';'FFC2h';'Fz';'F2';'F4';'F6';'F8';'F10';'AFF6h';'AFF2h';'AF4';'AF8';'RSSCM';'AFp2';'Oz';'OI2h';'O10';'PO10';'RITrap';'O2';'POO2';'PO4';'PO8';'PPO10h';'P10';'PPO6h';'PPO2h';'Pz';'P2';'P4';'P6';'P8';'RSTrap';'TPP8h';'CPP6h';'CPP4h';'CPP2h';'CP2';'CP4';'CP6';'TP8';'TP10';'TTP8h';'CCP6h';'CCP4h';'CCP2h';'CWR-x';'CWR-y';'CWR-z'};
                        refCh = 'CPz';
                    case 'noise_gy'
                        chNames = {'N-LSEOG';'N-AFp1';'N-AFz';'N-AF3';'N-AF7';'N-AFF5h';'N-AFF1h';'N-F1';'N-F3';'N-F5';'N-F7';'N-F9';'N-FFT9h';'N-FFT7h';'N-FFC5h';'N-FFC3h';'N-FFC1h';'N-FCz';'N-FC1';'N-FC3';'N-FC5';'N-FT7';'N-FT9';'N-LIEOG';'N-FTT7h';'N-FCC5h';'N-FCC3h';'N-FCC1h';'N-C1';'N-C3';'N-C5';'N-T7';'N-TTP7h';'N-CCP5h';'N-CCP3h';'N-CCP1h';'N-CP1';'N-CP3';'N-CP5';'N-TP7';'N-TP9';'N-LSJaw';'N-TPP7h';'N-CPP5h';'N-CPP3h';'N-CPP1h';'N-P1';'N-P3';'N-P5';'N-P7';'N-PPO5h';'N-PPO1h';'N-POz';'N-PO3';'N-PO7';'N-PPO9h';'N-P9';'N-PO9';'N-LIJaw';'N-O1';'N-POO1';'N-O9';'N-OI1h';'N-Iz';'NGY-x';'NGY-y';'NGY-z'};
                        refCh = 'N-CPz';
                    case 'noise_wr'
                        chNames = {'N-Cz';'N-C2';'N-C4';'N-C6';'N-T8';'N-RIEOG';'N-FTT8h';'N-FCC6h';'N-FCC4h';'N-FCC2h';'N-FC2';'N-FC4';'N-FC6';'N-FT8';'N-FT10';'N-FFT10h';'N-FFT8h';'N-FFC6h';'N-FFC4h';'N-FFC2h';'N-Fz';'N-F2';'N-F4';'N-F6';'N-F8';'N-F10';'N-AFF6h';'N-AFF2h';'N-AF4';'N-AF8';'N-RSEOG';'N-AFp2';'N-Oz';'N-OI2h';'N-O10';'N-PO10';'N-RIJaw';'N-O2';'N-POO2';'N-PO4';'N-PO8';'N-PPO10h';'N-P10';'N-PPO6h';'N-PPO2h';'N-Pz';'N-P2';'N-P4';'N-P6';'N-P8';'N-RSJaw';'N-TPP8h';'N-CPP6h';'N-CPP4h';'N-CPP2h';'N-CP2';'N-CP4';'N-CP6';'N-TP8';'N-TP10';'N-TTP8h';'N-CCP6h';'N-CCP4h';'N-CCP2h';'NWR-x';'NWR-y';'NWR-z'};
                        refCh = 'N-CPz';
                end
                
                
            end
    
            
    
    
            %## update set name
    %         EEG.setname = subfolderName; %e.g. cort_GY
            EEG.setname = LIVEAMP_SHORT_NAMES{LiveAmp_i}; %e.g. CGY
            
            %## update channel info
            chType = {};
            for ch_i = 1:EEG.nbchan %Note: EEG.nbchan accomodates 64 and 67 channel case (whether accelerometers recorded or not) whereas had we used length(chNames) then it might throw an error when there are fewer channels than expected
                tempLabel = chNames{ch_i};
                EEG.chanlocs(ch_i).labels = chNames{ch_i};
                
                if contains(tempLabel,{'N-LSEOG', 'N-LIEOG', 'N-RSEOG', 'N-RIEOG'},'IgnoreCase',true) 
                    chType{ch_i} = 'N-EOG'; %these are junk (idea to repurose noise sensors for extra EOG or EMG did not work out)
                elseif contains(tempLabel,{'N-LSJaw','N-LIJaw','N-RSJaw','N-RIJaw', ...
                                        'N-LSSCM', 'N-LISCM', 'N-LSTrap', 'N-LITrap', ...
                                        'N-RSSCM', 'N-RISCM', 'N-RSTrap', 'N-RITrap'},'IgnoreCase',true)
                    chType{ch_i} = 'N-EMG';%these are junk. 
                    %Note N-EOG and N-Jaw muscles we tried to record with V1 of
                    %cap but it did not work out. Meanwhile, N-NeckMuscles
                    %(SCM,Trap) are just electrodes connected to nothing at
                    %all. That's just how we chose to name them for v2 version
                    %of the cap to make them easy to delete later (N-LSSCM is
                    %the noise electrode that would otherwise be paired to
                    %whatever cortical electrode was repurposed to be a neck
                    %EMG sensor but N-LSSCM is not recording anything of
                    %value)
                elseif contains(tempLabel,'N-','IgnoreCase',true)
                    chType{ch_i} = 'Noise';
                elseif contains(tempLabel,{'-x','-y','-z'},'IgnoreCase',true)
                    chType{ch_i} = 'Acc';
                elseif contains(tempLabel,{'Trap','SCM'},'IgnoreCase',true)
                    chType{ch_i} = 'EMG';
                else
                    chType{ch_i} = 'EEG';
                end
            end
            
            EEG_chans = find(strcmpi(chType,'EEG')); %normal electrodes
            Noise_chans = find(strcmpi(chType,'Noise')); %noise electrodes
            NEOG_chans = find(strcmpi(chType,'N-EOG')); %EOG electrodes measured via noise system (sometimes works, sometimes doesn't)
            EMG_chans = find(strcmpi(chType,'EMG')); %neck EMG
            NEMG_chans = find(strcmpi(chType,'N-EMG')); %EMG electrodes measured via noise system (sometimes works, sometimes doesn't)
            Acc_chans = find(strcmpi(chType,'Acc')); %accelerometers
    %         Bad_chans = find(strcmpi(chType,'Bad')); %channels we know to never both with in the first place
    % %         EEG = pop_select(EEG,'nochannel',Bad_chans);
            
            EEG=pop_chanedit(EEG, 'settype',{num2str(EEG_chans), 'EEG'});
            EEG=pop_chanedit(EEG, 'settype',{num2str(Noise_chans), 'Noise'});
            EEG=pop_chanedit(EEG, 'settype',{num2str(NEOG_chans), 'N-EOG'});
            EEG=pop_chanedit(EEG, 'settype',{num2str(EMG_chans), 'EMG'});
            EEG=pop_chanedit(EEG, 'settype',{num2str(NEMG_chans), 'N-EMG'});
            EEG=pop_chanedit(EEG, 'settype',{num2str(Acc_chans), 'Acc'});
            
            
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, LiveAmp_i);
            
            %add generic location info
            [EEGLABDirectory,~,~] = fileparts(which('eeglab'));
            %     fileLoc_BEM = fullfile(EEGLABDirectory,'plugins\dipfit2.3\standard_BEM\elec\standard_1005.elc');
    %                 fileLoc_BEM = fullfile(EEGLABDirectory,'plugins\dipfit2.3\standard_BEM\elec\standard_1005_FixedO9andO10.elc'); %get modified file from Ryan
            
            %CHRISTINA BROOKS LOOK HERE 2021-10-13
            dipfitDir = dir(fullfile(EEGLABDirectory,'plugins','dipfit*')); %2021-10-13 not currently robust to there being .zip files alonside active dipfit directory
            dipfitDir = dipfitDir.name;
            
            try
                    fileLoc_BEM = fullfile(EEGLABDirectory,'plugins',dipfitDir,'standard_BEM\elec\standard_1005_FixedO9andO10.elc'); %get modified file from Ryan
                    EEG=pop_chanedit(EEG, 'nosedir','+Y', 'lookup', fileLoc_BEM);
                    disp('Successfully loaded custom-made file for default locations that Ryan made: standard_1005_FixedO9andO10.elc');
                    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            catch
                    disp('Could not find the custom .elc file Ryan made that contains default locations for channels O9 and O10. I am resorting to the premade version that comes with EEGLAB which lacks O9 and O10');
                    fileLoc_BEM = fullfile(EEGLABDirectory,'plugins',dipfitDir,'standard_BEM\elec\standard_1005.elc'); %original (not as good) file. You will probably be missing location info for a couple of channels but at least you have something
                    EEG=pop_chanedit(EEG, 'nosedir','+Y', 'lookup', fileLoc_BEM);
                    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end
    
            %set reference
            EEG=pop_chanedit(EEG, 'setref',{num2str([EEG_chans EMG_chans]) refCh}); %only EEG and neck EMG chans have CPz reference
            EEG=pop_chanedit(EEG, 'setref',{num2str([Noise_chans NEOG_chans NEMG_chans]) refCh}); %only EEG and neck EMG chans have CPz reference
            
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
           
            
        else %tell user if subdirectory doesn't exist
            warning(['Cannot find data for ', subfolderName])
        end
        
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, LiveAmp_i); 
    end
    
    %## Tell them you finished and that you hope they didn't have to wait too long
    disp('Done importing all individual live amp files.');
    % eeglab redraw;
    toc
end

%% ===================================================================== %%
function [EEG] = addGlobalTimesToEvents(EEG)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    %BEWARE OF CROPPING! FURTHER TESTING NEEDED
    for event_i = 1:length(EEG.event)
    
        EEG.event(event_i).datetime = EEG.recordingStart.datetime + seconds( (EEG.event(event_i).latency - 1)/EEG.srate);
        EEG.event(event_i).datestr = datestr(EEG.event(event_i).datetime);
    
    end
end

%% ===================================================================== %%
function [EEG_out] = AlignLiveAmps_new(ALLEEG,savDir)
    %ALIGNLIVEAMPS Summary of this function goes here
    %   Detailed explanation goes here
    
    if isempty(ALLEEG)
        disp('add functionality later');
    end
    
    nAmps = length(ALLEEG);
    
    %## State your business and clock in
    disp('Merging all LiveAmp Files');
    tic
    
    %## make sure all EEG structures have the same fields for their EVENT structure
    % find all possible fields
    f = fieldnames(ALLEEG(1).event);
    for amp_i = 2:nAmps
        f2 = fieldnames(ALLEEG(amp_i).event);
        f = union(f,f2);
    end
    %go back and make sure all eeg sets have those fields
    for amp_i = 1:nAmps
        for field_i = 1:length(f)
            FIELD = f{field_i};
            X = [];
            if ~isfield(ALLEEG(amp_i).event,FIELD)
                %                 tempEEG(:).f(field_i) = [];
                [ALLEEG(amp_i).event.(FIELD)] = deal(X);
            end
        end
    end
    
    %## make sure all EEG structures have the same fields for their UREVENT structure
    % find all possible fields
    f = fieldnames(ALLEEG(1).urevent);
    for amp_i = 2:nAmps
        f2 = fieldnames(ALLEEG(amp_i).urevent);
        f = union(f,f2);
    end
    %go back and make sure all eeg sets have those fields
    for amp_i = 1:nAmps
        for field_i = 1:length(f)
            FIELD = f{field_i};
            X = [];
            if ~isfield(ALLEEG(amp_i).urevent,FIELD)
                %                 tempEEG(:).f(field_i) = [];
                [ALLEEG(amp_i).urevent.(FIELD)] = deal(X);
            end
        end
    end
    
    %## make sure all EEG structures have the same fields for their CHANLOCS structure
    % find all possible fields
    f = fieldnames(ALLEEG(1).chanlocs);
    for amp_i = 2:nAmps
        f2 = fieldnames(ALLEEG(amp_i).chanlocs);
        f = union(f,f2);
    end
    %go back and make sure all eeg sets have those fields
    for amp_i = 1:nAmps
        for field_i = 1:length(f)
            FIELD = f{field_i};
            X = [];
            if ~isfield(ALLEEG(amp_i).chanlocs,FIELD)
                %                 tempEEG(:).f(field_i) = [];
                [ALLEEG(amp_i).chanlocs.(FIELD)] = deal(X);
            end
        end
    end
    
    %## Initialize "merged" EEG structure
    EEG_out = ALLEEG(1);
    EEG_out.data = [];
    EEG_out.event = [];
    EEG_out.urevent = [];
    EEG_out.chanlocs = [];
    EEG_out.nbchan = 0;
    
    %## Find relative shifting
    %Try automatic
    disp('We are now using the accelerometers built into the EEG amps to rough align with cross-correlation');
    initialGuess = accelAlign(ALLEEG, savDir); %ROEHL - edited to pass in save directory
    disp(['Rough shift based on accelerometers = ',num2str(seconds(initialGuess)),' seconds.']);
    disp('If your results look bad and you want to manually override the auto-determined shift, then uncomment a line in this code to use the userFindShift_ALLEEG_new function'); 
    
    %Allow user to further tweak relative shifting  %If auto results not great, we can ignore suggested shift?
    %       [ relativeShift_vector ] = userFindShift_ALLEEG_new(ALLEEG, initialGuess);   
    relativeShift_vector = initialGuess; %lazy way out
    
    %## Find sync buddies across all datasets
    for eeg_i = 1:nAmps
        % A1TrigEventInd = find(strcmpi('M  1',{ALLEEG(1).event.type}));
        A_eventLat{eeg_i} = [ALLEEG(eeg_i).event.latency]; % we can lump all events into one here (M  1 or T  1 who cares) since we later check for consistency across all amps
        A_event_localTime{eeg_i} = ALLEEG(eeg_i).times(A_eventLat{eeg_i})/1000; %seconds
        % A1_event_globalTime = ALLEEG(1).recordingStart.datetime + seconds(A1_event_localTime) + seconds(0); %datetime
        A_event_globalTime{eeg_i} = ALLEEG(eeg_i).recordingStart.datetime + seconds(A_event_localTime{eeg_i}) + relativeShift_vector(eeg_i); %datetime
    end
    
    Master = A_event_globalTime{1}; %assume master is 1 but it doesn't really matter
    AllOthers = A_event_globalTime(2:end);
    thres = milliseconds(75);
    systemMatch = [];
    
    [ A_ind, B_ind ] = findSyncBuddies_multi( Master, AllOthers, thres );
    if length(A_ind)~=length(Master)
        error('uh oh. problem with findSyncBuddies_multi. ask ryan'); %2020/10/29 rjd: this shouldn't actually happen now because I fixed code in findSyncBuddies_multi, but lets check just in case
    end
    
    nTrig_A = sum(A_ind);
    nTrig_B = cellfun(@sum,B_ind);
    nTrig_all = [nTrig_A, nTrig_B]
    if any(diff(nTrig_all)) %someone disagrees
        errordlg('we may have a sync mistmatch or maybe we need to program smarter');
    else
        disp(['Looks like all systems have ',num2str(nTrig_A),' sync events']);
        systemMatch = [find(A_ind);...
            cell2mat( cellfun(@find,B_ind','UniformOutput',false) )]
    end
    
    %## Figure out which system started recording first and which ended last (note we can reasonably assume a delay of at least a couple minutes between recordings)
    for eeg_i = 1:nAmps
        recordingStartTimesAdjusted(eeg_i) = ALLEEG(eeg_i).recordingStart.datetime + relativeShift_vector(eeg_i);
        recordingStopTimesAdjusted(eeg_i) = ALLEEG(eeg_i).recordingStart.datetime + relativeShift_vector(eeg_i) + seconds(ALLEEG(eeg_i).times(end)/1000);
    end
    [earliestStartTime firstAmp] = min(recordingStartTimesAdjusted);
    [lastStopTime lastAmp] = max(recordingStopTimesAdjusted);
    
    %## Define all time where any data was recorded
    Fs = ALLEEG(1).srate; %could update to max over all eeg sets for future proofing
    
    %Option 1: select all time
    plotALLEEGevents(ALLEEG, relativeShift_vector, savDir); % %ROEHL - edited to pass in save directory
    timeOfInterest = earliestStartTime-seconds(10):seconds(1/Fs):lastStopTime+seconds(10); %add a 10 second buffer to each end so we make sure we include everything
    
    npts = length(timeOfInterest);
    EEG_out.recordingStart.datetime = timeOfInterest(1);
    
    
    %## Interpolate continuous data and events from each system to a synced set
    %Sync continuous data
    disp('Proceding to sync continuous data');
    for eeg_i = 1:nAmps
        disp(['Syncing ',num2str(eeg_i)]);
        
        Fs_DAQa         = Fs;
        t_DAQa          = timeOfInterest; %time points we want all our data to show up at
        t_DAQb          = ALLEEG(eeg_i).recordingStart.datetime + seconds(ALLEEG(eeg_i).times/1000) + relativeShift_vector(eeg_i); %time points as recorded by system b (whichever is currently being warped, technically all get warped)
        masterAmp = 1; %arbitrarily choose system 1 as our "true" global reference
        tEvents_DAQa    = A_event_globalTime{masterAmp};
        tEvents_DAQa    = tEvents_DAQa(systemMatch(masterAmp,:));
        tEvents_DAQb    = A_event_globalTime{eeg_i};
        tEvents_DAQb    = tEvents_DAQb(systemMatch(eeg_i,:));
        data_DAQb       = ALLEEG(eeg_i).data;
        
        data_warped = zeros(size(data_DAQb,1),npts);
        
        if isempty(tEvents_DAQa) %added 2021-07-15
            error('We found no sync events, cannot proceed, possibly you need to widen the automatic search for the rough alignment with the accelerometers in accelAlign.m');
        end
        %2020-04-03 new smarter method to avoid extrapolation errors?
        if length(tEvents_DAQa)>1
            if tEvents_DAQa(end)-tEvents_DAQa(1) > hours(0.5) %if we have sync signals spread far enough apart we can extrapolate
                t_DAQb_in_DAQa_base =  interp1(tEvents_DAQb,tEvents_DAQa,t_DAQb,'linear');
                extrapind = isnat(t_DAQb_in_DAQa_base);
                if any(extrapind)
                    disp('We are doing some smart extrapolation for data warping');
                else
                    disp('We did not need to do any extrapolation for data warping');
                end
                t_DAQb_in_DAQa_base(extrapind) = interp1(tEvents_DAQb([1 end]),tEvents_DAQa([1 end]),t_DAQb(extrapind),'linear','extrap');
            else %constant shift
                meanShift = mean(tEvents_DAQa-tEvents_DAQb);
                meanShift.Format = 'mm:ss.SSS'; disp(meanShift); %show in minutes:seconds.fracsec
                t_DAQb_in_DAQa_base = t_DAQb + meanShift;
                disp('WARNING DOING CONSTANT SHIFT for data warping');
            end
        else %constant shift
            meanShift = mean(tEvents_DAQa-tEvents_DAQb); %mean is meaningless (hahaha?) here since it's a scalar
            meanShift.Format = 'mm:ss.SSS'; disp(meanShift); %show in minutes:seconds.fracsec
            t_DAQb_in_DAQa_base = t_DAQb + meanShift;
            disp('WARNING DOING CONSTANT SHIFT for data warping');
        end
    
        data_warped = interp1(t_DAQb_in_DAQa_base,data_DAQb',timeOfInterest,'linear');
        data_warped = data_warped';
        
        
        EEG_out.chanlocs = [EEG_out.chanlocs, ALLEEG(eeg_i).chanlocs];
        EEG_out.nbchan = length(EEG_out.chanlocs);
        EEG_out.data = [EEG_out.data; data_warped];
    end
    EEG_out.nbchan = size(EEG_out.data,1);
    EEG_out.pnts = size(EEG_out.data,2);
    EEG_out.xmin = 0;
    EEG_out.xmax = (EEG_out.pnts-1)/EEG_out.srate + EEG_out.xmin;
    disp('Synced continuous data');
    
    %sync events
    disp('Proceding to sync events');
    for eeg_i = 1:nAmps
        %rename newest eeg set's events
        for event_i = 1:length(ALLEEG(eeg_i).event)
            %check if it was a sync event
            if any( systemMatch(eeg_i,:)== event_i) %if was a sync event
                %rename as "ampname"-LASync
                ALLEEG(eeg_i).event(event_i).type = [ALLEEG(eeg_i).setname,'-','LASync']; %add amp prefix via the setname (e.g. CGY for cortical green yellow)
                
            else %not a sync event
                %rename as "ampname"-"event"
                ALLEEG(eeg_i).event(event_i).type = [ALLEEG(eeg_i).setname,'-',ALLEEG(eeg_i).event(event_i).type]; %add amp prefix via the setname (e.g. CGY for cortical green yellow)
                
            end
        end
        
        %keep track of original system warped events came from 
        for urevent_i = 1:length(ALLEEG(eeg_i).urevent)
            ALLEEG(eeg_i).urevent(urevent_i).ursys = ALLEEG(eeg_i).setname; 
        end
        
        
        %find warped event latencies
        nAllEvents_B = length(ALLEEG(eeg_i).event);
    
        tEvents_DAQa    = A_event_globalTime{masterAmp};
        tEvents_DAQa    = tEvents_DAQa(systemMatch(masterAmp,:)); %sync events
        tEvents_DAQb    = A_event_globalTime{eeg_i};
        tEvents_DAQb    = tEvents_DAQb(systemMatch(eeg_i,:)); %sync events
        t_DAQb = A_event_globalTime{eeg_i};
        
        %2020-04-03 new smarter method to avoid extrapolation errors?
        if length(tEvents_DAQa)>1
            if tEvents_DAQa(end)-tEvents_DAQa(1) > hours(0.5) %if we have sync signals spread far enough apart we can extrapolate
                allWarpedEventTimes_global =  interp1(tEvents_DAQb,tEvents_DAQa,t_DAQb,'linear'); %no extrapolation allowed at this pont
                extrapind = isnat(allWarpedEventTimes_global); %find parts to extrapolate
                if any(extrapind)
                    disp('We are doing some smart extrapolation for event warping');
                else
                    disp('We did not need to do any extrapolation for event warping');
                end
                allWarpedEventTimes_global(extrapind) = interp1(tEvents_DAQb([1 end]),tEvents_DAQa([1 end]),t_DAQb(extrapind),'linear','extrap'); %extrapolate in a safer manner
            else %constant shift
                meanShift = mean(tEvents_DAQa-tEvents_DAQb);
                meanShift.Format = 'mm:ss.SSS'; disp(meanShift); %show in minutes:seconds.fracsec
                allWarpedEventTimes_global = t_DAQb + meanShift;
                disp('WARNING DOING CONSTANT SHIFT for event warping');
            end
        else %constant shift
            meanShift = mean(tEvents_DAQa-tEvents_DAQb); %mean is meaningless (hahaha?) here since it's a scalar
            meanShift.Format = 'mm:ss.SSS'; disp(meanShift); %show in minutes:seconds.fracsec
            allWarpedEventTimes_global = t_DAQb + meanShift;
            disp('WARNING DOING CONSTANT SHIFT for event warping');
        end
        
        
        
        %find warped urevent latencies
        nAllUREvents_B = length(ALLEEG(eeg_i).urevent);
        %                     allUnwarpedUREventTimes = ALLEEG(eeg_i).times(cell2mat({ALLEEG(eeg_i).urevent.latency}))/1000;
        %                     allWarpedUREventTimes = interp1(B_syncEvent_localTime,A_event_localTime,allUnwarpedUREventTimes,'linear','extrap'); %linear interpolation done here to avoid breaking causality (sample i+1 should always take place after sample i even  if the timing in general is bad for DAQb). Extrapolation done so we can time warp data outside the event markers
        %                     allWarpedUREventTimes_global = interp1(B_syncEvent_localTime,A_event_globalTime{masterAmp}(systemMatch(masterAmp,:)),allUnwarpedUREventTimes,'linear','extrap'); %linear interpolation done here to avoid breaking causality (sample i+1 should always take place after sample i even  if the timing in general is bad for DAQb). Extrapolation done so we can time warp data outside the event markers
        if nAllUREvents_B == nAllEvents_B
            allWarpedUREventTimes_global = allWarpedEventTimes_global; %note line above assumes urevent and event the same (no events prev deleted
        else
            if nAllUREvents_B > nAllEvents_B
                disp('I am fixing your urevent structure because i think you did some cropping already');
                ALLEEG(eeg_i).urevent = ALLEEG(eeg_i).urevent( [ALLEEG(eeg_i).event.urevent] ); %select subset of urevents so that we only have the ones that match the cropping we already did
                nAllUREvents_B = length(ALLEEG(eeg_i).urevent);
                allWarpedUREventTimes_global = allWarpedEventTimes_global; %uncomment to run anyway if you are Ryan and you know what you are doing (like you had to manually crop some data for a couple of subjects who had really long files before merging amps)
            else
            disp('We have a problem over here. You are doing things I did not plan for');
            end
        end
        %update event and urevent structures
        prevNumUREvents = length(EEG_out.urevent);
        for urevent_i = 1:nAllUREvents_B
            %copy everything over
            if isstruct(EEG_out.urevent)
                EEG_out.urevent(end+1) = ALLEEG(eeg_i).urevent(urevent_i);
            else %stupid work around for initialization issues
                EEG_out.urevent = ALLEEG(eeg_i).urevent(urevent_i); %first entry so urevent just an empty set
            end
            
            %fix latency
            %                         newLat = round(1 + EEG_out.srate * allWarpedUREventTimes(urevent_i) ); %Ryan note: could/should we remove the "round" function?
            newLat = 1 + seconds( (allWarpedUREventTimes_global(urevent_i) - EEG_out.recordingStart.datetime) ) * EEG_out.srate; %Ryan note: could/should we remove the "round" function?
            EEG_out.urevent(end).latency = newLat;
            
            
        end
        
        for event_i = 1:nAllEvents_B
            %copy everything over
            if isstruct(EEG_out.event)
                EEG_out.event(end+1) = ALLEEG(eeg_i).event(event_i);
            else %stupid work around for initialization issues
                EEG_out.event = ALLEEG(eeg_i).event(event_i);
            end
            
            %fix latency
            %                         newLat = round(1+ EEG_out.srate * allWarpedEventTimes(event_i) ); %Ryan note: could/should we remove the "round" function?
            newLat =  1 + seconds( (allWarpedEventTimes_global(event_i) - EEG_out.recordingStart.datetime) ) * EEG_out.srate; %Ryan note: could/should we remove the "round" function?
            EEG_out.event(end).latency = newLat;
            
            %update urevent number reference
            if ~isempty(EEG_out.event(end).urevent)
                EEG_out.event(end).urevent = EEG_out.event(end).urevent + prevNumUREvents;
            end
        end
    end
    EEG_out = eeg_checkset(EEG_out,'eventconsistency'); %resort by latency
    disp('Synced events');
    
    %## other
    EEG_out.setname = 'Merged LiveAmp Files';
    
    %## Tell them you finished and that you hope they didn't have to wait too long
    disp('Done merging all individual LiveAmp files into one (oh boy oh boy oh boy!)');
    toc
end


%% Changelog
%2021-08-11: RJD cleaning up some old commented out section and updating
%display to reference "network drive" rather than M drive
%2022-06-01 Ryan replaced parfor with smarter method (warp matrix instead of row by row)    
%2021-01-18 RJD: adding more feedback to the user (disp function) about the
%use of the EEG accelerometers to auto find the rough shift and what the
%rough shift was that it calculated
%2021-07-01 lines above for setref may be confusing but still works
%out. Adding commentary to explain. Each amp (CGY, CWR, NGY, NWR)
%has a single reference. This function loops thru each amp one at a
%time. All channels on a given LiveAmp64 besides the accelerometer
%channels should have the same reference (either CPz or N-CPZ
%depending on whether we are using CGY/CWR or NGY/NWR). The main
%point here is to update the reference for all channels but the
%accelerometers. Noise_chans should be empty if the current amp was
%CGY or CWR so having two lines above do not conflict with each
%other, but could be written more clearly.

