function [outputArg1,outputArg2] = mim_ls_merge_eeg(subj_fpath,data_dir)
%MIM_LS_MERGE_EEG Summary of this function goes here
%   Detailed explanation goes here
% mergeLiveAmpsWithLoadsolOneSubj.m


%%
imu_fpath = [subj_fpath filesep 'IMU' ];
eeg_raw_fpath = [subj_fpath filesep 'EEG' filesep ];
%% params and local saving
subjStr = 'H2038_FU'; %Enter Subj Name here
saveLocal = true; %Save to local computer vs newtork
disp(subjStr);

%Note: if you need to, you can temporarily change the default sync buddy
%threshold parameter from 50ms to 100ms. 
%Simply: edit autoShiftLStoEEG_func.m 
%with thresList       =	milliseconds(   [100	] ); 
%Just make sure to change it back


%% other 
if isempty(MIMDataFolderLocal)
   disp('Could not find a local folder; add path to script!');
   return;
else
   disp(['We found local folder at: ',MIMDataFolderLocal]);
   clear MIMfolderGuess; 
end
% loadFileName = [subjStr,'.set']; %original expected fileName for merged EEG
loadFileName = [subjStr,'_EEG.set']; %new expected fileName (starting 9/24/2020)

%% Initialize Diary
if saveLocal
    EEGoutputFolder = fullfile(MIMDataFolderLocal,subjStr,'EEG','Merged');
else
    EEGoutputFolder = fullfile(MIMDataFolder,subjStr,'EEG','Merged');
end

%% load live amp file
if ~isempty( fullfile(MIMDataFolder,subjStr,'EEG','Merged',loadFileName) )
    EEG = pop_loadset('filename',loadFileName,'filepath',fullfile(MIMDataFolder,subjStr,'EEG','Merged'));
else
    error('%s) Unable to find EEG file...',);
end

%% find LS files
inputFolder = fullfile(MIMDataFolder,subjStr,'Loadsol','Imported');
tempFileList = dir(fullfile(inputFolder,'*.set'));
numFiles = length(tempFileList);
if numFiles == 0
    error('There are no imported loadsol files readily availabile for this subject. Exiting the script.');
end

%% process LS files
for LSfile_i = 1:numFiles
    fileName = tempFileList(LSfile_i).name;
    
    if numFiles > 1 %Check if need subfolders for Reports
        saveDir = fullfile(DiaryOutputFolder, ['LS_', num2str(LSfile_i)]);
        mkdir(saveDir); %Create folder if appending number (new folder)
    else
        saveDir = DiaryOutputFolder;
    end
    
    %load LS file i
    LS =  pop_loadset('filename',fileName,'filepath',fullfile(MIMDataFolder,subjStr,'Loadsol','Imported'));
      
    %find gait events if they don't already exist
    if  ~any( strcmpi('RHS',{LS.event.type}) ) 
        LS =  detect_gait_events(LS);
    end
    
    
    %for loadsol only?
    LS.urevent = LS.event;
        
    
    %Update events and urevents
    for event_i = 1:length(LS.event)
        LS.event(event_i).type = ['LS-',LS.event(event_i).type]; %e.g. 'A1-M  1'
    end
    for event_i = 1:length(LS.urevent)
        %             LS.urevent(event_i).type = ['LS-',LS.urevent(event_i).type]; %e.g. 'CGY-M  1'
        LS.urevent(event_i).ursys = 'LS';
    end
    
    %check for consistency with LiveAmp event structure
    % find all possible fields
    f = fieldnames(EEG.event);
    f2 = fieldnames(LS.event);
    f = union(f,f2);
    %go back and make sure all eeg sets have those fields
    for field_i = 1:length(f)
        FIELD = f{field_i};
        X = [];
        if ~isfield(LS.event,FIELD)
            [LS.event.(FIELD)] = deal(X);
        end
        if ~isfield(EEG.event,FIELD)
            [LS.event.(FIELD)] = deal(X);
        end
    end

    
    %check for consistency with LiveAmp urevent structure
    % find all possible fields
    f = fieldnames(EEG.urevent);
   if ~isempty( LS.urevent )
    f2 = fieldnames(LS.urevent);
    f = union(f,f2);
   end
    %go back and make sure all eeg sets have those fields
    for field_i = 1:length(f)
        FIELD = f{field_i};
        X = [];
        if ~isfield(LS.urevent,FIELD)
            [LS.urevent.(FIELD)] = deal(X);
        end
        if ~isfield(EEG.urevent,FIELD)
            [EEG.urevent.(FIELD)] = deal(X);
        end
    end
    
    %make sure LS and EEGs have the same fields for their CHANLOCS structure
    % find all possible fields
    f = fieldnames(EEG.chanlocs);
    f2 = fieldnames(LS.chanlocs);
    f = union(f,f2);
    %go back and make sure all eeg sets have those fields
    for field_i = 1:length(f)
        FIELD = f{field_i};
        X = [];
        if ~isfield(LS.chanlocs,FIELD)
            [LS.chanlocs.(FIELD)] = deal(X);
        end
        if ~isfield(EEG.chanlocs,FIELD)
            [EEG.chanlocs.(FIELD)] = deal(X);
        end
    end
    
    
    %generate info for events to use later
    uniqueEvents_A = unique({EEG.event.type});
    uniqueEvents_B = unique({LS.event.type});
    
    
    %ask user to pick sync signal on "master" eeg set
    disp(newline);
    for name_i = 1:length(uniqueEvents_A)
        disp([num2str(name_i),': ',uniqueEvents_A{name_i},' (',num2str( sum(strcmpi(uniqueEvents_A{name_i},{EEG.event.type})) ),' events)']);
    end
    disp('Please choose master sync signal from list above.');
    disp('Hint: it is probably CWR-M 1 or the one with the most events.');
    prompt = 'Master sync signal: ';
    eventTypeA = uniqueEvents_A{input(prompt)};

    %ask user to pick sync signal for newest eeg set we are joining to the master set
    %will automatically select LS if possible
    if length(find(strcmp(uniqueEvents_B, 'LS-Sync'))) == 1
        eventTypeB = uniqueEvents_B(find(strcmp(uniqueEvents_B, 'LS-Sync'))); %auto select LS
    else %have to do it manually
        disp(newline);
        prompt = ['Please choose sync signal for newest LS set: '];
        for name_i = 1:length(uniqueEvents_B)
            disp([num2str(name_i),': ',uniqueEvents_B{name_i},' (',num2str( sum(strcmpi(uniqueEvents_B{name_i},{LS.event.type})) ),' events).'])
        end
        eventTypeB = uniqueEvents_B{input(prompt)};
    end
    
    ATrigEventInd = find(strcmpi(eventTypeA,{EEG.event.type}));
    BTrigEventInd = find(strcmpi(eventTypeB,{LS.event.type}));

    A_eventLat = [EEG.event(ATrigEventInd).latency];
    B_eventLat = [LS.event(BTrigEventInd).latency];

    A_event_localTime = (A_eventLat-1)/EEG.srate; %seconds
    B_event_localTime = (B_eventLat-1)/LS.srate; %seconds
    
%     %Error structure
%     allSubjSyncLStoEEG(20).A_event_localTime = A_event_localTime; 
%     allSubjSyncLStoEEG(20).B_event_localTime = B_event_localTime; 
%     allSubjSyncLStoEEG(20).subj = subjStr;
%     save('allSubjSyncLStoEEG.mat', 'allSubjSyncLStoEEG')
%     
    trimPerc = 20;
    timeBetweenSyncEventsA =diff(A_event_localTime);
    timeBetweenSyncEventsB =diff(B_event_localTime);
    tempTrimA = trimmean(timeBetweenSyncEventsA,trimPerc);
    tempTrimB = trimmean(timeBetweenSyncEventsB,trimPerc);
%     figure; stem(diff(A_event_localTime));
    tempEquivLS_Fs = tempTrimB/tempTrimA*LS.srate; %estimate of what the correctec srate should be for the LS
   
    disp(['EEG trimmed mean diff between sync events = ',num2str(tempTrimA),' seconds']);
    disp(['Loadsol trimmed mean diff between sync events = ',num2str(tempTrimB),' seconds']);
%     disp(['Equiv sample rate you should set LS to help with alignment = ',num2str(tempEquivLS_Fs),' Hz']); %2021-01-21 RJD: this line isn't needed anymore. It was from a time when we thought LS might have time dilation issues but really it was because a bug in our EEG merging code caused time dilation in the merged EEG stream
    
    A_event_globalTime = EEG.recordingStart.datetime + seconds(A_event_localTime); %datetime
    B_event_globalTime = LS.recordingStart.datetime + seconds(B_event_localTime); %datetime
    
%     [tempTimeScaleFactor, tempShift] = testAutoStretchAndShiftEventsToLineUp(A_event_globalTime,B_event_globalTime);
    
     %Find rough shift to align systems %COLTON LOOK HERE%
     disp('You will now be asked to manually shift the LoadSol events relative to the EEG events to get them to roughly align.');
     disp('Manually enter non-zero numbers to shift the events by that amount (in seconds)');
     disp('Enter 0 when you are happy with the rough alignment');
    [ roughShift ] = userFindShift(A_event_globalTime, B_event_globalTime ) 
  %   roughShift = [];
    disp('Thank you for providing manual input. We are now using your estimate of the rough alignment and fine tuning it automatically');
    [ bestShift, bestThres ] = autoShiftLStoEEG_func( A_event_globalTime, B_event_globalTime, roughShift, saveDir)
%      [ A_ind, B_ind ] = findSyncBuddies( A_event_globalTime, B_event_globalTime+bestShift, milliseconds(50) );
        [ A_ind, B_ind ] = findSyncBuddies( A_event_globalTime, B_event_globalTime+bestShift, bestThres );
     if sum(A_ind) ~= sum(B_ind)
         error('Something has caused us to find more sync buddies in one system than another so I cannot perfectly pair them anymore');
     end
     %%%%%%%%%%%%%warp continuous data%%%%%%%%%%%%%%%%%%%%%%%%
      disp(['Syncing continuous data. LS file ',num2str(LSfile_i)]);
    
    Fs_DAQa         = EEG.srate;
    t_DAQa          = EEG.recordingStart.datetime + seconds(EEG.times/1000); %time points we want all our data to show up at
    t_DAQb          = LS.recordingStart.datetime + seconds(LS.times/1000); %time points as recorded by system b (whichever is currently being warped, technically all get warped)
    masterAmp = 1; %arbitrarily choose system 1 as our "true" global reference
    tEvents_DAQa    = A_event_globalTime(A_ind);
    tEvents_DAQb    = B_event_globalTime(B_ind);
    data_DAQb       = LS.data;
    timeOfInterest = t_DAQa;
    npts = length(timeOfInterest);
    data_warped = zeros(size(data_DAQb,1),npts);
    
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
        if isnan(meanShift)
            error('something bad happened')
        else
        t_DAQb_in_DAQa_base = t_DAQb + meanShift;
        disp('WARNING DOING CONSTANT SHIFT for data warping');
        end
    end
    
    parfor row_i = 1:size(data_DAQb,1)
        %    data_interp(row_i, :) =  interp1(t_DAQb_in_DAQa_base,data_DAQb(row_i,:),t_DAQb_match_DAQa_sampling,'spline');
        data_warped (row_i, :) =  interp1(t_DAQb_in_DAQa_base,data_DAQb(row_i,:),timeOfInterest,'linear');
    end
    
    
    if LSfile_i == 1
    EEG.chanlocs = [EEG.chanlocs, LS.chanlocs]; %append chanlocs
    EEG.nbchan = length(EEG.chanlocs); %updating number of channels since we appended new LS channels
    EEG.data = [EEG.data; data_warped]; %appending new LS channels
    else
        numLSch = length(LS.chanlocs); %assumes all loadsol files have same exact num of data chans
        for LSch_i = 1:numLSch
        equivEEGch =  EEG.nbchan - numLSch + LSch_i;  %2021-04-20 RJD tried to fix bug with indexing when mergining 2+ LS files
%         timeWindowToReplace = find( isnan(EEG.data(end-LSch_i+1,:)) & ~isnan(data_warped(LSch_i,:)) ); %find where previous iteration didn't warp any data and we have new data avail to enter this iteration (file num)
        timeWindowToReplace = find( isnan(EEG.data(equivEEGch,:)) & ~isnan(data_warped(LSch_i,:)) ); %2021-04-20 RJD tried to fix bug %find where previous iteration didn't warp any data and we have new data avail to enter this iteration (file num)

%         EEG.data(end-LSch_i+1,timeWindowToReplace) = data_warped(LSch_i,timeWindowToReplace);
        EEG.data(equivEEGch,timeWindowToReplace) = data_warped(LSch_i,timeWindowToReplace);%2021-04-20 RJD tried to fix bug
        end
    end
 
    disp('Synced continuous LS data');
    
    
    %%%%%%%%%%%%%%%%%warp events%%%%%%%%%%%%%%%%%%%%%%%%
    
    %update names in LiveAmp merged set
    for event_i = 1:length(EEG.event)
        %check if it was a sync event
        if any( ATrigEventInd(A_ind) == event_i) %if was a sync event
            %rename as EEG-LSSync
            EEG.event(event_i).type = ['EEG','-','LSSync']; %add amp prefix via the setname (e.g. CGY for cortical green yellow)
            
        end
    end
    
    
  
    for event_i = 1:length(LS.event)
        %check if it was a sync event
        if any( BTrigEventInd(B_ind) == event_i) %if was a sync event
            %rename as LS-LSSync
            LS.event(event_i).type = ['LS','-','LSSync']; %add amp prefix via the setname (e.g. CGY for cortical green yellow)
            
        else %not a sync event
%             %rename as "ampname"-"event"
%             LS.event(event_i).type = ['LS','-',LS.event(event_i).type]; %add system prefix (should be none for LS data unless it already includes RHS/LHS/etc
%             %consider removing line above if don't want certain events to
%             %be overwritten like LHS/RHS...
        end
    end
    

    %keep track of original system warped events came from 
    for urevent_i = 1:length(LS.urevent)
        LS.urevent(urevent_i).ursys = 'LoadSol'; %possibly differentiate file 1 vs file 2 etc?
    end
    
     
    
   
     
    %find warped event latencies
    nAllEvents_B = length(LS.event);
    
    tEvents_DAQa    = A_event_globalTime(A_ind); %sync events
    tEvents_DAQb    = B_event_globalTime(B_ind);%sync events

     allB_eventLat = [LS.event.latency];
     allB_event_localTime = (allB_eventLat-1)/LS.srate; 
     allB_event_globalTime = LS.recordingStart.datetime + seconds(allB_event_localTime); %datetime
    t_DAQb = allB_event_globalTime;
    
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
    nAllUREvents_B = length(LS.urevent);
    if nAllUREvents_B == nAllEvents_B
        allWarpedUREventTimes_global = allWarpedEventTimes_global; %note line above assumes urevent and event the same (no events prev deleted
    else
        disp('We have a problem over here. You are doing things I did not plan for');
    end
    %update event and urevent structures
    prevNumUREvents = length(EEG.urevent);
    for urevent_i = 1:nAllUREvents_B
        %copy everything over
        if isstruct(EEG.urevent)
            EEG.urevent(end+1) = LS.urevent(urevent_i);
        else %stupid work around for initialization issues
            EEG.urevent = LS.urevent(urevent_i); %first entry so urevent just an empty set
        end
        
        %fix latency
        %                         newLat = round(1 + EEG_out.srate * allWarpedUREventTimes(urevent_i) ); %Ryan note: could/should we remove the "round" function?
        newLat = 1 + seconds( (allWarpedUREventTimes_global(urevent_i) - EEG.recordingStart.datetime) ) * EEG.srate; %Ryan note: could/should we remove the "round" function?
        EEG.urevent(end).latency = newLat;
    end
    
    for event_i = 1:nAllEvents_B
        %copy everything over
        if isstruct(EEG.event)
            EEG.event(end+1) = LS.event(event_i);
        else %stupid work around for initialization issues
            EEG.event = LS.event(event_i);
        end
        
        %fix latency
        %                         newLat = round(1+ EEG_out.srate * allWarpedEventTimes(event_i) ); %Ryan note: could/should we remove the "round" function?
        newLat =  1 + seconds( (allWarpedEventTimes_global(event_i) - EEG.recordingStart.datetime) ) * EEG.srate; %Ryan note: could/should we remove the "round" function?
        EEG.event(end).latency = newLat;
        
        %update urevent number reference
        if ~isempty(EEG.event(end).urevent)
            EEG.event(end).urevent = EEG.event(end).urevent + prevNumUREvents;
        end
    end
     
     EEG = eeg_checkset(EEG,'eventconsistency'); %resort by latency
     disp('Synced events');
    
end


EEG.setname = 'Merged LiveAmp and Loadsol Files';


%Add global times to events to make it easier to look at manually
EEG = addGlobalTimesToEvents(EEG);

%% Checking accelerometers 
AccInd = find(strcmpi({EEG.chanlocs.type},'Acc')); %finding where the accelerometers are
tempEEG = pop_select(EEG, 'channel', AccInd); %selecting only these channels
%setting bounds of movement data from 0 to 1
Resultant = sqrt(sum(tempEEG.data.^2,1)); %taking magnitude of xyz vector
ScaleFactor = prctile(Resultant,99.5);
Resultant_Norm = Resultant/ScaleFactor; %normalizing by dividing my 99.5% of maximum value

globalEEGtimes = tempEEG.recordingStart.datetime + milliseconds(tempEEG.times); %grabbing global times


for event_i = 1:length(tempEEG.event)
    EEG.event(event_i).datetime = tempEEG.recordingStart.datetime + seconds( (tempEEG.event(event_i).latency - 1)/tempEEG.srate);
    EEG.event(event_i).datestr = datestr(tempEEG.event(event_i).datetime);
end


%finding events
RTOind = find(strcmpi({tempEEG.event.code},'RTO')); %might switch to event.type in the future
LTOind = find(strcmpi({tempEEG.event.code},'LTO')); %Left Toe Off indices
RHSind = find(strcmpi({tempEEG.event.code},'RHS')); %Right Heel Strike indices, etc.
LHSind = find(strcmpi({tempEEG.event.code},'LHS'));
 
myFig = figure;
stem([EEG.event(RTOind).datetime], ones(size(RTOind)),'b'); hold on; %marking when RTO happens
stem([EEG.event(LTOind).datetime], ones(size(LTOind)),'r'); %when LTO happens
stem([EEG.event(RHSind).datetime], ones(size(RHSind)),'g'); %when RHS happens
stem([EEG.event(LHSind).datetime], ones(size(LHSind)),'k'); %when LHS happens
plot(globalEEGtimes, Resultant_Norm) %then add the accelerometer resultant calculated before
legend('RTO', 'LTO', 'RHS', 'LHS', 'Resultant');
saveas(myFig,fullfile(DiaryOutputFolder, 'AccelGaitEventAlign.fig')); %saving figure showing resultant accel signal with gait events overlay

%% Save data

EEG.setname = [subjStr,'_EEGandLS'];
fileNameNoExt = [subjStr,'_EEGandLS'];
if saveLocal
    disp('Saving locally because hipergator is slow');
    outputFolder = fullfile(MIMDataFolderLocal,subjStr,'EEG','Merged');
else
    outputFolder = fullfile(MIMDataFolder,subjStr,'EEG','Merged');
end
mkdir(outputFolder);

% datetime('20190911-152201','InputFormat','yyyyMMdd-HHmmss')
outputFileName = [fileNameNoExt,'.set'];
disp(['Saving: ', fullfile(outputFolder,outputFileName)]);
EEG = pop_saveset( EEG, 'filepath', outputFolder, 'filename', fileNameNoExt);
disp('Done saving.');

if isempty(CURRENTSET)
    CURRENTSET = 0;
end
    
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab redraw; 

%% End diary
diary off


%% Tell user to copy files if needed
disp(newline);
disp('Finished merging loadsol to live amps.')
if saveLocal
    disp(newline)
    disp('You chose to save the data and figures to a local folder temporarily (to save time with slow M drive). Please cut and paste the files for this subject from the local folder to the M drive.');
    disp('Note if you are having trouble copying the diary file over to the M drive (error saying matlab has it open somehow), then please close matlab and retry). You can also let Ryan know you had issues with the diary so he can try to figure out a fix');
else
    disp(newline);
    disp('You saved to the M drive so there is nothing left for you to do.');
end
end

