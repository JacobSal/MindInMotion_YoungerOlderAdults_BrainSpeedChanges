function [outputArg1,outputArg2] = imu_plots(inputArg1,inputArg2)
%IMU_PLOTS Summary of this function goes here


%% Optionally Plot 3D IMU trajectory
% movie_samplePeriod  = .1;%.02 default %1/EEG.srate;

%plot in global frame (don't have to deal with trying to figure out how to
%rotate the quatPlot to show everything in the body frame of reference)
posPlot = pos_est; %5-14-21 we could crop here but might get discontinuity in the video and we don't care the video is just for fun
quatPlot = quaternConj(orient_data); %not sure why need to use conjugate but it looks correct compared to other way


% %2020-01-05 Note: temp crop (remove later) %RD commented out on
% 2021-03-24
% posPlot = posPlot(5000:end-5000,:);
% quatPlot = quatPlot(5000:end-5000,:);


% Extend final sample to delay end of animation
extraTime = 1;
onesVector = ones(extraTime*(1/movie_samplePeriod), 1);
posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

% Create 6 DOF animation
SamplePlotFreq = 10;%20 for 500 Hz to appear normal speed?

Spin =120*2;
spinMat = [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)]; %original
% spinMat = [0*ones(length(posPlot),1),0*ones(length(posPlot), 1)];
%0,0 x= horiz, z = vert (side or frontal)
%90,0 y = horiz, z = vert (side or fontal)
%90,-90 x=vert, y = horiz (top down view?)
%-90,90 north up, east right

%note trail could be 'Off' 'DotsOnly' or 'All'
if doMovie
    mySixDofAnimation(posPlot, quatern2rotMat(quatPlot), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'DotsOnly', ...
        'Position', [9 39 1280 768], 'View', spinMat, ...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
        'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/movie_samplePeriod) / SamplePlotFreq));
end

%% make static image of entire timecourse of trial from 3 orthogonal angles
figure;
if DoCrop
    ExactCropLatency = round(EEG.srate*ExactCrop+1);
    nPairsTemp = size(ExactCrop);
    nPairsReal = nPairsTemp(1); %roundabout way but length does work for 2x2 matrix as well
    
    for pair_i = 1:1:nPairsReal
        UsefulCropLatencies = ExactCropLatency(pair_i,1):ExactCropLatency(pair_i,2);
        subplot(1,3,1); plot(-pos_est_VAF(UsefulCropLatencies,2),pos_est_VAF(UsefulCropLatencies,1),':'); if pair_i < nPairsReal; hold on; end; xlabel('Right'); ylabel('Foward'); title('Transverse'); %top down 
        subplot(1,3,2); plot(pos_est_VAF(UsefulCropLatencies,1),pos_est_VAF(UsefulCropLatencies,3),':'); if pair_i < nPairsReal; hold on; end; xlabel('Forward'); ylabel('Up'); title('Sagittal'); %side view
        subplot(1,3,3); plot(-pos_est_VAF(UsefulCropLatencies,2),pos_est_VAF(UsefulCropLatencies,3),':'); if pair_i < nPairsReal; hold on; end; xlabel('Right'); ylabel('Up'); title('Coronal'); %from behind or in front
    end
else 
    subplot(1,3,1); plot(-pos_est_VAF(:,2),pos_est_VAF(:,1),':'); xlabel('Right'); ylabel('Foward'); title('Transverse'); %top down 
    subplot(1,3,2); plot(pos_est_VAF(:,1),pos_est_VAF(:,3),':'); xlabel('Forward'); ylabel('Up'); title('Sagittal'); %side view
    subplot(1,3,3); plot(-pos_est_VAF(:,2),pos_est_VAF(:,3),':'); xlabel('Right'); ylabel('Up'); title('Coronal'); %from behind or in front
end

if doSave
    tempSubFolder = 'Orthogonal_View_Figs';
    if ~exist(fullfile(OutputFolder,tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder,tempSubFolder));
    end
    saveas(gcf, fullfile(OutputFolder,tempSubFolder, saveName));
    close
end

%% make new EEG set with A/P M/L vert positions 
fakeEEG = pop_select(tempEEG,'channel',1:12); %easy way to start a new eeg structure
fakeEEG.data = pos_est'; fakeEEG.pnts = size(fakeEEG.data,2);


%no fancy pca. use estimate of average Front, Left, Up frame 
fakeEEG.data(1:3,:) = pos_est_VAF';
fakeEEG.data(4:6,:) = vel_est_VAF';
fakeEEG.data(7:9,:) = acc_est_VAF';
fakeEEG.data(10:12,:) = nan([3, size(orient_data,1)]); %fill with nan for now because I don't know what this would be. We probably don't need it anyway

fakeEEG.chanlocs(1).labels = 'xPos';
fakeEEG.chanlocs(2).labels = 'yPos';
fakeEEG.chanlocs(3).labels = 'zPos';
fakeEEG.chanlocs(4).labels = 'xVel';
fakeEEG.chanlocs(5).labels = 'yVel';
fakeEEG.chanlocs(6).labels = 'zVel';
fakeEEG.chanlocs(7).labels = 'xAcc';
fakeEEG.chanlocs(8).labels = 'yAcc';
fakeEEG.chanlocs(9).labels = 'zAcc';
fakeEEG.chanlocs(10).labels = 'xRot';
fakeEEG.chanlocs(11).labels = 'yRot';
fakeEEG.chanlocs(12).labels = 'zRot';

disp('taking a break');

%add jerk
jerk_VAF = [diff(acc_est_VAF); 0 0 0]/(1/fakeEEG.srate); %take time derivative of acceleration to get jerk (in VAF frame)
jerk_resultant = sqrt( sum( jerk_VAF.^2,2 )); %calc resultant jerk
fakeEEG.data([13 14 15],:) = jerk_VAF'; 
fakeEEG.data(16,:) = jerk_resultant';
fakeEEG.chanlocs(13).labels = 'xJerk';
fakeEEG.chanlocs(14).labels = 'yJerk';
fakeEEG.chanlocs(15).labels = 'zJerk';
fakeEEG.chanlocs(16).labels = 'rJerk'; %rJerk = resultant Jerk



%% Time warp data stride by stride keeping only good strides
%Copied from timewarpIMUERP (that older script probably no longer needed)
tempEEG = fakeEEG;

gaitEvents = tempEEG.event( strcmpi('RHS',{tempEEG.event.type}) | strcmpi('RTO',{tempEEG.event.type}) | strcmpi('LHS',{tempEEG.event.type}) | strcmpi('LTO',{tempEEG.event.type})); %grab just right heel strike info

RHSevents = gaitEvents( strcmpi('RHS',{gaitEvents.type}) ); %grab just right heel strike info
RHSind_relAllGaitEvents = find(strcmpi({gaitEvents.type},'RHS')); %finds event number relative to all other events for these RHS

tDes = 0:.1:100;%Desired 0 to 100 percent of gait cycle at a fixed rate
nStrides = length(RHSevents)-1;
nTpts = length(tDes);
nCh = size(tempEEG.data,1);
warpedEpochedData = NaN(nCh,nTpts,nStrides);
goodStrides = false(1,nStrides);


% rejecting bad strides (code had issue with h2002 0p25_2
for stride_i = 1:(nStrides)
    %Determine if this stride is a valid stride (gait events follow expected number
    %and ordering)
   % indGaitEventsInThisSTride = []; ContainedEventTypes = [];
    CurrRHSind = RHSind_relAllGaitEvents(stride_i); %current RHS
    NextRHSind = RHSind_relAllGaitEvents(stride_i + 1); %next RHS (does not exceed length of variable based on definition of stride_i)
    indGaitEventsInThisStride = find([gaitEvents.latency] > [gaitEvents(CurrRHSind).latency] & [gaitEvents.latency] < [gaitEvents(NextRHSind).latency]); %indices sandwiched between bounds
    ContainedEventTypes = {gaitEvents(indGaitEventsInThisStride).type}; %array of names that these indices correspond to
    if length(ContainedEventTypes) == 3 %Expect to find LTO, LHS, RTO in between each pair of RHS
        if all(strcmpi(ContainedEventTypes,{'LTO','LHS','RTO'})) 
            goodStrides(stride_i) = true; %this stride is good
        end
    end
    
    tStartStride = (RHSevents(stride_i).latency-1)/tempEEG.srate;
    tEndStride = (RHSevents(stride_i+1).latency-1)/tempEEG.srate;  
    if DoCrop
        MaxLatency = (RHSevents(end).latency+1); %creates an artificially larger timeset so we don't exceed dimensions
        CropBoolean = zeros(1,fix(MaxLatency)); %start the boolean at 0s for the size of this larger set
        ExactCropLatencies = EEG.srate*ExactCrop+1; %converting time to latencies 
        tempPairs = size(ExactCrop); %if we have stranger cropping this matters
        for set_i = 1:1:tempPairs(1) %grabbing the first num is really the number of pairs, they are N x 2, like [0 10] or [0 10; 30 50], etc.
            LatenciesOfInterest = ExactCropLatencies(set_i,:); %looking at pairs one by one. first [0 10] then [30 50]
            CropBoolean(LatenciesOfInterest(1):LatenciesOfInterest(2)) = 1; %setting those latencies true in the boolean
        end
     
        StrideLatencies = [tStartStride tEndStride]*EEG.srate+1; %convert to latencies
        StrideLatencySpan = StrideLatencies(1):StrideLatencies(2); %just changing to a start:stop format 
        %the logic is that the boolean tells us all the good things as 1. if
        %we look at the boolean during the latencies of this stride, and
        %they are all 1, then the stride is definitely good. if that's
        %true, then the "unique" output will give 1 as the only thing.
        %otherwise there will be at least one 0, so the unique part won't
        %be just 1.
        %--
        %sometimes matlab liked to put the latency in scientific notation
        %which annoyed the indexing part, so using fix "rounds it" to the
        %same thing. all the latencies are integers so this should not
        %cause problems.
        if all(CropBoolean(fix(StrideLatencySpan)))
            %we are happy since this falls in cropping range
            %might be smarter to not set true in case they do not satisfy
            %the above 
        else
            goodStrides(stride_i) = false; %kick stride out based on the Problematic IMU input
        end
    end
    
    %warp data (whether good or not. we will deal with bad strides later)
    tInd = find(tempEEG.times/1000 >= tStartStride & tempEEG.times/1000 <= tEndStride);
    if ~isempty(tInd)
        tempData = tempEEG.data(:,tInd); %grab data for this stride
        tAct = linspace(0,100,length(tInd)); %actual time our samples occured at when expressed as percentage of gait cycle
        tempWarped = interp1(tAct, tempData',tDes)';
        warpedEpochedData(:,:,stride_i) = tempWarped;
    else
        %do nothing. leave as NaN
    end
end
%keep only good strides
%disp(['Removing bad strides: ', find(~goodStrides)]); %this line glitches??
fprintf('Removing bad strides: '); disp(find(~goodStrides)); %random fix (:
warpedEpochedData = warpedEpochedData(:,:,goodStrides); %selecting good strides only

warpedEpochedData(2,:,:) = -warpedEpochedData(2,:,:); %flip left to right to match another paper
warpedEpochedData(5,:,:) = -warpedEpochedData(5,:,:); %flip left to right to match another paper

disp(newline);
disp('Finished warping IMU data stride by stride');

%% Calc Stride by Stride Measures 
%Note: 3 positions and 3 velocities for indexing
SDcutoff = 2.5; %2021-03-24 RD: changed from 3 to 2.5 to match what we are doing with LS analysis

%this method didn't work after looking at it for 1 hour for weird indexing
%stuff so I'm just hard coding for now.
% for nSet = 1:6 %for 6 conditions
%     [outlierBoolArrays] = MakeOutlierBoolAndStructs(warpedEpochedData,SDcutoff,nSet);
% end

APcomExcAllStrides = squeeze(range(warpedEpochedData(1,:,:))); %x dir. AP
tempMean = mean(APcomExcAllStrides);
tempSD = std(APcomExcAllStrides);
inputArray = APcomExcAllStrides;
outlierBoolArrayAPexc = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
figure; %sgtitle([SubjStr,' ',saveName]);
subplot(3,4,1); stem(1:length(APcomExcAllStrides),APcomExcAllStrides); hold on; stem(find(outlierBoolArrayAPexc),APcomExcAllStrides(outlierBoolArrayAPexc),'r'); title('P2P AP Exc (m)'); xlabel('Stride Num');
realMean = mean(APcomExcAllStrides(~outlierBoolArrayAPexc));
realSD = std(APcomExcAllStrides(~outlierBoolArrayAPexc));
realCOV = 100*realSD/realMean;
myStructure.APexc_mean = realMean;
myStructure.APexc_COV = realCOV;
StrideByStrideStruct.APexc = APcomExcAllStrides(~outlierBoolArrayAPexc);

MLcomExcAllStrides = squeeze(range(warpedEpochedData(2,:,:))); %y dir. ML
tempMean = mean(MLcomExcAllStrides);
tempSD = std(MLcomExcAllStrides);
inputArray = MLcomExcAllStrides;
outlierBoolArrayMLexc = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
subplot(3,4,2); stem(1:length(MLcomExcAllStrides),MLcomExcAllStrides); hold on; stem(find(outlierBoolArrayMLexc),MLcomExcAllStrides(outlierBoolArrayMLexc),'r'); title('P2P ML Exc (m)'); xlabel('Stride Num');
realMean = mean(MLcomExcAllStrides(~outlierBoolArrayMLexc));
realSD = std(MLcomExcAllStrides(~outlierBoolArrayMLexc));
realCOV = 100*realSD/realMean;
myStructure.MLexc_mean = realMean;
myStructure.MLexc_COV = realCOV;
StrideByStrideStruct.MLexc = MLcomExcAllStrides(~outlierBoolArrayMLexc);

UDcomExcAllStrides = squeeze(range(warpedEpochedData(3,:,:))); %z dir. UD
tempMean = mean(UDcomExcAllStrides);
tempSD = std(UDcomExcAllStrides);
inputArray = UDcomExcAllStrides;
outlierBoolArrayUDexc = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
subplot(3,4,3); stem(1:length(UDcomExcAllStrides),UDcomExcAllStrides); hold on; stem(find(outlierBoolArrayUDexc),UDcomExcAllStrides(outlierBoolArrayUDexc),'r'); title('P2P UpDown Exc (m)'); xlabel('Stride Num');
realMean = mean(UDcomExcAllStrides(~outlierBoolArrayUDexc));
realSD = std(UDcomExcAllStrides(~outlierBoolArrayUDexc));
realCOV = 100*realSD/realMean;
myStructure.UDexc_mean = realMean;
myStructure.UDexc_COV = realCOV;
StrideByStrideStruct.UDexc = UDcomExcAllStrides(~outlierBoolArrayUDexc);

PeakAntVelAllStrides = squeeze(range(warpedEpochedData(4,:,:))); %x dir. AP
tempMean = mean(PeakAntVelAllStrides);
tempSD = std(PeakAntVelAllStrides);
inputArray = PeakAntVelAllStrides;
outlierBoolArrayAPvel = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
subplot(3,4,5); stem(1:length(PeakAntVelAllStrides),PeakAntVelAllStrides); hold on; stem(find(outlierBoolArrayAPvel),PeakAntVelAllStrides(outlierBoolArrayAPvel),'r'); title('P2P AP Vel (m/s)'); xlabel('Stride Num');
realMean = mean(PeakAntVelAllStrides(~outlierBoolArrayAPvel));
realSD = std(PeakAntVelAllStrides(~outlierBoolArrayAPvel));
realCOV = 100*realSD/realMean;
myStructure.PeakAntVel_mean = realMean;
myStructure.PeakAntVel_COV = realCOV;
StrideByStrideStruct.PeakAntVel = PeakAntVelAllStrides(~outlierBoolArrayAPvel);

PeakLatVelAllStrides = squeeze(range(warpedEpochedData(5,:,:))); %y dir. ML
tempMean = mean(PeakLatVelAllStrides);
tempSD = std(PeakLatVelAllStrides);
inputArray = PeakLatVelAllStrides;
outlierBoolArrayMLvel = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
subplot(3,4,6); stem(1:length(PeakLatVelAllStrides),PeakLatVelAllStrides); hold on; stem(find(outlierBoolArrayMLvel),PeakLatVelAllStrides(outlierBoolArrayMLvel),'r'); title('P2P ML Vel (m/s)'); xlabel('Stride Num');
realMean = mean(PeakLatVelAllStrides(~outlierBoolArrayMLvel));
realSD = std(PeakLatVelAllStrides(~outlierBoolArrayMLvel));
realCOV = 100*realSD/realMean;
myStructure.PeakLatVel_mean = realMean;
myStructure.PeakLatVel_COV = realCOV;
StrideByStrideStruct.PeakLatVel = PeakLatVelAllStrides(~outlierBoolArrayMLvel);

PeakUDVelAllStrides = squeeze(range(warpedEpochedData(6,:,:))); %z dir. UD
tempMean = mean(PeakUDVelAllStrides);
tempSD = std(PeakUDVelAllStrides);
inputArray = PeakUDVelAllStrides;
outlierBoolArrayUDvel = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
subplot(3,4,7); stem(1:length(PeakUDVelAllStrides),PeakUDVelAllStrides); hold on; stem(find(outlierBoolArrayUDvel),PeakUDVelAllStrides(outlierBoolArrayUDvel),'r'); title('P2P UpDown Vel (m/s)'); xlabel('Stride Num');
realMean = mean(PeakUDVelAllStrides(~outlierBoolArrayUDvel));
realSD = std(PeakUDVelAllStrides(~outlierBoolArrayUDvel));
realCOV = 100*realSD/realMean;
myStructure.PeakUpDownVel_mean = realMean;
myStructure.PeakUpDownVel_COV = realCOV;
StrideByStrideStruct.PeakUpDownVel = PeakUDVelAllStrides(~outlierBoolArrayUDvel);


%Jerk
p2pAPJerkAllStrides = squeeze(range(warpedEpochedData(13,:,:))); %x dir. AP
    tempMean = mean(p2pAPJerkAllStrides);
    tempSD = std(p2pAPJerkAllStrides);
    inputArray = p2pAPJerkAllStrides;
    outlierBoolArrayAPjerk = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
    subplot(3,4,9); stem(1:length(p2pAPJerkAllStrides),p2pAPJerkAllStrides); hold on; stem(find(outlierBoolArrayAPjerk),p2pAPJerkAllStrides(outlierBoolArrayAPjerk),'r'); title('P2P AP Jerk (m/sec^3)'); xlabel('Stride Num');
    realMean = mean(p2pAPJerkAllStrides(~outlierBoolArrayAPjerk));
    realSD = std(p2pAPJerkAllStrides(~outlierBoolArrayAPjerk));
    realCOV = 100*realSD/realMean;
myStructure.APJerk_mean = realMean;
myStructure.APJerk_COV = realCOV;
StrideByStrideStruct.APJerk = p2pAPJerkAllStrides(~outlierBoolArrayAPjerk);

p2pMLJerkAllStrides = squeeze(range(warpedEpochedData(14,:,:))); %y dir. ML
    tempMean = mean(p2pMLJerkAllStrides);
    tempSD = std(p2pMLJerkAllStrides);
    inputArray = p2pMLJerkAllStrides;
    outlierBoolArrayMLjerk = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
    subplot(3,4,10); stem(1:length(p2pMLJerkAllStrides),p2pMLJerkAllStrides); hold on; stem(find(outlierBoolArrayMLjerk),p2pMLJerkAllStrides(outlierBoolArrayMLjerk),'r'); title('P2P ML Jerk (m/sec^3)'); xlabel('Stride Num');
    realMean = mean(p2pMLJerkAllStrides(~outlierBoolArrayMLjerk));
    realSD = std(p2pMLJerkAllStrides(~outlierBoolArrayMLjerk));
    realCOV = 100*realSD/realMean;
myStructure.MLJerk_mean = realMean;
myStructure.MLJerk_COV = realCOV;
StrideByStrideStruct.MLJerk = p2pMLJerkAllStrides(~outlierBoolArrayMLjerk);

p2pUDJerkAllStrides = squeeze(range(warpedEpochedData(15,:,:))); %x dir. AP
    tempMean = mean(p2pUDJerkAllStrides);
    tempSD = std(p2pUDJerkAllStrides);
    inputArray = p2pUDJerkAllStrides;
    outlierBoolArrayUDjerk = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
    subplot(3,4,11); stem(1:length(p2pUDJerkAllStrides),p2pUDJerkAllStrides); hold on; stem(find(outlierBoolArrayUDjerk),p2pUDJerkAllStrides(outlierBoolArrayUDjerk),'r'); title('P2P UpDown Jerk (m/sec^3)'); xlabel('Stride Num');
    realMean = mean(p2pUDJerkAllStrides(~outlierBoolArrayUDjerk));
    realSD = std(p2pUDJerkAllStrides(~outlierBoolArrayUDjerk));
    realCOV = 100*realSD/realMean;
myStructure.UDJerk_mean = realMean;
myStructure.UDJerk_COV = realCOV;
StrideByStrideStruct.UDJerk = p2pUDJerkAllStrides(~outlierBoolArrayUDjerk);


p2pResJerkAllStrides = squeeze(range(warpedEpochedData(16,:,:))); %x dir. AP
    tempMean = mean(p2pResJerkAllStrides);
    tempSD = std(p2pResJerkAllStrides);
    inputArray = p2pResJerkAllStrides;
    outlierBoolArrayResjerk = (inputArray < tempMean - SDcutoff*tempSD) | (inputArray > tempMean + SDcutoff*tempSD);
    subplot(3,4,12); stem(1:length(p2pResJerkAllStrides),p2pResJerkAllStrides); hold on; stem(find(outlierBoolArrayResjerk),p2pResJerkAllStrides(outlierBoolArrayResjerk),'r'); title('P2P Resultant Jerk (m/sec^3)'); xlabel('Stride Num');
    realMean = mean(p2pResJerkAllStrides(~outlierBoolArrayResjerk));
    realSD = std(p2pResJerkAllStrides(~outlierBoolArrayResjerk));
    realCOV = 100*realSD/realMean;
myStructure.ResJerk_mean = realMean;
myStructure.ResJerk_COV = realCOV;
StrideByStrideStruct.ResJerk = p2pResJerkAllStrides(~outlierBoolArrayResjerk);


if doSave
    tempSubFolder = 'Stride_by_Stride_Figs';
    if ~exist(fullfile(OutputFolder,tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder,tempSubFolder));
    end
    saveas(gcf, fullfile(OutputFolder,tempSubFolder, saveName));
    close
    
    tempSubFolder = 'Outcome_Measures';
    if ~exist(fullfile(OutputFolder,tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder,tempSubFolder));
    end
    save(fullfile(OutputFolder,tempSubFolder, [saveName, '.mat']), 'myStructure');
    
    tempSubFolder = 'Stride_by_Stride_Structs';
    if ~exist(fullfile(OutputFolder,tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder,tempSubFolder));
    end
    save(fullfile(OutputFolder,tempSubFolder, [saveName, '.mat']), 'StrideByStrideStruct');
end

%% grand average  plots of body frame positions and velocities as a function of the gait cycle (avg over all strides)
posERP = nanmean(warpedEpochedData(1:3,:,:),3)';
velERP = nanmean(warpedEpochedData(4:6,:,:),3)'; %nanmean not needed anymore?
jerkERP = nanmean(warpedEpochedData(13:16,:,:),3)';


% plot AP and ML and Up/down position and velcoity
figure;  %sgtitle([SubjStr,' ',saveName]);
subplot(3,4,1);
    plot(tDes,posERP(:,1));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Displacement (m)'); legend({'AP'});%front pos
    grid on;
    APcomExc= range(posERP(:,1)); %Anterior Posterior Center of Mass Excursion (m)
    text(0.05,0.1,['P-P = ', num2str(APcomExc)],'Units','normalized')
subplot(3,4,2);
    plot(tDes,posERP(:,2));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Displacement (m)'); legend({'ML'}); %right pos (in another paper)
    grid on;
    MLcomExc = range(posERP(:,2));
    text(0.05,0.1,['P-P = ', num2str(MLcomExc)],'Units','normalized')
subplot(3,4,3);
    plot(tDes,posERP(:,3));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Displacement (m)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    UDcomExc = range(posERP(:,3));
    text(0.05,0.1,['P-P = ', num2str(UDcomExc)],'Units','normalized')
subplot(3,4,5);
    plot(tDes,velERP(:,1));
    %             ylabel('velocity (m/s)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'AP'});
    grid on;
    PeakAntVel = range(velERP(:,1)); %really p2p a/p vel not peak ant vel
    text(0.05,0.1,['P-P = ', num2str(PeakAntVel)],'Units','normalized')
subplot(3,4,6);
    plot(tDes,velERP(:,2));
    %             ylabel('velocity (m/s)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'ML'});
    grid on;
    PeakLatVel = range(velERP(:,2)); %really p2p m/l velcoity not peak lateral vel
    text(0.05,0.1,['P-P = ', num2str(PeakLatVel)],'Units','normalized')
subplot(3,4,7);
    plot(tDes,velERP(:,3));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    PeakUDvel = range(velERP(:,3)); %really p2p up/down vel
    text(0.05,0.1,['P-P = ', num2str(PeakUDvel)],'Units','normalized')
subplot(3,4,9);
    plot(tDes,jerkERP(:,1));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    PeakAPJerk = range(jerkERP(:,1)); %really p2p up/down vel
    text(0.05,0.1,['P-P = ', num2str(PeakAPJerk)],'Units','normalized')
%jerk
subplot(3,4,10);
    plot(tDes,jerkERP(:,2));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    PeakMLJerk = range(jerkERP(:,2)); %really p2p up/down vel
    text(0.05,0.1,['P-P = ', num2str(PeakMLJerk)],'Units','normalized')    
subplot(3,4,11);
    plot(tDes,jerkERP(:,3));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    PeakUDJerk = range(jerkERP(:,3)); %really p2p up/down vel
    text(0.05,0.1,['P-P = ', num2str(PeakUDJerk)],'Units','normalized')
subplot(3,4,12);
    plot(tDes,jerkERP(:,4));
    %             ylabel('displacement (m)');
    xlabel('% RHS to RHS'); title('Velocity (m/s)'); legend({'UpDown'}); %right pos (in another paper)
    grid on;
    PeakResJerk = range(jerkERP(:,4)); %really p2p up/down vel
    text(0.05,0.1,['P-P = ', num2str(PeakResJerk)],'Units','normalized')

    
if doSave
    tempSubFolder = 'Grand_Avg_Pos_Vel_Plots';
    if ~exist(fullfile(OutputFolder, tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder, tempSubFolder))
    end
    saveas(gcf, fullfile(OutputFolder, tempSubFolder, saveName));
    close
end

%% stride by stride position velocity plots
% GoodStrides = false(nStrides,1); %weird lifehack for H2021 TM_Med_1
% GoodStrides(APcomExcAllStrides < 0.1) = true;

figure;    %sgtitle([SubjStr,' ',saveName]); %3/31/2021 change to make the bad strides appear as dashed lines.
subplot(3,4,1); plot(tDes,squeeze(warpedEpochedData(1,:,~outlierBoolArrayAPexc))); xlim([0 100]); hold on; if any(outlierBoolArrayAPexc); plot(tDes,squeeze(warpedEpochedData(1,:,outlierBoolArrayAPexc)), '--'); end; plot(tDes,posERP(:,1),'k','LineWidth',3); title('AP Displacement (m)'); xlabel('% RHS to RHS');
subplot(3,4,2); plot(tDes,squeeze(warpedEpochedData(2,:,~outlierBoolArrayMLexc))); xlim([0 100]); hold on; if any(outlierBoolArrayMLexc); plot(tDes,squeeze(warpedEpochedData(2,:,outlierBoolArrayMLexc)), '--'); end; plot(tDes,posERP(:,2),'k','LineWidth',3); title('ML Displacement (m)'); xlabel('% RHS to RHS');
subplot(3,4,3); plot(tDes,squeeze(warpedEpochedData(3,:,~outlierBoolArrayUDexc))); xlim([0 100]); hold on; if any(outlierBoolArrayUDexc); plot(tDes,squeeze(warpedEpochedData(3,:,outlierBoolArrayUDexc)), '--'); end; plot(tDes,posERP(:,3),'k','LineWidth',3); title('UpDown Displacement (m)'); xlabel('% RHS to RHS');
subplot(3,4,5); plot(tDes,squeeze(warpedEpochedData(4,:,~outlierBoolArrayAPvel))); xlim([0 100]); hold on; if any(outlierBoolArrayAPvel); plot(tDes,squeeze(warpedEpochedData(4,:,outlierBoolArrayAPvel)), '--'); end; plot(tDes,velERP(:,1),'k','LineWidth',3); title('AP Velocity (m/s)'); xlabel('% RHS to RHS');
subplot(3,4,6); plot(tDes,squeeze(warpedEpochedData(5,:,~outlierBoolArrayMLvel))); xlim([0 100]); hold on; if any(outlierBoolArrayMLvel); plot(tDes,squeeze(warpedEpochedData(5,:,outlierBoolArrayMLvel)), '--'); end; plot(tDes,velERP(:,2),'k','LineWidth',3); title('ML Velocity (m/s)'); xlabel('% RHS to RHS');
subplot(3,4,7); plot(tDes,squeeze(warpedEpochedData(6,:,~outlierBoolArrayUDvel))); xlim([0 100]); hold on; if any(outlierBoolArrayUDvel); plot(tDes,squeeze(warpedEpochedData(6,:,outlierBoolArrayUDvel)), '--'); end; plot(tDes,velERP(:,3),'k','LineWidth',3); title('UpDown Velocity (m/s)'); xlabel('% RHS to RHS');

subplot(3,4,9); plot(tDes,squeeze(warpedEpochedData(13,:,~outlierBoolArrayAPjerk))); xlim([0 100]); hold on; if any(outlierBoolArrayAPjerk); plot(tDes,squeeze(warpedEpochedData(13,:,outlierBoolArrayAPjerk)), '--'); end; plot(tDes,jerkERP(:,1),'k','LineWidth',3); title('AP Jerk (m/sec^3)'); xlabel('% RHS to RHS');
subplot(3,4,10); plot(tDes,squeeze(warpedEpochedData(14,:,~outlierBoolArrayMLjerk))); xlim([0 100]); hold on; if any(outlierBoolArrayMLjerk); plot(tDes,squeeze(warpedEpochedData(14,:,outlierBoolArrayMLjerk)), '--'); end; plot(tDes,jerkERP(:,2),'k','LineWidth',3); title('ML Jerk (m/sec^3)'); xlabel('% RHS to RHS');
subplot(3,4,11); plot(tDes,squeeze(warpedEpochedData(15,:,~outlierBoolArrayUDjerk))); xlim([0 100]); hold on; if any(outlierBoolArrayUDjerk); plot(tDes,squeeze(warpedEpochedData(15,:,outlierBoolArrayUDjerk)), '--'); end; plot(tDes,jerkERP(:,3),'k','LineWidth',3); title('UpDown Jerk (m/sec^3)'); xlabel('% RHS to RHS');
subplot(3,4,12); plot(tDes,squeeze(warpedEpochedData(16,:,~outlierBoolArrayResjerk))); xlim([0 100]); hold on; if any(outlierBoolArrayResjerk); plot(tDes,squeeze(warpedEpochedData(16,:,outlierBoolArrayResjerk)), '--'); end; plot(tDes,jerkERP(:,4),'k','LineWidth',3); title('Resultant Jerk (m/sec^3)'); xlabel('% RHS to RHS');




%spaghetti plot movie ONLY DO FOR TESTING SAYS COLTON (otherwise we overwrite the figure made above when saving below... ): )
% figure; 
% for tempstride_i = 1:nStrides
%     plot(tDes,squeeze(warpedEpochedData(1,:,tempstride_i))); hold on;
%     pause(.1);
% end

% figure;    %sgtitle([SubjStr,' ',saveName]); 
% subplot(2,3,1); plot(tDes,squeeze(warpedEpochedData(1,:,:))); xlim([0 100]); hold on; plot(tDes,posERP(:,1),'k','LineWidth',3); title('AP Displacement (m)'); xlabel('% RHS to RHS');
% subplot(2,3,2); plot(tDes,squeeze(warpedEpochedData(2,:,:))); xlim([0 100]); hold on; plot(tDes,posERP(:,2),'k','LineWidth',3); title('ML Displacement (m)'); xlabel('% RHS to RHS');
% subplot(2,3,3); plot(tDes,squeeze(warpedEpochedData(3,:,:))); xlim([0 100]); hold on; plot(tDes,posERP(:,3),'k','LineWidth',3); title('UpDown Displacement (m)'); xlabel('% RHS to RHS');
% subplot(2,3,4); plot(tDes,squeeze(warpedEpochedData(4,:,:))); xlim([0 100]); hold on; plot(tDes,velERP(:,1),'k','LineWidth',3); title('AP Velocity (m/s)'); xlabel('% RHS to RHS');
% subplot(2,3,5); plot(tDes,squeeze(warpedEpochedData(5,:,:))); xlim([0 100]); hold on; plot(tDes,velERP(:,2),'k','LineWidth',3); title('ML Velocity (m/s)'); xlabel('% RHS to RHS');
% subplot(2,3,6); plot(tDes,squeeze(warpedEpochedData(6,:,:))); xlim([0 100]); hold on; plot(tDes,velERP(:,3),'k','LineWidth',3); title('UpDown Velocity (m/s)'); xlabel('% RHS to RHS');

if doSave
    tempSubFolder = 'Stride_by_Stride_Pos_Vel_Plots';
    if ~exist(fullfile(OutputFolder, tempSubFolder),'dir')
        mkdir(fullfile(OutputFolder, tempSubFolder));
    end
%     saveas(gcf, 'R:\Ferris-Lab\share\MindInMotion\Scripts\testFigure2')
    saveas(gcf, fullfile(OutputFolder, tempSubFolder, saveName));
    close
end

end

