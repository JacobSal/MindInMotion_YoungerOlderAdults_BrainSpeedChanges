%   Project Title: Run a graph analysis for multiple subjects
%
%   Code Designer: Jacob salminen
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/2_GLOBAL_BATCH/alpha/run_GLOBAL.sh
%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% Initialization
%## TIME
tic
%% REQUIRED SETUP 4 ALL SCRIPTS
%- DATE TIME
dt = datetime;
dt.Format = 'ddMMyyyy';
%- VARS
USER_NAME = 'jsalminen'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% PATH TO YOUR GITHUB REPO
%- GLOBAL VARS
REPO_NAME = 'par_EEGProcessing';
%- determine OS
if strncmp(computer,'PC',2)
    DO_UNIX = false;
    PATH_EXT = 'M';
else  % isunix
    DO_UNIX = true;
    PATH_EXT = 'dferris';
end
%## DEBUG: PATHROOT OVERRIDE
if DO_UNIX
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
else
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'GitHub']; % path 2 your github folder
end
%% SETWORKSPACE
%- define the directory to the src folder
source_dir = [PATH_ROOT filesep REPO_NAME filesep 'src'];
run_dir = [source_dir filesep '_test' filesep '1_paper_MIM'];
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
%- set workspace
setWorkspace
%% PARPOOL SETUP
if DO_UNIX
%     eeg_options;
    pop_editoptions('option_parallel',1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    % pp.NumWorkers = POOL_SIZE-3;
    % pp.NumThreads = 1;
    fprintf('Number of workers: %i',pp.NumWorkers);
    fprintf('Number of threads: %i',pp.NumThreads);
    %- make meta data directory for slurm
    mkdir([BATCH_DIR filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat(run_dir, getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, POOL_SIZE, 'IdleTimeout', 1440);
end
%% ================================================================= %%
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
    'H1022','H1024','H1026','H1027','H1033','H1034'};
SUBJ_HMA = {'H2002', 'H2010', 'H2015', 'H2017', 'H2020', 'H2021', 'H2022', 'H2023',...
    'H2025', 'H2026', 'H2034', 'H2059', 'H2062', 'H2082', 'H2095'};
SUBJ_NMA = {'NH3008', 'NH3043', 'NH3055', 'NH3059', 'NH3069', ...
    'NH3070', 'NH3074', 'NH3086', 'NH3090', 'NH3104', 'NH3105', 'NH3106', 'NH3112', 'NH3114'};
DATASET_PROCESSNUM = 2;
TRIAL_TYPES = {'rest','0p5','0p25','0p75', '1p0','flat','low','med','high'};
%- Subject Picks
% SUBJ_PICS = {SUBJ_YNG,SUBJ_HMA,SUBJ_NMA};
SUBJ_PICS = {SUBJ_YNG};
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
%% ===================================================================== %%
%## PROCESSING PARAMS
%- subjstruct and saving
LOAD_STUDY = true;
%## POST-PROCESSING PARAMS
%- datetime override
dt = '26122022';
%- hard define
conn_meas = 'dDTF';
% load_trials = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
% load_trials = {'rest','0p25','0p5','0p75','1p0'}; %,'flat','low','med','high'};
% load_trials = {'flat','low','med','high'};
load_trials = {'rest','0p25','0p5','0p75','1p0'};
%- soft define
subjinfDir = [SUBJINF_DIR filesep sprintf('%s',dt)];
save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep '_figs'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%- create new subject directory
if ~exist(subjinfDir,'dir')
    mkdir(subjinfDir);
end
%% ==== Load Study ==== %%
STUDIES = cell(1,length(load_trials));
ALLEEGS = cell(1,length(load_trials));
%- Create STUDY & ALLEEG structs
for cond_i = 1:length(load_trials)
    study_fName = sprintf('%s_MIM_study',load_trials{cond_i});
    fprintf('');
    if ~exist([load_dir filesep study_fName '.study'],'file')
        error('ERROR. study file does not exist');
    else
        if DO_UNIX
            [STUDIES{cond_i},ALLEEGS{cond_i}] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',load_dir);
        else
            [STUDIES{cond_i},ALLEEGS{cond_i}] = pop_loadstudy('filename',[study_fName '.study'],'filepath',load_dir);
        end
    end
end
%- 
MAIN_STUDY = STUDIES{1};
MAIN_ALLEEG = ALLEEGS{1};
%- extract component array
comps_store = zeros(length(MAIN_ALLEEG),length(MAIN_STUDY.cluster));
for clus_i = 2:length(MAIN_STUDY.cluster)
    sets_i = MAIN_STUDY.cluster(clus_i).sets;
    for j = 1:length(sets_i)
        comps_store(sets_i(j),clus_i) = MAIN_STUDY.cluster(clus_i).comps(j);
    end
end
%% Load some data
subj_i = 5;
cond_i = 1;
%- loop
EEG = ALLEEGS{cond_i}(subj_i);
% EEG = cnctanl_loadCAT(EEG,SUBJSTRUCT,iterTrial,subj_i,'BootStrap');
EEG = cnctanl_loadCAT(EEG,'NonzeroTest');
%## ASSIGN TEMPORARY EEG STRUCTURE
tmpEEG = EEG;
%% VISUALIZE TIME FREQ GRID (components to component)
%## TIMEFREQ AVG PLOT
figs_save_dir = [save_dir filesep 'indvidual_TxF_conn' filesep load_trials{cond_i} filesep sprintf('%i',subj_i)];
if ~exist(figs_save_dir,'dir')
    mkdir(figs_save_dir)
end
connMethods = {'dDTF'}; %{'dDTF','GGC','ffDTF'};
tmpcl = squeeze(comps_store(subj_i,:));
[tmpcl,idxcl] = sort(tmpcl);
idxcl = idxcl(tmpcl ~= 0);
tmpcl = tmpcl(tmpcl ~= 0);
clusts = cell(1,length(idxcl));
comps = zeros(1,length(idxcl));
for i = 1:length(clusts)
    clusts{i} = sprintf('%i',idxcl(i));
    comps(i) = tmpcl(i);
end
ComponentNames = EEG.CAT.curComponentNames;
%- plot
% fprintf('%s) Plotting Time Frequency Plot...\n',SUBJSTRUCT(subj_i).subjStr)
cnctanl_visualize(tmpEEG,PATHS,comps,clusts,...
    'IndvCell_TxF',...
    'savePath',figs_save_dir,...
    'connMeasures',connMethods,...
    'CLIM',[0.0,0.002]);
%%
% est_fitMVAR(EEG)
% OUTPUT:
%  AR    multivariate autoregressive model parameter
%  RC    reflection coefficients (= -PARCOR coefficients)
%  PE    remaining error variance
%  [AR,RC,PE] = mvar(Data,g.morder)
%  e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR,1)),Data);
%-
% tmpSrcdata = EEG.CAT.srcdata;
% [AR, PE, argsout] = mvar_arfit(tmpSrcdata,modelOrder,'sbc',true);

%%
%- Note:
% -2 = lower epoch time limit in (s), (xmin)
% 2 = upper epoch time limit in (s), (xmax)
% 500 = sampling rate in (hz), (EEG.srate)
% 0.5 = window length for the model (s), (MODEL.winlen)
% 0.0240 = window step size in (s),(MODEL.winstep)
% 146 = numWins
% 7 = model order
% 14 = number of components being compared (multiple variables)
% then 7*14 = number of coefficient that will exist for every component in
% window t;
% (i.e., size(AR{t}) = (14x98) && size(AR) = 146 && numWins = 146)
% t = 1;
% see. https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4071888
% see. https://www3.stat.sinica.edu.tw/statistica/oldpdf/A15n112.pdf
% see. Conn = est_mvtransfer(varargin)
% see. [Model,cfg] = est_fitMVAR(varargin)
% see. [ARF,PE,argsout] = mvar_vieiramorf(varargin)
g = [];
%- model
AR = EEG.CAT.MODEL.AR;
PE = EEG.CAT.MODEL.PE;
%- params
epochN = size(EEG.CAT.srcdata,3);
timeN = size(EEG.CAT.srcdata,2);
tmpL = size(EEG.CAT.srcdata,2);
tmpXmin = 0;
tmpXmax = epochN*(timeN/EEG.srate);
tmpTimes = (tmpXmin:(1/EEG.srate):tmpXmax)*1000;
total_time = [tmpXmin tmpXmax];
epochDur = EEG.xmax-EEG.xmin;
compN = size(EEG.CAT.srcdata,1);
%- window & time
mod_order = EEG.CAT.MODEL.morder;
winLenPnts  = round(EEG.CAT.MODEL.winlen*EEG.srate);
winStepPnts = round(EEG.CAT.MODEL.winstep*EEG.srate);
numWins   = length(EEG.CAT.MODEL.winStartTimes);
%- initialize loop params
y = zeros(compN,length(EEG.CAT.MODEL.winStartTimes),size(EEG.CAT.srcdata,3));
ptrue = zeros(compN,length(EEG.CAT.MODEL.winStartTimes),size(EEG.CAT.srcdata,3));
ytrue = zeros(compN,length(EEG.CAT.MODEL.winStartTimes),size(EEG.CAT.srcdata,3));
timeConn = zeros(1,numWins);
timeDat = zeros(1,numWins);
timeX1 = zeros(length(EEG.CAT.MODEL.winStartTimes),size(EEG.CAT.srcdata,3));
%-
nchs = compN;
z = 2*pi*1i/EEG.srate;
I = eye(nchs);
% sqrtinvcov = diag(diag(PEt).^(-1/2));
offset = 0; %-ceil(1);
%- Connectivity Play Params
ConnEst = cell(1,numWins);

%- loop over epochs in data
for epoch_i = 1:size(EEG.CAT.srcdata,3)    
    %- loop over windows in epoch
    for t = 1:numWins
        %- extract model coefficients
        ARt = reshape(AR{t},nchs,nchs,mod_order);
        %- extract noise covariance matrix
        if size(PE{t},2)>nchs
            PEt = PE{t}(:,nchs*mod_order+1:nchs*(mod_order+1));
        end
%         PEt = reshape(PE{t},nchs,nchs,mod_order+1);
%         sqrtinvcov = diag(diag(PEt).^(-1/2));
        %- apply coefficients to window
        mod_start = ceil((EEG.CAT.MODEL.winStartTimes(t)+EEG.CAT.MODEL.winlen/2)*EEG.srate);
        %- 
%         winpnts = (mod_start-mod_order+1):(mod_start);
        winpnts = (mod_start-mod_order-1):(mod_start-2);
        tmpDat = squeeze(EEG.CAT.srcdata(:,winpnts,epoch_i));
%         tmpDat = squeeze(EEG.CAT.srcdata(:,winpnts,epoch_i));
        ut = diag(PEt);
%         ut = diag(sqrtinvcov);
        for k = 1:mod_order            
            ARk = squeeze(ARt(:,:,k));
            tmpARD = ARk*tmpDat(:,k);
            y(:,t,epoch_i) = (y(:,t,epoch_i) + tmpARD); %+ PEt(:,:,k)*tmpDat(:,k); % + PEt*tmpDat(:,k); %+ ut;
%             ptrue(:,t,epoch_i) = ptrue(:,t,epoch_i) + ut; %- PEt(:,:,k)*tmpDat(:,k); %
        end
        y(:,t,epoch_i) = y(:,t,epoch_i); % + ut;
        ytrue(:,t,epoch_i) = EEG.CAT.srcdata(:,mod_start,epoch_i);
        ptrue(:,t,epoch_i) = PEt*tmpDat(:,k); %ut;
        timeConn(t) = mod_start/EEG.srate - EEG.CAT.MODEL.winstep; % - (1/EEG.srate);
        timeDat(t) = mod_start/EEG.srate - EEG.CAT.MODEL.winstep + (1/EEG.srate);
        %- some ASS logic
        if epoch_i == 1
            timeX1(t,epoch_i) = mod_start/EEG.srate;
            epochConnEnd = mod_start/EEG.srate;
        else
            if t == 1
                % time for next epoch = last epoch time + time lag ; where
                % time lag = (total epoch time - conn last window epoch
                % time + conn first window epoch time)
                timeX1(t,epoch_i) = timeX1(end,epoch_i-1) + (epochDur-epochConnEnd) + mod_start/EEG.srate;
            else
                timeX1(t,epoch_i) =  timeX1(t-1,epoch_i) + EEG.CAT.MODEL.winstep;
            end
        end
        %{
        ConnEst{t} = est_mvtransfer('AR',ARt,...
                                    'C',PEt,...
                                    'freqs',(1:60),...
                                    'srate',EEG.srate,...
                                    'connmethods',{'dDTF'},...
                                    'arg_direct',true);
        %}
        
    end
end
%% Sum of Squares & covariance calculation
raw_t = EEG.times/1000+2;
diff_yx = zeros(size(y));
cov_yx = zeros(size(y));
for epoch_i = 1:size(EEG.CAT.srcdata,3)
    for comp_i = 1:size(y,1)
        for t = 2:numWins
            %- create line
            coefficients = polyfit([timeConn(t-1), timeConn(t)], [y(comp_i,t-1,epoch_i), y(comp_i,t,epoch_i)], 1);
            a = coefficients(1);
            b = coefficients(2);
            tt = raw_t(round(timeConn(t-1)*EEG.srate):round(timeConn(t)*EEG.srate));
            x = tt*a+b;
            %- extract MVAR
            raw_dat = squeeze(EEG.CAT.srcdata(comp_i,round(timeDat(t-1)*EEG.srate):round(timeDat(t)*EEG.srate),epoch_i));
            diff_yx(comp_i,t-1,epoch_i) = sum(sqrt((raw_dat - x).^2)); % root mean squared
            cov_yx(comp_i,t-1,epoch_i) = sum((x-mean(x)).*(raw_dat-mean(raw_dat)))/(length(raw_dat)-1); % covariance 
            %{
            figure;
            hold on;
            plot(tt,x);
            plot(tt,raw_dat);
            hold off;
            %}
        end
    end
end

%% PLOTS?
%- compare mvar to a particular epoch
comp_i = 1;
epoch_i = 11;
figure;
%- subplot 1
subplot(3,1,1);
% h = [];
hold on;
h = plot(timeConn,squeeze(y(comp_i,:,epoch_i)),'DisplayName',sprintf('MVAR comp %i',comp_i));
h.Color = color.Blue;
h = plot(EEG.times(round(timeConn(1)*EEG.srate):round(timeConn(end)*EEG.srate))/1000+2,squeeze(EEG.CAT.srcdata(comp_i,round(timeConn(1)*EEG.srate):round(timeConn(end)*EEG.srate),epoch_i)),'-','DisplayName',sprintf('RAW comp %i',comp_i));
h.Color = color.lightRed;
hold off;
title('MVAR Estimate @ t vs. RAW EEG for all times');
ylim([-5,5])
ylabel('Voltage (uV)')
xlabel('Epoch Time (s)');
legend();
%- subplot 2
% h = [];
subplot(3,1,2);
hold on;
h = plot(timeConn,squeeze(y(comp_i,:,epoch_i)),'DisplayName',sprintf('MVAR comp %i',comp_i));
h.Color = color.Blue;
h = plot(timeDat,ytrue(comp_i,:,epoch_i)','-','DisplayName',sprintf('y(t) comp %i',comp_i));
h.Color = color.lightRed;
hold off;
title('MVAR Estimate @ t vs. RAW EEG @ t');
ylim([-5,5])
ylabel('Voltage (uV)')
xlabel('Epoch Time (s)');
legend();
%- subplot 3
% h = [];
subplot(3,1,3);
%-
hold on;
disp_name = sprintf('sqrt(RAW(comp %i)-MVAR(comp %i))^2',comp_i,comp_i);
% disp_name = sprintf('sum((x-mean(x))*(y-mean(y))/N-1');
% tmp = cov_yx(comp_i,:,epoch_i); 
tmp = sqrt(ytrue(comp_i,:,epoch_i)-squeeze(y(comp_i,:,epoch_i))).^2;
h = plot(timeConn,tmp,'DisplayName',disp_name);
h.Color = color.Blue;
h = plot(timeDat,squeeze(ptrue(comp_i,:,epoch_i)),'DisplayName',sprintf('err coeffs'));
h.Color = color.lightRed;
hold off;
title('sqrt( MVAR - Raw @ t)^2 vs. Covariance Matrix Coefficients @ t');
ylim([-5,5])
ylabel('Voltage (uV)')
xlabel('Epoch Time (s)');
legend();
%% CONTINUOUS SCROLL
%## (RESTING STATE) TIME SCROLL PLOT
connMethods = {'dDTF','GGC'};
freqs = EEG.CAT.Conn.freqs;
comp_i = 1;
j = 2;
alpha = 0.05;
%- connectivity 
connScroll = squeeze(y(comp_i,:,:))';
connScroll = squeeze(connScroll(:));
% timeX1 = timeConn(1):EEG.CAT.MODEL.winstep:(EEG.CAT.MODEL.winstep*size(y,2)*size(y,3)+EEG.CAT.MODEL.winstep*size(y,1));
timeX1 = timeX1(:);

%- epochs
% tmpEpochsi = squeeze(EEG.icaact(comp_i,:,:))';
% tmpEpochsi = tmpEpochsi(:);
% timeN = size(EEG.CAT.srcdata,2);
% tmpXmin = 0;
% tmpXmax = epochN*(timeN/EEG.srate);
% timeX2 = (tmpXmin:(1/EEG.srate):(tmpXmax-(1/EEG.srate)));
compScroll = squeeze(ytrue(comp_i,:,:))';
compScroll = squeeze(compScroll(:));


%- plot
figure();
hold on;
% h = plot(timeX2,tmpEpochsi,'DisplayName',sprintf('RAW comp %i',comp_i));
% h.Color = color.Blue;
h = plot(timeX1,compScroll,'DisplayName',sprintf('RAW comp %i',comp_i));
h.Color = color.Blue;
h = plot(timeX1,connScroll,'DisplayName',sprintf('MVAR comp %i',comp_i));
h.Color = color.lightRed;
hold off;
legend()
%% Connectivity Play
%- (01/08/2022), JS: Doesnt work right now
%{
ARest = cell(1,numWins);
PEest = cell(1,numWins);
winStartIdx = [];
epochTimeLims = [EEG.xmin,EEG.xmax];
tidx = getindex(EEG.CAT.times,epochTimeLims*1000);
winLenPnts  = round(EEG.CAT.MODEL.winlen*EEG.srate); % window size in points
winStepPnts = round(EEG.CAT.MODEL.winstep*EEG.srate);
if isempty(winStartIdx)
    winStartIdx  = tidx(1):winStepPnts:(tidx(2)-winLenPnts)+1;
end

for t = 1:numWins
    mod_start = ceil((EEG.CAT.MODEL.winStartTimes(t)+EEG.CAT.MODEL.winlen/2)*EEG.srate);
%     winpnts = (mod_start-mod_order):(mod_start-1);
    winpnts = winStartIdx(t):winStartIdx(t)+winLenPnts-1;
    [ARest{t},PEest{t}] = mvar_vieiramorf('data',EEG.CAT.srcdata(:,winpnts,:));
end
%}
%% Version History
%{
v1.0.01132023.0 : Initializing versioning for future iterations. Previous
Params:

%}