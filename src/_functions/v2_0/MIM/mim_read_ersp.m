function [STUDY,data_ersp,xvals,yvals,events_ersp,params] = mim_read_ersp(STUDY,cluster_number,varargin)
%MIM_REJECT_ICS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Jacob Salminen
% Code Date: 06/09/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, 
%## TIME
tic
%## DEFINE DEFAULTS
%- find eeglab on path
% tmp = strsplit(path,';');
% b1 = regexp(tmp,'eeglab','end');
% b2 = tmp(~cellfun(@isempty,b1));
% PATH_EEGLAB = b2{1}(1:b1{1});
% fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- fileextension
FILE_EXT = '.icatimef';
%- designnumber
designnumber = STUDY.currentdesign;
%- timelimits
timelimits = STUDY.etc.erspparams.timerange;
%- freqlimits
freqlimits = STUDY.etc.erspparams.freqrange;
%- singletrials
singletrials = 'off'; %STUDY.etc.erspparams.singletrials;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'cluster_number',@isnumeric);
%## OPTIONAL
addOptional(p,'designnumber',designnumber,@isnumeric);
addOptional(p,'timelimits',timelimits,@isnumeric);
addOptional(p,'freqlimits',freqlimits,@isnumeric);
addOptional(p,'singletrials',singletrials,@ischar);
%## PARAMETER
parse(p,STUDY,cluster_number,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
designnumber = p.Results.designnumber;
timelimits = p.Results.timelimits;
freqlimits = p.Results.freqlimits;
singletrials = p.Results.singletrials;
%- PARAMETER
%% ===================================================================== %%
%## ASSIGN PARAMETERS FOR READING
%-
opts = {'timelimits',timelimits,'freqlimits',freqlimits,'singletrials',singletrials};
%- 
design_var = struct(STUDY.design(designnumber).variable);
%- assign session
allSessions = { STUDY.datasetinfo.session };
allSessions(cellfun(@isempty, allSessions)) = {1};
allSessions = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);
%- prepare loop
data_ersp = cell(1,length(STUDY.datasetinfo));
params = cell(1,length(STUDY.datasetinfo));
xvals = cell(1,length(STUDY.datasetinfo));
yvals = cell(1,length(STUDY.datasetinfo));
events_ersp = cell(1,length(STUDY.datasetinfo));
for subj_i = 1:length(STUDY.datasetinfo)
    fprintf('Loading subject %s...\n',STUDY.datasetinfo(subj_i).subject);
    sub_comps = (subj_i == STUDY.cluster(cluster_number).sets); 
    comp_list = STUDY.cluster(cluster_number).comps(sub_comps);
    fileName = getfilename({STUDY.datasetinfo(subj_i).filepath}, STUDY.datasetinfo(subj_i).subject, { STUDY.datasetinfo(subj_i).session }, FILE_EXT, length(uniqueSessions) == 1);
    [data_ersp{subj_i}, params{subj_i}, xvals{subj_i}, yvals{subj_i}, events_ersp{subj_i} ] = std_readfile( fileName, 'designvar', design_var, opts{:}, 'components', comp_list);
end
%## TIME
toc

end

% get file base name: filepath and sess are cell array (in case 2 files per subject)
% ----------------------------------------------------------------------------------
function filebase = getfilename(filepath, subj, sess, fileSuffix, onlyOneSession)
    if onlyOneSession
        filebase = fullfile(filepath{1}, [ subj fileSuffix ] );
    else
        if isempty(sess)
            sess = { '1' };
        end
        for iSess = 1:length(sess)
            if isnumeric(sess{iSess})
                sesStr   = [ '0' num2str(sess{iSess}) ];
            else
                sesStr   = [ '0' sess{iSess} ];
            end
            filebase{iSess} = fullfile(filepath{iSess}, [ subj '_ses-' sesStr(end-1:end) fileSuffix ] );
        end
        if length(unique(filebase)) < length(filebase)
            filebase = unique(filebase); % order is not important
        end
    end
end
