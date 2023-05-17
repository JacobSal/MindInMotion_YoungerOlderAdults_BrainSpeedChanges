function [STUDY,ALLEEG] = eeglab_cluster(STUDY,ALLEEG,varargin)
%SUBJ_I_CLUSTER Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
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
DO_RECOMPUTE = 'off';
%- spec params
SPEC_MODE = 'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
LOG_TRIALS = 'off'; %options: 'on'
FREQ_FAC = 4;
SPEC_CONTINUOUS = 'off';
PAD_RATIO = 2;
SPEC_EPOCH_LIMS = [0,4];
ICLABEL_BRAIN_CRIT = 0.5;
RV_BRAIN_CRIT = 0.15;
CLUS_NUM = 10;
%* 
FREQ_LIMITS = [1,100];
errorMsg = 'Value must be of format [INT1,INT2] where INT2 > INT1. Sets value limits for calculated spectral measure'; 
fl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%* 
CYCLE_LIMITS = [3,0.8];
errorMsg = 'Value must be of format [DOUBLE1,DOUBLE2]. Varies the length of the window used to determine the power of a certain frequency band';
cl_validFcn = @(x) assert((isnumeric(x) && length(x) == 2) || isempty(x),errorMsg);
%*
DO_ERSP_CALC = false;
%- cluster params
SET_CLUSTER_N = 10; %8;
DO_VISUALIZATION = false;
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'FREQ_LIMITS',FREQ_LIMITS,fl_validFcn); 
addParameter(p,'CYCLE_LIMITS',CYCLE_LIMITS,cl_validFcn); 
addParameter(p,'SPEC_MODE',SPEC_MODE,@ischar); 
addParameter(p,'LOG_TRIALS',LOG_TRIALS,@ischar); 
addParameter(p,'FREQ_FAC',FREQ_FAC,@isnumeric); 
addParameter(p,'SPEC_CONTINUOUS',SPEC_CONTINUOUS,@ischar); 
addParameter(p,'PAD_RATIO',PAD_RATIO,@isnumeric); 
addParameter(p,'SPEC_EPOCH_LIMS',SPEC_EPOCH_LIMS,@isnumeric); 
addParameter(p,'SET_CLUSTER_N',SET_CLUSTER_N,@isnumeric); 
addParameter(p,'DO_VISUALIZATION',DO_VISUALIZATION,@islogical); 
addParameter(p,'DO_ERSP_CALC',DO_ERSP_CALC,@islogical);
addParameter(p,'ICLABEL_BRAIN_CRIT',ICLABEL_BRAIN_CRIT,@isnumeric); 
addParameter(p,'RV_BRAIN_CRIT',RV_BRAIN_CRIT,@isnumeric);
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
FREQ_LIMITS         = p.Results.FREQ_LIMITS;
CYCLE_LIMITS        = p.Results.CYCLE_LIMITS; %#ok<NASGU>
SPEC_MODE           = p.Results.SPEC_MODE;
LOG_TRIALS          = p.Results.LOG_TRIALS; %#ok<NASGU>
FREQ_FAC            = p.Results.FREQ_FAC;
SPEC_CONTINUOUS     = p.Results.SPEC_CONTINUOUS; %#ok<NASGU>
PAD_RATIO           = p.Results.PAD_RATIO; %#ok<NASGU>
SET_CLUSTER_N       = p.Results.SET_CLUSTER_N; %#ok<NASGU>
SPEC_EPOCH_LIMS     = p.Results.SPEC_EPOCH_LIMS; %#ok<NASGU>
%% ===================================================================== %%
%## SAVE STUDY
for subj_i = 1:length(ALLEEG)
    %- check iclabel 
    if ~isfield(ALLEEG(subj_i).etc,'ic_classification')
       error('For some reason ICLABELS were not calculated for subject %s',ALLEEG(subj_i).subject)
    end
    %- check dipfit
    if ~isfield(ALLEEG(subj_i).dipfit,'model')
        disp(ALLEEG(subj_i).dipfit)
        error('For some reason DIPFITS were not calculated for subject %s',ALLEEG(subj_i).subject)
    end
end
%%
% make sure ALLEEG & STUDY structs are consistent
[STUDY, ALLEEG] = std_checkset(STUDY,ALLEEG);
% remove dipoles from analysis.
fprintf('\n');
% this computes power spectral density for designated study design
fprintf('==== Performing Precomputing ====\n');
if strcmp(SPEC_MODE,'psd')
%     SPEC_EPOCHLIM = [-1,1];
%     BOUNDARIES = (1:length(SPEC_EPOCHLIM):ALLEEG(1).xmax);
    BOUNDARIES = [];
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                    'components',... 
                                    'allcomps','on',...
                                    'recompute',DO_RECOMPUTE,...
                                    'spec','on',...
                                    'scalp','on',...
                                    'specparams',...
                                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                                    'freqrange',FREQ_LIMITS,'boundaries',BOUNDARIES});
else
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                    'components',...                               
                                    'recompute',DO_RECOMPUTE,...
                                    'spec','on',...
                                    'scalp','on',...
                                    'specparams',...
                                    {'specmode',SPEC_MODE,'freqfac',FREQ_FAC,...
                                    'freqrange',FREQ_LIMITS});
end
if DO_ERSP_CALC
    [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
                                'components',...
                                'ersp','on',...
                                'erspparams',...
                                {'cycles' CYCLE_LIMITS,'freqs',FREQ_LIMITS,...
                                'padratio',PAD_RATIO});
end
fprintf('done.\n');
%% PRECLUSTERING
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,...
%                                 'components',...
%                                 'rmicacomps','on');
% make sure ALLEEG & STUDY structs are consistent
[STUDY, ALLEEG] = std_checkset(STUDY,ALLEEG);
fprintf('\n');
% CLUST METHOD (MIM 04/23/2023)
[STUDY, ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
                        { 'dipoles', 'norm' 1, 'weight' 5 },...
                        { 'scalp', 'norm' 1, 'weight' 5 });

% this uses the designated algorithm to cluster components
% Try: sweeping from 5 to 15 cluster numbers and see how they compare
% Try: using average number of components for all subjects.
fprintf('\n');
[STUDY]         = pop_clust(STUDY,ALLEEG,...
                    'algorithm','kmeans',...
                    'clus_num',CLUS_NUM); %,'outliers',1);

%## ADD ANATOMICAL LABELS
% [~, atlas_cell] = add_anatomical_labels(STUDY,ALLEEG);
% STUDY.etc.add_anatomical_labels = atlas_cell;
% STUDY = std_movecomp(STUDY,ALLEEG, from_cluster, to_cluster, comps);
%## KEEP TRACK OF REJECTED COMPS
%% SAVE STUDY && CONVERT ALLEEG PATHS TO UNIX
STUDY.urcluster = STUDY.cluster;
% check the STUDY & ALLEEG for consistency, and save.
[STUDY,ALLEEG] = std_checkset(STUDY,ALLEEG);
[STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                        STUDY.filename,STUDY.filepath,...
                                        'RESAVE_DATASETS','on');
end
