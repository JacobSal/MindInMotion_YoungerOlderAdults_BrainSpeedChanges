function [STUDY,spec_savef,spec_subjcorr_savef] = mim_gen_spec_data(STUDY,ALLEEG,...
    des_i,cluster_i,design_char,save_dir,varargin)
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic1 = tic;
%## DEFINE DEFAULTS
STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000); 
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'comps','all');
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'des_i',@isnumeric);
addRequired(p,'cluster_i',@isnumeric);
addRequired(p,'design_char',@ischar);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'SPEC_PARAMS',SPEC_PARAMS,@isstruct);
addParameter(p,'STAT_PARAMS',STAT_PARAMS,@isstruct);
parse(p,STUDY,ALLEEG,des_i,cluster_i,...
    design_char,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
SPEC_PARAMS = p.Results.SPEC_PARAMS;
%% ===================================================================== %%
%## SPECTRUM CALCULATIONS
spec_savef = '';
spec_subjcorr_savef = '';
%##
fprintf('Computing specdata for cluster %i...\n',cluster_i)
tic
try
    [~,specdata,specfreqs,pgroup,pcond,pinter] = std_specplot(STUDY, ALLEEG, 'clusters',cluster_i,...
        'comps',SPEC_PARAMS.comps,'subject',SPEC_PARAMS.subject,'freqrange',SPEC_PARAMS.freqrange,'design',des_i);
    %- save dat
    spec_data = struct('specdata',specdata,'specfreqs',specfreqs,...
        'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
    par_save(spec_data,save_dir,sprintf('spec_data_cl%i_%s.mat',cluster_i,design_char));
    spec_savef = [save_dir filesep sprintf('spec_data_cl%i_%s.mat',cluster_i,design_char)];
    %- save fig
    fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
    saveas(fig_i,[save_dir filesep sprintf('spec_plot_cl%i_%s_allcomps',cluster_i,design_char)])
catch e
    fprintf(['error. identifier: %s\n',...
         'error. %s\n',...
         'stack. %s\n'],e.identifier,e.message,getReport(e));
end
toc
%##
fprintf('Computing specdata with substracted subject mean for cluster %i...\n',cluster_i)
tic
try
    STUDY.etc.specparams.subtractsubjectmean = 'on';
    [~,specdata,specfreqs,pgroup,pcond,pinter] = std_specplot(STUDY, ALLEEG, 'clusters',cluster_i,...
        'comps',SPEC_PARAMS.comps,'subject',SPEC_PARAMS.subject,'freqrange',SPEC_PARAMS.freqrange,...
        'subtractsubjectmean','on','design',des_i);
    %- save dat
    spec_data = struct('specdata',specdata,'specfreqs',specfreqs,...
        'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
    par_save(spec_data,save_dir,sprintf('spec_data_cl%i_%s_subtractmean.mat',cluster_i,design_char));
    spec_subjcorr_savef = [save_dir filesep sprintf('spec_data_cl%i_%s_subtractmean.mat',cluster_i,design_char)];
    %- save fig
    fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
    saveas(fig_i,[save_dir filesep sprintf('spec_plot_cl%i_%s_allcomps_subtractsubjectmean',cluster_i,design_char)])
catch e
    fprintf(['error. identifier: %s\n',...
         'error. %s\n',...
         'stack. %s\n'],e.identifier,e.message,getReport(e));
end
toc 
fprintf('done.\n')
close all
%## time
toc(tic1)
end

