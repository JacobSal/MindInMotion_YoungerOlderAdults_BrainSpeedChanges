function [STUDY,ersp_savef,ersp_subbase_savef,ersp_subbase_combase_savef,ersp_singletrial_subbase_savef] = mim_gen_ersp_data(STUDY,ALLEEG,warping_times,...
    ersp_load_params,des_i,cluster_i,design_char,save_dir,varargin)
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
ERSP_PARAMS = struct('subbaseline','on',...
    'timerange',[warping_times(1) warping_times(end)],...
    'ersplim',[-1.5,1.5],...
    'freqrange',[1,200]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'warping_times',@isnumeric);
addRequired(p,'ersp_load_params',@isstruct);
addRequired(p,'des_i',@isnumeric);
addRequired(p,'cluster_i',@isnumeric);
addRequired(p,'design_char',@ischar);
addRequired(p,'save_dir',@ischar);
%## OPTIONAL
%## PARAMETER
addParameter(p,'ERSP_PARAMS',ERSP_PARAMS,@isstruct);
addParameter(p,'STAT_PARAMS',STAT_PARAMS,@isstruct);
parse(p,STUDY,ALLEEG,warping_times,ersp_load_params,des_i,cluster_i,design_char,save_dir,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
%## MAKE DIRS
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
ERSP_PARAMS = p.Results.ERSP_PARAMS;
STAT_PARAMS = p.Results.STAT_PARAMS;
%% ===================================================================== %%
ersp_savef = {};
ersp_subbase_savef = {};
ersp_subbase_combase_savef = {};
ersp_singletrial_subbase_savef = {};
%- override
% STAT_PARAMS.singletrials = 'on';
%## PARAMS SETUP
tmpSTUDY = pop_statparams(STUDY, 'condstats', STAT_PARAMS.condstats,...
        'method',STAT_PARAMS.method,...
        'singletrials',STAT_PARAMS.singletrials,'mode',STAT_PARAMS.mode,...
        'fieldtripalpha',STAT_PARAMS.fieldtripalpha,'fieldtripmethod',STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',STAT_PARAMS.fieldtripnaccu);
tmpSTUDY_commonbase = pop_erspparams(tmpSTUDY, 'subbaseline','on',...
        'timerange',ERSP_PARAMS.timerange, 'ersplim',ERSP_PARAMS.ersplim);  % 'subbaseline' - ['on'|'off'] subtract the same baseline across conditions for ERSP     
%% ERSP CALCULATIONS
% NOTE (07/18/2023) std_ersplot_customParams.m has edits (find:'Chang') as
% well as a call to another edited function std_readata_customParams.m
%## ERSP calculation for no normalization 
fprintf('Gathering ERSP without any normalization for cluster %i\n',cluster_i)
tic
try
    [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot(STUDY,ALLEEG,'clusters',cluster_i,...
        'freqrange',ERSP_PARAMS.freqrange,'design',des_i); %PLOTS THE ENTIRE cluster_i'th CLUSTER ,% did I apply subtractsubjectmean???, no correction at all with the current setting
    %- save dat
    ersp_data = struct('allerspdata',{allerspdata},'alltimes',{alltimes},'allfreqs',{allfreqs},...
        'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
    par_save(ersp_data,save_dir,sprintf('ersp_data_cl%i_%s.mat',cluster_i,design_char));
    ersp_savef = [save_dir filesep sprintf('ersp_data_cl%i_%s.mat',cluster_i,design_char)];
    %- save fig
    fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
    saveas(fig_i,[save_dir filesep sprintf('allcomps_ersp_plot_cl%i_%s',cluster_i,design_char)])
    close(fig_i)
catch e
    fprintf(['error. code block 1\n',...
        'error. identifier: %s\n',...
        'error. %s\n',...
        'stack. %s\n'],e.identifier,e.message,getReport(e));
end
toc
%% ===================================================================== %%
%## ERSP calculation for normalization 
fprintf('Gathering ERSP after baseline correction using times {%0.2f,%0.2f] for cluster %i\n',warping_times(1),warping_times(5),cluster_i);
tic
try
    disp(tmpSTUDY.etc.statistics)
%                 disp(ersp_load_params.common_base);
    [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot_customParams(tmpSTUDY,ALLEEG,ersp_load_params,'clusters',cluster_i,...
        'freqrange',ERSP_PARAMS.freqrange,'design',des_i);
    %## DEBUG std_readdata
    %{
    params = tmpSTUDY.etc.erspparams;
    params.concatenate = 'off';
    stats = tmpSTUDY.etc.statistics;
    opt.datatype = 'ersp';  = []; opt.clusters = cluster_i; 
    [~, allersp, alltimes, allfreqs, events, paramsersp] = std_readdata_customParams(tmpSTUDY, ALLEEG, ersp_load_params, 'clusters', opt.clusters, 'datatype', opt.datatype, ...
            'component', opt.comps, 'singletrials', stats.singletrials, 'subbaseline', params.subbaseline, 'timerange', params.timerange, 'freqrange', params.freqrange, 'design', des_i, 'concatenate', params.concatenate);
    %}
    %- save dat
    ersp_data = struct('allerspdata',{allerspdata},'alltimes',{alltimes},'allfreqs',{allfreqs},...
        'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
    par_save(ersp_data,save_dir,sprintf('ersp_data_cl%i_%s_subbaselined.mat',cluster_i,design_char));
    ersp_subbase_savef = [save_dir filesep sprintf('ersp_data_cl%i_%s_subbaselined.mat',cluster_i,design_char)];
    %- save fig
    fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
    saveas(fig_i,[save_dir filesep sprintf('ersp_plot_cl%i_%s_allcomps_subbaselined',cluster_i,design_char)])
    close(fig_i)
catch e
    fprintf(['error. code block 2\n',...
        'error. identifier: %s\n',...
        'error. %s\n',...
        'stack. %s\n'],e.identifier,e.message,getReport(e));
end
toc 
%% ===================================================================== %%

fprintf('Gathering ERSP after baseline correction and common baseline using times {%0.2f,%0.2f] for cluster %i\n',warping_times(1),warping_times(5),cluster_i);
tic
try
    disp(tmpSTUDY_commonbase.etc.statistics)
    disp(tmpSTUDY_commonbase.etc.erspparams)
    [~,allerspdata,alltimes,allfreqs,pgroup,pcond,pinter] = std_erspplot_customParams(tmpSTUDY_commonbase,ALLEEG,ersp_load_params,'clusters',cluster_i,...
        'freqrange',ERSP_PARAMS.freqrange,'design',des_i); 
    %- save dat
    ersp_data = struct('allerspdata',{allerspdata},'alltimes',{alltimes},'allfreqs',{allfreqs},...
        'pgroup',{pgroup},'pcond',{pcond},'pinter',{pinter});
    par_save(ersp_data,save_dir,sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',cluster_i,design_char));
    ersp_subbase_combase_savef = [save_dir filesep sprintf('spec_data_cl%i_%s_subbaselined_commonbase.mat',cluster_i,design_char)];
    %- save fig
    fig_i = get(groot,'CurrentFigure');
%             fig_i.Position = [500 300 1480 920];
    saveas(fig_i,[save_dir filesep sprintf('ersp_plot_cl%i_%s_allcomps_subbaselined_commonbase',cluster_i,design_char)])
    close(fig_i)
catch e
    fprintf(['error. code block 3\n',...
        'error. identifier: %s\n',...
        'error. %s\n',...
        'stack. %s\n'],e.identifier,e.message,getReport(e));
end
toc

%% ===================================================================== %%
%- override
%{
STAT_PARAMS.singletrials = 'on';
tmpSTUDY = pop_statparams(STUDY, 'condstats', STAT_PARAMS.condstats,...
        'method',STAT_PARAMS.method,...
        'singletrials',STAT_PARAMS.singletrials,'mode',STAT_PARAMS.mode,...
        'fieldtripalpha',STAT_PARAMS.fieldtripalpha,'fieldtripmethod',STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',STAT_PARAMS.fieldtripnaccu);
params = tmpSTUDY.etc.erspparams;
params.concatenate = 'off';
stats = tmpSTUDY.etc.statistics;
opt.datatype = 'ersp';
opt.clusters = cluster_i; 
opt.comps = [];
[~, allersp, alltimes, allfreqs, ~, ~] = std_readdata_customParams(tmpSTUDY, ALLEEG, ersp_load_params, 'clusters', opt.clusters, 'datatype', opt.datatype, ...
        'component', opt.comps, 'singletrials', stats.singletrials, 'subbaseline', params.subbaseline, 'timerange', params.timerange, 'freqrange', params.freqrange, 'design', des_i, 'concatenate', params.concatenate);
[allerspdata,~,~] = mim_prune_ersp_trials(tmpSTUDY,allersp,cluster_i,des_i);
%- save dat
ersp_data = struct('allerspdata',{allerspdata},'alltimes',{alltimes},'allfreqs',{allfreqs},...
    'pgroup',{[]},'pcond',{[]},'pinter',{[]});
par_save(ersp_data,save_dir,sprintf('spec_data_cl%i_%s_singletrials_subbaselined.mat',cluster_i,design_char));
ersp_singletrial_subbase_savef = [save_dir filesep sprintf('spec_data_cl%i_%s_singletrials_subbaselined.mat',cluster_i,design_char)];
%}
%## time
toc(tic1)
end

