function [stat,cfg] = eeglab_ft_stats(data)
%EEGLAB_FT_STATS Summary of this function goes here
%   This is an annotated wrapper for ft_fieldtripstatistics that bipases
%   some eeglab mechanisms that can cause bugs (USE WITH CAUTION)

%-
g.chandim = 0;
g.chanlocs = struct.empty;
% alter: design, uvar, tail, correcttail, ivar, clustercritval,
% effect for various statistical tests
cfg = [];
% 'indepsamplesT'           independent samples T-statistic,
% 'indepsamplesF'           independent samples F-statistic,
% 'indepsamplesregrT'       independent samples regression coefficient T-statistic,
% 'indepsamplesZcoh'        independent samples Z-statistic for coherence,
% 'depsamplesT'             dependent samples T-statistic,
% 'depsamplesFmultivariate' dependent samples F-statistic MANOVA,
% 'depsamplesregrT'         dependent samples regression coefficient T-statistic,
% 'actvsblT'                activation versus baseline T-statistic.
%## SPECIFIC TO FT_STATISTICS_MONTECARLO
% cfg.statistic = 'depsamplesT';
cfg.statistic = 'depsamplesFunivariate';
cfg.correcttail = 'alpha';
%- how to combine single samples that belong to a cluster:
%'maxsum','maxsize','wcm' (i.e., weighted cluster mass). see.
%Hayasaka & Nichols (2004) NeuroImage.
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.clusteralpha = 0.05;
% cfg.computecritval = 'yes';
% cfg.clustercritval = [];
% cfg.clustertail = 0;
cfg.correctm = 'cluster';
cfg.tail        = 1;
cfg.correcttail = 'no';
cfg.feedback = 'no';
cfg.alpha = 0.05;
%- independent variables
cfg.ivar = 1;
%- unit variables
cfg.uvar = 2;
% cfg.wvar = [];
% cfg.cvar = [];
% cfg.randomseed = 'yes'
% cfg.resampling = 'permutation';
% cfg.computestat = 'yes';
cfg.mcorrect = 'cluster'; 
cfg.numrandomization = 2000;
cfg.precondition = [];
% cfg.dim = [1 122 60];
% cfg.dimord = 'chan_freq_time';
%## SPECIFIC TO FT_FREQSTATISTICS
%- determines the ft function called (e.g.,
%ft_statistics_montecarlo)
cfg.channel = {'c1'};
cfg.latency = [1,60];
cfg.frequency = [1,122];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter = 'powspctrm';
cfg.method = 'montecarlo';
cfg.design = [];
cfg.correcttail = 'no';
%- odd logic for selecting the right statistic based on the availability of
%the ft_statfun_depsamplesFmultivariate function
% tmpP = fileparts(which('ft_freqstatistics'));
% if exist([tmpP filesep 'statfun' filesep 'ft_statfun_depsamplesFmultivariate.m'],'file')
%     cfg.statistic   = 'depsamplesFunivariate';
% else
%     cfg.statistic   = 'depsamplesF';
% end
%## 
%- 
%- newdata is a cell (groups x conds) of structs with fields: dimord, powspctrmdimord, powspctrm, label, freq, time
% powspctrm is size subjs x 1? x freqs x time
[newdata,design1,design2,design3] = makefieldtripdata(data, g.chandim, g.chanlocs);
cfg.design = [design1; design3];
% design1 is a 1xN vector containing designs for the conditions
% (e.g., 4 conditions with 3 subjects would require a vector
% size 1x12 with numbers 1-4.) (i.e., condition stats, paired and unpaired)
% design2 is a 1xN vector containing designs for the group (i.e.,
% group stats)
% design3 is a 1xN vector containing designs for the subjects
% (i.e., paired)
% see. statcondfieldtrip for more configuration examples.
stat = ft_freqstatistics(cfg, newdata{:});
%- 
% data: Samples x Vars matrix (e.g., 4 conds, 18 subjs, 60 latencies, 122 freqs
% (uvar = 2, ivar = 1 (speed)) = [60*122=7320 x 4*18=72]
% stat = ft_statistics_montecarlo();
end
%% =========================================== %%
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1
            val = 2;
        elseif size(a,2) == 1
            val = 1;
        else
            val = 2;
        end
    end
end
function [newdata, design1, design2, design3] = makefieldtripdata(data, chandim, chanlocs)
    newdata = {};
    swapdim = [];
    for i = 1:length(data(:))
        newdata{i}.dimord    = 'rpt_chan_freq_time';
        newdata{i}.powspctrmdimord = 'rpt_chan_freq_time';
        switch myndims(data{1})
          case 1
            newdata{i}.powspctrm = data{i};
            
          case 2
            if chandim
                 newdata{i}.powspctrm = transpose(data{i});
            else 
                newdata{i}.powspctrm = reshape(transpose(data{i}), size(data{i},2), 1, size(data{i},1));
            end
            
          case 3
            if chandim == 2 % chandim can be 1 or 2
                swapdim = [2 1];
            end
            if chandim
                 newdata{i}.powspctrm = permute(data{i}, [3 1 2]);
            else 
                newdata{i}.powspctrm = permute(data{i}, [3 4 1 2]); % 4 is a singleton dimension
            end
            
          case 4
            newdata{i}.powspctrm = permute(data{i}, [4 3 1 2 ]); % Fixed dimension from [4 1 2 3]
        end
        
        newdata{i}.label     = cell(1,size(newdata{i}.powspctrm,2));
        newdata{i}.label(:)  = { 'cz' };
        for ic = 1:length(newdata{i}.label)
            newdata{i}.label{ic} = [ 'c' num2str(ic) ];
        end
        newdata{i}.freq      = [1:size(newdata{i}.powspctrm,3)];
        newdata{i}.time      = [1:size(newdata{i}.powspctrm,4)];
        
        % below in case channels are specified
        % not that statistics are done on time x frequencies or channels
        % so time x frequency x channels do not work yet here
        if ~isempty(chanlocs)
            newdata{i}.powspctrm = squeeze(newdata{i}.powspctrm);
            newdata{i}.label     = { chanlocs.labels };
            newdata{i}.freq      = 1;
            newdata{i}.time      = 1;
        end
        if isempty(chanlocs) && size(newdata{i}.powspctrm,2) ~= 1
            newdata{i}.dimord    = 'rpt_freq_time';
            newdata{i}.powspctrmdimord = 'rpt_freq_time';
        end
    end
    design1 = [];
    design2 = [];
    design3 = [];
    for i = 1:size(data,2)
        for j = 1:size(data,1)
            nrepeat = size(data{i}, ndims(data{i}));
            ij = j+(i-1)*size(data,1);
            design1 = [ design1 ones(1, nrepeat)*i ];
            design2 = [ design2 ones(1, nrepeat)*j ];
            design3 = [ design3 [1:nrepeat] ];
        end
    end
end
