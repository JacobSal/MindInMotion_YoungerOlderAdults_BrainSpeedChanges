% exerpt from std_erspplot 

function [pcond, pgroup, pinter] = erspStats(STUDY,allersp)

        %get stats parameters
        stats = STUDY.etc.statistics;
        stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
        if isempty(STUDY.design(STUDY.currentdesign).variable)
            stats.paired = { };
        else
            stats.paired = { STUDY.design(STUDY.currentdesign).variable(:).pairing };
        end

        %get ersp params
        params = STUDY.etc.erspparams;
        params.plottf =[];
        % select specific time and freq
        % -----------------------------

        if ~isempty(params.plottf)
            if length(params.plottf) < 3
                params.plottf(3:4) = params.plottf(2);
                params.plottf(2)   = params.plottf(1);
            end
            [~, fi1] = min(abs(allfreqs-params.plottf(1)));
            [~, fi2] = min(abs(allfreqs-params.plottf(2)));
            [~, ti1] = min(abs(alltimes-params.plottf(3)));
            [~, ti2] = min(abs(alltimes-params.plottf(4)));
            for index = 1:length(allersp(:))
                allersp{index} = mean(mean(allersp{index}(fi1:fi2,ti1:ti2,:,:),1),2);
                allersp{index} = reshape(allersp{index}, [1 size(allersp{index},3) size(allersp{index},4) ]);
            end

            % prepare channel neighbor matrix for Fieldtrip
            statstruct = std_prepare_neighbors(STUDY, ALLEEG);
            stats.fieldtrip.channelneighbor = statstruct.etc.statistics.fieldtrip.channelneighbor;

            params.plottf = { params.plottf(1:2) params.plottf(3:4) };
            [pcond, pgroup, pinter] = std_stat(allersp, stats);
            if (~isempty(pcond) && length(pcond{1}) == 1) || (~isempty(pgroup) && length(pgroup{1}) == 1), pcond = {}; pgroup = {}; pinter = {}; end % single subject STUDY
        else
            [pcond, pgroup, pinter] = std_stat(allersp, stats);
            if (~isempty(pcond ) && (size( pcond{1},1) == 1 || size( pcond{1},2) == 1)) || ...
                    (~isempty(pgroup) && (size(pgroup{1},1) == 1 || size(pgroup{1},2) == 1))
                pcond = {}; pgroup = {}; pinter = {};
                disp('No statistics possible for single subject STUDY');
            end % single subject STUDY
        end
    end