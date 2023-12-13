function [ALLEEG] = format_alleeg_cell(ALLEEG)
%format_alleeg_cell Summary of this function goes here
%   Detailed explanation goes here
% Code Designer: Jacob Salminen (11/25/2022)
%## TIME
tic
%- ALLEEG
errorMsg = 'Value must be of format {EEG_STRUCT,EEG_STRUCT,...}. EEG strucutres you''d like to concatenate and ensure formatting on.';
al_validFcn = @(x) assert(iscell(x) && isstruct(x{1}), errorMsg);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',al_validFcn);
parse(p,ALLEEG);
%% ===================================================================== %%
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
%- get current fields
fss = cell(1,length(ALLEEG));
for subj_i = 1:length(ALLEEG)
    fss{subj_i} = fields(ALLEEG{subj_i})';
    disp(size(fields(ALLEEG{subj_i})'));
end
fss = unique([fss{:}]);
fsPrev = fss;
%- add fields and concatenate
for subj_i = 1:length(ALLEEG)
    EEG = ALLEEG{subj_i};
    fs = fields(EEG);
    %- add fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fs)),fsPrev,'UniformOutput',false);
    out = [out{:}];
    addFs = fsPrev(~out);
    if any(~out)
        for j = 1:length(addFs)
            EEG.(addFs{j}) = [];
            fprintf('%s) Adding fields %s\n',EEG.subject,addFs{j});
        end
    end
    ALLEEG{subj_i} = orderfields(EEG);
end
%- CONCATENATE ALLEEG
ALLEEG = cellfun(@(x) [[]; x], ALLEEG);
%##
toc
end

