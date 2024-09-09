function [args] = eeglab_struct2args(struct)
%EEGLAB_STRUCT2ARGS Summary of this function goes here
%   Detailed explanation goes here
%% Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'struct',@isstruct);
parse(p,struct);
%%
fn = fieldnames(struct);
sc = struct2cell(struct);
args = cell(length(fn)+length(sc),1);
cnt = 1;
for i = 1:length(fn)
    args{cnt} = fn{i};
    args{cnt+1} = sc{i};
    cnt = cnt + 2;
end
end

