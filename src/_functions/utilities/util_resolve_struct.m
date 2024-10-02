function [struct_out] = util_resolve_struct(struct_in,varargin)
%EEGLAB_RESOLVE_STRUCT Summary of this function goes here
%
%   IN:
%       struct_in, CELL ARRAY of STRUCTS
%           cell array of structures that you want to concatentate, but
%           they have mismatched fields. The algorithm will add missing
%           fields to indices lacking them.
%       Optional; struct_ids, CELL ARRAY of CHARS
%           cell array of characters that can be used for verbose debugging
%           purposes. (e.g., cell array of subject names).
%   OUT:
%       struct_out, STRUCT ARRAY
%
%   Version History --> See details at the end of the script.
%   Previous Version: n/a
%   Summary:  
%
struct_ids = {};
struct_ids_vfcn = @(x) iscell(x) && length(x) == length(struct_in);
cat_logo();
%## TIME
t = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'struct_in',@iscell);
%## OPTIONAL
addOptional(p,'struct_ids',struct_ids,struct_ids_vfcn)
%## PARSE
parse(p, struct_in, varargin{:});
%## SET DEFAULTS
struct_ids = p.Results.struct_ids;
%% ===================================================================== %%
struct_fields = cell(1,length(struct_in));
for i = 1:length(struct_in)
    struct_fields{i} = fields(struct_in{i})';
    disp(size(fields(struct_in{i})'));
end
struct_fields = unique([struct_fields{:}]);
iter_fields = struct_fields;
for i = 1:length(struct_in)
    tmp_struct = struct_in{i};
    tmp_fields = fields(tmp_struct);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,tmp_fields)),iter_fields,'UniformOutput',false); 
    out = [out{:}];
    field_add = iter_fields(~out);
    if any(~out)
        for j = 1:length(field_add)
            tmp_struct.(field_add{j}) = [];
            if isempty(struct_ids)
                fprintf('Index %i) Adding fields %s\n',i,field_add{j})
            else
                fprintf('%s) Adding fields %s\n',struct_ids{i},field_add{j})
            end
        end
    end 
%     struct_in{subj_i} = EEG;
    struct_in{i} = orderfields(tmp_struct);
end
%- CONCATENATE struct_in
struct_out = cellfun(@(x) [[]; x], struct_in);
fprintf('util_resolve_struct processing time: %0.2f\n\n',toc(t));
end

