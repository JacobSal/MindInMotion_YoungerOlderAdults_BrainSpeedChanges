function [cell_of_struct] = chk_struct_fields(cell_of_struct,varargin)
%CHK_STRUCT_FIELDS Summary of this function goes here
%   Detailed explanation goes here
%## BOOKKEEPING (i.e., ADD fields not similar across EEG structures)
fss = cell(1,length(cell_of_struct));
for subj_i = 1:length(cell_of_struct)
    fss{subj_i} = fields(cell_of_struct{subj_i});
end
fss = unique([fss{:}]);
fsPrev = fss;
for subj_i = 1:length(cell_of_struct)
    tmp_s = cell_of_struct{subj_i};
    fs = fields(tmp_s);
    % delete fields not present in other structs.
    out = cellfun(@(x) any(strcmp(x,fsPrev)),fs,'UniformOutput',false); 
    out = [out{:}];
    addFs = fs(~out);
    if any(~out)
        for j = 1:length(addFs)
            tmp_s.(addFs{j}) = [];
            fprintf('%s) Adding %s %s\n',tmp_s.subject,addFs{j})
        end
    end 
    cell_of_struct{subj_i} = tmp_s;
end
end

