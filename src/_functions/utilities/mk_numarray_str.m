function [numarray_str] = mk_numarray_str(numarray)
%MK_NUMARRAY_STR Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
%## REQUIRED
addRequired(p,'numarray',@isnumeric);
%## CHECK
parse(p,numarray);
%% ===================================================================== %%
if any(size(numarray) > 1)
    numarray_str = [sprintf('[%s,',string(numarray(1:end-1))),sprintf('%s]',string(numarray(end)))];
else
    numarray_str = [sprintf('[%s]',string(numarray(1)))];
end
end

