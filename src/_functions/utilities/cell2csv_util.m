function [char] = cell2csv_util(cell_arr)
%CELL2CSV Summary of this function goes here
%   Detailed explanation goes here
sz = size(cell_arr);
char = cell(sz(1)*sz(2),1);
if iscell(cell_arr)
    chkn = cellfun(@isnumeric,cell_arr);
    if any(chkn,'all')
        cell_arr = cellfun(@string,cell_arr);
        cell_arr = cellstr(cell_arr);
    end
elseif isnumeric(cell_arr)
    cell_arr = string(cell_arr);
    cell_arr = cellstr(cell_arr);
end

if sz(1) > 1 && sz(2) == 1
    for i = 1:sz(1)
        char{i,1} = [sprintf('%s,',cell_arr{i,1:sz(2)-1}), sprintf('%s;',cell_arr{i,sz(2)})];
    end
elseif sz(2) > 1 && sz(1) == 1
    char{1,1} = [sprintf('%s,',cell_arr{1,1:sz(2)-1}), sprintf('%s;',cell_arr{1,sz(2)})];
elseif sz(1) > 1 && sz(2) > 1
    char{1,1} = [sprintf('%s,',cell_arr{1,1:sz(2)-1}), sprintf('%s;',cell_arr{1,sz(2)})];
else
    char{1,1} = [sprintf('%s',cell_arr{1})];
end
char = [char{:}];
end

