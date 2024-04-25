function [t] = sprintf_table(table_in)
%SPRINTF_TABLE Summary of this function goes here
%   Detailed explanation goes here
%
% see. PrintTable.m
t = PrintTable;
t.addRow('',table_in.Properties.VariableNames{:});
for i = 1:size(table_in,1)
    rn = table_in.Properties.RowNames;
    row = num2cell([table_in{i,:}]);
    t.addRow(rn{i},row{:});
end
% t.display;
t.HasHeader = true;
% t.Format = 'tex';
% t.Caption = '';
% t.print;
end

