function [b] = validate_struct(x,DEFAULT_STRUCT)
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT: 
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Chang Liu, Jacob Salminen
% Code Date: 04/28/2023, MATLAB R2020b
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
% Copyright (C) Chang Liu, liu.chang1@ufl.edu    
b = false;
struct_name = inputname(2);
%##
fs1 = fields(x);
fs2 = fields(DEFAULT_STRUCT);
vals1 = struct2cell(x);
vals2 = struct2cell(DEFAULT_STRUCT);
%- check field names
fprintf('Running structure checks...\n');
chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
if ~all(chk)
    fprintf(2,'Fields for struct do not match for %s\n',struct_name);
    return
end
%- check field value's class type
for f = 1:length(fs2)
    ind = strcmp(fs2{f},fs1);
    if any(ind)
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind}));
            return
        end
    else
        fprintf(2,'Struct.%s is not present in input structure. Setting to default.\n',fs2{f});
    end
end
b = true;
end

