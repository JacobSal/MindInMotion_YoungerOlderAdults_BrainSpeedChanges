function [SUBJSTRUCT] = subj_loadStruct(fPath,fName,varargin)
%SUBJ_RESAVESTRUCT Summary of this function goes here
%   
%   IN:
%       SUBJSTRUCT
% 
%   OUT: 
%   IMPORTANT: 
%

%## TIME
tic
%## DEFINE DEFAULTS
%- validation
validPath = (@(x) ischar(x) || isempty(x));
Defaults = {[]};
p = inputParser;
%## REQUIRED
addRequired(p,'fPath',validPath)
addRequired(p,'fName',validPath)
%## OPTIONAL
addOptional(p,'path_prefix',Defaults{1},validPath)
parse(p,fPath,fName, varargin{:});

path_prefix = p.Results.path_prefix;
%## SET DEFAULTS
%% ===================================================================== %%
%## MAIN PROCESS
if strncmp(computer,'PC',2)
    if isempty(path_prefix)
        path_prefix = 'M';
    end
    SUBJSTRUCT = par_load(fPath,fName,path_prefix);
else  % isunix
    if isempty(path_prefix)
        path_prefix = 'dferris';
    end
    SUBJSTRUCT = par_load(fPath,fName,path_prefix);
end
fprintf('Loaded SUBJSTRUCT from %s\n',fullfile(fPath,fName));

end

