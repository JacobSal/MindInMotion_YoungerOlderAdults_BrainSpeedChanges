function [allSubj,stat,idxval] = getDirInfMIM(PATHS)
%GETSUBJDIRDATA Summary of this function goes here
%   Detailed explanation goes here
%   
%   IN: 
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct)
%## OPTIONAL
%## PARAMETER
parse(p, PATHS);
%## SET DEFAULTS
%## SAVE FOLDER HANDLER (HIGHEST ORDER FUNCTION ONLY)
% ----------------------------------------------------------------------- %
%## Get all subject path data and organize into health and non-healthy
% parsing character arrays for regexp
patH = '[H]+\d*';
patNH = '[N][H]+\d*';
pat1 = '[H]+1+\d*';
pat2 = '[H]+2+\d*';
pat3 = '[H]+3+\d*';
subjInf = dir(PATHS.path4MIMData);
sN = {subjInf.name};
%## parse strings to get age and health information
% Healthy Subjects
Hidx = cellfun(@(x) regexp(x,patH),sN,'UniformOutput',false); Hidx = cellfun(@(x) isequal(x,1), Hidx,'UniformOutput',false); Hidx = [Hidx{:}];
% Not Healthy Subjects
NHidx = cellfun(@(x) regexp(x,patNH),sN,'UniformOutput',false); NHidx = cellfun(@(x) isequal(x,1), NHidx,'UniformOutput',false); NHidx = [NHidx{:}];
% Young Subjects
idxA1 = cellfun(@(x) regexp(x,pat1),sN,'UniformOutput',false); idxA1 = cellfun(@(x) isequal(x,1), idxA1,'UniformOutput',false); idxA1 = [idxA1{:}];
% Middle Aged Subjects
idxA2 = cellfun(@(x) regexp(x,pat2),sN,'UniformOutput',false); idxA2 = cellfun(@(x) isequal(x,1), idxA2,'UniformOutput',false); idxA2 = [idxA2{:}];
% Older Subjects
idxA3 = cellfun(@(x) regexp(x,pat3),sN,'UniformOutput',false); idxA3 = cellfun(@(x) isequal(x,1), idxA3,'UniformOutput',false); idxA3 = [idxA3{:}];

subjInfH = subjInf(Hidx);
subjInfNH = subjInf(NHidx);
allSubj = cat(1,subjInfH,subjInfNH);

fprintf('there are %i Health (H) patients, and %i Not Health (NH) patients\n',length(subjInfH),length(subjInfNH));
fprintf('there are %i ''Age 1'' patients, %i ''Age 2'' patients, %i ''Age 3'' patients\n',sum(idxA1),sum(idxA2),sum(idxA3));

stat = struct('nHealthy', length(subjInfH), ...
              'nNotHealthy', length(subjInfNH),...
              'nAge1', sum(idxA1),...
              'nAge2', sum(idxA2),...
              'nAge3', sum(idxA3));
idxval = {Hidx, NHidx, idxA1, idxA2, idxA3};
            
%## TIME
toc
end

