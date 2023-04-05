function [compsOut,clustOut,distOut,posxyz,anatTable] = get_Dist2Anatomy(STUDY,ALLEEG,INDEX,anatomyCoords,anatomyName,subjectString,varargin)
%GETDIST2ANATOMY Summary of this function goes here
%   Detailed explanation goes here
%## TIME
tic
%## DEFINE DEFAULTS
Defaults = {25};
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY')
addRequired(p,'ALLEEG')
addRequired(p,'INDEX')
addRequired(p,'anatomyCoords')
addRequired(p,'anatomyName')
addRequired(p,'subjectString');
%## OPTIONAL
%## PARAMETER
addParameter(p,'threshold',Defaults{1},@isnumeric)
parse(p, STUDY, ALLEEG, INDEX, anatomyCoords, anatomyName, subjectString, varargin{:});
%## SET DEFAULTS
threshold = p.Results.threshold;
%% ===================================================================== %%
%- params
subj_i = INDEX;
components = {};
clustNum = {};
distOut = {};
anatBool = {};
%- loop through clusters and get clusters and components for subj_i
for j = 2:length(STUDY.cluster) % skip parentcluster    
    compI = (STUDY.cluster(j).sets == subj_i);
    tmpcomps = STUDY.cluster(j).comps(compI);
    components = [components, tmpcomps];
    clustNum = [clustNum, repmat(j,1,length(tmpcomps))];
end
%- extract information
components = [components{:}];
clustNum = [clustNum{:}];
%- perform distance calculation between each dipole and anatomy Coordinate
posxyz = {ALLEEG(subj_i).dipfit.model(components).posxyz};
for k = 1:length(posxyz)
    tmpD = sqrt(sum((anatomyCoords-posxyz{k}).^2));
    distOut = [distOut, tmpD];
end
distOut = [distOut{:}];
%- create indexing array
for k = 1:length(distOut)
    anatBool = [anatBool (distOut(k) < threshold)];
end
%- create table
anatTable = table(components',clustNum',distOut',anatBool','VariableNames',{'ComponentNumber','ClusterNumber','Distance2Anatomy','ThreshBool'});
%- extract thresheld distances
anatBool = [anatBool{:}];
idx = logical(anatBool);
compsOut = components(idx);
clustOut = clustNum(idx);
distOut = distOut(idx);
t = table(compsOut',clustOut',distOut','VariableNames',{'ComponentNumber','ClusterNumber','Distance2Anatomy'});
fprintf('Subject %s) Anatomy %s Assignment:\n',subjectString,anatomyName);
disp(t)
% fprintf('Subject %s) Anatomy Assignment:\n',subjectString);
% fprintf('Components | '); fprintf('%i |',compsOut); fprintf('\n');
% fprintf('Cluster number | '); fprintf('%i |',clustOut); fprintf('\n');
% fprintf('Distance to anatomy (mm) | '); fprintf('%0.1f |',distOut{idx}); fprintf('\n');
end

