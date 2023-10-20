function [STUDY,atlas_name,centroids] = add_anatomical_labels(STUDY,ALLEEG,varargin)
%ADD_ANATOMICAL_LABELS Summary of this function goes here
%   Detailed explanation goes here
%   IN: 
%   OUT: 
%   IMPORTANT
%
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designer: Noelle Jacobsen, Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Noelle Jacobsen

%## TIME
tic
%## DEFINE DEFAULTS
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%-
%## Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'FIELDTRIP_PATH',FIELDTRIP_PATH,@ischar);
parse(p,STUDY,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS

%% ===================================================================== %%
[STUDY, ~] = std_centroid(STUDY,ALLEEG,1:length(STUDY.cluster),'dipole');
atlas = ft_read_atlas([FIELDTRIP_PATH filesep 'template' filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii']);
%- extract centroid locations
centroids = [STUDY.cluster(1:end).centroid];
centroids = [centroids.dipole];
tmp_i = find(isnan([centroids.rv]));
for i = tmp_i
    centroids(i).posxyz = [0,0,0];
    centroids(i).momxyz = [0,0,0];
end
dipfit_roi = cat(1,centroids.posxyz);
%- params
atlas_name = cell(size(dipfit_roi,1),2);
for i = 2:length(STUDY.cluster)
    cfg              = [];
    %- (01/19/2023), JS: Perhaps have the function loop through each dipoel
    %in the cluster, determine a location of interest, and do it that way?
    %Doesn't seem like passing in a set of points helps...
%     cfg.roi = STUDY.cluster(i).preclust.preclustdata;
    cfg.roi        = dipfit_roi(i,:);
    cfg.output     = 'multiple';
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
    cfg.sphere = 3;
    labels = ft_volumelookup(cfg, atlas);
    if ~isempty(labels)
        % disp(labels.count(labels.count ~= 0))
        [val, indx] = max(labels.count);
        if strcmp(labels.name(indx),'no_label_found')
            sub_indx = find(labels.count ~= 0 & labels.count < val);
            if isempty(sub_indx)
                atlas_name{i,1} = ['CLs ',num2str(i)];
                atlas_name{i,2} = labels.name{indx};
                continue;
            end
            atlas_name{i,1} = ['CLs ',num2str(i)];
            tmp = sprintf('(N=%i) %s',labels.count(sub_indx(1)),labels.name{sub_indx(1)});
            for ic_i = sub_indx(2:end)
                tmp = [tmp,' & ', sprintf('(N=%i) %s',labels.count(ic_i),labels.name{ic_i})];
            end
            atlas_name{i,2} = tmp;
        else
            atlas_name{i,1} = ['CLs ',num2str(i)];
            atlas_name{i,2} = labels.name{indx};
        end
    else
        warning('No label found for cluster %i',i);
    end
end

fprintf('Cluster \t\t Label\n');
fprintf('________________________\n');
for i = 2:size(dipfit_roi,1)
    label = cellstr(atlas_name{i,2});
    cl =  cellstr(atlas_name{i,1});
    fprintf('%s\t\t%s\n',cl{:},label{:})
end
end

