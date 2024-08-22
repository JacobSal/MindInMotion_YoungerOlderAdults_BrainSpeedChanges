function [anatomy_struct] = eeglab_get_anatomy(STUDY,varargin)
%EEGLAB_GET_ANATOMY Summary of this function goes here
%
%   IN:
%       STUDY, EEGLAB STUDY STRUCT
%           
%       Optional; ALLEEG,EEGLAB ALLEEG STRUCT
%           use this option if given error that dipole locations are needed
%           to produce centroid data or cluster dipole data.
%   OUT:
%       struct_out, STRUCT ARRAY
%
%   Version History --> See details at the end of the script.
%   Previous Version: n/a
%   Summary:  
%
%## PARAMS
%- find fieldtrip on path
if ~ispc
    tmp = strsplit(path,':');
else
    tmp = strsplit(path,';');
end
b1 = regexp(tmp,'fieldtrip','end');
b2 = tmp(~cellfun(@isempty,b1));
path_fieldtrip = b2{1}(1:b1{1});
fprintf('fieldtrip path: %s\n',path_fieldtrip);
%- atlas paths
ATLAS_FPATHS = {[path_fieldtrip filesep 'template',...
        filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
    [path_fieldtrip filesep 'template' filesep 'atlas'...
        filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
    [path_fieldtrip filesep 'template' filesep 'atlas'...
        filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
% see. https://www.fieldtriptoolbox.org/template/atlas/
atlas_fpath_vfcn = @(x) 
%- alleeg default
ALLEEG = struct();
%- logo
cat_logo();
%## TIME
t = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
%## OPTIONAL
addOptional(p,'ALLEEG',ALLEEG,@isstruct)
%## PARAMETER
addParameter(p,'ATLAS_FPATHS',ATLAS_FPATHS)
%## PARSE
parse(p, struct_in, varargin{:});
%## SET DEFAULTS
ALLEEG = p.Results.ALLEEG;
% [STUDY,centroid] = std_centroid(STUDY,ALLEEG,double(string(clusters)),'dipole');
%% ===================================================================== %%
atlas_i = 1;
atlas_name_store = cell(length(main_cl_inds),1);
titles_store = cell(length(main_cl_inds),1);
%- cfg
cfg            = [];
cfg.roi        = dip1;
cfg.output     = 'single';
cfg.atlas      = atlas;
cfg.inputcoord = 'mni';
cfg.verbose = 0;
%- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
cfg.sphere = 3;
%- loop
for k_i = 1:length(main_cl_inds)
    cl_i = main_cl_inds(k_i);
    %## ANATOMY
    STUDY.cluster(cl_i).centroid.dipole.posxyz = mean(STUDY.cluster(cl_i).all_diplocs);
    % dip1 = STUDY.cluster(cl_i).all_diplocs;
    atlas = ft_read_atlas(ATLAS_FPATHS{atlas_i});
    atlas_name = 'error';
    label_i = ft_volumelookup(cfg, atlas);
    if ~isempty(label_i)
        counts = sum([label_i.count],2);
        [val, indx] = max(counts);
        names = label_i(1).name;
        if strcmp(names(indx),'no_label_found')
            sub_indx = find(counts ~= 0 & counts < val);
            if ~isempty(sub_indx)
                atlas_name = names{sub_indx};
            end
        else
            atlas_name = names{indx};
        end
    end
    atlas_name_store{cl_i} = atlas_name;
    titles_store{cl_i} = sprintf('CL%i: %s\n[%0.1f,%0.1f,%0.1f]\n',cl_i,atlas_name,STUDY.cluster(cl_i).centroid.dipole.posxyz);
end
end

