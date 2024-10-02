function [STUDY,anatomy_struct,dipfit_structs,topo_cells,txt_store] = eeglab_get_anatomy(STUDY,varargin)
%EEGLAB_GET_ANATOMY Summary of this function goes here
%
%   IN:
%       STUDY, EEGLAB STUDY STRUCT
%           
%       Parameter; ALLEEG,EEGLAB ALLEEG STRUCT
%           use this option if given error that dipole locations are needed
%           to produce centroid data or cluster dipole data.
%
%       Parameter; ANATOMY_STRUCT, STRUCT
%           struct containing fields defining atlases to use for anatomy
%           labeling, group characters for output, clusters to identify
%           anatomy on, whether to save information internal to this
%           function, and data availble if you have ran this function
%           before.
%   OUT:
%       struct_out, STRUCT ARRAY
%
%   Version History --> See details at the end of the script.
%   Previous Version: n/a
%   Summary:  
%
%## PARAMS
%-
anatomy_struct = struct.empty; %#ok<NASGU>
txt_store = {};
dipfit_structs = [];
topo_cells = [];
%- find fieldtrip on path
if ~ispc
    tmp = strsplit(path,':');
else
    tmp = strsplit(path,';');
end
%*
b1 = regexp(tmp,'fieldtrip','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_fieldtrip = b2{1}(1:b1{1});
    fprintf('fieldtrip path: %s\n',path_fieldtrip);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('fieldtrip path not found.\n');
    end
end
%*
b1 = regexp(tmp,'AAL3','end');
b2 = tmp(~cellfun(@isempty,b1));
try
    path_aal3 = b2{1}(1:b1{1});
    fprintf('ALL3 path: %s\n',path_aal3);
catch ME
    switch ME.identifier
        case 'MATLAB:badsubscript'
            fprintf('AAL3 path not found.\n');
    end
end
%- atlas paths test
% ANATOMY_STRUCT.atlas_fpath = {[path_fieldtrip filesep 'template',...
%         filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii'],... % MNI atalst
%     [path_fieldtrip filesep 'template' filesep 'atlas'...
%         filesep 'afni' filesep 'TTatlas+tlrc.HEAD'],... % tailarch atlas
%     [path_fieldtrip filesep 'template' filesep 'atlas'...
%         filesep 'spm_anatomy' filesep 'AllAreas_v18_MPM.mat']}; % also a discrete version of this
% see. https://www.fieldtriptoolbox.org/template/atlas/
% ANATOMY_STRUCT.atlas_fpath = {[path_aal3 filesep 'AAL3v1.nii'],...
%     [path_aal3 filesep 'AAL3v1_1mm.nii'],...
%     [path_aal3 filesep 'ROI_MNI_V7.nii'],...
%     [path_aal3 filesep 'ROI_MNI_V7_1mm.nii']}; % also a discrete version of this
%- test AAL3 spm
% addpath('M:\jsalminen\GitHub\par_EEGProcessing\submodules\spm12');
%## ANATOMY STRUCT
inds = cellfun(@(x) contains(x,'Outlier','IgnoreCase',true) || contains(x,'Parentcluster'),{STUDY.cluster.name});
% DEF_ANATOMY_STRUCT = struct('atlas_fpath',{{[path_aal3 filesep 'AAL3v1.nii'],...
%                                             [path_aal3 filesep 'AAL3v1_1mm.nii'],...
%                                             [path_aal3 filesep 'ROI_MNI_V7.nii'],...
%                                             [path_aal3 filesep 'ROI_MNI_V7_1mm.nii']}},...
%     'group_chars',{unique({STUDY.datasetinfo.group})},...
%     'cluster_inds',find(~inds),...
%     'save_inf',true,...
%     'topo_cells',{{}},...
%     'dipfit_structs',struct.empty);
DEF_ANATOMY_STRUCT = struct('atlas_fpath',{{[path_aal3 filesep 'AAL3v1.nii']}},...
    'group_chars',{unique({STUDY.datasetinfo.group})},...
    'cluster_inds',find(~inds),...
    'anatomy_calcs',{{'all centroid','all aggregate'}},... % ('all calcs','group centroid','all centroid','group aggregate','all aggregate')
    'save_inf',true,...
    'save_dir',STUDY.filepath,...
    'topo_cells',{{}},...
    'dipfit_structs',struct.empty);
%## ALLEEG DEFAULT
ALLEEG = struct.empty;
%- logo
cat_logo();
%## TIME
tt = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'ALLEEG',ALLEEG,@isstruct)
addParameter(p,'ANATOMY_STRUCT',DEF_ANATOMY_STRUCT,@(x) validate_struct(x,DEF_ANATOMY_STRUCT));
%## PARSE
parse(p, STUDY, varargin{:});
%## SET DEFAULTS
ANATOMY_STRUCT = p.Results.ANATOMY_STRUCT;
ALLEEG = p.Results.ALLEEG;
%-
ANATOMY_STRUCT = set_defaults_struct(ANATOMY_STRUCT,DEF_ANATOMY_STRUCT);
%- extract params for readability
cluster_inds = ANATOMY_STRUCT.cluster_inds;
group_chars = ANATOMY_STRUCT.group_chars;
%% PRELOAD ALLEEG SET INFORMATION ====================================== %%
CALC_STRUCT = struct('cluster_inds',ANATOMY_STRUCT.cluster_inds,...
    'save_inf',ANATOMY_STRUCT.save_inf,...
    'recalculate',false);
[STUDY,~,~] = eeglab_get_topodip(STUDY,...
    'CALC_STRUCT',CALC_STRUCT,...
    'ALLEEG',ALLEEG);
%% ATLAS CALCULATION =================================================== %%
def_anatomy_struct = struct('cluster',[],...
    'dips',[],...
    'subjects',{''},...
    'group',{''},...
    'group_n',[],...
    'calculation',{''},...
    'anatomy_label',{''},...
    'atlas',{''},...
    'atlas_label',{''});
anatomy_struct = def_anatomy_struct;
cnt = 1;
% atlas_name_store = cell(length(cluster_inds),1);
if ANATOMY_STRUCT.save_inf
    txt_store = cell(length(cluster_inds),1);
    f = fopen([ANATOMY_STRUCT.save_dir filesep 'anatomy_output.txt'],'w');
end
cnttxt = 1;
chk_a = any(strcmp(ANATOMY_STRUCT.anatomy_calcs,'all calcs'));
chk_aa = any(strcmp(ANATOMY_STRUCT.anatomy_calcs,'all aggregate'));
chk_ac = any(strcmp(ANATOMY_STRUCT.anatomy_calcs,'all centroid'));
chk_gc = any(strcmp(ANATOMY_STRUCT.anatomy_calcs,'group centroid'));
chk_ga = any(strcmp(ANATOMY_STRUCT.anatomy_calcs,'group aggregate'));
for atlas_i = 1:length(ANATOMY_STRUCT.atlas_fpath)
    atlas_fpath = ANATOMY_STRUCT.atlas_fpath{atlas_i};
    tmp = strsplit(atlas_fpath,filesep);
    tmpn = strsplit(tmp{end},'.');
    tmpp = strjoin(tmp(1:end-1),filesep);
    %## READ XML DATA FOR ATLAS LABELS
    tree = xmlread([tmpp filesep tmpn{1},'.xml']);
    label_xml = tree.getElementsByTagName('label');
    coord_xml = tree.getElementsByTagName('coordinate_system');
    ver_xml = tree.getElementsByTagName('version');
    nn_xml = tree.getElementsByTagName('name');
    def_atlas_xml_struct = struct('index',[],...
        'name',{''},...
        'atlas_name',{char(nn_xml.item(0).getFirstChild.getData)},...
        'atlas_version',{char(ver_xml.item(0).getFirstChild.getData)},...
        'coord_system',{char(coord_xml.item(0).getFirstChild.getData)});
    atlas_xml_struct = def_atlas_xml_struct;
    cntas = 1;
    for i=0:label_xml.getLength-1
        thisLabelItem = label_xml.item(i);
        indexNode = thisLabelItem.getElementsByTagName('index');
        indexNode = indexNode.item(0);
        nameNode = thisLabelItem.getElementsByTagName('name');
        nameNode = nameNode.item(0);
        atlas_xml_struct(cntas).index = double(string(indexNode.getFirstChild.getData));
        atlas_xml_struct(cntas).name = char(nameNode.getFirstChild.getData);
        cntas = cntas + 1;
        atlas_xml_struct(cntas) = def_atlas_xml_struct;
    end
    %## GET LABELS
    for cl_i = 1:length(cluster_inds)
        %## AGGREGATE ANATOMY FOR ALL DIPS IN CL
        if chk_a || chk_aa
            dip_in = STUDY.cluster(cluster_inds(cl_i)).all_diplocs;
            STUDY.cluster(cluster_inds(cl_i)).centroid.dipole.posxyz = mean(dip_in);
            atlas = ft_read_atlas(atlas_fpath);
            anatomy_out = 'error';
            cfg              = [];
            cfg.roi        = dip_in;
            cfg.output     = 'single';
            cfg.atlas      = atlas;
            cfg.verbose = 0;
            %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
            cfg.sphere = 3;
            label_i = ft_volumelookup(cfg, atlas);
            if ~isempty(label_i)
                counts = sum([label_i.count],2);
                [val, indx] = max(counts);
                names = label_i(1).name;
                if strcmp(names(indx),'no_label_found')
                    sub_indx = find(counts ~= 0 & counts < val);
                    if ~isempty(sub_indx)
                        anatomy_out = names{sub_indx};
                    end
                else
                    anatomy_out = names{indx};
                end
            end
            %- out
            if ~strcmp(anatomy_out,'error')
                tmp = strsplit(anatomy_out,' ');
                ind = [atlas_xml_struct.index] == double(string(tmp{2}));
                anatomy_out = atlas_xml_struct(ind).name;
            end
            anatomy_struct(cnt).cluster = cluster_inds(cl_i);
            anatomy_struct(cnt).dips = dip_in;
            anatomy_struct(cnt).subjects = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).subject};
            anatomy_struct(cnt).group = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group};
            anatomy_struct(cnt).group_n = length(STUDY.cluster(cluster_inds(cl_i)).sets);
            anatomy_struct(cnt).calculation = 'aggregate label for all';
            anatomy_struct(cnt).anatomy_label = anatomy_out;
            anatomy_struct(cnt).atlas = atlas;
            tmp = strsplit(atlas_fpath,filesep);
            anatomy_struct(cnt).atlas_label = tmp{end};
            cnt = cnt + 1;
            anatomy_struct(cnt) = def_anatomy_struct;
        end
        %## CENTROID ANATOMY FOR MEAN DIP
        if chk_a || chk_ac
            dip_in = STUDY.cluster(cluster_inds(cl_i)).centroid.dipole.posxyz;
            atlas = ft_read_atlas(atlas_fpath);
            atlas_name_ct = 'error';
            cfg              = [];
            cfg.roi        = dip_in;
            cfg.output     = 'multiple';
            cfg.atlas      = atlas;
            cfg.verbose = 0;
            %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
            cfg.sphere = 3;
            label_i = ft_volumelookup(cfg, atlas);
            if ~isempty(label_i)
                counts = sum([label_i.count],2);
                [val, indx] = max(counts);
                names = label_i(1).name;
                if strcmp(names(indx),'no_label_found')
                    sub_indx = find(counts ~= 0 & counts < val);
                    if ~isempty(sub_indx)
                        atlas_name_ct = names{sub_indx};
                    end
                else
                    atlas_name_ct = names{indx};
                end
            end
            %- out
            if ~strcmp(atlas_name_ct,'error')
                tmp = strsplit(atlas_name_ct,' ');
                ind = [atlas_xml_struct.index] == double(string(tmp{2}));
                atlas_name_ct = atlas_xml_struct(ind).name;
            end
            anatomy_struct(cnt).cluster = cluster_inds(cl_i);
            anatomy_struct(cnt).dips = dip_in;
            anatomy_struct(cnt).subjects = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).subject};
            anatomy_struct(cnt).group = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group};
            anatomy_struct(cnt).group_n = length(STUDY.cluster(cluster_inds(cl_i)).sets);
            anatomy_struct(cnt).calculation = 'centroid label for all';
            anatomy_struct(cnt).anatomy_label = atlas_name_ct;
            anatomy_struct(cnt).atlas = atlas;
            tmp = strsplit(atlas_fpath,filesep);
            anatomy_struct(cnt).atlas_label = tmp{end};
            cnt = cnt + 1;
            anatomy_struct(cnt) = def_anatomy_struct;
        end
        %## GROUP ANATOMY FOR AGGREGATE GROUP IN CL
        if chk_a || chk_ga
            atlas_name_gag = cell(1,length(group_chars));
            centroid_gct = cell(1,length(group_chars));
            for g_i = 1:length(group_chars)
                g_inds = cellfun(@(x) strcmp(x,group_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group});                   
                % s_inds = STUDY.cluster(clusters(cl_i)).sets(g_inds);        
                % dip1 = STUDY.cluster(clusters(k_i)).centroid.dipole.posxyz;
                dip_in = STUDY.cluster(cluster_inds(cl_i)).all_diplocs(g_inds,:);
                centroid_gct{g_i} = dip_in;
                atlas = ft_read_atlas(atlas_fpath);
                atlas_name_gag{g_i} = 'error';
                cfg              = [];
                cfg.roi        = dip_in;
                cfg.output     = 'multiple';
                cfg.atlas      = atlas;
                cfg.verbose = 0;
                %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
                cfg.sphere = 3;
                label_i = ft_volumelookup(cfg, atlas);
                if ~isempty(label_i)
                    counts = sum([label_i.count],2);
                    [val, indx] = max(counts);
                    names = label_i(1).name;
                    if strcmp(names(indx),'no_label_found')
                        sub_indx = find(counts ~= 0 & counts < val);
                        if ~isempty(sub_indx)
                            atlas_name_gag{g_i} = names{sub_indx};
                        end
                    else
                        atlas_name_gag{g_i} = names{indx};
                    end
                end
                %- out
                if ~strcmp(atlas_name_gag{g_i},'error')
                    tmp = strsplit(atlas_name_gag{g_i},' ');
                    ind = [atlas_xml_struct.index] == double(string(tmp{2}));
                    atlas_name_gag{g_i} = atlas_xml_struct(ind).name;
                end
                anatomy_struct(cnt).cluster = cluster_inds(cl_i);
                anatomy_struct(cnt).dips = dip_in;
                anatomy_struct(cnt).subjects = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets(g_inds)).subject};
                anatomy_struct(cnt).group = {group_chars{g_i}};
                anatomy_struct(cnt).group_n = sum(g_inds);
                anatomy_struct(cnt).calculation = 'aggregate label for group';
                anatomy_struct(cnt).anatomy_label = atlas_name_gag{g_i};
                anatomy_struct(cnt).atlas = atlas;
                tmp = strsplit(atlas_fpath,filesep);
                anatomy_struct(cnt).atlas_label = tmp{end};
                cnt = cnt + 1;
                anatomy_struct(cnt) = def_anatomy_struct;
            end
        end
        %## GROUP CENTROID ANATOMY FOR GROUP IN CL
        if chk_a || chk_gc
            atlas_name_gct = cell(1,length(group_chars));
            centroid_gct = cell(1,length(group_chars));
            for g_i = 1:length(group_chars)
                g_inds = cellfun(@(x) strcmp(x,group_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group});                
                % s_inds = STUDY.cluster(cluster_inds(cl_i)).sets(g_inds);        
                % dip1 = STUDY.cluster(clusters(k_i)).centroid.dipole.posxyz;
                dip_in = mean(STUDY.cluster(cluster_inds(cl_i)).all_diplocs(g_inds,:),1);
                centroid_gct{g_i} = dip_in;
                atlas = ft_read_atlas(atlas_fpath);
                atlas_name_gct{g_i} = 'error';
                cfg              = [];
                cfg.roi        = dip_in;
                cfg.output     = 'multiple';
                cfg.atlas      = atlas;
                cfg.verbose = 0;
                %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
                cfg.sphere = 3;
                label_i = ft_volumelookup(cfg, atlas);
                if ~isempty(label_i)
                    counts = sum([label_i.count],2);
                    [val, indx] = max(counts);
                    names = label_i(1).name;
                    if strcmp(names(indx),'no_label_found')
                        sub_indx = find(counts ~= 0 & counts < val);
                        if ~isempty(sub_indx)
                            atlas_name_gct{g_i} = names{sub_indx};
                        end
                    else
                        atlas_name_gct{g_i} = names{indx};
                    end
                end
                %- out
                if ~strcmp(atlas_name_gct{g_i},'error')
                    tmp = strsplit(atlas_name_gct{g_i},' ');
                    ind = [atlas_xml_struct.index] == double(string(tmp{2}));
                    atlas_name_gct{g_i} = atlas_xml_struct(ind).name;
                end
                anatomy_struct(cnt).cluster = cluster_inds(cl_i);
                anatomy_struct(cnt).dips = dip_in;
                anatomy_struct(cnt).subjects = {STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets(g_inds)).subject};
                anatomy_struct(cnt).group = {group_chars{g_i}};
                anatomy_struct(cnt).group_n = sum(g_inds);
                anatomy_struct(cnt).calculation = 'centroid label for group';
                anatomy_struct(cnt).anatomy_label = atlas_name_gct{g_i};
                anatomy_struct(cnt).atlas = atlas;
                tmp = strsplit(atlas_fpath,filesep);
                anatomy_struct(cnt).atlas_label = tmp{end};
                cnt = cnt + 1;
                anatomy_struct(cnt) = def_anatomy_struct;
            end
        end
        %## TEXT OUTPUT
        str_ctg = [];
        for g_i = 1:length(group_chars)
            g_inds = cellfun(@(x) strcmp(x,group_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group});
            str_ctg = [str_ctg, sprintf('Group %s N=%i\n',group_chars{g_i},sum(g_inds))];
            %-
            if chk_a || chk_ga || chk_gc
                str_ctg = [str_ctg, sprintf('Group %s Centroid Dip: [%0.1f,%0.1f,%0.1f]\n',group_chars{g_i},centroid_gct{g_i})];
            end
            %-
            if chk_a || chk_ga
                str_ctg = [str_ctg, sprintf('Group %s Aggregate Label: CL%i: %s\n',group_chars{g_i},cluster_inds(cl_i),atlas_name_gag{g_i})];
            end
            %-
            if chk_a || chk_gc
                str_ctg = [str_ctg, sprintf('Group %s Centroid: CL%i: %s\n',group_chars{g_i},cluster_inds(cl_i),atlas_name_gct{g_i})];
            end
        end
        str_cta = [];
        if chk_a || chk_ac || chk_aa
            str_cta = [str_cta,sprintf('\n\nAtlas %s; CL%i: N=%i\n',char(nn_xml.item(0).getFirstChild.getData),cluster_inds(cl_i),length(STUDY.cluster(cluster_inds(cl_i)).sets)),...
                sprintf('All Centroid Dip: [%0.1f,%0.1f,%0.1f]\n',STUDY.cluster(cluster_inds(cl_i)).centroid.dipole.posxyz)];
        end
        if chk_a || chk_ac
             str_cta = [str_cta, sprintf('ALL Centroid Label: %s\n',atlas_name_ct)];
        end
        if chk_a || chk_aa
            str_cta = [str_cta, sprintf('ALL Aggregate Label: %s\n',anatomy_out)];
        end
        txt_store{cnttxt} = [str_cta,...
            str_ctg];
        cnttxt = cnttxt + 1;
    end
end
%## SAVE
if ANATOMY_STRUCT.save_inf
    txt_store = txt_store(~cellfun(@isempty,txt_store));
    cellfun(@(x) fprintf(f,x),txt_store);
    fclose(f);
    save([ANATOMY_STRUCT.save_dir filesep 'anatomy_struct_out.mat'],'anatomy_struct');
end
fprintf('eeglab_get_anatomy.m done: %0.1f\n\n',toc(tt));
end
%% ===================================================================== %%
function chk = chk_field(fcn,struct_in,field_in)
    chk = isfield(struct_in,field_in);
    if chk
        chk = chk && feval(fcn,struct_in.(field_in));
    end
end
