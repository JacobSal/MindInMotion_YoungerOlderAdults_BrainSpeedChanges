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
    'save_inf',true,...
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
atlas_fpath = ANATOMY_STRUCT.atlas_fpath;
group_chars = ANATOMY_STRUCT.group_chars;
%% PRELOAD ALLEEG SET INFORMATION ====================================== %%
CALC_STRUCT = struct('cluster_inds',ANATOMY_STRUCT.cluster_inds,...
    'save_inf',ANATOMY_STRUCT.save_inf,...
    'recalculate',false);
[STUDY,dipfit_structs,topo_cells] = eeglab_get_topodip(STUDY,...
    'CALC_STRUCT',CALC_STRUCT,...
    'ALLEEG',ALLEEG);
% chk = zeros(length(CALC_STRUCT.cluster_inds),7);
% for i = 1:length(CALC_STRUCT.cluster_inds)
%     cl_i = CALC_STRUCT.cluster_inds(i);
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topox');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topoy');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topoall');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topo');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topopol');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'all_diplocs');
%     chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'dipole');
%     % fprintf('%i\n',chk);
%     % STUDY.cluster(cl_i).centroid.dipole.posxyz
% end
% chk = all(chk,[2,1]);
% if ~isempty(ALLEEG) && ~chk
%     [STUDY,~] = std_centroid(STUDY,ALLEEG,double(string(cluster_inds)),'dipole');
%     tmp_study = STUDY;
%     if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') ||...
%             isfield(tmp_study.cluster,'topopol') 
%         tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
%         tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
%         tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
%         tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
%         tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
%     end
%     if ~isfield(tmp_study.cluster,'topo')
%         tmp_study.cluster(1).topo = [];
%     end
%     for j = 1:length(cluster_inds) % For each cluster requested
%         cl_i = cluster_inds(j);
%         tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
%         STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
%         STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
%         STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
%         STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
%         STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
%     end
% elseif ~chk
%     %%
%     % (08/27/2024) JS, Little bit more memory efficient than loading whole
%     % ALLEEG structure.
%     if isempty(ANATOMY_STRUCT.topo_cells) && isempty(ANATOMY_STRUCT.dipfit_structs)
%         def_tmp_topo = struct('topotmp',[],...
%                 'topo',{{}},...
%                 'topox',[],...
%                 'topoy',[],...
%                 'comp',[],...
%                 'subj',[]);
%         dipfit_structs = cell(1,length(STUDY.datasetinfo));
%         topo_cells = cell(1,length(STUDY.datasetinfo));
%         % centroid_topo = cell(1,length(CLUSTER_INDS));
%         tt = tic;
%         for s_i = 1:length(STUDY.datasetinfo)
%             %## LOAD .SET
%             fprintf('Loading %s .set file\n',STUDY.datasetinfo(s_i).subject);
%             tmp = load([STUDY.datasetinfo(s_i).filepath filesep STUDY.datasetinfo(s_i).filename],'-mat');
%             %## LOAD DIPOLE INFORMATION
%             if ~isfield(tmp, 'dipfit')
%                warndlg2(['No dipole information available in dataset ' num2str(s_i) ], 'Aborting compute centroid dipole');
%                return;
%             end
%             dipfit_structs{s_i} = tmp.dipfit;
%             %## LOAD TOPO INFORMATION
%             chanlocs = tmp.chanlocs(tmp.icachansind);
%             tmp_topo = def_tmp_topo;
%             %## (08/27/2024) JS, leaving off here for now, need to optimize to make
%             %sure its only loading in the component topos that are necessary for
%             %the averaging calculation. Would be silly to load in all comp data.
%             %UPDATE: using size(icawinv,2) for now.
%             for ch_i = 1:size(tmp.icawinv,2)
%                 [~, grid, ~, Xi, Yi] = topoplot(tmp.icawinv(:,ch_i), chanlocs,...
%                     'verbose', 'off',...
%                     'electrodes', 'on' ,'style','both',...
%                     'plotrad',0.55,'intrad',0.55,...
%                     'noplot', 'on', 'chaninfo', tmp.chaninfo);
%                 tmp_topo(ch_i).topo = grid;
%                 tmp_topo(ch_i).topox = Xi;
%                 tmp_topo(ch_i).topoy = Yi;
%                 tmp_topo(ch_i).topotmp = grid(1:4:end);
%                 tmp_topo(ch_i).comp = ch_i;
%                 tmp_topo(ch_i).subj = s_i;
%                 if ch_i < size(tmp.icawinv,2)
%                     tmp_topo(ch_i+1) = def_tmp_topo;
%                 end
%             end
%             topo_cells{s_i} = tmp_topo;
%         end
%         dipfit_structs = cat(2,dipfit_structs{:});
%         if ANATOMY_STRUCT.save_inf
%             par_save(dipfit_structs,STUDY.filepath,'dipfit_structs.mat');
%             par_save(topo_cells,STUDY.filepath,'topo_cells.mat');
%         end
%         fprintf('Done loading dipfit structures: %0.2gs\n',toc(tt));
%     end
%     %% GET TOPO PLOT CENTROIDS ============================================= %%
%     % (08/27/2024) CL, Notes added by CL: this chunk code grabs scalp topographs from each
%     % subject and each component
%     centroid_topo = cell(1,length(cluster_inds));
%     centroid_avg = cell(1,length(cluster_inds));
%     centroid_new = cell(1,length(cluster_inds));
%     for cl_i = 1:length(cluster_inds) %go over all requested clusters
%         if isempty( STUDY.cluster(cluster_inds(cl_i)).topo)
%             numitems = length(STUDY.cluster(cluster_inds(cl_i)).comps);
%             for k = 1:numitems % go through all components
%                 comp  = STUDY.cluster(cluster_inds(cl_i)).comps(k);
%                 abset = STUDY.cluster(cluster_inds(cl_i)).sets(cond,k);
%                 if ~isnan(comp) && ~isnan(abset)
%                     % [grid yi xi] = std_readtopo(ALLEEG, abset, comp);
%                     if ~isfield(centroid_topo{cl_i}, 'topotmp') || isempty(centroid_topo{cl_i}.topotmp)
%                         centroid_topo{cl_i}.topotmp = zeros([ size(grid(1:4:end),2) numitems ]);
%                     end
%                     tmp = topo_cells{abset};
%                     centroid_topo{cl_i}.topotmp(:,k) = tmp(comp(k)).topotmp; %grid(1:4:end); % for inversion
%                     centroid_topo{cl_i}.topo{k} = tmp(comp(k)).topo; %grid;
%                     centroid_topo{cl_i}.topox = tmp(comp(k)).topox; %xi;
%                     centroid_topo{cl_i}.topoy = tmp(comp(k)).topoy; %yi;
%                 end
%             end
%             fprintf('\n');
%             %update STUDY
%             % tmpinds = find(isnan(centroid_topo{cl_i}.topotmp(:,1)));
%             for cond  = 1
%                 % if CLUSTER_INDS(1) > 0
%                 %     ncomp = length(STUDY.cluster(CLUSTER_INDS(cl_i)).comps);
%                 % end
%                 [tmp, pol] = std_comppol(centroid_topo{cl_i}.topotmp);
%                 fprintf('%d/%d polarities inverted while reading component scalp maps\n', ...
%                         length(find(pol == -1)), length(pol));
%                 % centroid_orig = centroid_topo;%make a copy just in case
%                 subj_in_cluster = unique(STUDY.cluster(cluster_inds(cl_i)).sets);%Subjects in this cluster
%                 for j = 1:length(unique(STUDY.cluster(cluster_inds(cl_i)).sets)) 
%                     comp_ind = find(STUDY.cluster(cluster_inds(cl_i)).sets == subj_in_cluster(j));
%                     nICs = length(comp_ind);
%                     for k = 1:nICs
%                         centroid_new{cl_i}.topo{comp_ind(k)} =  pol(comp_ind(k))*centroid_topo{cl_i}.topo{comp_ind(k)};
%                         if k == 1
%                             avgScalp =  centroid_new{cl_i}.topo{comp_ind(k)} /nICs;
%                             centroid_avg{cl_i}.topo{j} = avgScalp;
%                         else
%                             avgScalp =  centroid_new{cl_i}.topo{comp_ind(k)} /nICs + avgScalp;
%                             centroid_avg{cl_i}.topo{j} = avgScalp;
%                         end
%                     end
%                 end
%                 % Then take average across all participants
%                 nitems = length(subj_in_cluster);
%                 for k = 1:nitems
%                     if k == 1
%                         allscalp = centroid_avg{cl_i}.topo{k}/nitems;
%                     else
%                         allscalp = centroid_avg{cl_i}.topo{k}/nitems + allscalp;
%                     end
%                 end
%                 STUDY.cluster(cluster_inds(cl_i)).topox   = centroid_topo{cl_i}.topox;
%                 STUDY.cluster(cluster_inds(cl_i)).topoy   = centroid_topo{cl_i}.topoy;
%                 STUDY.cluster(cluster_inds(cl_i)).topoall = centroid_avg{cl_i}.topo;
%                 STUDY.cluster(cluster_inds(cl_i)).topo    = allscalp;
%                 STUDY.cluster(cluster_inds(cl_i)).topopol = pol;
%             end
%         else
%             centroid_topo{cl_i}.topox = STUDY.cluster(cluster_inds(cl_i)).topox;
%             centroid_topo{cl_i}.topoy = STUDY.cluster(cluster_inds(cl_i)).topoy;
%             centroid_topo{cl_i}.topo  = STUDY.cluster(cluster_inds(cl_i)).topoall;
%         end
%     end
%     %% CENTROID CALCULATION ================================================ %%
%     centroid_dip = cell(1,length(STUDY.cluster));
%     chk = zeros(length(STUDY.cluster),2);
%     for i = 1:length(STUDY.cluster)
%         cl_i = cluster_inds(i);
%         chk(i,1) = isfield(STUDY.cluster(cl_i),'centroid');
%         chk(i,2) = isempty(STUDY.cluster(cntcl_i).centroid);
%     end
%     if ~chk
%         error('All cluster centroids calculated, no need to continue.\n');
%     end
%     %##
%     for cl_i = 1:length(cluster_inds)
%         max_r = 0;
%         len = length(STUDY.cluster(cluster_inds(cl_i)).comps);
%         tmppos = 0;
%         tmpmom = 0;
%         tmprv = 0;
%         ndip = 0;
%         for k = 1:len 
%             fprintf('.');
%             comp  = STUDY.cluster(cluster_inds(cl_i)).comps(k);
%             abset = STUDY.cluster(cluster_inds(cl_i)).sets(1,k);
%             if ~isempty(dipfit_structs(abset).model(comp).posxyz)
%                 ndip = ndip +1;
%                 tmppos = tmppos + dipfit_structs(abset).model(comp).posxyz;
%                 tmpmom = tmpmom + dipfit_structs(abset).model(comp).momxyz;
%                 tmprv = tmprv + dipfit_structs(abset).model(comp).rv;
%                 if strcmpi(dipfit_structs(abset).coordformat, 'spherical')
%                    if isfield(dipfit_structs(abset), 'hdmfile') %dipfit 2 spherical model
%                        load('-mat', dipfit_structs(abset).hdmfile);
%                        max_r = max(max_r, max(vol.r));
%                    else % old version of dipfit
%                        max_r = max(max_r,max(dipfit_structs(abset).vol.r));
%                    end
%                end
%             end
%         end
%         centroid_dip{cl_i}.dipole.posxyz =  tmppos/ndip;
%         centroid_dip{cl_i}.dipole.momxyz =  tmpmom/ndip;
%         centroid_dip{cl_i}.dipole.rv =  tmprv/ndip;
%         if strcmpi(dipfit_structs(abset).coordformat, 'spherical') && (~isfield(dipfit_structs(abset), 'hdmfile')) %old dipfit
%             centroid_dip{cl_i}.dipole.maxr = max_r;
%         end
%         STUDY.cluster(cluster_inds(cl_i)).centroid.dipole = centroid_dip{cl_i}.dipole;
%     end
% else
%     fprintf('Information Calculated... Proceeding to Anatomy Calculation\n');
% end
%% ATLAS CALCULATION =================================================== %%
def_anatomy_struct = struct('cluster',[],...
    'dips',[],...
    'subjects',{''},...
    'group',{''},...
    'calculation',{''},...
    'anatomy_label',{''},...
    'atlas',{''},...
    'atlas_label',{''});
anatomy_struct = def_anatomy_struct;
cnt = 1;
atlas_name_store = cell(length(cluster_inds),1);
txt_store = cell(length(cluster_inds),1);
if ANATOMY_STRUCT.save_inf
    f = fopen([STUDY.filepath filesep 'anatomy_output.txt'],'w');
end
cnttxt = 1;
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
        anatomy_struct(cnt).calculation = 'aggregate label for all';
        anatomy_struct(cnt).anatomy_label = anatomy_out;
        anatomy_struct(cnt).atlas = atlas;
        tmp = strsplit(atlas_fpath,filesep);
        anatomy_struct(cnt).atlas_label = tmp{end};
        cnt = cnt + 1;
        anatomy_struct(cnt) = def_anatomy_struct;
        %## CENTROID ANATOMY FOR MEAN DIP
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
        anatomy_struct(cnt).calculation = 'centroid label for all';
        anatomy_struct(cnt).anatomy_label = atlas_name_ct;
        anatomy_struct(cnt).atlas = atlas;
        tmp = strsplit(atlas_fpath,filesep);
        anatomy_struct(cnt).atlas_label = tmp{end};
        cnt = cnt + 1;
        anatomy_struct(cnt) = def_anatomy_struct;
        %## GROUP ANATOMY FOR AGGREGATE GROUP IN CL
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
            anatomy_struct(cnt).calculation = 'aggregate label for group';
            anatomy_struct(cnt).anatomy_label = atlas_name_gag{g_i};
            anatomy_struct(cnt).atlas = atlas;
            tmp = strsplit(atlas_fpath,filesep);
            anatomy_struct(cnt).atlas_label = tmp{end};
            cnt = cnt + 1;
            anatomy_struct(cnt) = def_anatomy_struct;
        end
        
        
        %## GROUP CENTROID ANATOMY FOR GROUP IN CL
        atlas_name_gct = cell(1,length(group_chars));
        centroid_gct = cell(1,length(group_chars));
        for g_i = 1:length(group_chars)
            g_inds = cellfun(@(x) strcmp(x,group_chars{g_i}),{STUDY.datasetinfo(STUDY.cluster(cluster_inds(cl_i)).sets).group});                   
            s_inds = STUDY.cluster(cluster_inds(cl_i)).sets(g_inds);        
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
            anatomy_struct(cnt).calculation = 'centroid label for group';
            anatomy_struct(cnt).anatomy_label = atlas_name_gct{g_i};
            anatomy_struct(cnt).atlas = atlas;
            tmp = strsplit(atlas_fpath,filesep);
            anatomy_struct(cnt).atlas_label = tmp{end};
            cnt = cnt + 1;
            anatomy_struct(cnt) = def_anatomy_struct;
        end
        %- create group string
        str_ct = [];
        for g_i = 1:length(group_chars)
            str_ct = [str_ct, sprintf('Group %s Centroid: CL%i: %s\n',group_chars{g_i},cluster_inds(cl_i),atlas_name_gct{g_i}),...
                        sprintf('Group %s Centroid Dip: [%0.1f,%0.1f,%0.1f]\n',group_chars{g_i},centroid_gct{g_i}),...
                        sprintf('Group %s Aggregate Label: CL%i: %s\n\n',group_chars{g_i},cluster_inds(cl_i),atlas_name_gag{g_i})];
        end
        txt_store{cnttxt} = [sprintf('\nAtlas %s; CL%i: N=%i\n',char(nn_xml.item(0).getFirstChild.getData),cluster_inds(cl_i),length(STUDY.cluster(cluster_inds(cl_i)).sets)),...
        sprintf('ALL Aggregate Label: %s\n',anatomy_out),...
        sprintf('ALL Centroid Label: %s\n',atlas_name_ct),...
        sprintf('All Centroid Dip: [%0.1f,%0.1f,%0.1f]\n',STUDY.cluster(cluster_inds(cl_i)).centroid.dipole.posxyz),...
        str_ct];
        % atlas_name_store{cl_i} = sprintf('CL%i: %s\n',cluster_inds(cl_i),anatomy_out);
        cnttxt = cnttxt + 1;
    end
end
%## SAVE
if ANATOMY_STRUCT.save_inf
    txt_store = txt_store(~cellfun(@isempty,txt_store));
    cellfun(@(x) fprintf(f,x),txt_store);
    fclose(f);
    save([STUDY.filepath filesep 'anatomy_struct_out.mat'],'anatomy_struct');
end
fprintf('eeglab_get_anatomy.m done: %0.1f\n',toc(tt));
end
%% ===================================================================== %%
function chk = chk_field(fcn,struct_in,field_in)
    chk = isfield(struct_in,field_in);
    if chk
        chk = chk && feval(fcn,struct_in.(field_in));
    end
end
