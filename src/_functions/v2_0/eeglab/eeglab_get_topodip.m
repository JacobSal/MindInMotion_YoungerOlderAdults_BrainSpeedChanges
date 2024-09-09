function [STUDY,dipfit_structs,topo_cells] = eeglab_get_topodip(STUDY,varargin)
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
TOPO_PLOTRAD = 0.55;
TOPO_INTRAD = 0.55;
dipfit_structs = [];
topo_cells = [];
%## ANATOMY STRUCT
inds = cellfun(@(x) contains(x,'Outlier','IgnoreCase',true) || contains(x,'Parentcluster'),{STUDY.cluster.name});
DEF_CALC_STRUCT = struct('cluster_inds',find(~inds),...
    'save_inf',true,...
    'recalculate',true,...
    'topo_cells',{{}},...
    'dipfit_structs',struct.empty);
%## ALLEEG DEFAULT
ALLEEG = struct.empty;
%- logo
cat_logo();
%## TIME
tt0 = tic;
%## DEFINE DEFAULTS
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'ALLEEG',ALLEEG,@isstruct)
addParameter(p,'CALC_STRUCT',DEF_CALC_STRUCT,@(x) validate_struct(x,DEF_CALC_STRUCT));
%## PARSE
parse(p, STUDY, varargin{:});
%## SET DEFAULTS
CALC_STRUCT = p.Results.CALC_STRUCT;
ALLEEG = p.Results.ALLEEG;
%-
CALC_STRUCT = set_defaults_struct(CALC_STRUCT,DEF_CALC_STRUCT);
%% PRELOAD ALLEEG SET INFORMATION ====================================== %%
chk = zeros(length(CALC_STRUCT.cluster_inds),7);
for i = 1:length(CALC_STRUCT.cluster_inds)
    cl_i = CALC_STRUCT.cluster_inds(i);
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topox');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topoy');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topoall');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topo');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'topopol');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'all_diplocs');
    chk(i,1) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'dipole');
    % fprintf('%i\n',chk);
    % STUDY.cluster(cl_i).centroid.dipole.posxyz
end
chk = all(chk,[2,1]);
chkk = (isempty(CALC_STRUCT.topo_cells) && isempty(CALC_STRUCT.dipfit_structs)) || CALC_STRUCT.recalculate;
if ~isempty(ALLEEG)
    if ~chk || CALC_STRUCT.recalculate;
        [STUDY,~] = std_centroid(STUDY,ALLEEG,double(string(CALC_STRUCT.cluster_inds)),'dipole');
        tmp_study = STUDY;
        if isfield(tmp_study.cluster,'topox') || isfield(tmp_study.cluster,'topoall') ||...
                isfield(tmp_study.cluster,'topopol') 
            tmp_study.cluster = rmfield(tmp_study.cluster,'topox');
            tmp_study.cluster = rmfield(tmp_study.cluster,'topoy');
            tmp_study.cluster = rmfield(tmp_study.cluster,'topoall');
            tmp_study.cluster = rmfield(tmp_study.cluster,'topo');
            tmp_study.cluster = rmfield(tmp_study.cluster,'topopol');
        end
        if ~isfield(tmp_study.cluster,'topo')
            tmp_study.cluster(1).topo = [];
        end
        for j = 1:length(CALC_STRUCT.cluster_inds) % For each cluster requested
            cl_i = CALC_STRUCT.cluster_inds(j);
            tmp_study = std_readtopoclust_CL(tmp_study,ALLEEG,cl_i);% Using this custom modified code to allow taking average within participant for each cluster
            STUDY.cluster(cl_i).topox = tmp_study.cluster(cl_i).topox;
            STUDY.cluster(cl_i).topoy = tmp_study.cluster(cl_i).topoy;
            STUDY.cluster(cl_i).topoall = tmp_study.cluster(cl_i).topoall;
            STUDY.cluster(cl_i).topo = tmp_study.cluster(cl_i).topo;
            STUDY.cluster(cl_i).topopol = tmp_study.cluster(cl_i).topopol;
        end
    end
elseif ~chk
    if chkk
        % (08/27/2024) JS, Little bit more memory efficient than loading whole
        % ALLEEG structure.
        def_tmp_topo = struct('topotmp',[],...
                'topo',{{}},...
                'topox',[],...
                'topoy',[],...
                'comp',[],...
                'subj',[]);
        dipfit_structs = cell(1,length(STUDY.datasetinfo));
        topo_cells = cell(1,length(STUDY.datasetinfo));
        % centroid_topo = cell(1,length(CALC_STRUCT.cluster_inds));
        tt1 = tic;
        for s_i = 1:length(STUDY.datasetinfo)
            %## LOAD .SET
            fprintf('Loading %s .set file\n',STUDY.datasetinfo(s_i).subject);
            % tt2 = tic;
            tmp = load([STUDY.datasetinfo(s_i).filepath filesep STUDY.datasetinfo(s_i).filename],'-mat',...
                'dipfit','chanlocs','icawinv','icachansind','chaninfo');
            % fprintf('Done loading %s .set file: %0.3g\n',STUDY.datasetinfo(s_i).subject,toc(tt2)/60)
            %## LOAD DIPOLE INFORMATION
            if ~isfield(tmp, 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(s_i) ], 'Aborting compute centroid dipole');
               return;
            end
            dipfit_structs{s_i} = tmp.dipfit;
            %## LOAD TOPO INFORMATION
            chanlocs = tmp.chanlocs(tmp.icachansind);
            tmp_topo = def_tmp_topo;
            %## (08/27/2024) JS, leaving off here for now, need to optimize to make
            %sure its only loading in the component topos that are necessary for
            %the averaging calculation. Would be silly to load in all comp data.
            %UPDATE: using size(icawinv,2) for now.
            for ch_i = 1:size(tmp.icawinv,2)
                [~, grid, ~, Xi, Yi] = topoplot(tmp.icawinv(:,ch_i), chanlocs,...
                    'verbose', 'off',...
                    'electrodes','on','style','both',...
                    'plotrad',TOPO_PLOTRAD,'intrad',TOPO_INTRAD,...
                    'noplot', 'on', 'chaninfo', tmp.chaninfo);
                tmp_topo(ch_i).topo = grid;
                tmp_topo(ch_i).topox = Xi;
                tmp_topo(ch_i).topoy = Yi;
                if isempty(grid)
                    tmp_topo(ch_i).topotmp  = zeros([size(grid(1:4:end),2)]);
                end
                tmp_topo(ch_i).topotmp = grid(1:4:end);
                tmp_topo(ch_i).comp = ch_i;
                tmp_topo(ch_i).subj = s_i;
                if ch_i < size(tmp.icawinv,2)
                    tmp_topo(ch_i+1) = def_tmp_topo;
                end
            end
            topo_cells{s_i} = tmp_topo;
        end
        dipfit_structs = cat(2,dipfit_structs{:});
        if CALC_STRUCT.save_inf
            par_save(dipfit_structs,STUDY.filepath,'dipfit_structs.mat');
            par_save(topo_cells,STUDY.filepath,'topo_cells.mat');
        end
        fprintf('Done loading dipfit structures: %0.2g min\n',toc(tt1)/60);
    else
        fprintf('Using inputted information...\n');
    end
    %% GET TOPO PLOT CENTROIDS ============================================= %%
    % (08/27/2024) CL, Notes added by CL: this chunk code grabs scalp topographs from each
    % subject and each component
    centroid_topo = cell(1,length(CALC_STRUCT.cluster_inds));
    centroid_avg = cell(1,length(CALC_STRUCT.cluster_inds));
    centroid_new = cell(1,length(CALC_STRUCT.cluster_inds));
    for cl_i = 1:length(CALC_STRUCT.cluster_inds) %go over all requested clusters
        chk = chk_field(@isempty,STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)),'topo');
        %-
        if ~chk
            numitems = length(STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).comps);
            for k = 1:numitems % go through all components
                comp  = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).comps(k);
                abset = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).sets(k);
                if ~isnan(comp) && ~isnan(abset)
                    % [grid yi xi] = std_readtopo(ALLEEG, abset, comp);
                    
                    tmp = topo_cells{abset};
                    centroid_topo{cl_i}.topotmp(:,k) = tmp(comp(k)).topotmp; %grid(1:4:end); % for inversion
                    centroid_topo{cl_i}.topo{k} = tmp(comp(k)).topo; %grid;
                    centroid_topo{cl_i}.topox = tmp(comp(k)).topox; %xi;
                    centroid_topo{cl_i}.topoy = tmp(comp(k)).topoy; %yi;
                end
            end
            fprintf('\n');
            %## UPDATE STUDY
            [~, pol] = std_comppol(centroid_topo{cl_i}.topotmp);
            fprintf('%d/%d polarities inverted while reading component scalp maps\n', ...
                    length(find(pol == -1)), length(pol));
            % centroid_orig = centroid_topo;%make a copy just in case
            subj_in_cluster = unique(STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).sets);%Subjects in this cluster
            for j = 1:length(unique(STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).sets)) 
                comp_ind = find(STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).sets == subj_in_cluster(j));
                nICs = length(comp_ind);
                for k = 1:nICs
                    centroid_new{cl_i}.topo{comp_ind(k)} =  pol(comp_ind(k))*centroid_topo{cl_i}.topo{comp_ind(k)};
                    if k == 1
                        avgScalp =  centroid_new{cl_i}.topo{comp_ind(k)} /nICs;
                        centroid_avg{cl_i}.topo{j} = avgScalp;
                    else
                        avgScalp =  centroid_new{cl_i}.topo{comp_ind(k)} /nICs + avgScalp;
                        centroid_avg{cl_i}.topo{j} = avgScalp;
                    end
                end
            end
            % Then take average across all participants
            nitems = length(subj_in_cluster);
            for k = 1:nitems
                if k == 1
                    allscalp = centroid_avg{cl_i}.topo{k}/nitems;
                else
                    allscalp = centroid_avg{cl_i}.topo{k}/nitems + allscalp;
                end
            end
            STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topox   = centroid_topo{cl_i}.topox;
            STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topoy   = centroid_topo{cl_i}.topoy;
            STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topoall = centroid_avg{cl_i}.topo;
            STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topo    = allscalp;
            STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topopol = pol;
        else
            centroid_topo{cl_i}.topox = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topox;
            centroid_topo{cl_i}.topoy = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topoy;
            centroid_topo{cl_i}.topo  = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).topoall;
        end
    end
    %% CENTROID CALCULATION ================================================ %%
    centroid_dip = cell(1,length(STUDY.cluster));
    chk = zeros(length(STUDY.cluster));
    for i = 1:length(STUDY.cluster)
        cl_i = CALC_STRUCT.cluster_inds(i);
        chk(i) = chk_field((@(x) ~isempty(x)),STUDY.cluster(cl_i),'centroid');
    end
    chk = all(chk);
    if ~chk
        error('All cluster centroids calculated, no need to continue.\n');
    end
    %##
    for cl_i = 1:length(CALC_STRUCT.cluster_inds)
        max_r = 0;
        len = length(STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).comps);
        tmppos = 0;
        tmpmom = 0;
        tmprv = 0;
        ndip = 0;
        for k = 1:len 
            fprintf('.');
            comp  = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).comps(k);
            abset = STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).sets(1,k);
            if ~isempty(dipfit_structs(abset).model(comp).posxyz)
                ndip = ndip +1;
                tmppos = tmppos + dipfit_structs(abset).model(comp).posxyz;
                tmpmom = tmpmom + dipfit_structs(abset).model(comp).momxyz;
                tmprv = tmprv + dipfit_structs(abset).model(comp).rv;
                if strcmpi(dipfit_structs(abset).coordformat, 'spherical')
                   if isfield(dipfit_structs(abset), 'hdmfile') %dipfit 2 spherical model
                       load('-mat', dipfit_structs(abset).hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(dipfit_structs(abset).vol.r));
                   end
               end
            end
        end
        centroid_dip{cl_i}.dipole.posxyz =  tmppos/ndip;
        centroid_dip{cl_i}.dipole.momxyz =  tmpmom/ndip;
        centroid_dip{cl_i}.dipole.rv =  tmprv/ndip;
        if strcmpi(dipfit_structs(abset).coordformat, 'spherical') && (~isfield(dipfit_structs(abset), 'hdmfile')) %old dipfit
            centroid_dip{cl_i}.dipole.maxr = max_r;
        end
        STUDY.cluster(CALC_STRUCT.cluster_inds(cl_i)).centroid.dipole = centroid_dip{cl_i}.dipole;
    end
else
    fprintf('Topo & Dipole Information Calculated... Change CALC_STRUCT.recalculate if you want to recalculate information.\n');
end
fprintf('Done eeglab_get_topodip.m: %0.3g',toc(tt0));
end
%% ===================================================================== %%
function chk = chk_field(fcn,struct_in,field_in)
    chk = isfield(struct_in,field_in);
    if chk
        chk = chk && feval(fcn,struct_in.(field_in));
    end
end

