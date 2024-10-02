function [tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(fig_dir,cond_order,varargin)
%FIGS_TO_EXCEL Summary of this function goes here
%   Detailed explanation goes here
ORDER_CHAR = 'order';
VALIDATION_CHAR = 'validation';
VALID_WHITE_CHARS = {'Ljung-Box','ACF','Box-Pierce','Li-Mcleod'};
% HNDL_IND_ORDER_LINE = 1;
% HNDL_IND_ORDER_HIST = 3;
p = inputParser;
%## REQUIRED
addRequired(p,'fig_dir',@ischar)
addRequired(p,'cond_order',@iscell)
%## PARAMETER
addParameter(p,'ORDER_CHAR',ORDER_CHAR,@ischar);
addParameter(p,'VALIDATION_CHAR',VALIDATION_CHAR,@ischar);
parse(p,fig_dir,cond_order,varargin{:});
%-
ORDER_CHAR = p.Results.ORDER_CHAR;
VALIDATION_CHAR = p.Results.VALIDATION_CHAR;

%%
fdir = dir([fig_dir filesep '*.fig']);
%- define cells for table construction
subj_name = cell(length(fdir),1);
cond_name = cell(length(fdir),1);
window_num = cell(length(fdir),1);
white_sig_ljb = cell(length(fdir),1);
white_sig_acf = cell(length(fdir),1);
white_sig_boxp = cell(length(fdir),1);
white_sig_limcl = cell(length(fdir),1);
perc_cons = cell(length(fdir),1);
stability_ind = cell(length(fdir),1);
info_crit_modord_x = cell(length(fdir),1);
info_crit_hq_line = cell(length(fdir),1);
info_crit_aic_line = cell(length(fdir),1);
hist_counts_x = cell(length(fdir),1);
hist_aic_amnts = cell(length(fdir),1);
hist_hq_amnts = cell(length(fdir),1);
for i = 1:length(fdir)
    file = load([fdir(i).folder filesep fdir(i).name], '-mat');
    fig = get(groot,'CurrentFigure');
    close(fig)
    cond_n = feval(@cond_subj_regexp,fdir(i).name);
    cond_name{i} = cond_order{cond_n};
    subj_name(i) = feval(@subj_name_regexp,fdir(i).name);
    if contains(fdir(i).name,ORDER_CHAR)
        for ch_i = 1:length(file.hgS_070000.children)
            if ~isempty(file.hgS_070000.children(ch_i).children)
                chk1 = any(contains({file.hgS_070000.children(ch_i).children.type},'graph2d.lineseries'));
                chk2 = any(contains({file.hgS_070000.children(ch_i).children.type},'specgraph.baseline'));
                chk3 = any(contains({file.hgS_070000.children(ch_i).children.type},'patch'));
                if chk2 && chk3 && chk1
                    %## histogram
                    for ch_j = 1:length(file.hgS_070000.children(ch_i).children)
                        chk = ~cellfun(@isempty,{file.hgS_070000.children(ch_i).children.children});
                        try
                            if contains(file.hgS_070000.children(ch_i).children(ch_j).properties.String,'aic')
                                aic_hist_Y = file.hgS_070000.children(ch_i).children(chk).properties.YData;
                                hist_X     = file.hgS_070000.children(ch_i).children(chk).properties.XData;
                            elseif contains(file.hgS_070000.children(ch_i).children(ch_j).properties.String,'hq')
                                hq_hist_Y  = file.hgS_070000.children(ch_i).children(chk).properties.YData;
                            end
                        catch
                            % fprintf('Didn''t find valid child\n')
                        end
                    end
                elseif chk1 && ~chk2 && ~chk3
                    %## Line Graph
                    for ch_j = 1:length(file.hgS_070000.children(ch_i).children)
                        try
                            if contains(file.hgS_070000.children(ch_i).children(ch_j).properties.DisplayName,'aic')
                                aic_line_Y = file.hgS_070000.children(ch_i).children(ch_j).properties.YData;
                                %- XData only needs extracting once
                                line_X =     file.hgS_070000.children(ch_i).children(ch_j).properties.XData;
                            elseif contains(file.hgS_070000.children(ch_i).children(ch_j).properties.DisplayName,'hq')
                                hq_line_Y  = file.hgS_070000.children(ch_i).children(ch_j).properties.YData;
                            end
                        catch
                            % fprintf('Didn''t find valid child\n')
                        end
                    end
                else 
                    warning('No extractable handle for order plot');
                end
            end
        end
        info_crit_modord_x{i} = squeeze(line_X).';
        info_crit_aic_line{i} = squeeze(aic_line_Y).';
        info_crit_hq_line{i} = squeeze(hq_line_Y).';
        hist_counts_x{i} = squeeze(hist_X).';
        hist_aic_amnts{i} = squeeze(aic_hist_Y).';
        hist_hq_amnts{i} = squeeze(hq_hist_Y).';
    elseif contains(fdir(i).name,VALIDATION_CHAR)
        for ch_i = 1:length(file.hgS_070000.children)
            try
                chk1 = cellfun(@(x) find_hndl_prop(x,'String','Whiteness Significance'),{file.hgS_070000.children(ch_i).children.properties});
                if any(chk1)
                    chk_ljb = cellfun(@(x) find_hndl_prop(x,'DisplayName','Ljung-Box'),{file.hgS_070000.children(ch_i).children.properties});
                    chk_acf = cellfun(@(x) find_hndl_prop(x,'DisplayName','ACF'),{file.hgS_070000.children(ch_i).children.properties});
                    chk_boxp = cellfun(@(x) find_hndl_prop(x,'DisplayName','Box-Pierce'),{file.hgS_070000.children(ch_i).children.properties});
                    chk_limcl = cellfun(@(x) find_hndl_prop(x,'DisplayName','Li-McLeod'),{file.hgS_070000.children(ch_i).children.properties});
                    white_sig_ljb{i} = file.hgS_070000.children(ch_i).children(chk_ljb).properties.YData;
                    white_sig_acf{i} = file.hgS_070000.children(ch_i).children(chk_acf).properties.YData;
                    white_sig_boxp{i} = file.hgS_070000.children(ch_i).children(chk_boxp).properties.YData;
                    white_sig_limcl{i} = file.hgS_070000.children(ch_i).children(chk_limcl).properties.YData;
                    %- Extract XData Once
                    window_num{i} = file.hgS_070000.children(ch_i).children(chk_ljb).properties.XData;
                end
                chk1 = cellfun(@(x) find_hndl_prop(x,'String','Percent Consistency'),{file.hgS_070000.children(ch_i).children.properties});
                if any(chk1)
                    chk2 = cellfun(@(x) contains(x,'graph2d.lineseries'),{file.hgS_070000.children(ch_i).children.type});
                    perc_cons{i} = file.hgS_070000.children(ch_i).children(chk2).properties.YData;
                end
                chk1 = cellfun(@(x) find_hndl_prop(x,'String','Stability Index'),{file.hgS_070000.children(ch_i).children.properties});
                if any(chk1)
                    % chk2 = cellfun(@(x) contains(x,'graph2d.lineseries'),{file.hgS_070000.children(ch_i).children.type});
                    chk3 = cellfun(@(x) find_hndl_prop(x,'Marker','.'),{file.hgS_070000.children(ch_i).children.properties});
                    perc_cons{i} = file.hgS_070000.children(ch_i).children(chk3).properties.YData;
                end
            catch
                warning('No extractable handle for validation plot');
            end
        end
    else
        fprintf('file %s is not an order or validation figure...\n',fdir(i).name)
    end
end
%%
tbl_out = table(subj_name,cond_name,info_crit_modord_x,info_crit_aic_line,info_crit_hq_line,...
    hist_counts_x,hist_aic_amnts,hist_hq_amnts,window_num,white_sig_ljb,white_sig_acf,white_sig_boxp,...
    white_sig_limcl,perc_cons,stability_ind);
%-
subj_names = [tbl_out.subj_name(:)];
usubj = unique(subj_names);
cond_names = [tbl_out.cond_name(:)];
ucond = unique(cond_names);
nn = length(usubj)*length(ucond);
%-
subj_name = cell(nn,1);
cond_name = cell(nn,1);
mean_white_sig_ljb = zeros(nn,1);
mean_white_sig_acf = zeros(nn,1);
mean_white_sig_boxp = zeros(nn,1);
mean_white_sig_limcl = zeros(nn,1);
mean_perc_cons = zeros(nn,1);
mean_stability_ind = zeros(nn,1);
mean_info_crit_hq_line = zeros(nn,1);
min_info_crit_hq_line = zeros(nn,1);
min_modorder_info_crit_hq_line = zeros(nn,1);
mean_info_crit_aic_line = zeros(nn,1);
min_info_crit_aic_line = zeros(nn,1);
min_modorder_info_crit_aic_line = zeros(nn,1);
mean_hist_aic_amnts = zeros(nn,1);
mean_hist_hq_amnts = zeros(nn,1);
%-
cnt = 1;
for subj_i = 1:length(usubj)
    indss = strcmp(subj_names,usubj{subj_i});
    for cond_i = 1:length(ucond)
        indsc = strcmp(cond_names,ucond{cond_i});
        inds = find(indss & indsc);
        subj_name{cnt} = usubj{subj_i};
        cond_name{cnt} = ucond{cond_i};
        mean_white_sig_ljb(cnt) = mean([tbl_out.white_sig_ljb{inds}]);
        mean_white_sig_acf(cnt) = mean([tbl_out.white_sig_acf{inds}]);
        mean_white_sig_boxp(cnt) = mean([tbl_out.white_sig_boxp{inds}]);
        mean_white_sig_limcl(cnt) = mean([tbl_out.white_sig_limcl{inds}]);
        mean_perc_cons(cnt) = mean([tbl_out.perc_cons{inds}]);
        mean_stability_ind(cnt) = mean([tbl_out.stability_ind{inds}]);
        %-
        xx = [tbl_out.info_crit_modord_x{inds}];
        mean_info_crit_hq_line(cnt) = mean([tbl_out.info_crit_hq_line{inds}]);
        min_info_crit_hq_line(cnt) = min([tbl_out.info_crit_hq_line{inds}]);
        mean_info_crit_aic_line(cnt) = mean([tbl_out.info_crit_aic_line{inds}]);
        min_info_crit_aic_line(cnt) = min([tbl_out.info_crit_aic_line{inds}]);
        i1 = [tbl_out.info_crit_hq_line{inds}] == min_info_crit_hq_line(cnt);
        i2 = [tbl_out.info_crit_aic_line{inds}] == min_info_crit_aic_line(cnt);
        min_modorder_info_crit_hq_line(cnt) = xx(i1);
        min_modorder_info_crit_aic_line(cnt) = xx(i2);
        %-
        tmp1 = [];
        tmp2 = [];
        xx = [tbl_out.hist_counts_x{inds}];
        yy1 = [tbl_out.hist_aic_amnts{inds}];
        yy2 = [tbl_out.hist_hq_amnts{inds}];
        for i = 1:length([tbl_out.hist_counts_x{inds}])
            tmp1 = [tmp1; repmat(xx(i),yy1(i),1)];
            tmp2 = [tmp2; repmat(xx(i),yy2(i),1)];
        end
        mean_hist_aic_amnts(cnt) = mean(tmp1);
        mean_hist_hq_amnts(cnt) = mean(tmp2);
        cnt = cnt + 1;
    end
end
tbl_summary_out = table(mean_white_sig_ljb, mean_white_sig_acf, mean_white_sig_boxp,...
    mean_white_sig_limcl, mean_perc_cons, mean_stability_ind, mean_info_crit_hq_line,...
    min_info_crit_hq_line, min_modorder_info_crit_hq_line, mean_info_crit_aic_line, ...
    min_info_crit_aic_line, min_modorder_info_crit_aic_line, mean_hist_aic_amnts,...
    mean_hist_hq_amnts);
end
%% ===================================================================== %%
function [cond_n] = cond_subj_regexp(str)
    cond_n = regexp(str,'_(\d*)_','match');
    cond_n = strsplit(cond_n{1},'_');
    cond_n = str2double(cond_n{~cellfun(@isempty,cond_n)});
end
%##
function [subj_name] = subj_name_regexp(str)
    % subj_name = regexp(str,'Pilot[\d+]*','match');
    subj_name = regexp(str,'H[\d+]*','match');
end
%##
function [bool] = find_hndl_prop(obj,prop_name,prop_val)
    try
        bool = any(contains(obj.(prop_name),prop_val));
        % if length(obj.(prop_name))
        %     bool = any(bool);
        % end
    catch
        bool = false;
    end
end
