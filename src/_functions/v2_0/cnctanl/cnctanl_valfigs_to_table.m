function [tbl_out,tbl_summary_out] = cnctanl_valfigs_to_table(fig_dir,cond_order,varargin)
%FIGS_TO_EXCEL Summary of this function goes here
%   Detailed explanation goes here
ORDER_CHAR = 'order';
VALIDATION_CHAR = 'validation';
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
        %% Line Graph
        line_X =     file.hgS_070000.children(2).children(1).properties.XData;
        aic_line_Y = file.hgS_070000.children(2).children(1).properties.YData;
        hq_line_Y  = file.hgS_070000.children(2).children(2).properties.YData;
        %% Histograms
        hist_X     = file.hgS_070000.children(4).children(1).properties.XData;
        aic_hist_Y = file.hgS_070000.children(4).children(1).properties.YData;
        hq_hist_Y  = file.hgS_070000.children(5).children(1).properties.YData;
%         header = ["Model Order", "AIC_BITS", "HQ_BITS", "AIC_HIST", "HQ_HIST"];
        info_crit_modord_x{i} = squeeze(line_X).';
        info_crit_aic_line{i} = squeeze(aic_line_Y).';
        info_crit_hq_line{i} = squeeze(hq_line_Y).';
        hist_counts_x{i} = squeeze(hist_X).';
        hist_aic_amnts{i} = squeeze(aic_hist_Y).';
        hist_hq_amnts{i} = squeeze(hq_hist_Y).';
    elseif contains(fdir(i).name,VALIDATION_CHAR)
        window_num{i} = file.hgS_070000.children(2).children(1).properties.XData;
        white_sig_ljb{i} = file.hgS_070000.children(2).children(1).properties.YData;
        white_sig_acf{i} = file.hgS_070000.children(2).children(2).properties.YData;
        white_sig_boxp{i} = file.hgS_070000.children(2).children(3).properties.YData;
        white_sig_limcl{i} = file.hgS_070000.children(2).children(4).properties.YData;
        perc_cons{i} = file.hgS_070000.children(4).children(1).properties.YData;
        stability_ind{i} = file.hgS_070000.children(5).children(1).properties.YData;
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