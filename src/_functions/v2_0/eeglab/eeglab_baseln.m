function [allersp_full,allersp_crop,baseln_com,baseln_trial] = eeglab_baseln(allersp,alltimes,allfreqs,base_tlims,base_flims,...
    varargin)
%## INPUTS
%-
% IN:
%       allersp, CELL
%           should be log-transformed time-frequency data (freq x time x subjs)
%       alltimes, DOUBLE/SINGLE          
%           vector of times in milliseconds with size(allersp{1},2)
%       allfreqs, DOUBLE/SINGLE
%           vector of frequencies in hertz with size(allersp{1},1)
%       basetimes, DOUBLE/SINGLE
%           a pair of times(ms) [BEGINNING,END] to baseline allersp to
%       basefreqs, DOUBLE/SINGLE
%           a pair of freqs(hz) [BEGINNING,END] to baseline allersp to.
%           This is more for plotting and doesn't have a major effect on
%           baselining outcome.
%-
DO_COMMON_BASE = false;
DO_SUBJ_BASE = true;
GROUPWISE_PARAMS = struct('single_grp_treatment',false,...
    'group_wts',[1,1,1]);
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'allersp',@iscell);
addRequired(p,'alltimes',@isnumeric);
addRequired(p,'allfreqs',@isnumeric);
addRequired(p,'base_tlims',@isnumeric);
addRequired(p,'base_flims',@isnumeric);
%## OPTIONAL
%## PARAMETER
% addParameter(p,'PLOT_STRUCT',PLOT_STRUCT,@(x) validate_struct(x,PLOT_STRUCT));
addParameter(p,'DO_COMMON_BASE',DO_COMMON_BASE,@islogical);
addParameter(p,'DO_SUBJ_BASE',DO_SUBJ_BASE,@islogical);
addParameter(p,'GROUPWISE_PARAMS',GROUPWISE_PARAMS,@struct);
parse(p,allersp,alltimes,allfreqs,base_tlims,base_flims,varargin{:});
%## SET DEFAULTS
DO_COMMON_BASE = p.Results.DO_COMMON_BASE;
DO_SUBJ_BASE = p.Results.DO_SUBJ_BASE;
GROUPWISE_PARAMS = p.Results.GROUPWISE_PARAMS;
%%
baseln_com = zeros(size(allersp{1},1),size(allersp{1},3));
baseln_trial = cell(size(allersp));
base_tinds = find(alltimes>=base_tlims(1) & alltimes<=base_tlims(end));
base_finds = find(allfreqs>=base_flims(1) & allfreqs<=base_flims(end));
tmp_nolog = zeros(size(allersp{1},1),size(allersp{1},2),size(allersp,1));
allersp_bcom1 = cell(size(allersp));
allersp_bcom2 = cell(size(allersp));
% allersp_log_subBase = cell(size(allersp));
allersp_btrial = cell(size(allersp));
tmp_allersp = allersp;
%## WITHIN CONDITION & ?GROUP BASELINE
if DO_SUBJ_BASE
    for i = 1:size(tmp_allersp,2)
        for j = 1:size(tmp_allersp,1)
            erspdata = tmp_allersp{j,i}(:,:,:);
            %- mean power for each person
            tmp_base_subj = mean(erspdata(:,base_tinds,:),2);
            %- subtract out each subjects baseline for each condition
            allersp_btrial{j,i} = tmp_allersp{j,i}(:,:,:)-repmat(tmp_base_subj,[1,length(alltimes),1]);
            %- store base
            baseln_trial{j,i} = mean(tmp_base_subj,3);
        end
    end
    tmp_allersp = allersp_btrial;
end
%## BETWEEN CONDITION BASELINE
if DO_COMMON_BASE
%     tmp = tmp_allersp;
    %- loop over groups
    for i = 1:size(tmp_allersp,2)
        %- loop over subjects
        for n = 1:size(tmp_allersp{1,i},3)
            %- convert log-spectrum to non-log
            for j = 1:size(tmp_allersp,1)
                tmp_nolog(:,:,j) = 10.^(tmp_allersp{j,i}(:,:,n)/20);
            end
            %- mean across conditions then mean across times 
            tmp_base_com = mean(mean(tmp_nolog(:,base_tinds,:),3),2);

            %- log-subtraction of baseline for a particular subject
            for j = 1:size(tmp_allersp,1)
                allersp_bcom1{j,i}(:,:,n) = tmp_nolog(:,:,j)./(repmat(tmp_base_com,1,size(tmp_nolog,2)));
                allersp_bcom2{j,i}(:,:,n) = 20*log10(allersp_bcom1{j,i}(:,:,n)); 
            end
            %- store
            baseln_com(:,n) = tmp_base_com;
        end
        %- store
%         tmp{:,i} = tmp_allersp
    end
    %- store
    tmp_allersp = allersp_bcom2;
end
allersp_full = tmp_allersp;
allersp_crop = cellfun(@(x) x(base_finds,base_tinds,:),tmp_allersp,'uniformoutput',false);
end