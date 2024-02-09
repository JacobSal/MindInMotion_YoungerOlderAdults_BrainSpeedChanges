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
parse(p,allersp,alltimes,allfreqs,base_tlims,base_flims,varargin{:});
%## SET DEFAULTS
DO_COMMON_BASE = p.Results.DO_COMMON_BASE;
DO_SUBJ_BASE = p.Results.DO_SUBJ_BASE;
%%
baseln_com = zeros(size(allersp{1},1),size(allersp{1},3));
baseln_trial = cell(size(allersp));
base_tinds = find(alltimes>=base_tlims(1) & alltimes<=base_tlims(end));
base_finds = find(allfreqs>=base_flims(1) & allfreqs<=base_flims(end));
tmp_nolog = zeros(size(allersp{1},1),size(allersp{1},2));
allersp_bcom1 = cell(size(allersp));
allersp_bcom2 = cell(size(allersp));
% allersp_log_subBase = cell(size(allersp));
allersp_btrial = cell(size(allersp));
tmp_allersp = allersp;
%## WITHIN CONDITION & ?GROUP BASELINE
if DO_SUBJ_BASE
    for j = 1: length(tmp_allersp)
        erspdata = tmp_allersp{j}(:,:,:);
        %- mean power for each person
        tmp_base_subj = mean(erspdata(:,base_tinds,:),2);
        %- subtract out each subjects baseline for each condition
        allersp_btrial{j,1} = tmp_allersp{j}(:,:,:)-repmat(tmp_base_subj,[1,length(alltimes),1]);
        %- store base
        baseln_trial{j,1} = mean(tmp_base_subj,3);
    end
    tmp_allersp = allersp_btrial;
end
%## BETWEEN CONDITION BASELINE
if DO_COMMON_BASE
    %- loop over subjects
    for n = 1:size(tmp_allersp{1},3)
        %- convert log-spectrum to non-log
        for j = 1:4
            tmp_nolog(:,:,j) = 10.^(tmp_allersp{j}(:,:,n)/20);
        end
        %- mean across conditions then mean across times 
        tmp_base_com = mean(mean(tmp_nolog(:,base_tinds,:),3),2);
        
        %- log-subtraction of baseline for a particular subject
        for j = 1:4
            allersp_bcom1{j,1}(:,:,n) = tmp_nolog(:,:,j)./(repmat(tmp_base_com,1,size(tmp_nolog,2)));
            allersp_bcom2{j,1}(:,:,n) = 20*log(allersp_bcom1{j,1}(:,:,n)); 
        end
        %- store
        baseln_com(:,n) = tmp_base_com;
    end
    %- store
    tmp_allersp = allersp_bcom2;
end
allersp_full = tmp_allersp;
allersp_crop = cellfun(@(x) x(base_finds,base_tinds,:),tmp_allersp,'uniformoutput',false);
end
%% ===================================================================== %%
function [b] = validate_struct(x,DEFAULT_STRUCT)
    b = false;
    struct_name = inputname(2);
    %##
    fs1 = fields(x);
    fs2 = fields(DEFAULT_STRUCT);
    vals1 = struct2cell(x);
    vals2 = struct2cell(DEFAULT_STRUCT);
    %- check field names
    chk = cellfun(@(x) any(strcmp(x,fs2)),fs1);
    if ~all(chk)
        fprintf(2,'\nFields for struct do not match for %s\n',struct_name);
        return
    end
    %- check field value's class type
    for f = 1:length(fs2)
        ind = strcmp(fs2{f},fs1);
        chk = strcmp(class(vals2{f}),class(vals1{ind}));
        if ~chk
            fprintf(2,'\nStruct.%s must be type %s, but is type %s\n',fs2{f},class(vals2{f}),class(vals1{ind}));
            return
        end
    end
    b = true;
end