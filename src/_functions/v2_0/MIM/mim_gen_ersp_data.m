function [] = mim_gen_ersp_data(STUDY,ALLEEG,tmpSTUDY,tmpSTUDY_commonbase, ...
    myErspParams,clusternum,outputdir,file_keyword,varargin)
% CAT CODE
%  _._     _,-'""`-._
% (,-.`._,'(       |\`-/|
%     `-.-' \ )-`( , o o)
%           `-    \`_`"'-
% Code Designers: Chang Liu, Jacob Salminen
% Code Date: 07/17/2023, MATLAB 2020b
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Copyright (C) Chang Liu,
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
%## Define Parser
p = inputParser;
%## REQUIRED
addRequired(p,'STUDY',@isstruct);
addRequired(p,'ALLEEG',@isstruct);
addRequired(p,'study_fName',@ischar);
addRequired(p,'study_fPath',@ischar);
%## OPTIONAL
%## PARAMETER
parse(p,ALLEEG,study_fName,study_fPath,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%##
% DIP_NUM = 1;
% DIP_PLOT = 'off';
%## MAKE DIRS
if ~exist(study_fPath,'dir')
    mkdir(study_fPath);
end
%% ===================================================================== %%
%
%{
        disp('Gathering esrp WITHOUT single trial normalization or any normalization')
        tic
        [STUDY,allerspdata, alltimes, allfreqs, pgroup,pcon,pinter] = std_erspplot(STUDY, ALLEEG, 'clusters', clusternum,...
            'freqrange', [0 200]); %PLOTS THE ENTIRE cluster_i'th CLUSTER ,% did I apply subtractsubjectmean???, no correction at all with the current setting
        save(fullfile(outputdir,['readESRP_', num2str(STUDY.currentdesign),file_keyword,'.mat']), 'allerspdata', 'alltimes', 'allfreqs', 'pgroup', 'pcon', 'pinter');
        saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword, '_All_Comps_ERSP.fig']))
        close
        toc
%}
        % ----------------------------------------
        % In the std_erspplot_customParams, it also had the
        % std_readata_customParams. 
        disp('Gathering esrp after baseline correction using warp(1) to warp(5)')
        tic
        disp(tmpSTUDY.etc.statistics)
        disp(myErspParams);
        [tmpSTUDY,allerspdata2, alltimes2, allfreqs2, pgroup2,pcon2,pinter2] = std_erspplot_customParams(tmpSTUDY, ALLEEG, myErspParams, 'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_', num2str(STUDY.currentdesign),file_keyword,'_subBase.mat']), 'allerspdata2', 'alltimes2', 'allfreqs2', 'pgroup2', 'pcon2', 'pinter2');
        saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword, '_All_Comps_ERSP_subBase.fig']))
% 
%             toc

        % ---------------------------------------- Commonbase!!!
        disp('Gathering esrp after baseline correction using warp(1) to warp(5) and apply commonbase')
        tic
        disp(tmpSTUDY_commonbase.etc.statistics)
        disp(tmpSTUDY_commonbase.etc.erspparams)
        [tmpSTUDY_commonbase_out,allerspdata3, alltimes3, allfreqs3, pgroup3,pcon3,pinter3] = std_erspplot_customParams(tmpSTUDY_commonbase, ALLEEG, myErspParams, 'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_', num2str(STUDY.currentdesign),file_keyword,'_subBase_commonBase.mat']), 'allerspdata3', 'alltimes3', 'allfreqs3', 'pgroup3', 'pcon3', 'pinter3');
        saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword, '_All_Comps_ERSP_subBase_commonBase.fig']))

                        % ---------------------------------------- Commonbase!!!
%}
%{
        disp('Gathering esrp after fulltrial correction and baseline correction using warp(1) to warp(5) and apply commonbase')
        tic
        disp(tmpSTUDY_commonbase.etc.statistics)
        disp(tmpSTUDY_commonbase.etc.erspparams)
        % if trialbase = 'on', no commonbaseline subtraction is
        % performed in newtimefbaseln
        [tmpSTUDY_commonbase_out,allerspdata4, alltimes4, allfreqs4, pgroup4,pcon4,pinter4] = std_erspplot_customParams(tmpSTUDY_commonbase, ALLEEG, myErspParams_trialbasefull,...
            'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_', num2str(STUDY.currentdesign),file_keyword,'_subBase_fullTrialbase_commonBase.mat']), 'allerspdata4', 'alltimes4', 'allfreqs4', 'pgroup4', 'pcon4', 'pinter4');
        saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword, '_All_Comps_ERSP_subBase_fullTrialbase_commonBase.fig']))
        %}
        % 
%}
end

