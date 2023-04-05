function [ALLEEG] = wrapper_cnctanl_stats(ALLEEG,varargin)
%WRAPPER_CNCTANL_STATS Summary of this function goes here
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
% Code Designer: Jacob Salminen
% Code Date: 12/30/2022, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%## TIME
tic
%## DEFINE DEFAULTS
DO_PHASERAND = true;
DO_BOOTSTRAP = false;
p = inputParser;
%## REQUIRED
addRequired(p,'ALLEEG',@isstruct);
%## OPTIONAL
%## PARAMETER
addParameter(p,'DO_PHASERAND',DO_PHASERAND,@islogical);
addParameter(p,'DO_BOOTSTRAP',DO_BOOTSTRAP,@islogical);
parse(p,ALLEEG,varargin{:});
%## SET DEFAULTS
%- OPTIONALS
%- PARAMETER
DO_BOOTSTRAP = p.Results.DO_BOOTSTRAP;
DO_PHASERAND = p.Results.DO_PHASERAND;
%- Define Defaults
%% ===================================================================== %%
%## PHASE RANDOMIZATION && NONZERO TEST
if DO_PHASERAND
    for cond_i=1:length(ALLEEG)
        %- Generate Phase Randomized Distribution
        fprintf('\n==== PHASE RANDOMIZING CONNECTIVITY MEASURES ====\n')
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        fPath = ALLEEG(cond_i).filepath;
        if ~exist([fPath filesep fName],'file')
            [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'PhaseRnd');
            %- Save Phase randomized distribution
            tmp = ALLEEG(cond_i).CAT.PConn;
            par_save(tmp,fPath,fName,'_PhaseRnd');
            fprintf('done.\n')
            %- Conduct Nonzero Test (Need Phase Randomized Distribution)
            fprintf('\n==== NONZERO STATISTICS ====\n')
            [ALLEEG(cond_i),~] = cnctanl_groupStats(ALLEEG(cond_i),'NonZero');
            %- Save Nonzero Stats
            tmp = ALLEEG(cond_i).CAT.Stats;
            fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
            par_save(tmp,ALLEEG(cond_i).filepath,fName,'_NonZero');
            fprintf('done.\n')
        else
            ALLEEG(cond_i) = cnctanl_loadCAT(ALLEEG(cond_i),'NonzeroTest');
        end
        
        
    end
end
%## BOOTSTRAP
if DO_BOOTSTRAP
    for cond_i=1:length(ALLEEG)
        fprintf('\n==== CALCULATING BOOTSTRAP MEASURES ====\n')
        %- calculate BootStrap distribution
        ALLEEG(cond_i) = cnctanl_groupStats(ALLEEG(cond_i),'BootStrap');    
        %- save BootStrap distribution=
        tmp = ALLEEG(cond_i).CAT.PConn;
        fName = strsplit(ALLEEG(cond_i).filename,'.'); fName = [fName{1} '.mat'];
        par_save(tmp,ALLEEG(cond_i).filepath,fName,'_BootStrap');
    end
end
%## TIME
toc
end

