function [SUBJSTRUCT] = subj_createParStruct(PATHS, studyName, studyPath, subjNames, subjFolders, subFolderNames, subFolderCodes, varargin)
%GETSUBJDIRDATA Summary of this function goes here
%   Detailed explanation goes here
%   
%   IN: 
%       PATHS, STRUCT
%          
%       studyName, CHAR
%          
%       studyPath, CHAR
%          
%       subjNames, CELL of CHARS
%           
%       subjFolders, CELL of CHARS
%           
%       subFolderNames, CELL of CHARS
%           
%       subFolderCodes, CELL of CHARS
%           '1' ::: RAW DIRECTORY; 
%           '2' ::: ICA DIRECTORY; 
%           '3' ::: EPOCHED ICA/PROCESSED DIRECTORY; 
%           '4' ::: HEADMODEL DIRECTORY; 
%           '5' ::: HYPERCOMPUTER OUTPUTS; 
%           '6' ::: 
%
%   OUT: 
%   IMPORTANT: 

%## TIME
tic
%## DEFINE DEFAULTS

Defaults = {};
p = inputParser;
%## REQUIRED
addRequired(p,'PATHS',@isstruct)
addRequired(p,'studyName',@ischar);
addRequired(p,'studyPath',@ischar);
addRequired(p,'subjNames',@iscell)
addRequired(p,'subjFolders',@iscell)
addRequired(p,'subFolderNames',@iscell)
addRequired(p,'subFolderCodes',@iscell)
%## OPTIONAL
%## PARAMETER
parse(p, PATHS, studyName, studyPath, subjNames, subjFolders, subFolderNames, subFolderCodes, varargin{:});
%## SET DEFAULTS
PATH_PREFIX = [];
LOAD_DATA = false;
SAVE_DATA = true;
%## SAVE FOLDER HANDLER (HIGHEST ORDER FUNCTION ONLY)
% ----------------------------------------------------------------------- %
ln1 = length(subjNames);
ln2 = length(subjFolders);
ln3 = length(subFolderNames);
ln4 = length(subFolderCodes);

disp('Assigning patient data to subject strucutre');
path_i = struct('label',[],'filepath',[],'filename',[]);
path_i_s = [path_i; path_i];
pathsSubj = struct('conditions',path_i_s,...
                   'groups',path_i_s,...
                   'sessions',path_i_s);
pathsEEG = struct('raw',path_i_s,...
                  'epoched',path_i_s,...
                  'headmodel',path_i_s, ...
                  'cnctAnl',path_i_s,...
                  'ica',path_i_s,...
                  'sourcelocalize',path_i_s);
sourceLocalize = struct('mrifile',path_i_s,...
                        'hasAcpc',0,...
                        'hasFid',0,...
                        'acpcNames',{{'ac','pc','xzpoint','rpoint'}},...
                        'fidNames',{{'nasion','lhj','rhj','xzpoint'}},...
                        'acpcCoords',{{[0,0,0],[0,0,0],[0,0,0],[0,0,0]}},...
                        'fidCoords',{{[0,0,0],[0,0,0],[0,0,0],[0,0,0]}});
META = struct('epoched',[],...
              'cnctanl',[]);
cnctAnl = struct('current',path_i_s,...
                 'epoched',path_i_s);
SUBJSTRUCT = struct('subjStr','',...
                       'folder',[],...
                       'subfolders',[],...
                       'subfolderRoles',[],...
                       'pathsSubj',pathsSubj,...
                       'pathsEEG',pathsEEG,...
                       'sourceLocalize',sourceLocalize,...
                       'cnctAnl',cnctAnl,...
                       'META',META,...
                       'SAVE',path_i,...
                       'PATHS',PATHS);


if (ln1 ~= ln2)
    error('SUBJSTRUCT:subjects','the length of ''subjNames'' must equal ''subjFolders''');
end
if (ln3 ~= ln4)
    error('SUBJSTRUCT:subjects','the length of ''subFolderNames'' must equal ''subFolderCodes''');
end 
tmpS = SUBJSTRUCT;
%## LOOP & MAKE ENTRIES
for cnd = 1:ln1
    % append another entry
    if cnd > 1
        SUBJSTRUCT = horzcat(SUBJSTRUCT,tmpS);
    end
    % assign inputs 
    SUBJSTRUCT(cnd).SAVE.filepath = studyPath;
    SUBJSTRUCT(cnd).SAVE.filename = sprintf('%s.mat',studyName);
    SUBJSTRUCT(cnd).SAVE.label = studyName;
    SUBJSTRUCT(cnd).folder = subjFolders{cnd};
    SUBJSTRUCT(cnd).subjStr = subjNames{cnd};
    SUBJSTRUCT(cnd).subfolders = subFolderNames;
    SUBJSTRUCT(cnd).subfolderRoles = subFolderCodes;
end

%## WRAP UP & SAVE
if SAVE_DATA
    SUBJSTRUCT = subj_saveStruct(SUBJSTRUCT,PATH_PREFIX);
elseif LOAD_DATA && isfile([SUBJSTRUCT(1).SAVE.filepath filesep SUBJSTRUCT(1).SAVE.filename])
    fprintf('Study ''%s'' already exists. Loading...\n',studyName);
    SUBJSTRUCT = par_load(SUBJSTRUCT(1).SAVE.filepath,SUBJSTRUCT(1).SAVE.filename,PATH_PREFIX);
else
    warning('Study ''%s'' does not exist, and saving was not specified.\nNo file saved.\n',studyName)
end

%## TIME
toc
end

