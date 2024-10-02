classdef workspace_init
    %WORKSPACE_INIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        src_dir (1,:) char {exist} = ''
        submods_dir (1,:) char {exist} = ''
    end
    
    methods
        function obj = workspace_init()
            %WORKSPACE Construct an instance of this class
            %   Detailed explanation goes here
            tmp_pwd = dir(['.' filesep]);
            tmp_pwd = tmp_pwd(1).folder;
            [src_dir,submods_dir] = obj.get_dirs(tmp_pwd);
            obj.src_dir = src_dir;
            obj.submods_dir = submods_dir;
        end
        
        function [src_dir,submod_dir] = get_dirs(~,tmppwd,src_char,submod_char)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            tmp = strsplit(tmppwd,filesep);
            src_ind = find(strcmp(tmp,src_char));
            %- Add source directory where setup functions are
            if ~ispc
                src_dir = [filesep strjoin(tmp(1:src_ind),filesep)];
            else
                src_dir = strjoin(tmp(1:src_ind),filesep);
            end
            %- Add submodules directory where packages are 
            if ~ispc
                submod_dir = [filesep strjoin(tmp(1:src_ind-1),filesep) filesep submod_char];
            else
                submod_dir = [strjoin(tmp(1:src_ind-1),filesep) filesep submod_char];
            end
        end

        function set_paths(~)
            
        end
    end
end

