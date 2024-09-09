function [EEG,cmd_out] = mim_prep_hpg_amica(EEG,float_fPath,amica_out_fPath,email_char,avg_ref_pca_reduction,varargin)
%       EEG is the EEG struct you want to save on M drive and run AMICA on via the
%       hipergator
%
%   IN: 
%       EEG, STRUCT (EEGLAB)
%           is the EEG struct you want to save on M drive and run AMICA on via the
%           hipergator           
%       float_fPath, CHAR
%           is the path where your EEG float (.fdt) file is located. (hint:
%           its usually where your .set file is.
%       amica_out_fPath, CHAR
%           place where you would like your AMICA .params & .sh files to be
%           saved. (recommended: save them to the same folder as your
%           cleaned .set and .fdt file.
%       email_char, CHAR
%           is a string with your email if you want to be notified when the
%           hipergator ran your data
%       avg_ref_pca_reduction, DOUBLE
%           is a number. It can be 1 or 2 or maybe even more.
%           Example: set to 3 if using EEG, EMG, and Noise recordings all separately
%           avg ref to each other.
%
%   OUT: 
%       EEG, struct
%           
%       cmd_out, CHAR
%           a command line entry that runs the .sh file use srun & pmix_v3  
%   IMPORTANT: 
%fileNameNoExt,amicaOutputFolder_local,...
%     amicaOutputFolder_unix,avgRefPCAReduction,emailStr
%## TIME
tic
%## DEFINE DEFAULTS
%- HARD DEFINES
NUM_MODELS          = 1;
MAX_ITERS           = 2500;
MAX_DURATION        = hours(3); %minutes, days, etc.
NUM_NODES           = 64; % How many nodes to request
NUM_TASKS           = NUM_NODES; % Number of MPI jobs
NUM_MEM             = ceil(512*NUM_MODELS*1.5); % memory per cpu used (in MB) %Changed to 2*512=1024 (06/16/2022)
NUM_TASKS_PER_NODE  = 1; % ryan edit: noticed the code ran best when you only used a single CPU from each node so i would generally set this to 1
MAX_DURATION.Format = 'hh:mm:ss';
NUM_TIME            = char(MAX_DURATION); % Wall time hh:mm:ss
QOS_NAME            = 'dferris-b'; % Use 'dferris' or 'dferris-b' (burst allocation)
%- DEFAULTS
AMICA_DLL_PATH = [filesep 'blue' filesep 'dferris' filesep,...
                'share' filesep 's.peterson' filesep 'test' filesep,...
                'AMICA_15' filesep 'amica15ub'];
disp('Converting data to double...');
tmpdata = double(EEG.data);
PCA_KEEP = min([rank(tmpdata), EEG.nbchan-avg_ref_pca_reduction]); % PCA value to use

p = inputParser;
%## REQUIRED
addRequired(p,'EEG',@isstruct)
addRequired(p,'float_fPath',@ischar);
addRequired(p,'amica_out_fPath',@ischar);
addRequired(p,'email_char',@ischar);
addRequired(p,'avg_ref_pca_reduction',@isnumeric);
%## OPTIONAL
%## PARAMETER
addParameter(p,'PCA_KEEP',PCA_KEEP,@isnumeric);
addParameter(p,'AMICA_DLL_PATH',AMICA_DLL_PATH,@ischar);
%## SET DEFAULTS
parse(p,EEG,float_fPath,amica_out_fPath,email_char,avg_ref_pca_reduction,varargin{:});
AMICA_DLL_PATH = p.Results.AMICA_DLL_PATH;
PCA_KEEP = p.Results.PCA_KEEP;
%% ===================================================================== %%
% Number of PCA components should be the number of channels minus the
% reduction value: see. getRank()
if ~exist(float_fPath,'file')
    error('ERROR. %s does not exist',float_fPath);
end
if ~ispc
    amica_out_fPath = convertPath2UNIX(amica_out_fPath);
else
    amica_out_fPath = convertPath2Drive(amica_out_fPath);
end
if ~exist(amica_out_fPath,'file')
    fprintf('Making directory ''amica_out_fPath'' %s...\n',amica_out_fPath);
    mkdir(amica_out_fPath);
end
%## STORE AMICA/BASH PARAMS
EEG.etc.amica_run.pcaKeep   = PCA_KEEP;
EEG.etc.amica_run.numChans  = EEG.nbchan;
EEG.etc.amica_run.numFrames = length(EEG.times);
EEG.etc.amica_run.amica_out_fPath = amica_out_fPath;

%## CREATE PARAM FILE
[~,param_out_fPath] = make_param_amica( float_fPath, amica_out_fPath, PCA_KEEP,...
    EEG.nbchan, length(EEG.times), NUM_MODELS, MAX_ITERS);

%## CREATE HIPERGATOR BASH FILE
% bash_out_fPath = [amica_out_fPath filesep 'run_amica_hipergator.sh'];
[~,bash_out_fPath] = make_amica_bash(EEG, amica_out_fPath, param_out_fPath, AMICA_DLL_PATH,...
                    email_char, NUM_NODES, NUM_TASKS, NUM_TASKS_PER_NODE, NUM_MEM, NUM_TIME, QOS_NAME);
fprintf('Done creating .param and .sh files for AMICA.\n');
% disp('Done creating AMICA stuff');
% disp('Use MobaXTerm to submit your job (use right click to paste)');
% disp(newline);
% disp(['cd ' convertPath2UNIX(scriptFileDir,'dferris')]);
% disp(['sbatch runAMICA_hipergator.sh']);
% disp('squeue -A dferris');

cmd_out = {['cd ' convertPath2UNIX(amica_out_fPath)];...
              sprintf('sbatch  %s',convertPath2UNIX(bash_out_fPath))};
%% ===================================================================== %%
function [fid,out_fPath] = make_param_amica(float_fPath, out_fPath, pcaKeep, numChans, numFrames, numModels, maxIter)
    % MAKE FILE AND ADD PARAMS
    PARAM_FILE_FNAME = 'input_hipergator.param';
    fid = fopen([out_fPath filesep PARAM_FILE_FNAME],'w+'); %changed from 'w'
    fprintf(fid,'files %s\n',convertPath2UNIX(float_fPath));
    fprintf(fid,'outdir %s\n',convertPath2UNIX(out_fPath));
    fprintf(fid,'num_models %d\n',numModels); %ryan temp switched to 3 
    fprintf(fid,'num_mix_comps 3\n');
    fprintf(fid,'pdftype 0\n');
    fprintf(fid,'block_size 128\n');
    fprintf(fid,'max_iter %d\n',maxIter); 
    fprintf(fid,'num_samples 1\n');
    fprintf(fid,'data_dim %d\n',numChans);
    fprintf(fid,'field_dim %d\n',numFrames);
    fprintf(fid,'field_blocksize 1\n');
    fprintf(fid,'share_comps 0\n'); %ryan temp switched to 1 %SHARE COMPS HERE
    fprintf(fid,'share_start 100\n');
    fprintf(fid,'comp_thresh 0.990000\n');
    fprintf(fid,'share_iter 100\n');
    fprintf(fid,'lrate 0.100000\n');
    fprintf(fid,'minlrate 1.000000e-08\n');
    fprintf(fid,'lratefact 0.500000\n');
    fprintf(fid,'rholrate 0.050000\n');
    fprintf(fid,'rho0 1.500000\n');
    fprintf(fid,'minrho 1.000000\n');
    fprintf(fid,'maxrho 2.000000\n');
    fprintf(fid,'rholratefact 0.500000\n');
    fprintf(fid,'kurt_start 3\n');
    fprintf(fid,'num_kurt 5\n');
    fprintf(fid,'kurt_int 1\n');
    fprintf(fid,'do_newton 1\n');
    fprintf(fid,'newt_start 50\n');
    fprintf(fid,'newt_ramp 10\n');
    fprintf(fid,'newtrate 1.000000\n');
    fprintf(fid,'do_reject 1\n'); %2021-12-10 RJD turned on do_reject based on makoto processing pipeline recommendation
    %https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code#What_do_those_AMICA_parameters_mean.3F_.2804.2F06.2F2018_updated.29
    %"If you want to know which data points were rejected, you can check
    %EEG.etc.amica.Lht. If any datapoint shows 0, it means the datapoints were rejected by AMICA. 
    %Note also that you might want to use this info for further data cleaning if you don't mind rejecting randomly 
    %distributing single datapoints."
    fprintf(fid,'numrej 15\n');
    fprintf(fid,'rejsig 3.000000\n');
    fprintf(fid,'rejstart 1\n');
    fprintf(fid,'rejint 1\n');
    %note: there can be memory issues if your dataset is large and num_tasks is
    %greater than 1. You can either 1) use PCA reduction on the data or 2) set 
    %max_threads to 1 and set num_tasks to be a large number with one task per 
    %node. AMICA runs quickly this way
    % fprintf(fid,['max_threads %d\n'],num_tasks); 
    fprintf(fid,'max_threads %d\n',1); 
    fprintf(fid,'writestep 250\n'); %ryan switched from 10 to 250 since performance was slowing down a lot during writing
    fprintf(fid,'write_nd 0\n');
    fprintf(fid,'write_LLt 1\n'); %ryan switched from 0 to 1
    fprintf(fid,'decwindow 1\n');
    fprintf(fid,'max_decs 3\n');
    fprintf(fid,'update_A 1\n');
    fprintf(fid,'update_c 1\n');
    fprintf(fid,'update_gm 1\n');
    fprintf(fid,'update_alpha 1\n');
    fprintf(fid,'update_mu 1\n');
    fprintf(fid,'update_beta 1\n');
    fprintf(fid,'invsigmax 100.000000\n');
    fprintf(fid,'invsigmin 0.000000\n');
    fprintf(fid,'do_rho 1\n');
    fprintf(fid,'load_rej 0\n');
    fprintf(fid,'load_W 0\n');
    fprintf(fid,'load_c 0\n');
    fprintf(fid,'load_gm 0\n');
    fprintf(fid,'load_alpha 0\n');
    fprintf(fid,'load_mu 0\n');
    fprintf(fid,'load_beta 0\n');
    fprintf(fid,'load_rho 0\n');
    fprintf(fid,'load_comp_list 0\n');
    fprintf(fid,'do_mean 1\n');
    fprintf(fid,'do_sphere 1\n');
    fprintf(fid,'doPCA 1\n');
    fprintf(fid,'pcakeep %d\n',pcaKeep);
    fprintf(fid,'pcadb 30.000000\n');
    fprintf(fid,'byte_size 4\n');
    fprintf(fid,'doscaling 1\n');
    fprintf(fid,'scalestep 1\n');
    fclose(fid);
    out_fPath = [out_fPath filesep PARAM_FILE_FNAME];
end


function [fid,sh_out_fPath] = make_amica_bash(EEG, out_fPath, param_fPath, amica_dll_fPath,...
        emailStr, numNodes, numTasks, numTasksPerNode, numMem, numTime, qosName)
    AMICA_BASH_FNAME = 'run_amica_hipergator.sh';
    sh_out_fPath = [out_fPath filesep AMICA_BASH_FNAME];
    %## CREATE BASH FILE
    fid = fopen(sh_out_fPath,'w+'); %changed from 'w'
    %- SHEBANG
    fprintf(fid,'#!/bin/bash\n'); % change from 'sh' to 'bash', supposedly better features?
    %- JOB NAME
    fprintf(fid,'#SBATCH --job-name=%s_AMICA\n',EEG.subject); % Job name
    %- MAIL AND STANDARD OUTPUT: Mail events (NONE, BEGIN, END, FAIL, ALL)
    fprintf(fid,'#SBATCH --mail-type=ALL\n');
    fprintf(fid,'#SBATCH --mail-user=%s\n',emailStr); % Where to send mail
    fprintf(fid,'#SBATCH --output=%s/%%j_amica_out.log\n', convertPath2UNIX(out_fPath)); % Standard output and error log
    %- NODES, NTASKS, & MEMORY MGMT
    fprintf(fid,'#SBATCH --nodes=%d\n',numNodes); % Use one node
    fprintf(fid,'#SBATCH --ntasks=%d\n',numTasks); % Run a single task
    fprintf(fid,'#SBATCH --ntasks-per-node=%d\n',numTasksPerNode); % number of tasks per node; ryan edit: added this to have more control over how tasks are split across nodes
    fprintf(fid,'#SBATCH --mem-per-cpu=%dmb\n',numMem); % Total memory limit
    fprintf(fid,'#SBATCH --distribution=cyclic:cyclic\n'); % Distribute tasks cyclically first among nodes and then among sockets within a node
    %- ALLOCATABLE TIME
    fprintf(fid,'#SBATCH --time=%s\n',numTime); % Time limit hrs:min:sec
    fprintf(fid,'#SBATCH --account=dferris\n'); % Account name
    fprintf(fid,'#SBATCH --qos=%s\n',qosName); % Quality of service name
    fprintf(fid,'#SBATCH --partition=hpg2-compute\n'); % added 12/15/2020 because recommended by hipergator IT when random nodes on older part of cluster would cause AMICA to crash. This limits selection to hipergator 2
    %- ECHO DETAILS & MODULE LOAD
    fprintf(fid,'cd %s\n',convertPath2UNIX(out_fPath));
    fprintf(fid,['echo "Date              = $(date)"\n',...
                 'echo "Hostname          = $(hostname -s)"\n',...
                 'echo "Working Directory = $(pwd)"\n',...
                 'echo ""\n',...
                 'echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"\n',...
                 'echo "Number of Tasks Allocated      = $SLURM_NTASKS"\n',...
                 'echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"\n',...
                 'echo "slurm_mem_per_cpu $SLURM_MEM_PER_CPU"\n',...
                 'echo "slurm_mem_per_gpu $SLURM_MEM_PER_GPU"\n',...
                 'echo "slurm_mem_per_node $SLURM_MEM_PER_NODE"\n']);
%     fprintf(fid,'module purge\n'); % (03/07/2023) JS, not sue if this is needed, but its use is encouraged on the HiperGator Wiki.
    
    fprintf(fid,'module load ufrc\n');
    fprintf(fid,'module load intel/2020 openmpi/4.1.5\n');
    fprintf(fid,'srun --mpi=pmix_v3 %s %s\n',convertPath2UNIX(amica_dll_fPath), convertPath2UNIX(param_fPath));
    fclose(fid);
end

end