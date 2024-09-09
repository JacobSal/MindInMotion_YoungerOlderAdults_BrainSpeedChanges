% ==============================================================
%
%    GENERAL
%
%
%    INPUT/S
%   
%        
%    OUTPUT/S
%
%      -
%
%    PENDING WORK
%
%      -
%
%    KNOWN BUG/S
%
%      -
%
%    COMMENT/S
%
%      -
%
%    RELATED FUNCTION/S
%
%      
%
%    ABOUT
%
%      -Created:     November 2003
%      -Last update: 
%      -Revision:    0.3.1
%      -Author:      R. S. Schestowitz, University of Manchester
% Copyright (c) 2016, Roy Schestowitz
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% ==============================================================

% Instructions: follow the three simple steps below -- (1), (2) and (3)

max=100;
          % (1) Set this to the total number of iterations

progress_bar_position = 0;

time_for_this_iteration = 0.01;
          % (2) Provide initial time estimate for one iteration

for i=1:max,
	   tic;
	   
	   
	   % (3) Place all computations here
	   
	   
	   progress_bar_position = progress_bar_position + 1 / max;
           clc;
           disp(['|=================================================|']);
           progress_string='|';       
           for counter = 1:floor(progress_bar_position * 100 / 2),
               progress_string = [progress_string, '#'];
           end
           disp(progress_string);
           disp(['|================= ',num2str(floor(progress_bar_position * 100)),'% completed =================|']);
                          % display progress per cent
           steps_remaining = max - i;
           minutes = floor(time_for_this_iteration * steps_remaining / 60);
           seconds = rem(floor(time_for_this_iteration *  steps_remaining), 60);
           disp(' ');
           if (seconds > 9),
             disp(['            Estimated remaining time: ', num2str(minutes), ':', num2str(seconds)]);
                          % show time indicators
           else
             disp(['            Estimated remaining time: ', num2str(minutes), ':0', num2str(seconds)]);
           end
           time_for_this_iteration = toc;
end
