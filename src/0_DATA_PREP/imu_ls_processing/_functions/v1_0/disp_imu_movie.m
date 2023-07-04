function [] = disp_imu_movie()
%DISP_IMU_MOVIE Summary of this function goes here
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
% Code Date: 02/13/2023, MATLAB 2019a
% Copyright (C) Jacob Salminen, jsalminen@ufl.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

%% VALIDATION PLOTS
% movie_samplePeriod  = .1;%.02 default %1/EEG.srate;

%plot in global frame (don't have to deal with trying to figure out how to
%rotate the quatPlot to show everything in the body frame of reference)
posPlot = world_frame_pos; %5-14-21 we could crop here but might get discontinuity in the video and we don't care the video is just for fun
quatPlot = quaternConj(sensor_frame_ori); %not sure why need to use conjugate but it looks correct compared to other way

% %2020-01-05 Note: temp crop (remove later) %RD commented out on
% 2021-03-24
% posPlot = posPlot(5000:end-5000,:);
% quatPlot = quatPlot(5000:end-5000,:);

% Extend final sample to delay end of animation
extraTime = 1;
onesVector = ones(extraTime*(1/movie_samplePeriod), 1);
posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];

% Create 6 DOF animation
SamplePlotFreq = 20; %20 for 500 Hz to appear normal speed?

Spin = 120*2;
spinMat = [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)]; %original
% spinMat = [0*ones(length(posPlot),1),0*ones(length(posPlot), 1)];
%0,0 x= horiz, z = vert (side or frontal)
%90,0 y = horiz, z = vert (side or fontal)
%90,-90 x=vert, y = horiz (top down view?)
%-90,90 north up, east right

%note trail could be 'Off' 'DotsOnly' or 'All'
if DO_MOVIE
    SixDofAnimation(posPlot, quatern2rotMat(quatPlot), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'DotsOnly', ...
        'Position', [9 39 1280 768], 'View', spinMat, ...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
        'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/movie_samplePeriod) / SamplePlotFreq));
end

end

