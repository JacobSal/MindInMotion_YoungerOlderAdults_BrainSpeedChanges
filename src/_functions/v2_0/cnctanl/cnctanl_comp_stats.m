function [StatsMatrix] = cnctanl_comp_stats(varargin)
%CNCTANL_COMP_STATS Summary of this function goes here
% Create a Time-Frequency Grid from an EEGLAB/SIFT dataset. 
% For details on the Interactive Time-Frequency Grid see [1].
%
% ----------------------------------------------------------------------------------------------------------------------------------------
% Input                             Information                                                                                           
% ----------------------------------------------------------------------------------------------------------------------------------------
% ...| EEG:                         EEGLAB dataset(s)                                                                                     
%                                   This is an array of at most two EEGLAB structures.                                                    
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'MANDATORY INPUT'                                                                    
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| Conn:                        SIFT Conn object                                                                                      
%                                   This is typically stored in EEG.CAT.Conn                                                              
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'MANDATORY INPUT'                                                                    
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% [+] PlotConditionDifference:      Plot difference between selected conditions                                                           
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% .......| ConditionOrder:          Order in which to take difference                                                                     
%                                   Possible values: {''}                                                                                 
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| Stats:                       A structure containing statistics                                                                     
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| VisualizationMode:           Visualization Modes                                                                                   
%                                   Create Time-Frequency imageplots, Causality x Frequency plots (collapsing across time), Causality x   
%                                   Time plots (collapsing across frequency)                                                              
%                                   Possible values: {'TimeXFrequency', 'TimeXCausality', 'FrequencyXCausality'}                          
%                                   Default value  : 'TimeXFrequency'                                                                     
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| msubset:                     Subset of the full matrix to keep                                                                     
%                                   Lower/upper triangle ('tril'/'triu'), diagonals ('diag'), everything except diagonal ('nodiag'),      
%                                   everything ('all').                                                                                   
%                                   Possible values: {'tril', 'triu', 'diag', 'nodiag', 'all'}                                            
%                                   Default value  : 'all'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% [+] MatrixLayout:                 Select the measure and layout                                                                         
%                                   Possible values: {'Full', 'Partial'}                                                                  
%                                   Default value  : 'Full'                                                                               
%                                   Input Data Type: string                                                                               
% ....[+] Full:                                                                                                                           
%                                                                                                                                         
% .......| Estimator:               Estimator to visualize                                                                                
%                                   Possible values: {''}                                                                                 
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| ColorLimits:             Color/Y-axis scaling limits                                                                           
%                                   If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                   scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                   prctile(abs(Conn),scalar)                                                                             
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 100                                                                                  
%                                   Input Data Type: real number (double)                                                                 
% ....[+] Partial:                                                                                                                        
%                                                                                                                                         
% .......| UpperTriangle:           Estimator to render on upper triangle                                                                 
%                                   Possible values: {'none', ''}                                                                         
%                                   Default value  : 'none'                                                                               
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| UT_ColorLimits:          Color/Y-axis scaling limits for upper triangle                                                        
%                                   If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                   scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                   prctile(abs(Conn),scalar)                                                                             
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 100                                                                                  
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| LowerTriangle:           Estimator to render on upper triangle                                                                 
%                                   Possible values: {'none', ''}                                                                         
%                                   Default value  : 'none'                                                                               
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| LT_ColorLimits:          Color/Y-axis scaling limits for lower triangle                                                        
%                                   If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                   scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                   prctile(abs(Conn),scalar)                                                                             
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 100                                                                                  
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| Diagonal:                Estimator to render on diagonal                                                                       
%                                   Possible values: {'none', ''}                                                                         
%                                   Default value  : 'none'                                                                               
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| D_ColorLimits:           Color/Y-axis scaling limits for diagonal                                                              
%                                   If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                   scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                   prctile(abs(Conn),scalar)                                                                             
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 100                                                                                  
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| AllColorLimits:          Color/Y-axis scaling limits for all subplots                                                          
%                                   If set, overrides all other colorlimits options. If [min max], scale by [min max]. If scalar, and     
%                                   all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc     
%                                   maxprc] where maxprc is prctile(abs(Conn),scalar)                                                     
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| ColorLimits:                 Color/Y-axis scaling limits                                                                           
%                                   If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                   scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                   prctile(abs(Conn),scalar)                                                                             
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 100                                                                                  
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| TimesToPlot:                 [Min Max] Time range to image (sec)                                                                   
%                                   Leave blank to use all timewindows                                                                    
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| FrequenciesToPlot:           Vector of frequencies (Hz) to image                                                                   
%                                   Leave blank to use all frequencies                                                                    
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| TimeWindowsToPlot:           Time window centers (sec)                                                                             
%                                   If a vector of times, will plot a separate curve for each specified time                              
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% [+] PlotContour:                  Plot contours around significant regions                                                              
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% .......| ContourColor:            Contour Color                                                                                         
%                                   Can use any allowable Matlab color specification (see 'help ColorSpec').                              
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [0 0 0]                                                                              
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% [+] Thresholding:                 Thresholding options                                                                                  
%                                   You can choose to use statistics (passed in as 'stats' structure), or simple percentile or absolute   
%                                   thresholds.                                                                                           
%                                   Possible values: {'None', 'Statistics', 'Simple'}                                                     
%                                   Default value  : 'None'                                                                               
%                                   Input Data Type: string                                                                               
% ....[+] None:                                                                                                                                                                                                                                                              
% ....[+] Statistics:                                                                                                                     
%                                                                                                                                         
% .......| PlotConfidenceIntervals: Plot confidence intervals (if available)                                                              
%                                   Does not apply to for time-frequency images.                                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : logical(false)                                                                       
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% .......| ThresholdingMethod:      Method to use for significance masking                                                                
%                                   Possible values: {'none'}                                                                             
%                                   Default value  : 'none'                                                                               
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| AlphaSignificance:       P-value threshold for significance. e.g., 0.05 for p<0.05                                             
%                                   Possible values: [0 1]                                                                                
%                                   Default value  : 0.05                                                                                 
%                                   Input Data Type: real number (double)                                                                 
% ....[+] Simple:                     f                                                                                                    
%                                                                                                                                         
% .......| PercentileThreshold:     Percentile threshold                                                                                  
%                                   If of form [percentile, dimension], percentile is applied elementwise across the specified            
%                                   dimension.                                                                                            
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 0                                                                                    
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| AbsoluteThreshold:       Exact threshold                                                                                       
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| Baseline:                    Time range of baseline [Min Max] (sec)                                                                
%                                   Will subtract baseline from each point. Leave blank for no baseline.                                  
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| FigureHandles:               Vector of figure handles to superimpose new graph onto                                                
%                                   New figures and grid will *not* be created. Old grid will be used and new subplots overlaid           
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| Smooth2D:                    Smooth time-freq image                                                                                
%                                   This will apply nearest-neighbor interpolation.                                                       
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : logical(false)                                                                       
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% ...| XTickLabels:                 Labels for X-Tickmarks                                                                                
%                                   Must equal number of time windows                                                                     
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| YTickLabels:                 Labels for Y-Tickmarks                                                                                
%                                   Must equal number of time windows                                                                     
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| VariablesToKeep:             List of indices of channels to keep                                                                   
%                                   Can be [vector], a subset of [1:nbchan]                                                               
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| PlottingOrder:               Specify index order                                                                                   
%                                   Subset of [1:nbchan] in which to arrange columns/rows. Useful for grouping channels.                  
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| SourceMarginPlot:            What to plot on margins                                                                               
%                                   Options: 'Topoplot': plot source scalp projection. 'Dipole': plot dipole                              
%                                   Possible values: {'none', 'topoplot', 'dipole', 'customtopo'}                                         
%                                   Default value  : 'customtopo'                                                                         
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| TopoplotOptions:             Additional options (name,value) for topoplot                                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: cell array of strings (cellstr)                                                      
%                                                                                                                                         
% ...| CustomTopoMatrix:            Custom topoplot matrix                                                                                
%                                   For N channels/sources, this is a 1 X N cell array of symmetric matrices comprised the topoplot       
%                                   *surface* (not a component vector) for each channel/source. This is provided as input to              
%                                   toporeplot() if 'SourceMarginPlot' is chosen to be 'customtopo'.                                      
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% [+] DipolePlottingOptions:        Options for dipole plotting                                                                           
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| mri:                     Dipplot MRI structure                                                                                 
%                                   Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a    
%                                   path to a Matlab file containing MRI structure. Default uses MNI brain.                               
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| DipoleCoordinateFormat:  Coordinate format for dipplot                                                                         
%                                   Possible values: {'spherical', 'mni'}                                                                 
%                                   Default value  : 'mni'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% .......| ShowCortexMesh:          Show cortex surface instead of MRI volume                                                             
%                                   Only valid if EEG.dipfit.surfmesh and EEG.dipfit.reducedMesh are present. These are structures        
%                                   containing fields .faces and .vertices                                                                
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : logical(false)                                                                       
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% .......| ColorROIs:               Color ROIs                                                                                            
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : logical(false)                                                                       
%                                   Input Data Type: boolean                                                                              
%                                                                                                                                         
% .......| DipoleSize:              Dipole sphere size                                                                                    
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 80                                                                                   
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| DipplotOptions:          Additional dipplot options                                                                            
%                                   Cell array of <'name',value> pairs of additional options for dipplot (see 'doc dipplot')              
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : '{}'                                                                                 
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% .......| row_view:                View angle for row marginals                                                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [1 0 0]                                                                              
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% .......| col_view:                View angle for column marginals                                                                       
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [0 0 1]                                                                              
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| NodeLabels:                  List of labels for each node. e.g., {'Node1','Node2',...}                                             
%                                   Leave blank to use defaults.                                                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| FrequencyMarkers:            Vector of frequencies (Hz) at which to draw horizontal lines                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| FrequencyMarkerColor:        Coloring for frequency markers                                                                        
%                                   If an [1 x 3] array of RBG values, then color all lines using this color. If an [N x 3] matrix of     
%                                   RBG values, then color the kth line with the colorspec from the kth row. If empty then cycle          
%                                   through colorlist                                                                                     
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| EventMarkers:                Event marker time and style                                                                           
%                                   Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. {    
%                                   { 0.2 'y' ':' 2} { 1.5 'r' ':' 2}} will render two dotted-line event makers, yellow at 200 ms and     
%                                   red at 1500 ms                                                                                        
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : {{0, 'r', ':', 2}}                                                                   
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| FrequencyScale:              Make the y-scale logarithmic or linear                                                                
%                                   Possible values: {'linear', 'log'}                                                                    
%                                   Default value  : 'linear'                                                                             
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| Transform:                   transform the data (logarithmically or other)                                                         
%                                   Possible values: {'log', 'linear', ''}                                                                
%                                   Default value  : 'linear'                                                                             
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| YTickLabelLoc:               Y-tick label location                                                                                 
%                                   Possible values: {'left', 'right', 'both'}                                                            
%                                   Default value  : 'right'                                                                              
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| TitleString:                 Figure title string                                                                                   
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'n/a'                                                                                
%                                   Input Data Type: string                                                                               
%                                                                                                                                         
% ...| TitleFontSize:               Title Font Size                                                                                       
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 12                                                                                   
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| AxesFontSize:                Axes Font Size                                                                                        
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 11                                                                                   
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| TextColor:                   Text color                                                                                            
%                                   See 'doc ColorSpec'.                                                                                  
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [1 1 1]                                                                              
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| LineColor:                   Linecolor for lineplots                                                                               
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [1 1 1]                                                                              
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| PatchColor:                  FaceColor for shaded regions                                                                          
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [1 1 1]                                                                              
%                                   Input Data Type: real number (double)                                                                 
%                                                                                                                                         
% ...| Colormap:                    Colormap                                                                                              
%                                   Matlab expression denoting colormap to use (e.g., 'jet(64)'). See 'help colormap'.                    
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : 'jet(300)'                                                                           
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| BackgroundColor:             Background Color                                                                                      
%                                   See 'doc ColorSpec'.                                                                                  
%                                   Possible values: 'Unrestricted'                                                                       
%                                   Default value  : [0 0 0]                                                                              
%                                   Input Data Type: any evaluable Matlab expression.                                                     
%                                                                                                                                         
% ...| TFCellColorScheme:           Color scheme for TimeFreqCell popout                                                                  
%                                   Possible values: {'black', 'white', 'eeglab'}                                                         
%                                   Default value  : 'black'                                                                              
%                                   Input Data Type: string  
%
% ----------------------------------------------------------------------------------------------------------------------------------------
% Output                             Information                                                                                           
% ----------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                   
%
%       StatsMatrix:                
%
%
% See Also: pop_vis_TimeFreqGrid()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen, 2010, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% extract some stuff from inputs for arg defaults
Conn = arg_extract(varargin,'Conn',2);

if ~ischar(Conn) && ~isempty(Conn)
    Conn = Conn(1);
    ConnNames   = hlp_getConnMethodNames(Conn);
    conndef     = ConnNames{1};
    freqrange   = [Conn.freqs(1) Conn.freqs(end)];
    freqdef     = Conn.freqs; %['[' num2str(freqrange(1)) ':' num2str(Conn.freqs(2)-Conn.freqs(1)) ':' num2str(freqrange(end)) ']'];
    timerange   = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    timedef     = [timerange(1) timerange(end)];
    clear Conn;
else
    ConnNames = {''};
    conndef = '';
    [freqrange, freqdef, timerange, timedef] = deal([]);
end

% get some defaults from ALLEEG
ALLEEG = arg_extract(varargin,{'ALLEEG','EEG'},1);
[MyComponentNames MyChannelNames] = deal([]);
if ~ischar(ALLEEG) && ~isempty(ALLEEG)
    if isfield(ALLEEG(1).CAT,'curComponentNames') && ~isempty(ALLEEG(1).CAT.curComponentNames)
        MyComponentNames = ALLEEG(1).CAT.curComponentNames;
    else
        MyComponentNames = ALLEEG(1).CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    if isfield(ALLEEG(1),'chanlocs') && ~isempty(ALLEEG(1).chanlocs)
        MyChannelNames = {ALLEEG(1).chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:ALLEEG(1).nbchan)'))');
    end
    
    % set allowable options for marginal plots
    sourceMarginOptions = {'none'};
    if isfield(ALLEEG(1),'chanlocs') && isfield(ALLEEG(1).chanlocs,'X')
        sourceMarginOptions = [sourceMarginOptions 'topoplot'];
    end
    if isfield(ALLEEG(1),'dipfit') && ~isempty(ALLEEG(1).dipfit)
        sourceMarginOptions = [sourceMarginOptions 'dipole'];
    end
    if ~isempty(arg_extract(varargin,{'customTopoMatrix'},[],[]))
    	sourceMarginOptions = [sourceMarginOptions 'customtopo'];
    end
    
    % get the condition names
    for cnd=1:length(ALLEEG)
        conditionNames{cnd} = ALLEEG(cnd).condition;
        setNames{cnd}       = ALLEEG(cnd).setname;
        fileNames{cnd}      = ALLEEG(cnd).filename;
    end
    
    % set up condition difference order defaults
    if length(ALLEEG)==2
        if ~any(cellfun(@isempty,conditionNames))
            % use condition names for labeling
            CondDiffOrderDefaults = ...
                {sprintf('%s-%s',conditionNames{1},conditionNames{2}), ...
                sprintf('%s-%s',conditionNames{2},conditionNames{1})};
            CondLabels = conditionNames;
        elseif ~any(cellfun(@isempty,setNames))
            % use set names for labeling
            CondDiffOrderDefaults = ...
                {sprintf('%s-%s',setNames{1},setNames{2}), ...
                sprintf('%s-%s',setNames{2},setNames{1})};
            CondLabels = setNames;
        elseif ~any(cellfun(@isempty,fileNames))
            % use set filenames for labeling
            CondDiffOrderDefaults = ...
                {sprintf('%s-%s',fileNames{1},fileNames{2}), ...
                sprintf('%s-%s',fileNames{2},fileNames{1})};
            CondLabels = fileNames;
        else
            % use set numbers for labeling
            CondDiffOrderDefaults = {'Set 1 - Set 2','Set 2 - Set 1'};
            CondLabels = {'Set 1','Set 2'};
        end
    else
        CondDiffOrderDefaults = {''};
        if ~isempty(conditionNames{1})
            CondLabels = conditionNames;
        elseif ~isempty(setNames{1})
            CondLabels = setNames;
        elseif ~isempty(fileNames{1})
            CondLabels = fileNames;
        else
            CondLabels = {'Set 1'};
        end
    end
    
else
    CondDiffOrderDefaults = {''};
    sourceMarginOptions = {'none','topoplot','dipole','customtopo'};
end


% determine whether statistics are present
% (PROBLEM 'Stats' argument is case-sensitive)
if isfield(varargin{1},'icaact') && length(varargin)==2
    stats = [];
elseif isfield(varargin{1},'icaact')
    stats = arg_extract(varargin(3:end),'Stats',[],[]);
else
    stats = arg_extract(varargin,'Stats',[],[]);
end

if isempty(stats)
    usestatsdef = [];  % false
    StatThreshMethods = {'none'};
    def_alpha = 0.05;
else
    usestatsdef = {};  % true
    methods = intersect_bc(fieldnames(stats),ConnNames);
    StatThreshMethods = intersect_bc(fieldnames(stats.(methods{1})),{'pval','thresh','logical'});
    StatThreshMethods = [StatThreshMethods 'none'];
    def_alpha = stats.alpha;
end

% clear stats ALLEEG;

% setup the argument list
% -----------------------------------------------------
g = arg_define([0 2],varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory,[],'EEGLAB dataset(s). This is an array of at most two EEGLAB structures.','type','expression'),...
    arg_norep({'Conn'},mandatory,[],'SIFT Conn object. This is typically stored in EEG.CAT.Conn','type','expression'),...
    arg_subtoggle({'plotCondDiff','PlotConditionDifference'},{}, ...
    {...
    arg({'condOrder','ConditionOrder'},CondDiffOrderDefaults{1},CondDiffOrderDefaults,'Order in which to take difference.') ...
    }, 'Plot difference between selected conditions','cat','DisplayProperties'), ...
    arg_norep({'stats','Stats'},[],[],'A structure containing statistics.','type','expression'), ...
    arg({'vismode','VisualizationMode'},'',{'','TimeXFrequency','TimeXCausality','FrequencyXCausality'},'Visualization Modes. Create Time-Frequency imageplots, Causality x Frequency plots (collapsing across time), Causality x Time plots (collapsing across frequency)'), ...
    arg_nogui({'msubset'},'all',{'tril','triu','diag','nodiag','all'},'Subset of the full matrix to keep. Lower/upper triangle (''tril''/''triu''), diagonals (''diag''), everything except diagonal (''nodiag''), everything (''all'').'), ...
    arg_subswitch({'MatrixLayout'},'Full', ...
    {'Full', ...
    { ...
    arg({'estimator','Estimator'},ConnNames{1},ConnNames,'Estimator to visualize','shape','row') ...
    arg({'clim','ColorLimits'},100,[],'Color/Y-axis scaling limits. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    } ...
    'Partial', ...
    {...
    arg({'triu','UpperTriangle'},ConnNames{1},['none' ConnNames],'Estimator to render on upper triangle.','shape','row'), ...
    arg({'ut_clim','UT_ColorLimits'},100,[],'Color/Y-axis scaling limits for upper triangle. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    arg({'tril','LowerTriangle'},ConnNames{1},['none' ConnNames],'Estimator to render on upper triangle.','shape','row'), ...
    arg({'lt_clim','LT_ColorLimits'},100,[],'Color/Y-axis scaling limits for lower triangle. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    arg({'diag','Diagonal'},ConnNames{1},['none' ConnNames],'Estimator to render on diagonal.','shape','row') ...
    arg({'d_clim','D_ColorLimits'},100,[],'Color/Y-axis scaling limits for diagonal. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    arg({'clim','AllColorLimits','ColorLimits'},[],[],'Color/Y-axis scaling limits for all subplots. If set, overrides all other colorlimits options. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    }, ...
    },'Select the measure and layout','cat','DisplayProperties'), ...
    arg_nogui({'clim','ColorLimits'},100,[],'Color/Y-axis scaling limits. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    arg({'timeRange','TimesToPlot'},timedef,[],'[Min Max] Time range to image (sec). Leave blank to use all timewindows','shape','row','type','denserealdouble','cat','DisplayProperties'), ...
    arg({'freqValues','FrequenciesToPlot'},freqdef,[],'Vector of frequencies (Hz) to image. Leave blank to use all frequencies','type','expression','shape','row','cat','DisplayProperties'), ...
    arg_nogui({'windows','TimeWindowsToPlot'},[],[],'Time window centers (sec). If a vector of times, will plot a separate curve for each specified time','shape','row','cat','DisplayProperties'),...
    arg_subtoggle({'pcontour','PlotContour'},[], ...
    {...
    arg({'contourcolor','ContourColor'},[0 0 0],[],'Contour Color. Can use any allowable Matlab color specification (see ''help ColorSpec'').','shape','row','type','expression','cat','DisplayProperties') ...
    }, 'Plot contours around significant regions','cat','DisplayProperties'), ...
    arg_subswitch({'thresholding','Thresholding'},'None', ...
    {'None' ...
    { ...
    arg_norep({'dummy1'},[],[],'dummy') ...
    }, ...
    'Statistics' ...
    {...
    arg({'plotci','PlotConfidenceIntervals'},false,[],'Plot confidence intervals (if available). Does not apply to for time-frequency images.'), ...
    arg({'sigthreshmethod','ThresholdingMethod'},StatThreshMethods{1},StatThreshMethods,'Method to use for significance masking') ...
    arg({'alpha','AlphaSignificance'},def_alpha,[0 1],'P-value threshold for significance. e.g., 0.05 for p<0.05') ...
    }, ...
    'Simple' ...
    {...
    arg({'prcthresh','PercentileThreshold'},0,[],'Percentile threshold. If of form [percentile, dimension], percentile is applied elementwise across the specified dimension.','type','denserealdouble','shape','row','cat','Thresholding'), ...
    arg({'absthresh','AbsoluteThreshold'},[],[],'Exact threshold.','cat','Thresholding') ...
    } ...
    }, 'Thresholding options. You can choose to use statistics (passed in as ''stats'' structure), or simple percentile or absolute thresholds.','cat','Thresholding'), ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point. Leave blank for no baseline.','shape','row','type','denserealdouble','cat','DataProcessing'), ...
    arg_nogui({'fighandles','FigureHandles'},[],[],'Vector of figure handles to superimpose new graph onto. New figures and grid will *not* be created. Old grid will be used and new subplots overlaid'), ...
    arg({'smooth','Smooth2D'},false,[],'Smooth time-freq image. This will apply nearest-neighbor interpolation.','cat','DataProcessing'), ...
    arg_nogui({'xord','XTickLabels'},[],[],'Labels for X-Tickmarks. Must equal number of time windows or XTickLoc','type','expression','cat','DisplayProperties'), ...
    arg_nogui({'yord','YTickLabels'},[],[],'Labels for Y-Tickmarks. Must equal number of time windows or YTickLoc','type','expression','cat','DisplayProperties'), ...
    arg_nogui({'xloc','XTickLoc'},[],[],'Locations for X-Tickmarks.','type','expression','cat','DisplayProperties'), ...
    arg_nogui({'yloc','YTickLoc'},[],[],'Locations for Y-Tickmarks.','type','expression','cat','DisplayProperties'), ...
    arg_norep({'channels','VariablesToKeep'},[],[],'*deprecated* List of indices of channels to keep. Can be [vector], a subset of [1:nbchan]'), ...
    arg({'plotorder','PlottingOrder'},[],[],'Specify index order. Subset of [1:nbchan] in which to arrange columns/rows. Useful for grouping channels.','cat','DisplayProperties'), ...
    arg({'topoplot','SourceMarginPlot'},sourceMarginOptions{end},sourceMarginOptions,'What to plot on margins. Options: ''Topoplot'': plot source scalp projection. ''Dipole'': plot dipole','cat','DisplayProperties'), ...
    arg_nogui({'topoplot_opts','TopoplotOptions'},{},[],'Additional options (name,value) for topoplot','type','cellstr'), ...    
    arg_nogui({'customTopoMatrix','CustomTopoMatrix'},[],[],'Custom topoplot matrix. For N channels/sources, this is a 1 X N cell array of symmetric matrices comprised the topoplot *surface* (not a component vector) for each channel/source. This is provided as input to toporeplot() if ''SourceMarginPlot'' is chosen to be ''customtopo''.','shape','matrix','cat','DisplayProperties'), ...
    arg_sub({'dipplot','DipolePlottingOptions'},[], ...
    { ...
    arg_nogui({'mri'},'',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
    arg({'coordformat','DipoleCoordinateFormat'},'mni',{'spherical','mni'},'Coordinate format for dipplot','type','char','shape','row'), ...
    arg({'showCortexMesh','ShowCortexMesh'},isstruct(ALLEEG(1)) && ~isempty(ALLEEG(1).dipfit) && isfield(ALLEEG(1).dipfit,'model') && isfield(ALLEEG(1).dipfit.model,'meshVertices'),[],'Show cortex surface instead of MRI volume. Only valid if EEG.dipfit.surfmesh and EEG.dipfit.reducedMesh are present. These are structures containing fields .faces and .vertices'), ...
    arg({'colorROIs','ColorROIs'},isstruct(ALLEEG(1)) && isfield(ALLEEG(1).dipfit,'surfmesh'),[],'Color ROIs.'), ...
    arg({'dipsize','DipoleSize'},80,[],'Dipole sphere size'), ...
    arg_nogui({'dipplotopt','DipplotOptions'},'{}','','Additional dipplot options. Cell array of <''name'',value> pairs of additional options for dipplot (see ''doc dipplot'')','type','expression','shape','row') ...
    arg({'row_view'},[1 0 0],[],'View angle for row marginals'), ...
    arg({'col_view'},[0 0 1],[],'View angle for column marginals'), ...
    },'Options for dipole plotting'), ...
    arg({'nodelabels','NodeLabels'},MyComponentNames,[],'List of labels for each node. e.g., {''Node1'',''Node2'',...}. Leave blank to use defaults.','type','expression','cat','DisplayProperties'),...
    arg({'foilines','FrequencyMarkers'},[],[],'Vector of frequencies (Hz) at which to draw horizontal lines','cat','FrequencyMarkers'), ...
    arg({'foilinecolor','FrequencyMarkerColor'},[],[],'Coloring for frequency markers. If an [1 x 3] array of RBG values, then color all lines using this color. If an [N x 3] matrix of RBG values, then color the kth line with the colorspec from the kth row. If empty then cycle through colorlist','shape','matrix','cat','FrequencyMarkers'), ...
    arg({'events','EventMarkers'},{{0 'r' ':' 2}},[],'Event marker time and style. Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. { { 0.2 ''y'' '':'' 2} { 1.5 ''r'' '':'' 2}} will render two dotted-line event makers, yellow at 200 ms and red at 1500 ms','type','expression','shape','row','cat','DisplayProperties'), ...
    arg({'freqscale','FrequencyScale'},'linear',{'linear','log'},'Make the y-scale logarithmic or linear','cat','DisplayProperties'), ...
    arg_nogui({'transform','Transform'},'linear',{'log','linear',''},'transform the data (logarithmically or other)'), ...
    arg({'yTickLoc','YTickLabelLoc'},'right',{'left','right','both'},'Y-tick label location.','cat','TextAndFont'), ...
    arg({'titleString','TitleString'},[],[],'Figure title string','type','char','shape','row','cat','TextAndFont'), ...
    arg({'titleFontSize','TitleFontSize'},12,[],'Title Font Size','cat','TextAndFont'), ...
    arg({'axesFontSize','AxesFontSize'},11,[],'Axes Font Size','cat','TextAndFont'), ...
    arg({'textColor','TextColor'},[1 1 1],[],'Text color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg({'linecolor','LineColor'},[1 1 1],[],'Linecolor for lineplots','shape','row','type','denserealdouble','cat','TextAndFont'), ...
    arg({'patchcolor','PatchColor'},[1 1 1],[],'FaceColor for shaded regions','shape','row','type','denserealdouble','cat','TextAndFont'), ...
    arg({'colormap','Colormap'},'jet(300)',[],'Colormap. Matlab expression denoting colormap to use (e.g., ''jet(64)''). See ''help colormap''.','type','expression','cat','DisplayProperties'), ...
    arg({'backgroundColor','BackgroundColor'},[0 0 0],[],'Background Color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg({'colorscheme','TFCellColorScheme','ColorScheme'},'eeglab',{'black','white','eeglab'},'Color scheme for TimeFreqCell popout'), ...
    arg_norep({'report_args'},[],[],'Need this to allow recursive calls...') ...
    );

%     arg_sub({'subplotargs','SubplotExpansionProperties'},[],@vis_TimeFreqCell,'Additional arguments for subplot callback function.','cat','SubplotExpansion'), ...

% Commit ALLEEG and Conn variables to workspace
[data g] = hlp_splitstruct(g,{'ALLEEG','Conn'});
arg_toworkspace(data);
clear data;

if length(Conn)>1
    % we are creating difference plot
    ConnMatrix  = Conn(1).(CEstimator) - Conn(2).(CEstimator);
    TwoSidedThresholding = true;
else
    ConnMatrix  = Conn.(CEstimator);
end

if ~isempty(g.baseline)
    % subtract baseline from connectivity matrix
    ConnMatrix = hlp_rmbaseline(ConnMatrix,g.baseline,erWinCenterTimes);
    TwoSidedThresholding = true;
end

if ~strcmpi(g.msubset,'all')
    % only plot a subset of the full TF Grid matrix
    switch lower(g.msubset)
        case 'nodiag'
            % NaN diagonal elements of nch x nch grid
            Dd = single(~eye(nch));
            Dd(Dd==0)=nan;
        case 'diag'
            % NaN off-diagonal elements of nch x nch grid
            Dd = eye(nch);
            Dd(Dd==0)=nan;
        case 'tril'
            % NaN everything except lower triangle of nch x nch grid
            Dd = tril(ones(nch),-1);
            Dd(Dd==0)=nan;
        case 'triu'
            % NaN everything except upper triangle of nch x nch grid
            Dd = triu(ones(nch),1);
            Dd(Dd==0)=nan;
    end
    ConnMatrix=ConnMatrix.*(repmat(Dd,[1,1,size(ConnMatrix,3),size(ConnMatrix,4)]));
%     if willPlotStatCI(g,CEstimator)
%         for k=1:2
%             g.stats.(CEstimator).ci(k,:,:,:,:) ...
%                 = squeeze(g.stats.(CEstimator).ci(k,:,:,:,:) ...
%                 .*(repmat(Dd,[1,size(g.stats.(CEstimator).ci,4),size(g.stats.(CEstimator).ci,5)])));
%         end
%     end
end

% ----------------------------------------------------------------
% | Apply statistics and thresholding
% ----------------------------------------------------------------

if strcmpi(g.thresholding.arg_selection,'simple')
    if any(ConnMatrix(:)<0), TwoSidedThresholding = true; end
    
    if ~isempty(g.thresholding.prcthresh)
        % percentile Thresholding
        
        if length(g.thresholding.prcthresh)>1
            % get percentiles across a specified dimension
            dim=g.thresholding.prcthresh(2);
            sz=size(ConnMatrix); nd=length(sz);
            odims = setdiff_bc(1:nd,dim);
            Cntmp = permute(ConnMatrix,[odims(1) dim odims(2:end)]);    % put dim in second dim (for reshape)
            Cntmp = reshape(Cntmp,[prod(sz(odims)) sz(dim)]);           % we'll take percentiles for ea. col
            
            if TwoSidedThresholding
                % apply two-sided thresholding
                StatsMatrix(2,:)= prctile(Cntmp,g.thresholding.prcthresh(1),1);             % upper limit
                StatsMatrix(1,:)= prctile(Cntmp,100-g.thresholding.prcthresh(1),1);         % lower limit
                
                % expand StatsMatrix to same size as ConnMatrix
                SMtmp(2,:,:,:,:,:) = repmat(StatsMatrix(2,:),[sz(odims(1)) 1 sz(odims(2:end))]);
                SMtmp(1,:,:,:,:,:) = repmat(StatsMatrix(1,:),[sz(odims(1)) 1 sz(odims(2:end))]);
                StatsMatrix = ipermute(SMtmp,[1 1+[odims(1) dim odims(2:end)]]);
                clear SMtmp
            else
                % apply single-sided thresholding
                StatsMatrix = prctile(Cntmp,g.thresholding.prcthresh(1),1);
                
                % expand StatsMatrix to same size as ConnMatrix
                StatsMatrix = repmat(StatsMatrix,[sz(odims(1)) 1 sz(odims(2:end))]);
                StatsMatrix = ipermute(StatsMatrix,[odims(1) dim odims(2:end)]);
            end
            
            clear Cntmp
            
        else
            % get percentiles of complete data matrix
            if TwoSidedThresholding
                StatsMatrix(2,:)= prctile(ConnMatrix(:),g.thresholding.prcthresh(1),1);     % upper limit
                StatsMatrix(1,:)= prctile(ConnMatrix(:),100-g.thresholding.prcthresh(1),1); % lower limit
            else
                StatsMatrix = prctile(ConnMatrix(:),g.thresholding.prcthresh);
            end
        end
    end
    
    % additonally, use absolute thresholding
    if ~isempty(g.thresholding.absthresh)
        % use scalar threshold
        StatsMatrix = g.thresholding.absthresh; %*ones(size(ConnMatrix));
    end
    
elseif willPlotStats(g,CEstimator)
    
    if ~isfield(g.stats,'tail')
        if any(ConnMatrix(:)<0), TwoSidedThresholding = true; end
    else
        switch g.stats.tail
            case {'both'}
                TwoSidedThresholding = true;
            otherwise
                TwoSidedThresholding = false;
        end
    end
    
    % Use Statistics Structure for thresholding
    if nargin>3 && isfield(g.stats,CEstimator)
        if length(g.stats)>1
            fprintf('WARNING: Stats contains more than one structure. Taking the first one...\n');
            g.stats = g.stats(1);
        end
        if isstruct(g.stats.(CEstimator))
            if ~isfield(g.stats.(CEstimator),g.thresholding.sigthreshmethod)
                fprintf('ERROR: %s is not a field of g.stats.%s\n',g.thresholding.sigthreshmethod,CEstimator);
                return;
            end
            StatsMatrix = g.stats.(CEstimator).(g.thresholding.sigthreshmethod);
        else
            StatsMatrix = g.stats.(CEstimator);
        end
    else
        StatsMatrix = [];
    end
    
    if strcmpi(g.thresholding.sigthreshmethod,'pval')
        StatsMatrix = StatsMatrix <= g.thresholding.alpha;
    end
else
    TwoSidedThresholding = any(ConnMatrix(:)<0);
end

% ---------------------------------------------------------------
% | Apply significance mask
% ---------------------------------------------------------------

% preserve the original connectivity matrix
OrigConnMatrix = ConnMatrix;

if ~isempty(StatsMatrix) && isempty(g.windows) && ~g.pcontour.arg_selection
    
    if TwoSidedThresholding % two-sided thresholds (x < lothresh | x > hithresh = nan)
        if isscalar(StatsMatrix) % && ~g.pcontour.arg_selection
            % uniform two-sided threshold
            ConnMatrix(abs(ConnMatrix) < abs(StatsMatrix)) = 0;
        elseif length(StatsMatrix)==2  %&& ~g.pcontour.arg_selection
            % uniform two-sided threshold
            ConnMatrix(ConnMatrix > StatsMatrix(1) & ConnMatrix < StatsMatrix(2)) = 0;
        elseif isequal(size(StatsMatrix),size(ConnMatrix)) && islogical(StatsMatrix)
            % logical thresholding (e.g., p-value)
            ConnMatrix(~StatsMatrix) = 0;
        elseif isequal(size(StatsMatrix),[2 size(ConnMatrix)])
            % two-sided numeric thresholding
            ConnMatrix(ConnMatrix > squeeze(StatsMatrix(1,:,:,:,:))  ...
                & ConnMatrix < squeeze(StatsMatrix(2,:,:,:,:))) = 0;
        else
            error('unknown statistical thresholding paradigm');
        end
        
    else % single-sided thresholding (x < thresh = 0)
        if isscalar(StatsMatrix) % && ~g.pcontour.arg_selection
            ConnMatrix(ConnMatrix < StatsMatrix) = 0;
        elseif isvector(StatsMatrix) && ndims(ConnMatrix)>3
            %                 sz=size(StatsMatrix);
            %                 % expand vector thresh to dims of ConnMatrix
            %                 StatsMatrix = repmat(StatsMatrix,fastif(sz(1)==1,[length(freqs) 1],[1 length(erWinCenterTimes)]));
        elseif isequal(size(StatsMatrix),size(ConnMatrix))
            if islogical(StatsMatrix)
                % logical thresholding (e.g., p-value)
                ConnMatrix(~StatsMatrix) = 0;
            else
                % single-sided numeric thresholding
                ConnMatrix(ConnMatrix < StatsMatrix) = 0;
            end
        else
            error('unknown statistical thresholding paradigm');
        end
    end
    
end


if ~isempty(g.windows)
    % Instead of Time-Frequency plots, user wishes to plot causal
    % spectra for individual selected window(s)
    g.windows = unique_bc(g.windows);
    windowIndex = getindex(erWinCenterTimes,g.windows);
    ConnMatrix=ConnMatrix(:,:,:,windowIndex);
    if ~isempty(StatsMatrix) && ~isscalar(StatsMatrix)
        StatsMatrix = StatsMatrix(:,:,:,windowIndex);
    end
    
    OrigConnMatrix = OrigConnMatrix(:,:,:,windowIndex);
    
    % select window for confidence interval
    if ~isempty(g.stats) && isfield(g.stats.(CEstimator),'ci') && g.thresholding.plotci
        g.stats.(CEstimator).ci = g.stats.(CEstimator).ci(:,:,:,:,windowIndex);
    end
end
end

