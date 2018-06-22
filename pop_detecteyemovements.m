% detecteyemovements() - detect saccades & fixations in eye tracking data. 
%                Saccade detection is based on the algorithm by 
%                Engbert & Mergenthaler (2006). Saccades are defined as 
%                (monocular or binocular) outliers in 2D velocity space.
%                Velocity thresholds for saccade detection are determined 
%                adaptively as a multiple of the (median-based) SD of all 
%                data samples in the epoch. Fixations are defined as 
%                intervals in-between saccades. Eye movements can be added 
%                as new events to EEGLAB's event structure. For various 
%                other options, see below.
%
% Usage:
%   >> EEG = pop_detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,...
%            degperpixel,smooth,globalthresh,clusterdist,clustermode,
%            plotfig,writesac,writefix)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data (see pop_importeyetracker)
%   left_eye_xy  - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of left eye. Leave empty []
%                  if the left eye was not recorded.
%   right_eye_xy - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of right eye. Leave empty []
%                  if the right eye was not recorded.
%   vfac         - [double] velocity factor ("lambda") to determine 
%                  the velocity threshold for saccade detection 
%                  (cf. Engbert & Mergenthaler, 2006)
%   mindur       - [integer] minimum saccade duration (in samples)
%                  (cf. Engbert & Mergenthaler, 2006)
%   degperpixel  - [double] visual angle of one screen pixel
%                  if this value is left empty [], saccade characteristics 
%                  are reported in the original data metric (pixel?) 
%                  instead of in degrees of visual angle
%   smooth       - [0/1] if set to 1, the raw data is smoothed to suppress
%                  noise. Recommended for high native ET sampling rates.
%   globalthresh - [0/1]. Use the same thresholds for all epochs? 
%                  0: Adaptive velocity thresholds are computed 
%                  individually for each data epoch. 
%                  1: Adaptive velocity thresholds are first computed for 
%                  each epoch, but then the mean thresholds are applied to 
%                  each epochs (i.e. same detection parameters are used for
%                  all epochs). Setting is irrelevant if the input data is 
%                  still continuous (= only one data epoch).
%   clusterdist  - [integer] value in sampling points that defines the
%                  minimum allowed fixation duration between two saccades. 
%                  If the off- and onsets of two temp. adjacent sacc. are 
%                  closer together than 'clusterdist' samples, these 
%                  saccades are regarded as a "cluster" and treated 
%                  according to the 'clustermode' setting (see below).
%                  clusterdist is irrelevant if clustermode == 1.
%   clustermode  - [1,2,3,4]. Integer between 1 and 4.
%                  1: keep all saccades, do nothing
%                  2: keep only first saccade of each cluster
%                  3: keep only largest sacc. of each cluster
%                  4: combine all movements into one (longer) saccade
%                     this new saccade is defined as the movement that 
%                     occurs between the onset of the 1st saccade in the
%                     cluster and the offset of the last sacc. in cluster.
%                     WARNING: CLUSTERMODE 4 is experimental and untested!  
%   plotfig      - [0/1] Show a figure with eye movement properties?
%                  0: do not plot a figure. 
%                  1: plot a figure displaying properties of detected 
%                  saccades & fixations
%   writesac     - [0/1]: Add saccades to EEG.event?
%                  0: detect saccades, but do not store them in EEG.event. 
%                  1: add detected saccades as new events to EEG.event.
%   writefix     - [0/1]: Add fixations to EEG.event?
%                  0: detect fixations, but do not add them to EEG.event.
%                  1: add detected fixations as new events to EEG.event.
%                  Note: It is recommended to first test the parameters of
%                  saccade/fixation detection without adding events.
%                  For this, set writesac and writefix to 0.
%
% Outputs:
%   EEG         - EEG structure. If writesac or writefix were set to 1, 
%                 the EEG structure (EEG.event/EEG.urevent/EEG.epoch)
%                 will contain additional "saccade" and "fixation" events
%                 with their respective properties
%
% See also: vecvel, velthresh, microsacc_plugin, binsacc, saccpar, 
%           mergesacc, addevents
%
%
% An example call of the function might look like this: 
% >> EEG = pop_detecteyemovements(EEG,[],[33 34],6,4,0.037,1,0,25,4,1,1,0)
%
% In this example, the eye position data for the right eye is stored in 
% channels 33 (horiz.) and 34 (vertical). The left eye was not recorded. 
% The velocity threshold is set to 6 times the (median-based) 
% SD of all velocity samples in the epoch. The minimum duration of
% saccades to be detected is 4 samples. In the experiment, one screen 
% pixel corresponded to 0.037 degrees of visual angle. 
% The raw data is smoothed prior to saccade detection (smooth: 1). 
% Adaptive velocity thresholds (X and Y-threshold for each eye) are 
% determined individually for each data epoch (globalthresh: 0). For saccades 
% separated by fixations of less than 25 samples, only the first saccade 
% is kept (clusterdist: 25, clustermode: 2). A figure with the 
% results is plotted. Detected saccades are stored as new events in 
% EEG.event, but fixations are not stored.
% 
% The eye movement detection is based on:
%
% Engbert, R., & Kliegl, R. (2003). Microsaccades uncover the orientation
% of covert attention. Vision Research, Vol. 43, 1035-1045
%
% ...as well as...
%
% Engbert, R., & Mergenthaler, K. (2006). Microsaccades are triggered by 
% low retinal image slip, PNAS, Vol. 103 (18), 7192-7197
%
% Author: od
% Copyright (C) 2009-2018 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de 

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
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

function [EEG, com] = pop_detecteyemovements(EEG, left_eye_xy, right_eye_xy, vfac, mindur, degperpixel, smooth, globalthresh, clusterdist, clustermode, plotfig, writesac, writefix)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 11
        % pop up dialogue
        [left_eye_xy, right_eye_xy, vfac, mindur, degperpixel, smooth, globalthresh, clusterdist, clustermode, plotfig, writesac, writefix] = dlg_detecteyemovements(mfilename, EEG.chanlocs, EEG.srate);       
    end
    
    % detect eye movements
    EEG = detecteyemovements(EEG, left_eye_xy, right_eye_xy, vfac, mindur, degperpixel, smooth, globalthresh, clusterdist, clustermode, plotfig, writesac, writefix);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% return history command string
allArgs = vararg2str({left_eye_xy, right_eye_xy, vfac, mindur, degperpixel, smooth, globalthresh, clusterdist, clustermode, plotfig, writesac, writefix});
com = sprintf('EEG = %s(EEG,%s)',mfilename,allArgs);
return