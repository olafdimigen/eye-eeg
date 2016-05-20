% pop_rej_eyecontin() - reject intervals of continuous data during which
%                the eye tracker recorded extreme (out-of-range) values.
%                This function is useful for several purposes: to reject
%                data with eye blinks, to control fixation position, and to
%                reject intervals bad eye tracking data prior to saccade or
%                fixation detection.
%                EEGLAB function pop_select() is used to do the actual
%                rejection. A new boundary event is added at the position
%                of each of the resulting data breaks. To perform a similar
%                rejection for epoched data, use pop_rej_eyeepoch.
%
% Usage:
%   >> EEG = pop_rej_eyecontin(EEG,chns,minvals,maxvals,windowsize)
%
% Inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye
%                  tracking data
%   chns         - [vector of channel indices] Indices of ET channels to 
%                  check for out-of-range values
%   minvals      - [vector of length(chns)] minimum allowed values for the
%                  corresponding data channel in 'chns'. Data points
%                  with values < minvals will be rejected. minvals needs to
%                  have the same length as chns and maxvals. If you only 
%                  want to test for maxvals, minvals can be left empty [].
%   maxvals      - [vector of length(chns)] maximum allowed values for the
%                  corresponding data channel in 'chns'. Data points with
%                  with values > maxvals will be rejected. maxvals needs to
%                  have the same length as chns and minvals.  If you only 
%                  want to test for minvals, maxvals can be left empty [].
%   windowsize   - [integer] if windowsize is set to a value > 0, an 
%                  additional plusminus 'windowsize' data points will be 
%                  removed before and after each interval of out-of-range 
%                  data. We recommended to set this to a reasonable value 
%                  (e.g., 50 samples in case of 500 Hz data) because eye 
%                  blinks (usually characterized by zeros or negative 
%                  values in the eye track) are often flanked by 
%                  additional bad samples (recorded during 
%                  the partial occlusion of the pupil while the eye is 
%                  closing and re-opening) that may otherwise not exceed 
%                  the minvals/maxvals thresholds.
%
% Outputs:
%   EEG         - EEG structure with bad intervals removed. For each
%                 removed interval of bad data, a new boundary event is
%                 inserted into EEG.event.
%
%   An example call of the function might look like this:
%   >> EEG = pop_rej_eyecontin(EEG,[33 34],[1 1],[1024 768],50)
%
%   In this example, the data of channels 33 and 34, containing the
%   horizontal (33) and vertical (34) position of the left eye, are tested
%   for pixel values smaller than 1 pixel or larger than 1024 pixels
%   (channel 33) and smaller than 1 pixel or larger than 768 pixels
%   (channel 34). Screen resolution during the experiment was 1024 x 768
%   pixel. Values outside this range likely reflect eye blinks
%   or intervals during which the eye was not properly tracked. Around each
%   interval of out-of-range data an addtional plusminus 50 samples of
%   data are rejected as well. This is useful to exclude bad data samples
%   at the beginning and end of eye blinks.
%
% See also: rej_eyecontin, pop_rej_eyeepoch, pop_select
%
% Author: od
% Copyright (C) 2009-2013 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de / ulrich.reinacher.1@hu-berlin.de

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

function [EEG, com] = pop_rej_eyecontin(EEG,chans,minvals,maxvals,windowsize)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

if isempty(EEG.data)
    fprintf('\npop_rej_eyeepoched(): No dataset loaded!\n')
    return;
end

try
    if nargin < 5       
        % pop up dialogue
        [chans minvals maxvals windowsize] = dlg_rej_eyecontin(mfilename,EEG);
    end
    
    [EEG, badblocks] = rej_eyecontin(EEG,chans,minvals,maxvals,windowsize);
       
    com = sprintf('%s = pop_rej_eyecontin(%s,%s,%s,%s,%s)',inputname(1),inputname(1),vararg2str(chans),vararg2str(minvals),vararg2str(maxvals),vararg2str(windowsize));
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        % disp('User decided to cancel the selection of input.')
        % to many ways to end up here: Invalid source file (name )provided, insufficent input
        return
    else
        rethrow(err);
    end
end

return