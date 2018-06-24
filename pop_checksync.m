% pop_checksync() - estimate quality of the data synchronization between
% eye-tracker and EEG by computing the crosscorrelation between the
% (horizontal) gaze position signal and the electro-oculogram (or other
% lateral frontal EEG channels that receive strong corneoretinal artifacts). 
% Since both signals reflect the rotation of the eye ball, there should be
% virtually no time lag between both signals, that is, the cross
% correlation should peak at a lag of zero. This function computes the
% crosscorrelation between gaze and EOG, plots it, and stores the results
% in the field EEG.misc.xcorrelation_EEG_ET of the EEG structure.
%
% Usage:
%   >> EEG = pop_checksync(EEG,eye_x,heog_l,heog_r,plotfig)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data (see pop_importeyetracker)
%   gaze_x        - [one channel index], 
%                  specify channel index of the X-component (horizontal)
%                  gaze position signal of the eye tracker. If you recorded
%                  from both eyes (binocular), just pick one of the eyes 
%                  (left or right eye)
%   heog_l       - [channel index], index of EOG channel on LEFT side of 
%                  the head/face. If no proper EOG electrode was placed
%                  near the eye (bad idea!), use a lateral EEG channel that 
%                  is as close to the eye as possible
%   heog_r       - [channel index], index of EOG channel on RIGHT side of 
%                  the head/face. If no proper EOG electrode was placed
%                  near the eye (bad idea!), use a lateral EEG channel that 
%                  is as close to the eye as possible
%  plotfig       - show a plot with result of the crosscorrelation function
%
% Outputs:
%   EEG         - EEG structure with cross-correlation info added to the
%                 field EEG.misc.xcorr_EEG_ET
%
% See also: pop_importeyetracker, synchronize
%
%
% An example call of the function might look like this: 
% >> EEG = pop_checksync(EEG,65,1,2,1)
%
% In this example, the horizontal gaze position of the left eye of the 
% eye tracker was stored in channel 65 of the synchronized EEG/ET dataset. 
% The signal of the horizontal EOG electrodes, placed on the left and right 
% canthus of each eye, is stored in channels 1 and channels 2, respectively. 
% The last input indicates that EYE-EEG plots a figure that shows the 
% cross-correlation function between the ET channels and the bipolar-referenced 
% (left minus right electrode) EOG. If synchronization is good, the time 
% lag of the maximum cross-correlation should be near zero.
% 
% IMPORTANT NOTE: In order to get clean results, you should first remove
% bad or missing data in the eye tracking channels, for example due to
% eye blinks, using method pop_rej_eyecontin(). Otherwise, the cross-corr.
% function is distorted by this data. In some cases, it might also be
% necessary to remove intervals with very bad EOG recordings, but this is
% usually less of a problem.
%
% Note that this method was first applied in Dimigen et al., 2011, JEP:GEN. 
% Please cite this paper when using this method. Thanks!
%
% Author: od
% Copyright (C) 2009-2018 Olaf Dimigen, HU Berlin
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

function [EEG, com] = pop_checksync(EEG, chan_gaze_x, chan_heog_l, chan_heog_r, plotfig)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 4
        % pop up dialogue
        [chan_gaze_x, chan_heog_l, chan_heog_r, plotfig] = dlg_checksync(mfilename, EEG.chanlocs);       
    end
    
    % compute the cross-correlation function to check the synchronization
    EEG = checksync( EEG, chan_gaze_x, chan_heog_l, chan_heog_r, plotfig);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% return history command string
allArgs = vararg2str({chan_gaze_x, chan_heog_l, chan_heog_r, plotfig});
com = sprintf('EEG = %s(EEG,%s)',mfilename,allArgs);
return