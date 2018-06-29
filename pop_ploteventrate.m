% pop_ploteventrate() - plot rate of (micro)saccades or fixations relative to
%                     the epoch time-locking event (usually stimulus onset)
%
% Usage:
%   >> [times rate] = pop_ploteventrate(EEG)
%
% Inputs:
%  EEG         - [string] EEG struct. EEG must already contain detected
%                saccades and/or fixations in EEG.event, either imported from
%                the (Eyelink) raw data (>> pop_importeyetracker()) or
%                detected later using function pop_detecteyemovements()
%
% rate_event    - [string] name of event (in EEG.event). This
%                could be "saccade" (default) or "fixation"
%
% Outputs:
%            - figure showing average rate (events per second) of the
%              "rate_event" relative to time-locking event of the epoch
%
% Note: this function is for epoched data, but will also produce a plot for
%       continuous data, although a rate plot will usually not make much
%       sense for continuous data (= 1 long data epoch).
%       The functions requires that saccade or fixation events were already
%       imported or detected using pop_detecteyemovements()
%
% See also: pop_detecteyemovements, pop_ploteyemovments
%
% Example: A call to this function might look like this:
%
% >> pop_ploteventrate(EEG) % --> pops dialogue window
% >> pop_ploteventrate(EEG,'saccade')
% >> pop_ploteventrate(EEG,'fixation')
% >> [times, rate] = pop_ploteventrate(EEG,'saccade') % to make your own plots
%
% Note: if you plot the rate of fixations in epoched data, it will be
% very high in the first sample of the epoch, since this is usually
% automatically the start of a fixation, if fixations were detected in
% epoched data (and if the epoch does not start with a saccade)
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

function com = pop_ploteventrate(EEG,rate_event)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 2
        % pop dialogue
        rate_event = dlg_ploteventrate(mfilename,EEG);
    end
    
    ploteventrate(EEG,rate_event);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end


%% Show pop_resample() EEG.event.duration warning if EEGLAB version < X
try
    vers1 = eeg_getversion;
    firstdot = strfind(vers1,'.');
    vers2 = str2double(vers1(1:firstdot+1)); % take .X version
    if vers2 < 14.1 % pop_resample bug fixed in EEGLAB version 14.1
        fprintf('\n\n%s(): WARNING! This function only provides correct results if you\ndid *NOT* change the SAMPLING RATE of your data. Before version 14.1 of EEGLAB,\n\tfunction pop_resample did not update EEG.event.duration after resampling leading to\nwrong saccade and fixation durations.',mfilename);
    end
catch
end

% return history command string
allArgs = vararg2str({rate_event});
com = sprintf('%s(EEG,%s);',mfilename,allArgs);
return