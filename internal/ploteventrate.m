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
%
% Outputs:
%  timebins    - time points for plotting the rate
%  ratepersec  - frequency of events for plotting the rate
%  all_lats    - [double] latencies of all events of type "rate_event"
%
% Plots figure showing average rate (events per second) of the "rate_event"
% relative to time-locking event of the epoch
%
% This functions requires that saccade or fixation events were already
% imported or detected using pop_detecteyemovements()
%
% Note: this function is for epoched data only, because it plots the rate of
% event relative to the time-locking event(s) of the epochs. A rate plot
% for continuous data would not make too much sense.
%
% See also: pop_ploteventrate, pop_detecteyemovements, pop_ploteyemovments
%
% Example: A call to this function might look like this:
% >> [times4plot, rate4plot, all_lats] = ploteventrate(EEG,'saccade') % to make your own plots
%
% Note: if you plot the rate of fixations in epoched data, it will be
% very high in the first sample of the epoch, since the start of an
% epoch is also usually (by definition) the start of a fixation, if
% fixations were detected in epoched data rather than continuous data
% (and if the epoch does not start with a saccade, which is less likely).
%
% Author: od
% Copyright (C) 2009-2020 Olaf Dimigen, HU Berlin
% olaf.dimigen@hu-berlin.de

% Plans for future versions
% -- Allow to enter multiple event types that will be plotted on top of
%    each other
% -- Fix problem that "fixation" events are numerous at the first sample
%    of an epoch, which is completely logical but distracting

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

function [timebins, ratepersec, all_lats] = ploteventrate(EEG,rate_event)

if nargin < 2
    help(mfilename);
    return;
end

if size(EEG.data,3) == 1
    
    fprintf('\n%s: Your dataset is continuous.\n',mfilename)
    warning(sprintf('pop_ploteventrate() plots an event rate relative to the epoch time-locking event.\nIt therefore only works on epoched datasets.'))
    return;
    
end


%% plotting constants (change as you like)
AVG_EVENTS_PER_BIN = 30;
BINS_FOR_SMOOTHING = 5;  % smooth with moving avg. across X bins

try
    %% check whether any "rate_event" is present in EEG.event
    if ~isempty(rate_event)
        ix_re = find(cellfun(@(x) strcmp(x,rate_event),{EEG.event.type}), 1);
        if isempty(ix_re)
            warning('%s(): No events of type %s were found in EEG.event. Cannot plot their rate.', mfilename, rate_event)
            return
        end
    else
        error('%s rate_event not defined!',mfilename)
    end
    
    %% compute latencies of the event-of-choice
    [dummy, lats] = eeg_getepochevent(EEG,rate_event,[],'latency'); % output in ms
    all_lats = cell2mat(lats);
    
    %% estimate a reasonable number of "bins" for the rate histogram
    binsteps = round(1000 ./ (length(all_lats)./AVG_EVENTS_PER_BIN));
    edges    = EEG.times(1):binsteps:EEG.times(end);
    
    %% get histogram of events
    n_abs    = histc(all_lats,edges); % instead of histcounts() [backw. compatibility]
    % n_abs    = histcounts(all_lats,edges);
    
    %% compute rate per second (normalize histogram)
    epochduration    = EEG.pnts./EEG.srate;
    n_per_sec        = (n_abs / EEG.trials) / (epochduration/length(edges));
    % n_per_sec_smooth = smooth(n_per_sec,BINS_FOR_SMOOTHING);  % moving average
    n_per_sec_smooth = movingAverage(n_per_sec, BINS_FOR_SMOOTHING); % avoids commerical toolbox function smooth()
    
    %% plot rate figure
    figure; hold on;
    figtitle = sprintf('Rate of %s events',rate_event);
    title(figtitle);
    bar(edges,n_per_sec,'k');
    p2 = plot(edges,n_per_sec_smooth,'r','linewidth',2);
    l = legend(p2,sprintf('Rate smoothed over %i bins',BINS_FOR_SMOOTHING));
    set(l,'box','off');
    box on
    xlim([EEG.times(1) EEG.times(end)])
    xlabel('Time after time-locking event [ms]')
    ylabel(sprintf('%s events per second',rate_event))
    
    timebins   = edges;
    ratepersec = n_per_sec;
    
catch err
    
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
    
end
end % function ploteventrate

% function movingAverage
% helper function for smoothing: convolve with kernel to do moving average
function y = movingAverage(xx, w)
    k = ones(1, w)/w; % create kernel for convolution
    y = conv(xx,k,'same');
end












