% findmatchingevents() - find corresponding events in EEG and eye tracker.
%
% Go through events in the EEG, identify those events where an ET event 
% of the same type (e.g., type "123") is found within a search radius of
% "sampleradius" samples.
%
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

function [ hasPartner ] = findmatchingevents( events_eeg, events_et, sampleradius )

hasPartner = zeros(size(events_eeg,1),2);

% go tru EEG events (events_eeg)
% is there an event_et of same type within sampleradius samples?
for e = 1:size(events_eeg,1)
    try
        % find events of the same type (e.g., "123") in ET
        sameType = find(events_et(:,2) == events_eeg(e,1));
        
        % get temporal distances to all matching events
        distances = (events_et(sameType,3) - events_eeg(e,3));
        
        % mark ET events within a reasonable temporal distance
        closeEnough = abs(distances) <= sampleradius;
        
        % multiple matching events nearby?
        if sum(closeEnough)>1
            fprintf('\n%s(): Warning: For event type \"%i\", multiple (%i) events were found within radius of %i samples!',mfilename,events_eeg(e,1),sum(closeEnough),sampleradius)
            [val_min,ix_min] = min(abs(distances(closeEnough)));
            % if multiple matching events are found, take closest one
            % prefer earlier one in case of same distance, i.e.,
            % prefer -1 over +1
            hasPartner(e,:) = [sameType(ix_min),val_min];
        else % default: one matching event found
            hasPartner(e,:) = [sameType(closeEnough), distances(closeEnough)];
        end
    catch % no matching event found
        hasPartner(e,:) = [NaN,NaN];
    end
end

% Note: the following exotic scenario is not addressed: two EEG events reference
% the same nearby ET event, because the event was only transmitted to the ET once (and lost once)

end