% synchronize() - subfunction called by pop_importeyetracker()
%          for help type >> help pop_importeyetracker
%
% Usage:
%   >> eyetracker = synchronize(ET, startEndEvent, eegEvents, eegrate, ...
%                   n_eegsmp, doRegression, filterEyetrack, plotFig)
%
% Authors: ur & od
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

function [ET, eegEvents, syncQuality] = synchronize(ET, startEndEvent, eegEvents, eegrate, n_eegsmp, doRegression, filterEyetrack, plotFig)

if nargin < 8
    help(mfilename)
    return;
end

%% define search radius to find matching (shared) events in ET & EEG
% Search for RADIUS samples around each EEG event to find an ET event of 
% the same type (e.g. '123'). Hard-coded. It might be sensible to increase 
% the search RADIUS in case of very high sampling rates (e.g. 2000 Hz).
RADIUS = 5;

%% get start-event & end-event for synchronization
startevent = startEndEvent(1);
endevent   = startEndEvent(2);
clear startEndEvent

%% find those events in ET.event: take first "start-event", last "end-event"
ix_startevent = find(ET.event(:,2) == startevent,1);
ix_endevent   = find(ET.event(:,2) == endevent,1,'last');
assert(ix_endevent > ix_startevent,'\n%s(): Synchronization not possible. There is a problem with your chosen synchronization events in the ET data.\nThe last instance of your end-event %i does not occur after the first instance of your start-event %i.\n',mfilename,endevent,startevent);
% ...and get their time in samples
starteventTime = ET.event(ix_startevent,1);
endeventTime   = ET.event(ix_endevent,  1);

%% for EEG events, keep only those between start-event and end-event
eegEvents      = eegEvents(find(eegEvents(:,1) == startevent,1):find(eegEvents(:,1) == endevent,1,'last'),:);
eegEvents(:,3) = eegEvents(:,2)-eegEvents(1,2)+1; % subtract offset

% remember no. of EEG samples recorded between start-event and end-event
n_eegsmp_range       = eegEvents(end,2)-eegEvents(1,2)+1;
assert(n_eegsmp_range > 0,'\n%s(): Synchronization not possible. There is a problem with your chosen synchronization events in the EEG.\nThe last instance of your end-event %i does not happen after the first instance of your start-event %i.\n',mfilename,endevent,startevent);


%% estimate ET sampling rate
% How many samples were recorded between start- and end-event?
n_eyesmp_range = sum(ET.data(:,1) >= starteventTime & ET.data(:,1) <= endeventTime);

% Count number of ET samples recorded and determine median 
% inter-sample-interval. Rather than just counting samples, this takes into 
% consideration that the ET recording may have occasionally been paused 
% causing large forward jumps in the timestamp. Note that the ETr sampling 
% rate is estimated based on the assumption that "EEG.srate" is correct
ix_starteventSample = find(ET.data(:,1) == starteventTime);
ix_endeventSample   = find(ET.data(:,1) == endeventTime);
% get inter-sample-intervals
et_df      = diff(ET.data(ix_starteventSample:ix_endeventSample,1));
% divide by median inter-sample interval
stepFactor = et_df/median(et_df);
eyerate    = (sum(stepFactor) / n_eegsmp_range) * eegrate;

%% feedback: sampling rates
fprintf('\n\n-- %i recorded EEG samples in synchronization range',n_eegsmp_range);       
fprintf('\n-- %i recorded ET samples in synchronization range',n_eyesmp_range); 
fprintf('\n-- %.2f Hz EEG sampling rate',eegrate);
fprintf('\n-- %.2f Hz estimated ET sampling rate',eyerate);
if eyerate > eegrate
    fprintf('\n-- Eye track will be downsampled to match the EEG sampling rate\n');
    fprintf('\n-- Note: the plugin does not yet implement a low-pass filter to prevent aliasing.\n');
else
    fprintf('\n-- Eye track will be upsampled to match the EEG sampling rate\n');
end

%% produce new eye tracker time
% create [n_eegsmp_range] linear spaced new sample times between
% starteventTime and endeventTime. Creating this new, regularly spaced
% timevector for interpolation of ET data is based on the assumptions that
% the EEG a has truly regular sampling interval and that startevent and
% endevent were transmitted to each system without signicant jitter/delay
ET.newtime = linspace(starteventTime, endeventTime, n_eegsmp_range)';

% for each existing ET event, search for closest timestamp in new time
new_ix = zeros(length(ET.event),1);
for k = ix_startevent:ix_endevent
    new_ix(k) = searchclosest(ET.newtime,ET.event(k,1));
end
ET.event(:,3) = new_ix; % assign updated sample index to ET events

% identify "shared" events in ET & EEG, i.e. find EEG events that have an 
% ET event of same type no more than RADIUS samples away
eegEventHasPartner = findmatchingevents(eegEvents,ET.event,RADIUS);

ix_keepEEG = find(~isnan(eegEventHasPartner(:,1))); % EEG events with partner
ix_keepET  = eegEventHasPartner(ix_keepEEG,1);      % ET events that are partners of EEG events
x = eegEvents(ix_keepEEG,2); % sample of EEG events with a partner
y = ET.event(ix_keepET,1);   % timestamp of ET events with a partner

%% optional: linear regression
if doRegression
    
    % do linear regression based on all shared events to optimize
    % the recorded latencies of start-event and end-event. Helps to reduce 
    % error if start-event and/or end-event were transmitted with jitter.
    % Won't help in case of constant transmission delays for all events.
    
    % linear regression based on shared events
    b  = [ones(length(x),1) x] \ y;   
    % for plotting only
    yhat = b(1)+b(2).*x;
    % sum of squares of data residuals
    st = sum((y-mean(y)).^2);
    % sum of squares of estimate residuals
    sr = sum((y-yhat).^2);
    % coefficient of determination (R^2)
    r2 = (st-sr) / st;
    
    % correct latency of start-/end-event using regression function
    ET.event(ix_startevent,1) = round( b(1)+b(2).*eegEvents(1,2)   );
    ET.event(ix_endevent,1)   = round( b(1)+b(2).*eegEvents(end,2) );
    starteventTime = ET.event(ix_startevent,1);
    endeventTime   = ET.event(ix_endevent,1);
    
    % repeat linspace() with "refined" latencies of start-/end-event
    ET.newtime = linspace(starteventTime, endeventTime, n_eegsmp_range)';
end


%% linear interpolation
if ~isempty(ET.data)
    
    %% check for backward jumps in time
    % Note: occasionally ET timestamps can jump backwards during
    % continuous recording. This cannot be repaired since it remains
    % unknown how much time did really pass during a backward jump. 
    % Problem observed with SMI trackers & older versions of the IView X 
    % recording software
    assert(all(diff(ET.data(:,1)) > 0),'\n%s(): Sample time in the eye tracking data is not monotonically increasing.\nCheck for backward jumps in the timestamp column of your raw data and report this to your ET manufacturer.',mfilename);

    %% check whether filters are installed
    if filterEyetrack
        fprintf('\nWarning: Filtering to prevent aliasing is not yet implemented.')
        fprintf('\nIt is planned for future versions.')
    end    
    
    %% ET is downsampled --> filter to prevent aliasing
    if filterEyetrack
        % hicutoff = min([eegrate eyerate])/2;
        % if eegrate <= eyerate
        %   fprintf('\nET data is being low-pass filtered (IIR) to prevent aliasing...\n');
        %   ...
        %   ...
        % end
    end
    
    %% linear piecewise interpolation by interp1: generate new ET data
    % note: minimal extrapolation may be necessary if doRegression == true
    et_new = interp1(ET.data(:,1),ET.data(:,2:end),ET.newtime,'linear','extrap');

    %% ET is upsampled --> filter to prevent images in spectrum
    if filterEyetrack
        % if eegrate > eyerate
        %   fprintf('\nET data is being low-pass filtered (IIR) to prevent image artifacts in spectrum...\n');
        %   ...
        %   ...
        % end
    end
    
    %% store new ET data columns
    ET.syncdata = [ET.newtime et_new];
else
    ET.syncdata = [];
end


%% update the sample index of events in the ET
% for existing ET events, search for closest timestamp in the new time
new_ix = zeros(length(ET.event),1);
for k = ix_startevent:ix_endevent
    new_ix(k) = searchclosest(ET.newtime,ET.event(k,1));
end
ET.event(:,3) = new_ix;


%% update the sample index of imported eye movement events in ET
if isfield(ET,'eyeevent')
    
    eventTypes = fieldnames(ET.eyeevent)';
    
    % go tru eye movement event types in Eyelink data (e.g., 'L_fixation', 'R_saccade',...)
    for t = 1:length(eventTypes)
        eventType = eventTypes{t};
        
        % delete events not within synchronized time range
        inRange = (ET.eyeevent.(eventType).data(:,1) > starteventTime) & (ET.eyeevent.(eventType).data(:,2) < endeventTime) ;
        ET.eyeevent.(eventType).data = ET.eyeevent.(eventType).data(inRange,:);
        ET.eyeevent.(eventType).eye  = ET.eyeevent.(eventType).eye(inRange,:);
        n_events = size(ET.eyeevent.(eventType).data,1);
        
        % find corresponding sample in new ET time
        new_ix   = zeros(n_events,3);
        new_time = zeros(n_events,3);
        for e = 1:2*n_events
            [new_ix(e) new_time(e)] = searchclosest(ET.newtime,ET.eyeevent.(eventType).data(e));
        end
        
        % update event "duration"
        new_ix(:,3) = new_ix(:,2)-new_ix(:,1)+1;
        
        % translate to event indices in continuous EEG (add offset)
        new_ix(:,1:2) = new_ix(:,1:2)-1 + eegEvents(1,2);
        
        % store event
        ET.eyeevent.(eventType).data(:,1:3) = new_ix;
    end
end

%% find matching events again to update estimate of synchronization error
eegEventHasPartner = findmatchingevents(eegEvents,ET.event,RADIUS);


%% feedback: estimate sync error based on all "shared" events
n_total     = size(eegEventHasPartner,1);
n_nopartner = sum(isnan(eegEventHasPartner(:,2)));
[count bin] = hist(eegEventHasPartner(:,end),-RADIUS:RADIUS);
syncQuality = [bin;count]';


%% feedback on synchronization quality: MATLAB console
fprintf('\nSynchronization quality:')
if n_nopartner == 0
    fprintf('\nFor each EEG event, an ET event of the same type was found within plusminus %i samples',RADIUS);
    fprintf('\nSynch. error is the latency difference between matching ET/EEG events after synchronization.');
    fprintf('\nIt is distributed as follows:\n');
else
    fprintf('\nWarning: For %i of %i EEG events (%.1f%%), no ET event of same name was found within plusminus %i samples',n_nopartner,n_total,100*n_nopartner/n_total,RADIUS);
    fprintf('\nPossibly, some events were not transmitted to the eye tracker?');
    fprintf('\nSynch. error is the latency difference between matching ET/EEG events after synchronization.');
    fprintf('\nFor the remaining %i events, it is distributed as follows:\n',n_total-n_nopartner);
end
fprintf('\n%s\t%s','Error [smp]','Events');
fprintf('\n-------------------------\n');
disp(syncQuality);
fprintf('-------------------------');
fprintf('\nShared events less than +/- %i samples apart: %i\n',RADIUS,n_total-n_nopartner);


%% figure with feedback on synchronizatino quality
% three subplots:
% - event latencies visualized after synchronization (all events)
% - scatterplot of event latencies in original time ('shared' events only)
% - histogram of sync. error ('shared' events only)
if plotFig
    figure('Name','synchronize(): Synchronization results');
    
    %% latency comparisons after synchronization (all events)
    subplot(2,2,1:2); hold on; box on;
    title('Overview: Event latencies in synchronized dataset','fontweight','bold')
    
    % get random set of colors (hsv randomized)
    colors_rgb   = colormap('hsv');
    colors_rgb   = [colors_rgb rand(size(colors_rgb,1),1)];
    colors_rgb   = sortrows(colors_rgb,4);
    
    % trigger types
    eeg_type     = eegEvents(:,1);
    et_type      = ET.event(ix_startevent:ix_endevent,2);
    types_unique = unique([eeg_type; et_type]); % unique trigger types
    
    % get trigger latencies (after interpolation) relative to start of
    % the EEG recording (rather than relative to the start-event).
    % Reason: also display time intervals outside of synchronization range
    eeg_abssmp = eegEvents(:,2);
    offset     = eeg_abssmp(1)-1; % samples until start-event
    et_abssmp  = ET.event(ix_startevent:ix_endevent,3)+offset;
    
    % randomized color table over unique types
    rgb_type = colors_rgb( mod( (1:size(types_unique,1))-1, size(colors_rgb,1)) +1, :); % pick color
    % generate brackground
    h_bgGr = hggroup;
    % plot areas outside sync. range grey
    fill([0 eeg_abssmp(1) eeg_abssmp(1) 0],[0 0 1 1],[.5 .5 .5], 'Parent',h_bgGr);
    fill([eeg_abssmp(end) n_eegsmp n_eegsmp eeg_abssmp(end)],[0 0 1 1],[.5 .5 .5], 'Parent',h_bgGr);
    % horizontal separator line
    plot([0 n_eegsmp],[0.5 0.5],'linestyle','-','color',[0 0 0],'linewidth',1.0, 'Parent',h_bgGr)
    % handle to feed legend()
    hleg = zeros(size(types_unique));
    % loop through unique event types
    for loop = 1 : length(types_unique);
        % 1. plot latencies of EEG triggers between start-event & end-event
        ix1 = find(types_unique(loop) == eeg_type);
        if ix1
            plot(eeg_abssmp(ix1),0.499,'linestyle','none','marker','.','color',rgb_type(loop,1:3),'handleVisibility','off');
            hline =  plot(eeg_abssmp([ix1 ix1])', [0 0.499],'linestyle','-','color',rgb_type(loop,1:3),'linewidth',1,'handleVisibility','off');
            hleg(loop) = hline(1);
        end
        % 2. plot latencies of ET triggers between start-event & end-event
        ix2 = find(types_unique(loop) == et_type);
        if ix2
            plot(et_abssmp(ix2),0.501,'linestyle','none','marker','.','color',rgb_type(loop,1:3),'handleVisibility','off');
            hline =    plot(et_abssmp([ix2 ix2])', [0.501 1],'linestyle','-','color',rgb_type(loop,1:3),'linewidth',1,'handleVisibility','off');
            if isempty(ix1) % unique in et
                hleg(loop) = hline(1);
            end
        end
    end
    ylim([0 1])
    xlim([0 n_eegsmp])
    xlabel('Time after start of EEG recording [samples]')
    set(gca,'yTick',[0.25 0.75]);
    set(gca,'yTickLabel',{'EEG events','ET events'});
    % show legend
    try % hotfix for v.0.337
        get(hleg,{'DisplayName'});
        set(hleg,'handleVisibility','on');
        set(hleg,{'DisplayName'}, cellstr(num2str(types_unique)) );
        legend('show')
    catch
        fprintf('\nSubplot legend not shown.')
    end
        
    %% scatterplot of event latencies ("shared" events only)
    subplot(2,2,3); hold on; box on;
    title('Event latencies of shared events','fontweight','bold')
    if doRegression
        plot(x,yhat,'r-');
        plot(x,y,'k.');
        plot(x(1),yhat(1),'ro');     % updated estimate: start-event latency
        plot(x(end),yhat(end),'ro'); % updated estimate: end-event latency
        l = legend('regression: ET on EEG','shared events','Location','NorthWest');
        % coeff. of determination (R2)
        xl = xlim; yl = ylim;
        text(xl(2)-0.4*(xl(2)-xl(1)),yl(2)-0.6*(yl(2)-yl(1)),sprintf('R^2: %.3f',r2),'fontsize',8);
    else
        plot(x,y,'k.');
        l = legend('shared events','Location','NorthWest');
    end
    set(l,'box','off','fontsize',8);
    xlabel('EEG latency (samples)')
    ylabel('Original ET latency (samples or time)')
    
    
    %% histogram of synchronization error (based on "shared" events only)
    % error = sample dist. between EEG and ET events after interpolation
    subplot(2,2,4); hold on; box on;
    title('Quality of synchronization','fontweight','bold')
    bar(bin,count,'k')
    set(gca,'xTick',-RADIUS:1:RADIUS);
    xlim([-RADIUS-0.5 RADIUS+0.5])
    xlabel('Offset between shared events (samples)')
    ylabel('Number of events')
end