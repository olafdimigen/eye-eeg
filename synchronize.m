% synchronize() - subfunction called by pop_importeyetracker()
%          for help type >> help pop_importeyetracker
%
% Usage:
%   >> eyetracker = synchronize(ET, startEndEvent, eegEvents, eegrate, ...
%                   n_eegsmp, doRegression, filterEyetrack, plotFig)
%
% Authors: ur & od
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

function [ET, eegEvents, syncQuality] = synchronize(ET, startEndEvent, eegEvents, eegrate, n_eegsmp, doRegression, filterEyetrack, plotFig, searchRadius)

if nargin < 9
    help(mfilename)
    return;
end

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

%% remember no. of EEG samples recorded between start-event and end-event
n_eegsmp_range       = eegEvents(end,2)-eegEvents(1,2)+1;
assert(n_eegsmp_range > 0,'\n%s(): Synchronization not possible. There is a problem with your chosen synchronization events in the EEG.\nThe last instance of your end-event %i does not happen after the first instance of your start-event %i.\n',mfilename,endevent,startevent);


%% estimate ET sampling rate (for user feedback/possible filtering)
% The ET sampling rate is estimated based on the assumption that
% the "EEG.srate" is correct, i.e. the clock of the EEG is treated as the
% master clock.
%
% How many samples were recorded between start- and end-event?
n_eyesmp_range = sum(ET.data(:,1) >= starteventTime & ET.data(:,1) <= endeventTime);
ix_starteventSample = find(ET.data(:,1) == starteventTime);
ix_endeventSample   = find(ET.data(:,1) == endeventTime);
% Bugfix, OD, 2017-02-15:
% special case for estimating the ET sampling rate in cases, where the
% parallel port inputs wer inserted as separate lines into the
% ET raw data (e.g. "INPUT" lines for Eyelink trackers). These lines have
% a timestamp that is slightly *different* from that of any data sample.
% Therefore, we need to search for the data sample that is closest to the
% input
if isempty(ix_starteventSample)
    %[ix_starteventSample, dummy] = searchclosest(ET.data(:,1),starteventTime);
    [ix_starteventSample] = nearestpoint(starteventTime,ET.data(:,1));
end
if isempty(ix_endeventSample)
    %[ix_endeventSample, dummy] = searchclosest(ET.data(:,1),endeventTime);
    [ix_endeventSample] = nearestpoint(endeventTime,ET.data(:,1));
end
n_etsmp_range = length(ix_starteventSample:ix_endeventSample);
eyerate = (n_etsmp_range / n_eegsmp_range) * eegrate; % estimated ET rate
% To do for future versions:
% Rather than just counting the ET samples to estimate the ET sampling rate,
% we need to considerer that the ET recording may have occasionally been
% *paused* by the user (causing large forward jumps in the ET timestamp).
% Therefore, we could use the median (robust against occasional recording
% pauses) interval between sucessive ET samples to estimate ET rate...
%
% % get median inter-sample-intervals for ET...
% et_df = median(diff(ET.data(ix_starteventSample:ix_endeventSample,1)));

%% feedback: sampling rates
fprintf('\n\n-- %i EEG samples in synchronization range',n_eegsmp_range);
fprintf('\n-- %i ET samples in synchronization range',n_eyesmp_range);
fprintf('\n-- %.2f Hz EEG sampling rate',eegrate);
fprintf('\n-- %.2f Hz estimated ET sampling rate [*]',eyerate);
fprintf('\n   ([*] taking EEG as master clock & assuming the ET recording was never paused)');
if eyerate > eegrate
    fprintf('\n-- Eye track will be downsampled to match the EEG sampling rate\n');
    fprintf('\n-- Note: the toolbox does not yet implement a low-pass filter to prevent aliasing.\n');
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

%% for each existing ET event, search for closest timestamp in new time
new_ix = zeros(length(ET.event),1);
% for k = ix_startevent:ix_endevent
%     new_ix(k) = searchclosest(ET.newtime,ET.event(k,1));
% end
new_ix(ix_startevent:ix_endevent) = nearestpoint(ET.event(ix_startevent:ix_endevent,1),ET.newtime);

ET.event(:,3) = new_ix; % assign updated sample index to ET events

%% identify "shared" events in ET & EEG
% find EEG events that have an ET event of same type no more than
% searchRadius samples away
eegEventHasPartner = findmatchingevents(eegEvents,ET.event,searchRadius);

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
    yhat = b(1)+b(2).*x; % for plotting of regression line only
    % sum of squares of data residuals
    st = sum((y-mean(y)).^2);
    % sum of squares of estimate residuals
    sr = sum((y-yhat).^2);
    % coefficient of determination (R^2)
    r2 = (st-sr) / st;
    % % root mean square error (RMSE)
    % rmse = sqrt(sr);
    
    % correct latency of start-/end-event based on this linear model fit
    ET.event(ix_startevent,1) = round( b(1)+b(2).*eegEvents(1,2)   );
    ET.event(ix_endevent,1)   = round( b(1)+b(2).*eegEvents(end,2) );
    starteventTime = ET.event(ix_startevent,1);
    endeventTime   = ET.event(ix_endevent,1);
    
    % now repeat linspace() with "refined" latencies of start-/end-event
    ET.newtime = linspace(starteventTime, endeventTime, n_eegsmp_range)';
end


%% linear interpolation
if ~isempty(ET.data)
    
    %% check for backward jumps in timestamps
    % Note: with some older eye trackers, ET timestamps can sometimes jump
    % backwards during the recording. This cannot be repaired since it remains
    % unknown how much time really passed during such a backward jump.
    % This problem was observed, in particular, with SMI trackers in combination
    % with older versions of SMI's "IView X" recording software
    assert(all(diff(ET.data(:,1)) > 0),'\n%s(): Sample time in the eye tracking data is not monotonically increasing.\nCheck for backward jumps in the timestamp column of your raw data and report this to your ET manufacturer.',mfilename);
    
    %% HANDLE RECORDING PAUSES IN EYE-TRACKING DATA 
    % If the ET recording was paused (generally not recommended) there will 
    % be jumps in the timestamp. During the linear interpolation (see below) 
    % the ET signal during these gaps will be also be linearly interpolated 
    % (i.e. there will be a linear trend between the last sample before and 
    % first sample after the pause. This is not helpful, because these 
    % intervals will not be easily detected as missing data (e.g. by saccade
    % detection algorithm or by rej_eyecontin.m), because they contain 
    % plausible non-zero values. 
    % This next step therefore detects jumps in the ET time stamp and sets 
    % the two samples adjacent to the recording pause to zero, so the pause 
    % will just contain zeros, event after interpolation and will be easyly
    % OD, 20180509
    %
    % Note: the threshold vor a "pause" was arbitraily set to the duration
    % of three median inter-sample intervals. If your eye-tracker produces
    % larger irregularites in the time stamps even without a pause, this 
    % values would need to be set higher
    SET_REC_PAUSES_TO_ZERO = true;
    PAUSE_THRESHOLD        = 3; % jump of 3 median inter-sample-intervals
    
    if SET_REC_PAUSES_TO_ZERO
                
        d               = diff(ET.data(:,1));
        ix_StartOfPause = find(d > PAUSE_THRESHOLD * median(d)); % find beginning of pauses
        % set first sample before pause and last sample after pause to zero
        zerosamplesA = ix_StartOfPause;   % start of pause
        zerosamplesB = ix_StartOfPause+1; % end of pause
        ET.data(unique([zerosamplesA zerosamplesB]), 2:end) = 0;
        
        % give feedback
        if ~isempty(ix_StartOfPause)
            fprintf('\n\n%s(): I detected %i PAUSES in the eye-tracking recording!',mfilename,length(ix_StartOfPause))
            fprintf('\nPauses were defined as jumps in the ET time stamp > %i samples',PAUSE_THRESHOLD)
            fprintf('\nThe last data sample before and first sample after each pause')
            fprintf('\nwill be set to zero to avoid interpolation artifacts during the pause!\n')
        end
    end
    
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
        % end
    end
    
    
    %% Linear piecewise interpolation (using interp1): generate new ET data
    % note: minimal extrapolation may be necessary if doRegression == true
    et_new = interp1(ET.data(:,1),ET.data(:,2:end),ET.newtime,'linear','extrap');
    
    %% ET is upsampled --> filter to prevent images in spectrum
    if filterEyetrack
        % if eegrate > eyerate
        %   fprintf('\nET data is being low-pass filtered (IIR) to prevent image artifacts in spectrum...\n');
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
% for k = ix_startevent:ix_endevent
%     new_ix(k) = searchclosest(ET.newtime,ET.event(k,1));
% end
new_ix(ix_startevent:ix_endevent) = nearestpoint(ET.event(ix_startevent:ix_endevent,1),ET.newtime);

ET.event(:,3) = new_ix;

%% update sample index of imported eye movement events in ET
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
        
        % find corresponding sample in new ('rescaled') ET time
        new_ix   = zeros(n_events,3);
        %new_time = zeros(n_events,3);
        
        % with searchclosest, somewhat obfuscated code
        % for e = 1:2*n_events
        %   [new_ix(e) new_time(e)] = searchclosest(ET.newtime,ET.eyeevent.(eventType).data(e));
        % end
        
        e = 1:n_events;
        % find event start sample in rescaled time
        new_ix(:,1)   = nearestpoint(ET.eyeevent.(eventType).data(e,1), ET.newtime);
        %new_time(:,1) = ET.newtime(new_ix(:,1));
        
        % find event end sample in rescaled time
        new_ix(:,2)   = nearestpoint(ET.eyeevent.(eventType).data(e,2), ET.newtime);
        %new_time(:,2) = ET.newtime(new_ix(:,2)); 
        
        % update event "duration"
        new_ix(:,3) = new_ix(:,2)-new_ix(:,1)+1;
        
        % translate to event indices in continuous EEG (add offset)
        new_ix(:,1:2) = new_ix(:,1:2)-1 + eegEvents(1,2);
        
        % store event
        ET.eyeevent.(eventType).data(:,1:3) = new_ix;
    end
    
    clear new_* n_events inRange
end


%% update timestamp of 'other' ET messages (new in Jan-2017, OD)
% = all messages starting with 'MSG', which are *not* keyword messages
% (used for synchronization) and *not* eyeevents (e.g. saccades).
if isfield(ET,'othermessages')
    
    % delete messages with timestamps not within synchronizeable time range
    NotInRange = [ET.othermessages(:).timestamp] < starteventTime | [ET.othermessages(:).timestamp] > endeventTime;
    ET.othermessages(NotInRange) = [];
    %     % go tru 'other' messages
    %     for m = 1:length(ET.othermessages)
    %         % find corresponding sample in new ('rescaled') ET time
    %         [new_ix, new_time] = searchclosest(ET.newtime,ET.othermessages(m).timestamp);
    %         new_ix             = new_ix-1 + eegEvents(1,2);
    %         % add the corresponding EEG time (= sample of EEG corresponding to
    %         % time that the message was send) to ET.othermessages
    %         ET.othermessages(m).EEGsample = new_ix;
    %     end    
    m = 1:length(ET.othermessages);
    new_ix = nearestpoint([ET.othermessages(m).timestamp], ET.newtime);
    new_ix = new_ix-1 + eegEvents(1,2);
    for j = m
        ET.othermessages(j).EEGsample = new_ix(j);
    end
end

%% find matching events again to update estimate of synchronization error
eegEventHasPartner = findmatchingevents(eegEvents,ET.event,searchRadius);

%% feedback: estimate sync error based on all "shared" events
n_total      = size(eegEventHasPartner,1);
n_nopartner  = sum(isnan(eegEventHasPartner(:,2)));
[count, bin] = hist(eegEventHasPartner(:,end),-searchRadius:searchRadius);
syncQuality  = [bin;count]';

%% feedback on synch. quality
fprintf('\nShared events of same type less than +/- %i samples apart: %i',searchRadius,n_total-n_nopartner);
fprintf('\n\nSynchronization quality:')
if n_nopartner == 0
    fprintf('\nFor all events in the EEG, an ET event of the same type was found within plusminus %i samples.',searchRadius);
    fprintf('\nSynch. error is the latency difference between matching ET/EEG events after synchronization.');
    fprintf('\nIt is distributed as follows:\n');
else
    fprintf('\nWarning: For %i of %i EEG events (%.1f%%), no ET event of same name was found within plusminus %i samples',n_nopartner,n_total,100*n_nopartner/n_total,searchRadius);
    fprintf('\nThis can occur, for example if some events were not transmitted to the ET (or EEG)');
    fprintf('\nIf only a few events do not have a \"partner\" event, and sync quality for the remaining events is good, you are probably safe.');
    fprintf('\nSynch. error is the latency difference between matching ET/EEG events after synchronization.');
    fprintf('\nFor the remaining %i events, it is distributed as follows:\n',n_total-n_nopartner);
end
fprintf('\n%s\t%s','Error [smp]','Events');
fprintf('\n-------------------------\n');
disp(syncQuality);
fprintf('-------------------------\n');
avg_abs_error = sum(abs(syncQuality(:,1)).*syncQuality(:,2)) ./ sum(syncQuality(:,2));
fprintf('Mean abs. sync. error (estimated from \"shared\" events): %.3f ms',avg_abs_error.*1000/eegrate); % note this is *not* RMSE
fprintf('\n-------------------------\n\n');

%% figure with feedback on synch quality
% three subplots:
% - event latencies visualized after synchronization (all events)
% - scatterplot of event latencies in original time ('shared' events only)
% - histogram of sync. error ('shared' events only)
if plotFig
    fprintf('\nMaking plot visualizing synchronization results...')
    
    syncfig = figure('Name','synchronize(): Synchronization results');

    % hotfix for v.085 due to changed legend behavior in R2017a 
    % avoid that legend auto-updates with new generic entries (e.g. "data 1")
    % this crashs scripts for large legend (Thanks, TheMathworks!)
    % 20180624 - OD
    try
        set(syncfig,'defaultLegendAutoUpdate','off')
    catch legend_autoupdate_err
    end
        
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
    % generate background
    h_bgGr = hggroup;
    % plot areas outside sync. range grey
    fill([0 eeg_abssmp(1) eeg_abssmp(1) 0],[0 0 1 1],[.5 .5 .5], 'Parent',h_bgGr);
    fill([eeg_abssmp(end) n_eegsmp n_eegsmp eeg_abssmp(end)],[0 0 1 1],[.5 .5 .5], 'Parent',h_bgGr);
    % horizontal separator line
    plot([0 n_eegsmp],[0.5 0.5],'linestyle','-','color',[0 0 0],'linewidth',1.0, 'Parent',h_bgGr)
    % handle to feed legend()
    
    hleg = zeros(size(types_unique));
    
    % loop through unique event types
    for loop = 1:length(types_unique)
        fprintf('\nPlotting events of type %i...',types_unique(loop))
        % 1. plot latencies of EEG triggers between start-event & end-event
        ix1 = find(types_unique(loop) == eeg_type);
        if ~isempty(ix1)
            plot(eeg_abssmp(ix1),0.499,'linestyle','none','marker','.','color',rgb_type(loop,1:3),'handleVisibility','off'); % "needle hats"
            hline = plot(eeg_abssmp([ix1 ix1])', [0 0.499],'linestyle','-','color',rgb_type(loop,1:3),'linewidth',1,'handleVisibility','off');
            hleg(loop) = hline(1);
        end
        % 2. plot latencies of ET triggers between start-event & end-event
        ix2 = find(types_unique(loop) == et_type);
        if ~isempty(ix2)
            plot(et_abssmp(ix2),0.501,'linestyle','none','marker','.','color',rgb_type(loop,1:3),'handleVisibility','off');
            hline = plot(et_abssmp([ix2 ix2])', [0.501 1],'linestyle','-','color',rgb_type(loop,1:3),'linewidth',1,'handleVisibility','off');
            if isempty(ix1) % unique in ET
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
    title('Comparison: Event latencies','fontweight','bold')
    if doRegression
        plot(x,yhat,'r-');
        plot(x,y,'k.');
        plot(x(1),yhat(1),'bo');     % updated estimate: start-event latency
        plot(x(end),yhat(end),'bo'); % updated estimate: end-event latency
        l = legend('regression line (ET on EEG)','shared events','start-/endEvent','Location','NorthWest');
        % coeff. of determination (R2)
        xl = xlim; yl = ylim;
        text(xl(2)-0.4*(xl(2)-xl(1)),yl(2)-0.6*(yl(2)-yl(1)),sprintf('R^2 = %.3f',r2),'fontsize',10);
    else
        plot(x,y,'k.');
        l = legend('shared events','Location','NorthWest');
    end
    set(l,'box','off','fontsize',8);
    xlabel('EEG latency (in samples)')
    ylabel('Original ET latency (timestamp)')
    
    %% histogram of synchronization error (based on "shared" events only)
    % error = sample dist. between EEG and ET events after interpolation
    subplot(2,2,4); hold on; box on;
    title('Quality of synchronization','fontweight','bold')
    bar(bin,count,'k')
    set(gca,'xTick',-searchRadius:1:searchRadius);
    xlim([-searchRadius-0.5 searchRadius+0.5])
    xlabel('Time diff. between shared events (in samples)')
    ylabel('Number of events')
end