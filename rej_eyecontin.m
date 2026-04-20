% rej_eyecontin() - reject intervals of continuous data during which
%                the eye tracker recorded extreme (out-of-range) values.
%                This function is useful for several purposes: to reject
%                data with eye blinks, to control fixation position, and to
%                reject intervals bad eye tracking data prior to saccade or
%                fixation detection.
%
%                There are two possibilities:
%                METHOD 1: reject bad intervals and insert boundary events
%                (=data breaks) at that position. This actually removes the
%                intervals with bad data from the dataset.
%                EEGLAB function pop_select() is used to do the actual
%                rejection. A new boundary event is added at the position
%                of each of the resulting data breaks. To perform a similar
%                rejection for epoched data, use pop_rej_eyeepoch.
%
%                METHOD 2: detect bad intervals in the eye tracking data,
%                and insert a "bad_ET" marker in the EEG.event structure for
%                every bad interval. These events will be considered by 
%                pop_detecteyemovements(), so that the bad ET data will not
%                distort the detection of saccades and fixations
%
%
%                Update November 2024: The function can now receive as
%                input continuous or epoched data. If run on epoched data
%                with rejectionmethod == 2, it will add events to the EEG
%                structure with an additional column 'epoch'. This allows 
%                to mask bad ET samples when using detecteyemovements.m on 
%                epoched data, as epoching does not update event duration, 
%                and thus 'bad_ET' derived from continous data creates 
%                wrong bad_ET masks on epoched data.
%                Note that the epoch added to EEG.event is the epoch an
%                event started at; if an event lasts several epochs, this
%                will not appear on the 'epoch' column, but will still
%                work correctly in detecteyemovements.m. The new optional
%                third output, seqs_bad, is a struct with information
%                regrading bad ET sequences per epoch. 
%                This implementation currently does not allow to use
%                rejectionmethod == 1 with epoched data. 
%                Rotem Krispil, rotem.krispil@mail.huji.ac.il. 
%
% Usage:
%   >> EEG = rej_eyecontin(EEG,chns,minvals,maxvals,windowsize,rejectionmethod)
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
%   rejectionnmethod  -  [1 or 2]
%                  1: reject bad data intervals (cuts intervals from data!)
%                     and inserts boundary events at the breaks
%                  2: add special marker "bad_ET" to EEG.event for each bad 
%                     data interval. These markers are considered by the
%                     eye movement detection function (pop_detecteyemovements.m)
%                     The dataset is not cut.
%
% Outputs:
%   EEG         - EEG structure with bad intervals removed or detected.
%                 For each removed interval of bad data, a new boundary
%                 event is inserted into EEG.event.
%
%   An example call of the function might look like this:
%   >> EEG = rej_eyecontin(EEG,[33 34],[1 1],[1024 768],50, 2)
%
%   In this example, the data of channels 33 and 34, containing the
%   horizontal (33) and vertical (34) position of the left eye, are tested
%   for pixel values smaller than 1 pixel or larger than 1024 pixels
%   (channel 33) and smaller than 1 pixel or larger than 768 pixels
%   (channel 34). Screen resolution during the experiment was 1024 x 768
%   pixel. Values outside this range likely reflect eye blinks
%   or intervals during which the eye was not properly tracked. Around each
%   interval of out-of-range data an addtional plusminus 50 samples of
%   data are marked as well. For each interval of bad data, a "bad_ET"
%   event is added to EEG.event. These bad intervals will not be considered 
%   by the saccade detection algorithm (pop_detecteyemovements.m)
%
% See also: pop_rej_eyecontin, pop_rej_eyeepoch, pop_select
%
% Author: od, rk (2024 update)
% Copyright (C) 2009-2021 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf@dimigen.de
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

function [EEG,seq_bad_overall, seqs_bad] = rej_eyecontin(EEG,chns,minvals,maxvals,windowsize,rejectionmethod)


%% input checks
if isempty(chns)
    error('\n%s(): Please provide indices of channels to test for out-of-range values.',mfilename)
end
if isempty(minvals) && isempty(maxvals)
    error('\n%s(): Please provide either minimum (minvals) or maximum (maxvals) values (or both) to test for out-of-range samples.',mfilename)
end
% Either minvals or maxvals can be left empty; in this case only the other
% criterion is used. If both minvals and maxvals are specified, both
% vectors need to have the same length (length(chns))
checkminima = true;
checkmaxima = true;
if length(minvals) ~= length(maxvals)
    if isempty(minvals), checkminima = false; end % only test for minvals
    if isempty(maxvals), checkmaxima = false; end % only test for maxvals
    if checkminima && checkmaxima
        error('%s(): If minimum and maximum values are specified, both vectors need to have the same length!',mfilename)
    end
elseif length(chns) ~= length(minvals)
    error('%s(): For each channel to test for out-or-range data, please specify one value in minvals and maxvals!',mfilename)
end

%% collect indices of out-of-range samples across all channels
fprintf('\n%s(): Rejecting intervals with out-of-range eye track',mfilename)
ix_bad = {};
nEpochs = size(EEG.data, 3);
nSamples = size(EEG.data, 2);
for e = 1:nEpochs
    ix_bad_e = [];
    for c = 1:length(chns)
        [chnid chntxt] = eeg_decodechan(EEG.chanlocs,chns(c)); % get channel names

        if checkminima && checkmaxima
            fprintf('\nChannel %i \"%s\": rejecting values < %i or > %i',chnid,chntxt{:},minvals(c),maxvals(c));
            ix = find(EEG.data(chns(c),:, e) < minvals(c) | EEG.data(chns(c),:, e) > maxvals(c));
        elseif checkminima
            fprintf('\nChannel %i \"%s\": rejecting values < %i',chnid,chntxt{:},minvals(c));
            ix = find(EEG.data(chns(c),:, e) > maxvals(c));
        else
            fprintf('\nChannel %i \"%s\": rejecting values > %i',chnid,chntxt{:},maxvals(c));
            ix = find(EEG.data(chns(c),:, e) < minvals(c) );
        end
        ix_bad_e = [ix_bad_e ix];
        
    end
    ix_bad{e} = unique(ix_bad_e);
end
%ix_bad = unique(ix_bad);
if windowsize > 0
    fprintf('\nRejecting an extra plusminus %i samples (%.0f ms) around each out-of-range interval.',windowsize,windowsize*1000/EEG.srate)
end

%% bad samples found?
if isempty(ix_bad)
    fprintf('\n\nNo out-of-range samples found.\n')
else
    seqs_bad = {};
    for e = 1:nEpochs
        % Check if epoch has bad data
        if ~isempty(ix_bad{e})

            % identify start and end indices of contiguous blocks of bad data
            seq_bad = findsequence2(ix_bad{e}');

            % add extra rejection zone around bad data
            % e.g., to remove sub-threshold bad samples at beginning/end of blinks
            seq_bad(:,1) = seq_bad(:,1) - windowsize;
            seq_bad(:,2) = seq_bad(:,2) + windowsize;

            % remove indices outside data boundaries
            seq_bad(seq_bad(:,1) < 1,1) = 1;
            seq_bad(seq_bad(:,2) > nSamples,2) = nSamples;
            
            seqs_bad{e} = seq_bad;
        end
        
    end
    % to handle overlapping blocks of bad data
    % translate again into vector of bad samples (0 = good, 1 = bad)
    badvectoroverall = zeros(nSamples * nEpochs,1);
    for e = 1:nEpochs
        if ~isempty(seqs_bad{e})
            seq_bad = seqs_bad{e};
            badvector = zeros(nSamples, 1);
            for n = 1:size(seq_bad,1)
                badvectoroverall((seq_bad(n,1) + nSamples*(e-1)) : (seq_bad(n,2) + nSamples*(e-1))) = 1;
                badvector(seq_bad(n,1):seq_bad(n,2)) = 1;
            end

            % get start/end of bad blocks again
            seqs_bad{e} = findsequence2(find(badvector));
        end
    end
    
    % get start/end of bad blocks again
    seq_bad_overall = findsequence2(find(badvectoroverall));
    
    if nEpochs > 1
        % Add epoch column
        epochs = floor(seq_bad_overall(:,1)/nSamples) + 1;
        seq_bad_overall(:,4) = epochs;
    else
        % If single epoch, extract from struct. seq_bad_overall is a double
        % already if single epoch - maintaining compatability with previous
        % version. 
        seqs_bad = seqs_bad{1};
    end
    
    fprintf('\n\nFound %i out-of-range samples in %i continuous intervals of bad data',sum(badvectoroverall),size(seq_bad_overall,1));
    
    % reject bad data
    switch rejectionmethod
        
        case 1 % remove using pop_select
            if nEpochs == 1
                fprintf('\nRemoving bad intervals from continuous data...\n\n')
                EEG = pop_select(EEG,'nopoint',seq_bad_overall(:,1:2));
            else
                error('%s(): Attempting to remove bad intervals from epoched data', mfilename)
            end
        case 2 % add "bad_ET" marker for each bad interval

            % Update march 2018 (OD)
            % Do not remove data, but introduce new bad_ET events in EEG.event
            fprintf('\nAdding %i bad_ET markers to the EEG.event structure...\nNot removing any data.\n\n',size(seq_bad,1))
            % Are there already bad_ET events (e.g., from runnning the 
            % function previously on continuous data)? If so, remove them.
            old_badET_events = contains({EEG.event.type}, {'bad_ET'});
            EEG = pop_editeventvals(EEG,'delete', find(old_badET_events));
            % add bad ET markers
            if nEpochs == 1  
                EEG = addevents(EEG,[seq_bad_overall(:,1) seq_bad_overall(:,3)],{'latency','duration'},'bad_ET');
            else
                EEG = addevents(EEG,[seq_bad_overall(:,1) seq_bad_overall(:,3:4)],{'latency','duration', 'epoch'},'bad_ET');
            end
            
        otherwise
            error('%s(): rejection method input not recognized',mfilename)
    end
end

end