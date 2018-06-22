%% Create optimized ICA training data, in which the samples around
%% certain events (e.g. saccade onsets) are overweighted
%
% EEG = overweightevents(EEG,event2overweight,sec_beforeevent,sec_afterevent,ow_proportion,removemean)
%
% INPUTS
%  EEG  - EEG dataset in EEGLAB format
%         can be epoched or continuous but must contain the event of
%         interest (e.g. 'saccade' events detected using EYE-EEG)
%
%  event2overweight - name of event in EEG.event to copy (e.g. 'saccade')
%
%  sec_beforeevent - time before the event onset to include (in seconds)
%
%  sec_afterevent  - time after the event onset to include (in seconds)
%
%  ow_proportion: proportion of the appendend overweighted samples relative
%                 to the current length of the dateset (e.g. 0.5)
% 
%  removemean - [boolean], for the appended epochs, remove mean value from each channel
%               across the whole epoch? (recommended)
%
% OUTPUTS
% - EEG with additional appended samples (=optimized ICA training dataset)
%
%  Please note
%    - The output EEG.data will be 2D (even if the input EEG.data was 3D)
%    - No additional events are added to EEG.event for the appended brief epochs
%
% Example call to this function:
% EEG_training = overweightevents(EEG,'saccade',-0.02, 0.01, 0.5, 1)
%
% This example call would create a dataset that consists of all the samples 
% of the original EEG.data (in 2D) plus appended samples around saccade onsets. 
% These appended samples comprise the time interval from -20 ms to +10 ms
% relative to all 'saccade' events found in EEG.event. Samples around 
% saccades will be repeatedly re-appended to the end of the dataset until the
% dataset is 0.5 times [ow_proportion] longer than before . That is, length of
% the newly created dataset will be 150% of its original length. The mean
% channel voltage will be removed from each appended short epoch.
%
% Note 1: Removing the epoch mean from the appended epochs is usually the best
% choice for ICA
%
% Note 2: Following recommendations in Dimigen (2018), the data should also 
% be optimally (high pass-) filtered for ICA *before* running this function
%
% -------------------------------------------------------------------------
%% The overweighting procedure was proposed and empirically evaluated in:
%
% Dimigen, O. (2018). Optimizing ICA-based ocular artifact correction 
% of EEG data recorded during free viewing. BioArXiv
%
% Please cite this reference paper if you use the method.
% -------------------------------------------------------------------------

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

function [EEG_overweighted] = overweightevents(EEG,event2overweight,timelim,ow_proportion,removemean)

fprintf('\n\n%s(): Creating data with overweighted artifacts for ICA training',mfilename);
fprintf('\n%s(): This may take a moment...', mfilename);
fprintf('\n%s(): Creating dataset with event *%s* overweighted.\n', mfilename, event2overweight);
fprintf('%s(): Appended samples will make up %.f percent of original length of dataset.\n\n', mfilename, event2overweight, ow_proportion*100);

if ~any(ismember({EEG.event.type},event2overweight))
    error('\n\n%s(): The overweighting event \"%s\" was not found in EEG.event.type.\nPlease check.',mfilename,event2overweight)
end


%% reshape original EEG.data to 2D format
EEG_tmp          = eeg_emptyset;
EEG_tmp.data     = EEG.data(:,:); % reshape data to 2D
EEG_tmp.pnts     = size(EEG_tmp.data,2);
EEG_tmp.srate    = EEG.srate;
EEG_tmp.chanlocs = EEG.chanlocs;
EEG_tmp.event    = EEG.event;
EEG_tmp.nbchan   = size(EEG_tmp.data,1);
EEG_tmp          = eeg_checkset(EEG_tmp);

%% how many overweight event samples to add?
npoints    = size(EEG_tmp.data,2); % number of EEG samples
nsacpoints = round(ow_proportion * npoints); % desired number of (redundant) saccade samples

%% create event-locked epochs to overweight
sac = pop_epoch(EEG,{event2overweight},timelim);

if removemean
    sac = pop_rmbase(sac,[]); % baseline subtracted across whole epoch
end
%sac = applytochannels(sac,1:NCHANS_EEG,'pop_rmbase( EEG,[]);');

%% overweight (=copy & re-append) overweight-event-locked epochs
tmpsacdata = sac.data(:,:); % 2D
tmpsacdata = repmat(tmpsacdata,1, ceil(nsacpoints./size(tmpsacdata,2)));
tmpsacdata = tmpsacdata(:,1:nsacpoints);
clear sac

%% create EEG dataset for overweight-event-locked epochs
EEG_sac          = eeg_emptyset;
EEG_sac.data     = tmpsacdata;
EEG_sac.pnts     = size(EEG_sac.data,2);
EEG_sac.srate    = EEG.srate;
EEG_sac.chanlocs = EEG.chanlocs;
EEG_sac.nbchan   = size(EEG_sac.data,1);
EEG_sac.event    = []; % overweighted events not in EEG.event structure
EEG_sac          = eeg_checkset(EEG_sac);

%% merge original EEG dastset with overweighted event-locked epochs
EEG_overweighted = pop_mergeset(EEG_tmp,EEG_sac,1);

fprintf('%s(): Done.\n', mfilename);