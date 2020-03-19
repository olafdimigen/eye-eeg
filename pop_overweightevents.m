%% Create optimized ICA training data, in which the samples around
%% certain events (e.g. saccade onsets) are overweighted
%
% EEG = pop_overweightevents(EEG,event2overweight,sec_beforeevent,sec_afterevent,ow_proportion,removemean)
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
%  ow_proportion: length of the appendend overweighted samples expressed as a 
%                 proportion of the original length of the dateset 
%                (e.g. 0.5 means that if the orig. dataset was 1000 samples long,
%                we will append 500 samples of overweighted data, creating a final
%                dataset that is 1500 samples long)
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
% EEG_training = pop_overweightevents(EEG,'saccade',-0.02, 0.01, 0.5, 1)
%
% This example call would create a dataset that consists of all the samples 
% of the original EEG.data (in 2D) plus appended samples around saccade onsets. 
% These appended samples comprise the time interval from -20 ms to +10 ms
% relative to all 'saccade' events found in EEG.event. Samples around 
% saccades will be repeatedly re-appended to the end of the dataset until the
% appended dat corresponds to 50% (factor 0.5) of the original length of the dataset.
% The mean % channel voltage will be removed from each appended short epoch (recommended!).
%
% -------------------------------------------------------------------------
% The overweighting procedure was described and evaluated in:
%
% Dimigen, O. (2020). Optimizing ICA-based ocular artifact correction 
% of EEG data recorded during free viewing. NeuroImage, 207, 116117
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

function [EEG_overweighted, com] = pop_overweightevents(EEG,event2overweight,timelim,ow_proportion,removemean)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 5
        % pop up dialogue
        [event2overweight,timelim,ow_proportion,removemean] = dlg_overweightevents(mfilename,EEG);
    end
  
    [EEG_overweighted] = overweightevents(EEG,event2overweight,timelim,ow_proportion,removemean);
       
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% fprintf('\n\n---------------------------------------------------------------\n')
% fprintf('\nA variable \"EEG_overweighted\" was created in the MATLAB workspace\n')
% fprintf('\n---------------------------------------------------------------\n')
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_overweighted, size(ALLEEG,2),'setname','Overweighted training data','gui','off'); 

% return history command string
allArgs = vararg2str({event2overweight, timelim, ow_proportion, removemean});
com = sprintf('[EEG_overweighted] = %s(EEG,%s)',mfilename,allArgs);
return