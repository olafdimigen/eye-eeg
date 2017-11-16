% applytochannels() - "wrapper" function that allows to apply
%                existing EEGLAB functions to only a subset of channels.
%                Only some EEGLAB functions (e.g. runICA) can be applied to
%                subsets of channels, while others (filtering, baseline
%                subtraction) are applied to all channels if called via the 
%                GUI. 
%                Many operations (e.g. baseline subtraction) make 
%                sense for EEG/MEG data, but not necessarily for peripheral 
%                data channels (e.g. gaze position, pupil, ECG, GSR...).
%
% Usage:
% >> [EEG inner_com] = applytochannels(EEG,targetchans,func_call)
%
% Inputs:
%   EEG          - [string] EEG struct
%   targetchans  - indices of channels to which function shall be applied
%   func_call    - EEGLAB function call, see below for examples
%
% Outputs:
%   EEG         - EEG struct
%   inner_com   - command string for EEGLAB history functionality (eegh)
%
%
% Example usage of the function might look like this:
%
%   >> EEG = applytochannels(EEG,[1:32],'pop_rmbase(EEG,[-100 0]);');
%
%   In this example, channels 1 to 32 are baseline corrected using EEGLAB
%   function pop_rmbase(), but not the other (non-EEG) channels
%
%   >> EEG = applytochannels(EEG,[1:32],'pop_eegfilt(EEG,0,80,[],[0],0,0,''fir1'',0);');
%
%   In this example, function pop_eegfilt() is applied only to channels
%   1 to 32 (the EEG channels) of a dataset, the other non-EEG
%   channels are not filtered. Make sure to properly espace special
%   characters like (')
%
% Author: od
% Copyright (C) 2009-2017 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [EEG inner_com] = applytochannels(EEG,targetchans,func_call)

fprintf('\n%s(): applying function ''%s'' to channels: %s\n',mfilename,func_call,vararg2str(targetchans))

inner_com = '';

size(EEG.data)

% remember original data
EEGorg = EEG;

% delete non-target channels
otherchans = ~ismember(1:EEG.nbchan,targetchans);

if size(EEG.data,3) == 1 % continuous data
    EEG.data(otherchans,:) = [];
else % epoched data
    EEG.data(otherchans,:,:) = [];
end

% update EEG.nbchan
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs(otherchans) = []; % avoid warning in console

%% execute the "inner" (wrapped) EEGLAB function & update "inner_com"
if ~isempty(strfind(func_call,'timtopo'))
    % functions with 1 output
    func_call_full = sprintf('[inner_com] = %s',func_call);
    eval(func_call_full);
else
    % functions with 2 outputs
    func_call_full = sprintf('[EEG inner_com] = %s',func_call);
    eval(func_call_full);
end

% eegfilt() returns 3D data (chan x time x epoch) as 2D (chan x (time/epoch))
% EEGLAB restores three-dimensionality by running eeg_checkset()
EEG = eeg_checkset(EEG);

%% join modified targetchannels with the original (untouched) channels
% continuous or epoched data?
if size(EEGorg.data,3) == 1
    EEGorg.data(targetchans,:) = EEG.data(:,:);
else
    EEGorg.data(targetchans,:,:) = EEG.data(:,:,:);
end

EEG = EEGorg;
EEG.nbchan = size(EEG.data,1);