% pop_applytochannels() - "wrapper" function that serves to apply existing 
%                EEGLAB functions to a subset of channels only.
%                Only some EEGLAB functions (e.g. runICA) can be applied to 
%                subsets of channels, while others (filtering, baseline 
%                subtraction) are not. However, many operations 
%                (e.g. baseline subtraction) make sense for EEG/MEG data, 
%                but not necessarily for peripheral data channels (e.g. 
%                gaze position, pupil diameter, ECG, GSR).
%
% Usage:
% >> [EEG com] = pop_applytochannels(EEG,targetchans,func_call)
%
% Required inputs:
%   EEG          - [string] EEG struct
%   targetchans  - indices of channels to which function shall be applied
%   func_call    - EEGLAB function call (omit anything before actual 
%                  function name, e.g. omit "EEG = ")
% Outputs:
%   EEG         - EEG struct
%   com         - command string for EEGLAB history functionality (eegh)
%                 this strings includes both outer and inner function call
%
% An example call of the function might look like this: 
%
%   >> EEG = pop_applytochannels(EEG,[1:32],'pop_rmbase(EEG);');
%
%   Channels 1 to 32 are baseline corrected using EEGLAB
%   function pop_rmbase(), the other (non-EEG) channels are not touched.
%
%   >> EEG = pop_applytochannels(EEG,[1:32],'pop_eegfilt(EEG);');
%
%   Function pop_eegfilt() is applied only to channels 1 to 32 (the EEG 
%   channels) of a dataset, the other (non-EEG) channels are not filtered
%
% Author: od
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

function [EEG, com] = pop_applytochannels(EEG,targetchannels,func_call)

com = '';

if nargin < 3
    help(mfilename);
    return;
end

% is this the official way to check whether dataset is loaded?
if isempty(EEG.data)
    fprintf('\npop_applyToSomeChannels(): No dataset loaded!\n')
    return;
end

try
    
    %% ask for channels to include
    if ~exist('targetchannels','var') | isempty(targetchannels)
        
        % pop up the targetchannel dialogue
        [targetchannels] = dlg_applytochannels(mfilename,EEG,func_call);
    end
    
    % run wrapper function
    [EEG inner_com] = applytochannels(EEG,targetchannels,func_call);
    
    %% return com string for command history
    % (including the com string of the embedded 'inner' EEGLAB function)
    
    if ~isempty(inner_com) % inner_com is empty if user canceled menu of embedded EEGLAB function

        % replace single quotes (') in inner command (e.g. 'firls') with double quotes ('')
        inner_com = strrep(inner_com,'''',''''''); % seriously crazy syntax!
        
        % remove part of assignment before equal sign (i.e.,: "[EEG =")
        temp = strfind(inner_com,'=');
        if ~isempty(temp)
            inner_com = inner_com(temp(1)+1:end);
        end        
        com = sprintf('%s = applytochannels(%s, %s,''%s'');',inputname(1),inputname(1),vararg2str(targetchannels),inner_com);
    end
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

return