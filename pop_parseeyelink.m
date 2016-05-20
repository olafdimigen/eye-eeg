% pop_parseeyelink() - parses text-converted (ascii) Eyelink raw data and 
%       saves it in a matlab structure. Will import messages and events, 
%       if detected. Builds table of events from trigger pulses sent 
%       during the experiment. If the user provides a keyword that marks 
%       the keyword-containing messages as synchronisation events, these
%       messages will be used to build the event table instead.
%       For Eyelink datasets, eye movement events detected online by the 
%       eye tracker (blinks, saccades, fixations) are also saved into the 
%       output MATLAB structure "ET".
%
% Usage:
%   >> ET = pop_parseeyelink
%
% Inputs:
%           - no inputs allowed
% Outputs:
%   ET      -  Eyetracker data structure. Note: This gets deleted
%           immediately, when used via EEGLAB menu. You will only need the
%           mat-file location in the following pop_importeyetracker.
%
%   An example use of the function might look like this: 
%   >> ET = parseeyelink
%
% See also:
%   parseeyelink, pop_importeyetracker, pop_parsesmi
%
% Author: ur
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

function [et, com] = pop_parseeyelink

com = '';

try
    
    %% pop up dialogue
    [fileToLoad, matToSave, useMessageTrigger, triggerKeyword] = dlg_parser(mfilename);
    
    %% check for identical load/save filenames
    if strcmp(fileToLoad,matToSave)
        error('%s(): Input and output file names should not be identical.',mfilename)
    end
    
    %% call parser
    if ~useMessageTrigger
        et = parseeyelink( fileToLoad, matToSave);
    else
        et = parseeyelink( fileToLoad, matToSave, triggerKeyword);
    end
  
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        et = [];
        return
    else
        rethrow(err);
    end
end

%% return the string command for history
if ~useMessageTrigger
    allArgs = vararg2str({fileToLoad, matToSave});
else
    allArgs = vararg2str({fileToLoad, matToSave, triggerKeyword});
end
com = sprintf('ET = %s(%s);','parseeyelink',allArgs);
return