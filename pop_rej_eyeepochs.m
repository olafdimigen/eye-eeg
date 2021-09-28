% pop_rej_eyeepoched() - remove epochs that contain out-of-range values
%                  in (eye tracking) channels. Function pops a dialogue to
%                  enter min-max values and channels to check for
%                  out-of-range data. The actual exclusion of epochs is
%                  peformed by simply calling EEGLAB function
%                  pop_eegthresh(). Note: This function does not store
%                  rejection information (without rejecting), but 
%                  immediately rejects all bad epochs from the dataset
%
% Usage:
%   >> EEG = pop_rej_eyeepoched(EEG)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data
%
% Outputs:
%   EEG         - EEG structure
%
% Author: od
% Copyright (C) 2009-2021 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [EEG, com] = pop_rej_eyeepochs(EEG)

com = '';

% display help if not enough arguments
if nargin < 1
    help(mfilename);
    return;
end

% is epoched dataset loaded?
if isempty(EEG.data) || ~isfield(EEG,'epoch')
    error('\n%s(): This function works on epoched data. No epochated dataset is loaded!\n',mfilename)
end

try
    %% ask for channels to include
    if ~exist('targetchannels','var') || isempty(targetchannels) 
        
        % pop up target channel dialogue
        %[chans minvals maxvals whatever rejectnow] = dlg_rej_eyeepochs(EEG);
        [chans minvals maxvals] = dlg_rej_eyeepochs(mfilename, EEG);
    end
        
    chans
    minvals
    maxvals
    
    % execute actual function: use EEGLAB 
    rejectnow = 1;
    [EEG,rejectindex,com] = pop_eegthresh(EEG,1,chans,minvals,maxvals,EEG.times(1)/1000,EEG.times(end)/1000,0,rejectnow);

    % alternative:
    % call pop_eegthreshold with last input set to 0 but then execute pop_select()
    % EEG = pop_select(EEG,'notrials',find(rejectindex));
        
    % note: there is a serious bug in: EEG = pop_rejepoch( EEG, [1 1 1 1 1 1 1 1 ....] ,0);
    % when *all* epochs are BAD, pop_rejepoch() only rejects one trial (the first)
    % because binary indices are not converted to numeric indices via find() 
    % before the rejection-vector is handed to pop_select(EEG,'notrial',rejectionvector)
    
    % two options:
    % 1. EEG = pop_rejepoch( EEG, rejectindex, 0);
    % 2. EEG = eeg_rejsuperpose( EEG, 0, 1, 0, 0, 0, 0, 0, 0);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

return