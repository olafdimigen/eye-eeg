% pop_eyetrackerica() - reject independent components based on their 
%                       covariance with simultaneously recorded eye 
%                       tracking data
%
% Usage:
%   >> EEG = pop_eyetrackerica(EEG, sacstring, fixstring, sactol, ...
%              threshratio, flagmode, plotfig, topomode)
%
% Inputs:
%  EEG        - [string] EEG struct. EEG must already contain extra 
%               channels with synchronized eye tracking data
%               (>> pop_importeyetracker()) and both saccade and fixation
%               events (imported or by using pop_detecteyemovements)
%  sacstring  - [string] name of saccade event (in EEG.event).
%               This is 'saccade' per default, but may differ if saccades 
%               were imported from Eyelink raw data (e.g. 'L_saccade')
%  fixstring  - [string] name of fixation event (in EEG.event).
%               This is 'fixation' per default, but may differ if fixations 
%               were imported from Eyelink raw data (e.g. 'R_fixation')
%  sactol     - [vector of two integers], e.g. [5 10]. 
%               Extra temporal tolerance around saccade onset and offset.
%               Saccade intervals will range from saccade-onset minus 
%               window(1) samples until saccade-offset plus window(2) 
%               samples. Correspondingly, the tolerance will be substracted
%               from the fixation intervals.
%  threshratio- critical variance ratio threshold for component rejection.
%               If the variance of IC activity within saccades is more than
%               threshratio larger than the IC activity during fixations, the
%               component will be flagged for rejection.
%  flagmode   - [integer: either 1,2, or 3]
%               1: do not flag any components for rejection (test mode)
%               2: flag components for rejection in EEG.reject.gcompreject, 
%               keep existing rejection flags (e.g. those flagged manually)
%               3: flag components for rejection in EEG.reject.gcompreject,
%               overwrite all existing flags in EEG.reject.gcompreject, i.e.
%               components formely flagged as "bad" may be changed to "good"
%  plotfig    - <boolean> [0/1] - show figure to evaluate the rejection
%               threshold
%  topomode   - <integer: 1,2,3, or 4> - determines whether after 
%               eyetrackerica a figure is shown with the 2D topographies of 
%               flagged (bad and/or good) components (using pop_selectcomp). 
%               Figure allows to toggle rejection flags manually.
%               Use one of four options:
%               1: plot topos of "bad" and "good" components (2 figures)
%               2: plot topos of "bad" components only (1 figure)
%               3: plot topos of "good" components only (1 figure)
%               4: do not plot topographies
%               Note: Option 1:3 requires electrode coordinates.
%
% Outputs:
%   EEG       - EEG structure with updated field EEG.reject.gcompreject. 
%               If rejectmode 1 or 2 was used, EEG.reject.gcompreject will 
%               contain updated flags flagging components with a variance
%               ratio larger than threshratio 
%  vartable   - [matrix of size ncomp x 3]
%               Columns of [vartable] contain the following information:
%               vartable(:,1): mean IC variance during saccade intervals
%               vartable(:,2): mean IC variance during fixation intervals
%               vartable(:,3): ratio, i.e., vartable(:,1)/vartable(:,2)
%               High ratios indicate that the component is more active
%               during saccades than during fixations, and may therefore
%               reflect corneoretinal or myogenic ocular artifact
%
% To actually remove the flagged components from your data, use EEGLAB 
% menu "Tools" > "Remove Components" or function pop_subcomp() afterwards
%
% See also: eyetrackerica, geticvariance
%
% Example: A call of the function might look like this:
%
% >> [EEG vartable] = pop_eyetrackerica(EEG,'saccade','fixation',[10 5],1.1,3,1,1)
%
% Explanation: compute the variance of the activity time course of each IC
% within saccade events (of type 'saccade' in EEG.event) and fixation
% events (of type 'fixation' in EEG.event). For this analysis, saccade 
% intervals are defined as lasting from 10 samples before fixation onset 
% until 5 samples after saccade offset [sactol]. Conversely, variance for
% fixation intervals is computed from 5 samples after fixation onset until 
% 10 samples before fixation offset. If the variance ratio 
% [var(sac)/var(fix)] for a component is larger than 1.1, set this 
% component to "reject". Overwrite existing rejection flags already present
% in EEG.reject.gcompreject. Show an extra figure to evaluate the setting
% for the threshold (maxcrit). Plot 2D topographies of components flagged
% as good and bad.
%
% The saccade/fixation "variance ratio" criterion was proposed by:
%
% Plöchl, M., Ossandon, J.P., & König, P. (2012). Combining EEG and 
% eye tracking: identification, characterization, and correction of eye 
% movement artifacts in electroencephalographic data. Frontiers in Human 
% Neuroscience, doi: 10.3389/fnhum.2012.00278.
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

function [EEG, com] = pop_eyetrackerica(EEG,sacstring,fixstring,sactolerance,threshratio,flagmode,plotfig,topomode)

com = '';

% check whether ICA was already computed
if isempty(EEG.icaweights)
    fprintf('%s(): EEG.icaweights is empty. You need to run ICA before calling this function.',mfilename);
    return
end

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 8
        % pop up dialogue
        [sacstring,fixstring,sactolerance,threshratio,flagmode,plotfig,topomode] = dlg_eyetrackerica(mfilename,EEG);       
    end
    
    EEG = eyetrackerica(EEG,sacstring,fixstring,sactolerance,threshratio,flagmode,plotfig,topomode);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% return history command string
allArgs = vararg2str({sacstring,fixstring,sactolerance,threshratio,flagmode,plotfig,topomode});
com = sprintf('[EEG vartable] = %s(EEG,%s)',mfilename,allArgs);
return