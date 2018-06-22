% pop_ploteyemovements() - plot basic propertiese of saccades and fixations
%
% Usage:
%   >> EEG = pop_ploteyemovments(EEG)
%
% Inputs:
%  EEG         - [string] EEG struct. EEG must already contain detected
%                saccades and/or fixations in EEG.event, either imported from
%                the raw data (>> pop_importeyetracker()) or detected later
%                using pop_detecteyemovements
%
%  sacstring   - [string] name of saccade event (in EEG.event).
%                This is 'saccade' per default, but may differ if saccades
%                were imported from Eyelink raw data (e.g. 'L_saccade')
%  fixstring   - [string] name of fixation event (in EEG.event).
%                This is 'fixation' per default, but may differ if fixations
%                were imported from Eyelink raw data (e.g. 'R_fixation')
%
%  metric      - [string] name of unit in which saccade amplitudes are stored
%                in EEG.event. Only used for plot labels.
%                (e.g., 'degree', 'pixels').
%
% Optional Inputs (command line only):
%
%  polar_flipy - [boolean: 0/1]: controls whether y-axis of the polar plot 
%                 showing saccade angles is flipped or not. It the ET 
%                 coordinate system has its origin in the upper left screen 
%                 corner, then an upward saccade has a *negative* vertical 
%                 movement component. Flipping the y-axis will then make 
%                 angles in the plot correspond to "real space". Plot labels
%                 with the angles (e.g. "90°") are unaffected by this.
%                 0: do not flip y-axis in polar plot
%                 1: flip y-axis in polar plot {default}
%
% Outputs:
%            - figure with saccade & fixation properties
%              
% See also: ploteyemovements
%
% Example: Calling this function might look like this:
% >> pop_ploteyemovements(EEG,'saccade','fixation',0)
%
% Author: od
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

function [EEG,com] = pop_ploteyemovements(EEG,sacstring,fixstring,metric,polar_flipy)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

if ~exist('polar_flipy','var')
    polar_flipy = true;
end


try
    if nargin < 4
        % ask how saccade/fixation events are named
        [sacstring fixstring metric] = dlg_ploteyemovements(mfilename,EEG);
    end
    
    %% check whether eye movement events are present in EEG.event
    if ~isempty(sacstring)
        ix_sac = find(cellfun(@(x) strcmp(x,sacstring),{EEG.event.type})); % bugfix v.0.337
    else
        ix_sac = [];
    end
    if ~isempty(fixstring)
        ix_fix = find(cellfun(@(x) strcmp(x,fixstring),{EEG.event.type})); % bugfix v.0.337
    else
        ix_fix = [];
    end
       
    %% collect eye movement properties from EEG.event    
    % for saccades
    if any(ix_sac)
        sac_amp   = [EEG.event(ix_sac).sac_amplitude]';
        sac_vmax  = [EEG.event(ix_sac).sac_vmax]';    
        
        % special case for saccade angle
        try % EM detected by plugin
            sac_angle = [EEG.event(ix_sac).sac_angle]';
        catch
            % EM imported from Eyelink (no angle info)
            fprintf('\n%s(): Computing saccade angles for Eyelink data...',mfilename);
            x_delta   = ([EEG.event(ix_sac).sac_endpos_x]-[EEG.event(ix_sac).sac_startpos_x])';
            y_delta   = ([EEG.event(ix_sac).sac_endpos_y]-[EEG.event(ix_sac).sac_startpos_y])';
            sac_angle = atan2(y_delta,x_delta) * 180/pi;
        end
    end
    
    % for fixations
    if any(ix_fix)
        % Bugfix, May, 2017, od: EEG.event.durations are in samples, not ms:
        fix_dur   = [EEG.event(ix_fix).duration]' .* (1000/EEG.srate); 
        fix_posx  = [EEG.event(ix_fix).fix_avgpos_x]';
        fix_posy  = [EEG.event(ix_fix).fix_avgpos_y]';
    end
    
    %% plot figure with saccades and/or fixations
    if any(ix_sac) && any(ix_fix)
        ploteyemovements(sac_amp,sac_vmax,sac_angle,fix_dur,fix_posx,fix_posy,metric,polar_flipy)
    elseif any(ix_sac)
        ploteyemovements(sac_amp,sac_vmax,sac_angle,[],[],[],metric,polar_flipy);
    elseif any(ix_fix)
        ploteyemovements([],[],[],fix_dur,fix_posx,fix_posy,metric,polar_flipy);
    else
        error('\n%s(): No saccades/fixations of specified type(s) found in EEG.event.\nDid you import or detect eye movements? (\"Eyetracker\" > \"Detect saccades & fixations\").',mfilename);
        return
    end
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% return history command string
allArgs = vararg2str({sacstring,fixstring,metric});
com = sprintf('%s(EEG,%s)',mfilename,allArgs);
return