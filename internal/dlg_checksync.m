% dlg_checksync - pops dialogue to enter channel names of horizontal ET
%                 gaze channel as well as channels containing the signal of
%                 the left and right horizontal EOG electrode
%                 called by pop_checksync()
%                 for help, type >> help pop_checksync
%
% Copyright (C) 2009-2016 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [gaze_x, heog_l, heog_r, plotfig] = dlg_checksync(callingFcn,chanlocs)

geometry = { 1 [2 0.8 0.5] [2 0.8 0.5] [2 0.8 0.5] [2 1.3] };

%% main menu
uilist = {...
    {'Style', 'text', 'string', 'Specify channels containing horiz. eye-track and horiz. EOG electrodes', 'fontweight', 'bold'},...
    {'Style', 'text', 'string', 'Eye-tracker channel with horiz. gaze position (X):'},...
    {'Style', 'edit', 'string', '', 'tag', 'chan_gaze_x', 'TooltipString','Pick channel that contains the horizontal gaze position of one of the eyes' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chan_gaze_x''), ''string'',tmp); clear tmp tmpchanlocs tmpval'},...
    ...
    {'Style', 'text', 'string', 'Horizontal EOG: LEFT electrode:'},...
    {'Style', 'edit', 'string', '', 'tag', 'chan_eog_l', 'TooltipString','Pick EOG (or frontolateral EEG) channels place left of the eyes' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chan_eog_l''), ''string'',tmp); clear tmp tmpchanlocs tmpval' },...
    ...
    {'Style', 'text', 'string', 'Horizontal EOG: RIGHT electrode:'},...
    {'Style', 'edit', 'string', '', 'tag', 'chan_eog_r', 'TooltipString','Pick EOG (or frontolateral EEG) channels place right of the eyes' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chan_eog_r''), ''string'',tmp); clear tmp tmpchanlocs tmpval' },...
    ...
    {'Style', 'text', 'string', 'Plot figure with cross-correlation function?'},...
    {'Style', 'checkbox','value', 1,'tag','plotfig'},...
    };

%% make GUI
[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Check synchronization -- ', callingFcn]);

%% process user input (cancel)
if isempty(results)
    return
end

gaze_x  = eeg_decodechan(chanlocs,outstruct.chan_gaze_x);
heog_l  = eeg_decodechan(chanlocs,outstruct.chan_eog_l);
heog_r  = eeg_decodechan(chanlocs,outstruct.chan_eog_r);
plotfig = outstruct.plotfig;

end