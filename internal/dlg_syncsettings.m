% dlg_syncsettings - pops dialogue, called by pop_importeyetracker
%                    see >> help pop_importeyetracker
% 
% Copyright (C) 2009-2020 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [syncEvents, importColumns, channelLabels, importEyeEvents, plotfig] = ...
    dlg_syncsettings( callingFcn, commonEvents, nEEG, nET, selectInPullDown, nColumns, importColumns, channelLabels, exampleData, hasEyeEvents, importEyeEvents)

% geometry depends on number of data columns in raw data
geometry            = cell(1,9+nColumns+2);
geometry(1:9)       = { 1 [2 1.5] [2 1.5] [2 1.5] 1 [2 1.5] 1 1 [1 2.5 1] };
geometry(10:end-2)  = {[1 2.5 1]};
geometry(end-1:end) = {1 [2 1.5]};

if hasEyeEvents
    importEnabled = 'on';
    importMessage = ' (check to import events)';
else
    importEnabled = 'off';
    importMessage = ' (no such events found)';
end

% strings for pulldown menu (synchronization events)
eventString = [num2str(commonEvents) repmat('  (EEG: ',length(commonEvents),1) num2str(nEEG') repmat('x ET: ',length(commonEvents),1) num2str(nET') repmat('x)',length(commonEvents),1)];

uilist = cell(1,15 + 3*nColumns + 3);
uilist(1:15) = {...
    { 'Style', 'text', 'string', 'Select events used for synchronization','fontweight','bold',},...
    {},...
    { 'Style', 'text', 'string', 'Event type (# in EEG, # in ET)'},...
    { 'Style', 'text', 'string', 'Start-event (begin synch. at first of these):'},...
    { 'Style', 'popupmenu', 'string', eventString, 'tag','startTrigger','value',selectInPullDown(1)},...
    { 'Style', 'text', 'string', 'End-event (end synch. at last of these):'},...
    { 'Style', 'popupmenu', 'string', eventString,'tag','endTrigger','value',selectInPullDown(2)},...
    {},...
    { 'Style', 'text', 'string', 'Show plot with synchronization results?'}, ...
    { 'Style', 'checkbox','value',1}, ...
    {},...
    { 'Style', 'text', 'string', 'Choose data columns to import','fontweight', 'bold',},...
    { 'Style', 'text', 'string', 'Example data'}, ...
    { 'Style', 'text', 'string', 'Name of data column'}, ...
    { 'Style', 'text', 'string', 'Import? (y/n)'}, ...
    };

for col = 1:nColumns
    uilist((16:18) + 3*(col-1)) = {
        {'Style', 'text', 'string', sprintf('%g',exampleData{col} ) }
        {'Style', 'edit', 'string', channelLabels(col) }
        {'Style', 'checkbox','value', ismember(col,importColumns)}
        };
end

uilist(end-2:end) = {...
    {},...
    {'Style', 'text', 'string', 'Import eye movement events from raw data?','HorizontalAlignment','center','enable',importEnabled },...
    {'Style', 'checkbox','string',importMessage,'HorizontalAlignment','center','value',importEyeEvents,'enable',importEnabled},...
    };

results = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Import and synchronize eye tracking data -- ', callingFcn]);

if isempty(results)
    return
end

% get results
syncEvents       = commonEvents([results{[1,2]}])';
plotfig          = results{3};
importColumns    = find( [ results{5:2:end-1} ] );
channelLabels    = [results{4:2:end-1}];
channelLabels    = channelLabels(importColumns);
importEyeEvents  = results{end};