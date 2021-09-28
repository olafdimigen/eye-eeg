% dlg_detecteogspikes - pops dialogue to enter parameters for the saccade
%                onset estimation based on spikes in the radial EOG
%                for help, type >> help pop_detecteogspikes
%
% Copyright (C) 2009-2021 Olaf Dimigen, HU Berlin
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

function [chans_eog,chans_parietal,threshfactor,mindist,clustermode,plotfig,writesac] = dlg_detecteogspikes(callingFcn,chanlocs,srate)

geometry = { 1 1 1 [2 0.8 0.5] [2 0.8 0.5] 1 1 [2 0.8 0.5] [2 0.8 0.5] [2 1.3] 1 1 [2 1.3] [2 1.3] };

%% strings for pulldown menu (saccade clustering methods, see mergesacc())
clustermethods = {'1. keep all spikes (do nothing)'; '2. keep first spike'; '3. keep largest spike'};

%% main menu
uilist = {...
    {'Style', 'text', 'string', 'Specify channels for computation of a ''radial'' EOG montage', 'fontweight', 'bold'},...
    {'Style', 'text', 'string', 'Recommended definition: mean of all facial EOG channels minus Pz'},...
    {},...
    {'Style', 'text', 'string', 'EOG channels:'},...
    {'Style', 'edit', 'string', '', 'tag', 'mychan_eog', 'TooltipString','do not add EOG channels recorded with a bipolar montage!' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); num2str(tmp); set(findobj(gcbf, ''tag'', ''mychan_eog''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' },...
    ...
    {'Style', 'text', 'string', 'Parietal electrode(s) (e.g., ''Pz''):'},...
    {'Style', 'edit', 'string', '', 'tag', 'mychan_parietal', 'TooltipString','you can also use P1/P2, POz, Oz' }, ...
    {'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); num2str(tmp); set(findobj(gcbf, ''tag'', ''mychan_parietal''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' },...
    {},...
    {'Style', 'text', 'string', 'Parameters for spike detection','fontweight', 'bold'},...
    ...
    {'Style', 'text', 'string', 'Threshold multiplier (in standard deviat. of filtered rEOG)'},...
    {'Style', 'edit', 'string', '2', 'tag', 'threshfactor'},...
    {'Style','text', 'string', 'SD'},...
    ...
    {'Style', 'text', 'string', 'For clusters of spikes separated by less than:'},...
    {'Style', 'edit','string', '50', 'tag', 'mindist_ms'},...
    {'Style', 'text', 'string', 'ms'},...
    ...
    {'Style', 'text', 'string', '...do the following:'},...
    {'Style', 'popupmenu', 'string',clustermethods,'tag','clustermode','value',3},...
    {},...
    {'Style', 'text', 'string', 'Add spikes to EEG.event?','fontweight', 'bold'},...
    ...
    {'Style', 'text', 'string', 'Add spikes as ''reog_spike'' events to EEG.event?'},...
    {'Style', 'checkbox','value', 0, 'string' '(uncheck to test settings first)','tag','writesac'},...
    ...
    {'Style', 'text', 'string', 'Plot figure showing detection results?'},...
    {'Style', 'checkbox','value', 1,'tag','plotfig'},...
    };


%% make GUI
[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Detect spikes in radial EOG -- ', callingFcn]);

%% process user input (cancel)
if isempty(results)
    disp('results is empty!')
    return
end

% get/transform outputs
chans_eog      = eeg_decodechan(chanlocs,outstruct.mychan_eog);
chans_parietal = eeg_decodechan(chanlocs,outstruct.mychan_parietal);
threshfactor   = str2num(outstruct.threshfactor);
mindist        = round(str2num(outstruct.mindist_ms)/(1000/srate));
clustermode    = outstruct.clustermode;
plotfig        = outstruct.plotfig;
writesac       = outstruct.writesac;

end