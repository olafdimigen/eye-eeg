% dlg_applytochannels() - pops dialogue to enter channel numbers, called by 
%                         pop_applytochannels()
%                         for help type >> help pop_applytochannels
%
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

function [targetChannels] = dlg_applytochannels(callingFcn,EEG,func_call)

chanlocs = EEG.chanlocs;

try
    % default: select all channels that are not of type 'EYE'
    non_et_chans = find(~ismember({EEG.chanlocs.type},'EYE'));        
catch
    % select all channels
    non_et_chans = 1:EEG.nbchan;
end

brief_func_call = func_call(1:strfind(func_call,'(')-1);

geometry      = cell(1,3);
geometry(1:3) = { 1 1 [1 1 0.5]};

uilist = cell(1,5);
uilist(1:5) = {...
    { 'Style', 'text', 'string', sprintf('Apply function ''%s'' only to following channels:',brief_func_call),'fontweight', 'bold'},...
    {},...
    { 'Style', 'text', 'string', 'Channel number(s):'},...
    { 'Style', 'edit', 'string', sprintf('%s',vararg2str(non_et_chans)),'tag','chans'},...
    { 'Style', 'pushbutton', 'string', '...', 'enable' fastif(isempty(chanlocs), 'off', 'on') 'callback' 'tmpchanlocs = EEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); num2str(tmp), set(findobj(gcbf, ''tag'', ''chans''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' }
    };

results = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Apply EEGLAB function to selected channels -- ', callingFcn]);

if isempty(results)
    return
else
    [chaninds chanlist] = eeg_decodechan(EEG.chanlocs,results{1});
    targetChannels = chaninds;
end