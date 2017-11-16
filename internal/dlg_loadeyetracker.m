% dlg_loadeyetracker - pops dialogue to enter file location of the (parsed)
%                      eye tracker raw data file (.mat)
%                      called by pop_importeyetracker
%                      for help, type >> help pop_importeyetracker
%
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

function [matFileToLoad] = dlg_loadeyetracker(callingFcn)

cb_selectSource = ['[filename filepath] = uigetfile2(''*.mat'', ''Select eyetracker mat-file - ',callingFcn,'()'');'...
    'if filename(1) ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', tagToUpdate), ''string'', [ filepath filename ]);' ...
    'end;' ...
    'clear filename filepath tagToUpdate' ...
    ];

geometry = {[3 5 2]};

uilist = { ...
    { 'Style', 'text', 'string', 'Eyetracker .mat file', 'horizontalalignment', 'right'}, ... % 'fontweight', 'bold'
    { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'fileToLoad' }, ...
    { 'Style', 'pushbutton', 'string', 'Browse...', 'callback', ['tagToUpdate = ''fileToLoad'';' cb_selectSource] },...
    };

[results newcomments err outstruct] = inputgui( geometry, uilist, ['pophelp(''',callingFcn,''');'], ['Load eyetracking data -- ',callingFcn,'()']);

if isempty(results) || isempty(results{1})
    return
end

matFileToLoad = outstruct.fileToLoad;