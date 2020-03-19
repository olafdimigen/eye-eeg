% dlg_calVisAngle - pops dialogue to compute visual angle per screen pixel
%            called by pop_detecteyemovements()
%            for help, type >> help pop_detecteyemovements
%
% Author: od
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

function alpha_per_pix = dlg_calcvisangle(callingFcn)

geometry = {1 1 [3 1 1] 1 [3 1 1] 1 [3 1 1]};

%% menu
uilist = {...
    {'Style', 'text', 'string', 'Compute visual angle of one screen pixel:','fontweight', 'bold'},...
    {},...
    {'Style', 'text', 'string', 'Horizontal width of screen  [in millimeters]'},...
    {'Style', 'edit', 'string', '400','tag','xmm'},...
    {'Style', 'text', 'string', 'mm'},...
    {},...
    {'Style', 'text', 'string', 'Horizontal resolution of screen [in pixels]'},...
    {'Style', 'edit', 'string', '1024','tag','xres' },...
    {'Style', 'text', 'string', 'pixel'},...
    {},...
    {'Style', 'text', 'string', 'Viewing distance [in millimeters]'},...
    {'Style', 'edit', 'string', '700','tag','distmm'},...
    {'Style', 'text', 'string', 'mm'},...
    };

[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Compute visual angle per pixel -- ', callingFcn] );

%% compute visual angle per screen pixel
if ~isempty(results)
    screenwidth = str2num(outstruct.xmm);
    resolution  = str2num(outstruct.xres);
    viewingdist = str2num(outstruct.distmm);
    
    try
        % calculate angle
        mm_per_pix    = screenwidth/resolution;
        alpha_per_pix = 180/pi*(2*atan(mm_per_pix/viewingdist/2));
        % fprintf('\n\nVisual angle per screen pixel is %f degrees.\n',alpha_per_pix);
    catch
        alpha_per_pix = [];
    end
else
    alpha_per_pix = [];
    return
end


end