% ploteyemovements - generate a figure showing basic properties of detected
%                    saccades and fixations
%
% Usage:
%
% >> ploteyemovements(sac_amp,sac_vmax,sac_angle,...
%                     fix_dur,fix_posx,fix_posy,metric,polar_flipy)
%
% Required inputs:
%
%   sac_amp     - [colum vector] saccade amplitude [pixel or degree]
%   sac_vmax    - [colum vector] saccade peak velocity [pixel or degree per second]
%   sac_angle   - [colum vector] saccade direction [angle in degree]
%   fix_dur     - [colum vector] fixation duration [in ms]
%   fix_posx    - [colum vector] horiz. fixation location (monoc. or binocular avg.)[pixel]
%   fix_posy    - [colum vector] vert. fixation location (monoc. or binocular avg.) [pixel]
%
%
%   Note: all input vectors need to have the same length
%
% Optional inputs:
%
%   metric      - [string] name of spatial unit in which sac_amp and sac_vmax
%                 are provided (e.g., 'degree' or 'pixel'). Fixation locations
%                 are always assumed to be in pixels. Only used for plot
%                 labels.
%   polar_flipy - [boolean: 0/1]: controls whether y-axis of the polar plot 
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
%             - figure with saccade & fixation properties
%
% Warning: the directional histogram ("rose") plot may be confusing.
% With most eye trackers the origin of the X-Y coordinate system is in the
% upper left corner of the computer screen, rather than in the lower left
% screen corner (Cartesian). In this case, an upwards-saccade has a
% negative movement component (going from a larger to a smaller Y-values
% in pixels). This also affects the saccade angle.
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

function ploteyemovements(sac_amp,sac_vmax,sac_angle,fix_dur,fix_posx,fix_posy,metric,polar_flipy)

if ~exist('metric','var')
    fprintf('\n%s(): Assuming that values are provided in degrees of visual angle',mfilename)
    metric = 'degree';
end

if ~exist('polar_flipy','var')
    polar_flipy = true;
    fprintf('\n%s(): Assuming that eye tracker coordinate system was in UPPER left screen corner.\nFlipping vertical axis of polar plot showing saccade angles!',mfilename)
end

%% display settings
POINTSPERBIN_1D = 20; % bin width for bar plots (avg. cases per bin) 
POINTSPERBIN_2D =  5; % bin width for "heatmap" (avg. cases per bin) 
MARKER          =  3; % point size in scatterplots

%% plot saccades and/or fixations?
 
% all saccade properties entered?
plotsac = false;
if ~isempty(sac_amp) && ~isempty(sac_vmax) && ~isempty(sac_angle)
    plotsac = true; 
    n_obs = length(sac_amp);
end
% all fixation properties entered?
plotfix = false;
if ~isempty(fix_dur) && ~isempty(fix_posx) && ~isempty(fix_posy)
    plotfix = true; 
    n_obs = length(fix_dur);
end
if ~plotsac && ~plotfix
    error('\n%s(): Function expects three saccade properties and/or three fixation properties as inputs.',mfilename)
end

%% get number of bins for 1D and 2D histograms
if n_obs > POINTSPERBIN_1D
    nbin_1d = round(n_obs/POINTSPERBIN_1D);
    nbin_2d = round(sqrt(n_obs/POINTSPERBIN_2D));
else
    nbin_1d = n_obs;
    nbin_2d = round(sqrt(nbin_1d));
    if nbin_2d < 1, nbin_2d = 1; end
end

%% plot figure

if plotsac && plotfix
    figure('Name','Properties of saccades and fixations'); rows = 2; j = 3;
elseif plotsac
    figure('Name','Properties of saccades'); rows = 1; j = 0;
else
    figure('Name','Properties of fixations'); rows = 1; j = 0;
end

if plotsac    
    
    fprintf('\n%s(): Plotting properties of %i saccades',mfilename,length(sac_amp));
    %% 1: saccade amplitude distribution (histogram)
    subplot(rows,3,1); hold on; title('Saccades: Amplitude','fontweight', 'bold')
    [n edges] = hist(sac_amp,nbin_1d);
    bar(edges,n,'b')
    xlabel(sprintf('Saccade amplitude [%s]',metric));
    ylabel('Cases')
    box on
    % xlim([ prctile(sac_amp,2) prctile(sac_amp,98) ]) % do not plot extreme outliers
    
    %% 2: main sequence (scatterplot)
    subplot(rows,3,2);
    loglog(sac_amp,sac_vmax,'b.','markersize',MARKER)
    xlabel(sprintf('Sacc. amplitude [%s]',metric));
    ylabel(sprintf('Sacc. peak velocity [%s/s]',metric));
    box on
    title('Saccades: Main sequence','fontweight', 'bold')
    
    %% 3: saccade orientation (directional histogram)
    subplot(rows,3,3);
    [t,r] = rose(sac_angle*pi/180,36); % angle in radians, plot 10° bins
    h = polar(t,r,'b-');
    hline = findobj(gca,'Type','line');
    set(hline,'LineWidth',1.2); % make line thicker
    title('Saccades: Angular histogram','fontweight', 'bold')
    if polar_flipy
        set(gca,'ydir','reverse');
        fprintf('\n%s(): Vertical axis of polar plot (showing saccade angles) is flipped.',mfilename)        
    else
        fprintf('\n%s(): Vertical axis of polar plot (showing saccade angles) is not flipped.',mfilename)
    end
end

if plotfix
    fprintf('\n%s(): Plotting properties of %i fixations',mfilename,length(fix_dur));
    
    %% 4: fixation durations (histogram)
    subplot(rows,3,1+j); hold on; title('Fixations: Durations','fontweight', 'bold')
    [n edges] = hist(fix_dur,nbin_1d);
    bar(edges,n,'b')
    xlabel('Fixation duration [ms]');
    ylabel('Cases')
    box on
    %xlim([ prctile(fix(:,3),2) prctile(fix(:,3),98) ]) % do not show extreme outliers
    
    %% 6: fixation locations (heatmap)
    subplot(rows,3,2+j);
    % histnd replaces hist3 for users without the statistics toolbox
    [pointCount, dimCenters] = hist2d([fix_posx fix_posy],[nbin_2d nbin_2d]);
    imagesc(dimCenters{1},dimCenters{2},pointCount');
    title('Fixations: Heatmap','fontweight', 'bold')
    xlabel('Horizontal position [pix]');
    ylabel('Vertical position [pix]')
    colormap('hot')
    if polar_flipy
        % default: assume that origin is in upper left screen corner
        set(gca,'YDir','reverse'); 
    else
        % origin of coordinate system in lower left screen corner
        fprintf('\n%s(): Y-axis of fixation location plot is not reversed.',mfilename) 
    end
    set(gca,'color',[0 0 0]); % set "canvas" of subplot to black
    xlimits = get(gca,'xlim');
    ylimits = get(gca,'ylim');
    axis equal % x-y proportion
    axis tight
    % colorbar('SouthOutside');
    
    %% 5: fixation locations (scatterplot)
    subplot(rows,3,3+j); hold on; title('Fixations: Locations','fontweight', 'bold')
    % axis equal % correct x-y proportion
    plot(fix_posx,fix_posy,'b.','markersize',MARKER)
    xlabel('Horizontal position [pix]');
    ylabel('Vertical position [pix]')
    axis([xlimits ylimits]) % use same axis limits as in heatmap
    set(gca,'YDir', 'reverse');
    axis equal % x-y proportion
    axis tight
    box on
end