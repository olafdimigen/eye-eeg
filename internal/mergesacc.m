% mergesacc() - handle temporal saccade clusters
%
% Usage:
%    >> sac = mergesacc(sac,x,mindist,mode)
%
% Inputs:
%    sac      output of saccpar (size: nsac * 9 columns)
%    x        eye position of corresponding data epoch
%    mindist  minimum required distance (in samples) between detected
%             movements in order to treat them as separate saccades
%             (i.e. minimum allowed fixation duration)
%    mode     [1,2,3, or 4] Specifies way in which saccade clusters
%             are treated:
%             1: do nothing, keep all saccades
%             2: only first saccade of each cluster is kept
%             3: only largest saccade of each cluster is kept
%             4: saccades of the cluster are merged into one (longer) sacc.
%
% Outputs:
%    sac      modified sac, changed according to "mode"
%             % columns of [sac]:
%             1: saccade onset (sample)
%             2: saccade offset (sample)
%             3: duration (samples)
%             4: NaN (originally: delay between eyes in samples)
%             5: vpeak (peak velocity)
%             6: saccade distance
%             7: saccade angle (based on distance)
%             8: saccade amplitude
%             9: saccade angle (based on amplitude)
%
% Note: This function is not part of the original algorithms proposed by
% Engbert et al. (2003, 2006). It intends to solve the problem that
% single eye movements may be detected as multiple saccades by microsacc,
% binsacc, and saccpar. For example, overshoots towards the end of saccades
% ("glissades") are frequently detected as a separate saccade events, even
% though they are probably better regarded as part of the same movement. 
% One consequence of this can be extremely short fixation durations 
% (e.g. < 50 ms).
%
% Note: This function is experimental.
%
% See also:
% smoothdata, vecvel, velthresh, microsacc_plugin, binsacc, saccpar
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

function sac = mergesacc(sac,x,mindist,mode)

if nargin < 4
    error('\n%s() expects four inputs',mfilename)
end

% column with sacc. size information
% col. 6: distance  (sacc. endpoint minus startpoint)
% col. 8: amplitude (max minus min in full sacc. trajectory)
AMPCOL = 6;

nsac  = size(sac,1);
keep  = true(nsac,1); % vector of saccades to keep

% determine temporal distance to preceding saccade
if nsac > 1
    tempdist = [NaN; sac(2:end,1)-sac(1:end-1,2)];
else
    tempdist = NaN;
end

% get index of saccades that form temporal clusters
ix_close = find(tempdist < mindist);

if ~isempty(ix_close)
    
    % get beginning, end, length of saccade clusters
    clustlist      = findsequence2(ix_close);
    clustlist(:,1) = clustlist(:,1)-1; % include first saccade
    clustlist(:,3) = clustlist(:,3)+1; % update cluster length
    
    % go tru saccade clusters
    for n = 1:size(clustlist,1)
        
        ix = clustlist(n,1):clustlist(n,2);
        
        % how to treat saccade clusters?
        switch mode
            case 1
                % do nothing
                return
                
            case 2
                % keep first saccade of each cluster
                keep(ix(2):ix(end)) = 0;
                
            case 3
                % keep largest saccade of each cluster
                keep(ix) = 0;
                [tmp,ix_max] = max(sac(ix,AMPCOL));
                keep(ix(ix_max)) = 1;
                
            case 4
                % merge all saccades of cluster into one movement
                newsac(1) = sac(ix(1),1);           % onset
                newsac(2) = sac(ix(end),2);         % offset
                newsac(3) = newsac(2)-newsac(1)+1;  % duration
                newsac(4) = NaN;                    % originally: delay between eyes
                
                % saccade peak velocity (vmax), here: max of individual vmax'es of cluster saccades
                newsac(5) = max(sac(ix,5));
                
                % alternative: go back to raw data to get overall vmax
                % newsac(5) = max( sqrt( vel(newsac(1):newsac(2),1).^2 + vel(newsac(1):newsac(2),2).^2 ) );
                
                % saccade distance (dx,dy) & angle (based on distance)
                dx = x(newsac(2),1)-x(newsac(1),1);
                dy = x(newsac(2),2)-x(newsac(1),2);
                newsac(6) = sqrt(dx.^2+dy.^2);
                newsac(7) = atan2(dy,dx); % angle1
                
                % saccade amplitude (dX,dY) & angle (based on amplitude)
                i = newsac(1):newsac(2);
                [minx, ix1] = min(x(i,1));
                [maxx, ix2] = max(x(i,1));
                [miny, iy1] = min(x(i,2));
                [maxy, iy2] = max(x(i,2));
                dX = sign(ix2-ix1)*(maxx-minx);
                dY = sign(iy2-iy1)*(maxy-miny);
                newsac(8) = sqrt(dX.^2+dY.^2);
                newsac(9) = atan2(dX,dY); % angle2
                
                % insert new saccade into table
                sac(ix(1),1:9)      = newsac;
                keep(ix(1))         = 1;
                keep(ix(2):ix(end)) = 0; % delete other sacc. of cluster
                
            otherwise
                error('\n%s(): Cluster mode not recognized. 1, 2, 3, or 4 expected.',mfilename)
        end
    end % cluster loop
    
    % remove other saccades of clusters
    sac = sac(keep,:);
end