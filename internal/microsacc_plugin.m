function [sac, radius] = microsacc_plugin(x,vel,VFAC,MINDUR,MSDX,MSDY)
%-------------------------------------------------------------------
%
%  FUNCTION microsacc.m
%  Detection of monocular candidates for microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 2.1, 03 OCT 05)
%
%-------------------------------------------------------------------
%
%  INPUT:
%
%  x(:,1:2)         position vector
%  vel(:,1:2)       velocity vector
%  VFAC             relative velocity threshold
%  MINDUR           minimal saccade duration

%  Optional inputs (modification by OD for plugin:)
%  MSDX             threshold for x-component (horizontal)
%  MSDY             threshold for y-component (vertical)

%  OUTPUT:
%
%  sac(1:num,1)   onset of saccade
%  sac(1:num,2)   end of saccade
%  sac(1:num,3)   peak velocity of saccade (vpeak)
%  sac(1:num,4)   horizontal component     (dx)
%  sac(1:num,5)   vertical component       (dy)
%  sac(1:num,6)   horizontal amplitude     (dX)
%  sac(1:num,7)   vertical amplitude       (dY)
%
%
%-------------------------------------------------------------------
% Note by olaf.dimigen@hu-berlin.de (2012):
% This version of microsacc() has been modified for use with the EYE-EEG toolbox. 
% msdx and msdy can now be provided as function inputs. Their computation 
% is outsourced into new function velthresh(), which contains code that was 
% formerly part of microsacc(). Helpful to compute threholds globally 
% (over all data epochs of a subject) or to apply fixed thresholds.
%-------------------------------------------------------------------

% // changes to original microsacc()
% // olaf.dimigen@hu-berlin.de
if nargin < 5
    [msdx msdy] = velthresh(vel);
else
    msdx = MSDX;
    msdy = MSDY;
end
% // end of changes

radiusx = VFAC*msdx;
radiusy = VFAC*msdy;
radius = [radiusx radiusy];

% compute test criterion: ellipse equation
test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;
indx = find(test>1);

% determine saccades
N = length(indx); 
sac = [];
nsac = 0;
dur = 1;
a = 1;
k = 1;
while k<N
    if indx(k+1)-indx(k)==1
        dur = dur + 1;
    else
        if dur>=MINDUR
            nsac = nsac + 1;
            b = k;
            sac(nsac,:) = [indx(a) indx(b)];
        end
        a = k+1;
        dur = 1;
    end
    k = k + 1;
end

% check for minimum duration
if dur>=MINDUR
    nsac = nsac + 1;
    b = k;
    sac(nsac,:) = [indx(a) indx(b)];
end

% compute peak velocity, horizonal and vertical components
for s=1:nsac
    % onset and offset
    a = sac(s,1); 
    b = sac(s,2); 
    % saccade peak velocity (vpeak)
    vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
    sac(s,3) = vpeak;
    % saccade vector (dx,dy)
    dx = x(b,1)-x(a,1); 
    dy = x(b,2)-x(a,2); 
    sac(s,4) = dx;
    sac(s,5) = dy;
    % saccade amplitude (dX,dY)
    i = sac(s,1):sac(s,2);
    [minx, ix1] = min(x(i,1));
    [maxx, ix2] = max(x(i,1));
    [miny, iy1] = min(x(i,2));
    [maxy, iy2] = max(x(i,2));
    dX = sign(ix2-ix1)*(maxx-minx);
    dY = sign(iy2-iy1)*(maxy-miny);
    sac(s,6:7) = [dX dY];
end