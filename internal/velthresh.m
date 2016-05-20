function [msdx msdy] = velthresh(vel)
% VELTHRESH - compute median-based standard deviation (SD) estimators
% of velocity time series
%
%-------------------------------------------------------------------
% INPUT:
%
% vel: n-by-2 matrix with hor. [vel(:,1)] and vert. [vel(:,2)] eye velocity
%
% OUTPUT:
%
% msdx: median-based SD estimator for horizontal eye velocity
% msdy: median-based SD estimator for vertical eye velocity
%-------------------------------------------------------------------
%
% Note by olaf.dimigen@hu-berlin.de (OD): The code below was originally 
% part of the microsacc() function by Engbert et al. but was outsourced 
% into a stand-alone function for usage with the plugin. This allows to
% compute velocity thresholds globally (based on all epochs of a given 
% participant) or to work with fixed rather than relative thresholds.

% compute threshold
msdx = sqrt( median(vel(:,1).^2) - (median(vel(:,1)))^2 );
msdy = sqrt( median(vel(:,2).^2) - (median(vel(:,2)))^2 );
if msdx<realmin
    msdx = sqrt( mean(vel(:,1).^2) - (mean(vel(:,1)))^2 );
    if msdx<realmin
        %error('msdx<realmin in microsacc.m');
        error('msdx<realmin in vecthresh.m. Did you exclude blinks/missing data before saccade detection?'); % // added by O.D.
    end
end
if msdy<realmin
    msdy = sqrt( mean(vel(:,2).^2) - (mean(vel(:,2)))^2 );
    if msdy<realmin
        %error('msdy<realmin in microsacc.m');
        error('msdy<realmin in vecthresh.m. Did you exclude blinks/missing data before saccade detection?'); % // added by O.D.
    end
end

end