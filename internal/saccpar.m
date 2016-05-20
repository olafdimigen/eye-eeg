function sac = saccpar(bsac)
%-------------------------------------------------------------------
%
%  FUNCTION saccpar.m
%  Calculation of binocular saccade parameters;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 1.0, 01 AUG 05)
%
%-------------------------------------------------------------------
%
%  INPUT: binocular saccade matrix from FUNCTION binsacc.m
%
%  sac(:,1:14)       binocular microsaccades 
%
%  OUTPUT:
%
%  sac(:,1:9)        parameters avergaed over left and right eye data
%
%---------------------------------------------------------------------
if size(bsac,1)>0
    sacr = bsac(:,1:7);
    sacl = bsac(:,8:14);

    % 1. Onset
    a = min([sacr(:,1)'; sacl(:,1)'])';

    % 2. Offset
    b = max([sacr(:,2)'; sacl(:,2)'])';

    % 3. Duration
    DR = sacr(:,2)-sacr(:,1)+1;
    DL = sacl(:,2)-sacl(:,1)+1;
    D = (DR+DL)/2;

    % 4. Delay between eyes
    delay = sacr(:,1) - sacl(:,1);

    % 5. Peak velocity
    vpeak = (sacr(:,3)+sacl(:,3))/2;

    % 6. Saccade distance
    dist = (sqrt(sacr(:,4).^2+sacr(:,5).^2)+sqrt(sacl(:,4).^2+sacl(:,5).^2))/2;
    angle1 = atan2((sacr(:,5)+sacl(:,5))/2,(sacr(:,4)+sacl(:,4))/2);
    
    % Note added by olaf.dimigen@hu-berlin.de:
    % The origin of the eye tracking coordinate system is usually in the 
    % upper left screen corner, not in the lower left screen corner as in a
    % Cartesian coordinate system. Remember this when interpreting saccade 
    % angles (angle1 and angle2) that are based on the vertical movement
    % component.
    
    % 7. Saccade amplitude
    ampl = (sqrt(sacr(:,6).^2+sacr(:,7).^2)+sqrt(sacl(:,6).^2+sacl(:,7).^2))/2;
    angle2 = atan2((sacr(:,7)+sacl(:,7))/2,(sacr(:,6)+sacl(:,6))/2);
    
    sac = [a b D delay vpeak dist angle1 ampl angle2];
else
    sac = [];
end