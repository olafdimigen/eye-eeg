function [sac, monol, monor] = binsacc(sacl,sacr);
%-------------------------------------------------------------------
%
%  FUNCTION binsacc.m
%  Detection of binocular microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 2.0, 18 JUL 05)
%
%-------------------------------------------------------------------
%
%  INPUT: saccade matrices from FUNCTION microsacc.m
%
%  sacl(:,1:7)       microsaccades detected from left eye
%  sacr(:,1:7)       microsaccades detected from right eye
%
%  OUTPUT:
%
%  sac(:,1:14)       binocular microsaccades (right eye/left eye)
%  monol(:,1:7)      monocular microsaccades of the left eye
%  monor(:,1:7)      monocular microsaccades of the right eye
%
%---------------------------------------------------------------------

if size(sacr,1)*size(sacl,1)>0
    
    % determine saccade clusters
    TR = max(sacr(:,2));
    TL = max(sacl(:,2));
    T = max([TL TR]);
    s = zeros(1,T+1);
    for i=1:size(sacl,1)
        s(sacl(i,1)+1:sacl(i,2)) = 1;
    end
    for i=1:size(sacr,1)
        s(sacr(i,1)+1:sacr(i,2)) = 1;
    end
    s(1) = 0;
    s(end) = 0;
    m = find(diff(s~=0));
    N = length(m)/2;
    m = reshape(m,2,N)';
    
    % determine binocular saccades
    NB = 0;
    NR = 0;
    NL = 0;
    sac = [];
    monol = [];
    monor = [];
    for i=1:N
        l = find( m(i,1)<=sacl(:,1) & sacl(:,2)<=m(i,2) );
        r = find( m(i,1)<=sacr(:,1) & sacr(:,2)<=m(i,2) );
        if length(l)*length(r)>0
            ampr = sqrt(sacr(r,6).^2+sacr(r,7).^2);
            ampl = sqrt(sacl(l,6).^2+sacl(l,7).^2);
            [h, ir] = max(ampr);
            [h, il] = max(ampl);
            NB = NB + 1;
            sac(NB,:) = [sacr(r(ir),:) sacl(l(il),:)];
        else
            % determine monocular saccades
            if length(l)==0
                NR = NR + 1;
                monor(NR,:) = sacr(r,:);
            end
            if length(r)==0
                NL = NL + 1;
                monol(NL,:) = sacl(l,:);
            end
        end
    end
else
    % special cases of exclusively monocular saccades
    if size(sacr,1)==0
        sac = [];
        monor = [];
        monol = sacl;
    end
    if size(sacl,1)==0
        sac = [];
        monol = [];
        monor = sacr;
    end
end