function [safeVal,t,u] = checkSafeIntersection(p1,p2,q1,q2)
%p1,p2 are endpoints of line P.  Likewise for q1/q2/Q.
%returns TRUE if NO INTERSECTION
%returns FALSE if YES INTERSECTION

%simpler than orientation checking:
% https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
s=q2-q1;
r=p2-p1;
u=0; t=0;

uDen=cross2(r,s);
uNum=cross2(q1-p1,r);

if abs(uDen) <= 1e-4
    if abs(uNum)<1e-4
        safeVal=0; %colinear
        fprintf('Colinear')
    else
        safeVal=1;
    end
else
    tNum=cross2(q1-p1,s);
    t=tNum/uDen;
    u=uNum/uDen;
    if t<=1 && t>0 && u<=1 && u>0
        safeVal=0;
%         t_alarm=t
%         u_alarm=u
    else
        safeVal=1;
    end
end


end

