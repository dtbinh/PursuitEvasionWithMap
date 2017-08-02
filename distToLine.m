function dist = distToLine(x_in,p1,p2)
%modified from javascript code, see https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

x=x_in(1); y=x_in(2); x1=p1(1); y1=p1(2); x2=p2(1); y2=p2(2);
A=x-x1; B=y-y1; C=x2-x1; D=y2-y1;

dd=A*C+B*D; len_sq=C^2+D^2;
param=-1;
if len_sq>1e-6 %length nonzero
    param=dd/len_sq;
end

if param<0
    xx=x1; yy=y1;
elseif param>1
    xx=x2;yy=y2;
else
    xx=x1+param*C;
    yy=y1+param*D;
end
dx=x-xx;
dy=y-yy;
dist=sqrt(dx^2+dy^2);


% %Distance to FULL line, not line segment
%dist=abs( (p2(2)-p1(2))*x(1)-(p2(1)-p1(1))*x(2)+p2(1)*p1(2)-p2(2)*p1(1))/sqrt((p2(2)-p1(2))^2+(p2(1)-p1(1))^2);

end

