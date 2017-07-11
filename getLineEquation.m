function [a,b,c] = getLineEquation(p1,p2)
%expressed in aX+bY+c=0 form
if p1(1)~=p2(1)
    m=(p2(2)-p1(2))/(p2(1)-p1(1));
    a=m; b=1;
    c=p1(2)-m*p1(1);
else
    a=0;b=1;c=p1(2); %error condition
end

end

