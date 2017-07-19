function [x1,P1] = linearKFStep(x0,z,A_tr,B_tr,Gammak,P0,Q0,u,H,R0,xbarPRIOR)

nX=length(x0);
if nargin<=10
    x1bar=A_tr*x0+B_tr*u;
else
    x1bar=xbarPRIOR;  %just a simplification
end
P1bar=A_tr*P0*A_tr'+Gammak*Q0*Gammak';
Sk=H*P1bar*H'+R0;
Kk=P1bar*H'*inv(Sk);
x1=x1bar+Kk*(z-H*x1bar);
P1=(eye(nX)-Kk*H)*P1bar;


end

