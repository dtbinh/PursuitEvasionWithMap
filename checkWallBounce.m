function [xF,bounceFlag] = checkWallBounce(x0,u_control,cdd,dt,nPvec,wallPoints,numObj)
bounceFlag=0;  %bounceFlag=1 if it hits the wall

%TO DO: HANDLE CASE OF SLIDING OFF OF WALL

nX=length(x0);
eyeHalfNX=eye(nX/2); zerosHalfNX=zeros(nX/2,nX/2);
A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

%propagate expected dynamics
xF_withoutCollision=A_tr*x0+B_tr*u_control+Gammak*nPvec;

if isSafeFromFutureCollision(1,x0,A_tr,B_tr,u_control,wallPoints,numObj)==1
    xF=xF_withoutCollision;
    bounceFlag=0;
else
    %(1) find point/time of collision, get V there
    %(2) set V|_ =0, V|| unchanged
    %(3) see if it hits another wall
    %(4) see if it slides off
    contOuterLoop=1; iter=0; thisWallInd=9001; thisPPInd=9001;
    tt=0; %internal time var
    xS=x0;
    ucon=u_control;
    while contOuterLoop==1
        uMin=9001;
        iter=iter+1;
        for i2=1:length(numObj)
            wpt=wallPoints{i2};
            nl=length(wpt);
            for i3=1:nl
                if i3<nl; p2=wpt(:,i3+1); else; p2=wpt(:,1);end
                [~,~,u_control]=checkSafeIntersection(x0,xF_withoutCollision(1:2),wpt(:,i3),p2);
                if u_control>0 && u_control<1 && u_control<uMin
                    thisWallInd=i2; thisPPInd=i3; uMin=u_control;
                end
            end
        end
        dumvar=wallPoints{thisWallInd};
        nl=length(dumvar);
        p1=dumvar(:,thisPPInd);
        if thisPPInd<nl; p2=dumvar(:,thisPPInd+1); else; p2=dumvar(:,1);end
        wallVec=unit_vector(p2-p1);
        
        %iterate to find time to next collision
        dttemp=0;
        contInnerLoop=1;
        while contInnerLoop==1
            dttemp=dttemp+dt*.01; %.01% increments
            A_temp=[eyeHalfNX (1*dttemp-cdd*dttemp^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dttemp)*eyeHalfNX];
            B_temp=[dttemp^2/2*eyeHalfNX; dttemp*eyeHalfNX];
            Gamma_temp=[dttemp^2/2*eyeHalfNX; dttemp*eyeHalfNX];
            xFtemp=A_temp*xS+B_temp*ucon+Gamma_temp*nP;
            doesNotCollide=checkSafeIntersection(x0,xFtemp,p1,p2);
            
            %iterate until collision, check time
            if tt+dttemp>dt || doesNotCollide==0
                contInnerLoop=0;
            end
            
        end
        
        tr=dttemp-dt*.01; %time remaining to next collision
        A_temp=[eyeHalfNX (1*tr-cdd*tr^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*tr)*eyeHalfNX];
        B_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
        Gamma_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
        
        xS=A_temp*xS+B_temp*ucon+Gamma_temp*nP;
        
        ucon=dot(u_control,wallVec)*sign(dot(ucontrol,wallVec))*wallVec;
        
        if iter>=50 || tt>=dt
            contOuterLoop=0;
        end
    end
    
end


end

