function [xF,bounceFlag] = checkWallBounceAndPropagateDynamics(x0,u_control,cdd,dt,nPvec,wallPoints,numObj)
bounceFlag=0;  %bounceFlag=1 if it hits the wall

%TO DO: HANDLE CASE OF SLIDING OFF OF WALL WHEN WALL ENDS

nX=length(x0);
eyeHalfNX=eye(nX/2); zerosHalfNX=zeros(nX/2,nX/2);
A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

%propagate expected dynamics
xF_withoutCollision=A_tr*x0+B_tr*u_control+Gammak*nPvec;
%xfstoch=A_tr*x0+B_tr*u_control;

if isSafeFromFutureCollision(1,x0,A_tr,B_tr,u_control,wallPoints,numObj,0)==1
    xF=xF_withoutCollision;
    bounceFlag=0;
else
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
                [~,~,udum]=checkSafeIntersection(x0(1:nX/2),xF_withoutCollision(1:nX/2),wpt(:,i3),p2);
                if udum>0 && udum<1 && udum<uMin
                    thisWallInd=i2; thisPPInd=i3; uMin=udum;
                end
            end
        end
        
        if uMin==9001 %no collision found, propagate to end of dt
            tr=dt-tt;
            A_temp=[eyeHalfNX (1*tr-cdd*tr^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*tr)*eyeHalfNX];
            B_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
            Gamma_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
            xS=A_temp*xS+B_temp*ucon+Gamma_temp*nPvec;
        else
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
                xFtemp=A_temp*xS+B_temp*ucon+Gamma_temp*nPvec;
                doesNotCollide=checkSafeIntersection(x0(1:nX/2),xFtemp(1:nX/2),p1,p2);
                
                %iterate until collision, check time
                if tt+dttemp>dt || doesNotCollide==0
                    contInnerLoop=0;
                end
                
            end
            
            tr=dttemp-dt*.01; %time remaining to next collision
            A_temp=[eyeHalfNX (1*tr-cdd*tr^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*tr)*eyeHalfNX];
            B_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
            Gamma_temp=[tr^2/2*eyeHalfNX; tr*eyeHalfNX];
            xS=A_temp*xS+B_temp*ucon+Gamma_temp*nPvec;
            
            ucon=dot(u_control,wallVec)*sign(dot(u_control,wallVec))*wallVec;
            xS(nX/2+1:nX)=dot(xS(nX/2+1:nX),wallVec)*sign(dot(xS(nX/2+1:nX),wallVec))*wallVec
        end
        
        if iter>=50 || tt>=dt
            contOuterLoop=0;
        end
    end
    
    xF=xS;
    
end


end

