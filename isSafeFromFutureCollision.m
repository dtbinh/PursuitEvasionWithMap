function [isSafe,nWI,nPPI,nWD] = isSafeFromFutureCollision(ttf,xp1,A_tr,B_tr,...
            uhist,wallPoints,numObj,virtualWallWidth)
%if xp1 is within vWW of a wall, only consider strategies that move away.
%  Set==0 to ignore this action.
if nargin<=7
    virtualWallWidth=0;
end
nWI=0; nPPI=0;
nX=length(xp1);
contloop=true;
ncount=0;
isSafe=1;  %returns 1 for safe
if virtualWallWidth==0
    while contloop %break loop immediately if collision flag seen
        tvec=[]; uvec=[];
        ncount=ncount+1;
        if ncount>=ttf||ncount>=20; contloop=false; end
        xp2=A_tr*xp1+B_tr*uhist(:,:,ncount);
        xp2Shift=xp2*1; %work out better safe shifting metric
        for i2=1:numObj %iterate through objects
            wpt=wallPoints{i2};
            nl=length(wpt);
            
            for i3=1:nl %iterate through connected point pairs
                p1=wpt(:,i3);
                
                %more efficient wraparound when i3->endpt
                if i3<nl; p2=wpt(:,i3+1); else; p2=wpt(:,1);end
                
                [isSafeThis,ttt,uuu]=checkSafeIntersection(xp1(1:nX/2),xp2Shift(1:nX/2),p1,p2); %#ok<ASGLU>
                if isSafeThis==0
                    contloop=false;
                    isSafe=0;
                end
            end
        end
        xp1=xp2;
    end
else
    %find nearest wall somehow
    for j2=1:numObj
        wpt=wallPoints{j2};
        nl=length(wpt);
        nearestWallIndex=0;  %contains the index number into wallPoints j2
        nearestPointPairIndex=0; %nPPI is j3 where the nearest wall is
        %along j3 to j3+1. If j3==len(wpt), the pair is (j3,1)
        nearestWallDist=9001;
        
        xpEnd=(A_tr)^ttf*xp1; %Consider where the endpoint would be if the
        %current trajectory was maintained without
        %control.  Addressing coasting behavior is
        %the most important thing to do.
        
        for j3=1:nl
            p1=wpt(:,j3);
            if j3<nl; p2=wpt(:,j3+1); else; p2=wpt(:,1);end
            
            %get distance to line, see wikipedia article on distance to
            %a line from a point
            thisDist=abs( (p2(2)-p1(2))*xpEnd(1)-(p2(1)-p1(1))*xpEnd(2)+p2(1)*p1(2)-p2(2)*p1(1))...
                /sqrt((p2(2)-p1(2))^2+(p2(1)-p1(1))^2);
            
            %if distance is decreasing, store this
            if thisDist<nearestWallDist
                nearestWallDist=thisDist;
                nearestWallIndex=j2;
                nearestPointPairIndex=j3;
            end
        end
    end
    nWI=nearestWallIndex; nPPI=nearestPointPairIndex; nWD=nearestWallDist;
    
    
    if nWD<=virtualWallWidth  %only move THROUGH virtual wall if player is too close to it
        %same loop as before but only consider strategies that "collide" with
        %shifted virtual wall
        doesPassThroughVW=1;
        while contloop %break loop immediately if collision flag seen OR if the chosen strategy doesn't go through the wall
            ncount=ncount+1;
            if ncount>=ttf||ncount>=20; contloop=false; end
            xp2=A_tr*xp1+B_tr*uhist(:,:,ncount);
            xp2Shift=xp2*1; %work out better safe shifting metric
            for i2=1:numObj
                wpt=wallPoints{i2};
                nl=length(wpt);
                for i3=1:nl
                    p1=wpt(:,i3);
                    if i3<nl; p2=wpt(:,i3+1); else; p2=wpt(:,1);end %more efficient wraparound
                    
                    if i2==nearestWallIndex && i3==nearestPointPairIndex %shift to nearest wall
                        
                        
                        %DO STUFF HERE
                        
                        
                        p1shift=p1;
                        p2shift=p2;
                    else
                        p1shift=p1;
                        p2shift=p2;
                    end
                    %see if chosen strategy breaks through wall
                    isSafeThis_temp=checkSafeIntersection(xp1(1:nX/2),xp2Shift(1:nX/2),p1shift,p2shift);
                    
                    %if correct strategy DOES break through virtual wall
                    if ncount==ttf && i2==nearestWallIndex && i3==nearestPointPairIndex
                        if isSafeThis_temp==1;doesPassThroughVW=0;else;doesPassThroughVW=1;end; %we DO want to pass through VW
                    else
                        isSafeThis=isSafeThis_temp;
                    end
                    
                    if isSafeThis==0 || doesPassThroughVW==0
                        contloop=false;
                        isSafe=0;
                    end
                end
            end
            xp1=xp2;
        end
    else
        %call recursively without using the virtual wall
        isSafe=isSafeFromFutureCollision(ttf,xp1,A_tr,B_tr,uhist,wallPoints,numObj,0); 
    end
end

end

