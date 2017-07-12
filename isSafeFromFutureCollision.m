function [isSafe,nWI,nPPI,nWD] = isSafeFromFutureCollision(ttf,xp1,A_tr,B_tr,...
            uhist,wallPoints,numObj,virtualWallWidth)
%Only consider strategies that stay outside of virtualWallWidth of a
%physical wall.  If xp1 is already within vWW of a wall, then this will
%only present strategies that go through it.
%NOTE:  Only handles a *single* wall, if within vWW of two walls, this will
%move away from the nearest while trying to avoid running into the other.
%  Set virtualWallWidth==0 to ignore this action.

%The return argument isSafe is thus:
%(a) if it is NOT already behind a virtual wall, then isSafe=1 IFF it does not
%break through a virtual wall
%(b) if it is already behidn a virtual wall, then isSafe=1 IFF if it moves
%through the virtual wall to safety and does NOT break through another
%virtual wall

if nargin<=7
    virtualWallWidth=0;
end
nWI=0; nPPI=0; nWD=9001;
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
    
    %see if it is already between a virtual wall and a real wall
    [nWI_temp,nPPI_temp,nWD_temp]=nearestWallInformation(xp1,wallPoints,numObj);
    
    %output nearest wall to initial point
    nWI=nWI_temp;
    nWD=nWD_temp;
    nPPI=nPPI_temp;
    
    if nWD_temp<virtualWallWidth
        
        %DO STUFF HERE
                %same loop as before but only consider strategies that "collide" with
        %shifted virtual wall
        countPassagesThroughVW=0;
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
                    
                    
                    %DO STUFF HERE
                    p1shift=p1;
                    p2shift=p2;
                    
                    
                    
                    %see if chosen strategy breaks through wall
                    isSafeThis_temp=checkSafeIntersection(xp1(1:nX/2),xp2Shift(1:nX/2),p1shift,p2shift);
                    
                    %if correct strategy DOES break through virtual wall
                    if i2==nWI_temp && i3==nPPI_temp
                        if isSafeThis_temp==0 %see if we've broken through the nearest VW
                            countPassagesThroughVW=countPassagesThroughVW+1;
                        end %we DO want to pass through VW that we are behind
                        %
                        %now, check to ensure that the "real" wall that you
                        %are sandwiched between is not broken
                        isSafeThis=checkSafeIntersection(xp1(1:nX/2),xp2Shift(1:nX/2),p1,p2);
                        %
                        %With this order of calls, we only count how many
                        %times the nearest VW has been broken (odd number
                        %means that, in the end, it escapes).  To continue
                        %the loop (corresp isSafeThis), we ONLY need to
                        %guarantee that we have not broken through the real
                        %wall.
                        %
                    else %for virtual walls that we do NOT want to break through
                        isSafeThis=isSafeThis_temp; 
                    end
                    
                    if isSafeThis==0
                        contloop=false;
                        isSafe=0;
                    end
                end
            end
            xp1=xp2;
        end
        if mod(countPassagesThroughVW,2)==1 && isSafe==1
            isSafe=1; %if it passes through the VW an odd number of times AND
                      %it is safe from other collisions, then this control
                      %is safe.
        else
            isSafe=0;
        end
        
    else %the goal is to avoid future collisions
        %find nearest wall on trajectory
        xpEnd=xp1;
        continueFindWall=1;
        k=0;
        while continueFindWall==1
        	k=k+1;
            xpEnd=A_tr*xpEnd+B_tr*uhist(:,k);
            [nWI_temp,nPPI_temp,nWD_temp]=nearestWallInformation(xpEnd,wallPoints,numObj);
            if nWD_temp<nWD
                nWI=nWI_temp;
                nPPI=nPPI_temp;
                nWD=nWD_temp;
            end
            if k>=ttf
                continueFindWall=0;
            end
            if nWD<=virtualWallWidth
                isSafe=0;
                continueFindWall=0;
            end
        end
        
    end
    
    
%     if nWD<=virtualWallWidth  %only move THROUGH virtual wall if player is too close to it
%     else
%         %call recursively without using the virtual wall
%         isSafe=isSafeFromFutureCollision(ttf,xp1,A_tr,B_tr,uhist,wallPoints,numObj,0); 
%     end
end

end

