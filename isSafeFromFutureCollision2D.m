function isSafe = isSafeFromFutureCollision2D(ttf,xp1,A_tr,B_tr,uhist,wallPoints,numObj)
%old generic versoin, only handles 2D collisions.  Models objects as
%straight-line point pairs
contloop=true;
ncount=0;
isSafe=true;
while contloop %break loop immediately if collision flag seen
    ncount=ncount+1;
    if ncount>=ttf||ncount>=20; contloop=false; end
    xp2=A_tr*xp1+B_tr*uhist(:,:,ncount);
    xp2Shift=xp2*1; %work out better safe shifting metric
    for i2=1:numObj
        wpt=wallPoints{i2};
        nl=length(wpt);
        for i3=1:nl
            p1=wpt(:,i3);
            if i3<nl; p2=wpt(:,i3+1); else; p2=wpt(:,1);end
            isSafeThis=checkSafeIntersection(xp1(1:2),xp2Shift(1:2),p1,p2);
            if ~isSafeThis
                contloop=false;
                isSafe=false;
            end
        end
    end
    xp1=xp2;
end

end

