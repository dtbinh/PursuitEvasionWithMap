function [nearestWallIndex, nearestPointPairIndex, nearestDistanceToWall] = nearestWallInformation(x,wallPoints,numObj)
%returns nearest wall index and distance to the wall
%nearestPointPair returns the *first* index of the point in the appropriate
%wall point pair matrix.  If the correct line is from wpt(:,end) to
%wpt(:,1), nPPI=END instead of 1 (to avoid confusion with 1-2 connection)

nearestWallIndex=0; nearestPointPairIndex=0; nearestDistanceToWall=9001;

for j2=1:numObj
    wpt=wallPoints{j2};
    nl=length(wpt);
    for j3=1:nl
        p1=wpt(:,j3);
        if j3<nl; p2=wpt(:,j3+1); else; p2=wpt(:,1);end
        
        %get distance to line, see wikipedia article on distance to
        %a line from a point
        thisDist=abs( (p2(2)-p1(2))*x(1)-(p2(1)-p1(1))*x(2)+p2(1)*p1(2)-p2(2)*p1(1))/sqrt((p2(2)-p1(2))^2+(p2(1)-p1(1))^2);
        
        %if distance is shortest so far, store it
        if thisDist<nearestDistanceToWall
            nearestDistanceToWall=thisDist;
            nearestWallIndex=j2;
            nearestPointPairIndex=j3;
        end
    end
end

end

