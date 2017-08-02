function outstruct = genControlsFromNashmatAndWaypoints(nashmat,x0,pointlist)
%pointlist is the list of maze node points
%x0 is initial pose
cc=0;
for i=1:length(pointlist)
    wpt=pointlist{i};
    for jk=1:length(wpt)
        cc=cc+1;
        outstruct.probability{cc}=nashmat(cc);
        outstruct.waypoints{cc}=[x0 wpt(:,1:jk)];
    end
end


end

