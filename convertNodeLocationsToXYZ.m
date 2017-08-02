function [xyzlist] = convertNodeLocationsToXYZ(nodelist,maze, nodeToMazeDists)
%moves list of array locations to 0,0,0 at origin
%nodeToMazeDists are the distances corresponding to size of X/Y/Z blocks in the maze
%the returned waypoints are at the centers of blocks

if nargin<=2 %use 
    nodeToMazeDists=[1;1;1];
end

xyzlist={};
[a,~,~]=size(maze);
for ik=1:length(nodelist)
    thisPoints=nodelist{ik};
    thisPointsFixed=zeros(3,length(thisPoints));
    thisPointsFixed(1,:)=nodeToMazeDists(1)*(thisPoints(2,:)-1)+nodeToMazeDists(1)/2; %x~columns
    thisPointsFixed(2,:)=nodeToMazeDists(2)*(a-thisPoints(1,:))+nodeToMazeDists(2)/2;
    thisPointsFixed(3,:)=nodeToMazeDists(3)*(thisPoints(3,:)-1)+nodeToMazeDists(3)/2;
    xyzlist{ik}=thisPointsFixed;
end




end

