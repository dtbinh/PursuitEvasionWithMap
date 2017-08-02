function [count,allValidPaths,pathHist,visitedMat] =...
    countPaths3D(mazeMat,visitedMat,xyz,xyzMazeEnd,count,pathHist,allValidPaths)
%pathhist is 2xN matrix of current path history
%visitedMat is AxB matrix of visited points.  1=visited
%mazeMat is AxB matrix of full/empty coordinates. 0~empty, 1~blocked
%when called from main, pathHist=[], allValidPaths={}
%count is number of valid path
%NOTE: COORDINATES RUN FROM TOP LEFT OF MATRIX AND DENOTE ELEMENTWISE
%   POSITION IN MAP MATRIX.  [1,1] IS TOP LEFT CORNER, [A,B] IS BOTTOM
%   RIGHT CORNER

pathHist=[pathHist xyz];

if xyz==xyzMazeEnd %check if you've arrived.  If so, store and stop
    count=count+1;
    nn=length(allValidPaths);
    visitedMat(xyz(1),xyz(2))=1; %prepare to backtrack
    allValidPaths{nn+1}=pathHist;
else
    visitedMat(xyz(1),xyz(2),xyz(3))=1; %store location
    if mazeMat(xyz(1),xyz(2),xyz(3))==0 %current square is free
        %pathHist=[pathHist xyCurrent];
        
        %check to see if the current point is in the maze map
        if isValidCell3D(xyz,mazeMat)
            
            %break into 6 search directions
            %east
            xyzT=xyz+[1;0;0];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %move right, recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
            %west
            xyzT=xyz-[1;0;0];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
            %north
            xyzT=xyz+[0;1;0];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
            %south
            xyzT=xyz-[0;1;0];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
            %up
            xyzT=xyz+[0;0;1];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
            %down
            xyzT=xyz-[0;0;1];
            if isValidCell3D(xyzT,mazeMat)
                if visitedMat(xyzT(1),xyzT(2),xyzT(3))==0 %recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths3D(mazeMat,visitedMat,...
                        xyzT,xyzMazeEnd,count,pathHist,allValidPaths);
                end
            end
        end
    end
end

%backtrack http://www.techiedelight.com/find-total-number-unique-paths-maze-source-destination/
lPH=length(pathHist);
if lPH>=2
    pathHist=pathHist(:,1:lPH-1);
elseif lPH==1
    pathHist=[];
end
visitedMat(xyz(1),xyz(2),xyz(3))=0;

end

