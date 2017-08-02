function [count,allValidPaths,pathHist,visitedMat] =...
    countPaths(mazeMat,visitedMat,xyCurrent,xyMazeEnd,count,pathHist,allValidPaths)
%pathhist is 2xN matrix of current path history
%visitedMat is AxB matrix of visited points.  1=visited
%mazeMat is AxB matrix of full/empty coordinates. 0~empty, 1~blocked
%when called from main, pathHist=[], allValidPaths={}
%NOTE: COORDINATES RUN FROM TOP LEFT OF MATRIX AND DENOTE ELEMENTWISE
%   POSITION IN MAP MATRIX.  [1,1] IS TOP LEFT CORNER, [A,B] IS BOTTOM
%   RIGHT CORNER


pathHist=[pathHist xyCurrent];

if xyCurrent==xyMazeEnd %check if you've arrived.  If so, store and stop
    count=count+1;
    nn=length(allValidPaths);
    visitedMat(xyCurrent(1),xyCurrent(2))=1; %prepare to backtrack
    allValidPaths{nn+1}=pathHist;
else
    visitedMat(xyCurrent(1),xyCurrent(2))=1; %store location
    if mazeMat(xyCurrent(1),xyCurrent(2))==0 %current square is free
        %pathHist=[pathHist xyCurrent];
        %check to see if the current point is in the maze map
        if isValidCell(xyCurrent,mazeMat)
            %break into 4 directions to search.  (6 if 3D)
            if isValidCell(xyCurrent+[1;0],mazeMat)
                if visitedMat(xyCurrent(1)+1,xyCurrent(2))==0 %move right, recursive call
                    [count,allValidPaths,pathHist,visitedMat]=countPaths(mazeMat,visitedMat,...
                        xyCurrent+[1;0],xyMazeEnd,count,pathHist,allValidPaths);
                end
            end
            if isValidCell(xyCurrent-[1;0],mazeMat)
                if visitedMat(xyCurrent(1)-1,xyCurrent(2))==0
                    [count,allValidPaths,pathHist,visitedMat]=countPaths(mazeMat,visitedMat,...
                        xyCurrent-[1;0],xyMazeEnd,count,pathHist,allValidPaths);
                end
            end
            if isValidCell(xyCurrent+[0;1],mazeMat)
                if visitedMat(xyCurrent(1),xyCurrent(2)+1)==0
                    [count,allValidPaths,pathHist,visitedMat]=countPaths(mazeMat,visitedMat,...
                        xyCurrent+[0;1],xyMazeEnd,count,pathHist,allValidPaths);
                end
            end
            if isValidCell(xyCurrent-[0;1],mazeMat)
                if visitedMat(xyCurrent(1),xyCurrent(2)-1)==0
                    [count,allValidPaths,pathHist,visitedMat]=countPaths(mazeMat,visitedMat,...
                        xyCurrent-[0;1],xyMazeEnd,count,pathHist,allValidPaths);
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
visitedMat(xyCurrent(1),xyCurrent(2))=0;

end

