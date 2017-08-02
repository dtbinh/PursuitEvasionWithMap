clear;clc;

% xvajPpurEst; xvajEpurEst;
% xvajPevaEst; xvajEevaEst;
xvajPpurEst=zeros(12,1); xvajPpurEst(1:3)=[-3;3.5;1.5];
xvajEpurEst=zeros(12,1); xvajEpurEst(1:3)=[.5;3.5;1.5];

tf=10;

mapFloor=[1 1 1 1
          1 1 1 1
          1 0 0 0
          1 1 1 1];
mapGame=[0 0 0 0
         1 0 1 0
         1 0 1 0
         1 0 0 0];
mapCeil=[1 1 1 1
         1 1 1 1
         1 1 1 1
         1 1 1 1];
map(:,:,1)=mapFloor;
map(:,:,2)=mapGame;
mapStart=[1;1;2]; mapEnd=[4;4;2];

[nmod,allValidPaths]=countPaths3D(map,zeros(size(map)),mapStart,mapEnd,0,[],{});
pointlist=convertNodeLocationsToXYZ(allValidPaths,map);
%solve game for pursuer
[nashPpur,nashEpur,eqLoc]=solveMazeGame(tf,xvajPpurEst,xvajEpurEst,pointlist,pointlist);
%postprocess into waypoints lists AND expected weighting of each set
utemp=genControlsFromNashmatAndWaypoints(nashPpur,xvajPpurEst(1:3),pointlist);
uPurTrueProb=utemp.probability; uPurTrueTraj=utemp.waypoints;
utemp=genControlsFromNashmatAndWaypoints(nashPpur,xvajEpurEst(1:3),pointlist);
uEvaExpectedProb=utemp.probability; uEvaExpectedTraj=utemp.waypoints;
%choose a set of waypoints (for case with max(nashP)!=1)
uPurTrueIndex=randsample(length(nashPpur),1,true,nashPpur);
uPurTrueWaypoints=uPurTrueTraj{uPurTrueIndex};
uPurTruePossible=find(nashPpur>0);
uEvaExpectedPossible=find(nashEpur>0);

%[uP_expected,uE_true]=solveMazeGame(xvajPpurEst,xvajEpurEst,pointlist,pointlist);




