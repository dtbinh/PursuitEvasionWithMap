function [uPur,uEva] = generateCullSolveProcess(ttf,uPconst,uPtheta,uEconst,uEtheta,xPurEst,xEvaEst,vWW)
%generates uconMat, culls for collisions, solves the game, and then processes into controls
%designed for *polar* controls only
%set vWW=-1 for free space action

global cddGG nXGG dtGG umaxPurGG umaxEvaGG
nX=nXGG;
cdd=cddGG;
dt=dtGG;
umaxPur=umaxPurGG;
umaxEva=umaxEvaGG;

eyeHalfNX=eye(nX/2);
zerosHalfNX=zeros(nX/2,nX/2);

A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

uconPtemp=uConMat_Rt(ttf,uPconst,uPtheta);
if max(uconPtemp)>umaxPur
    fprintf('Max thrust exceeded\n')
end
uconEtemp=uConMat_Rt(ttf,uEconst,uEtheta);
if max(uconEtemp)>umaxEva
    fprintf('Max thrust exceeded\n')
end

[a1,a2,a3,~]=size(uconPtemp);
uconPcull=zeros(a1,a2,a3,2);
[a1,a2,a3,~]=size(uconEtemp);
uconEcull=zeros(a1,a2,a3,2);
nmod_pur_temp=length(uconPtemp); nmod_eva_temp=length(uconEtemp);
%remove collision points from map; technically, each party should
%create their own uconCull mats but whatever

mapEdges; %load wall locations into memory

if vWW==-1
    uconPcull=uconPtemp;
    uconEcull=uconEtemp;
    nmod_pur=nmod_pur_temp;
    nmod_eva=nmod_eva_temp;
else
    numSafeControls=0;
    for i1=1:nmod_pur_temp
        uhist=uconPtemp(:,:,:,i1);  %Should check for collisions from x0=left of quad and x0=right of quad
        isControlSafe=isSafeFromFutureCollision(ttf,xPurEst,A_tr,B_tr,uhist,wallPoints,numObj,vWW);
        if isControlSafe==1
            numSafeControls=numSafeControls+1;
            uconPcull(:,:,:,numSafeControls)=uconPtemp(:,:,:,i1);
        end
    end
    uconPcull=uconPcull(:,:,:,1:numSafeControls); %removes [0;0] controls added when initializing
    nmod_pur=length(uconPcull);
    numSafeControls=0;
    for i1=1:nmod_eva_temp
        uhist=uconEtemp(:,:,:,i1);
        isControlSafe=isSafeFromFutureCollision(ttf,xEvaEst,A_tr,B_tr,uhist,wallPoints,numObj,vWW);
        if isControlSafe
            %fprintf('foundsafeE\n')
            numSafeControls=numSafeControls+1;
            uconEcull(:,:,:,numSafeControls)=uconEtemp(:,:,:,i1);
        end
    end
    uconEcull=uconEcull(:,:,:,1:numSafeControls); %removes [0;0] controls added when initializing
    nmod_eva=length(uconEcull);
end

hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
hE=ttf*ones(nmod_eva,1);

strategiesE=struct('constant',uconEcull,'horizon',hE);
strategiesP=struct('constant',uconPcull,'horizon',hP);

[Jpur,Jeva]=generateCostMatrices(strategiesP,strategiesE,xPurEst,xEvaEst,0);
[eqLoc,nashReturnFlag,~]=findRDEq(-Jpur,-Jeva); %pass value matrices rather than cost matrices
if nashReturnFlag>=1 %if there IS an RDEq
    [uPur,uEva] = processNashType1(eqLoc,umaxPur,umaxEva,uconPcull,uconEcull);
elseif nashReturnFlag==0 %suboptimal
    fprintf('Running LH2\n')
    %split for efficiency
    [uPur,uEva] = processNashType0(-Jpur,-Jeva,umaxPur,umaxEva,uconPcull,uconEcull);
else
    uPur=zeros(nX/2,1);
    uEva=zeros(nX/2,1);
end




end

