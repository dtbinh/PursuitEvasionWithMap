clear;clc;
%closes the loop for control switching in cluttered maps

%NOTES
% work out better metric for safe shifting to account for covariance matrix
% shift wall by virtual wall amount in isSafeFromFutureCollision
% vWW handling would be better if it provided a "reward" for moving away
%   from the wall rather than culling those that don't move away
%
%BUGS
%If a point is exactly vWW away from a wall, then all controls are unsafe
%vWW is still rejecting all controls if too close

tmax=25;
lookAheadTime=1; %number of steps to look ahead in optimization; set 0 to use t_remaining instead
dt=1;
MCmax=1;
nSim=ceil(tmax/dt); %simulation time within MC
fineificationSteps=1; %no refinement if ==0
scaleRed=8; %each fineification step has scaleRed+1 elements in each state. DO NOT USE <2
%nSim=1;
virtualWallWidth=0; %if within vWW of wall, only consider strategies that
                     %move away from the wall.  Set ==0 to ignore.  Set ==-1
                     %for free space

flagDiminishHorizonNearNSim=1;  %if ==1, consider lookAheadTime
         %if TRemain>lookAheadTime, use TRemain if TRemain<lookAheadTime
flagUseDetm=1;  %zero noise if flag==1
flagUseFeedback=0; %generate control history via feedback
flagUseModeControlInsteadOfMean=1;

muAllP=zeros(MCmax,nSim,2);
muAllE=zeros(MCmax,nSim,2);
nsimTrunc=zeros(MCmax,1);
captureIndex=zeros(MCmax,1);
flagPlotAsGo=0;  %=1 to plot motion
flagPlotFinal=0; %=1 to plot final dist/speed history
flagPlotArena=1; %=1 to plot arena edges
plotint=.1;

nX=4;   %full size of state space (even number)
nU=2;
nV=nX/2;
nZ=2*nX;

eStore=zeros(nX,nSim,MCmax);

for MCL=1:MCmax
        
    uphist=[];
    uehist=[];
    
    cdd=.25; %drag coefficient, =0 to ignore
    eyeNX=eye(nX); eyeHalfNX=eye(nX/2);
    zerosNX=zeros(nX,nX); zerosHalfNX=zeros(nX/2,nX/2);
    A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
    B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    thetaPpur=0; thetaEeva=0;
    Hstack=eye(nZ); %each player measures each player's positions
    P0=.1*eyeNX;
    P0stack=[P0 zerosNX; zerosNX P0];
    P_pur=P0stack; P_eva=P0stack;
    Q0p=.00001*eye(nV); %process noise, pursuer
    Q0e=.00001*eye(nV);
    Q0stack=[Q0p zerosHalfNX; zerosHalfNX Q0e];
    R0P=.00005*eye(nZ); %measurement noise
    R0E=R0P;
    cholQ0p_T=chol(Q0p)';
    cholQ0e_T=chol(Q0e)';
    cholR0P_T=chol(R0P)';
    cholR0E_T=chol(R0E)';
    
    %cost matrices
    QfPur=8*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepPur=zeros(nX,nX);
    RPur=25*eye(nU);  %NOTE: discrete R = continuous R*dt^2
    RPurScaled=RPur*dt^2;
    
    QfEva=5*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepEva=zeros(nX,nX);
    REva=.01*eye(nU);
    REvaScaled=REva*dt^2;
    global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG QnoisePGG QnoiseEGG...
        RPurGG umaxPurGG umaxEvaGG nXGG QnoiseGG dtGG cddGG nUGG%#ok
    QfEvaGG=QfEva; QfPurGG=QfPur; QstepEvaGG=QstepEva; QnoisePGG=Q0p; cddGG=cdd; QnoiseEGG=Q0e;
    QstepPurGG=QstepPur; REvaGG=REvaScaled; RPurGG=RPurScaled; nXGG=nX; dtGG=dt; nUGG=nU;
    
    
    nmod_pur=2; %number of models considered
    MijE = 1/nmod_pur*ones(nmod_pur,nmod_pur); %probability of mode switching
    
    nmod_eva=2; %number of models considered
    MijP = 1/nmod_pur*ones(nmod_eva,nmod_eva); %probability of mode switching

    flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
    flagUseXPrevInsteadOfXModel=0;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.    
    
    mu_min=10^-2;
    normpdf_diag_DEBUG=.01; %with high confidence, diagonal elements of normpdf matrix go to zero
    
    
    xEva=[1;.3;0;0];
    xPur=[.11;.45;.1;0];
    
    umaxPur=.1;
    umaxEva=.1;
    umaxPurGG=umaxPur; umaxEvaGG=umaxEva;
    JpurRunning=0;  %FOR COMPARISON WITH KUMARV2 ONLY
    JevaRunning=0;
    
    captureThresh=.5;
    flagSign1Flip=0;
    flagSign2Flip=0;
    
    %set randn's seed RNG
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    ehat0=xEva-xPur+(chol(P0))'*randn(nX,1);
    if flagUseDetm==1
        ehat0=xEva-xPur; %deterministic behavior
    end
    ehat_prev_pur=ehat0;
    
    Pstore_pur=zeros(nX,nX,nmod_pur,nSim);
    Pstore_pur(:,:,1,1)=P0;
    Pstore_pur(:,:,2,1)=P0;
    xhathist_pur=zeros(nX,nmod_pur,nSim);
    
    mukhist_pur=[];
    mu0=1/nmod_pur*ones(nmod_pur,1);
    LambdaVec_pur=zeros(nmod_pur,nSim);
    mukhist_pur=[mukhist_pur mu0];
    muPrev_pur=mu0;
    xhathist_weighted_pur=ehat0;
    
    if flagPlotArena==1
        figure(1); clf;
        mapEdges;
        for i1=1:numObj
            wpt=wallPoints{i1};
            nn=length(wpt);
            for i2=1:nn-1
                hold on
                plot([wpt(1,i2) wpt(1,i2+1)],[wpt(2,i2) wpt(2,i2+1)],'k');
            end
            hold on
            plot([wpt(1,1) wpt(1,nn)],[wpt(2,1) wpt(2,nn)],'k');
        end
        purPlot=plot(xPur(1),xPur(2),'g*');
        evaPlot=plot(xEva(1),xEva(2),'r*');
        axis(axisSet);
    end

    
    Pstore_eva=zeros(nX,nX,nmod_eva,nSim);
    Pstore_eva(:,:,1,1)=P0;
    Pstore_eva(:,:,2,1)=P0;
    xhathist_eva=zeros(nX,nmod_eva,nSim);
    
    mukhist_eva=[];
    mu0=1/nmod_eva*ones(nmod_eva,1);
    LambdaVec_eva=zeros(nmod_eva,nSim);
    mukhist_eva=[mukhist_eva mu0];
    muPrev_eva=mu0;
    xhathist_weighted_eva=ehat0;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    global umax_GG %#ok<TLEV,NUSED>
    
    e_true=zeros(nX,nSim);
    
    LambdaTemp=zeros(nmod_pur,1);
    
    optionsForFMC=optimset('Display','off');
    numSimRuns=0;
    i=1;
    
    if flagPlotAsGo==1
        figure(42);clf;
        gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red');
    end
    
    for i=1:nSim
        numSimRuns=numSimRuns+1;
        i %#ok<NOPTS>
        
        if i==1
            xPpur=xPur;
            xEpur=xEva;
            xPeva=xPur;
            xEeva=xEva;
        end
        x0p=[xPpur;xEpur]; x0e=[xPeva;xEeva]; %a priori estimates for xP,xE for EKF
        
        %generate greedy open-loop controls
        if flagUseFeedback ~= 1
            if lookAheadTime==0
                ttf=ceil(tmax/dt + 1 - i);
            elseif flagDiminishHorizonNearNSim==1
                tRemain=nSim-i+1;
                if tRemain>=lookAheadTime
                    ttf=lookAheadTime;
                else
                    ttf=tRemain;
                end
            else
                ttf=lookAheadTime;
            end
%             
%             uconPtemp=uConMat_Rt(ttf,uPconst,uPtheta);
%             uconEtemp=uConMat_Rt(ttf,uEconst,uEtheta);
        end
        
        uPthetaskipP=30*pi/180; uEthetaskipP=uPthetaskipP;
        uPconstskipP=umaxPur/5; uEconstskipP=umaxEva/5;
        uPthetaskipE=uPthetaskipP; uEthetaskipE=uEthetaskipP;
        uPconstskipE=uPconstskipP; uEconstskipE=uEconstskipP;
        uPthetaP = 0:uPthetaskipP:2*pi-uPthetaskipP; %.177
        uPconstP = .01:uPconstskipP:umaxPur;
        uEthetaP = 0:uEthetaskipP:2*pi-uEthetaskipP; %.097
        uEconstP = .01:uEconstskipP:umaxEva;
        uPthetaE = 0:uPthetaskipE:2*pi-uPthetaskipP; %.177
        uPconstE = .01:uPconstskipE:umaxPur;
        uEthetaE = 0:uEthetaskipE:2*pi-uEthetaskipP; %.097
        uEconstE = .01:uEconstskipE:umaxEva;
        
        for jk=1:fineificationSteps+1
            %get NE controls
            [uPurTrue,uEvaExpected]=generateCullSolveProcess(ttf,uPconstP,uPthetaP,uEconstP,uEthetaP,xPpur,xEpur,virtualWallWidth);
            [uPurExpected,uEvaTrue]=generateCullSolveProcess(ttf,uPconstE,uPthetaE,uEconstE,uEthetaE,xPeva,xEeva,virtualWallWidth);
            
            %The next block finds the local NE neighborhood, set up finer discretization (xAyB~x of P from perspective of B)
            %
            %extracting true parameters from optimal controls
            uPconstCenterP=norm(uPurTrue); uEconstCenterP=norm(uEvaExpected);
            uEconstCenterE=norm(uEvaTrue); uPconstCenterE=norm(uPurExpected);
            uPthetaCenterP=atan2(uPurTrue(2),uPurTrue(1)); uEthetaCenterP=atan2(uEvaExpected(2),uEvaExpected(1));
            uEthetaCenterE=atan2(uEvaTrue(2),uEvaTrue(1)); uPthetaCenterE=atan2(uPurExpected(2),uPurExpected(1));
            %uconP,uconE from pursuer perspective.  Mesh is refined by
            %scaleRed but the reduced scale is not stored
            uPconstP=max(0,uPconstCenterP-uPconstskipP):uPconstskipP*2/scaleRed:min(uPconstCenterP+uPconstskipP,umaxPur);
            uEconstP=max(0,uEconstCenterP-uEconstskipP):uEconstskipP*2/scaleRed:min(uEconstCenterP+uEconstskipP,umaxEva);
            uPthetaP=uPthetaCenterP-uPthetaskipP : uPthetaskipP*2/scaleRed : uPthetaCenterP+uPthetaskipP;
            uEthetaP=uEthetaCenterP-uEthetaskipP : uEthetaskipP*2/scaleRed : uEthetaCenterP+uEthetaskipP;
            %uconP, uconE from evader perspective
            uPconstE=max(0,uPconstCenterE-uPconstskipE):uPconstskipE*2/scaleRed:min(uPconstCenterE+uPconstskipE,umaxPur);
            uEconstE=max(0,uEconstCenterE-uEconstskipE):uEconstskipE*2/scaleRed:min(uEconstCenterE+uEconstskipE,umaxEva);
            uPthetaE=uPthetaCenterE-uPthetaskipE : uPthetaskipE*2/scaleRed : uPthetaCenterE+uPthetaskipE;
            uEthetaE=uEthetaCenterE-uEthetaskipE : uEthetaskipE*2/scaleRed : uEthetaCenterE+uEthetaskipE;
            %stores new mesh size
            uPconstskipP=uPconstskipP*2/scaleRed; uPthetaskipP=uPthetaskipP*2/scaleRed;
            uEconstskipP=uEconstskipP*2/scaleRed; uEthetaskipP=uEthetaskipP*2/scaleRed;
            uPconstskipE=uPconstskipE*2/scaleRed; uPthetaskipE=uPthetaskipE*2/scaleRed;
            uEconstskipE=uEconstskipE*2/scaleRed; uEthetaskipE=uEthetaskipE*2/scaleRed;
        end
        
        uphist=[uphist uPurTrue]; %#ok
        uehist=[uehist uEvaTrue]; %#ok
        
        %process noise
        if flagUseDetm==1 %deterministic behavior
            nP=zeros(nV,1);
            nE=zeros(nV,1);
        else
            nP=cholQ0p_T*randn(nV,1);
            nE=cholQ0e_T*randn(nV,1);
        end
        %Propagate states
%        xPur=A_tr*xPur+B_tr*uPurTrue+Gammak*nP;
%        xEva=A_tr*xEva+B_tr*uEvaTrue+Gammak*nE;
%         
%         %Does it bounce off of a wall?
        wallbounceP=0; wallbounceE=0;
        [xPur,wallbounceP]=checkWallBounceAndPropagateDynamics(xPur,uPurTrue,cdd,dt,nP,wallPoints,numObj);
        [xEva,wallbounceE]=checkWallBounceAndPropagateDynamics(xEva,uEvaTrue,cdd,dt,nE,wallPoints,numObj);
        
        
        if flagPlotArena==1
            pause(.1);
            delete(purPlot); delete(evaPlot);
            mapEdges;
            for i1=1:numObj
                wpt=wallPoints{i1};
                nn=length(wpt);
                for i2=1:nn-1
                    hold on
                    plot([wpt(1,i2) wpt(1,i2+1)],[wpt(2,i2) wpt(2,i2+1)],'k');
                end
                hold on
                plot([wpt(1,1) wpt(1,nn)],[wpt(2,1) wpt(2,nn)],'k');
            end
            purPlot=plot(xPur(1),xPur(2),'g*');
            evaPlot=plot(xEva(1),xEva(2),'r*');
            axis(axisSet);
            legend(strcat('t=',num2str(i)));
        end
        
        %Measurement
        if flagUseDetm==1
            zE=randn(nZ,1);
            zP=randn(nZ,1);
        else
            zE=zeros(nZ,1);
            zP=randn(nZ,1);
        end
        zEva=Hstack*[xPur;xEva]+cholR0E_T*zE;
        zPur=Hstack*[xPur;xEva]+cholR0E_T*zP;
        Astack=[A_tr zerosNX; zerosNX A_tr]; Bstack=[B_tr zeros(4,2); zeros(4,2) B_tr];
        uCombPur=[uPurTrue;uEvaExpected];
        uCombEva=[uPurExpected;uEvaTrue];
        GammaLinStackP=[Gammak zeros(4,2); zeros(4,2) Gammak];
        GammaLinStackE=[Gammak zeros(4,2); zeros(4,2) Gammak];
        if wallbounceE==0 && wallbounceP==0 %use standard KF if none hit the wall
            [xstackP,P_pur]=linearKFStep(x0p,zPur,Astack,Bstack,GammaLinStackP,P_pur,Q0stack,uCombPur,Hstack,R0P);
            [xstackE,P_eva]=linearKFStep(x0e,zEva,Astack,Bstack,GammaLinStackE,P_eva,Q0stack,uCombEva,Hstack,R0E);
        else  %use modified KF (with prior) to handle the bounce off of the wall
            [xstackP,P_pur]=linearKFStep(x0p,zPur,Astack,Bstack,GammaLinStackP,P_pur,Q0stack,uCombPur,Hstack,R0P,[xPur;xEva]);
            [xstackE,P_eva]=linearKFStep(x0e,zEva,Astack,Bstack,GammaLinStackE,P_eva,Q0stack,uCombEva,Hstack,R0E,[xPur;xEva]);
        end
        
        
        xPpur=xstackP(1:nX); xEpur=xstackP(nX+1:end); xPeva=xstackE(1:nX); xEeva=xstackE(nX+1:end);
        
        %Various outputs
        uP=uPurTrue
        uE=uEvaTrue
        xPurLoc=xPur
        xEvaLoc=xEva
        
        if flagPlotAsGo==1
            figure(42);delete(gameStatePlot)
            pause(plotint)
            gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red')
        end
        
        
    end
    nsimTrunc(MCL)=numSimRuns;
end


% 
%         if max(uPconst)>umaxPur || max(uEconst)>umaxEva
%             fprintf('Warning: Exceeding maximum thrust\n');
%         end
%         
%         [a1,a2,a3,~]=size(uconPtemp);
%         uconPcull=zeros(a1,a2,a3,2);
%         [a1,a2,a3,~]=size(uconEtemp);
%         uconEcull=zeros(a1,a2,a3,2);
%         nmod_pur_temp=length(uconPtemp); nmod_eva_temp=length(uconEtemp);
%         %remove collision points from map; technically, each party should
%         %create their own uconCull mats but whatever
%         mapEdges; %load wall locations into memory
%         nearestWallIndex=0; nearestWallDist=0; nearestWallPointPair=0;
%         numSafeControls=0;
%         for i1=1:nmod_pur_temp
%             uhist=uconPtemp(:,:,:,i1);  %Should check for collisions from x0=left of quad and x0=right of quad
%             isControlSafe=isSafeFromFutureCollision(ttf,xPpur,A_tr,B_tr,uhist,wallPoints,numObj,virtualWallWidth);
%             if isControlSafe==1
%                 numSafeControls=numSafeControls+1;
%                 uconPcull(:,:,:,numSafeControls)=uconPtemp(:,:,:,i1);
%             end
%         end
%         uconPcull=uconPcull(:,:,:,1:numSafeControls); %removes [0;0] controls added when initializing
%         nmod_pur=length(uconPcull);
%         numSafeControls=0;
%         nWDstore=[];
%         for i1=1:nmod_eva_temp
%             uhist=uconEtemp(:,:,:,i1);
%             [isControlSafe,nWI,nPPI,nWD]=isSafeFromFutureCollision(ttf,xEeva,A_tr,B_tr,uhist,wallPoints,numObj,virtualWallWidth);
%             nWDstore=[nWDstore nWD];
%             if isControlSafe
%                 %fprintf('foundsafeE\n')
%                 numSafeControls=numSafeControls+1;
%                 uconEcull(:,:,:,numSafeControls)=uconEtemp(:,:,:,i1);
%             end
%         end
%         uconEcull=uconEcull(:,:,:,1:numSafeControls); %removes [0;0] controls added when initializing
%         nmod_eva=length(uconEcull);
%         
%         hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
%         hE=ttf*ones(nmod_eva,1);
%         KmatP=[];
%         KmatE=[];
%         
%         uclassVecPpur=zeros(nmod_pur,1); uclassVecEpur=zeros(nmod_eva,1);
%         uclassVecPeva=zeros(nmod_pur,1); uclassVecEeva=zeros(nmod_eva,1);
%         
%         strategiesE=struct('constant',uconEcull,'horizon',hE);
%         strategiesP=struct('constant',uconPcull,'horizon',hP);
%         %KmatE, KmatP are nu x nx x T x nmodE/nmodP feedback matrices,
%         %where T is the maximum time horizon considered
%         [JpPur,JePur]=generateCostMatrices(strategiesP,strategiesE,xPpur,xEpur,0);
%         [JpEva,JeEva]=generateCostMatrices(strategiesP,strategiesE,xPeva,xEeva,0);
%         
%         %convert cost matrices to payoff matrices
%         VpPur = -JpPur;
%         VePur = -JePur;
%         VpEva = -JpEva;
%         VeEva = -JeEva;
%         [eqLocP,nashReturnFlagP,~]=findRDEq(VpPur,VePur);
%         [eqLocE,nashReturnFlagE,~]=findRDEq(VpEva,VeEva);
%         %Process equilibria into controls
%         %NOTE: For readability, this could be made into a function to be
%         %called twice but this would require passing large matrices to a
%         %function that may not use them.  The choice is deliberate here.
%         if nashReturnFlagP>=1 %if there IS an RDEq
%             uClassPp=eqLocP(1,1);
%             uClassEp=eqLocP(2,1);
%             [uPurTrue,uEvaExpected] = processNashType1(eqLocP,umaxPur,umaxEva,uconPcull,uconEcull);
%         elseif nashReturnFlagP==0 %suboptimal
%             fprintf('Running LH2\n')
%             %split for efficiency
%             [uPurTrue,uEvaExpected] = processNashType0(VpPur,VePur,umaxPur,umaxEva,uconPcull,uconEcull);
%         end
%         if nashReturnFlagE>=1
%             [uPurExpected,uEvaTrue] = processNashType1(eqLocE,umaxPur,umaxEva,uconPcull,uconEcull);
%         elseif nashReturnFlagE==0 %suboptimal
%             fprintf('Running LH2')
%             [uPurExpected,uEvaTrue] = processNashType0(VpEva,VeEva,umaxPur,umaxEva,uconPcull,uconEcull);
%         end





