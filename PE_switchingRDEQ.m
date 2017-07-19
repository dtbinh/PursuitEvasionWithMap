clear;clc;
%closes the loop for control switching by risk-dominance criterion

%NOTES
% Swap so as to use MC to generate OPTIMAL play and non-MC to generate ACTUAL play

tmax=1;
lookAheadTime=1; %number of steps to look ahead in optimization; set 0 to use t_remaining instead
dt=1;
MCmax=1;
nSim=ceil(tmax/dt); %simulation time within MC
%nSim=1;

flagDiminishHorizonNearNSim=1;  %if ==1, consider lookAheadTime
         %if TRemain>lookAheadTime, use TRemain if TRemain<lookAheadTime
flagUseDetm=0;  %zero noise if flag==1
flagUseFeedback=0; %generate control history via feedback
flagUseModeControlInsteadOfMean=1;

muAllP=zeros(MCmax,nSim,2);
muAllE=zeros(MCmax,nSim,2);
nsimTrunc=zeros(MCmax,1);
captureIndex=zeros(MCmax,1);
flagPlotAsGo=0;  %=1 to plot motion
flagPlotFinal=0; %=1 to plot final dist/speed history
plotint=.1;

nX=2;   %full size of state space (even number)
nU=nX/2;
nV=nX/2;
nZ=nX/2;

MC_rdeq_Max=100; %length of MC used to generate payoffs
if flagUseDetm==1
    MC_rdeq_Max=1;
end

eStore=zeros(nX,nSim,MCmax);

for MCL=1:MCmax
    
    
    uphist=[];
    uehist=[];
    upDethist=[];
    ueDethist=[];
    
    C=300; %cost offset to guarantee positivity of V=C-J
    cdd=0;
    eyeHalfNX=eye(nX/2);
    zerosHalfNX=zeros(nX/2,nX/2);
    A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
    B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    H=[eyeHalfNX zeros(nX/2)];
    P0=diag([.1*ones(1,nX/2) .1*ones(1,nX/2)]);
    Q0=.01*eyeHalfNX; %process noise
    Q0s=Q0*sqrt(2);
    R0P=.05*eye(nZ); %measurement noise
    R0E=R0P;
    cholQ0_T=chol(Q0)';
    cholR0P_T=chol(R0P)';
    cholR0E_T=chol(R0E)';
    P_pur=P0; P_eva=P0;
    
    %cost matrices
    QfPur=8*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepPur=zeros(nX,nX);
    RPur=110*eye(nU);  %NOTE: discrete R = continuous R*dt^2 %it IS dt^2, NOT dt^2/2
    RPurScaled=RPur*dt^2;
    
    QfEva=5*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepEva=zeros(nX,nX);
    REva=125*eye(nU);
    REvaScaled=REva*dt^2;
    global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG ...
        RPurGG umaxPurGG umaxEvaGG nXGG QnoiseGG dtGG cddGG nUGG%#ok
    QfEvaGG=QfEva; QfPurGG=QfPur; QstepEvaGG=QstepEva; QnoiseGG=Q0; cddGG=cdd;
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
    
    if nX==2
        xEva=[10;0];
        xPur=[0;.1];
    else
        xEva=[5.000;3;0;0];
        xPur=[0;0;0;.1];
    end
    umaxPur=.5;
    umaxEva=.5;
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
    
    
    %set randn's seed RNG
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    ehat0=xEva-xPur+(chol(P0))'*randn(nX,1);
    if flagUseDetm==1
        ehat0=xEva-xPur; %deterministic behavior
    end
    ehat_prev_eva=ehat0;
    
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
    
    %----- Random number seed
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
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
            ehatP=ehat0;
            for j=1:nmod_pur
                xhathist_pur(:,j,1)=ehat0;
            end
            ehatE=ehat0;
            for j=1:nmod_eva
                xhathist_eva(:,j,1)=ehat0;
            end
        end
        
       
        
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
            
%             k1 = 1.0e-04 *[-0.4066   -0.0004]*ehat_prev_pur;
%             k2 = 1.0e-04 *[-0.3578   -0.0003]*ehat_prev_pur;
%             u1p=[(k1-.0005):.0001:(k1+.005) (k2-.0005):.0001:(k2+.005) 0:.001:.2]; %.177
%             u1e=u1p; %.097
            uPtheta=0:15*pi/180:2*pi; %.177
            uPconst=0:.05:.25;
            uEtheta=0:.01:.25; %.097
            uEconst=0:.05:.25;
            
            uconP=uConMat_Rt(ttf,umaxPur,uPtheta,uPconst);
            uconE=uConMat_Rt(ttf,umaxEva,uEtheta,uEconst);
        end
        
        
        nmod_pur=length(uconP);
        nmod_eva=length(uconE);
        
        %remove collisions
        for i1=1:nmod_pur
            uhist=uconP(:,:,:,i1);
            
            
        end
            
        hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
        hE=ttf*ones(nmod_eva,1);
        KmatP=[];
        KmatE=[];
        typesE=ones(ttf,nmod_eva);
        typesP=ones(ttf,nmod_pur);
        
        uPurTrue=zeros(nU,1); uEvaTrue=uPurTrue; uEvaExpected=uPurTrue; uPurExpected=uPurTrue;
        uclassVecPpur=zeros(nmod_pur,1); uclassVecEpur=zeros(nmod_eva,1);
        uclassVecPeva=zeros(nmod_pur,1); uclassVecEeva=zeros(nmod_eva,1);
        
        %counter variables for Truth/eXpected for mode control instead of
        %mean control
        uPTct=zeros(nmod_pur,1); uPXct=zeros(nmod_pur,1);
        uETct=zeros(nmod_eva,1); uEXct=zeros(nmod_eva,1);
        
        cctemp=1/MC_rdeq_Max; %normalization constant
        for mcpayoffdumvar=1:MC_rdeq_Max+1 %MC for payoffs
            %if noise cost is flat across all time steps then ignore noise
            %find true optimal play on last iteration
            minTE=min(typesE); minTP=min(typesP); maxTE=max(typesE); maxTP=max(typesP);
            if max(hP)==min(hP) && max(hE)==min(hE) && max(maxTE)==min(minTE)...
                    && max(maxTP)==min(minTP) && min(minTP)==min(minTE)
                noiseIsSunkCost=0; %if all strategies for all times are same type AND same horizon
            else
                noiseIsSunkCost=1;
            end
            
            strategiesE=struct('matrices',KmatE,'constant',uconEcull,'horizon',hE,'types',typesE);
            strategiesP=struct('matrices',KmatP,'constant',uconPcull,'horizon',hP,'types',typesP);
            %KmatE, KmatP are nu x nx x T x nmodE/nmodP feedback matrices,
            %where T is the maximum time horizon considered
            if flagUseDetm==1 || mcpayoffdumvar==MC_rdeq_Max+1
                noiseP=zeros(nX,1); noiseE=zeros(nX,1);
            else
                noiseP=chol(P_pur)*randn(nX,1); noiseE=chol(P_eva)*randn(nX,1);
            end
            ehatPsim=ehatP+noiseP;
            ehatEsim=ehatE+noiseE;
            [JpPur,JePur]=generateCostMatrices(strategiesP,strategiesE,xPpur,xEpur);
            [JpEva,JeEva]=generateCostMatrices(strategiesP,strategiesE,xPeva,xEeva);
            
            %convert cost matrices to payoff matrices
            VpPur=C-JpPur;
            VePur=C-JePur;
            VpEva=C-JpEva;
            VeEva=C-JeEva;
            [eqLocP,nashReturnFlagP,~]=findRDEq(VpPur,VePur);
            [eqLocE,nashReturnFlagE,~]=findRDEq(VpEva,VeEva);
            if nashReturnFlagP>=1 %if there IS an RDEq
                uClassPp=eqLocP(1,1);
                uClassEp=eqLocP(2,1);
                if flagUseModeControlInsteadOfMean==1
                    uptemp=zeros(nmod_pur,1); uptemp(uClassPp)=1;
                    uetemp=zeros(nmod_eva,1); uetemp(uClassEp)=1;
                    uPTct=uPTct+uptemp; uEXct=uEXct+uetemp;
                else
                    if typesP(uClassPp)==1
                        uPurTrueTemp = vectorSaturationF(uconP(:,:,1,uClassPp),0,umaxPur);
                    else
                        uPurTrueTemp = vectorSaturationF(KmatP(:,:,1,uClassPp)*ehatP,0,umaxPur);
                    end
                    if typesE(uClassEp)==1
                        uEvaExpectedTemp = vectorSaturationF(uconE(:,:,1,uClassEp),0,umaxEva);
                    else
                        uEvaExpectedTemp = vectorSaturationF(KmatE(:,:,1,uClassEp)*ehatE,0,umaxEva);
                    end
                end
            elseif nashReturnFlagP==0 %suboptimal
                fprintf('Running LH2')
                nashMixedP=LH2(VpPur,VePur);
                KpM=zeros(size(KmatP(:,:,1,1)));
                KeM=zeros(size(KmatE(:,:,1,1)));
                nashP=nashMixedP{1};
                nashE=nashMixedP{2};
                u0p=zeros(nU,1);
                u0e=zeros(nU,1);
                if flagUseModeControlInsteadOfMean==1
                    uptemp=nashP;
                    uetemp=nashE;
                    uPXct=uPXct+uptemp;
                    uETct=uETct+uetemp;
                else
                    for kk=1:length(nashP)
                        if nashP(kk) > 0
                            if typesP(kk)==0
                                u0p=u0p+nashP(kk)*KmatP(:,:,1,kk)*ehatP;
                            else
                                u0p=u0p+nashP(kk)*uconP(:,:,1,kk);
                            end
                        end
                    end
                    for kk=1:length(nashE)
                        if nashE(kk) > 0
                            if typesE(kk)==0
                                u0e=u0e+nashE(kk)*KmatE(:,:,1,kk)*ehatP;
                            else
                                u0e=u0e+nashE(kk)*uconE(:,:,1,kk);
                            end
                        end
                    end
                    uPurTrueTemp=vectorSaturationF(u0p,0,umaxPur);
                    uEvaExpectedTemp=vectorSaturationF(u0e,0,umaxEva);
                end
            end
            
            if nashReturnFlagE>=1
                uClassPe=eqLocE(1,1);
                uClassEe=eqLocE(2,1);
                if flagUseModeControlInsteadOfMean==1
                    uptemp=zeros(nmod_pur,1); uptemp(uClassPe)=1;
                    uetemp=zeros(nmod_eva,1); uetemp(uClassEe)=1;
                    uETct=uETct+uetemp; uPXct=uPXct+uptemp;
                else
                    if typesP(uClassPe)==1
                        uPurExpectedTemp = vectorSaturationF(uconP(:,:,1,uClassPe),0,umaxPur);
                    else
                        uPurExpectedTemp = vectorSaturationF(KmatP(:,:,1,uClassPe)*ehatE,0,umaxPur);
                    end
                    if typesE(uClassEe)==1
                        uEvaTrueTemp = vectorSaturationF(uconE(:,:,1,uClassEe),0,umaxEva);
                    else
                        uEvaTrueTemp = vectorSaturationF(KmatE(:,:,1,uClassEe)*ehatE,0,umaxEva);
                    end
                end
            elseif nashReturnFlagE==0 %suboptimal
                fprintf('Running LH2')
                nashMixedE=LH2(VpEva,VeEva);
                KpM=zeros(size(KmatP(:,:,1,1)));
                KeM=zeros(size(KmatE(:,:,1,1)));
                nashP=nashMixedE{1};
                nashE=nashMixedE{2};
                u0p=zeros(nU,1);
                u0e=zeros(nU,1);
                if flagUseModeControlInsteadOfMean==1
                    uptemp=nashP;
                    uetemp=nashE;
                    uPXct=uPXct+uptemp;
                    uETct=uETct+uetemp;
                else
                    for kk=1:length(nashP)
                        if nashP(kk) > 0
                            if typesP(kk)==0
                                u0p=u0p+nashP(kk)*KmatP(:,:,1,kk)*ehatE;
                            else
                                u0p=u0p+nashP(kk)*uconP(:,:,1,kk);
                            end
                        end
                    end
                    for kk=1:length(nashE)
                        if nashE(kk) > 0
                            if typesE(kk)==0
                                u0e=u0e+nashE(kk)*KmatE(:,:,1,kk)*ehatE;
                            else
                                u0e=u0e+nashE(kk)*uconE(:,:,1,kk);
                            end
                        end
                    end
                    uPurExpectedTemp=vectorSaturationF(u0p,0,umaxPur);
                    uEvaTrueTemp=vectorSaturationF(u0e,0,umaxEva);
                end
                
            end
            if mcpayoffdumvar<=MC_rdeq_Max
                if flagUseModeControlInsteadOfMean~=1
                    uPurExpected=uPurExpected+cctemp*uPurExpectedTemp;
                    uEvaExpected=uEvaExpected+cctemp*uEvaExpectedTemp;
                    uPurTrue=uPurTrue+cctemp*uPurTrueTemp;
                    uEvaTrue=uEvaTrue+cctemp*uEvaTrueTemp;
                end
            else
                if flagUseModeControlInsteadOfMean==1
                    uPurTrueTemp=uconP(:,1,1,uClassPp);
                    uEvaTrueTemp=uconE(:,1,1,uClassEe);
                end
                upDethist=[upDethist uPurTrueTemp];
                ueDethist=[ueDethist uEvaTrueTemp];
            end
        end
        
        
        if flagUseModeControlInsteadOfMean==1
            [~,indPT]=max(uPTct);
            [~,indPX]=max(uPXct);
            [~,indET]=max(uETct);
            [~,indEX]=max(uEXct);
            uPurTrue=uconP(:,1,1,indPT);
            uPurExpected=uconP(:,1,1,indPX);
            uEvaTrue=uconE(:,1,1,indET);
            uEvaExpected=uconE(:,1,1,indEX);
        end
        
        uphist=[uphist uPurTrue];
        uehist=[uehist uEvaTrue];
        
        %process noise
        if flagUseDetm==1 %deterministic behavior
            nP=zeros(nX,1);
            nE=zeros(nX,1);
        else
            nP=Gammak*cholQ0_T*randn(nV,1);
            nE=Gammak*cholQ0_T*randn(nV,1);
        end
        
        
        %Propagate states
        xPur=A_tr*xPur+B_tr*uPurTrue+nP;
        xEva=A_tr*xEva+B_tr*uEvaTrue+nE;
        eTrue=xEva-xPur;
        
        eStore(:,i,MCL)=eTrue;
        
        JpurRunning=JpurRunning+uPurTrue'*RPurScaled*uPurTrue;
        JevaRunning=JevaRunning+uEvaTrue'*REvaScaled*uEvaTrue;

        %Measurement
        ze=H*(xEva-xPur)+cholR0E_T*randn(nZ,1);
        zp=H*(xEva-xPur)+cholR0P_T*randn(nZ,1);
        if flagUseDetm==1
            ze=H*(xEva-xPur);
            zp=H*(xEva-xPur); %deterministic behavior
        end
        
        %KF
%         ehat_prev_eva=eTrue+dt^2/2*chol(P0)*randn(nX,1);
%         ehat_prev_pur=eTrue+dt^2/2*chol(P0)*randn(nX,1);
        [ehatP,P_pur]=pointMassKFStep(ehatP,P_pur,Q0s,uEvaExpected,uPurTrue,zp,H,R0P,dt);
        [ehatE,P_eva]=pointMassKFStep(ehatE,P_eva,Q0s,uEvaTrue,uPurExpected,ze,H,R0E,dt);
        
        %Various outputs
%        kHist=KmatE(:,:,:,uClassEe)
%        kE=KmatE(:,:,1,uClassEe)
%        kP=KmatP(:,:,1,uClassPp)
        uP=uPurTrue
        uE=uEvaTrue
%        xPurLoc=xPur
%        xEvaLoc=xEva
        e=eTrue;
        vp=VpPur;
        ve=VePur;
        
        
        if flagPlotAsGo==1
            figure(42);delete(gameStatePlot)
            pause(plotint)
            gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red')
        end
        
        
        
        %capture times
        captureIndex(MCL)=i; %will be overwritten continuously until
        %the player is captured (the loop breaks at capture, after
        %incrementing captureIndex).
        
        %If captured at a time step or thresholds have flipped
        if norm(eTrue(1:nX/2))<=captureThresh
            fprintf('Captured \n')
            break
        end
        %If captured at a time step or thresholds have flipped
        if norm(eTrue(1:nX/2))<=captureThresh
            fprintf('Captured \n')
            break
        end
    end
    nsimTrunc(MCL)=numSimRuns;
    
    %for comparison with KumarV2 ONLY
    JpurFinal_minimized=JpurRunning+eTrue'*QfPur*eTrue %minimized
    JevaFinal_maximized=-JevaRunning+eTrue'*QfEva*eTrue %maximized
    
    
end
% 
% fastestEnd=min(captureIndex);
% 
% finalEndPlot=finalEnd-1;
% 
% if MCmax>1 %plotting rules differ for multiple vs. single simulation lengths
%     figure(1);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_pur
%         plot(tkhist(1:finalEndPlot),meanMuCalcP(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Pursuer strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('q1','q2')
%     
%     figure(2);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_eva
%         plot(tkhist(1:finalEndPlot),meanMuCalcE(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Evader strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('r1','r2')
% else
%     finalEndPlot=captureIndex(1)
%     figure(1);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_pur
%         plot(tkhist(1:finalEndPlot),meanMuCalcP(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Pursuer strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('q1','q2')
%     
%     figure(2);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_eva
%         plot(tkhist(1:finalEndPlot),meanMuCalcE(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Evader strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('r1','r2')
% end
if flagPlotFinal==1
    mState=mean(eStore,3);
    figure(2)
    subplot(2,1,1)
    plot(1:1:nSim,mState(1,:));
    title('Monte Carlo simulation')
    xlabel('Time step')
    ylabel('Distance')
    subplot(2,1,2)
    plot(1:1:nSim,mState(2,:));
    xlabel('Time step')
    ylabel('Relative speed')
end


figure(3);clf
subplot(2,1,1)
plot(dt*(1:1:numSimRuns),eStore(1,1:numSimRuns),'r')
hold on
plot(dt*(1:1:numSimRuns),eStore(2,1:numSimRuns),'b')
hLeg = legend('$\|e\|$','$\|\dot{e}\|$');
set(hLeg,'Interpreter','latex');
set(hLeg,'Interpreter','latex');
title('Relative pos and velocity')
subplot(2,1,2)
stairs(dt*(1:1:numSimRuns),uphist,'r')
hold on
stairs(dt*(1:1:numSimRuns),upDethist,'k')
hold on
stairs(dt*(1:1:numSimRuns),uehist,'b')
hold on
stairs(dt*(1:1:numSimRuns),ueDethist,'g')
hLeg = legend('$\|u_p\|$','$\|u_p\|_{opt}$','$\|u_e\|$','$\|u_e\|_{opt}$');
set(hLeg,'Interpreter','latex');
set(hLeg,'Interpreter','latex');
title('Applied control and optimal control')


% CODE GRAVEYARD
% generating/concatenating feedback controls
for indentFeedback=1:1
    if flagUseFeedback==1
        %generate greedy feedback matrices for matrix case
        ttf=ceil(tmax/dt - i);
        k1p = .05:.1:.25;
        k2p = [0 .1];
        possibleCombAtOneTimeStep=combvec(k1p,k2p);
        possibleCombosPrev=possibleCombAtOneTimeStep;
        for i0=1:ttf
            possibleCombosPrev=combvec(possibleCombosPrev, possibleCombAtOneTimeStep);
        end
        nnp=length(possibleCombosPrev);
        KpmatList=zeros(nU,nX,ttf,nnp);
        for i1=1:nnp
            for i2=1:ttf
                %grab the i2'th block of possibleCombosPrev
                if nX==2
                    nblock=nU*(i2-1);
                    KpmatList(:,:,i2,i1)=[possibleCombosPrev(nblock+1,i1)...
                        possibleCombosPrev(nblock+2,i1)];
                elseif nX==4
                    nblock=nU*(i2-1);
                    k1Ind=possibleCombosPrev(nblock+1,i1);
                    k2Ind=possibleCombosPrev(nblock+2,i1);
                    KpmatList(:,:,i2,i1)=[k1Ind 0 k2Ind 0; 0 k1Ind 0 k2Ind];
                end
            end
        end
        k1e=k1p;
        k2e=k2p;
        possibleCombAtOneTimeStep=combvec(k1e,k2e);
        possibleCombosPrev=possibleCombAtOneTimeStep;
        for i0=1:ttf
            possibleCombosPrev=combvec(possibleCombosPrev,possibleCombAtOneTimeStep);
        end
        nne=length(possibleCombosPrev);
        KematList=zeros(nU,nX,ttf,nne);
        for i1=1:nne
            for i2=1:ttf
                %grab the i2'th block of possibleCombosPrev
                if nX==2
                    nblock=nU*(i2-1);
                    KematList(:,:,i2,i1)=[possibleCombosPrev(nblock+1,i1)...
                        possibleCombosPrev(nblock+2,i1)];
                elseif nX==4
                    nblock=nU*(i2-1);
                    k1Ind=possibleCombosPrev(nblock+1,i1);
                    k2Ind=possibleCombosPrev(nblock+2,i1);
                    KematList(:,:,i2,i1)=[k1Ind 0 k2Ind 0; 0 k1Ind 0 k2Ind];
                end
            end
        end
        sizeKmat=size(KpmatList);
        nmod_pur=nnp;
        nmod_eva=nne;
        hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
        hE=ttf*ones(nmod_eva,1);
        KmatP=KpmatList;
        KmatE=KematList;
        typesE=zeros(ttf,nmod_eva);
        typesP=zeros(ttf,nmod_pur);
        uconE=[]; %negation is handled in cost generation code
        uconP=[];
    end
end



