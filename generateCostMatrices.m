function [Jpur,Jeva] = generateCostMatrices(strategiesP,strategiesE,xP,xE,noiseTypeFlag)
%inputs:
%     strategiesP (pursuer) is a struct containing
%         'constant': T x nmodP  constant max thrust vectors such that
%             u=uconst*unit(x)
%         'horizon': nmodP x 1  vector of integers containing finite-time horizon
%             length of cost function to evaluate, from 1 to T
%         for all pursuit strategies
%     strategiesE is the same for all evader strategies
%     xhat0 is the ehat of the player computing J matrices
%     noiseTypeFlag refers to how the cost due to noise is to be
%         calculated.
%             noiseTypeFlag=0 means ignore noise terms in J (sunk cost)
%             noiseTypeFlag=1 means to calculate them
%             This may also be done with an internal Monte Carlo but I prefer to be
%             efficient with my code.
% outputs:
%     Jpur, Jeva are cost matrices (NOT payoff matrices) for PE players
%         row player is pursuer, column player is evader
%         size nmodP x nmodE.
%         These matrices are E(J) (as the game is stochastic).  Thus, the
%         noise cost is higher for higher time horizons.

%notes:
%When calculating J, if the two time horizons for different strategies do
%not match, then J=9001^2;
%nU=nX/2 is assumed

nmodP=length(strategiesP.horizon);
nmodE=length(strategiesE.horizon);

global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG RPurGG nXGG QnoisePGG dtGG cddGG nUGG QnoiseEGG;
QfEva=QfEvaGG; QfPur=QfPurGG; QstepEva=QstepEvaGG; dt=dtGG;
QstepPur=QstepPurGG; REva=REvaGG; RPur=RPurGG; QnoiseP=QnoisePGG; QnoiseE=QnoiseEGG;
nX=nXGG; cdd=cddGG; nU=nUGG;

Jpur=9001*ones(nmodP,nmodE);
Jeva=42*ones(nmodP,nmodE);

eyeHalfNX=eye(nX/2);
zerosHalfNX=zeros(nX/2,nX/2);

A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

%Precomputing sunk noise cost for speed.  This saves about 30% of time
%spent in-file.  Other parts of time cost involve feedback matrices and
%cannot be precomputed.
noise1p=trace(QfPur*Gammak*QnoiseP*Gammak');
noise2p=trace(QstepPur*Gammak*QnoiseP*Gammak');
noise1e=trace(QfEva*Gammak*QnoiseE*Gammak');
noise2e=trace(QstepEva*Gammak*QnoiseE*Gammak');

noiseCostPur=3409;
noiseCostEva=3409;
for i=1:nmodP
    for j=1:nmodE
        e=xE-xP;
        Jpurloc=0;
        Jevaloc=0;
        
        TT=min(strategiesP.horizon(i),strategiesE.horizon(j));
        if TT ~= max(strategiesP.horizon(i),strategiesE.horizon(j)) %only compare strategies of equal length
            Jpur(i,j)=9001;
            Jeva(i,j)=9001;
        else
            TT=max(TT,1); %ensure that TT>=1
            
            for k=1:TT
                %use terminal vs. stepping cost
                if k==TT
                    QuseP=QfPur;
                    QuseE=QfEva;
                else
                    QuseP=QstepPur;
                    QuseE=QstepEva;
                end
                
                uconP=strategiesP.constant(:,:,:,i);                
                uconE=strategiesE.constant(:,:,:,j);

                %generate the noise contribution of the cost matrix
                if noiseTypeFlag==0 %noise contribution for pursuer
                    noiseCostPur=0;
                elseif noiseTypeFlag==1
                    if k==TT
                        noiseCostPur=noise1p;
                    else
                        noiseCostPur=noise2p;
                    end
                end
                if noiseTypeFlag==0
                    noiseCostEva=0;
                elseif noiseTypeFlag==1
                    if k==TT
                        noiseCostEva=noise1e;
                    else
                        noiseCostEva=noise2e;
                    end
                end
                
                %generate and saturate u
                up=uconP(:,:,k);
                ue=uconE(:,1,k);
                %up=vectorSaturationF(up,0,umaxPur);
                %ue=vectorSaturationF(ue,0,umaxEva);
                
                %advance dynamics
                e=A_tr*e+B_tr*ue-B_tr*up;
                
                %calculate costs
                Jpurloc=Jpurloc + noiseCostPur + e'*QuseP*e + up'*RPur*up;
                Jevaloc=Jevaloc + noiseCostEva - e'*QuseE*e + ue'*REva*ue;
            end
            Jpur(i,j)=Jpurloc;
            Jeva(i,j)=Jevaloc;
        end
        
        
    end
end






end

