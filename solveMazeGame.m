function [nashP,nashE] = solveMazeGame(tf,xvajP,xvajE,pointlistP,pointlistE)

[Jpur,Jeva] = generateCostMatrixMazeGame(tf,xvajP,xvajE,pointlistP,pointlistE);
[nmod_pur,nmod_eva]=size(Jpur);

[eqLoc,nashReturnFlag,~]=findRDEq(-Jpur,-Jeva); %pass value matrices rather than cost matrices
[~,b]=size(eqLoc);
if b>=2
    numPoss=length(eqLoc);
    jpurTruncated=zeros(numPoss,numPoss); jevaTruncated=jpurTruncated;
    for j=1:numPoss
        for k=1:numPoss
            jpurTruncated(j,k)=-Jpur(eqLoc(1,j),eqLoc(2,k));
            jevaTruncated(j,k)=-Jeva(eqLoc(1,j),eqLoc(2,k));
        end
    end
    jpurMini=jpurTruncated;
    jevaMini=jevaTruncated;
end
if nashReturnFlag>=1 %if there IS an RDEq
    nashP=zeros(nmod_pur,1); nashP(eqLoc(1,1))=1;
    nashE=zeros(nmod_eva,1); nashE(eqLoc(2,1))=1;
elseif nashReturnFlag==0 %suboptimal
    fprintf('Running LH2\n')
    %split for efficiency
    nashMixedP=LH2(-Jpur,-Jeva);
    nashP=nashMixedP{1};
    nashE=nashMixedP{2};
end



end

