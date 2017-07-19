function uconP = uConMat_Rt(ttf,uconst,uthetha)
%ASSUMES PRE-SATURATED FOR COMPUTATIONAL EFFICIENCY
%generates uCon for nU=2 when u=R+theta
%uconP = nU x 1 x ttf x nMod  matrix of nMod possible control combinations

nU=2;

uconP=zeros(nU,1,ttf,(length(uconst)*length(uthetha))^ttf);
comb1vec=combvec(uconst,uthetha); %uconst becomes top row, utheta becomes bottom row

%generate stacked uconP matrix
if ttf==1
    uconP(:,:,1,:)=comb1vec(1,:).*[cos(comb1vec(2,:));sin(comb1vec(2,:))];
elseif ttf>=2
    possibleCombPrev=comb1vec;
    for i1=2:ttf
        possibleCombPrev=combvec(comb1vec,possibleCombPrev);
    end
    for i1=1:length(possibleCombPrev)
        for i2=0:floor(length(possibleCombPrev(:,i1))/nU)-1
            nblock=i2*nU+1:(i2+1)*nU;
            ttemp=possibleCombPrev(nblock,i1);
            uconP(:,1,i2+1,i1)=ttemp(1)*[cos(ttemp(2));sin(ttemp(2))];
        end
    end
end



end

