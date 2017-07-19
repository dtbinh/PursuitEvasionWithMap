clear;clc;
%ros pose messages come in, controls go out

rosinit
global cddGG nXGG dtGG umaxPurGG umaxEvaGG
cddGG=0.4;
nXGG=2;
dtGG=1;
umaxPurGG=.5;
umaxEvaGG=.5;
%game variables
global xPpur xEpur xPeva xEeva
purPoseSub=rossubscriber('/posestackP',@callbackPurPose);  %must initialize ros publisher pose
evaPoseSub=rossubscriber('/posestackE',@callbackEvaPose);  %callback sets global pose vars

purControllerPub=rospublisher('controllerP','/CUSTOMMESSAGETYPE_CONTROL');
evaControllerPub=rospublisher('controllerE','/CUSTOMMESSAGETYPE_CONTROL');

cont=true;
n=0;
while cont
    n=n+1;
    
    xPpurTemp=xPpur;
    xEpurTemp=xEpur;
    xPevaTemp=xPeva;
    xEevaTemp=xEeva;
    
    uPc=.1:.075:.4; uEc=uPc;
    uPt=0:30*pi/180:2*pi-30*pi/180; uEt=uPt;
    [uPtempP,uEtempP]=generateCullSolveProcess(1,uPc,uPt,uEc,uTc,xPpurTemp,xEpurTemp,0);
    [uPtempE,uTtempE]=generateCullSolveProcess(1,uPc,uPt,uEc,uTc,xPevaTemp,xEevaTemp,0);
    
    %note: Header is empty, may be filled later for utility
    purControlMessage=rosmessage(controllerP);
    purControlMessage.uSelf=[uPtempP;0];
    purControlMessage.uOpponent=[uEtempP;0];
    evaControlMessage=rosmessage(controllerE);
    evaControlMessage.uSelf=[uEtempE;0];
    evaControlMessage.uOpponent=[uPtempE;0];
    
    send(purControllerPub,purControlMessage);
    send(evaControllerPub,evaControlMessage);
    
    if n>=1200
        cont=false;
    end
    pause(1);  %corresp spin()
end


%posestackP,posestackE are topics publishing a custom message type that
%includes estimated pose for both E and P




