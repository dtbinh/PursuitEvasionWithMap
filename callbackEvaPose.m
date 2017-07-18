function callbackEvaPose(~,msg)
global xPeva xEeva
xPeva=msg.xPur.position;
xEeva=msg.xEva.position;


end

