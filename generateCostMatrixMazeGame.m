function [JP,JE] = generateCostMatrixMazeGame(tf,xvaj0p,xvaj0e,pointlistP,pointlistE)

%game constants
m_pur=1;
m_eva=1;
Rpur=1*m_pur; %control cost will be a multiple of weight
Reva=2*m_eva; %minsnap costs are often O(6)
Qfpur=5000;
Qfeva=4500;
Qpur=.01;
Qeva=0;

numPathsP=length(pointlistP);
numPathsE=length(pointlistE);

countpathsP=0; countpathsE=0;
for i1=1:numPathsP
    wpt=pointlistP{i1};
    countpathsP=countpathsP+length(wpt);
end
for i2=1:numPathsE
    wpt=pointlistE{i2};
    countpathsE=countpathsE+length(wpt);
end

JP=zeros(countpathsP,countpathsE);
JE=zeros(countpathsP,countpathsE);

tEhistmat={};uEhistmat={};

JuEmat=zeros(numPathsE,1);
cc=0;
xvajTemp=xvaj0e;
xyzEend=zeros(3,countpathsE);
p_xElist=cell(countpathsE,1); p_yElist=cell(countpathsE,1); p_zElist=cell(countpathsE,1);
dt = 0.01;
for iE=1:numPathsE %precompute one player for efficiency
    wpt=pointlistE{iE};
    for ii=1:length(wpt)
        cc=cc+1;
        %set up waypoint vector
        x_w=[xvajTemp(1) wpt(1,1:ii)];
        y_w=[xvajTemp(2) wpt(2,1:ii)];
        z_w=[xvajTemp(3) wpt(3,1:ii)];
        vx_w = nan(1,ii+1); vy_w = nan(1,ii+1); vz_w = nan(1,ii+1);
        ax_w = nan(1,ii+1); ay_w = nan(1,ii+1); az_w = nan(1,ii+1);
        jx_w = nan(1,ii+1); jy_w = nan(1,ii+1); jz_w = nan(1,ii+1);
        vx_w(1)=xvajTemp(4); vy_w(1)=xvajTemp(5); vz_w(1)=xvajTemp(6);
        ax_w(1)=xvajTemp(7); ay_w(1)=xvajTemp(7); az_w(1)=xvajTemp(9);
        jx_w(1)=xvajTemp(10); jy_w(1)=xvajTemp(11); jz_w(1)=xvajTemp(12);
        Waypoints.x = [x_w; vx_w; ax_w; jx_w];
        Waypoints.y = [y_w; vy_w; ay_w; jy_w];
        Waypoints.z = [z_w; vz_w; az_w; jz_w];
        Waypoints.t = linspace(0,tf,ii+1);
        
        %solve minsnap
        [SolCoeff,ju,t_w]=minSnapGradientDescent(Waypoints);
        JuEmat(cc)=Reva*ju;
        xyzEend(:,cc)=wpt(:,ii);
        
        %process into xyz coordinates
        if Qpur>0 || Qeva>0
            [p_x,p_y,p_z]=processSolCoeff(t_w,SolCoeff,dt);
            p_xElist{cc}=p_x; p_yElist{cc}=p_y; p_zElist{cc}=p_z;
        end
    end
end

cc=0;
xvajTemp=xvaj0p;
for iP=1:numPathsP
    wpt=pointlistP{iP};
    for ii=1:length(wpt)
        cc=cc+1;
        x_w=[xvajTemp(1) wpt(1,1:ii)];
        y_w=[xvajTemp(2) wpt(2,1:ii)];
        z_w=[xvajTemp(3) wpt(3,1:ii)];
        vx_w = nan(1,ii+1); vy_w = nan(1,ii+1); vz_w = nan(1,ii+1);
        ax_w = nan(1,ii+1); ay_w = nan(1,ii+1); az_w = nan(1,ii+1);
        jx_w = nan(1,ii+1); jy_w = nan(1,ii+1); jz_w = nan(1,ii+1);
        vx_w(1)=xvajTemp(4); vy_w(1)=xvajTemp(5); vz_w(1)=xvajTemp(6);
        ax_w(1)=xvajTemp(7); ay_w(1)=xvajTemp(7); az_w(1)=xvajTemp(9);
        jx_w(1)=xvajTemp(10); jy_w(1)=xvajTemp(11); jz_w(1)=xvajTemp(12);
        Waypoints.x = [x_w; vx_w; ax_w; jx_w];
        Waypoints.y = [y_w; vy_w; ay_w; jy_w];
        Waypoints.z = [z_w; vz_w; az_w; jz_w];
        Waypoints.t = linspace(0,tf,ii+1);
        
        [SolCoeff,JuP,t_w]=minSnapGradientDescent(Waypoints);
        xyzPend=wpt(:,ii);
        
        if Qpur>0 || Qeva>0
            [p_xP,p_yP,p_zP]=processSolCoeff(t_w,SolCoeff,dt);
        end
        
        for iE=1:countpathsE %can also do shortest maze-end solution for xPur to xEva instead of through-object linear distance
            JPrun=0; JErun=0;
            if Qpur>0 || Qeva>0
                e_xdiff=p_xP-p_xElist{iE};
                e_ydiff=p_yP-p_yElist{iE};
                e_zdiff=p_zP-p_zElist{iE};
                ex2=dot(e_xdiff,e_xdiff); ey2=dot(e_ydiff,e_ydiff); ez2=dot(e_zdiff,e_zdiff);
            end
            if Qpur>0
                JPrun=-Qpur*(ex2+ey2+ez2);
            end
            if Qeva>0
                JErun=Qeva*(ex2+ey2+ez2);
            end
            eDiff=xyzPend-xyzEend(:,iE);
            JP(cc,iE)=Rpur*JuP-eDiff'*Qfpur*eDiff+JPrun;
            JE(cc,iE)=JuEmat(iE)+eDiff'*Qfeva*eDiff+JErun;
        end
    end
end


end





