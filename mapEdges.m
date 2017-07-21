%reads 2D maps of points containing connected vertices of objects.  The
%pursuer/evader will generate actions that do not pass through walls.

%read in matrics as a row of points (ie 2x4 for a rectangle in R2)
a=[[0;0] [5;0] [5;2] [0;2]];
%b=[[1;.9] [4;.9] [4;1.1] [1;1.1]];
axisSet=[-1 6 -1 3]; %axes for use in wall plot
%wallPoints={a,b};
wallPoints={a};
numObj=length(wallPoints);



