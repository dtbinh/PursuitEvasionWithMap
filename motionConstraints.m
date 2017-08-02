function [Aineq,bineq] = motionConstraints(n,t_w,n_w,constraints,n_points)

pmin = constraints(1);
pmax = constraints(2);
vmin = constraints(3);
vmax = constraints(4);
amin = constraints(5);
amax = constraints(6);
jmin = constraints(7);
jmax = constraints(8);

time = linspace(t_w(1), t_w(end), n_points);

Aineq = [];
bineq = [];
k = 1;
for i = 1:numel(time)
    if ~(time(i) <= t_w(k+1))
        k = k+1;
    end
    [Tx,Tv,Ta,Tj] = TimeVectors(time(i),n);
    Aineq = [Aineq;
             zeros(1,(k-1)*(n+1)) -Tx' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1))  Tx' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1)) -Tv' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1))  Tv' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1)) -Ta' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1))  Ta' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1)) -Tj' zeros(1,(n_w-k-1)*(n+1));
             zeros(1,(k-1)*(n+1))  Tj' zeros(1,(n_w-k-1)*(n+1));];
         
    bineq = [bineq;
             -pmin;
             pmax;
             -vmin
             vmax
             -amin
             amax
             -jmin
             jmax];
end