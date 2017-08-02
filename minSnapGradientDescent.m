function [SolCoeff,Cost,t_w] = minSnapGradientDescent(Waypoints)

%Solve min snap with initial guess of times
[SolCoeff,Cost] = solveMinSnap(Waypoints);

%Backtracking gradient descent
epsilon = 0.05;
step = inf;
tf = Waypoints.t(end);
dt_w = diff(Waypoints.t);
Waypoints_new = Waypoints;
Waypoints_g = Waypoints;
while(step > epsilon)
    %%Calculate gradient vector
    m = numel(Waypoints.t)-1;
    g = (-1/(m-1))*ones(1,m);
    h = 1e-5;
    for i = 1:m
        gi = g;
        gi(i) = 1;
        dt_wi = dt_w + h*gi;
        Waypoints_new = Waypoints;
        Waypoints_new.t = [0 cumsum(dt_wi)];
        [~,Cost_i] = solveMinSnap(Waypoints_new);
        gradf(1,i) = (Cost_i-Cost)/h;
    end
    
    %Normalize gradient vector
    gradf = min(abs(dt_w))*gradf/max(abs(gradf));
    
    %Backtracking line search
    powers = 1:6;
    alpha = 0.9*(0.5.^powers);
    CurCost = inf*ones(1,numel(alpha));
    dt_wi_alpha = zeros(numel(alpha),m);
    for i = 1:numel(alpha)

        dt_wi_alpha(i,:) = dt_w - alpha(i)*gradf;
        
        %Renormalize dt vector such that final time (tf) doesn't change
        dt_wi_alpha(i,:) = tf*dt_wi_alpha(i,:)/(sum(dt_wi_alpha(i,:)));
        
        %Solve for the current cost
        Waypoints_new = Waypoints;
        Waypoints_new.t = [0 cumsum(dt_wi_alpha(i,:))];
% tic
        [~,CurCost(i)] = solveMinSnap(Waypoints_new);
%  toc
        %If solution is better, save best solution and exit loop
        if((CurCost(i) < Cost))
            break;
        end
    end

    if(min(CurCost) < Cost)
        index = find(CurCost == min(CurCost));
        dt_w = dt_wi_alpha(index,:);
        step = (Cost - CurCost(index))/Cost;
        Cost = CurCost(index);
    else
        break;
    end

end

% dt_w = tf*dt_w/sum(dt_w);
Waypoints.t = [0 cumsum(dt_w)];
[SolCoeff,FinalCost] = solveMinSnap(Waypoints);
t_w = Waypoints.t;

