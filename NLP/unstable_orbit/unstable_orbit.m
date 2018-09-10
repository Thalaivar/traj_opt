N = 100;

y0 = [1; 0];
tspan = [0, 3];
[D, cheb_x] = cheb_diff(N);
cheb_t = flip(((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2);

[t, state] = ode45(@system, tspan, y0);

ode_samples = interp1(t, state, cheb_t);
xy_samp = get_xy(ode_samples);
plot(xy_samp(:,1), xy_samp(:,2), '-b','Linewidth', 2.0)
hold on
x0 = zeros(2*N+2+1, 1); lb = ones(2*N+2+1, 1); ub = ones(2*N+2+1, 1);
x0(1:N+1, 1) = ode_samples(:,1); x0(N+2:2*N+2, 1) = ode_samples(:,2); x0(end,1) = cheb_t(end) - cheb_t(1);
lb(1:N+1, 1) = 0*lb(1:N+1, 1); lb(N+2:2*N+2, 1) = -3*lb(N+2:2*N+2, 1); lb(end,1) = 2;
ub(1:N+1, 1) = 2*ub(1:N+1, 1); ub(N+2:2*N+2, 1) = 7*ub(N+2:2*N+2, 1); ub(end,1) = 3;

options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 200000, 'StepTolerance', 1e-15, 'MaxIterations', 1000);
[x, fval] = fmincon(@objfun, x0, [], [], [], [], lb, ub, @state_const, options);

r = x(1:N+1,1); theta = x(N+2:2*N+2,1); T = x(end,1);
xy = get_xy([r,theta]);
plot(xy(:,1), xy(:,2), '.-r')

function [c, ceq] = state_const(x)
    N = (length(x) - 3)/2;
    rk = x(1:N+1, 1); thetak = x(N+2:2*N+2,1); T = x(end,1);
    
    [D, cheb_x] = cheb_diff(N);
    %tspan = [0, 3];
    %cheb_t = flip(((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2);

    diffmat = [D, zeros(N+1, N+1); zeros(N+1, N+1), D];
    xdot_cap = (2/T)*diffmat*[rk;thetak];
    
    fr = zeros(N+1,1); ftheta = zeros(N+1,1);
    for i = 1:N+1
        temp = system(0, [rk(i);thetak(i)]);
        fr(i, 1) = temp(1,1); ftheta(i, 1) = temp(2,1);
    end
    xdot = [fr;ftheta];
    
    c = [];
    ceq = xdot_cap - xdot;
    %[M, index] = min(abs(cheb_t - T));
    %index = N+1;
    %ceq(end+1) = x(1,1) - x(index,1);
    %ceq(end+1) = abs(x(N+2,1) - x(N+1+index,1)) - 2*pi;
    
    
end

function ydot = system(t, y)
    r = y(1,1); theta = y(2,1);
    a = -3; b = 5; c = pi; d = -2; nu = -pi; omega = pi;
    rdot = d*nu*r + a*r^3;
    thetadot = omega + c*nu + b*r^2;
    ydot = [rdot;thetadot];
end

function theta_c = constrain_theta(theta)
    theta_c = theta;
    if theta > 2*pi
        n = floor(theta/(2*pi));
        theta_c = theta - n*2*pi;
    end
    
    if theta < 0
        n = floor(theta/(2*pi));
        theta_c = -n*2*pi + theta;
    end
end

function xy = get_xy(state)
    r = state(:,1); theta = state(:,2);
    xy = zeros(size(state));
    for i = 1:length(r)
      xy(i,1) = r(i)*cos(theta(i)); xy(i,2) = r(i)*sin(theta(i));
    end
end

function f = objfun(x)
    N = (length(x)-3)/2;
    f = (x(1,1) - x(N+1,1))^2 + (abs(x(N+2,1) - x(2*N+2,1)) - 2*pi)^2 + x(end)^2;
end