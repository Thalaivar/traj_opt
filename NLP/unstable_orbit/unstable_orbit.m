N = 100;

y0 = [0.1;0.1];
tspan = [0, 2];
[D, cheb_x] = cheb_diff(N);
cheb_t = flip(((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2);

[t, state] = ode45(@system, cheb_t, y0);

%r = state(:,1); theta = state(:,2);
%x = zeros(size(r)); y = zeros(size(r));
%for i = 1:length(r)
%   x(i,1) = r(i, 1)*cos(theta(i, 1)); y(i,1) = r(i, 1)*sin(theta(i, 1));
%end
%plot(t, r, '.-b', t, theta, '-r')

ode_samples = state;
for i = 1:length(ode_samples)
    theta = ode_samples(i,2);
    if theta > 2*pi, theta = theta - 2*pi;, end
    if theta < 0, theta = 2*pi + theta;, end
    ode_samples(i,2) = theta;
end
x0 = zeros(2*N + 2 + 1, 1);
x0(1:N+1, 1) = ode_samples(:,1); x0(N+2:2*N+2, 1) = ode_samples(:,2);
x0(end,1) = -1*(cheb_t(1) - cheb_t(end));
lb = ones(2*N + 2 + 1, 1); ub = ones(2*N + 2 + 1, 1);
lb(1:N+1, 1) = 0*lb(1:N+1, 1); ub(1:N+1, 1) = 2*ub(1:N+1, 1);
lb(N+2:2*N+2, 1) = 0*lb(N+2:2*N+2, 1); ub(N+2:2*N+2, 1) = 2*pi*ub(N+2:2*N+2, 1);
lb(end,1) = 2; ub(end,1) = Inf;

options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 200000, 'StepTolerance', 1e-15, 'MaxIterations', 1000);
[x, fval] = fmincon(@objfun, x0, [], [], [], [], lb, ub, @state_const, options);

r = x(1:N+1,1); theta = x(N+2:2*N+2,1);
x1 = zeros(size(r)); x2 = zeros(size(r));
for i = 1:length(r)
   x1(i,1) = r(i, 1)*cos(theta(i, 1)); x2(i,1) = r(i, 1)*sin(theta(i, 1));
end
plot(x1,x2)

function ydot = system(t, y)
    r = y(1,1); theta = y(2,1);
    d = -2.0; nu = -pi; a = -3.0; c = pi; b = 5.0; omega = pi;
    
    rdot = d*nu*r + a*(r^3);
    thetadot = omega + c*nu + b*(r^2);
    ydot = [rdot; thetadot];
end

function [c, ceq] = state_const(x)
    N = (length(x) - 3)/2;
    rk = x(1:N+1, 1); thetak = x(N+2:2*N+2,1); T = x(end,1);
    
    [D, ~] = cheb_diff(N);
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
    
    %ceq(end+1) = x(1)-x(N+1);
    %ceq(end+1) = x(N+2)-x(2*N+2);
end

function f =  objfun(x)
    N = (length(x) - 3)/2;
    f = 10*(x(N+2,1)-x(2*N+2,1))^2 + 10*(x(1,1)-x(N+1,1))^2 + x(end)^2;
end
