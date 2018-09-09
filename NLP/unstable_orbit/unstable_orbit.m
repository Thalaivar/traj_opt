N = 100;

y0 = [0.1;0.1];
tspan = [0, 3];
[D, cheb_x] = cheb_diff(N);
cheb_t = ((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2;

[t, ode_samples] = ode45(@system, tspan, y0);

r = ode_samples(:,1); theta = ode_samples(:,2);
x = zeros(size(r)); y = zeros(size(r));
for i = 1:length(r)
   x(i,1) = r(i, 1)*cos(theta(i, 1)); y(i,1) = r(i, 1)*sin(theta(i, 1));
end
plot(x, y)

x0 = zeros(2*N + 2 + 1, 1);
x0(1:N+1, 1) = ode_samples(:,1); x0(N+2:2*N+2, 1) = ode_samples(:,2);
x0(end,1) = cheb_t(1) - cheb_t(end);
lb = ones(2*N + 2 + 1, 1); ub = ones(2*N + 2 + 1, 1);

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
    
    fx = zeros(N+1,1); fy = zeros(N+1,1);
    for i = 1:N+1
        temp = system(0, [xk(i);yk(i)]);
        fx(i, 1) = temp(1,1); fy(i, 1) = temp(2,1);
    end
    xdot = [fx;fy];
    
    c = [];
    ceq = xdot_cap - xdot;
    
    ceq(end+1) = x(1)-x(N+1);
    ceq(end+1) = x(N+2)-x(2*N+2)
end

function f =  objfun(x)
    f = x(end)^2;
end
