
clear all
clc

N = 100;

% stable orbit
d = -2;
a = -3;
b = 5;
c = pi;
nu = -pi;
omega = pi;
p = [a, b, c, d ,nu, omega];

[D, cheb_x] = cheb_diff(N-1);
tspan = [0, 2];
Tguess = tspan(2) - tspan(1);
cheb_t = 0.5*Tguess*(1 - cheb_x);
y0 = [0.1;0];
odeopts = odeset('AbsTol',1e-11,'RelTol',1e-11,'MaxStep',0.001);
[t, xguess] = ode45(@(t, y) sys_model(t, y, p), cheb_t, y0, odeopts);

% first check for initial guess
%xsamp = xguess(:,1).*cos(xguess(:,2));
%ysamp = xguess(:,1).*sin(xguess(:,2));
%plot(xsamp, ysamp, 'b', 'Linewidth', 2.0 )
%hold on

x0 = [xguess(:,1);xguess(:,2);Tguess];
lb = ones(size(x0)); ub = Inf*ones(size(x0)); 
lb(1:N,1) = 0*lb(1:N,1); % r should always be > 0
lb(N+1:2*N,1) = -Inf*lb(N+1:2*N,1); lb(end,1) = 0;

% check if inital guess satisfies constraints
%[c, ceq] = state_const(x0, p, N, D);
%norm(ceq)/length(ceq)
alg = 'sqp';

options = optimoptions(@fmincon,'Algorithm',alg,'Display','iter');
options.MaxIterations = 200;
options.MaxFunctionEvaluations = 1000000;
x = fmincon(@(x) objfun(x, N), x0, [], [], [], [], lb, ub,@(x) state_const(x, p, N, D), options);

T = x(end,1);
fprintf('Time period = %.2f\n',T);

xpos = x(1:N).*cos(x(N+1:2*N));
ypos = x(1:N).*sin(x(N+1:2*N));
plot(xpos,ypos,'b', 'LineWidth', 2.0);
xlabel('x');
ylabel('y');
axis equal
grid on
hold on

r = x(1:N,1); theta = x(N+1:2*N,1);
index = floor(rand*100);
y0 = [r(index);theta(index)];
tspan = [0, 10];
odeopts = odeset('AbsTol',1e-11,'RelTol',1e-11,'MaxStep',0.001);
[t, x_prop] = ode45(@(t, y) sys_model(t, y, p), tspan, y0, odeopts);

xpos1 = x_prop(:,1).*cos(x_prop(:,2));
ypos1 = x_prop(:,1).*sin(x_prop(:,2));
plot(xpos1,ypos1,'.-r');
plot(y0(1,1)*cos(y0(2,1)), y0(1,1)*sin(y0(2,1)), 'sm')
xlabel('x');
ylabel('y');
axis equal
grid on


