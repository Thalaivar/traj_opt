clear all

N = 100;

% simulate VDP system

% y0 = [0; 1];
y0 = [1 ; 1];

global nu
nu = 0.5; 

tspan = [0, 6.5];
% tspan = [6.961, 13.49];

[D, cheb_x] = cheb_diff(N);

% cheb_t = tspan;
% cheb_t = - ((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2;
cheb_t = 0.5*(1-cheb_x)*tspan(2);

[t, y] = ode45(@VDP, cheb_t, y0);

ode_samples = y;

% ode_samples = interp1(t, y, cheb_t);

% figure
% subplot(1,2,1)
% plot(cheb_t, ode_samples(:,1), '-r', cheb_t, ode_samples(:,2), '-b');
% xlim([0,cheb_t(end)]);
% subplot(1,2,2)
% plot(ode_samples(:,1),ode_samples(:,2));

%%
% alg = 'sqp';
alg = 'interior-point';

% setup NLP
x0 = ones(2*N + 2 + 1, 1);
x0(1:N+1, 1) = ode_samples(:, 1); x0(N+2:2*N + 2, 1) = ode_samples(:, 2);

x0(end) = cheb_t(end);

% x0(1,1) =  x0(N+1,1); x0(N+2,1) =  x0(2*N+2, 1);

lb = ones(2*N + 2 + 1, 1); lb(1:N+1, 1) = -3*lb(1:N+1, 1); lb(N+2:2*N + 2, 1) = -3*lb(N+2:2*N + 2, 1);
lb(end) = 0;
ub = ones(2*N + 2 + 1, 1); ub(1:N+1, 1) = 3*ub(1:N+1, 1); ub(N+2:2*N + 2, 1) = 3*ub(N+2:2*N + 2, 1);
ub(end) = Inf;

options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', alg,...
                       'MaxFunctionEvaluations', 200000, 'StepTolerance', 1e-15,...
                       'MaxIterations', 10,'ScaleProblem',true);
[x, fval] = fmincon(@objfun, x0, [], [], [], [], lb, ub, @state_const, options);

% visualzing the solution
T = x(2*N+2+1,1);
tvec = 0.5*(1-cheb_x)*T; 
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
fprintf('Time period (T) = %.2f s\n',T);

figure
subplot(1,2,1)
hold on
plot(tvec,x1);
plot(tvec,x2);
hold off
xlim([0,T])
legend('x_1','x_2');

subplot(1,2,2)
plot(x1,x2)
xlabel('x_1');
ylabel('x_2');
title('Phase Portrait');

% ------------------------------------------------ %
% ------------ required functions ---------------- %
% ------------------------------------------------ %
% to simulate system
function ydot = VDP(~, y)
global nu
%     nu = 1;
    y1 = y(1, 1); y2 = y(2, 1);
    y1dot = y2;
    y2dot = nu*(1 - y1^2)*y2 - y1;
    ydot = [y1dot;y2dot];
end

% state constraints
function [c, ceq] = state_const(x)
    N = (length(x)-3)/2;

    xk = x(1:N+1, 1); yk = x(N+2:2*N+2, 1); T = x(end);
    % get chebyshev matrix
    [D, ~] = cheb_diff(N);
    
    % get interpolated derivative
    diffmat = [D, zeros(N+1, N+1); zeros(N+1, N+1), D];
    xdot_cap = -(2/T)*diffmat*[xk;yk];
    
    % get actual derivative
    fx = zeros(N+1,1); fy = zeros(N+1,1);
    for i = 1:N+1
        temp = VDP(0, [xk(i);yk(i)]);
        fx(i, 1) = temp(1,1); fy(i, 1) = temp(2,1);
    end
    xdot = [fx;fy];
    
    c = [];
    ceq = xdot - xdot_cap;
end

function f = objfun(x)
    N = (length(x)-3)/2;
    scaler = 1;
    f = scaler*(x(1,1) - x(N+1, 1))^2 + (x(N+2, 1) - x(2*N+2, 1))^2;
end