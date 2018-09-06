N = 100;
n_p = 4;

% simulate LV system
y0 = [0.9; 0.9]; % starting at 900 of prey and predator
tspan = [6.301, 22.099];   % time in years
[D, cheb_x] = cheb_diff(N);
cheb_t = ((tspan(2)-tspan(1))/2)*cheb_x + (sum(tspan))/2;

[t, y] = ode45(@LVS, tspan, y0);

% ode_samples = sample_ode(t, cheb_t, y);
ode_samples = interp1(t, y, cheb_t, 'pchip');
plot(cheb_t, ode_samples(:,1), '-r', cheb_t, ode_samples(:,2), '-b');

% optimisation algorithm
% alg = 'interior-point';
alg = 'sqp';

% setup NLP
x0 = zeros(2*N + 2 + n_p + 1, 1);
x0(1:N+1,1) = ode_samples(:,1); x0(N+2:2*N+2,1) = ode_samples(:,2); % sampled ode solution
x0(2*N + 2 + 1:2*N + 2 + n_p, 1) = [2/3; 4/3; 1; 1];    % parameters
x0(end, 1) = cheb_t(1) - cheb_t(end); % time period
lb = ones(2*N + 2 + n_p + 1, 1);
lb(1:N+1, 1) = 0.4*lb(1:N+1, 1); lb(N+2:2*N+2, 1) = 0.2*lb(N+2:2*N+2, 1);
lb(2*N+2+1:2*N+2+n_p, 1) = 0*lb(2*N+2+1:2*N+2+n_p, 1);
lb(end,1) = x0(end,1);
ub = ones(2*N + 2 + n_p + 1, 1);
ub(1:N+1, 1) = 1.8*ub(1:N+1, 1); ub(N+2:2*N+2, 1) = 1.0*ub(N+2:2*N+2, 1);
ub(2*N+2+1:end, 1) = Inf*ub(2*N+2+1:end, 1);

options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', alg, 'MaxFunctionEvaluations', 200000, 'StepTolerance', 1e-15);
[x, fval] = fmincon(@objfun, x0, [], [], [], [], lb, ub, @state_const, options);

% visualzing the solution
T = x(end);
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

% get ode samples (REDUNDANT)
% function samples = sample_ode(t, cheb_t, y)
%     % holds indices of time points in t closest to those in cheb_t
%     indices = zeros('like', cheb_t);
%     for i = 1:length(cheb_t)
%         curr_time = cheb_t(i);
%         temp = zeros('like', t);
%         for j = 1:length(t)
%             temp(j) = (t(j) - curr_time)^2;
%         end
%         [M, I] = min(temp);
%         indices(i) = I(1);
%     end
%     samples = zeros(length(indices),2);
%     for i = 1:length(indices)
%         samples(i,:) = y(indices(i), :);
%     end
% end

% for simulating the system
function ydot = LVS(~, y)
    y1 = y(1,1); y2 = y(2,1);
    p = [2/3, 4/3, 1, 1];
    ydot1 = p(1)*y1 - p(2)*y1*y2;
    ydot2 = p(3)*y1*y2 - p(4)*y2;
    ydot = [ydot1; ydot2];
end

% when using constants as dec vars
function ydot = LVS_opt(p, y)
    y1 = y(1,1); y2 = y(2,1);
    p = [2/3, 4/3, 1, 1];
    ydot1 = p(1)*y1 - p(2)*y1*y2;
    ydot2 = p(3)*y1*y2 - p(4)*y2;
    ydot = [ydot1; ydot2];
end

% state constraints
function [c, ceq] = state_const(x)
    N = (length(x)-7)/2;
    n_p = 4;
    xk = x(1:N+1, 1); yk = x(N+2:2*N+2, 1); T = x(2*N + 2 + 1, 1);
    p = x(2*N+2+1:2*N+2 + n_p, 1);
    % get chebyshev matrix
    [D, ~] = cheb_diff(N);
    
    % get interpolated derivative
    diffmat = [D, zeros(N+1, N+1); zeros(N+1, N+1), D];
    xdot_cap = (2/T)*diffmat*[xk;yk];
    
    % get actual derivative
    fx = zeros(N+1,1); fy = zeros(N+1,1);
    for i = 1:N+1
        temp = LVS_opt(p, [xk(i);yk(i)]);
        fx(i, 1) = temp(1,1); fy(i, 1) = temp(2,1);
    end
    xdot = [fx;fy];
    
    c = [];
    ceq = xdot - xdot_cap;
end

% cost function
function f = objfun(x)
    scaler = 10.0;
%     N = 100;
    N = (length(x)-7)/2;
    f = (x(1,1) - x(N+1, 1))^2 + (x(N+2, 1) - x(2*N+2, 1))^2;
end