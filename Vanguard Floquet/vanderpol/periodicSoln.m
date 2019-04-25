clearvars
% retrieve periodic solution for Vanderpol oscillator
mu = -4;
% choice = 'levenberg';
choice = 'fourier';
if strcmp(choice, 'fourier')
    N = 300; d = 2;
%     X0 = init_guess(N);
    load('mu_m4.mat')
    mu = -4;
%     NN = 0.5*(length(X)-1);
%     [~,x] = fourierdiff(NN);
%     tt = x*X(2*NN+1,1)/(2*pi);
%     [~,x] = fourierdiff(N);
%     t = x*X(2*NN+1,1)/(2*pi);
%     X0 = [interp_sinc(tt, X(1:NN,1), t)'; interp_sinc(tt, X(1+NN:2*NN,1), t)'];
%     X0(end+1,1) = X(2*NN+1,1);
    X0 = X;
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000000, 'MaxIterations', 10000);
    A = zeros(3,2*N+1); A(1,1) = -1; A(2,N+1) = -1; A(3,2*N+1) = 1;
    b = [-1; -1; 25];
    X = fmincon(@(X) objfun(X, mu), X0, A, b, [], [], [], [], [], options);
    plot(X(1:N,1), X(N+1:2*N,1), 'b', 'LineWidth', 1.25);
    xlabel('$x$', 'Interpreter', 'latex');
    xlabel('$\dot{x}$', 'Interpreter', 'latex');
    title(strcat('Time histories for $\mu =$ ', num2str(mu)), 'Interpreter', 'latex');
elseif strcmp(choice, 'levenberg')
    X0 = [2; 2; 15];
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'MaxFunctionEvaluations', 1000000, 'MaxIterations', 10000);
    X = lsqnonlin(@(X) resvec(X,mu), X0, [], [], options);
    X0 = [X(1,1); X(2,1)];
    tspan = [0, X(3,1)];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    [t,y] = ode45(@(t,y) vanderpol(t, y, mu), tspan, X0, options);
    plot(t,y(:,1),'r','LineWidth', 1.25)
    hold on
    plot(t,y(:,2),'b','LineWidth', 1.25)
    legend('x', 'xdot')
    title(strcat('Time histories for $\mu =$ ', num2str(mu)), 'Interpreter', 'latex');
end
%% functions
function X0 = init_guess(N)
    T0 = 10;
    [D,x] = fourierdiff(N);
    t = x*T0/(2*pi);
    x1 = 2.5*(sin(2*pi*t/T0))';
    x2 = 2*pi*D*x1/T0;
    X0 = [x1;x2;T0];
end

function dy = vanderpol(t, y, mu)
    dy(1,1) = y(2);
    dy(2,1) = mu*(1-(y(1))^2)*y(2) - y(1);
end

function f = objfun(X, mu)
    N = (length(X)-1)/2;
    [D,~] = fourierdiff(N);
    x1 = X(1:N,1); x2 = X(N+1:2*N); T = X(2*N+1,1);
    xdot_cap = [2*pi*D*x1/T; 2*pi*D*x2/T];
    xdot = zeros(size(xdot_cap));
    for i = 1:N
        dy = vanderpol([], [x1(i); x2(i)], mu);
        xdot(i,1) = dy(1);
        xdot(i+N,1) = dy(2);
    end
    f = norm(xdot_cap - xdot)/(2*N);
end

function f = resvec(X, mu)
    X0 = [X(1,1); X(2,1)];
    tspan = [0, X(3,1)];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    [~,y] = ode45(@(t,y) vanderpol(t, y, mu), tspan, X0, options);
    f = y(1,:)' - y(end,:)';
end