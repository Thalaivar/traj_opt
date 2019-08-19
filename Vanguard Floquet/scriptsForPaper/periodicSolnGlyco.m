clearvars

a = 0.08; B = 0.6;
N = 300; d = 2;

% X0 = init_guess(N);

% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000000, 'MaxIterations', 10000, 'UseParallel', true);
% 
% % A = zeros(3, 2*N+1);
% % A(1,1) = -1; 
% % A(2,N+1) = -1; 
% % A(3,2*N+1) = 1;
% % b = [-1; -1; 5];
% 
% X = fmincon(@(X) objfun(X, a, B), X0, [], [], [], [], [], [], [], options);
% plot(X(1:N,1), X(N+1:2*N,1), 'b', 'LineWidth', 1.25);
% xlabel('$x$', 'Interpreter', 'latex');
% xlabel('$y$', 'Interpreter', 'latex');

X0 = [0.75; 0.5; 10];
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'MaxFunctionEvaluations', 1000000, 'MaxIterations', 10000);
X = lsqnonlin(@(X) resvec(X, a, B), X0, [], [], options);
X0 = [X(1,1); X(2,1)];
N = 300;
[~,tspan] = fourierdiff(N);
tspan = tspan*X(3,1)/(2*pi);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t,y] = ode45(@(t,y) glycolytic(t, y, a, B), tspan, X0, options);
plot(y(:,1), y(:,2))

function X0 = init_guess(N)
    T0 = 1;
    [~,x] = fourierdiff(N);
    t = x*T0/(2*pi);
    x1 = -2.5*(sin(2*pi*t/T0))';
    x2 = x1;
    X0 = [x1;x2;T0];
end

function f = objfun(X, a, b)
    N = (length(X)-1)/2;
    [D,~] = fourierdiff(N);
    x1 = X(1:N,1); x2 = X(N+1:2*N); T = X(2*N+1,1);
    xdot_cap = [2*pi*D*x1/T; 2*pi*D*x2/T];
    xdot = zeros(size(xdot_cap));
    for i = 1:N
        dy = glycolytic([x1(i), x2(i)], a, b);
        xdot(i,1) = dy(1);
        xdot(i+N,1) = dy(2);
    end
    f = norm(xdot_cap - xdot);
end

function dy = glycolytic(t, Y, a, b)
    x = Y(1); y = Y(2);
    dy(1,1) = -x + a*y + y*(x^2);
    dy(2,1) = b - a*y - (x^2)*y;
end

function f = resvec(X, a, b)
    X0 = [X(1,1); X(2,1)];
    tspan = [0, X(3,1)];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    [~,y] = ode45(@(t,y) glycolytic(t, y, a, b), tspan, X0, options);
    f = y(1,:)' - y(end,:)';
end