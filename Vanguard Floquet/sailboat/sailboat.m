N = 100;
X0 = initguess(N);
% X0 = xguess;
b = 0.18;
% [c , ceq] = constFun(X0);
lb = -1*ones(3*N+1,1); ub = 1*ones(3*N+1,1);
lb(end,1) = 0; ub(end,1) = 50;
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'UseParallel', true);
options.MaxFunctionEvaluations = 10000000;
options.MaxIterations = 10000;
sol = fmincon(@(X) objfun(X, b), X0, [], [], [], [], lb, ub, @constFun, options);
plot(sol(1:N),'r')
hold on
plot(sol(1+N:2*N),'b')
plot(sol(1+2*N:3*N),'g')
legend('x1', 'x2', 'u')

[~,chebx] = chebdiff(N-1);
chebt = 0.5*sol(3*N+1,1)*(1-chebx);
save('sailboat_data.mat', 'sol', 'N', 'chebt');

function X0 = initguess(N)
    [D, x] = chebdiff(N-1);
    T = 10;
    t = (0.5*T*(1-x))';
    x1 = 0.09*sin(2*pi*t/T);
    x2 = -2*D*x1/T;
    u = -2*D*x2/T;
    X0 = [x1;x2;u;T];
end

function f = objfun(X, b)
    N = (length(X)-1)/3;
    X_th = [X(1:N,1), X(N+1:2*N,1), X(2*N+1:3*N,1)];
    [~,w] = clencurt(N-1);
    g = zeros(N,1);
    for i = 1:N
        g(i,1) = 0.25*(X_th(i,1)^2 + X_th(i,2)^4) - 0.5*(X_th(i,2)^2 - b*X_th(i,3)^2);
    end
    f = (w*g)/2;
end

function [c, ceq] = constFun(X)
    N = (length(X)-1)/3;
    X_th = [X(1:N,1), X(N+1:2*N,1), X(2*N+1:3*N,1)];
    T = X(3*N+1,1);
    Xdot_cap = zeros(2*N,1);
    [D,~] = chebdiff(N-1);
    % get spectral derivatives
    for i = 1:2
        j = (i-1)*N+1;
        Xdot_cap(j:j+N-1,1) = -(2/T)*D*X_th(:,i); 
    end
    % get actual derivatives
    Xdot = zeros(2*N,1);
    for i = 1:N
        xdot = sailboat_model([], [X_th(i,1); X_th(i,2)], X_th(i,3));
        Xdot(i,1) = xdot(1); Xdot(i+N,1) = xdot(2);
    end
    
    ceq(1:2*N,1) = Xdot - Xdot_cap;
    ceq(2*N+1,1) = X_th(1,1) - X_th(N,1);
    ceq(2*N+2,1) = X_th(1,2) - X_th(N,2);
    ceq(2*N+3,1) = X_th(1,3) - X_th(N,3);
    c = [];
end

function dy = sailboat_model(t, y, u)
    dy = [y(2); u];
end