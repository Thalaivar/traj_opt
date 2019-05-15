addpath('../lib')
addpath('../diagnostics')

clearvars

load('solutions\EE50.mat')
N = 50; p = 0.25;
% X0 = interpolateSol(sol, N);
% X0 = sol;
X0 = initGuess(N, 'eight');
% plot3(X0(1:N), -X0(N+1:2*N), -X0(2*N+1:3*N));

lb = [-500*ones(2*N,1); -Inf*ones(N,1);   -0.001;   0];
ub = [ 500*ones(2*N,1);     zeros(N,1); 100; 100]; 

[D, ~] = fourierdiff(N);
column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'UseParallel', true);
options.MaxIterations = 10000;
options.MaxFunctionEvaluations = 100000000;
options.StepTolerance = 1e-12;

sol = fmincon(@(X) objectiveFunctionDF(X), X0, [], [], [], [], lb, ub, @(X) constraintFunctionDF(X, p, D, DD), options);

rmpath('../lib')
rmpath('../diagnostics')

function X0 = initGuess(N, type)
    if strcmp(type, 'circle')
        VR_0 = 0.1; tf_0 = 10;
        [~, x] = fourierdiff(N);
        t = tf_0*x/(2*pi);
        x = 20*(-1 + cos(2*pi*t/tf_0));
        y = -20*(2^0.5)*sin(2*pi*t/tf_0);
        z = -20*(1 - cos(2*pi*t/tf_0)) - 0.2;    
    elseif strcmp(type, 'eight')
        VR_0 = 0.2; tf_0 = 20;
        [~, x] = fourierdiff(N);
        t = tf_0*x/(2*pi);
        x = 20*sin(4*pi*t/tf_0);
        y = 40*sqrt(2)*sin(2*pi*t/tf_0);
        z = 20*(sin(4*pi*t/tf_0) - 1) - 0.15; 
    end
    X0 = [x';y';z';tf_0;VR_0];
end

function X = interpolateSol(sol, N)
    n = (length(sol)-2)/3;
    T = sol(3*n+1);
    [~,t] = fourierdiff(n); t = T*t/(2*pi);
    [~,tt] = fourierdiff(N); tt = T*tt/(2*pi);
    X(1:N,1) = (interp_sinc(t, sol(1:n), tt))';
    X(N+1:2*N,1) = (interp_sinc(t, sol(n+1:2*n), tt))';
    X(2*N+1:3*N,1) = (interp_sinc(t, sol(2*n+1:3*n), tt))';
    X(3*N+1,1) = sol(3*n+1); X(3*N+2,1) = sol(3*n+2);
end