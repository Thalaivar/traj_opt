clearvars
addpath('floquet');
addpath('trajectory');
load('solutions\same_wind_shear\EOS01_01_50.mat');

N = 10 0;

[D, cheb_x] = cheb_diff(N-1);
t = 0.5*ac.tf*(1-cheb_x);

FTM_expo = get_FTM(ac, 'friedmann');
FTM_expo = [FTM_expo(1:3,1:3),FTM_expo(1:3,6);FTM_expo(6,1:3),FTM_expo(6,6)];
[V, eig] = eig(FTM_expo);
eigvals = zeros(4,1);
for i = 1:4
    eigvals(i,1) = eig(i,i);
end
eigvals = log(eigvals)/ac.tf;
[~,max_eig] = max(real(eigvals));
v = complex(real(V(:,max_eig)), imag(V(:,max_eig))); 
lam = complex(real(eigvals(max_eig)), imag(eigvals(max_eig)));
[t, u] = ode45(@(t, u) odefunc(t, u, ac, lam), t, [real(v);imag(v)]);

X0 = zeros(8*N + 2,1);
for i = 1:4
    j = (i-1)*N + 1;
    X0(j:j+(N-1),1) = u(:,i);
    X0(4*N+j:4*N+j+(N-1),1) = u(:,i+4);
end

X0(end-1,1) = real(lam);
X0(end,1) = imag(lam);

[~, ceq] = constFun_colloc(X0, ac);

lb = -10*ones(8*N + 2,1); ub = 10*ones(8*N + 2,1);
lb(end-1:end,1) = [-0.1; -1]; ub(end-1:end,1) = [Inf; 1];
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'UseParallel', true);
options.MaxFunctionEvaluations = 100000;
options.MaxIterations = 1000;
options.StepTolerance = 1e-12;
x = fmincon(@(x) objfun_colloc(x), X0, [], [], [], [], lb, ub, @(x) constFun_colloc(x, ac), options);

rmpath('floquet');
rmpath('trajectory');

function udot = odefunc(t, u, ac, eigval)
    u_re = u(1:4,1); u_im = u(5:8,1);
    A = ac.get_jac(t, 'FD');
    A = [A(1:3,1:3),A(1:3,6);A(6,1:3),A(6,6)];
    udot(1:4,:) = (A - real(eigval)*eye(4))*u_re + imag(eigval)*u_im;
    udot(5:8,:) = (A - real(eigval)*eye(4))*u_im - imag(eigval)*u_re;
end