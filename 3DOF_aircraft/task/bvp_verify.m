clearvars

addpath('..\trajectory')
addpath('..\floquet')
addpath('..')

global ac;
load('..\solutions\same_wind_shear\EOS1_01_50.mat')

FTM = get_FTM(ac, 'friedmann');
FTM = [FTM(1:3,1:3),FTM(1:3,6);FTM(6,1:3),FTM(6,6)];
[V, D] = eig(FTM);
D = log(D)/ac.tf;
eigval = zeros(1,4);
for i = 1:4
    eigval(i) = D(i,i);
end
[lam, i] = max(real(eigval));
v = V(:,i);
solinit = bvpinit(linspace(0, ac.tf, 1000), v, lam);
sol = bvp4c(@odefunc, @bcfun, solinit);

y0 = v;
tspan = linspace(0, ac.tf, 1000);
[t, y] = ode45(@(t,y) odefunc(t, y, lam), tspan, y0);

error = zeros(1, length(t));
for i = 1:length(t)
    error(i) = norm(sol.y(:,i) - y(i,:)');
end
plot(error)

function dydx = odefunc(t, y, lam)
    global ac;
    A = ac.get_jac(t, 'FD');
    A = [A(1:3,1:3),A(1:3,6);A(6,1:3),A(6,6)];
    dydx = (A - lam*eye(4))*y;
end

function res = bcfun(ya, yb, lam)
    res = ya - yb;
    res(end+1,1) = norm(ya) - 1;
end