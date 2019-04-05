addpath('..\trajectory')
addpath('..\floquet')
addpath('..')
global ac;

load('..\solutions\same_wind_shear\EOS1_01_50.mat')

% initial guess
lam = 0;
FTM = get_FTM(ac, 'friedmann');
[V, D] = eig(FTM);
D = log([D(3,3), D(4,4), D(5,5), D(6,6)])/ac.tf;
[~,i] = max(real(D));
global lam;
lam = D(i);
V = V(:,3:end);
v = V(:,i);
solinit = bvpinit(linspace(0,ac.tf,1000), v);
sol = bvp4c(@odefunc, @bfunc, solinit);

rmpath('..\trajectory')
rmpath('..')

function dydx = odefunc(t, y)
    global ac;
    global lam;
    A = ac.get_jac(t, 'FD');
    dydx = (A - lam*eye(6))*y;
end

function res = bfunc(ya, yb)
    res = ya - yb;
end