clearvars

addpath('../trajectory')
addpath('../constraint_funcs')
addpath('../floquet')
addpath('../solutions')
addpath('..')

load('solutions/lin_O.mat')

t = linspace(0, ac.tf, 1000);
error = zeros(1,length(t)-1);

for k = 1:length(t)
    A_FD = ac.get_jac(t(k), 'FD')
    A_AN = ac.get_jac(t(k), 'analytic')
    error(k) = norm(A_FD(1:3,1:3) - A_AN(1:3,1:3))/norm(A_FD(1:3,1:3));
end

rmpath('..')
rmpath('../solutions')
rmpath('../trajectory')
rmpath('../constraint_funcs')
rmpath('../floquet')


