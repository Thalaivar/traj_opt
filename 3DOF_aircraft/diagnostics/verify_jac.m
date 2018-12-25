%clearvars

addpath('trajectory')
addpath('constraint_funcs')
addpath('floquet')
addpath('solutions')

%load('solutions/lin_O.mat')

J_cap = zeros(6);

t = linspace(0, ac.tf, 1000);
error = zeros(1,length(t)-1);

for k = 1:length(t)-1
    sigma = get_traj(t(k), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sigma);
    for i = 1:6
        x0 = [ac.x(1);ac.x(3);ac.x(2);sigma(1);sigma(2);sigma(3)]; I = eye(6);
        h = t(k+1) - t(k);
        f1 = ac.non_flat_model(t(k), x0-h*I(:,i)); f2 = ac.non_flat_model(t(k), x0+h*I(:,i));
        J_cap(:,i) = (f2 - f1)/(2*h);
    end
    A = ac.get_jac(t(k))
    J_cap(1:3,1:3)
   % error(k) = norm(J_cap - A)/norm(A);
end

rmpath('solutions')
rmpath('trajectory')
rmpath('constraint_funcs')
rmpath('floquet')

