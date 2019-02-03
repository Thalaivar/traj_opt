%clearvars 

addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');
 
%load('solutions/trajectory_opt/expo_O_shaped.mat')

%[ac, sol] = optimize_stability(ac, [sol(1:end-2,1);sol(end,1)], p);
[ac, sol] = optimize_stability(ac, sol, p);

% options = optimoptions('fminunc', 'Display', 'iter-detailed', 'StepTolerance', 1e-15);
% sol = fminunc(@(X) objfun_eig(X, ac), ac.tf, options);

rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');