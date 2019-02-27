%clearvars 

addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');
 
%load('solutions/trajectory_opt/lin_O.mat')

%[ac, sol] = optimize_stability(ac, [sol(1:end-2,1);sol(end,1)], p);
%xguess = [sol(1:end-2); sol(end,1)];
% [ac, sol] = optimize_stability(ac, sol, p);

M = [250, 500, 1000];
for i = 1:3
    [ac, sol] = optimize_stability(ac, sol, p, M(i));
    save(strcat('better_sol_', int2str(M(i))));
end

%ac = aircraft(); p = 0.25;
%ac.N = 5;
%p = 0.25;
%[ac, sol] = optimize_traj(ac, [], p);

rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');