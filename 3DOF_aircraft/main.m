 %clearvars 

addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');
 
%load('solutions/trajectory_opt/lin_O.mat')

%[ac, sol] = optimize_stability(ac, [sol(1:end-2,1);sol(end,1)], p);
%xguess = [sol(1:end-2); sol(end,1)];
clearvars
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\trajectory_opt\LOT_XX_YY_50.mat')
ac.p = 1;
[ac, sol] = optimize_stability(ac, sol, p, 50, [1, 0.1]);
save('LOS1_01_50.mat');

clearvars
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\trajectory_opt\LOT_XX_YY_50.mat')
ac.p = 1;
[ac, sol] = optimize_stability(ac, sol, p, 50, [0.1, 0.1]);
save('LOS01_01_50.mat');

clearvars
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\trajectory_opt\LOT_XX_YY_50.mat')
ac.p = 1;
[ac, sol] = optimize_stability(ac, sol, p, 50, [10, 1]);
save('LOS10_1_50.mat');

clearvars
load('D:\ranjth_mohan\traj_opt\3DOF_aircraft\solutions\trajectory_opt\LOT_XX_YY_50.mat')
ac.p = 1;
[ac, sol] = optimize_stability(ac, sol, p, 50, [10, 10]);
save('LOS10_10_50.mat');

% xguess = [sol(1:end-2); sol(end,1)];
% p = ac.p;
% [ac, sol] = optimize_stability(ac, xguess, p, 250);

% M = [250, 500, 1000];
% for i = 1:3
%     [ac, sol] = optimize_stability(ac, sol, p, M(i));
%     save(strcat('better_sol_', int2str(M(i))));
% end


% ac.p = 0.25; p = 0.25;
% [ac, sol] = optimize_traj(ac, sol, p, 1000);

rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');