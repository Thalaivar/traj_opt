clearvars

addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');

%load('solutions/trajectory_opt/lin_O.mat');

ac = aircraft();
ac.N = 5; % no. of harmonics
p = 1; ac.p = p; % wind model param

[ac, sol] = optimize_traj(ac, [], p);
visualisation('traj-3d', ac);
axis equal

rmpath('trajectory')
rmpath('constraint_funcs')
rmpath('floquet')
