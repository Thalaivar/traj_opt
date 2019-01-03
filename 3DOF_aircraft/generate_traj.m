clc

addpath('trajectory')
addpath('constraint_funcs')
addpath('floquet')

ac = aircraft();
ac.N = 8; % no. of harmonics
p = 1; % wind model param

[ac, sol] = optimize_traj(ac, [], p);
visualisation('traj-3d', ac);
axis equal

rmpath('trajectory')
rmpath('constraint_funcs')
rmpath('floquet')
