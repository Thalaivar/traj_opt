addpath('solutions');
addpath('floquet');
addpath('trajectory');
addpath('task');
addpath('constraint_funcs');
 
clearvars
%load('solutions\trajectory_opt\EOTXX_YY_50.mat');
load('EOS10_01_50.mat')
ac.p = 0.25;
%ac.VR = ac.VR + 0.5;
[ac, sol] = optimize_stability(ac, [sol(1:end-2,1);sol(end,1)], ac.p, 50, [10,0.1], {'stability', 'stability'}, 'unstable'); 
sol = [sol(1:end-1,1);ac.VR;sol(end,1)];
save('EOS10_01_50.mat');

rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');