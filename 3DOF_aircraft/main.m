addpath('solutions');
addpath('floquet');
addpath('trajectory');  
addpath('task');
addpath('constraint_funcs');

% clearvars
% load('solutions\trajectory_opt\EOTXX_YY_50.mat')
% load('EOS01_01_50.mat')

sol = 
ac.p = 0.25;
ac.VR = ac.VR + 0.5;

tic
[ac, sol] = optimize_stability(ac, X0, ac.p, 50, [1, 0.1], {'stability', 'stability'}, 'stable'); 
optTime = toc;

sol = [sol(1:end-1,1);ac.VR;sol(end,1)];
% save('solutions/max_E.mat')

rmpath('solutions');
rmpath('floquet');
rmpath('trajectory');
rmpath('task');
rmpath('constraint_funcs');