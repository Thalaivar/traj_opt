% clear all
addpath('..\lib\')

ac = create_acmod();
% load('../full_state_discrete/solutions/trajectoryOptimized/E50.mat')
% sol = inclinedCircleGuess(50, 'loiter', ac);
% sol = convert3DoF(sol, N, 'loiter', ac);
tic
solStruct = stabilityOptimization(sol, 'loiter', 'not-same', 'VR');
optTime = toc;

rmpath('..\lib\')