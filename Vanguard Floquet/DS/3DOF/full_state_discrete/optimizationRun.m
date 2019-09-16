clearvars
addpath('../../lib/')

% saveFile = 'solutions/trajectoryOptimized/VROpt_E70.mat';
% saveFile = 'solutions/IGs/IG_circle_expo_50.mat';
saveFile = 'VROpt_E60.mat';
% saveFile = 'temp.mat';
tic;
solStruct = stabilityOptimization(saveFile, 'circle', 'not-same');
optTime = toc;

sol = solStruct.sol;
N = solStruct.N;
p = solStruct.p;

rmpath('../../lib/')
