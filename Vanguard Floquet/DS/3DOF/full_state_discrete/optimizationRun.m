clearvars
addpath('../../lib/')

% saveFile = 'solutions\trajectoryOptimized\E50.mat';
saveFile = 'lineSep1.mat';

tic
solStruct = stabilityOptimization(saveFile, 'circle', 'not-same', 'stability');
optTime = toc;

% sol = solStruct.sol;
% N = solStruct.N;
% p = solStruct.p;
% save('solutions\maxObjFun\SEE80_diffWind.mat')
rmpath('../lib/')