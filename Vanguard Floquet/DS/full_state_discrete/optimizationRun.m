clearvars
addpath('../lib/')

saveFile = 'solutions\trajectoryOptimized\EE80.mat';

tic
solStruct = stabilityOptimization(saveFile, 'eight', 'not-same', 'stability');
optTime = toc;

sol = solStruct.sol;
N = solStruct.N;
p = solStruct.p;
save('solutions\maxObjFun\SEE80_diffWind.mat')
rmpath('../lib/')