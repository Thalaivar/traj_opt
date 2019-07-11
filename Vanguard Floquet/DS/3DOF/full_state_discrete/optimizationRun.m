clearvars
addpath('../../lib/')

saveFile = 'solutions/trajectoryOptimized/L50.mat';
% saveFile = 'solutions\finalResults\diffWind\LS_100_AM_GS_E_S_L.mat'; 
% saveFile = 'solutions\IGs\IG_eight_linear_70.mat';

tic;
solStruct = stabilityOptimization(saveFile, 'circle', 'same');
optTime = toc;

sol = solStruct.sol;
N = solStruct.N;
p = solStruct.p;

save("LS_100_O_L.mat", 'sol', 'N', 'p');

rmpath('../../lib/')
