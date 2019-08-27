clearvars
addpath('../../lib/')

% saveFile = 'solutions/trajectoryOptimized/VROpt_E70.mat';
% saveFile = 'solutions/IGs/IG_circle_expo_50.mat';

tic;
solStruct = stabilityOptimization('VROpt_E50.mat', 'circle', 'not-same');
optTime = toc;

sol = solStruct.sol;
N = solStruct.N;
p = solStruct.p;

save("LS1_sameWind_E50.mat", 'sol', 'N', 'p', 'optTime');

clearvars

saveFile = 'solutions/trajectoryOptimized/E50.mat';

tic;
solStruct = stabilityOptimization(saveFile, 'circle', 'not-same');
optTime = toc;

sol = solStruct.sol;
N = solStruct.N;
p = solStruct.p;

save("LS1_diffWind_E50.mat", 'sol', 'N', 'p', 'optTime');

rmpath('../../lib/')
