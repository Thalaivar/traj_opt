clearvars

addpath('..\lib')

saveFile = 'solutions/trajectoryOptimized/E200.mat';
% saveFile = 'solutions/IG_circle_expo_50.mat';
% saveFile = 'solutions/sumObjFun/SE50_diffWind3.mat';
tic
solStruct = stabilityOptimization(saveFile, 'not-same', 'VR');
optTime = toc;

rmpath('..\lib')