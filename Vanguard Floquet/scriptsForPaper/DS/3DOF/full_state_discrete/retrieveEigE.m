clearvars

addpath('../../lib')

%load('solutions/trajectoryOptimized/E50.mat')
load('solutions/finalResults/diffWind/LS_1_AM_GS_O_S_E.mat')

[D, fGrid] = fourierdiff(N) ;

[X,U,T,VR] = stateControlMat3DOF(sol,N,'not-same');    
chiLinearTerm = nan;
if length(sol) == 8*N+3, chiLinearTerm = sol(8*N+3); end

if ~isnan(chiLinearTerm)
    t = T*fGrid/(2*pi);
    for i = 1:N
        X(i,2) = X(i,2) + chiLinearTerm*t(i);
    end
end

trajData.X = X; trajData.U = U; trajData.N = N; 
trajData.VR = VR; trajData.T = T; trajData.D = D;
trajData.p = p; trajData.fourierGrid = fGrid;
trajData.type = 'full-state';

[spectralFE, eigE, ~, ~, eigVec] = spectralMethod(trajData);

% [FTM, freidmannFE] = freidmannMethod(trajData, fGrid);
% 
% [~, tmFE] = timeMarchMethod(trajData);

rmpath('../../lib')