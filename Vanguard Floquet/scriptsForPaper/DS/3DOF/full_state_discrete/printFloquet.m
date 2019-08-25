addpath('../../lib/')
% load('solutions\trajectoryOptimized\EE70.mat')
load('resultsForPaper\diffWind_EE70.mat')

[D, fGrid] = fourierdiff(N);
[x,u,T,VR] = stateControlMat3DOF(sol,N,'not-same');

chiLinearTerm = nan;
if length(sol) == 8*N+3, chiLinearTerm = sol(8*N+3); end

if ~isnan(chiLinearTerm)
    t = T*fGrid/(2*pi);
    for i = 1:N
        x(i,2) = x(i,2) + chiLinearTerm*t(i);
    end
end

trajData.N = N; trajData.X = x; trajData.U = u;
trajData.T = T; trajData.VR = VR; trajData.D = D;
trajData.p = p;

[FE, eigE, AM, groupSizes, eigVec] = spectralMethod(trajData);
display(max(real(FE)));

rmpath('../../lib/')