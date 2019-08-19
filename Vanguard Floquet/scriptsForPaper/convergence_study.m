clearvars
addpath('../DS/lib/')

load('E500.mat')
d = 4;
trajData.shape = 'circle'; trajData.type = 'full-state';
trajData = stateControlMat(sol, N, p, trajData);

[~, timeMarchFE] = timeMarchMethod(trajData);
timeMarchFE = sortrows(timeMarchFE, 'descend');

%%
N = linspace(40, 500, 47);
FEData = zeros(d, length(N));
discpData = zeros(d, length(N));
for i = 1:length(N)
    [~, FE] = spectralMethodMod(trajData, N(i));
    FE = sortrows(FE, 'descend');
    for j = 1:d
        discpData(j,i) = abs(FE(j) - timeMarchFE(j))/abs(timeMarchFE(j));
    end
    FEData(:,i) = FE';
end

function [eigE, FE] = spectralMethodMod(trajData, N)
    d = 4;
    T = trajData.T;
    [D,t] = fourierdiff(N); t = t*T/(2*pi);
    Dmat = zeros(d*N); Mmat = zeros(d*N);
    
    tt = trajData.T*trajData.fourierGrid/(2*pi);
    X = zeros(N,6); U = zeros(N,3);
    for i = 1:6
        X(:,i) = interp_sinc(tt, trajData.X(:,i), t);
    end
    for i = 1:N
        X(i,2) = X(i,2) + trajData.chiLinearTerm*t(i);
    end

    U(:,1) = interp_sinc(tt, trajData.U(:,1), t);
    U(:,2) = interp_sinc(tt, trajData.U(:,2), t);
    
    for i = 1:N
        A = sysModel(trajData.p, X(i,:), [U(i,:), trajData.VR]);
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
           end
           Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
        end
    end
    
    eigE = eig(Dmat - Mmat);
    eigE = -1*eigE;
    FE = identifyFloquet(eigE, N, trajData.T);
end
function A = sysModel(p, Z, U)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = p;
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
%     A = A(1:3,1:3);
end
function trajData = stateControlMat(sol, N, p, trajData)
    x = zeros(N, 6);
    u = zeros(N, 3);
    for i = 1:8
        j = (i-1)*N + 1;
        if i <= 6
            x(:,i) = sol(j:j+N-1,1);
        else
            u(:,i-6) = sol(j:j+N-1,1);
        end
    end
    T = sol(8*N+1,1); VR = sol(8*N+2,1);
    [D,fourierGrid] = fourierdiff(N); t = T*fourierGrid/(2*pi);
    if strcmp(trajData.shape, 'circle')
        trajData.chiLinearTerm = sol(8*N+3);
%         for i = 1:N
%             x(i,2) = x(i,2) + chiLinearTerm*t(i);
%         end
    end
    trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
    trajData.X = x; trajData.U =  u; trajData.N = N; trajData.p = p; trajData.D = D;
end