clearvars

load('../3DOF/full_state_discrete/solutions/trajectoryOptimized/E500.mat')
trajData.shape = 'circle'; trajData.type = 'full-state';
trajData = stateControlMat(sol, N, p, trajData);

% get time marched estimates for truth
[~,timeMarchFE] = timeMarchMethod(trajData);

display("Time marching successful!")

% FE = spectralMethod(trajData);
% val = maxDiscrepancy(timeMarchFE, FE);
%%
[eigE,FE] = spectralMethodMod(trajData, 300);
 
% N = linspace(50, 500, 10); 
% maxDiscp = zeros(1, 10);
% for i = 1:length(N)
%     [~,FE] = spectralMethodMod(trajData, N(i));
%     maxDiscp(i) = maxDiscrepancy(timeMarchFE, FE);
% end 

% plot(N, maxDiscp, '--om')

function val = maxDiscrepancy(timeMarchFE, FE)
    timeMarchFE = sort(timeMarchFE, 'ComparisonMethod', 'real');
    FE = sort(FE, 'ComparisonMethod', 'real');
    i = 1; val = 0;
    while i <= length(FE)
        if imag(FE(i)) ~= 0
            cmplx1 = complex(real(timeMarchFE(i)), abs(imag(timeMarchFE(i))));
            cmplx2 = complex(real(FE(i)), abs(imag(FE(i))));
            currDiscp = abs(cmplx1 - cmplx2)/abs(cmplx1);
            val = max(val, currDiscp);
            i = i + 2;
        else
            currDiscp = abs(FE(i)-timeMarchFE(i))/abs(timeMarchFE(i));
            val = max(val, currDiscp);
            i = i + 1;
        end
    end
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
end
function trajData = stateControlMat(sol, N, p, trajData)
    x = zeros(N, 6);
    u = zeros(N, 3);
    if strcmp(trajData.type, 'full-state')
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
            chiLinearTerm = sol(8*N+3);
            for i = 1:N
                x(i,2) = x(i,2) + chiLinearTerm*t(i);
            end
        end
    elseif strcmp(trajData.type, 'diff-flat')
        T = sol(3*N+1); VR = sol(3*N+2);
        diffFac = 2*pi/T;

        [D,fourierGrid] = fourierdiff(N);
        column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
        DD = toeplitz(column2,column2([1, N:-1:2])); 
        x = [sol(1:N), sol(N+1:2*N), sol(2*N+1:3*N)];
        dx = [diffFac*D*x(:,1), diffFac*D*x(:,2), diffFac*D*x(:,3)];
        ddx = [(diffFac^2)*DD*x(:,1), (diffFac^2)*DD*x(:,2), (diffFac^2)*DD*x(:,3)];

        state = zeros(N,3); u = zeros(N,3);
        for i = 1:N
            [xx, uu] = diffFlatModel([x(i,:);dx(i,:);ddx(i,:)], p, VR);
            state(i,:) = xx; u(i,:) = uu;
        end
        x = [state, x];
    end
    trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
    trajData.X = x; trajData.U =  u; trajData.N = N; trajData.p = p; trajData.D = D;
end