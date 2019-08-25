function solStruct = stabilityOptimization(initialGuessSaveFile, type, windShear)
    load(initialGuessSaveFile)
%     p = 0.25;
    
%     N = 80; 
%     sol = interpolateSolution('full-state', sol, N, type);
%     
    [D, fourierGrid] = fourierdiff(N);
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
    
    trajData.p = p; trajData.fourierGrid = fourierGrid;
    trajData.D = D; trajData.DD = DD; trajData.N = N;
    
%     freidmannGridN = 500;
%     [~,freidmannGrid] = fourierdiff(freidmannGridN);
%     trajData.freidmannGrid = freidmannGrid;
    
    if strcmp(windShear, 'same')
        if strcmp(type, 'circle')
            X0 = [sol(1:8*N+1); sol(8*N+3)];
        else
            X0 = sol(1:8*N+1);
        end
        trajData.VR = sol(8*N+2);
    else
        X0 = sol;
    end
    
    [lb, ub] = optimbounds(N, type, windShear, sol(8*N+2));
    
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'UseParallel', true);
    options.MaxIterations = 10000;
    options.MaxFunctionEvaluations = 100000000;
%     options.StepTolerance = 1e-12;
    
    sol = fmincon(@(X) optimizeStabilityFSD(X, trajData, windShear), X0, [], [], [], [], lb, ub, @(X) constrainFunctionFSD(X, trajData, windShear), options);
    
    if strcmp(windShear, 'same')
        if strcmp(type, 'circle')
            sol = [sol(1:8*N+1); trajData.VR; sol(8*N+2)];
        else
            sol = [sol; trajData.VR];
        end
    end
    solStruct.sol = sol; solStruct.N = N; solStruct.p = trajData.p;
%     solStruct.shape = trajData.shape;
end

function [lb, ub] = optimbounds(N, type, windShear, VR0)
    if strcmp(type, 'circle')
        if strcmp(windShear, 'same')
            lb = zeros(8*N+2,1); ub = zeros(8*N+2,1);
        else
            lb = zeros(8*N+3,1); ub = zeros(8*N+3,1);
        end
    else
        if strcmp(windShear, 'same')
            lb = zeros(8*N+1,1); ub = zeros(8*N+1,1);
        else
            lb = zeros(8*N+2,1); ub = zeros(8*N+2,1);
        end
    end
    lbounds = [10, -2*pi, -pi/4, -200, -200, -200,  -0.2, -pi/3];
    ubounds = [80,  2*pi,  pi/4,  200,  200,-0.01,  1.17,  pi/3];
    for i = 1:8
        j = (i-1)*N; 
        lb(j+1:j+N,1) = lbounds(i)*ones(N,1);
        ub(j+1:j+N,1) = ubounds(i)*ones(N,1);
    end
    lb(8*N+1,1) = 0; ub(8*N+1,1) = 150;
    if strcmp(windShear, 'same')
        if strcmp(type, 'circle')
            lb(8*N+2) = -Inf; ub(8*N+2) = Inf;
        end
    else
        lb(8*N+2) = 0; 
        ub(8*N+2) = 1.2*VR0;
%         ub(8*N+2) = 100;
        if strcmp(type, 'circle')
            lb(8*N+3) = -Inf; ub(8*N+3) = Inf;
        end
    end
end