function solStruct = stabilityOptimization(initialGuessSaveFile, windShear, objType)
    load(initialGuessSaveFile)

%     N = 500;
%     sol = interpolateSolution('diff-flat', sol, N);
    
    [D, fourierGrid] = fourierdiff(N);
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix

    trajParams.p = p; trajParams.t = fourierGrid;
    trajParams.D = D; trajParams.DD = DD;
    
    [~,freidmannGrid] = fourierdiff(500);
    trajParams.freidmannGrid = freidmannGrid;
    
    if strcmp(windShear, 'same')
        X0 = sol(1:3*N+1);
        trajParams.VR = sol(3*N+2);
    else
        X0 = sol;
    end
    
    [lb, ub] = optimbounds(N, windShear, sol(3*N+2));
%     options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'UseParallel', true);
%     options.MaxIterations = 10000;
%     options.MaxFunctionEvaluations = 100000000;
%     options.StepTolerance = 1e-12;
    options = optimoptions('ga', 'Display', 'iter', 'UseParallel', true, 'PopulationSize', 1000, 'MaxGenerations', 100000, 'ConstraintTolerance', 1e-5);
    sol = ga(@(X) optimizeStabilityDFD(X, trajParams, windShear, objType), length(X0), [], [], [], [], lb, ub, @(X) constraintFunctionDFD(X, trajParams, windShear), options);  
%     sol = fmincon(@(X) optimizeStabilityDFD(X, trajParams, windShear, objType), X0, [], [], [], [], lb, ub, @(X) constraintFunctionDFD(X, trajParams, windShear), options);
    
    if strcmp(windShear, 'same')
        sol = [sol; trajParams.VR];
    end
    solStruct.sol = sol; solStruct.N = N; solStruct.p = trajParams.p;
end

function [lb, ub] = optimbounds(N, windShear, VR0)
    if strcmp(windShear, 'same')
        lb = zeros(3*N+1,1); ub = zeros(3*N+1,1);
    else
        lb = zeros(3*N+2,1); ub = zeros(3*N+2,1);
    end
    lbounds = [-500, -500, -500];
    ubounds = [500, 500, -0.01];
    for i = 1:3
        j = (i-1)*N; 
        lb(j+1:j+N,1) = lbounds(i)*ones(N,1);
        ub(j+1:j+N,1) = ubounds(i)*ones(N,1);
    end
    lb(3*N+1,1) = 5; ub(3*N+1,1) = 150;
    if strcmp(windShear, 'not-same')
        lb(3*N+2,1) = 0.05; ub(3*N+2,1) = 1.2*VR0;
    end
end