function solStruct = stabilityOptimization(sol, shape, windShear, objType)
    nState = 12; nControl = 4; nOptim = nState + nControl;
    ac = create_acmod();
    trajData.shape = shape;
    trajData.type = windShear;
    
    if strcmp(shape, 'loiter'), N = (length(sol)-3)/nOptim; else, N = (length(sol)-2)/nOptim; end
    
    if strcmp(windShear, 'same')
        if strcmp(shape, 'loiter'), X0 = [sol(1:nOptim*N+1);sol(nOptim*N+3)]; else, X0 = sol(1:nOptim*N+1); end
        trajData.VR = sol(nOptim*N+2);
    else
        X0 = sol;
    end
    trajData.T = sol(nOptim*N+1);
    
    [D, fourierGrid] = fourierdiff(N);
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
    
    trajData.p = ac.p_exp; trajData.D = D; trajData.DD = DD; trajData.fourierGrid = fourierGrid;
    trajData.ac = ac; trajData.N = N;
    
    [lb, ub] = getBounds(nOptim, N, shape, windShear, sol(nOptim*N+2));
    options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', true, 'Algorithm', 'sqp');
    options.MaxIterations = 10000;
    options.MaxFunctionEvaluations = 100000000;
    
    sol = fmincon(@(X) objectiveFunc(X,trajData,objType),X0,[],[],[],[],lb,ub,@(X) constraintFunc(X, trajData), options);
    
    if strcmp(windShear, 'same')
        if strcmp(shape, 'loiter')
            sol = [sol(1:nOptim*N+1);trajData.VR;sol(nOptim*N+2)];
        else
            sol = [sol;trajData.VR];
        end
    end
    
    solStruct.sol = sol; solStruct.N = trajData.N; solStruct.p = ac.p_exp;
end

function [lb, ub] = getBounds(nOptim, N, shape, windShear, VR0)
    %         [  u,    v,    w,    p,    q,    r,   phi,  theta,  psi,    x,    y,    z,    df,    da,    de,    dr ]
    lbounds = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -pi/3, -pi/3, -2*pi, -200, -200, -200, -pi/4, -pi/4, -pi/4, -pi/4];
    ubounds = [ Inf,  Inf,  Inf,  Inf,  Inf,  Inf,  pi/3,  pi/3,  2*pi,  200,  200,-0.01,  pi/4,  pi/4,  pi/4,  pi/4];
    lb = zeros(nOptim*N,1); ub = zeros(nOptim*N,1);
    
    if nOptim ~= length(lbounds), error("Insufficient number of bounds."); end
    
    for i = 1:length(lbounds)
        j = (i-1)*N;
        lb(j+1:j+N) = lbounds(i)*ones(N,1);
        ub(j+1:j+N) = ubounds(i)*ones(N,1);
    end
    % bound on T
    lb(end+1) = 0; ub(end+1) = 150;
    if strcmp(windShear, 'not-same')
        % bound on VR
        lb(end+1) = 0; ub(end+1) = 100;
    end
    if strcmp(shape, 'loiter')
        % bound on psi linear term
        lb(end+1) = -Inf; ub(end+1) = Inf;
    end
end