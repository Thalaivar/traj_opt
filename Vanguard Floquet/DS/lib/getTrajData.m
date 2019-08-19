function trajData = getTrajData(sol, N, prm)
    if(strcmp(prm.model, '6dof')
        X = zeros(N,6);
    elseif(strcmp(prm.model, '3dof'))
        X = zeros(N, 6);
        U = zeros(N, 3);
        for i = 1:8
            j = (i-1)*N + 1;
            if i <= 6
                X(:,i) = sol(j:j+N-1,1);
            else
                U(:,i-6) = sol(j:j+N-1,1);
            end
        end
        T = sol(8*N+1); VR = sol(8*N+2);
        [D,fourierGrid] = fourierdiff(N); t = T*fourierGrid/(2*pi);
        if(prm.loiter), trajData.chiLinearTerm = sol(8*N+3); end
        if(prm.correctChi)
            for i = 1:N
                X(i,2) = X(i,2) + trajData.chiLinearTerm*t(i);
            end
        end
        trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
        trajData.X = X; trajData.U =  U; trajData.N = N; trajData.p = prm.p; trajData.D = D;
    elseif(strcmp(prm.model, 'vanderpol'))
        X = zeros(N,2);
        X(:,1) = sol(1:N,1);
        X(:,2) = sol(N+1:2*N,1);
        T = sol(2*N+1);
        [D, fourierGrid] = fourierdiff(N);
        trajData.X = X; trajData.N = N; trajData.T = T;
        trajData.D = D; trajData.fourierGrid = fourierGrid;
        trajData.mu = prm.mu;
    elseif(strcmp(prm.model, 'glycolytic'))
        X = zeros(N,2);
        X(:,1) = sol(1:N,1);
        X(:,2) = sol(N+1:2*N,1);
        T = sol(2*N+1);
        [D, fourierGrid] = fourierdiff(N);
        trajData.X = X; trajData.N = N; trajData.T = T;
        trajData.D = D; trajData.fourierGrid = fourierGrid;
        trajData.a = prm.a; trajData.b = prm.b;
    end
end