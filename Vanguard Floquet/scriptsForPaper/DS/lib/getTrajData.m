function trajData = getTrajData(sol, N, prm)
    if(strcmp(prm.model, '6dof'))
        X = zeros(N,12); 
        U = zeros(N,4);
        for i = 1:16
            j = (i-1)*N + 1;
            if i <= 12
                X(:,i) = sol(j:j+N-1,1);
            else
                U(:,i-12) = sol(j:j+N-1,1);
            end
        end
        T = sol(16*N+1); VR = sol(16*N+2);
        [D,fourierGrid] = fourierdiff(N); t = T*fourierGrid/(2*pi);
        if(prm.loiter), trajData.chiLinearTerm = sol(16*N+3); end
        if(prm.correctChi)
            for i = 1:N
                X(i,9) = X(i,9) + trajData.chiLinearTerm*t(i);
            end
        end
        trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
        trajData.X = X; trajData.U =  U; trajData.N = N; trajData.p = prm.p; trajData.D = D;
        load('rawMaterial\acmod.mat')
        trajData.ac = ac;
        
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
        column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
        DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
        if(prm.loiter), trajData.chiLinearTerm = sol(8*N+3); end
        if(prm.correctChi)
            for i = 1:N
                X(i,2) = X(i,2) + trajData.chiLinearTerm*t(i);
            end
        end
        trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
        trajData.X = X; trajData.U =  U; trajData.N = N; trajData.p = prm.p; trajData.D = D;
        trajData.DD = DD;
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