function [FE, eigE, eigVec] = spectralMethodGen(trajData, jacEval, d)
    % spectral method runs on same number of grid points as trajectory
    D = trajData.D; N = trajData.N; T = trajData.T;
    
    Dmat = zeros(d*N); Mmat = zeros(d*N);
    for i = 1:N
        if(isfield(trajData, 'U'))
            A = jacEval(trajData.X(i,:), trajData.U(i,:), trajData);
        else
            A = jacEval(trajData.X(i,:), trajData);
        end
        
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
           end
           Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
        end
    end
    
    [eigVec, eigE] = eig(Dmat - Mmat);
    eigE = -1*diag(eigE);
    
    FE = identifyFloquet(eigE, N, trajData.T);
end