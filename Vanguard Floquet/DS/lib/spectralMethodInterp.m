function [FE, eigE] = spectralMethodInterp(trajData, jacEval, N, d)
    T = trajData.T;
    [D,t] = fourierdiff(N); t = t*T/(2*pi);
    Dmat = zeros(d*N); Mmat = zeros(d*N);
    
    tt = trajData.T*trajData.fourierGrid/(2*pi);
    Xdim = size(trajData.X); Xdim = Xdim(2);
    X = zeros(N, Xdim);
    for i = 1:Xdim
        X(:,i) = interp_sinc(tt, trajData.X(:,i), t);
    end
    
    if(isfield(trajData, 'chiLinearTerm'))
        for i = 1:N
            X(i,2) = X(i,2) + trajData.chiLinearTerm*t(i);
        end
    end
    
    if(isfield(trajData, 'U'))
        Udim = size(trajData.U); Udim = Udim(2);
        U = zeros(N, Udim);
        for i = 1:Udim
            U(:,i) = interp_sinc(tt, trajData.U(:,i), t);
        end
    end
    
    for i = 1:N
        if(isfield(trajData, 'U'))
            A = jacEval(X(i,:), U(i,:), trajData);
        else
            A = jacEval(X(i,:), trajData);
        end
        
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