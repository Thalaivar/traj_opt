function [trajData, chiLinearTerm] = stateControlMat(X, trajData)
    nState =12; nControl = 4; nOptim = nState + nControl;
    if strcmp(trajData.shape, 'loiter')
        if strcmp(trajData.type, 'not-same')
            N = (length(X)-3)/nOptim;
            chiLinearTerm = X(N*nOptim+3);
            VR = X(nOptim*N+2); trajData.VR = VR;
        else
            N = (length(X)-2)/nOptim;
            chiLinearTerm = X(N*nOptim+2);
        end
    else
        if strcmp(trajData.type, 'not-same')
            N = (length(X)-2)/nOptim;
            VR = X(nOptim*N+2); trajData.VR = VR;
        else
            N = (length(X)-1)/nOptim;
        end
    end
    
    trajData.N = N;
    T = X(nOptim*N+1); trajData.T = T;
    Z = zeros(N,nState); U = zeros(N,nControl); 
    for i = 1:nOptim
        j = (i-1)*N;
        if i <= nState
            Z(:,i) = X(j+1:j+N);
        else
            U(:,i-nState) = X(j+1:j+N);
        end
    end
    
    trajData.X = Z; trajData.U = U;
end