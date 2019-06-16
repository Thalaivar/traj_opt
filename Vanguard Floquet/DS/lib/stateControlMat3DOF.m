function [X,U,T,VR] = stateControlMat3DOF(sol, N, windShear)
% returns state and control matrices for 3DoF aircraft
% X = [V,chi,gamma,x,y,z] : nState = 6;
% U = [CL, mu, CT] : nControl = 3;
% for full state discretization, nControl = 2 (CT = 0);
    nControl = 2; nState = 6;
    X = zeros(N, nState); U = zeros(N,3);
    for i = 1:nState
        j = (i-1)*N;
        X(:,i) = sol(j+1:j+N);
    end
    for i = nState+1:nState+nControl
        j = (i-1)*N;
        U(:,i-nState) = sol(j+1:j+N);
    end
    
    T = sol(8*N+1); 
    if strcmp(windShear, 'not-same')
        VR = sol(8*N+2);
    else
        VR = [];
    end
end