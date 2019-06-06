function [c, ceq] = constraintFunc(X, trajData)
    nState = 12; nControl = 4; nAeroCoeffs = 6;
    if strcmp(trajData.shape, 'loiter')
        [trajData, psiLinearTerm] = stateControlMat(X, trajData);
    else
        trajData = stateControlMat(X, trajData);
    end
    
    D = trajData.D;
    t = trajData.fourierGrid; shape = trajData.shape;
    T = trajData.T; N = trajData.N;
    
    % spectral derivatives
    diffFac = 2*pi/T;
    spectralXdot = zeros(N*nState,1);
    for i = 1:nState
        j = (i-1)*N;
        if i == 9
            if strcmp(shape, 'loiter')
                spectralXdot(j+1:j+N,1) = diffFac*D*trajData.X(:,i) + psiLinearTerm*ones(N,1);
            else
                spectralXdot(j+1:j+N,1) = diffFac*D*trajData.X(:,i);
            end 
        else
            spectralXdot(j+1:j+N,1) = diffFac*D*trajData.X(:,i);
        end
    end
    
    % for loiter trajectories, add linear term to psi
    t = t/diffFac;
    if strcmp(shape, 'loiter')
        for i = 1:N
            trajData.X(i,9) = trajData.X(i,9) + psiLinearTerm*t(i);
        end
    end
    
    % derivatives from dynamics
    Xdot = zeros(N*nState,1); 
    ac = trajData.ac;
    aeroCoeffs = zeros(N, nAeroCoeffs);
    for i = 1:N
        Z = trajData.X(i,:); U = trajData.U(i,:);
        [dZ, coeffs] = ac.stateDervs([Z,U], trajData.VR);
        aeroCoeffs(i,:) = coeffs;
        for j = 1:nState
            Xdot(i+(j-1)*N,1) = dZ(j,1);
        end
    end
    
    % dynamics constraints
    ceq = spectralXdot - Xdot;
    
    % control input derivatives
    [d_df, dd_df] = dftDerv(trajData.U(:,1), T);
    [d_da, dd_da] = dftDerv(trajData.U(:,2), T);
    [d_de, dd_de] = dftDerv(trajData.U(:,3), T);
    [d_dr, dd_dr] = dftDerv(trajData.U(:,4), T);
    
    % get quanitites for constraints
    V = sqrt(trajData.X(:,1).^2 + trajData.X(:,2).^2 + trajData.X(:,3).^2);
    alpha = atan(trajData.X(:,3)./trajData.X(:,1));
    beta  = asin(trajData.X(:,2)./V);
    CL = aeroCoeffs(:,3);
    
    % constraints on states and controls
    c(1:N,1) = trajData.X(:,12) + 0.5*ac.b*abs(sin(trajData.X(:,7))); % wing tip clearance
    % CL constraint
    c(end+1:end+N) = CL - 1.17;
    c(end+1:end+N) = -0.2 - CL;
    % V constraint
    c(end+1:end+N) = V - 80;
    c(end+1:end+N) = 10 - V;
    % beta constraint
    c(end+1:end+N) = beta - 5*pi/180;
    c(end+1:end+N) = -5*pi/180 - beta;
    % alpha constraint
    c(end+1:end+N) = alpha - 15*pi/180;
    c(end+1:end+N) = -15*pi/180 - alpha;
    % control rate^2 constraints
    limSet = pi/60;
    c(end+1:end+N) = dd_df - limSet;
    c(end+1:end+N) = -dd_df - limSet;
    c(end+1:end+N) = dd_de - limSet;
    c(end+1:end+N) = -dd_de - limSet;
    c(end+1:end+N) = dd_da - limSet;
    c(end+1:end+N) = -dd_da - limSet;
    c(end+1:end+N) = dd_dr - limSet;
    c(end+1:end+N) = -dd_dr - limSet;
    
    % control rate constraints
%     limSet = pi/20;
%     c(end+1:end+N) = d_df - limSet;
%     c(end+1:end+N) = -d_df - limSet;
%     c(end+1:end+N) = d_de - limSet;
%     c(end+1:end+N) = -d_de - limSet;
%     c(end+1:end+N) = d_da - limSet;
%     c(end+1:end+N) = -d_da - limSet;
%     c(end+1:end+N) = d_dr - limSet;
%     c(end+1:end+N) = -d_dr - limSet;
    
    c(end+1:end+N) = 0.5*trajData.X(:,4).*ac.b./V - 0.1;
    c(end+1:end+N) = -0.1 - 0.5*trajData.X(:,4).*ac.b./V;
    
    % periodicity constraint for psi
    if strcmp(shape, 'loiter')
        ceq(end+1) = psiLinearTerm*T + 2*pi;
    end
end